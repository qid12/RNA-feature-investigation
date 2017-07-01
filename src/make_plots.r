#' @param exp.dat expression data, with each row a gene and each column a sample
#' @param c.mat a matrix recoding genes of each set
#' @param sample_label a named vector, specifying samples needed and their label
#' 
#' @return a matrix, with each row a geneset and each column a sample
#'   
calculate_geneset_mean <- function(exp.dat, c.mat, sample_label) {
  new.dat = exp.dat[,sample_label]
  mean.mat <- matrix(NA, nrow = dim(c.mat)[1], ncol = dim(new.dat)[2])
  for (i in 1:dim(c.mat)[1]) {
    tmp.mat = new.dat[which(c.mat[i,]==1), ]
    mean.mat[i, ] = unlist(apply(tmp.mat, FUN = mean, MARGIN = 2, na.rm = T))
  }
  dimnames(mean.mat) = list(rownames(c.mat), names(sample_label))
  mean.mat
}

#---------------------------------------------

draw_geneset_info <- function(geneset_name, c.mat, exp.dat, output_dir, size, sample_label) {
  tmp.dat = exp.dat[c.mat[geneset_name, ]==1, ]
  cor.mat = cor(t(tmp.dat), method = "spearman", use = "pairwise.complete.obs")
  pheatmap(cor.mat, filename = paste(output_dir, geneset_name, "_cor.png", sep = ""))
  
  pheatmap(tmp.dat, annotation_col = sample_label, show_colnames = F, 
           filename = paste(output_dir, geneset_name, "_exp.png", sep = ""))
  
  random.mat <- generate_random_sets(dim(exp.dat)[2], sum(c.mat[geneset_name, ]==1), size)
  random.dat = calculate_geneset_mean(exp.dat, random.mat, colnames(exp.dat))
  geneset.dat = unlist(apply(exp.dat[which(c.mat[geneset_name,]==1), ], FUN = mean, MARGIN = 2, na.rm = T))
  background.df = melt(random.dat)
  background.df$label = paste("random sample", background.df$Var2, sep = "")
  sample.df = data.frame(value = geneset.dat, label = unique(background.df$label))
  p = ggplot(background.df, aes(label, value))
  p + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6)) + 
    geom_point(aes(x = label, colour = "red"), data = sample.df) 
  ggsave(paste(output_dir, geneset_name, "_mean.png", sep = ""), plot = last_plot(),  width = 10, height = 4.5)
  
}

#----------------------------------------------------------

generate_coexpression_network <- function(mean.mat, cutoff) {
  nodes = rownames(mean.mat)
  cor.mat = cor(t(mean.mat))
  cor.filter.mat = (cor.mat >= cutoff)
  linkage_table = c()
  for (i in 1:(dim(cor.mat)[2]-1)) {
    for (j in (i+1):dim(cor.mat)[2]) {
      if (cor.filter.mat[i,j]) {
        linkage_table = rbind(linkage_table, c(nodes[i],nodes[j],cor.mat[i,j]))
      }
    }
  }
  linkage_table
}


