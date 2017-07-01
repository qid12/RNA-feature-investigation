#' @param genesets list of genesets
#' @param all_genes vector of gene names, corresponding to genes in the expression data
#' 
#' @return a binary matrix recoding genes of each set; there might be missing gene due to data specification 
#' 
prepare_Cmat <- function(genesets, all_genes) {
  C.matrix = matrix(NA, nrow = length(genesets), ncol = length(all_genes))
  for (i in 1:length(genesets)) {
    C.matrix[i,] = all_genes %in% genesets[[i]]
  }
  C.matrix
}

#--------------------------------------------------------
#' @param N number of total genes
#' @param k size of required genesets
#' @param pool_size number of samples
#' 
#' @return a binary matrix recoding random gene sets
#' 
generate_random_sets <- function(N, k, pool_size) {
  randomC.mat <- matrix(0, nrow = pool_size, ncol = N) 
  for (i in 1:pool_size) {
    randomC.mat[i, sample(1:N, size = k, replace = F)] = 1
  }
  randomC.mat
}

#-----------------------------------------------------

local_stats_mean <- function(exp.dat, c.mat) {
  mean.mat <- matrix(NA, nrow = dim(c.mat)[1], ncol = dim(exp.dat)[2])
  for (i in 1:dim(c.mat)[1]) {
    tmp.mat = exp.dat[which(c.mat[i,]==1), ]
    mean.mat[i, ] = unlist(apply(tmp.mat, FUN = mean, MARGIN = 2, na.rm = T))
  }
  unlist(apply(mean.mat, FUN = sd, MARGIN = 1, na.rm = T))
}

#-----------------------------------------------------

local_stats_cor <- function(exp.dat, c.mat, ref_matrix) {
  cor.mat <- matrix(NA, nrow = dim(c.mat)[1], ncol = 1)
  n = sum(c.mat[1,]==1)
  for (i in 1:dim(c.mat)[1]) {
    tmp.mat = exp.dat[which(c.mat[i,]==1), ]
    tmp.cor = cor(t(tmp.mat), method = "spearman", use = "pairwise.complete.obs")
    cor.mat[i] = sum(abs(as.vector(tmp.cor * ref_matrix[1:dim(tmp.cor)[1], 1:dim(tmp.cor)[2]])))/(n*n-n)
  }
  
  as.vector(cor.mat)
}

#----------------------------------------------------

generate_lower_matrix <- function(k) {
  new.mat = rep(0, times= k*k); dim(new.mat) = c(k,k)
  for (i in 2:k) {
    j = 1
    while (j < i) {
      new.mat[i,j] = 1
      j = j + 1
    }
  }
  new.mat
}

#----------------------------------------------------

pattern_filter <- function(df, method, alpha, p.adjust.method) {
  random_mean = df$mean_stat[df$label == "random sample"]
  random_cor = df$cor_stat[df$label == "random sample"]
  df = df[df$label == "real sample",]
  df$mean_less = NA
  df$mean_greater = NA
  df$cor_less = NA
  df$cor_greater = NA
  
  if (method == "wilcox") {
    for (m in 1:(dim(expC.mat)[1]-pop_size)) {
      test1 = wilcox.test(df$mean_stat[m], random_mean, alternative = "less")
      test2 = wilcox.test(df$mean_stat[m], random_mean, alternative = "greater")
      test3 = wilcox.test(df$cor_stat[m], random_cor, alternative = "less")
      test4 = wilcox.test(df$cor_stat[m], random_cor, alternative = "greater")
      df$mean_less[m] = test1$p.value
      df$mean_greater[m] = test2$p.value
      df$cor_less[m] = test3$p.value
      df$cor_greater[m] = test4$p.value
    }
    df$mean_less = p.adjust(df$mean_less, p.adjust.method)
    df$mean_greater = p.adjust(df$mean_greater, p.adjust.method)
    df$cor_less = p.adjust(df$cor_less, p.adjust.method)
    df$cor_greater = p.adjust(df$cor_greater, p.adjust.method)
    filter1 = df$mean_less <= alpha & df$cor_less <= alpha
    filter2 = df$mean_greater <= alpha & df$cor_greater <= alpha
  } else if (method == "t.test") {
    for (m in 1:(dim(expC.mat)[1]-pop_size)) {
      test1 = t.test(log(random_mean), mu = log(df$mean_stat[m]), alternative = "less")
      test2 = t.test(log(random_mean), mu = log(df$mean_stat[m]), alternative = "greater")
      test3 = t.test(log(random_cor), mu = log(df$cor_stat[m]), alternative = "less")
      test4 = t.test(log(random_cor), mu = log(df$cor_stat[m]), alternative = "greater")
      df$mean_less[m] = test1$p.value
      df$mean_greater[m] = test2$p.value
      df$cor_less[m] = test3$p.value
      df$cor_greater[m] = test4$p.value
    }
    df$mean_less = p.adjust(df$mean_less, p.adjust.method)
    df$mean_greater = p.adjust(df$mean_greater, p.adjust.method)
    df$cor_less = p.adjust(df$cor_less, p.adjust.method)
    df$cor_greater = p.adjust(df$cor_greater, p.adjust.method)
    
    filter1 = df$mean_less <= alpha & df$cor_less <= alpha
    filter2 = df$mean_greater <= alpha & df$cor_greater <= alpha
  }
  as.numeric(filter2)*2 + as.numeric(filter1)
}



