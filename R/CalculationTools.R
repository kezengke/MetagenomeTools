
#' Function to normalize counts table
normFun <- function(table) {
  n<-colSums(table)
  sumx<-sum(table)
  for (j in 1:ncol(table)) {
    table[,j]<-table[,j]/n[j]
  }
  table<-log10(table*(sumx/ncol(table))+1)
  table<-data.frame(table, check.names = F)
  return(table)
}

#' Function to calculate t-test results
calcTtest <- function(table, meta) {
  t_stats<-apply(table, 1, function(x){t.test(unlist(x)~meta$conditions)$stat})
  t_test_p<-apply(table, 1, function(x){t.test(unlist(x)~meta$conditions)$p.value})

  t_results<-cbind(t_stats, t_test_p)
  rownames(t_results)<-rownames(table)
  colnames(t_results)<-c("stats", "pval")
  t_results<-data.frame(t_results, check.names = F)
  return(t_results)
}

#' Function to calculate Wilcoxon results
calcWilcox <- function(table, meta) {
  wilcox_stats<-apply(table, 1, function(x){wilcox.test(unlist(x)~meta$conditions)$stat})
  wilcox_p<-apply(table, 1, function(x){wilcox.test(unlist(x)~meta$conditions)$p.value})

  wilcox_results<-cbind(wilcox_stats, wilcox_p)
  rownames(wilcox_results)<-rownames(table)
  colnames(wilcox_results)<-c("stats", "pval")
  wilcox_results<-data.frame(wilcox_results, check.names = F)
  return(wilcox_results)
}

#' Function to calculate DESeq2 results
calcDESeq2 <- function(table, meta) {
  #solve deseq2 all 0 issue
  table<-table+1

  meta$conditions<-factor(meta$conditions)
  dds1 <- DESeqDataSetFromMatrix(countData=table,
                                 colData=meta,
                                 design=~conditions)
  dds2 <- DESeq(dds1)
  res <- results(dds2, cooksCutoff=FALSE, independentFiltering=FALSE)
  deseq_results<-cbind(res$stat, res$pvalue)

  rownames(deseq_results)<-rownames(table)
  colnames(deseq_results)<-c("stats", "pval")
  deseq_results<-data.frame(deseq_results, check.names = F)
  return(deseq_results)
}

#' Function to generate Log10 p-values and assign test statistics directions
directPFun<-function(ttest, deseq){
  ttest$direction<-ttest$stats/abs(ttest$stats)
  ttest$logP<-log10(ttest$pval)
  ttest$p_directed<-ttest$direction*ttest$logP
  deseq$direction<-deseq$stats/abs(deseq$stats)
  deseq$logP<-log10(deseq$pval)
  deseq$p_directed<-deseq$direction*deseq$logP
  deseq$p_directed<-deseq$p_directed*-1
  directP<-cbind(ttest$p_directed, deseq$p_directed)
  directP<-data.frame(directP, check.names = F)
  # directP<-directP[!apply(directP, 1, function(row) any(is.infinite(row))), ]
  return(directP)
}
