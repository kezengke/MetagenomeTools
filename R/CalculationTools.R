library(DESeq2)
library(edgeR)
library(coin)

#' Load in counts table
LoadCountsT <- function(filePath) {
  if (grepl("\\.txt$", filePath, ignore.case = TRUE)) {
    countsT<-read.table(filePath, sep = "\t", header = T, row.names = 1, check.names = F)
    return(countsT)
  } else if (grepl("\\.csv$", filename, ignore.case = TRUE)) {
    countsT<-read.csv(file, header = T, row.names = 1, check.names = F)
    return(countsT)
  } else {
    stop("unsupported filetype")
  }
}

#' Load in metadata
LoadMeta <- function(filePath) {
  if (grepl("\\.txt$", filePath, ignore.case = TRUE)) {
    meta<-read.table(filePath, sep = "\t", header = T, row.names = 1, check.names = F)
    return(meta)
  } else if (grepl("\\.csv$", filename, ignore.case = TRUE)) {
    meta<-read.csv(file, header = T, row.names = 1, check.names = F)
    return(meta)
  } else {
    stop("unsupported filetype")
  }
}

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
safe_t_test <- function(x, y) {
  if (sd(x) == 0 || sd(y) == 0) {
    return(list(p.value = NA, estimate = NA, warning = "Data are essentially constant"))
  } else {
    return(t.test(x, y))
  }
}

calcTtest <- function(newT, meta) {
  unique_groups <- unique(meta)
  group1 <- newT[, meta == unique_groups[1], drop = FALSE]
  group2 <- newT[, meta == unique_groups[2], drop = FALSE]

  t_results <- apply(newT, 1, function(row) {
    x <- row[meta == unique_groups[1]]
    y <- row[meta == unique_groups[2]]
    test_result <- safe_t_test(x, y)
    return(c(p.value = test_result$p.value, estimate = test_result$estimate, warning = test_result$warning))
  })

  return(t_results)
}

# calcTtest <- function(table, meta) {
#   t_stats<-apply(table, 1, function(x){t.test(unlist(x)~meta$conditions)$stat})
#   t_test_p<-apply(table, 1, function(x){t.test(unlist(x)~meta$conditions)$p.value})
#
#   t_results<-cbind(t_stats, t_test_p)
#   rownames(t_results)<-rownames(table)
#   colnames(t_results)<-c("stats", "pval")
#   t_results<-data.frame(t_results, check.names = F)
#   return(t_results)
# }

#' Function to calculate Wilcoxon results
calcWilcox <- function(table, meta) {
  wilcox_stats<-apply(table, 1, function(x){statistic(wilcox_test(unlist(x)~factor(meta$conditions)))})
  wilcox_p<-apply(table, 1, function(x){pvalue(wilcox_test(unlist(x)~factor(meta$conditions)))})

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

#' Function to calculate edgeR results
calcEdgeR <- function(table, meta) {
  group <- meta$condition
  dgList <- DGEList(counts=table, group = group)
  dgList <- calcNormFactors(dgList, method="TMM")
  dgList <- estimateDisp(dgList)
  et <- exactTest(dgList)
  res <- et$table
  edger_results<-cbind(res$logFC, res$PValue)

  rownames(edger_results)<-rownames(table)
  colnames(edger_results)<-c("stats", "pval")
  edger_results<-data.frame(edger_results, check.names = F)
  return(edger_results)
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

#' Function to resample counts table with multiple (Rnorm)
resampleRNORM <- function(table, meta, multiple) {
  if (nrow(table) == 0 || nrow(meta) == 0) {
    stop("Input table or meta data frame is empty.")
  }

  groups <- unique(meta$conditions)
  if (length(groups) != 2) {
    stop("The function currently supports exactly two groups.")
  }

  group1 <- groups[1]
  group2 <- groups[2]

  # Function to calculate mean and sd for each group
  calculateMeanSd <- function(z) {
    g1 <- unlist(z[meta$conditions == group1])
    g2 <- unlist(z[meta$conditions == group2])
    c(mean1 = mean(g1), mean2 = mean(g2), sd1 = sd(g1), sd2 = sd(g2))
  }

  # Apply the function to each row of the table and combine results into a data frame
  MeanSd_table <- t(apply(table, 1, calculateMeanSd))
  MeanSd_table <- as.data.frame(MeanSd_table)
  rownames(MeanSd_table) <- rownames(table)

  # Function to generate resampled counts for each row
  resample_counts <- function(row_index) {
    z <- table[row_index,]
    g1 <- unlist(z[meta$conditions == group1])
    g2 <- unlist(z[meta$conditions == group2])

    ng1 <- rnorm(n = length(g1), mean = MeanSd_table[row_index, "mean1"], sd = sqrt(multiple * (MeanSd_table[row_index, "sd1"])^2))
    ng2 <- rnorm(n = length(g2), mean = MeanSd_table[row_index, "mean2"], sd = sqrt(multiple * (MeanSd_table[row_index, "sd2"])^2))

    c(ng1, ng2)
  }

  # Apply the resampling function to each row
  newT <- t(sapply(seq_len(nrow(table)), resample_counts))

  # Set column and row names
  colnames(newT) <- c(rownames(meta)[meta$conditions == group1], rownames(meta)[meta$conditions == group2])
  newT <- newT[, colnames(table)]
  rownames(newT) <- rownames(table)

  # Replace negative values with 0 and round to integer
  newT[newT < 0] <- 0
  newT <- data.frame(round(newT, digits = 0), check.names = F)

  return(newT)
}

#' Function to put original taxon back in resampled countsT and calculate new stats (t-test)
BackTestTtest <- function(oldT, newT, meta) {
  tRES<-c()
  for (i in 1:nrow(oldT)) {
    testing_taxa<-oldT[i, ]
    testing_frame<-newT
    #placing the taxon needs to be tested into the old counts table
    testing_frame[i, ]<-testing_taxa

    testing_frame<-normFun(testing_frame)
    res<-calcTtest(testing_frame[i, ], meta)

    tRES<-rbind(tRES, res)
  }

  rownames(tRES)<-rownames(oldT)
  tRES<-data.frame(tRES)
  return(tRES)
}

#' Function to put original taxon back in resampled countsT and calculate new stats (DESeq2)
BackTestDESeq2 <- function(oldT, newT, meta) {
  dRES<-c()
  for (i in 1:nrow(oldT)) {
    testing_taxa<-oldT[i, ]
    testing_frame<-newT
    #placing the taxon needs to be tested into the old counts table
    testing_frame[i, ]<-testing_taxa

    out<-calcDESeq2(testing_frame, meta)
    res<-out[i, ]
    dRES<-rbind(dRES, res)
  }

  rownames(dRES)<-rownames(oldT)
  dRES<-data.frame(dRES)
  return(dRES)
}
