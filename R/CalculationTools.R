library(DESeq2)
library(edgeR)
library(coin)
library(ALDEx2)

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

#' Function to calculate t-test results (for resampling down)
calcTtest2 <- function(table, meta) {
  t_results <- apply(table, 1, function(x) {
    group1 <- x[meta$conditions == unique(meta$conditions)[1]]
    group2 <- x[meta$conditions == unique(meta$conditions)[2]]

    sd_group1 <- sd(group1)
    sd_group2 <- sd(group2)

    if (sd_group1 == 0 | sd_group2 == 0) {
      return(c(NA, NA))
    } else {
      test_result <- t.test(unlist(x) ~ meta$conditions)
      return(c(test_result$statistic, test_result$p.value))
    }
  })

  t_results <- t(t_results)
  rownames(t_results) <- rownames(table)
  colnames(t_results) <- c("stats", "pval")
  t_results <- data.frame(t_results, check.names = FALSE)

  return(t_results)
}

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

  # Wrap DESeq in tryCatch to handle errors
  dds2 <- tryCatch({
    DESeq(dds1)
  }, error = function(e) {
    # If error occurs, return NA instead of running DESeq
    return(NULL)
  })

  # If dds2 is NULL, skip the remaining lines and create deseq_results as NA
  if (is.null(dds2)) {
    deseq_results <- matrix(NA, nrow = nrow(table), ncol = 2)
  } else {
    res <- results(dds2, cooksCutoff=FALSE, independentFiltering=FALSE)
    deseq_results<-cbind(res$stat, res$pvalue)
  }

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

  # Wrap estimateDisp in tryCatch to handle errors
  dgList <- tryCatch({
    estimateDisp(dgList)
  }, error = function(e) {
    # If error occurs, return NA instead of running estimateDisp
    return(NULL)
  })

  # If dgList is NULL, skip the remaining lines and create edger_results as NA
  if (is.null(dgList)) {
    edger_results <- matrix(NA, nrow = nrow(table), ncol = 2)
  } else {
    et <- exactTest(dgList)
    res <- et$table
    edger_results <- cbind(res$logFC, res$PValue)
  }

  rownames(edger_results) <- rownames(table)
  colnames(edger_results) <- c("stats", "pval")
  edger_results <- data.frame(edger_results, check.names = F)

  return(edger_results)
}

#' Function to calculate ALDEx2 t-test results
calcALDEx2Ttest <- function(table, meta) {
  CLR<-aldex.clr(table, meta$conditions, mc.samples = 128, denom = "all", verbose = F, gamma = 0.5)
  results<-aldex.ttest(CLR, hist.plot = F, paired.test = F, verbose = F)
  effectSizes<-aldex.effect(CLR)

  aldex2T_results<-cbind(effectSizes$diff.btw, results$we.ep)

  rownames(aldex2T_results)<-rownames(table)
  colnames(aldex2T_results)<-c("stats", "pval")
  aldex2T_results<-data.frame(aldex2T_results, check.names = F)
  return(aldex2T_results)
}

#' Function to calculate ALDEx2 wilcoxon results
calcALDEx2Wilcoxon <- function(table, meta) {
  CLR<-aldex.clr(table, meta$conditions, mc.samples = 128, denom = "all", verbose = F, gamma = 0.5)
  results<-aldex.ttest(CLR, hist.plot = F, paired.test = F, verbose = F)
  effectSizes<-aldex.effect(CLR)

  aldex2W_results<-cbind(effectSizes$diff.btw, results$wi.ep)

  rownames(aldex2W_results)<-rownames(table)
  colnames(aldex2W_results)<-c("stats", "pval")
  aldex2W_results<-data.frame(aldex2W_results, check.names = F)
  return(aldex2W_results)
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

  # for deseq2 0 p-value outputs
  maxPosIndx<-which(deseq$direction>0 & deseq$pval == 0)
  minNegIndx<-which(deseq$direction<0 & deseq$pval == 0)

  maxPosP<-max(deseq$p_directed[deseq$p_directed < Inf])
  minNegP<-min(deseq$p_directed[deseq$p_directed > -Inf])

  deseq$p_directed[maxPosIndx]<-maxPosP + (maxPosP * 1e-6)
  deseq$p_directed[minNegIndx]<-minNegP + (minNegP * 1e-6)

  directP<-cbind(ttest$p_directed, deseq$p_directed)
  directP<-data.frame(directP, check.names = F)
  # directP<-directP[!apply(directP, 1, function(row) any(is.infinite(row))), ]
  return(directP)
}

#' Function to calculate Mean and Std
calcMeanStd<- function(table, meta) {
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

  return(MeanSd_table)
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

    newZ<-c(ng1, ng2)
    newZ<-newZ + abs(min(newZ))
    return(newZ)
  }

  # Apply the resampling function to each row
  newT <- t(sapply(seq_len(nrow(table)), resample_counts))

  # Set column and row names
  colnames(newT) <- c(rownames(meta)[meta$conditions == group1], rownames(meta)[meta$conditions == group2])
  newT <- newT[, colnames(table)]
  rownames(newT) <- rownames(table)

  # # Replace negative values with 0 and round to integer
  # newT[newT < 0] <- 0
  # Add smallest number to whole counts table
  # newT<-newT + abs(min(newT))
  newT <- data.frame(round(newT, digits = 0), check.names = F)

  return(newT)
}

#' Function to resample counts table (whole taxon) with multiple (Rnorm)
resampleWholeTaxonRNORM <- function(table, meta, multiple) {
  if (nrow(table) == 0 || nrow(meta) == 0) {
    stop("Input table or meta data frame is empty.")
  }

  # Function to calculate mean and sd for each group
  calculateMeanSd <- function(z) {
    c(meanTotal = mean(z), sdTotal = sd(z))
  }

  # Apply the function to each row of the table and combine results into a data frame
  MeanSd_table <- t(apply(table, 1, calculateMeanSd))
  MeanSd_table <- as.data.frame(MeanSd_table)
  rownames(MeanSd_table) <- rownames(table)

  # Function to generate resampled counts for each row
  resample_counts <- function(row_index) {
    z <- table[row_index,]

    nz <- rnorm(n = length(z), mean = MeanSd_table[row_index, "meanTotal"], sd = sqrt(multiple * (MeanSd_table[row_index, "sdTotal"])^2))

    nz<-nz+abs(min(nz))

    return(nz)
  }

  # Apply the resampling function to each row
  newT <- t(sapply(seq_len(nrow(table)), resample_counts))

  # Set column and row names
  colnames(newT) <- colnames(table)
  newT <- newT[, colnames(table)]
  rownames(newT) <- rownames(table)

  # # Replace negative values with 0 and round to integer
  # newT[newT < 0] <- 0
  # Add smallest number to whole counts table
  # newT<-newT + abs(min(newT))
  newT <- data.frame(round(newT, digits = 0), check.names = F)

  return(newT)
}

#' Function to resample counts table with multiple (Rnbinom)
resampleRNBINOM <- function(table, meta, multiple) {
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
  calculateMeanVar <- function(z) {
    g1 <- unlist(z[meta$conditions == group1])
    g2 <- unlist(z[meta$conditions == group2])
    c(mean1 = mean(g1), mean2 = mean(g2), var1 = var(g1), var2 = var(g2))
  }

  # Apply the function to each row of the table and combine results into a data frame
  MeanVar_table <- t(apply(table, 1, calculateMeanVar))
  MeanVar_table <- as.data.frame(MeanVar_table)
  rownames(MeanVar_table) <- rownames(table)

  exclude<-which(MeanVar_table$var1 <= MeanVar_table$mean1 | MeanVar_table$var2 <= MeanVar_table$mean2)

  # Function to generate resampled counts for each row
  resample_counts <- function(row_index) {
    z <- table[row_index,]
    g1 <- unlist(z[meta$conditions == group1])
    g2 <- unlist(z[meta$conditions == group2])

    n1<-length(g1)
    n2<-length(g2)

    r1<-((MeanVar_table[row_index, 1])^2)/((MeanVar_table[row_index, 3])-MeanVar_table[row_index, 1])
    p1<-r1/(r1+(MeanVar_table[row_index, 1]))

    r2<-((MeanVar_table[row_index, 2])^2)/((MeanVar_table[row_index, 4])-MeanVar_table[row_index, 2])
    p2<-r2/(r2+(MeanVar_table[row_index, 2]))

    ng1<-rnbinom(n = n1, size = r1,  prob = p1)
    ng2<-rnbinom(n = n2, size = r2,  prob = p2)

    newZ<-c(ng1, ng2)
    return(newZ)
  }

  if(length(exclude)>0){
    table<-table[-exclude, , drop = F]
    MeanVar_table<-MeanVar_table[-exclude, , drop = F]
  }

  # Apply the resampling function to each row
  newT <- t(sapply(seq_len(nrow(table)), resample_counts))

  newT[is.na(newT)] <- 0

  # Set column and row names
  colnames(newT) <- c(rownames(meta)[meta$conditions == group1], rownames(meta)[meta$conditions == group2])
  newT <- newT[, colnames(table)]
  rownames(newT) <- rownames(table)

  newT <- data.frame(round(newT, digits = 0), check.names = F)

  return(newT)
}

#' Function to resample counts table (whole taxon) with multiple (Rnbinom)
resampleWholeTaxonRNBINOM <- function(table, meta, multiple) {
  if (nrow(table) == 0 || nrow(meta) == 0) {
    stop("Input table or meta data frame is empty.")
  }

  # Function to calculate mean and sd for each group
  calculateMeanVar <- function(z) {
    c(meanTotal = mean(z), varTotal = var(z))
  }

  # Apply the function to each row of the table and combine results into a data frame
  MeanVar_table <- t(apply(table, 1, calculateMeanVar))
  MeanVar_table <- as.data.frame(MeanVar_table)
  rownames(MeanVar_table) <- rownames(table)

  exclude<-which(MeanVar_table$varTotal <= MeanVar_table$meanTotal)
  # Function to generate resampled counts for each row
  resample_counts <- function(row_index) {
    z <- table[row_index,]

    n<-length(z)

    r<-((MeanVar_table[row_index, 1])^2)/((MeanVar_table[row_index, 2])-MeanVar_table[row_index, 1])
    p<-r/(r+(MeanVar_table[row_index, 1]))

    nz<-rnbinom(n = n, size = r,  prob = p)

    return(nz)
  }

  if(length(exclude)>0){
    table<-table[-exclude, , drop = F]
    MeanVar_table<-MeanVar_table[-exclude, , drop = F]
  }

  # Apply the resampling function to each row
  newT <- t(sapply(seq_len(nrow(table)), resample_counts))

  # Set column and row names
  colnames(newT) <- colnames(table)
  newT <- newT[, colnames(table)]
  rownames(newT) <- rownames(table)

  # # Replace negative values with 0 and round to integer
  # newT[newT < 0] <- 0
  # Add smallest number to whole counts table
  # newT<-newT + abs(min(newT))
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

#' Function to calculate KS test results from uniform distribution
CalcKSpval <- function(pvalues) {
  ks_result<-ks.test(pvalues, "punif")

  format(ks_result$p.value, scientific = TRUE)
}

#' Function to calculate Chi-squared test results from uniform distribution
CalcChiSquaredPval <- function(pvalues) {
  numBins<-ceiling(1/0.05)
  breaks<-seq(0, 1, 0.05)
  observedCounts<-hist(pvalues, breaks, plot = F)$counts
  totalPvals<-length(pvalues)
  expectedCounts<-rep(totalPvals/numBins, numBins)
  pval<-chisq.test(observedCounts, p = rep(1/numBins, numBins))$p.value

  return(pval)
}
