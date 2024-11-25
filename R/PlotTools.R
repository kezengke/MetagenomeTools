library(ggplot2)

#' Function to assign directions to log10 pvalues of t-test
processTtestRes <- function(t_results) {
  #creating log10 pvalues for t-test results
  t_results$direction<-t_results$stats/abs(t_results$stats)
  #log10 pvals
  t_results$logP<-log10(t_results$pval)
  #pval with direction
  t_results$p_directed<-t_results$direction*t_results$logP

  return(t_results)
}

#' Function to assign directions to log10 pvalues of wilcoxon
processWilcoxonRes <- function(wilcox_results) {
  #creating log10 pvalues for wilcoxon results
  wilcox_results$direction<-wilcox_results$stats/abs(wilcox_results$stats)
  #log10 pvals
  wilcox_results$logP<-log10(wilcox_results$pval)
  #pval with direction
  wilcox_results$p_directed<-wilcox_results$direction*wilcox_results$logP

  return(wilcox_results)
}

#' Function to assign directions to log10 pvalues of DESeq2
processDESeq2Res <- function(deseq_results) {
  #creating log10 pvalues for DESeq2 results
  deseq_results$direction<-deseq_results$stats/abs(deseq_results$stats)
  #log10 pvals
  deseq_results$logP<-log10(deseq_results$pval)
  #pval with direction
  deseq_results$p_directed<-deseq_results$direction*deseq_results$logP
  deseq_results$p_directed<-deseq_results$p_directed*-1

  # for deseq2 0 p-value outputs
  maxPosIndx<-which(deseq_results$direction>0 & deseq_results$pval == 0)
  minNegIndx<-which(deseq_results$direction<0 & deseq_results$pval == 0)

  maxPosP<-max(deseq_results$p_directed[deseq_results$p_directed < Inf])
  minNegP<-min(deseq_results$p_directed[deseq_results$p_directed > -Inf])

  deseq_results$p_directed[maxPosIndx]<-maxPosP + (maxPosP * 1e-6)
  deseq_results$p_directed[minNegIndx]<-minNegP + (minNegP * 1e-6)

  return(deseq_results)
}

#' Function to assign directions to log10 pvalues of edgeR
processEdgeRRes <- function(edger_results) {
  #creating log10 pvalues for edgeR results
  edger_results$direction<-edger_results$stats/abs(edger_results$stats)
  #log10 pvals
  edger_results$logP<-log10(edger_results$pval)
  #pval with direction
  edger_results$p_directed<-edger_results$direction*edger_results$logP
  edger_results$p_directed<-edger_results$p_directed*-1

  return(edger_results)
}

#' Function to assign directions to log10 pvalues of edgeR
processALDEx2Res <- function(aldex2_results) {
  #creating log10 pvalues for edgeR results
  aldex2_results$direction<-aldex2_results$stats/abs(aldex2_results$stats)
  #log10 pvals
  aldex2_results$logP<-log10(aldex2_results$pval)
  #pval with direction
  aldex2_results$p_directed<-aldex2_results$direction*aldex2_results$logP
  aldex2_results$p_directed<-aldex2_results$p_directed*-1

  return(aldex2_results)
}


#' Function to plot Log10 p-values
Log10PvPPlot <- function(results1, results2, method1, method2, title = "Default title") {
  # Plotting log10 p-values with directions based on p-value significance by which algorithm
  df <- data.frame(results1$pval, results1$p_directed, results2$pval, results2$p_directed)
  rownames(df) <- rownames(results1)
  colnames(df) <- c("results1p", "results1Dp", "results2p", "results2Dp")

  df <- df %>%
    mutate(signif = case_when(
      results1p < 0.05 & results2p < 0.05 ~ "Both",
      results1p < 0.05 & results2p >= 0.05 ~ method1,
      results1p >= 0.05 & results2p < 0.05 ~ method2,
      TRUE ~ "Neither"
    ))

  counts <- df %>%
    group_by(signif) %>%
    summarise(count = n()) %>%
    mutate(label = paste(signif, " (n=", count, ")", sep = ""))

  df <- merge(df, counts, by = "signif")

  color_category<-c("t-test" = "#c44c52",
                    "Wilcoxon" = "#c8892b",
                    "DESeq2" = "#ebd374",
                    "edgeR" = "#94c6d4",
                    "ALDEx2t-test" = "#5d879e",
                    "ALDEx2Wilcoxon" = "#4863aa",
                    "Both" = "#a24c97",
                    "Neither" = "#cdda73",

                    "t-test-Before" = "#c44c52",
                    "t-test-After" = "#94c6d4",
                    "DESeq2-Before" = "#ebd374",
                    "DESeq2-After" = "#4863aa"
                    )

  df$plotColor<-color_category[df$signif]

  p<-ggplot(df, aes(x = results1Dp, y = results2Dp)) +
    geom_point(aes(color = label), alpha = 0.4, size = 3) +
    geom_abline(intercept = 0, slope = 1, col = "gray") +
    scale_color_manual(values = setNames(df$plotColor, df$label)) +
    labs(color = "Significant by") +
    theme_classic() +
    ggtitle(title) +
    xlab(paste0("log10(", method1, " p-value)")) +
    ylab(paste0("log10(", method2, " p-value)"))

  return(p)
}

#' Function to plot taxa correlation histograms
TaxaCorHistogram <- function(table, title = "Default title", plotCol) {
  spearCor<-cor(t(table), method = "spearman")
  corVal<-spearCor[lower.tri(spearCor, diag = F)]

  p<-ggplot(data = data.frame(corVal), aes(x = corVal)) +
    geom_histogram(binwidth = 0.05, fill = plotCol, color = "black", alpha = 0.7) +
    labs(x = "Spearman Correlation",
         y = "Frequency",
         title = title) +
    theme_minimal()

  return(p)
}

#' Function to plot p-value histograms
PvalHistogram <- function(table, histoCol, title = "Default title") {
  p<-ggplot(data = table, aes(x = pval)) +
    geom_histogram(binwidth = 0.05, fill = histoCol, color = "black", alpha = 0.7) +
    labs(x = "p-value",
         y = "Frequency",
         title = title) +
    theme_minimal()

  return(p)
}

