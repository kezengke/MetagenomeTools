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

#' Function to assign color code to different algorithm combo
get_combination <- function(method1, method2) {
  # Create a sorted vector to ignore the order of input methods
  methods <- sort(c(method1, method2))

  # Define all possible combinations and their respective results
  combinations <- list(
    `DESeq2 t-test` = c("turquoise1", "purple", "green", "red"),
    `t-test Wilcoxon` = c("turquoise1", "green", "red", "tan2"),
    `edgeR t-test` = c("turquoise1", "cornflowerblue", "green", "red"),
    `DESeq2 Wilcoxon` = c("turquoise1", "purple", "green", "tan2"),
    `DESeq2 edgeR` = c("turquoise1", "purple", "cornflowerblue", "green"),
    `edgeR Wilcoxon` = c("turquoise1", "cornflowerblue", "green", "tan2")
  )

  # Check for the combination and return the result
  for (combo in names(combinations)) {
    if (all(methods == unlist(strsplit(combo, " ")))) {
      return(combinations[[combo]])
    }
  }

  return("Invalid combination")
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
  df$signif <- factor(df$signif, levels = sort(unique(df$signif)))

  color_combo<-get_combination(method1, method2)

  p <- ggplot(df, aes(x = results1Dp, y = results2Dp, color = signif)) +
    geom_point(alpha = 0.3) +
    geom_abline(intercept = 0, slope = 1, col = "gray") +
    scale_color_manual(values = color_combo,
                       labels = sort(counts$label)) +
    labs(color = "Significant by") +
    theme_classic() +
    ggtitle(title) +
    xlab(paste0("log10(", method1, " p-value)")) +
    ylab(paste0("log10(", method2, " p-value)"))

  return(p)
}

