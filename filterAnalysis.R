#!/usr/bin/Rscript

marginalCounts <- read.csv("marginalCounts.csv", comment.char="#")

threshold <- function(total_freq_min, wildtype_freq_min, diff) {
  
  marginalCounts$total.freq.r1 = marginalCounts$Q1.t1.rep1.frequency + marginalCounts$Q2.t1.rep1.frequency + marginalCounts$Q3.t1.rep1.frequency + marginalCounts$Q4.t1.rep1.frequency + marginalCounts$Q5.t1.rep1.frequency
  marginalCounts$total.freq.r2 = marginalCounts$Q1.t1.rep2.frequency + marginalCounts$Q2.t1.rep2.frequency + marginalCounts$Q3.t1.rep2.frequency + marginalCounts$Q4.t1.rep2.frequency + marginalCounts$Q5.t1.rep2.frequency
  
  #filter out fs
  marginalCounts <- marginalCounts[!grepl("fs", marginalCounts$aaChange),]
  
  #filter analysis here
  marginalCountsFilter <- marginalCounts[marginalCounts$total.freq.r1 >= total_freq_min & marginalCounts$total.freq.r2 >= total_freq_min & marginalCounts$wildtype.t1.rep1.frequency < wildtype_freq_min, ]
  
  marginalCountsFilter$weight.score.r1 = ((0.2 * marginalCountsFilter$Q1.t1.rep1.frequency) + (0.4 * marginalCountsFilter$Q2.t1.rep1.frequency) + (0.6 * marginalCountsFilter$Q3.t1.rep1.frequency) + (0.8 * marginalCountsFilter$Q4.t1.rep1.frequency) + (1 * marginalCountsFilter$Q5.t1.rep1.frequency)) /
    (marginalCountsFilter$Q1.t1.rep1.frequency + marginalCountsFilter$Q2.t1.rep1.frequency + marginalCountsFilter$Q3.t1.rep1.frequency + marginalCountsFilter$Q4.t1.rep1.frequency + marginalCountsFilter$Q5.t1.rep1.frequency)
  
  marginalCountsFilter$weight.score.r2 = ((0.2 * marginalCountsFilter$Q1.t1.rep2.frequency) + (0.4 * marginalCountsFilter$Q2.t1.rep2.frequency) + (0.6 * marginalCountsFilter$Q3.t1.rep2.frequency) + (0.8 * marginalCountsFilter$Q4.t1.rep2.frequency) + (1 * marginalCountsFilter$Q5.t1.rep2.frequency)) /
    (marginalCountsFilter$Q1.t1.rep2.frequency + marginalCountsFilter$Q2.t1.rep2.frequency + marginalCountsFilter$Q3.t1.rep2.frequency + marginalCountsFilter$Q4.t1.rep2.frequency + marginalCountsFilter$Q5.t1.rep2.frequency)
  
  marginalCountsFilter$diff <- marginalCountsFilter$weight.score.r2 - marginalCountsFilter$weight.score.r1
  #view histogram of marginal counts filtered
  
  marginalCountsFilterStop <- marginalCountsFilter[grepl("Ter", marginalCountsFilter$hgvsp),]
  
  
  marginalCountsFilterOutlier <- marginalCountsFilter[marginalCountsFilter$diff < diff & marginalCountsFilter$diff > -1 * diff,]
  #filter out low quality nonsense
  marginalCountsFilterOutlierStop <- marginalCountsFilterOutlier[grepl("Ter", marginalCountsFilterOutlier$hgvsp),]
  
  marginalCountsFilterOutlierSyn <- marginalCountsFilterOutlier[grepl("=", marginalCountsFilterOutlier$hgvsp),]
  marginalCountsFilterOutlierNonSyn <- marginalCountsFilterOutlier[!grepl("=", marginalCountsFilterOutlier$hgvsp),]
  
  
  med_stop <- median((marginalCountsFilterOutlierStop$weight.score.r1 + marginalCountsFilterOutlierStop$weight.score.r2) / 2)
  med_syn <- median((marginalCountsFilterOutlierSyn$weight.score.r1 + marginalCountsFilterOutlierSyn$weight.score.r2) / 2)
  
  q75_stop <- as.numeric(summary((marginalCountsFilterOutlierStop$weight.score.r1 + marginalCountsFilterOutlierStop$weight.score.r2) / 2)[5])
  
  marginalCountsFilterOutlier$weight.score.norm.r1 <- (marginalCountsFilterOutlier$weight.score.r1 - med_stop) / (med_syn - med_stop)
  marginalCountsFilterOutlier$weight.score.norm.r2 <- (marginalCountsFilterOutlier$weight.score.r2 - med_stop) / (med_syn - med_stop)
  
  marginalCountsFilterOutlier$weight.score.avg <- (marginalCountsFilterOutlier$weight.score.norm.r1 + marginalCountsFilterOutlier$weight.score.norm.r2) / 2
  
  #marginalCountsFilterOutlier <- marginalCountsFilter[marginalCountsFilter$diff < 2 & marginalCountsFilter$diff > -2,]
  
  out <- tapply(marginalCountsFilterOutlier$weight.score.avg, marginalCountsFilterOutlier$hgvsp, mean)
  combinedCounts <- data.frame(x=names(out), y=out, row.names=NULL)
  colnames(combinedCounts) <- c("hgvsp", "score")
  
  #get n
  for (i in 1:length(combinedCounts$hgvsp)) {
    combinedCounts$n[i] <- 2 * length(grep(combinedCounts$hgvsp[i], marginalCountsFilterOutlier$hgvsp))
  }
  
  #get sd of counts
  
  for (i in 1:length(combinedCounts$hgvsp)) {
    combinedCounts$sd[i] <- sd(c(marginalCountsFilterOutlier$weight.score.norm.r1[grep(combinedCounts$hgvsp[i], marginalCountsFilterOutlier$hgvsp)], marginalCountsFilterOutlier$weight.score.norm.r2[grep(combinedCounts$hgvsp[i], marginalCountsFilterOutlier$hgvsp)]))
  }
  
  #get standard error
  combinedCounts$se <- combinedCounts$sd / sqrt(combinedCounts$n)
  
  maveVisOutput <- data.frame(c(combinedCounts$hgvsp), combinedCounts$score, combinedCounts$sd, combinedCounts$n, combinedCounts$se)
  colnames(maveVisOutput) <- c("hgvs_pro", "score", "sd", "df", "se")
  
  #get number of variants with NA(no sd or se)
  sum(is.na(maveVisOutput$sd))
  
  maveVisOutput <- maveVisOutput[-1,]#remove the p.=
  maveVisOutput <- maveVisOutput[maveVisOutput$se < 0.9,]
  
  weight_corr <- cor(marginalCountsFilterOutlier$weight.score.norm.r1, marginalCountsFilterOutlier$weight.score.norm.r2)
  freq_corr <- cor(marginalCountsFilterOutlier$total.freq.r1, marginalCountsFilterOutlier$total.freq.r2)
  
  all_coverage <- length(maveVisOutput$hgvs_pro)
  stop_coverage <- length(maveVisOutput$hgvs_pro[grepl("Ter", maveVisOutput$hgvs_pro)])
  syn_coverage <- length(maveVisOutput$hgvs_pro[grepl("=", maveVisOutput$hgvs_pro)])

  #return(c(total_freq_min, wildtype_freq_min, med_stop, med_syn, weight_corr, freq_corr, length(maveVisOutput$score), diff, q75_stop, all_coverage, stop_coverage, syn_coverage))
  cov_stop <- sd((marginalCountsFilterOutlierStop$weight.score.r1 + marginalCountsFilterOutlierStop$weight.score.r2) / 2) / mean((marginalCountsFilterOutlierStop$weight.score.r1 + marginalCountsFilterOutlierStop$weight.score.r2) / 2)
  return(c(total_freq_min, cov_stop, med_stop))
  
}

#freq_data <- data.frame(matrix(ncol=12, nrow=0))
freq_data <- data.frame(matrix(ncol=3, nrow=0))
colnames(freq_data) <- c("total_freq_min", "cov_stop", "med_stop")
#colnames(freq_data) <- c("total_freq_min", "wildtype_freq_min", "med_stop", "med_syn", "weight_corr", "freq_corr", "n", "diff", "q75_stop", "all_coverage", "stop_coverage", "syn_coverage")


total_freq <- c(1e-7, 2e-7, 5e-7, 1e-6, 2e-6, 5e-6, 1e-5, 2e-5, 5e-5, 1e-4, 2e-4, 5e-4, 1e-3, 2e-3, 5e-3, 1e-2, 2e-2, 5e-2, 1e-1, 2e-1, 5e-1)
wildtype_freq <- c(1e-4)
diff <- c(0.25)

for (i in total_freq) {
  for (j in wildtype_freq) {
    for(k in diff) {
      try(tmp <- threshold(i, j, k))
      try(freq_data <- rbind(freq_data, tmp))
      try(cat(i, "and", j, "and", k, "complete.", "\n"))
    }
  }
}

colnames(freq_data) <- c("total_freq_min", "cov_stop", "med_stop")
plot(freq_data$total_freq_min, freq_data$cov_stop, log = "x", xlab = "Total Frequency", ylab = "Coefficient of Variation", main = "Coefficient of Variation of Nonsense Distribution")
#colnames(freq_data) <- c("total_freq_min", "wildtype_freq_min", "med_stop", "med_syn", "weight_corr", "freq_corr", "n", "diff", "q75_stop", "all_coverage", "stop_coverage", "syn_coverage")
plot(freq_data$total_freq_min, freq_data$med_stop, log = "x", xlab = "Total Frequency", ylab = "Median", main = "Median of Nonsense Distribution", type = "l")

write.csv(freq_data, file = "freq_data.csv")

