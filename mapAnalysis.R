library(dplyr)
library(plotrix)
library(PerformanceAnalytics)

#total_freq_min = 1e-05
#wildtype_freq_min = 1e-04
#diff = 0.25

total_freq_min = 2e-5
wildtype_freq_min = 1e-4
diff = 0.04
#add frequencies 
marginalCounts$total.freq.r1 = marginalCounts$Q1.t1.rep1.frequency + marginalCounts$Q2.t1.rep1.frequency + marginalCounts$Q3.t1.rep1.frequency + marginalCounts$Q4.t1.rep1.frequency + marginalCounts$Q5.t1.rep1.frequency
marginalCounts$total.freq.r2 = marginalCounts$Q1.t1.rep2.frequency + marginalCounts$Q2.t1.rep2.frequency + marginalCounts$Q3.t1.rep2.frequency + marginalCounts$Q4.t1.rep2.frequency + marginalCounts$Q5.t1.rep2.frequency

#filter out fs
marginalCounts <- marginalCounts[!grepl("fs", marginalCounts$aaChange),]

marginalCountsFilter <- marginalCounts[marginalCounts$total.freq.r1 >= total_freq_min & marginalCounts$total.freq.r2 >= total_freq_min & marginalCounts$wildtype.t1.rep1.frequency < wildtype_freq_min, ]
marginalCountsFilter <- marginalCounts[marginalCounts$total.freq.r1 >= total_freq_min & marginalCounts$total.freq.r2 >= total_freq_min & marginalCounts$wildtype.t1.rep1.frequency < wildtype_freq_min, ]


#weighted average for rep1
marginalCountsFilter$weight.score.r1 = ((0.2 * marginalCountsFilter$Q1.t1.rep1.frequency) + (0.4 * marginalCountsFilter$Q2.t1.rep1.frequency) + (0.6 * marginalCountsFilter$Q3.t1.rep1.frequency) + (0.8 * marginalCountsFilter$Q4.t1.rep1.frequency) + (1 * marginalCountsFilter$Q5.t1.rep1.frequency)) /
  (marginalCountsFilter$Q1.t1.rep1.frequency + marginalCountsFilter$Q2.t1.rep1.frequency + marginalCountsFilter$Q3.t1.rep1.frequency + marginalCountsFilter$Q4.t1.rep1.frequency + marginalCountsFilter$Q5.t1.rep1.frequency)
  
marginalCountsFilter$weight.score.r2 = ((0.2 * marginalCountsFilter$Q1.t1.rep2.frequency) + (0.4 * marginalCountsFilter$Q2.t1.rep2.frequency) + (0.6 * marginalCountsFilter$Q3.t1.rep2.frequency) + (0.8 * marginalCountsFilter$Q4.t1.rep2.frequency) + (1 * marginalCountsFilter$Q5.t1.rep2.frequency)) /
  (marginalCountsFilter$Q1.t1.rep2.frequency + marginalCountsFilter$Q2.t1.rep2.frequency + marginalCountsFilter$Q3.t1.rep2.frequency + marginalCountsFilter$Q4.t1.rep2.frequency + marginalCountsFilter$Q5.t1.rep2.frequency)

marginalCountsFilter$diff <- marginalCountsFilter$weight.score.r2 - marginalCountsFilter$weight.score.r1

for(i in 1:length(marginalCountsFilter$hgvsp)) {
  marginalCountsFilter$sd[i] <- (sd(c(marginalCountsFilter$weight.score.r1[i], marginalCountsFilter$weight.score.r2[i]))) / (length(grep(marginalCountsFilter$hgvsp[i], marginalCountsFilter$hgvsp)) * 2)
}

hist(marginalCountsFilter$sd)
#view histogram of marginal counts filtered

marginalCountsFilterStop <- marginalCountsFilter[grepl("Ter", marginalCountsFilter$hgvsp),]

#marginalCountsFilterOutlier <- marginalCountsFilter[marginalCountsFilter$diff < diff & marginalCountsFilter$diff > -1 * diff,]

marginalCountsFilterOutlier <- marginalCountsFilter[marginalCountsFilter$sd < diff,]



#filter out low quality nonsense
marginalCountsFilterOutlierStop <- marginalCountsFilterOutlier[grepl("Ter", marginalCountsFilterOutlier$hgvsp),]

marginalCountsFilterOutlierSyn <- marginalCountsFilterOutlier[grepl("=", marginalCountsFilterOutlier$hgvsp),]
marginalCountsFilterOutlierNonSyn <- marginalCountsFilterOutlier[!grepl("=", marginalCountsFilterOutlier$hgvsp),]

#filter for counts in q1 and q2

med_stop <- median((marginalCountsFilterOutlierStop$weight.score.r1 + marginalCountsFilterOutlierStop$weight.score.r2) / 2)
med_syn <- median((marginalCountsFilterOutlierSyn$weight.score.r1 + marginalCountsFilterOutlierSyn$weight.score.r2) / 2)


marginalCountsFilterOutlier$weight.score.norm.r1 <- (marginalCountsFilterOutlier$weight.score.r1 - med_stop) / (med_syn - med_stop)
marginalCountsFilterOutlier$weight.score.norm.r2 <- (marginalCountsFilterOutlier$weight.score.r2 - med_stop) / (med_syn - med_stop)



#normalize weights

#A - B
#(A-median(A)) / (median(B) - median(A))

#histograms

#marginalCountsFilter$diff <- marginalCountsFilter$weight.score.norm.r2 - marginalCountsFilter$weight.score.norm.r1

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
#sum(is.na(maveVisOutput$sd))

maveVisOutput <- maveVisOutput[-1,]#remove the p.=
maveVisOutput <- maveVisOutput[maveVisOutput$se < 0.35,]
maveVisOutput <- maveVisOutput[!grepl("delins", maveVisOutput$hgvs_pro),]
write.csv(maveVisOutput, "ve_map_sa1.csv")

colnames(maveVisOutput)[1] = "hgvsp"
filtered <- merge(maveVisOutput, marginalCountsFilterOutlier, by = "hgvsp")
plot(filtered$weight.score.norm.r1, filtered$weight.score.norm.r2)
cor(filtered$weight.score.norm.r1, filtered$weight.score.norm.r2)


#corr
cat("weights:", cor(marginalCountsFilterOutlier$weight.score.norm.r1, marginalCountsFilterOutlier$weight.score.norm.r2))
plot(marginalCountsFilterOutlier$weight.score.norm.r1, marginalCountsFilterOutlier$weight.score.norm.r2, xlab = "Weighted Average Score (Rep 1)", ylab = "Weighted Average Score (Rep 2)", main = "Weighted Average Score Replicate Correlation")
cat("total freq:", cor(marginalCountsFilterOutlier$total.freq.r1, marginalCountsFilterOutlier$total.freq.r2))

colnames(vardata)[1] = "hgvs_pro"
positionMap <- merge(maveVisOutput, vardata, by.x = "hgvs_pro")

plot(positionMap$start, positionMap$score)




#Classfifying pathogenic and benign
#marginalCountsFilterOutlierSyn <- marginalCountsFilterOutlier[grepl("=", marginalCountsFilterOutlier$hgvsp),]
#marginalCountsFilterOutlierStop <- marginalCountsFilterOutlier[grepl("Ter", marginalCountsFilterOutlier$hgvsp),]

#correlations


cat("Q1:", cor(marginalCountsFilterOutlier$Q1.t1.rep1.frequency, marginalCountsFilterOutlier$Q1.t1.rep2.frequency))
cat("Q2:", cor(marginalCountsFilterOutlier$Q2.t1.rep1.frequency, marginalCountsFilterOutlier$Q2.t1.rep2.frequency))
cat("Q3:", cor(marginalCountsFilterOutlier$Q3.t1.rep1.frequency, marginalCountsFilterOutlier$Q3.t1.rep2.frequency))
cat("Q4:", cor(marginalCountsFilterOutlier$Q4.t1.rep1.frequency, marginalCountsFilterOutlier$Q4.t1.rep2.frequency))
cat("Q5:", cor(marginalCountsFilterOutlier$Q5.t1.rep1.frequency, marginalCountsFilterOutlier$Q5.t1.rep2.frequency))

#scoring variants
marginalCountsFilterOutlierSynNorm <- marginalCountsFilterOutlier[grepl("=", marginalCountsFilterOutlier$hgvsp),]
mean(marginalCountsFilterOutlierSynNorm$weight.score.avg) + 2 * std.error(marginalCountsFilterOutlierSynNorm$weight.score.avg)
mean(marginalCountsFilterOutlierSynNorm$weight.score.avg) - 2 * std.error(marginalCountsFilterOutlierSynNorm$weight.score.avg)
#lci of syn: 0.9582131
#uci of syn: 1.021833
#mean of syn 0.9900
#std error: 0.01623

marginalCountsFilterOutlierStopNorm <- marginalCountsFilterOutlier[grepl("Ter", marginalCountsFilterOutlier$hgvsp),]
marginalCountsFilterOutlierNonSynNorm <- marginalCountsFilterOutlier[!grepl("=", marginalCountsFilterOutlier$hgvsp),]
marginalCountsFilterOutlierSynNorm <- marginalCountsFilterOutlier[grepl("=", marginalCountsFilterOutlier$hgvsp),]

#marginalCountsFilterOutlierStopNorm <- marginalCountsFilterOutlierStopNorm[marginalCountsFilterOutlierStopNorm$Q1.t1.rep1.count > 100 & marginalCountsFilterOutlierStopNorm$Q1.t1.rep2.count > 100,]

hist(marginalCountsFilterOutlierSynNorm$weight.score.avg, xlim = c(-0.5,2), xlab = "Score", main = "Stop and Synonymous Distributions (Filter = 1e-5)", breaks = 30, col = c1, freq=T)
hist(marginalCountsFilterOutlierStopNorm$weight.score.avg[marginalCountsFilterOutlierStopNorm$weight.score.avg < 1.5], xlim = c(-0.5,2), main = "Stop and Synonymous Distributions (Filter = 1e-5)", xlab = "Score", breaks = 30, add = T, col = c2, freq = T)
hist(marginalCountsFilterOutlier$weight.score.avg, main = "Non-synonymous distribution", xlab = "Score", xlim = c(-0.5,2.5), breaks = 30, freq = F)

#Variant classification
totalvariants <- maveVisOutput
totalvariants$lci <- totalvariants$score - 1.96 * totalvariants$se
totalvariants$uci <- totalvariants$score + 1.96 * totalvariants$se

pathogenic_var <- filter(totalvariants, lci > 1.021833)
pathogenic_var2 <- filter(totalvariants, uci < 0.9582131)
path_total <- rbind(pathogenic_var, pathogenic_var2)
pathfinal <- parseHGVS(path_total$hgvs_pro, aacode = 1)
pathfinal <- pathfinal[c(4,5,6)]
pathfinal <- pathfinal[grepl("*", pathfinal$variant),]
write.csv(pathfinal, file = "path_sa1.csv", quote = F, row.names = F)

#split syn and stop and nonsyn region

data_reg2 <- grepl.sub(data=marginalCountsFilterOutlier, pattern = "(?:14[6-9]|1[5-9][0-9]|2[0-7][0-9]|28[0-7])", Var = "hgvsp")
data_reg3 <- grepl.sub(data=marginalCountsFilterOutlier, pattern = "(?:28[8-9]|29[0-9]|3[0-9]{2}|40[0-9]|41[0-8])", Var = "hgvsp")

stop_reg2 <- data_reg2[grepl("Ter", data_reg2$hgvsp),]
stop_reg3 <- data_reg3[grepl("Ter", data_reg3$hgvsp),]

syn_reg2 <- data_reg2[grepl("=", data_reg2$hgvsp),]
syn_reg3 <- data_reg3[grepl("=", data_reg3$hgvsp),]

#reg2-rep1
hist(stop_reg2$weight.score.norm.r1[stop_reg2$weight.score.norm.r1 < 0.9], ylim = c(0, 20), breaks = 30, col = c2, main = "Reg 2 Rep 1 Syn and Stop", xlab = "Score", xlim = c(-0.5,2))
hist(syn_reg2$weight.score.norm.r1, add = T, ylim = c(0,20), breaks = 30, col = c1)

#reg2-rep2
hist(stop_reg2$weight.score.norm.r2[stop_reg2$weight.score.norm.r2 < 0.9], ylim = c(0, 20), breaks = 30, col = c2, xlim = c(-0.5,2), main = "Reg 2 Rep 2 Syn and Stop", xlab = "Score")
hist(syn_reg2$weight.score.norm.r2, add = T, ylim = c(0,20), breaks = 30, col = c1, xlim = c(-0.5,2))

#reg3-rep1-
hist(stop_reg3$weight.score.norm.r1[stop_reg3$weight.score.norm.r1 < 0.9], ylim = c(0, 25), breaks = 30, col = c2, xlim = c(-0.5,2), main = "Reg 3 Rep 1 Syn and Stop", xlab = "Score")
hist(syn_reg3$weight.score.norm.r1, add = T, ylim = c(0,25), breaks = 30, col = c1, xlim = c(-0.5,2))

#reg3-rep2
hist(stop_reg3$weight.score.norm.r2[stop_reg3$weight.score.norm.r2 < 0.9], ylim = c(0, 25), breaks = 30, col = c2, xlim = c(-0.5,2), main = "Reg 3 Rep 2 Syn and Stop", xlab = "Score")
hist(syn_reg3$weight.score.norm.r2, add = T, ylim = c(0,25), breaks = 35, col = c1, xlim = c(-0.5,2))



#map stop codons to position and store
stop_pos <- data.frame(c(marginalCountsFilterOutlierStopNorm$hgvsp), marginalCountsFilterOutlierStopNorm$weight.score.avg)
colnames(stop_pos) <- c("hgvsp", "score")

stop_pos$hgvsp <- parseHGVS(stop_pos$hgvsp, aacode=1)[[4]]
plot(stop_pos$hgvsp, stop_pos$score)





#replicate filter

rep_data <- data.frame(unique(marginalCountsFilter$hgvsp))
colnames(rep_data) = "hgvsp"

for (i in 1:length(rep_data$hgvsp)) {
  data[i] <- c(marginalCountsFilter$weight.score.r1[marginalCountsFilter$hgvsp == rep_data$hgvsp[i]], marginalCountsFilter$weight.score.r2[marginalCountsFilter$hgvsp == rep_data$hgvsp[i]])
}








