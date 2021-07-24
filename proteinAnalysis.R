running_avg <- data.frame(matrix(ncol = 2))
colnames(running_avg) <- c("pos", "score")


maveVisOutputNoStop <- maveVisOutput[!grepl("=", maveVisOutput$hgvs_pro),]
maveVisOutputNoStop <- maveVisOutputNoStop[!grepl("Ter", maveVisOutputNoStop$hgvs_pro),]

colnames(maveVisOutputNoStop)[1] = "hgvsp"

scores <- c()
for (i in 1:273) {
  scores[i] <- print(mean(maveVisOutputNoStop$score[grepl(i + 145, maveVisOutputNoStop$hgvsp)]))
  #scores[i] <- print(summary(maveVisOutput$score[grepl(i + 145, maveVisOutput$hgvsp)])[5][[1]])
}

plot(1:273, scores, type = "l")


plot(rollmean(scores, k = 30), xlab = "Position", ylab = "Mean", main = "Moving Average of Variant Map Score (Window = 30)", type = "l")

aggregate <- c()

for (i in 1:273) {
  aggregate[i] <- print(sum((maveVisOutputNoStop$score[grepl(i + 145, maveVisOutputNoStop$hgvsp)] > 1.4)) / length ((maveVisOutputNoStop$score[grepl(i + 145, maveVisOutputNoStop$hgvsp)])))
}


misfold <- c()

for (i in 1:273) {
  misfold[i] <- print(sum((maveVisOutputNoStop$score[grepl(i + 145, maveVisOutputNoStop$hgvsp)] < 0.6)) / length ((maveVisOutputNoStop$score[grepl(i + 145, maveVisOutputNoStop$hgvsp)])))
}

plot(rollmean(aggregate, k = 40), xlab = "Position", type = "l", main = "Proportion of Aggregating (Red) and Misfolding (Blue) Variants", col = "red", ylim = c(0, 0.2), ylab = "Proportion of Variants of Each Type")
lines(rollmean(misfold, k = 40), type = "l", main = "Misfolding", col = "blue", ylim = c(0, 0.2))

#Analyze differences between hydrophobic and hydrophilic changes in aggregation region (235-315)

#235-315 (23[5-9]|2[4-9][0-9]|30[0-9]|31[0-5])
#355-418 (35[5-9]|3[6-9][0-9]|40[0-9]|41[0-8])

#misfolding: 315-355

agg_region <- maveVisOutputNoStop[grepl("(31[5-9]|3[2-4][0-9]|35[0-5])", maveVisOutputNoStop$hgvsp),]
#agg_region <- agg_region[agg_region$score > 1.2,]
agg_region <- agg_region[!grepl("=", agg_region$hgvsp),]
agg_region <- agg_region[!grepl("Ter", agg_region$hgvsp),]


agg_region$variant <- parseHGVS(agg_region$hgvsp)$variant
agg_region$ancestral <- parseHGVS(agg_region$hgvsp)$ancestral


charged <- c()

for (i in 1:length(agg_region$hgvsp)) {
  if (agg_region$variant[i] == "Lys" | agg_region$variant[i] == "Arg") {
    charged[i] <- agg_region$score[i]
  }
}

t.test(charged, agg_region$score)
boxplot(charged)
boxplot(agg_region$score)

aa_occurrence <- c()
#colnames(aa_occurrence) = c("variant", "n")
#list of amino acids

aa_ancestral <- c()
aa_variant <- c()

aa_code <- c("Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His", "Ile", "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val")

for(i in 1:20) {
  aa_ancestral[i] <- (length(grep(aa_code[i], agg_region$ancestral)) / length(agg_region$ancestral))
}

for(i in 1:20) {
  aa_variant[i] <- (length(grep(aa_code[i], agg_region$variant)) / length(agg_region$variant))
}


aa_diff = aa_variant - aa_ancestral
aa_diff <- data.frame(aa_code, aa_diff)
aa_diff <- aa_diff[order(-aa_diff$aa_diff),]

barplot(aa_diff$aa_diff, names.arg = aa_diff$aa_code, xlab = "Variant", ylab = "Frequency", main = "Difference in Occurence Between Variant and Ancestral (Position 315-355) - Misfolding")

#Varity analysisb
maveVisOutputNoStop$aa_pos <- parseHGVS(maveVisOutputNoStop$hgvsp, aacode = 1)$start
maveVisOutputNoStop$aa_ref <- (parseHGVS(maveVisOutputNoStop$hgvsp, aacode = 1))$ancestral
maveVisOutputNoStop$aa_alt <- parseHGVS(maveVisOutputNoStop$hgvsp, aacode = 1)$variant

aa_analysis <- full_join(maveVisOutputNoStop, varity_output, by = c("aa_ref", "aa_alt"))
aa_analysis <- distinct(aa_analysis, hgvsp, .keep_all = T)
aa_analysis <- aa_analysis[!is.na(aa_analysis$VARITY_ER),]

#235-315
#315-355
#355-418


aa_analysis_p1 <- aa_analysis[grepl("(31[5-9]|3[2-4][0-9]|35[0-5])", aa_analysis$pos),]
aa_analysis_p2 <- aa_analysis[grepl("(35[5-9]|3[6-9][0-9]|40[0-9]|41[0-8])", aa_analysis$pos),]

#significant
val1 <- (aa_analysis$pi_delta[aa_analysis$aa_pos.x > 315 & aa_analysis$aa_pos.x < 355 & aa_analysis$score > 1.4])
val2 <- ((aa_analysis$pi_delta[aa_analysis$aa_pos.x > 315 & aa_analysis$aa_pos.x < 355]))
t.test(val1, val2)
boxplot(val1, val2, main = "Change in pI of side chain (Pos. 315-355) of aggregating variants", names = c("Aggregating", "All"))

#trending towards significance
val1 <- (aa_analysis$mw_delta[aa_analysis$aa_pos.x > 235 & aa_analysis$aa_pos.x < 315 & aa_analysis$score > 1.4])
val2 <- ((aa_analysis$mw_delta[aa_analysis$aa_pos.x > 235 & aa_analysis$aa_pos.x < 315]))
wilcox.test(val1, val2)

#significant
val1 <- (aa_analysis$mw_delta[aa_analysis$aa_pos.x > 355 & aa_analysis$aa_pos.x < 418 & aa_analysis$score < 0.4])
val2 <- ((aa_analysis$mw_delta[aa_analysis$aa_pos.x > 355 & aa_analysis$aa_pos.x < 418]))
t.test(val1, val2)
boxplot(val1, val2, main = "Change in M.W. of side chain (Pos. 355-418) of misfolding variants", names = c("Misfolding", "All"))


val1 <- (aa_analysis$mw_delta[aa_analysis$aa_pos.x > 355 & aa_analysis$aa_pos.x < 418 & aa_analysis$score < 0.4])
val2 <- ((aa_analysis$mw_delta[aa_analysis$aa_pos.x > 355 & aa_analysis$aa_pos.x < 418]))
t.test(val1, val2)

#significant
val1 <- (aa_analysis$asa_delta[aa_analysis$aa_pos.x > 355 & aa_analysis$aa_pos.x < 418 & aa_analysis$score < 0.4])
val2 <- ((aa_analysis$asa_delta[aa_analysis$aa_pos.x > 355 & aa_analysis$aa_pos.x < 418]))
t.test(val1, val2)
boxplot(val1, val2, main = "Change in ASA of side chain (Pos. 355-418) of misfolding variants", names = c("Misfolding", "All"))


#trending towards significance (Wilcox)
val1 <- (aa_analysis$hydrophobic_delta[aa_analysis$aa_pos.x > 355 & aa_analysis$aa_pos.x < 418 & aa_analysis$score > 1.4])
val2 <- ((aa_analysis$hydrophobic_delta[aa_analysis$aa_pos.x > 355 & aa_analysis$aa_pos.x < 418]))
wilcox.test(val1, val2)

#significant
val1 <- (aa_analysis$vadw_delta[aa_analysis$aa_pos.x > 355 & aa_analysis$aa_pos.x < 418 & aa_analysis$score < 0.5])
val2 <- ((aa_analysis$vadw_delta[aa_analysis$aa_pos.x > 355 & aa_analysis$aa_pos.x < 418]))
t.test(val1, val2)


#almost near significant
val1 <- (aa_analysis$charge_delta[aa_analysis$aa_pos.x > 355 & aa_analysis$aa_pos.x < 418 & aa_analysis$score > 1.4])
val2 <- ((aa_analysis$charge_delta[aa_analysis$aa_pos.x > 355 & aa_analysis$aa_pos.x < 418]))

#t.test(val1, val2)
length(val1[val1 == 1])
length(val1[val1 == 0])
length(val1[val1 == -1])

length(val2[val2 == 1]) 
length(val2[val2 == 0]) 
length(val2[val2 == -1]) 
chisq.test(c(length(val1[val1 == 1]), length(val1[val1 == 0]), length(val1[val1 == -1])), p = c(length(val2[val2 == 1]) /329, length(val2[val2 == 0]) /329, length(val2[val2 == -1]) /329))

val1 <- (aa_analysis$charge_delta[aa_analysis$aa_pos.x > 235 & aa_analysis$aa_pos.x < 315 & aa_analysis$score < 0.4])
val2 <- ((aa_analysis$charge_delta[aa_analysis$aa_pos.x > 235 & aa_analysis$aa_pos.x < 315]))
wilcox.test(val1, val2)
chisq.test(c(length(val1[val1 == 1]), length(val1[val1 == 0]), length(val1[val1 == -1])), p = c(length(val2[val2 == 1]) /409, length(val2[val2 == 0]) /409, length(val2[val2 == -1]) /409))


hist(stk$score[grepl("Ter", stk$hgvs_pro)], breaks = 20, col = c2, freq = T, xlim = c(-1,4))
hist(stk$score[grepl("=", stk$hgvs_pro)], breaks = 30, col = c1, add = T, freq = T, xlim = c(-1,4))

