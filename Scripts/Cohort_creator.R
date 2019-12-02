library("stringr")
cohort_t = read.table("~/Koop_Klinghammer/Misc/RTV.tsv",sep ="\t", header = T, stringsAsFactors = F)

cohort_t$Copanlisib_val = log2(as.double(str_replace(cohort_t$Copanlisib_val, pattern = ",", "."))+1)
cohort_t$Cetuximab_val = log2(as.double(str_replace(cohort_t$Cetuximab_val, pattern = ",", "."))+1)
cohort_t$Copan_Cetux_val = log2(as.double(str_replace(cohort_t$Copan_Cetux_val, pattern = ",", "."))+1)
hist(cohort_t$Copanlisib_val)
hist(cohort_t$Cetuximab_val)
hist(cohort_t$Copan_Cetux_val)

plot(quantile(cohort_t$Copanlisib_val, seq(0,1,by = .1)))
lines(x=c(1,11),y=c(min(cohort_t$Copanlisib_val), max(cohort_t$Copanlisib_val)), col = "blue")
plot(quantile(cohort_t$Cetuximab_val, seq(0,1,by = .1)))
lines(x=c(1,11),y=c(0, max(cohort_t$Cetuximab_val)), col = "green")
plot(quantile(cohort_t$Copan_Cetux_val, seq(0,1,by = .1)))
lines(x=c(1,11),y=c(0, max(cohort_t$Copan_Cetux_val)), col = "darkred")

norm_cum_sim_copanlisib = round(cumsum(cohort_t$Copanlisib_val) / sum(cohort_t$Copanlisib_val) * 100,0)



threshold_Copanlisib = median(cohort_t$Copanlisib_val)
threshold_Cetuximab_val = median(cohort_t$Cetuximab_val)
threshold_Copan_Cetux_val = median(cohort_t$Copan_Cetux_val)
