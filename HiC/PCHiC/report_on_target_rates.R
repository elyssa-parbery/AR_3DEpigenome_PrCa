# report on target rates

setwd("/scratch/tr07/PCHiC_Level_2_LNCaP_DHT_SubSam/On_Target_Rates")

SAMPLE_NAMES <- list("DHT_10mins_Rep1","DHT_4h_1","DHT_05h_2","EtOH_30mins_Rep1","EtOH_16h_2",
                    "DHT_05h_1","DHT_10mins_Rep2","DHT_30mins_Rep2","DHT_16h_2","DHT_05mins_Rep2",
                    "DHT_20mins_Rep2","EtOH_30mins_Rep2","DHT_2h_1","EtOH_16h_1","DHT_05mins_Rep1",
                    "DHT_2h_2","DHT_30mins_Rep1","DHT_4h_2","DHT_20mins_Rep1","DHT_16h_1")
names(SAMPLE_NAMES) <- SAMPLE_NAMES

hicup.bam.length_list <- list()

on.target_list <- list()


for(i in 1:length(SAMPLE_NAMES)) {
  hicup.bam.length_list[[i]] <- read.table(paste0(SAMPLE_NAMES[[i]], ".hicup.bam.length.txt"), quote="\"", comment.char="")
  on.target_list[[i]] <- read.table(paste0(SAMPLE_NAMES[[i]], ".on-target.txt"), quote="\"", comment.char="")
}

names(hicup.bam.length_list) <- names(SAMPLE_NAMES)
names(on.target_list) <- names(SAMPLE_NAMES)

SAMPLE_NAMES.df <- as.data.frame(unlist(SAMPLE_NAMES, use.names=TRUE))
hicup.bam.length.df <- as.data.frame(unlist(hicup.bam.length_list, use.names=TRUE))
on.target_list.df <- as.data.frame(unlist(on.target_list, use.names=TRUE))

ontarget_calc.df <- merge(hicup.bam.length.df,on.target_list.df, by=0)
ontarget_calc.df$hicup.bam.div2 = ontarget_calc.df$`unlist(hicup.bam.length_list, use.names = TRUE)`/2
ontarget_calc.df$On.Target.Percentage =(ontarget_calc.df$`unlist(on.target_list, use.names = TRUE)`/ontarget_calc.df$hicup.bam.div2)*100

names(ontarget_calc.df) <- c("Sample","hicup.bam.length","on.target.length","hicup.bam.div2","On.Target.Percentage")

summary(ontarget_calc.df)
# on.target.length
#min = 12038949
#max = 26453857
#range = 14414908

### 16hour and Early timecourse only
list_16h <- c("DHT_05h_1","DHT_05h_2","DHT_16h_1","DHT_16h_2","DHT_2h_1","DHT_2h_2","DHT_4h_1","DHT_4h_2","EtOH_16h_1","EtOH_16h_2")
ontarget_calc_16h.df <- subset(ontarget_calc.df, ontarget_calc.df$Sample %in% list_16h)
ontarget_calc_Early.df <- subset(ontarget_calc.df, !(ontarget_calc.df$Sample %in% list_16h))


summary(ontarget_calc_16h.df)
# on.target.length
#min = 19624430
#max = 26453857
#range = 6829427

summary(ontarget_calc_Early.df)
# on.target.length
#min = 12038949
#max = 20575850
#range = 8536901

library(ggplot2)

ggplot(ontarget_calc.df, aes(x=as.factor(Sample), y=on.target.length)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

ggplot(ontarget_calc.df, aes(x=as.factor(Sample), y=On.Target.Percentage)) +
  geom_bar(stat = "identity") +
  ylim(c(0,100)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))


