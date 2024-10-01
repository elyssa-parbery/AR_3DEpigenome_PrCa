# split ChromHMM state bed files into one for each state per timepoint 
# and save nascent and steady states for deeptools plots

# load model per time
RNA_ChromHMM_5mins_DHT_states <- read.delim("~/external/RCLONE/nci_scratch/ec2963/PhD_Proj1_AR_Sig/MOFA/Initial_TimeCourse_2020_MOFA/Integrating_EarlyDHT/Deeptools_IN/RNA_ChromHMM_5mins_DHT_states.bed", header=FALSE)
RNA_ChromHMM_10mins_DHT_states <- read.delim("~/external/RCLONE/nci_scratch/ec2963/PhD_Proj1_AR_Sig/MOFA/Initial_TimeCourse_2020_MOFA/Integrating_EarlyDHT/Deeptools_IN/RNA_ChromHMM_10mins_DHT_states.bed", header=FALSE)
RNA_ChromHMM_20mins_DHT_states <- read.delim("~/external/RCLONE/nci_scratch/ec2963/PhD_Proj1_AR_Sig/MOFA/Initial_TimeCourse_2020_MOFA/Integrating_EarlyDHT/Deeptools_IN/RNA_ChromHMM_20mins_DHT_states.bed", header=FALSE)
RNA_ChromHMM_30mins_DHT_states <- read.delim("~/external/RCLONE/nci_scratch/ec2963/PhD_Proj1_AR_Sig/MOFA/Initial_TimeCourse_2020_MOFA/Integrating_EarlyDHT/Deeptools_IN/RNA_ChromHMM_30mins_DHT_states.bed", header=FALSE)

# subset each model per state
RNA_ChromHMM_5mins_DHT_state1 <- subset(RNA_ChromHMM_5mins_DHT_states, RNA_ChromHMM_5mins_DHT_states$V4 == 1)
RNA_ChromHMM_5mins_DHT_state2 <- subset(RNA_ChromHMM_5mins_DHT_states, RNA_ChromHMM_5mins_DHT_states$V4 == 2)
RNA_ChromHMM_5mins_DHT_state3 <- subset(RNA_ChromHMM_5mins_DHT_states, RNA_ChromHMM_5mins_DHT_states$V4 == 3)
RNA_ChromHMM_5mins_DHT_state4 <- subset(RNA_ChromHMM_5mins_DHT_states, RNA_ChromHMM_5mins_DHT_states$V4 == 4)

RNA_ChromHMM_10mins_DHT_state1 <- subset(RNA_ChromHMM_10mins_DHT_states, RNA_ChromHMM_10mins_DHT_states$V4 == 1)
RNA_ChromHMM_10mins_DHT_state2 <- subset(RNA_ChromHMM_10mins_DHT_states, RNA_ChromHMM_10mins_DHT_states$V4 == 2)
RNA_ChromHMM_10mins_DHT_state3 <- subset(RNA_ChromHMM_10mins_DHT_states, RNA_ChromHMM_10mins_DHT_states$V4 == 3)
RNA_ChromHMM_10mins_DHT_state4 <- subset(RNA_ChromHMM_10mins_DHT_states, RNA_ChromHMM_10mins_DHT_states$V4 == 4)

RNA_ChromHMM_20mins_DHT_state1 <- subset(RNA_ChromHMM_20mins_DHT_states, RNA_ChromHMM_20mins_DHT_states$V4 == 1)
RNA_ChromHMM_20mins_DHT_state2 <- subset(RNA_ChromHMM_20mins_DHT_states, RNA_ChromHMM_20mins_DHT_states$V4 == 2)
RNA_ChromHMM_20mins_DHT_state3 <- subset(RNA_ChromHMM_20mins_DHT_states, RNA_ChromHMM_20mins_DHT_states$V4 == 3)
RNA_ChromHMM_20mins_DHT_state4 <- subset(RNA_ChromHMM_20mins_DHT_states, RNA_ChromHMM_20mins_DHT_states$V4 == 4)

RNA_ChromHMM_30mins_DHT_state1 <- subset(RNA_ChromHMM_30mins_DHT_states, RNA_ChromHMM_30mins_DHT_states$V4 == 1)
RNA_ChromHMM_30mins_DHT_state2 <- subset(RNA_ChromHMM_30mins_DHT_states, RNA_ChromHMM_30mins_DHT_states$V4 == 2)
RNA_ChromHMM_30mins_DHT_state3 <- subset(RNA_ChromHMM_30mins_DHT_states, RNA_ChromHMM_30mins_DHT_states$V4 == 3)
RNA_ChromHMM_30mins_DHT_state4 <- subset(RNA_ChromHMM_30mins_DHT_states, RNA_ChromHMM_30mins_DHT_states$V4 == 4)

# assign Nascent and Steady States per model based on the ChromHMM outputs
##  Time   | Nascent State | Steady State  |
##  5min   |      4        |       1       |
##  10min  |      3        |       2       |
##  20min  |      1        |       3       |
##  30min  |      2        |       3       |

Nascent_ChromHMM_State_5mins.bed <- RNA_ChromHMM_5mins_DHT_state4[,c(1:3)]
Steady_ChromHMM_State_5mins.bed <- RNA_ChromHMM_5mins_DHT_state1[,c(1:3)]

Nascent_ChromHMM_State_10mins.bed <- RNA_ChromHMM_10mins_DHT_state3[,c(1:3)]
Steady_ChromHMM_State_10mins.bed <- RNA_ChromHMM_10mins_DHT_state2[,c(1:3)]

Nascent_ChromHMM_State_20mins.bed <- RNA_ChromHMM_20mins_DHT_state1[,c(1:3)]
Steady_ChromHMM_State_20mins.bed <- RNA_ChromHMM_20mins_DHT_state3[,c(1:3)]

Nascent_ChromHMM_State_30mins.bed <- RNA_ChromHMM_30mins_DHT_state2[,c(1:3)]
Steady_ChromHMM_State_30mins.bed <- RNA_ChromHMM_30mins_DHT_state3[,c(1:3)]

# draw van diagrams of overlaps in Nascent and Steady States
library(ChIPpeakAnno)
library(GenomicRanges)

names(Nascent_ChromHMM_State_5mins.bed) <- c("chrom", "start", "end")
names(Steady_ChromHMM_State_5mins.bed) <- c("chrom", "start", "end")

names(Nascent_ChromHMM_State_10mins.bed) <- c("chrom", "start", "end")
names(Steady_ChromHMM_State_10mins.bed) <- c("chrom", "start", "end")

names(Nascent_ChromHMM_State_20mins.bed) <- c("chrom", "start", "end")
names(Steady_ChromHMM_State_20mins.bed) <- c("chrom", "start", "end")

names(Nascent_ChromHMM_State_30mins.bed) <- c("chrom", "start", "end")
names(Steady_ChromHMM_State_30mins.bed) <- c("chrom", "start", "end")

Nascent_ChromHMM_State_5mins.gr <- makeGRangesFromDataFrame(Nascent_ChromHMM_State_5mins.bed)
Steady_ChromHMM_State_5mins.gr <- makeGRangesFromDataFrame(Steady_ChromHMM_State_5mins.bed)

Nascent_ChromHMM_State_10mins.gr <- makeGRangesFromDataFrame(Nascent_ChromHMM_State_10mins.bed)
Steady_ChromHMM_State_10mins.gr <- makeGRangesFromDataFrame(Steady_ChromHMM_State_10mins.bed)

Nascent_ChromHMM_State_20mins.gr <- makeGRangesFromDataFrame(Nascent_ChromHMM_State_20mins.bed)
Steady_ChromHMM_State_20mins.gr <- makeGRangesFromDataFrame(Steady_ChromHMM_State_20mins.bed)

Nascent_ChromHMM_State_30mins.gr <- makeGRangesFromDataFrame(Nascent_ChromHMM_State_30mins.bed)
Steady_ChromHMM_State_30mins.gr <- makeGRangesFromDataFrame(Steady_ChromHMM_State_30mins.bed)


Nascent_States <- makeVennDiagram(list(Nascent_ChromHMM_State_5mins.gr, Nascent_ChromHMM_State_10mins.gr, Nascent_ChromHMM_State_20mins.gr, Nascent_ChromHMM_State_30mins.gr), NameOfPeaks=c("DHT 5min", "DHT 10min", "DHT 20min", "DHT 30min"),
                                   totalTest=250000,
                                   fill=c("#d4bcfd", "#9adfd1", "#e3f1ae", "#f3e19e"), # circle fill color
                                   col=c("#866cb2", "#6bb5a6", "#a9c166", "#febb27"), #circle border color
                                   cat.col=c("#866cb2", "#6bb5a6", "#a9c166", "#febb27"))
Steady_States <- makeVennDiagram(list(Steady_ChromHMM_State_5mins.gr, Steady_ChromHMM_State_10mins.gr, Steady_ChromHMM_State_20mins.gr, Steady_ChromHMM_State_30mins.gr), NameOfPeaks=c("DHT 5min", "DHT 10min", "DHT 20min", "DHT 30min"),
                                       totalTest=250000,
                                       fill=c("#d4bcfd", "#9adfd1", "#e3f1ae", "#f3e19e"), # circle fill color
                                       col=c("#866cb2", "#6bb5a6", "#a9c166", "#febb27"), #circle border color
                                       cat.col=c("#866cb2", "#6bb5a6", "#a9c166", "#febb27"))

# write bedfiles out for deeptools:
setwd("~/external/RCLONE/nci_scratch/ec2963/PhD_Proj1_AR_Sig/MOFA/Initial_TimeCourse_2020_MOFA/Integrating_EarlyDHT")
write.table(Nascent_ChromHMM_State_5mins.bed, file = "Deeptools_IN/Nascent_ChromHMM_State_5mins.bed", quote = F, row.names = F, col.names = F, sep = "\t")
write.table(Steady_ChromHMM_State_5mins.bed, file = "Deeptools_IN/Steady_ChromHMM_State_5mins.bed", quote = F, row.names = F, col.names = F, sep = "\t")

write.table(Nascent_ChromHMM_State_10mins.bed, file = "Deeptools_IN/Nascent_ChromHMM_State_10mins.bed", quote = F, row.names = F, col.names = F, sep = "\t")
write.table(Steady_ChromHMM_State_10mins.bed, file = "Deeptools_IN/Steady_ChromHMM_State_10mins.bed", quote = F, row.names = F, col.names = F, sep = "\t")

write.table(Nascent_ChromHMM_State_20mins.bed, file = "Deeptools_IN/Nascent_ChromHMM_State_20mins.bed", quote = F, row.names = F, col.names = F, sep = "\t")
write.table(Steady_ChromHMM_State_20mins.bed, file = "Deeptools_IN/Steady_ChromHMM_State_20mins.bed", quote = F, row.names = F, col.names = F, sep = "\t")

write.table(Nascent_ChromHMM_State_30mins.bed, file = "Deeptools_IN/Nascent_ChromHMM_State_30mins.bed", quote = F, row.names = F, col.names = F, sep = "\t")
write.table(Steady_ChromHMM_State_30mins.bed, file = "Deeptools_IN/Steady_ChromHMM_State_30mins.bed", quote = F, row.names = F, col.names = F, sep = "\t")

# save
save.image("~/external/RCLONE/nci_scratch/ec2963/PhD_Proj1_AR_Sig/ChromHMM/Early_DHT_RNA/Assigning_Nascent_Steady_States.RData")