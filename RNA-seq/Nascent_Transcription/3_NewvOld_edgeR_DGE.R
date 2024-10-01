# New to old Differential Expression Analysis 
## TT-Seq vs Total RNA-Seq counts compared per sample

library(ggplot2)
library(GenomicRanges)
library(rtracklayer)
library(edgeR)
library(gridExtra)
library(data.table)
library(dplyr)
library(stringr)
library(matrixStats)
library(GO.db)

setwd("~/external/RCLONE/nci_scratch/ec2963/PhD_Proj1_AR_Sig/Total_RNA-Seq_Level_3/analysis/New_4sU_vs_unlabelled_RNA")

SAMPLES.NAME <- "GENE_TT_Total_RNASeq_LNCaP_DHT"
OUTPUT.FOLDER <- "DGE_"
TABLES.FOLDER <- "tables/"

raw.read.counts <- read.csv(paste0(TABLES.FOLDER, SAMPLES.NAME, "_expected_count_table.csv"), row.names=1)
tpm.data <- read.csv(paste0(TABLES.FOLDER, SAMPLES.NAME, "_TPM_table.csv"), row.names=1)

raw.read.counts <- round(raw.read.counts)
colnames(raw.read.counts) <- gsub("_TT_TT","_TT",colnames(raw.read.counts))

counts <- data.matrix(raw.read.counts) 
colnames(counts)<-colnames(raw.read.counts)
colnames(tpm.data) <- colnames(raw.read.counts)
tpm.data.matrix <- data.matrix(tpm.data)

sample.names <- colnames(raw.read.counts)

# range of library size/sequencing depth
library.size <- round(colSums(counts)/1e6, 1)
library.size

# The names argument tells the barplot to use the sample names on the x-axis
# The las argument rotates the axis names
barplot(library.size,names=colnames(library.size),las=2)
# Add a title to the plot
title("Barplot of TT- and Total RNA-Seq library sizes")


# perform PCA on log transformed count values
pca <- prcomp(t(log(counts+1)), center=TRUE)
# also see summary(pca)
pc1.var <- round(pca$sdev[1]^2/sum(pca$sdev^2)*100,2) 
pc2.var <- round(pca$sdev[2]^2/sum(pca$sdev^2)*100,2)
pc.data.frame <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], Name=colnames(counts), stringsAsFactors=F)

makeLab <- function(x,pc) paste0("PC",pc,": ", x, "% variance")

ggplot(pc.data.frame, aes(x=PC1, y=PC2, label=Name, fill=as.factor(Name)))+
  geom_point(aes(colour=factor(gsub("_Rep1|_Rep2", "", gsub("_Total|_TT", "", Name)))), show.legend = F) + 
  geom_text(aes(colour = factor(gsub("_Rep1|_Rep2", "", gsub("_Total|_TT", "", Name)))), show.legend = F) +
  xlab(makeLab(pc1.var,1)) + ylab(makeLab(pc2.var,2)) +
  expand_limits(x=c(-200, 220)) +
  expand_limits(y=c(-60,80)) +
  ggtitle("PC1 vs PC2 for log(count) of TT- and Total RNA-Seq in LNCaP DHT Time Course")

# hierarchical clustering for log(count) values
hc <- hclust(dist(t(log(counts+1))))       
plot(hc, labels=colnames(counts), main="Clustering samples on log(counts+1) of \nTT- and Total RNA-Seq in LNCaP DHT Time Course.")


# perform PCA on TPM values (use log(TPM) as TPM can be quite swayed by 
# differences in sample library size)
pca <- prcomp(t(log(tpm.data.matrix+1)), center=TRUE)
# also see summary(pca)
pc1.var <- round(pca$sdev[1]^2/sum(pca$sdev^2)*100,2) 
pc2.var <- round(pca$sdev[2]^2/sum(pca$sdev^2)*100,2)
pc.data.frame <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], Name=colnames(tpm.data.matrix))

ggplot(pc.data.frame, aes(x=PC1, y=PC2, label=Name, fill=as.factor(Name)))+
  geom_point(aes(colour=factor(gsub("_Rep1|_Rep2", "", gsub("_Total|_TT", "", Name)))), show.legend = F) + 
  geom_text(aes(colour = factor(gsub("_Rep1|_Rep2", "", gsub("_Total|_TT", "", Name)))), show.legend = F) +
  ggtitle("PC1 vs PC2 for log(TPM) values of TT- and Total RNA-Seq in LNCaP DHT Time Course") +
  xlab(makeLab(pc1.var,1)) + ylab(makeLab(pc2.var,2)) +
  expand_limits(x=c(-150,150)) +
  expand_limits(y=c(-50,50)) +
  theme(legend.position = "none")

# hierarchical clustering for log(TPM) vlaues
hc <- hclust(dist(t(log(tpm.data.matrix+1))))       
plot(hc, labels=colnames(tpm.data.matrix), main="Clustering samples on log(TPM+1) of \nTT- and Total RNA-Seq in LNCaP DHT Time Course.")


## Differential Gene Expression Analysis
### Design Matrix and Count Data
### Subset per treatment
Counts_EtOH <- counts[,c(1:2,11:12)]
Counts_5min <- counts[,c(3:4,13:14)]
Counts_10min <- counts[,c(5:6,15:16)]
Counts_20min <- counts[,c(7:8,17:18)]
Counts_30min <- counts[,c(9:10,19:20)]

tpm.data.matrix_EtOH <- tpm.data.matrix[,c(1:2,11:12)]
tpm.data.matrix_5min <- tpm.data.matrix[,c(3:4,13:14)]
tpm.data.matrix_10min <- tpm.data.matrix[,c(5:6,15:16)]
tpm.data.matrix_20min <- tpm.data.matrix[,c(7:8,17:18)]
tpm.data.matrix_30min <- tpm.data.matrix[,c(9:10,19:20)]


# set up design matrix
## EtOH
group.types_EtOH <- gsub("_Rep1|_Rep2", "", colnames(Counts_EtOH))
group_EtOH <- factor(group.types_EtOH)
design_EtOH <- model.matrix(~0+group_EtOH) 
colnames(design_EtOH) <- gsub("group", "", colnames(design_EtOH))
colnames(design_EtOH) <- c("LNCaP_EtOH_30min_TT", "LNCaP_EtOH_30min_Total")
design_EtOH
dge_gene_EtOH <- DGEList(counts=Counts_EtOH, group=group_EtOH)
dge_gene_EtOH$samples
keep <- rowSums(cpm(dge_gene_EtOH) > 0) >= 1
table(keep)
# FALSE  TRUE 
# 26623 34033 
dge_gene_EtOH <- dge_gene_EtOH[keep, , keep.lib.sizes=FALSE]
# TMM normalization
dge_gene_EtOH_norm1 <- calcNormFactors(dge_gene_EtOH, method="TMM")
dge_gene_EtOH_norm1$samples
### GLM estimates of variance (dispersion)
dge_gene_EtOH_norm1 <- estimateDisp(dge_gene_EtOH_norm1, design_EtOH)
dge_gene_EtOH_norm1$common.dispersion
#[1] 0.007651
dge_gene_EtOH_norm1 <- estimateGLMCommonDisp(dge_gene_EtOH_norm1, design_EtOH, verbose=TRUE)
#Disp = 0.0075 , BCV = 0.0866
dge_gene_EtOH_norm1 <- estimateGLMTrendedDisp(dge_gene_EtOH_norm1, design_EtOH)
dge_gene_EtOH_norm1 <- estimateGLMTagwiseDisp(dge_gene_EtOH_norm1, design_EtOH)

## 5min
group.types_5min <- gsub("_Rep1|_Rep2", "", colnames(Counts_5min))
group_5min <- factor(group.types_5min)
design_5min <- model.matrix(~0+group_5min) 
colnames(design_5min) <- gsub("group_5min", "", colnames(design_5min))
design_5min
dge_gene_5min <- DGEList(counts=Counts_5min, group=group_5min)
dge_gene_5min$samples
keep <- rowSums(cpm(dge_gene_5min) > 0) >= 1
table(keep)
# FALSE  TRUE 
# 26843 33813  
dge_gene_5min <- dge_gene_5min[keep, , keep.lib.sizes=FALSE]
# TMM normalization
dge_gene_5min_norm1 <- calcNormFactors(dge_gene_5min, method="TMM")
dge_gene_5min_norm1$samples
### GLM estimates of variance (dispersion)
dge_gene_5min_norm1 <- estimateDisp(dge_gene_5min_norm1, design_5min)
dge_gene_5min_norm1$common.dispersion
#[1] 0.01237
dge_gene_5min_norm1 <- estimateGLMCommonDisp(dge_gene_5min_norm1, design_5min, verbose=TRUE)
#Disp = 0.01252 , BCV = 0.1119 
dge_gene_5min_norm1 <- estimateGLMTrendedDisp(dge_gene_5min_norm1, design_5min)
dge_gene_5min_norm1 <- estimateGLMTagwiseDisp(dge_gene_5min_norm1, design_5min)

## 10min
group.types_10min <- gsub("_Rep1|_Rep2", "", colnames(Counts_10min))
group_10min <- factor(group.types_10min)
design_10min <- model.matrix(~0+group_10min) 
colnames(design_10min) <- gsub("group_10min", "", colnames(design_10min))
design_10min
dge_gene_10min <- DGEList(counts=Counts_10min, group=group_10min)
dge_gene_10min$samples
keep <- rowSums(cpm(dge_gene_10min) > 0) >= 1
table(keep)
# FALSE  TRUE 
# 26926 33730 
dge_gene_10min <- dge_gene_10min[keep, , keep.lib.sizes=FALSE]
# TMM normalization
dge_gene_10min_norm1 <- calcNormFactors(dge_gene_10min, method="TMM")
dge_gene_10min_norm1$samples
### GLM estimates of variance (dispersion)
dge_gene_10min_norm1 <- estimateDisp(dge_gene_10min_norm1, design_10min)
dge_gene_10min_norm1$common.dispersion
#[1] 0.01244
dge_gene_10min_norm1 <- estimateGLMCommonDisp(dge_gene_10min_norm1, design_10min, verbose=TRUE)
#Disp = 0.01225 , BCV = 0.1107 
dge_gene_10min_norm1 <- estimateGLMTrendedDisp(dge_gene_10min_norm1, design_10min)
dge_gene_10min_norm1 <- estimateGLMTagwiseDisp(dge_gene_10min_norm1, design_10min)

## 20min
group.types_20min <- gsub("_Rep1|_Rep2", "", colnames(Counts_20min))
group_20min <- factor(group.types_20min)
design_20min <- model.matrix(~0+group_20min) 
colnames(design_20min) <- gsub("group_20min", "", colnames(design_20min))
design_20min
dge_gene_20min <- DGEList(counts=Counts_20min, group=group_20min)
dge_gene_20min$samples
keep <- rowSums(cpm(dge_gene_20min) > 0) >= 1
table(keep)
# FALSE  TRUE 
# 27263 33393  
dge_gene_20min <- dge_gene_20min[keep, , keep.lib.sizes=FALSE]
# TMM normalization
dge_gene_20min_norm1 <- calcNormFactors(dge_gene_20min, method="TMM")
dge_gene_20min_norm1$samples
### GLM estimates of variance (dispersion)
dge_gene_20min_norm1 <- estimateDisp(dge_gene_20min_norm1, design_20min)
dge_gene_20min_norm1$common.dispersion
#[1] 0.08669
dge_gene_20min_norm1 <- estimateGLMCommonDisp(dge_gene_20min_norm1, design_20min, verbose=TRUE)
#Disp = 0.0847 , BCV = 0.291
dge_gene_20min_norm1 <- estimateGLMTrendedDisp(dge_gene_20min_norm1, design_20min)
dge_gene_20min_norm1 <- estimateGLMTagwiseDisp(dge_gene_20min_norm1, design_20min)

## 30min
group.types_30min <- gsub("_Rep1|_Rep2", "", colnames(Counts_30min))
group_30min <- factor(group.types_30min)
design_30min <- model.matrix(~0+group_30min) 
colnames(design_30min) <- gsub("group_30min", "", colnames(design_30min))
design_30min
dge_gene_30min <- DGEList(counts=Counts_30min, group=group_30min)
dge_gene_30min$samples
keep <- rowSums(cpm(dge_gene_30min) > 0) >= 1
table(keep)
# FALSE  TRUE 
# 26671 33985  
dge_gene_30min <- dge_gene_30min[keep, , keep.lib.sizes=FALSE]
# TMM normalization
dge_gene_30min_norm1 <- calcNormFactors(dge_gene_30min, method="TMM")
dge_gene_30min_norm1$samples
### GLM estimates of variance (dispersion)
dge_gene_30min_norm1 <- estimateDisp(dge_gene_30min_norm1, design_30min)
dge_gene_30min_norm1$common.dispersion
#[1] 0.01257
dge_gene_30min_norm1 <- estimateGLMCommonDisp(dge_gene_30min_norm1, design_30min, verbose=TRUE)
#Disp = 0.01313 , BCV = 0.1146 
dge_gene_30min_norm1 <- estimateGLMTrendedDisp(dge_gene_30min_norm1, design_30min)
dge_gene_30min_norm1 <- estimateGLMTagwiseDisp(dge_gene_30min_norm1, design_30min)

### Genewise biological coefficient of variation (BCV)
plotBCV(dge_gene_EtOH_norm1)
plotBCV(dge_gene_5min_norm1)
plotBCV(dge_gene_10min_norm1)
plotBCV(dge_gene_20min_norm1)
plotBCV(dge_gene_30min_norm1)

### Fit GLM QL and make contrasts
##GLM
library(statmod)

fit_glm_EtOH <- glmQLFit(dge_gene_EtOH_norm1, design_EtOH)
glm.EtOH.TT.vs.Total <- glmQLFTest(fit_glm_EtOH, contrast=c(1, -1))

fit_glm_5min <- glmQLFit(dge_gene_5min_norm1, design_5min)
glm.5min.TT.vs.Total <- glmQLFTest(fit_glm_5min, contrast=c(1, -1))

fit_glm_10min <- glmQLFit(dge_gene_10min_norm1, design_10min)
glm.10min.TT.vs.Total <- glmQLFTest(fit_glm_10min, contrast=c(1, -1))

fit_glm_20min <- glmQLFit(dge_gene_20min_norm1, design_20min)
glm.20min.TT.vs.Total <- glmQLFTest(fit_glm_20min, contrast=c(1, -1))

fit_glm_30min <- glmQLFit(dge_gene_30min_norm1, design_30min)
glm.30min.TT.vs.Total <- glmQLFTest(fit_glm_30min, contrast=c(1, -1))

# add all contrasts to a list for subsequent processing
glm.list <- list(glm.EtOH.TT.vs.Total=glm.EtOH.TT.vs.Total, 
                 glm.5min.TT.vs.Total=glm.5min.TT.vs.Total, 
                 glm.10min.TT.vs.Total=glm.10min.TT.vs.Total,
                 glm.20min.TT.vs.Total=glm.20min.TT.vs.Total,
                 glm.30min.TT.vs.Total=glm.30min.TT.vs.Total)

save.image("/scratch/tr07/ec2963/PhD_Proj1_AR_Sig/Total_RNA-Seq_Level_3/analysis/New_4sU_vs_unlabelled_RNA/GENE_TT_Total_RNASeq_EdgeR_TU_Normalised.RData")

# annotation
# FROM GTF FILE
gtf.annotation.file <- "/g/data/tr07/zoo/NCI_Epigenetics/NCI_Active_Jo_Projects_tr07/Users/elypar/Assets/DGE/data/hg38/gtf/gene_annotation.tsv"
gtf.anno <- read.table(gtf.annotation.file, sep= "\t", header=TRUE, stringsAsFactors=F)
colnames(gtf.anno)
# annotation file downloaded from http://www.ensembl.org/biomart/martview/ (July2017)
# has following columns: ensembl.gene.ID, chr, start, end, strand, 
# description, HGNC.symbol, entrez.gene.ID
annotation.file <- "/g/data/tr07/zoo/NCI_Epigenetics/NCI_Active_Jo_Projects_tr07/Users/elypar/Projects/CTCF_Knockdown/RNA_Seq_Level_3/data/Ensembl_87_GRCh38.p7_biomart_export.tsv"

if (!(file.exists(annotation.file))){
  untar(paste0("input/", annotation.file, ".tgz"))
}

DT <- fread(annotation.file)
setnames(DT, gsub(" ", ".", colnames(DT)))
setkey(DT, ensembl.gene.ID)
head(DT)

lengths.gene <- read.csv(paste0("tables/GENE_TT_Total_RNASeq_LNCaP_DHT_effective_length_table.csv"), row.names=1)
# get the average gene length for each row
lengths <- apply(lengths.gene, 1, mean)

FDR.CUTOFF <- 0.05
LOG.FC.CUTOFF <- 1

#### GLM Results

res <- lapply(names(glm.list), function(contrast){
  
  print (contrast)
  comparison <- gsub("glm.", "", contrast)
  
  partA <- strsplit(comparison, ".vs.")[[1]][1]
  partB <- strsplit(comparison, ".vs.")[[1]][2]
  
  SUB.FOLDER <- paste0(OUTPUT.FOLDER, gsub("glm.", "", contrast), "/")
  
  if (!(file.exists(SUB.FOLDER))){
    dir.create(SUB.FOLDER)
  }
  # Top table
  # the default method used to adjust p-values for multiple testing is BH.
  tt <- topTags(glm.list[[contrast]], n=nrow(counts))$table
  
  m <- match(rownames(tt), gtf.anno$gene.id)
  
  # assign the rownames(tt) as the gene_id as more specific with version 
  # number at end as originates from original gtf
  tt$gene.id <- rownames(tt)
  tt$gene.symbol <- gtf.anno$gene.name[m]
  tt$chr <- gtf.anno$chr[m]
  tt$start <- gtf.anno$start[m]
  tt$end <- gtf.anno$end[m]
  tt$strand <- gtf.anno$strand[m]
  tt$gene.type <- gtf.anno$gene.type[m]
  
  m <- match(gsub("\\.[0-9]*", "", rownames(tt)), DT$Gene.stable.ID)
  
  tt$description <- DT$Gene.description[m]
  tt$entrez.gene.id <- DT$Gene.name[m]
  
  # only keep chromosome names beginning with chr1..22, X, Y; 
  # remove patch chromosome assignments 
  # like JH806587.1, JH806587.1 etc
  tt <- tt[grep("chr*", tt$chr),]
  
  #Volcano plot
  #  plot(tt$logFC, -log10(tt$PValue), type="n", xlab="BRG1 KD <- -> scrambled logFC", ylab="-log10(p.value)", main="Volcano plot of Scrambled vs BRG1 KD")
  #plot(tt$logFC, -log10(tt$PValue), type="n", xlab=paste0(partB, " <- -> ", partA, " logFC"), ylab="-log10(p.value)", main=paste0("Volcano plot of ", partA, " vs ", partB))
  #text(tt$logFC, -log10(tt$PValue), labels = tt$gene.symbol, cex=0.5)
  #abline(h=-log10(tt$PValue[sum(tt$FDR < 0.05)]), col="red")
  
  vp.data <- tt[c("gene.symbol", "logFC", "PValue", "FDR")]
  vp.data = mutate(vp.data, sig=ifelse(((vp.data$FDR<FDR.CUTOFF)&(abs(vp.data$logFC)>LOG.FC.CUTOFF)),
                                       paste0("FDR<", FDR.CUTOFF), "Not Sig"))
  
  p <- ggplot(vp.data, aes(logFC, -log10(PValue))) +
    geom_point(aes(col=sig)) +
    scale_color_manual(values=c("red", "black")) +
    geom_text(data=filter(vp.data, ((FDR<FDR.CUTOFF)&(abs(logFC)>LOG.FC.CUTOFF))), aes(label=gene.symbol)) +
    ggtitle(paste0(comparison, " volcano plot")) +
    # edit the x label to signify contrast direction
    xlab(paste0(partB, " <- ", "logFC", " -> ", partA))
  
  print (p)
  
  # The function plotSmear generates a plot of the tagwise log-fold-changes against 
  # log-cpm (analogous to an MA-plot for microarray data). 
  # DE tags are highlighted on the plot
  de2 <- decideTestsDGE(glm.list[[contrast]], p.value=FDR.CUTOFF, lfc=LOG.FC.CUTOFF)
  de2tags <- rownames(glm.list[[contrast]])[as.logical(de2)]
  plotSmear(glm.list[[contrast]], de.tags=de2tags, 
            main=paste0("smear plot with GLM p.value < ", FDR.CUTOFF, " and LFC=", 
                        LOG.FC.CUTOFF, " cutoffs"))
  abline(h = c(-2, 2), col = "blue")
  
  # defining significant as FDR < FDR.CUTOFF and abs(logFC) > LOG.FC.CUTOFF 
  # add DGE.status column
  # UP, DOWN, NC
  tt$DGE.status <- "NC"
  # defining significant as FDR < FDR.CUTOFF and abs(logFC) > LOG.FC.CUTOFF 
  # for UP/DOWN
  if (nrow(tt[((tt$FDR < FDR.CUTOFF)&(tt$logFC > LOG.FC.CUTOFF)),]) > 0){
    tt[((tt$FDR < FDR.CUTOFF)&(tt$logFC > LOG.FC.CUTOFF)),]$DGE.status <- "UP"
  }
  if (nrow(tt[((tt$FDR < FDR.CUTOFF)&(tt$logFC < -LOG.FC.CUTOFF)),]) > 0){
    tt[((tt$FDR < FDR.CUTOFF)&(tt$logFC < -LOG.FC.CUTOFF)),]$DGE.status <- "DOWN"
  }
  
  # filtering/re-ordering
  column.order <- c(6:13, 1:5)
  tt <- tt[column.order]
  tt[is.na(tt)] <- ""
  
  write.table(tt, paste0(SUB.FOLDER, comparison, ".UpperQRT_EdgeR.GLM.DGE.FDR0.05.LFC1.tsv"), sep="\t", quote=F, row.names=F)
  
  print (nrow(tt[((tt$FDR < FDR.CUTOFF)&(abs(tt$logFC) > LOG.FC.CUTOFF)),]))
  print (paste0("UP: ", nrow(tt[(tt$DGE.status=="UP"),])))
  print (paste0("DOWN: ", nrow(tt[(tt$DGE.status=="DOWN"),])))
  
  # defining significant as FDR < FDR.CUTOFF and abs(logFC) > LOG.FC.CUTOFF
  sigtt <- tt[((tt$FDR < FDR.CUTOFF)&(abs(tt$logFC) > LOG.FC.CUTOFF)),]
  
})


