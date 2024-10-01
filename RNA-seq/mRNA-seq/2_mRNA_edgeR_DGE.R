library(edgeR)
library(ggplot2)
library(gridExtra)
library(data.table)
library(goseq)
library(GO.db)
library(dplyr)


SAMPLES.NAME <- "LNCaP_DHT"

OUTPUT.FOLDER <- "DGE_"

TABLES.FOLDER <- "tables/"

raw.read.counts <- read.csv(paste0(TABLES.FOLDER, SAMPLES.NAME, "_RUVr_norm_count_table.csv"), row.names=1)
tpm.data <- read.csv(paste0(TABLES.FOLDER, SAMPLES.NAME, "_TPM_table.csv"), row.names=1)

# round the rsem gene expected counts values to the nearest integer to input into edgeR
raw.read.counts <- round(raw.read.counts)
# matrix of counts with ENSGxxxxxxxx tags

counts <- data.matrix(raw.read.counts) 
colnames(counts)<-colnames(raw.read.counts)
counts<-counts[, c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15)]

colnames(tpm.data) <- colnames(raw.read.counts)
tpm.data<-tpm.data[, c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15)]
tpm.data.matrix <- data.matrix(tpm.data)

sample.names <- colnames(raw.read.counts)
# range of library size/sequencing depth
library.size <- round(colSums(counts)/1e6, 1)


# perform PCA on log transformed count values
pca <- prcomp(t(log(counts+1)), center=TRUE)
# also see summary(pca)
pc1.var <- round(pca$sdev[1]^2/sum(pca$sdev^2)*100,2) 
pc2.var <- round(pca$sdev[2]^2/sum(pca$sdev^2)*100,2)
pc.data.frame <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], Name=colnames(counts), stringsAsFactors=F)

makeLab <- function(x,pc) paste0("PC",pc,": ", x, "% variance")

ggplot(pc.data.frame, aes(x=PC1, y=PC2, label=Name))+
  geom_point() + geom_text(aes(colour = factor(gsub("_Rep1|_Rep2|_Rep3", "", pc.data.frame$Name)))) +
  xlab(makeLab(pc1.var,1)) + ylab(makeLab(pc2.var,2)) +
  ggtitle("PC1 vs PC2 for log(count) of LNCaP DHT") +
  expand_limits(x=c(-100, 100)) +
  theme(legend.position = "none")

# perform PCA on TPM values (use log(TPM) as TPM can be quite swayed by 
# differences in sample library size)
pca <- prcomp(t(log(tpm.data.matrix+1)), center=TRUE)
# also see summary(pca)
pc1.var <- round(pca$sdev[1]^2/sum(pca$sdev^2)*100,2) 
pc2.var <- round(pca$sdev[2]^2/sum(pca$sdev^2)*100,2)
pc.data.frame <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], Name=colnames(tpm.data.matrix))

ggplot(pc.data.frame, aes(x=PC1, y=PC2, label=Name))+
  geom_point() + geom_text(aes(colour=factor(gsub("_Rep1|_Rep2|_Rep3", "", pc.data.frame$Name)))) +
  xlab(makeLab(pc1.var,1)) + ylab(makeLab(pc2.var,2)) +
  ggtitle("PC1 vs PC2 for TPM values of LNCaP DHT") +
  expand_limits(x=c(-50,50)) +
  theme(legend.position = "none")


# hierarchical clustering for log(TPM) vlaues
hc <- hclust(dist(t(log(tpm.data.matrix+1))))       
plot(hc, labels=colnames(tpm.data.matrix), main="Clustering samples on log(TPM+1)")

hc <- hclust(dist(t(log(counts+1))))       
plot(hc, labels=colnames(counts), main="Clustering samples on log(counts+1)")

# Filter out ENSGxxxx tags whose coverage is so low that any group differences 
# aren't truly "real". 
# filter out tags whose rowcount <= degrees of freedom.
counts <- counts[rowSums(counts) >= 3,]
tpm.data.matrix <- tpm.data.matrix[rowSums(tpm.data.matrix) >= 3, ]

# set up design matrix
group.types <- gsub("_Rep1|_Rep2|_Rep3", "", colnames(counts))

group <- factor(group.types)

design <- model.matrix(~0+group) 
colnames(design) <- gsub("group", "", colnames(design))


dge_gene <- DGEList(counts=counts, group=group)
dge_gene$samples

###Filter out lowly expressed genes
cpm.cutoff <- 10/(min(dge_gene$samples$lib.size)/1000000)
cpm.cutoff
#0.3596098
##Use 0.6 for simplicity
keep <- rowSums(cpm(dge_gene) > 0.6) >= 3 # 3 is used here as each group as 3 replicates
table(keep)

# TMM normalization
dge_gene_norm1 <- calcNormFactors(dge_gene, method="TMM")
dge_gene_norm1$samples

# perform PCA on log(CPMS)
cpms <- cpm(dge_gene_norm1)
log.cpms <- log(cpms + 1)
# perform PCA on CPM values (use log(CPM) as CPM can be quite affected by differences in 
# sample library size)
pca <- prcomp(t(log.cpms), center=TRUE)
# also see summary(pca)
pc1.var <- round(pca$sdev[1]^2/sum(pca$sdev^2)*100,2) 
pc2.var <- round(pca$sdev[2]^2/sum(pca$sdev^2)*100,2)
pc.data.frame <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], Name=colnames(cpms))
ggplot(pc.data.frame, aes(x=PC1, y=PC2, label=Name))+
  geom_point() + geom_text(aes(colour = factor(gsub("_Rep1|_Rep2|_Rep3", "", pc.data.frame$Name)))) +
  xlab(makeLab(pc1.var,1)) + ylab(makeLab(pc2.var,2)) +
  ggtitle("PC1 vs PC2 for log(CPMS) of LNCaP DHT") +
  expand_limits(x=c(-100, 100)) +
  theme(legend.position = "none")


# GLM estimates of variance (dispersion)

# Dispersion
dge_gene_norm1 <- estimateDisp(dge_gene_norm1, design)
dge_gene_norm1$common.dispersion


# Fitting a model: common dispersion, fit a trended model, it the tagwise dispersion
dge_gene_norm1 <- estimateGLMCommonDisp(dge_gene_norm1, design, verbose=TRUE)
dge_gene_norm1 <- estimateGLMTrendedDisp(dge_gene_norm1, design)
dge_gene_norm1 <- estimateGLMTagwiseDisp(dge_gene_norm1, design)


# plot the genewise biological coefficient of variation (BCV) against gene abundance 
# (in log2 counts per million). 
# It displays the common, trended and tagwise BCV estimates.
plotBCV(dge_gene_norm1)

# Fitting linear model
fit <- glmFit(dge_gene_norm1, design)

##GLM
fit_glm <- glmQLFit(dge_gene_norm1, design)

glm.30mins.vs.EtOH <- glmQLFTest(fit_glm, contrast=c(0, 0, 1, 0, -1))
glm.2hrs.vs.EtOH <- glmQLFTest(fit_glm, contrast=c(0, 1, 0, 0, -1))
glm.4hrs.vs.EtOH <- glmQLFTest(fit_glm, contrast=c(0, 0, 0, 1, -1))
glm.16hrs.vs.EtOH <- glmQLFTest(fit_glm, contrast=c(1, 0, 0, 0, -1))



# add all contrasts to a list for subsequent processing
glm.list <- list(glm.30mins.vs.EtOH=glm.30mins.vs.EtOH, 
                 glm.2hrs.vs.EtOH=glm.2hrs.vs.EtOH,
                 glm.4hrs.vs.EtOH=glm.4hrs.vs.EtOH,
                 glm.16hrs.vs.EtOH=glm.16hrs.vs.EtOH)



# annotation
# FROM GTF FILE
setwd("/Volumes/griw/Cancer-Epigenetics/Projects/LNCaP_VCaP_DHT_Hi-C/Elyssa_2021_analysis_and_integration/RNA-Seq")
gtf.annotation.file <- "data/gene_annotation.tsv"
gtf.anno <- read.table(gtf.annotation.file, sep= "\t", header=TRUE, stringsAsFactors=F)

colnames(gtf.anno)

# annotation file downloaded from http://www.ensembl.org/biomart/martview/ (July2017)
# has following columns: ensembl.gene.ID, chr, start, end, strand, 
# description, HGNC.symbol, entrez.gene.ID
annotation.file <- "data/gene_annotation_mart_export.tsv"
if (!(file.exists(annotation.file))){
  untar(paste0("input/", annotation.file, ".tgz"))
}

DT <- fread(annotation.file)
setnames(DT, gsub(" ", ".", colnames(DT)))
setkey(DT, HGNC.symbol)
head(DT)

### annotations/DE results ###
library(ggplot2)
library(GenomicRanges)
library(rtracklayer)
library(edgeR)
library(gridExtra)
library(data.table)
library(dplyr)
library(stringr)
library(matrixStats)
options(stringsAsFactors=FALSE)
# Gene annotation - Load gtf
file <- "/Volumes/griw/Cancer-Epigenetics/Projects/LNCaP_VCaP_DHT_Hi-C/Elyssa_2021_analysis_and_integration/RNA-Seq/data/hg38/gtf/gencode.v35.annotation.gtf"
gtf <- fread(file, data.table = FALSE)
ex <- GRanges(gtf[,1], IRanges(gtf$V4, gtf$V5), strand = gtf$V7)
attribs <- c(gene_name = "gene_name", gene_type = "gene_type", 
             gene_id = "gene_id", tx_id = "transcript_id", tx_name = "transcript_name")
gtf.attribs <- data.frame(lapply(attribs, function(a) gsub(".*\"", "", gsub("\";", "", str_extract(gtf$V9, paste0(a, ".+?;"))))), stringsAsFactors = FALSE)
gtf.attribs <- cbind(as.data.frame(ex), gtf.attribs)
gtf.attribs <- gtf.attribs[!duplicated(gtf.attribs$gene_id),]
nrow(gtf.attribs)
head(gtf.attribs)
# Getting DE results
LOG.FC.CUTOFF <- 1.0
FDR.CUTOFF <- 0.05

##GLM
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
  de2tags <- rownames(lrt.list[[contrast]])[as.logical(de2)]
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
  
  write.table(tt, paste0(SUB.FOLDER, comparison, ".GLM.DGE.tsv"), sep="\t", quote=F, row.names=F)
  
  print (nrow(tt[((tt$FDR < FDR.CUTOFF)&(abs(tt$logFC) > LOG.FC.CUTOFF)),]))
  print (paste0("UP: ", nrow(tt[(tt$DGE.status=="UP"),])))
  print (paste0("DOWN: ", nrow(tt[(tt$DGE.status=="DOWN"),])))
  
  # defining significant as FDR < FDR.CUTOFF and abs(logFC) > LOG.FC.CUTOFF
  sigtt <- tt[((tt$FDR < FDR.CUTOFF)&(abs(tt$logFC) > LOG.FC.CUTOFF)),]
  
})

#[1] "glm.30mins.vs.EtOH"
#[1] 122
#[1] "UP: 74"
#[1] "DOWN: 48"
#[1] "glm.2hrs.vs.EtOH"
#[1] 793
#[1] "UP: 622"
#[1] "DOWN: 171"
#[1] "glm.4hrs.vs.EtOH"
#[1] 659
#[1] "UP: 432"
#[1] "DOWN: 227"
#[1] "glm.16hrs.vs.EtOH"
#[1] 806
#[1] "UP: 527"
#[1] "DOWN: 279"

# add annotation to TPM data and output
m <- match(rownames(tpm.data), gtf.anno$gene.id)

# assign the rownames(tt) as the gene_id as more specific with version number at end
# as originates from original gtf
tpm.data$gene.id <- rownames(tpm.data)
tpm.data$gene.symbol <- gtf.anno$gene.name[m]
tpm.data$chr <- gtf.anno$chr[m]
tpm.data$start <- gtf.anno$start[m]
tpm.data$end <- gtf.anno$end[m]
tpm.data$strand <- gtf.anno$strand[m]
tpm.data$gene.type <- gtf.anno$gene.type[m]

m <- match(gsub("\\.[0-9]*", "", rownames(tpm.data)), DT$Gene.stable.ID)

tpm.data$description <- DT$description[m]
tpm.data$entrez.gene.id <- DT$entrez.gene.ID[m]

# get averages for each set of duplicates/triplicates
tpm.data$mean_EtOH <- 
  round(rowMeans(tpm.data[,grep("EtOH", colnames(tpm.data))]), 2)
tpm.data$mean_30min <- 
  round(rowMeans(tpm.data[,grep("DHT_30min", colnames(tpm.data))]), 2)
tpm.data$mean_2h <- 
  round(rowMeans(tpm.data[,grep("DHT_2h", colnames(tpm.data))]), 2)
tpm.data$mean_4h <- 
  round(rowMeans(tpm.data[,grep("DHT_4h", colnames(tpm.data))]), 2)
tpm.data$mean_16h <- 
  round(rowMeans(tpm.data[,grep("DHT_16h", colnames(tpm.data))]), 2)


column.order <- c(17:27, 1:16)
tpm.data <- tpm.data[column.order]
tpm.data[is.na(tpm.data)] <- ""

write.table(tpm.data, paste0("annotated_LNCaP_DHT_TPM_table_GLM.tsv"), 
            sep="\t", quote=F, row.names=F)

#### NEXT

###Create top tables per comparison

# Top table  30mins vs EtOH
# the default method used to adjust p-values for multiple testing is BH.
tt1 <- topTags(lrt.30mins.vs.EtOH, n=nrow(counts))$table

m <- match(rownames(tt1), gtf.anno$gene.id)

# assign the rownames(tt) as the gene_id as more specific with version number at end
# as originates from original gtf
tt1$gene.id <- rownames(tt1)
tt1$gene.symbol <- gtf.anno$gene.name[m]
tt1$chr <- gtf.anno$chr[m]
tt1$start <- gtf.anno$start[m]
tt1$end <- gtf.anno$end[m]
tt1$strand <- gtf.anno$strand[m]
tt1$gene.type <- gtf.anno$gene.type[m]

m <- match(gsub("\\.[0-9]*", "", rownames(tt1)), DT$Ensembl.Gene.ID)

tt1$description <- DT$description[m]
tt1$entrez.gene.id <- DT$entrez.gene.ID[m]

# only keep chromosome names beginning with chr1..22, X, Y; remove patch chromosome assignments 
# like JH806587.1, JH806587.1 etc
tt1 <- tt1[grep("chr*", tt1$chr),]

# GOseq 30mins vs EtOH
lengths.gene <- read.csv(paste0("/Volumes/griw/Cancer-Epigenetics/Projects/LNCaP_VCaP_DHT_Hi-C/Elyssa_2021_analysis_and_integration/RNA-Seq/tables/LNCaP_DHT_effective_length_table.csv"), row.names=1)
# get the average gene length for each row
lengths <- apply(lengths.gene, 1, mean)
bias.data <- lengths[rownames(tt1)]
names(bias.data) <- tt1$gene.symbol
bias.data <- bias.data[!duplicated(names(bias.data))]
if (length(names(bias.data[(names(bias.data) == "")])) > 0){
  bias.data <- bias.data[-which(names(bias.data)=="")]
}
bias.data <- bias.data[-which(bias.data==0)]
if (length(names(bias.data[(is.na(names(bias.data)))])) > 0){
  bias.data <- bias.data[-which(is.na(names(bias.data)))]
}
sigtt1 <- tt1[((tt1$FDR < FDR.CUTOFF)&(abs(tt1$logFC) > LOG.FC.CUTOFF)),]
comparison.UP <- sigtt1$gene.symbol[sigtt1$logFC > 0]
comparison.DOWN <- sigtt1$gene.symbol[sigtt1$logFC < 0]

comparison.UP.DE <- sapply(names(bias.data), function (x) as.numeric(x %in% comparison.UP))
comparison.DOWN.DE <- sapply(names(bias.data), function (x) as.numeric(x %in% comparison.DOWN))

write.table(tt1, "/Volumes/griw/Cancer-Epigenetics/Projects/LNCaP_VCaP_DHT_Hi-C/Elyssa_2021_analysis_and_integration/RNA-Seq/LNCaP_DHT_30mins.vs.EtOH_DGE.tsv", sep = "\t", quote = F)


library(goseq)
library(magrittr)
library(dplyr)
library(ggplot2)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)


##Heatmap - DEG in 30mins vs EtOH ordered by PValue (top 300)
logCPM <- cpm(lrt.30mins.vs.EtOH,prior.count = 2,log = TRUE) 
o <- order(lrt.30mins.vs.EtOH$table$PValue)
logCPM <- logCPM[o[1:300],]
logCPM <- t(scale(t(logCPM)))
library(gplots)
col.pan <- colorpanel(100,"blue","white","red")
samples <- c("EtOH_R1", "EtOH_R2", "EtOH_R3", "30mins_R1", "30mins_R2", "30mins_R3", "2hrs_R1", "2hrs_R2", "2hrs_R3", "4hrs_R1", "4hrs_R2", "4hrs_R3", "16hrs_R1", "16hrs_R2", "16hrs_R3")
heatmap.2(logCPM,col=col.pan,labCol = samples, Rowv = TRUE,scale = "none",trace = "none",dendrogram = "both",cexRow = 1, cexCol = 1.4,margins = c(10,9), lhei = c(2,10),lwid = c(2,6))


###Heatmap - DEG in 2hrs vs EtOH ordered by PValue (top 300)
logCPM <- cpm(lrt.2hrs.vs.EtOH,prior.count = 2,log = TRUE) 
o <- order(lrt.2hrs.vs.EtOH$table$PValue)
logCPM <- logCPM[o[1:300],]
logCPM <- t(scale(t(logCPM)))
col.pan <- colorpanel(100,"blue","white","red")
heatmap.2(logCPM,col=col.pan, labCol = samples, Rowv = TRUE,scale = "none",trace = "none",dendrogram = "both",cexRow = 1, cexCol = 1.4,margins = c(10,9), lhei = c(2,10),lwid = c(2,6))

###Heatmap - DEG in 4hrs vs EtOH ordered by PValue (top 300)
logCPM <- cpm(lrt.4hrs.vs.EtOH,prior.count = 2,log = TRUE) 
o <- order(lrt.4hrs.vs.EtOH$table$PValue)
logCPM <- logCPM[o[1:300],]
logCPM <- t(scale(t(logCPM)))
col.pan <- colorpanel(100,"blue","white","red")
heatmap.2(logCPM,col=col.pan, labCol = samples, Rowv = TRUE,scale = "none",trace = "none",dendrogram = "both",cexRow = 1, cexCol = 1.4,margins = c(10,9), lhei = c(2,10),lwid = c(2,6))

###Heatmap - DEG in 16hrs vs EtOH ordered by PValue (top 300)
logCPM <- cpm(lrt.16hrs.vs.EtOH,prior.count = 2,log = TRUE) 
o <- order(lrt.16hrs.vs.EtOH$table$PValue)
logCPM <- logCPM[o[1:300],]
logCPM <- t(scale(t(logCPM)))
col.pan <- colorpanel(100,"blue","white","red")
heatmap.2(logCPM,col=col.pan, labCol = samples, Rowv = TRUE,scale = "none",trace = "none",dendrogram = "both",cexRow = 1, cexCol = 1.4,margins = c(10,9), lhei = c(2,10),lwid = c(2,6))

#### Top table 2hrs vs EtOH
# the default method used to adjust p-values for multiple testing is BH.
tt2 <- topTags(lrt.2hrs.vs.EtOH, n=nrow(counts))$table

m <- match(rownames(tt2), gtf.anno$gene.id)

# assign the rownames(tt) as the gene_id as more specific with version number at end
# as originates from original gtf
tt2$gene.id <- rownames(tt2)
tt2$gene.symbol <- gtf.anno$gene.name[m]
tt2$chr <- gtf.anno$chr[m]
tt2$start <- gtf.anno$start[m]
tt2$end <- gtf.anno$end[m]
tt2$strand <- gtf.anno$strand[m]
tt2$gene.type <- gtf.anno$gene.type[m]

m <- match(gsub("\\.[0-9]*", "", rownames(tt2)), DT$Ensembl.Gene.ID)

tt2$description <- DT$description[m]
tt2$entrez.gene.id <- DT$entrez.gene.ID[m]

# only keep chromosome names beginning with chr1..22, X, Y; remove patch chromosome assignments 
# like JH806587.1, JH806587.1 etc
tt2 <- tt2[grep("chr*", tt1$chr),]

# GOseq 2hrs vs EtOH
bias.data <- lengths[rownames(tt2)]
names(bias.data) <- tt2$gene.symbol
bias.data <- bias.data[!duplicated(names(bias.data))]
if (length(names(bias.data[(names(bias.data) == "")])) > 0){
  bias.data <- bias.data[-which(names(bias.data)=="")]
}
bias.data <- bias.data[-which(bias.data==0)]
if (length(names(bias.data[(is.na(names(bias.data)))])) > 0){
  bias.data <- bias.data[-which(is.na(names(bias.data)))]
}
sigtt2 <- tt2[((tt2$FDR < FDR.CUTOFF)&(abs(tt2$logFC) > LOG.FC.CUTOFF)),]
comparison.UP <- sigtt2$gene.symbol[sigtt2$logFC > 0]
comparison.DOWN <- sigtt2$gene.symbol[sigtt2$logFC < 0]

comparison.UP.DE <- sapply(names(bias.data), function (x) as.numeric(x %in% comparison.UP))
comparison.DOWN.DE <- sapply(names(bias.data), function (x) as.numeric(x %in% comparison.DOWN))

#### Top table 4hrs vs EtOH
# the default method used to adjust p-values for multiple testing is BH.
tt3 <- topTags(lrt.4hrs.vs.EtOH, n=nrow(counts))$table

m <- match(rownames(tt3), gtf.anno$gene.id)

# assign the rownames(tt) as the gene_id as more specific with version number at end
# as originates from original gtf
tt3$gene.id <- rownames(tt3)
tt3$gene.symbol <- gtf.anno$gene.name[m]
tt3$chr <- gtf.anno$chr[m]
tt3$start <- gtf.anno$start[m]
tt3$end <- gtf.anno$end[m]
tt3$strand <- gtf.anno$strand[m]
tt3$gene.type <- gtf.anno$gene.type[m]

m <- match(gsub("\\.[0-9]*", "", rownames(tt3)), DT$Ensembl.Gene.ID)

tt3$description <- DT$description[m]
tt3$entrez.gene.id <- DT$entrez.gene.ID[m]

# only keep chromosome names beginning with chr1..22, X, Y; remove patch chromosome assignments 
# like JH806587.1, JH806587.1 etc
tt3 <- tt3[grep("chr*", tt1$chr),]

# GOseq 4hrs vs EtOH
bias.data <- lengths[rownames(tt3)]
names(bias.data) <- tt3$gene.symbol
bias.data <- bias.data[!duplicated(names(bias.data))]
if (length(names(bias.data[(names(bias.data) == "")])) > 0){
  bias.data <- bias.data[-which(names(bias.data)=="")]
}
bias.data <- bias.data[-which(bias.data==0)]
if (length(names(bias.data[(is.na(names(bias.data)))])) > 0){
  bias.data <- bias.data[-which(is.na(names(bias.data)))]
}
sigtt3 <- tt3[((tt3$FDR < FDR.CUTOFF)&(abs(tt3$logFC) > LOG.FC.CUTOFF)),]
comparison.UP <- sigtt3$gene.symbol[sigtt3$logFC > 0]
comparison.DOWN <- sigtt3$gene.symbol[sigtt3$logFC < 0]

comparison.UP.DE <- sapply(names(bias.data), function (x) as.numeric(x %in% comparison.UP))
comparison.DOWN.DE <- sapply(names(bias.data), function (x) as.numeric(x %in% comparison.DOWN))


#### Top table 16hrs vs EtOH
# the default method used to adjust p-values for multiple testing is BH.
tt4 <- topTags(lrt.16hrs.vs.EtOH, n=nrow(counts))$table

m <- match(rownames(tt4), gtf.anno$gene.id)

# assign the rownames(tt) as the gene_id as more specific with version number at end
# as originates from original gtf
tt4$gene.id <- rownames(tt4)
tt4$gene.symbol <- gtf.anno$gene.name[m]
tt4$chr <- gtf.anno$chr[m]
tt4$start <- gtf.anno$start[m]
tt4$end <- gtf.anno$end[m]
tt4$strand <- gtf.anno$strand[m]
tt4$gene.type <- gtf.anno$gene.type[m]

m <- match(gsub("\\.[0-9]*", "", rownames(tt4)), DT$Ensembl.Gene.ID)

tt4$description <- DT$description[m]
tt4$entrez.gene.id <- DT$entrez.gene.ID[m]

# only keep chromosome names beginning with chr1..22, X, Y; remove patch chromosome assignments 
# like JH806587.1, JH806587.1 etc
tt4 <- tt4[grep("chr*", tt1$chr),]


## All genes from exp per cell line output
col.order <- c(3,4,5,6,1,2,7)
#column.order <- c(1:15, 1:5)
gtf.anno2 <- gtf.anno[col.order]
write.table(gtf.anno2, paste0("GTF_anno.bed"), sep="\t", quote=F, row.names=F)
