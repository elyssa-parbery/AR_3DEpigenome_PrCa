# EdgeR Size Factors calculation per sample group

## Use function to calculate per sample
Counts_EtOH <- raw.read.counts[,c(1:2,11:12)]
Counts_5min <- raw.read.counts[,c(3:4,13:14)]
Counts_10min <- raw.read.counts[,c(5:6,15:16)]
Counts_20min <- raw.read.counts[,c(7:8,17:18)]
Counts_30min <- raw.read.counts[,c(9:10,19:20)]

dge_gene_EtOH <- DGEList(counts=Counts_EtOH)
dge_gene_EtOH$samples
# TMM normalization
dge_gene_norm1_EtOH <- calcNormFactors(dge_gene_EtOH, method="TMM")
dge_gene_norm1_EtOH$samples
# group lib.size norm.factors
# LNCaP_EtOH_30min_Rep1_Total     1 72457205       0.7884
# LNCaP_EtOH_30min_Rep2_Total     1 96841575       0.7950
# LNCaP_EtOH_30min_Rep1_TT        1 24036982       1.2814
# LNCaP_EtOH_30min_Rep2_TT        1 30876781       1.2450

dge_gene_5min <- DGEList(counts=Counts_5min)
dge_gene_5min$samples
# TMM normalization
dge_gene_norm1_5min <- calcNormFactors(dge_gene_5min, method="TMM")
dge_gene_norm1_5min$samples
# group lib.size norm.factors
# LNCaP_DHT_05min_Rep1_Total     1 82349709       0.8566
# LNCaP_DHT_05min_Rep2_Total     1 95929073       0.8447
# LNCaP_DHT_05min_Rep1_TT        1 28479917       1.2236
# LNCaP_DHT_05min_Rep2_TT        1 38827864       1.1294

dge_gene_10min <- DGEList(counts=Counts_10min)
dge_gene_10min$samples
# TMM normalization
dge_gene_norm1_10min <- calcNormFactors(dge_gene_10min, method="TMM")
dge_gene_norm1_10min$samples
# group lib.size norm.factors
# LNCaP_DHT_10min_Rep1_Total     1 70628274       0.8010
# LNCaP_DHT_10min_Rep2_Total     1 72798096       0.7797
# LNCaP_DHT_10min_Rep1_TT        1 19851545       1.3262
# LNCaP_DHT_10min_Rep2_TT        1 26802190       1.2074

dge_gene_20min <- DGEList(counts=Counts_20min)
dge_gene_20min$samples
# TMM normalization
dge_gene_norm1_20min <- calcNormFactors(dge_gene_20min, method="TMM")
dge_gene_norm1_20min$samples
# group lib.size norm.factors
# LNCaP_DHT_20min_Rep1_Total     1 66838650       0.8407
# LNCaP_DHT_20min_Rep2_Total     1 73436576       0.8231
# LNCaP_DHT_20min_Rep1_TT        1 12455014       1.3697
# LNCaP_DHT_20min_Rep2_TT        1 42112948       1.0551


dge_gene_30min <- DGEList(counts=Counts_30min)
dge_gene_30min$samples
# TMM normalization
dge_gene_norm1_30min <- calcNormFactors(dge_gene_30min, method="TMM")
dge_gene_norm1_30min$samples
# group lib.size norm.factors
# LNCaP_DHT_30min_Rep1_Total     1 63087286       0.8304
# LNCaP_DHT_30min_Rep2_Total     1 71772197       0.7772
# LNCaP_DHT_30min_Rep1_TT        1 20291097       1.2462
# LNCaP_DHT_30min_Rep2_TT        1 20986088       1.2435
