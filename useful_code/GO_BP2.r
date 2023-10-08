region_peak <- readPeakFile(paste0("D:/ATAC_brain/mouse/GO_totaldiff5/", region, "_totaldiff.bed"), sep = "")
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
promoter <- getPromoters(TxDb = txdb, upstream = 1000, downstream = 1000)
tagMatrix <- getTagMatrix(region_peak, windows = promoter)
plotAvgProf(
tagMatrix,
xlim = c(-1000, 1000),
xlab = "Genomic Region (5'->3')",
ylab = "Peak Frequency"
)
region_peakAnnolist <- annotatePeak(
region_peak,
tssRegion = c(-1000, 1000),
TxDb = txdb,
annoDb = "org.Mm.eg.db"
)
region_peakAnno <- as.data.frame(region_peakAnnolist)
View(region_peakAnno)
region_biomart <- enrichGO(
gene = region_peakAnno$ENSEMBL,
keyType = "ENSEMBL",
OrgDb = org.Mm.eg.db,
ont = "BP",
pAdjustMethod = "BH",
qvalueCutoff = 0.05,
readable = TRUE
)
