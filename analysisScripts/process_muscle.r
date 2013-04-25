## do some raw data analysis on the muscle data from Wang et al.

# The original processing was done in R 2.13.1 and Bioconductor 2.8
# Be aware that slightly different results will be generated using newer versions of R
# and Bioconductor packages.

# download the data, and unTar / unzip it.

# download the supplemental data for GSE4411
# shell.exec("ftp://ftp.ncbi.nih.gov/pub/geo/DATA/supplementary/series/GSE4411/GSE4411%5FRAW%2Etar")

# load the raw data into an AffyBatch object to allow processing
require(affy)
require('mouse4302.db')
require(genefilter)
require(limma)

attFile <- 'muscleDescriptor.txt'
attData <- read.table(attFile,sep='\t',header=TRUE,stringsAsFactors=FALSE)
attData <- attData[,1:2]
attData$Control <- 'MCK'
attData$Control[grep('Control',attData$Type)] <- 'Control'
attData$status <- 'inn'
attData$status[grep('Denervated',attData$Type)] <- 'den'

cellDat <- ReadAffy(phenoData=attData,celfile.path='GSE4411_RAW')

# Perform RMA normalization
rmaDat <- rma(cellDat)
# saveDat('rmaDat',file='rma_data')

rmaFilt <- nsFilter(rmaDat, require.entrez=TRUE, remove.dupEntrez=FALSE, var.filter=FALSE, feature.exclude='^AFFX')
rmaDat <- rmaFilt$eset

sampDat <- pData(rmaDat)
rmaDat <- rmaDat[,(sampDat$Control == 'Control')]

f <- (factor(pData(rmaDat)[,'status']))
design <- model.matrix(~0 + f)
colnames(design) <- c("den", "inn")

fit <- lmFit(rmaDat, design)
contrast.matrix <- makeContrasts(den - inn, levels=design) 
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

outDat <- topTable(fit2, coef='den - inn', number=Inf, adjust.method='BH', sort.by='P', p.value=1)

save(outDat, file='muscleData.RData')