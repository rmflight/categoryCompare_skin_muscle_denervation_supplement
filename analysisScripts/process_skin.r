# Original processing was done using R 2.13.1 and Bioconductor packages from 
# v 2.8. This means that a newer versions of R will produce slightly different
# results, but they should not be significant.

# this processes the Skin data from Jeff Petruska. The GEO identifier is in the
# publication text.

require(affy)
require('rat2302.db')
require(genefilter)
require(limma)
library(affy)
attFile <- 'skinDescriptor.txt'

attData <- read.table(attFile,sep='\t',header=TRUE,stringsAsFactors=FALSE)

#setwd('raw_data')
attData$TimePoint <- 0
attData$TimePoint[attData$Time == '7 day'] <- 7
attData$TimePoint[attData$Time == '14 day'] <- 14

cellDat <- ReadAffy(phenoData=attData)
rmaDat <- rma(cellDat)

rmaFilt <- nsFilter(rmaDat, require.entrez=TRUE, remove.dupEntrez=FALSE, var.filter=FALSE, feature.exclude='^AFFX')
rmaDat <- rmaFilt$eset

f <- (factor(pData(rmaDat)[,'TimePoint']))
design <- model.matrix(~0 + f)
colnames(design) <- c("T0", "T7", "T14")

fit <- lmFit(rmaDat, design)
contrast.matrix <- makeContrasts(T7 - T0, T14 - T0, T14 - T7, levels=design) # not actually using T14-T7, but might later
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
fit2

difProbes_T7vT0 <- topTable(fit2, coef="T7 - T0", adjust.method="BH", p.value=1, number=Inf)
difProbes_T14vT0 <- topTable(fit2, coef="T14 - T0", adjust.method="BH", p.value=1, number=Inf)
difProbes_T14vT7 <- topTable(fit2, coef="T14 - T7", adjust.method="BH", p.value=1, number=Inf)

save('difProbes_T14vT0','difProbes_T14vT7','difProbes_T7vT0',file='skinData.RData')