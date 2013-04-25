# this function looks up all the homologs from one species list of entrez ids to another list of entrez ids from a different species
# note that this function only returns an ID if it has a homolog in the other species of interest. If there is no homolog in the other species, then it will not appear in the list at all.

homoLookup <- function(taxid1,taxid2){

# create a new table that organizes by homologene ID, with the two species in each to generate a table that can be easily used to get the data that one is interested in.
load('H:/datasets/homologene_data/homologeneData.RData')

isTax1 <- homoData$Taxon == taxid1
isTax2 <- homoData$Taxon == taxid2

homo1 <- unique(homoData$HomoGene[isTax1]) # need to worry about multiple homologs
homo2 <- unique(homoData$HomoGene[isTax2])

useHomo <- homo1[homo1 %in% homo2]

homoData <- homoData[((homoData$HomoGene %in% useHomo) & (isTax1 | isTax2)),] # this trims the data set so we don't have to search as many entries

# now start trying to pull stuff together to make an actual homology table.
nHomo <- length(useHomo)
homoGene <- vector('integer',length=0)
homoTaxid1 <- vector('integer',length=0)
homoTaxid2 <- vector('integer',length=0)
for (iHomo in 1:nHomo){
	tmpGene <- useHomo[iHomo]
	tmpHomo <- homoData[homoData$HomoGene == tmpGene,]
	tmpTaxid1 <- homoData$Entrez[(homoData$Taxon == taxid1) & (homoData$HomoGene == tmpGene)]
	tmpTaxid2 <- homoData$Entrez[(homoData$Taxon == taxid2) & (homoData$HomoGene == tmpGene)]
	nTax1 <- length(tmpTaxid1)
	nTax2 <- length(tmpTaxid2)
	
	# we have to have at least one Entrez entry for each Taxon, so if these are not the same length, then we must have multiples for at least one of them and need to do something accordingly (probably make a cross table so every one is cross listed for every other one)
	if (nTax2 != nTax1){
		lenReq <- nTax1 * nTax2
		outTax1 <- vector('integer',0)
		outTax1[1:lenReq] <- tmpTaxid1
		outTax2 <- vector('integer',0)
		outTax2[1:lenReq] <- tmpTaxid2
		outHomo <- vector('integer',0)
		outHomo[1:lenReq] <- tmpGene
		homoTaxid1 <- c(homoTaxid1,outTax1)
		homoTaxid2 <- c(homoTaxid2,outTax2)
		homoGene <- c(homoGene, outHomo)
	} else {
		outHomo <- vector('integer',0)
		outHomo[1:nTax1] <- tmpGene
		homoTaxid1 <- c(homoTaxid1, tmpTaxid1)
		homoTaxid2 <- c(homoTaxid2, tmpTaxid2)
		homoGene <- c(homoGene, outHomo)
	}
	
}

outDat <- data.frame(homoGene,homoTaxid1,homoTaxid2)
names(outDat) <- c('homoGene',as.character(taxid1),as.character(taxid2))
outDat
}