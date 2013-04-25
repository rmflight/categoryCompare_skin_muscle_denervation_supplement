#### 
# Explanatory Notes
####

# This file contains the series of commands used to generate the results in the 
# accompanying publication. This was run using R 2.14.1, and Bioconductor 2.9. If
# it is run using newer versions of R, the results will likely be slightly different,
# however the general results should be the same. Differences arise due to different
# mappings of probeset to Entrez, and new information about gene annotations.
# 
# Example scripts for generating the actual skin and muscle RData files are in 
# "skinProcessing.r" and "muscleProcessing.r", respectively. The original array
# processing was done using R 2.13.0 and Bioconductor 2.8.
# 
# The commands below assume processing was done in one session, without deleting and
# restarting. Therefore, it is also assumed that two instances of Cytoscape are started,
# with different ports used for each instance of CytoscapeRPC. The version of the
# CytoscapeRPC plugin used is also included in this zip file. 

####
# load up required packages
####

source("useful_utilities.r")
source("resaveImages.r")
library(org.Rn.eg.db)
library(org.Mm.eg.db)
library(rat2302.db)
library(GO.db)
library(KEGG.db)
library(categoryCompare)
library(graph)
library(RCytoscape)

####
# Load up the skin and muscle data sets
####
load('skinData.RData') # Rat skin data
load('muscleData.RData') # Mouse muscle data

# Convert affy probe IDs to Entrez Genes
difProbes_T7vT0$Entrez <- getAnnData(difProbes_T7vT0$ID,annPackage='rat2302.db',type='ENTREZID')
difProbes_T14vT0$Entrez <- getAnnData(difProbes_T14vT0$ID,annPackage='rat2302.db',type='ENTREZID')

# get rid of duplicate Entrez, keeping those with the lowest p-value, and also
# deleting entries that show different direction of FC
difProbes_T7vT0 <- dupCons(difProbes_T7vT0,'Entrez','logFC','adj.P.Val','min')
difProbes_T14vT0 <- dupCons(difProbes_T14vT0, 'Entrez', 'logFC', 'adj.P.Val', 'min')
outDat <- dupCons(outDat,'Entrez','denvinn.logFC','denvinn.Adj.P.Value','min')


# convert the rat stuff to mouse, and only keep that which is mapped to homologs
source('homologene_data/homoLookup.r')
homoDat <- homoLookup(10090,10116)

musEntrez <- outDat$Entrez
musHasHomo <- musEntrez %in% homoDat$'10090'
musEntrez <- musEntrez[musHasHomo]
outDat <- outDat[musHasHomo,]

# this function filters out those rat genes that don't have mouse homologs, and replaces the rat Entrez IDs with mouse Entrez IDs
homoFilterConver <- function(datTable,fromSpec,toSpec){
	tmpEntrez <- datTable$Entrez
	tmpFilt <- tmpEntrez %in% (homoDat[[fromSpec]])
	datTable <- datTable[tmpFilt,]
	tmpEntrez <- tmpEntrez[tmpFilt]
	fromTo <- homoDat[[toSpec]][match(tmpEntrez, homoDat[[fromSpec]])]
	datTable$Entrez <- as.character(fromTo)
	return(datTable)
}

difProbes_T7vT0 <- homoFilterConver(difProbes_T7vT0,'10116','10090')
difProbes_T14vT0 <- homoFilterConver(difProbes_T14vT0,'10116','10090')

# because we converted to mouse homologs, we need to filter out duplicates again in the rat data
difProbes_T7vT0 <- dupCons(difProbes_T7vT0,'Entrez','logFC','adj.P.Val','min')
difProbes_T14vT0 <- dupCons(difProbes_T14vT0, 'Entrez', 'logFC', 'adj.P.Val', 'min')

####
# begin analysis
####
# define our p-value cutoff
pCut <- 0.01

# for both data sets, just take a q-value (adj.p.value) cutoff and work with that
dif.dn.skinT7 <- unique(difProbes_T7vT0$Entrez[(difProbes_T7vT0$adj.P.Val <= pCut) & (difProbes_T7vT0$logFC < 0)])
dif.up.skinT7 <- unique(difProbes_T7vT0$Entrez[(difProbes_T7vT0$adj.P.Val <= pCut) & (difProbes_T7vT0$logFC >= 0)])

dif.dn.skinT14 <- unique(difProbes_T14vT0$Entrez[(difProbes_T14vT0$adj.P.Val <= pCut) & (difProbes_T14vT0$logFC < 0)])
dif.up.skinT14 <- unique(difProbes_T14vT0$Entrez[(difProbes_T14vT0$adj.P.Val <= pCut) & (difProbes_T14vT0$logFC >= 0)])

universe.skin <- unique(union(difProbes_T7vT0$Entrez, difProbes_T14vT0$Entrez))

dif.dn.musc <- unique(outDat$Entrez[(outDat$denvinn.Adj.P.Value <= pCut) & (outDat$denvinn.logFC < 0)])
dif.up.musc <- unique(outDat$Entrez[(outDat$denvinn.Adj.P.Value <= pCut) & (outDat$denvinn.logFC >= 0)])
universe.musc <- unique(outDat$Entrez)

# define the gene lists for use in categoryCompare
ccGenes <- list(skin.T7.DN=list(genes=dif.dn.skinT7,
																universe=universe.skin, 
																annotation="org.Mm.eg.db"), 
								skin.T7.UP=list(genes=dif.up.skinT7, 
																universe=universe.skin, 
																annotation="org.Mm.eg.db"), 
								skin.T14.DN=list(genes=dif.dn.skinT14, 
																 universe=universe.skin, 
																 annotation="org.Mm.eg.db"), 
								skin.T14.UP=list(genes=dif.up.skinT14, 
																 universe=universe.skin, 
																 annotation="org.Mm.eg.db"), 
								muscle.DN=list(genes=dif.dn.musc, 
															 universe=universe.musc, 
															 annotation="org.Mm.eg.db"), 
								muscle.UP=list(genes=dif.up.musc, 
															 universe=universe.musc, 
															 annotation="org.Mm.eg.db"))

ccGenes <- new("ccGeneList", ccGenes)
fdr(ccGenes) <- 0 # this doesn't help for this data, so speed it up by setting to 0
ccType(ccGenes) <- c("BP","CC","MF","KEGG")

ccOpts <- new("ccOptions", listNames=names(ccGenes), outType="none")
doComp <- readLines(con=file("listMemberships.txt", open="r+"))
compareNames(ccOpts) <- doComp

saveDat(inList=c("ccGenes", "difProbes_T7vT0", "difProbes_T14vT0", "outDat"), file="ccInData", compress="xz")

ccRes <- ccEnrich(ccGenes)
pvalueType(ccRes) <- 'pval'
minCount(ccRes) <- 5
pvalueCutoff(ccRes) <- 0.001

ccComp <- ccCompare(ccRes, ccOpts)

# due to the crazy number of possible list-memberships, we are going to go into all of the results, and find the actual list-memberships that exist, and redo the ccOptions object with only the ones that are actually present in the data
findAllComps <- function(ccCompCol){
	nType <- length(ccCompCol)
	allComp <- sapply(ccCompCol, function(x){
		tmpGraph <- x@mainGraph
		allNodes <- nodes(tmpGraph)
		unlist(unique(nodeData(tmpGraph, allNodes, "listMembership")))
	})
	allComp <- unique(unlist(allComp, use.names=F))
}

newComp <- findAllComps(ccComp)
newComp
# this is a rearrangement to get the colors in a better space
newComp_re <- c("skin.T7.UP", "skin.T14.UP", "muscle.UP", "skin.T14.DN", "muscle.DN", "skin.T7.UP,muscle.UP", "skin.T14.UP,muscle.UP", "skin.T14.DN,muscle.DN", "skin.T14.DN,muscle.UP")
newComp %in% newComp_re # make sure we got them all

ccOpts_new <- new("ccOptions", listNames=names(ccGenes), compareNames=newComp_re, outType="none")
ccComp_new <- ccCompare(ccRes, ccOpts_new)

# we break the edges to allow fewer edges having to be passed to Cytoscape
ccComp_new$BP <- breakEdges(ccComp_new$BP, 0.8)

####
# Visualize results in Cytoscape
####
# start Cytoscape and CytoscapeRPC plugin with port 9000 prior to the next command
ccBP <- ccOutCyt(ccComp_new$BP, ccOpts_new, postText="ccBP_001", rpcPort=9000)
breakEdges(ccBP, 1) # keep only edges between fully shared annotations

# go through the groups and save them to independent groups
##### 
# NOTE: only the first setup, and one other definition of a group is shown, but this was done individually for each group of nodes with more than two members
# there may be a way to find all of the groups automatically, but that was not done here
ccSavesBP <- cytOutNodes("lipidOxidation.muscleDN", ccBP)
ccSavesBP <- cytOutNodes("singleAnnot.skinT7UP", ccBP, ccSavesBP)

# save everything up to this point so that it can be reloaded and worked with
saveDat(inList=c("ccComp_new", "ccOpts_new", "ccComp", "ccOpts", "ccGenes", "ccBP", "ccSavesBP"))

saveImage.SVG(ccBP, file.path(getwd(), "ccBP_001_5.svg"))
deleteSVGText(svgFile=file.path(getwd(), "ccBP_001_5.svg"))
# and save a legend for it
ccBP_leg <- createLegend(ccOpts_new, "ccBP_001_5_legend")
saveImage.SVG(ccBP_leg, file.path(getwd(), "ccBP_001_5_leg.svg"))

# grab all of the descriptions, split them up at the first ., and then sort 
# them so that we can better group them for discussion
sortCytNodes <- function(cytNodeObj, regSrch){
	allDesc <- sapply(cytNodeObj, function(x){x$descStr})
	dotLoc <- regexpr(regSrch, allDesc)
	condDesc <- substr(allDesc, dotLoc+1, nchar(allDesc))
	condOrd <- order(condDesc)
	allDesc <- allDesc[condOrd]
	
	cytNodeObj <- cytNodeObj[condOrd]
	return(list(nodeObj=cytNodeObj, allDesc=allDesc))
}

sortedBP <- sortCytNodes(ccSavesBP, "\\.")
ccSavesBP <- sortedBP$nodeObj
ccAllDescBP <- sortedBP$allDesc

# now write this out to a file
mergDat <- mergeLists(ccGenes, ccOpts_new)
cytOutData(ccSavesBP, ccComp_new$BP, mergDat, orgType="header", fileName="ccSavesBP_001strict.txt", displayFile=F)

# for each group, lets go in, get the node coordinates, and place a new node at the middle
# of the group with the description so we can easily reference it when creating the legend for the figure
allNodes <- nodes(ccBP@graph)
nodePos <- getPosition(ccBP, allNodes)
xVals <- sapply(nodePos, function(x) x$x)
yVals <- sapply(nodePos, function(x) x$y)

sapply(seq(1, length(ccSavesBP)), function(x){
	grpNodes <- ccSavesBP[[x]]$nodes
	grpX <- mean(xVals[(names(xVals) %in% grpNodes)])
	grpY <- mean(yVals[(names(yVals) %in% grpNodes)])
	nodeName <- paste("descNode", x, sep="")
	addCyNode(ccBP, nodeName)
	ccBP@graph <- addNode(nodeName, ccBP@graph)
	setPosition(ccBP, nodeName, grpX, grpY)
	setNodeAttributesDirect(ccBP, "Desc", "character", nodeName, paste("DE.", allDesc[x], sep=""))
	return(NULL)
})
setNodeLabelRule(ccBP, "Desc")
saveImage.SVG(ccBP, file.path(getwd(), "ccBP_001_5_descStr.svg"))
deleteSVGText.notThis(file.path(getwd(), "ccBP_001_5_descStr.svg"), "DE\\.")

# now do a basic analysis of MF, CC, and KEGG to see if we get confirmation of what we already see
ccMF <- ccOutCyt(ccComp_new$MF, ccOpts_new, "MF.001")
breakEdges(ccMF, 0.8)
ccSavesMF <- cytOutNodes("transport.skinT14UP", ccMF)
ccSavesMF <- cytOutNodes("singleAnnot.skinT7UP", ccMF, ccSavesMF)
sortedMF <- sortCytNodes(ccSavesMF, "\\.")
cytOutData(sortedMF$nodeObj, ccComp_new$MF, mergDat, orgType="header", fileName="ccSavesMF_001strict.txt", displayFile=F)

ccCC <- ccOutCyt(ccComp_new$CC, ccOpts_new, "CC.001")
breakEdges(ccCC, 0.8)
ccSavesCC <- cytOutNodes("mitochondrialMembrane.skinT14DN.muscleDN", ccCC)
ccSavesCC <- cytOutNodes("singleAnnot.skinT7UP", ccCC, ccSavesCC)
sortedCC <- sortCytNodes(ccSavesCC, "\\.")
cytOutData(sortedCC$nodeObj, ccComp_new$CC, mergDat, orgType="header", fileName="ccSavesCC_001strict.txt", displayFile=F)

ccKEGG <- ccOutCyt(ccComp_new$KEGG, ccOpts_new, "KEGG.001")
ccSavesK <- cytOutNodes("oxidativePhosphorylation.skinT14DN.muscleDN", ccKEGG)
ccSavesK <- cytOutNodes("singleAnnot.skinT7UP", ccKEGG, ccSavesK)
sortedK <- sortCytNodes(ccSavesK, "\\.")
cytOutData(sortedK$nodeObj, ccComp_new$KEGG, mergDat, orgType="header", fileName="ccSavesKEGG_001strict.txt", displayFile=F)

#####
# relax constraints
#####

# now what if we relaxed the constraints a little bit
# ccRes2 <- ccRes
# minCount(ccRes2) <- 2
# pvalueCutoff(ccRes2) <- 0.01
# ccComp2 <- ccCompare(ccRes2, ccOpts)
# newComp2 <- findAllComps(ccComp2)
# newComp2 <- newComp2[!(newComp2 %in% "NA")]
# newComp2 <- newComp2[sample(length(newComp2), length(newComp2))]
# ccOpts2_new <- new("ccOptions", listNames=names(ccGenes), compareNames=newComp2, outType="none")
# ccComp2_new <- ccCompare(ccRes2, ccOpts2_new)
# ccComp2_new$BP <- breakEdges(ccComp2_new$BP, 0.8)
# 
# ccBP2 <- ccOutCyt(ccComp2_new$BP, ccOpts2_new, postText="ccBP2_01", rpcPort=9000)
# breakEdges(ccBP2, 1)

## this ends doing comparisons for CCM

#####
# FLA section
#####

####
# only do below if you have deleted the data from memory
####

# load("ccInData_100412.RData")
# pCut <- 0.01

# define the gene lists again if we need to
# dif.dn.skinT7 <- unique(difProbes_T7vT0$Entrez[(difProbes_T7vT0$adj.P.Val <= pCut) & (difProbes_T7vT0$logFC < 0)])
# dif.up.skinT7 <- unique(difProbes_T7vT0$Entrez[(difProbes_T7vT0$adj.P.Val <= pCut) & (difProbes_T7vT0$logFC >= 0)])
# 
# dif.dn.skinT14 <- unique(difProbes_T14vT0$Entrez[(difProbes_T14vT0$adj.P.Val <= pCut) & (difProbes_T14vT0$logFC < 0)])
# dif.up.skinT14 <- unique(difProbes_T14vT0$Entrez[(difProbes_T14vT0$adj.P.Val <= pCut) & (difProbes_T14vT0$logFC >= 0)])
# 
# universe.skin <- unique(union(difProbes_T7vT0$Entrez, difProbes_T14vT0$Entrez))
# 
# dif.dn.musc <- unique(outDat$Entrez[(outDat$denvinn.Adj.P.Value <= pCut) & (outDat$denvinn.logFC < 0)])
# dif.up.musc <- unique(outDat$Entrez[(outDat$denvinn.Adj.P.Value <= pCut) & (outDat$denvinn.logFC >= 0)])
# universe.musc <- unique(outDat$Entrez)


geneListMem <- list(skin.T7.DN=dif.dn.skinT7, 
										skin.T7.UP=dif.up.skinT7, 
										skin.T14.DN=dif.dn.skinT14, 
										skin.T14.UP=dif.up.skinT14, 
										muscle.DN=dif.dn.musc, 
										muscle.UP=dif.up.musc)

doComp <- readLines(con=file("listMemberships.txt", open="r+"))
ffOptions <- new("ccOptions", listNames=names(geneListMem), 
								 compareNames=doComp, outType="none")
# Assign the genes to lists based on the comparisons previously defined for CCM
nodeCompVec <- categoryCompare:::.compMem(geneListMem,ffOptions)
listIDs <- unique(nodeCompVec)
listIDs <- sort(listIDs)
listIDs <- listIDs[listIDs != 0]

ffUniverse <- union(universe.skin, universe.musc)

# create the geneList object for categoryCompare. This analysis is still using FLA,
# as we assigned genes to combinations of lists FIRST, but for the enrichment
# we are going to use the categoryCompare framework to do it, however
ffGenes <- lapply(listIDs, function(x){
	list(genes=names(nodeCompVec)[nodeCompVec == x], universe=ffUniverse, annotation="org.Mm.eg.db")
})

# this allows us to distinguish the combination lists from the comparison lists when we need them
names(ffGenes) <- gsub(",", "..", doComp[listIDs])
ffGenes <- new("ccGeneList", ffGenes)
ccType(ffGenes) <- c("BP","CC","MF","KEGG")
fdr(ffGenes) <- 0

ffRes <- ccEnrich(ffGenes)
minCount(ffRes) <- 5
pvalueType(ffRes) <- 'pval'
pvalueCutoff(ffRes) <- 0.001
pvalueCutoff(ffRes$KEGG) <- 0.05

## change the names of some of the comparisons
tmpCompNames <- names(ffGenes)[7:18]
splitNames <- strsplit(tmpCompNames, '..', fixed=TRUE)
tmpCompNames <- as.list(tmpCompNames)
allComps <- lapply(tmpCompNames, function(x){
	x2 <- strsplit(x, '..', fixed=TRUE)
	x2 <- x2[[1]]
	x3 <- vector("character", 0)
	for (iX in 1:length(x2)){
		x3[iX] <- paste(x,x2[iX],sep=",")
	}
	x3 
})

# these are essentially extra comparisons that need to be done, in the end we will lump them together with the first component
# The reason we are doing these, is because there are a lot of annotations that show up as NA when you only consider the defined
# contrasts using FLA. This is only necessary because we are using the categoryCompare framework for the analysis. If we weren't
# using categoryCompare, we would still see these, and make the same conclusions. The code that was used to determine 
# these were necessary is at the bottom of this script.
extraComps <- unlist(allComps)

# and set up some of the individual 2 way comparisons that should be done
comps2 <- strsplit(names(ffGenes[7:14]), "..", fixed=TRUE)
comps2 <- sapply(comps2, function(x){
	#browser(expr=TRUE)
	paste(x[1], x[2], sep=",")
})

compNames1 <- sample(names(ffGenes), length(ffGenes))
compNames2 <- c(names(ffGenes), comps2, extraComps)

ffOpts1 <- new("ccOptions", listNames=names(ffGenes), compareNames=compNames1, outType="none")
ffOpts2 <- new("ccOptions", listNames=names(ffGenes), compareNames=compNames2, outType="none")
ffComp1 <- ccCompare(ffRes, ffOpts1)
ffComp2 <- ccCompare(ffRes, ffOpts2)

ffListMem <- unique(as.character(findAllComps(ffComp2)))
ffListMem <- ffListMem[!(ffListMem %in% "NA")]
ffListMem <- sample(ffListMem, length(ffListMem))

ffOpts2 <- new("ccOptions", listNames=names(ffGenes), compareNames=ffListMem, outType="none")
ffComp2 <- ccCompare(ffRes, ffOpts2)

ffComp2$BP <- breakEdges(ffComp2$BP, 0.8)
ffBP <- ccOutCyt(ffComp2$BP, ffOpts2, "ffBP.001", rpcPort=9001)
breakEdges(ffBP, 1)

# big question is which ones in ff (FLA) are also in cc (CCM), and which ones are unique. Also, for those that are particularly interesting in cc, what are the p-values in ff
ffNodes <- nodes(ffComp2$BP@mainGraph)
ccNodes <- nodes(ccComp_new$BP@mainGraph)

commNodes <- intersect(ffNodes, ccNodes)
ccSpec <- ccNodes[!(ccNodes %in% commNodes)]
selectNodes(ccBP, ccSpec, preserve.current.selection=F)

# looks as though some of our groups are wholly contained. Can we calculate which ones ARE fully contained by the ccSpec list?
grpPerc <- sapply(ccSavesBP, function(x){
	grpNodes <- x$nodes
	sum(grpNodes %in% ccSpec) / length(grpNodes)
	})
names(grpPerc) <- ccAllDescBP

ccUniqGrps <- names(grpPerc[grpPerc == 1])
ccUniqIndx <- which(ccAllDesc %in% ccUniqGrps)
fileLoc <- file("ccVff_BP_001.txt", open="w+")
ccUniqNodes <- sapply(ccSavesBP[ccUniqIndx], function(x){
	# grab the nodes
	tmpNodes <- x$nodes
	
	# grab a table for those nodes from ccBP
	ccTable <- ccComp_new$BP@mainTable[(ccComp_new$BP@mainTable$ID %in% tmpNodes),]
	# grab the listMembership as well for these
	listMem <- unlist(nodeData(ccComp_new$BP@mainGraph, tmpNodes, "listMembership"))
	ccTable$ListMembership <- listMem
	
	ffTable <- ffComp2$BP@mainTable[(ffComp2$BP@mainTable$ID %in% tmpNodes),]
	
	annotGenes <- unique(unlist(ccComp_new$BP@allAnnotation[tmpNodes], recursive=T, use.names=F))
	mergTable <- mergDat[(mergDat$Entrez %in% annotGenes),]
	
	cat("\n\n", x$descStr, "\n", file=fileLoc)
	cat("ccTable\n", file=fileLoc)
	write.table(ccTable, file=fileLoc, sep="\t", row.names=F)
	cat("\nffTable\n", file=fileLoc)
	write.table(ffTable, file=fileLoc, sep="\t", row.names=F)
	cat("\nAnnotatedGenes\n", file=fileLoc)
	write.table(mergTable, file=fileLoc, sep="\t", row.names=F)
	
})
close(fileLoc)

## get the FLA groups for BP, MF, CC and KEGG
ffSavesBP <- cytOutNodes("carbohydrateMetabolism.muscleDN", ffBP)
ffSavesBP <- cytOutNodes("singleAnnot.skin7UP", ffBP, ffSavesBP)


ffMF <- ccOutCyt(ffComp2$MF, ffOpts2, "ffMF.001", rpcPort=9001)
breakEdges(ffMF, 1)

ffSavesMF <- cytOutNodes("nadBinding.skin14DN.muscleDN", ffMF)
ffSavesMF <- cytOutNodes("singleAnnot.skin7UP", ffMF, ffSavesMF)

ffCC <- ccOutCyt(ffComp2$CC, ffOpts2, "ffCC.001", rpcPort=9001)
breakEdges(ffCC, 1)

ffSavesCC <- cytOutNodes("mitochondrialMembrane.skin14DN.muscleDN", ffCC)
ffSavesCC <- cytOutNodes("singleAnnot.skin7UP.muscleUP", ffCC, ffSavesCC)

ffK <- ccOutCyt(ffComp2$KEGG, ffOpts2, "ffK.001", rpcPort=9001)
breakEdges(ffK, 0.5)
# just select all the nodes, because there are pretty much no groups
ffSavesK <- cytOutNodes("allNodes", ffK)

# and spit everything out into text files
ffMerg <- mergeLists(ffGenes, ffOpts2)

ffSortBP <- sortCytNodes(ffSavesBP, "\\.")
cytOutData(ffSortBP$nodeObj, ffComp2$BP, ffMerg, orgType="header", fileName="ffSavesBP_001strict.txt", displayFile=F)

ffSortMF <- sortCytNodes(ffSavesMF, "\\.")
cytOutData(ffSortMF$nodeObj, ffComp2$MF, ffMerg, orgType="header", fileName="ffSavesMF_001strict.txt", displayFile=F)

ffSortCC <- sortCytNodes(ffSavesCC, "\\.")
cytOutData(ffSortCC$nodeObj, ffComp2$CC, ffMerg, orgType="header", fileName="ffSavesCC_001strict.txt", displayFile=F)

ffSortK <- sortCytNodes(ffSavesK, "\\.")
cytOutData(ffSortK$nodeObj, ffComp2$KEGG, ffMerg, orgType="header", fileName="ffSavesKEGG_001strict.txt", displayFile=F)

## This code was used to determine what to do to get rid of a bunch of "NA" nodes
## that were appearing in the FF results.
## figure out which combinations of lists these NA nodes are significant in, and therefore what list combinations we need to keep them, if we should keep them
## this is what led to figuring out that we needed these extra nodes. 
# tmpNodes <- nodes(ffComp1$BP@mainGraph)
# tmpMember1 <- unlist(nodeData(ffComp1$BP@mainGraph, tmpNodes, 'listMembership'))
# tmpMember2 <- unlist(nodeData(ffComp2$BP@mainGraph, tmpNodes, 'listMembership'))
# sum(tmpMember1 == "NA")
# sum(tmpMember2 == "NA") # so we went from 110 to 87, that's not that impressive. So what the heck is going on
# isFFna <- names(tmpMember2[tmpMember2 == "NA"])
# allTable <- ffComp1$BP@mainTable
# isCount <- grep(".Count", names(allTable), fixed=TRUE)
# isPval <- grep(".Pvalue", names(allTable), fixed=TRUE)
# naTable <- allTable[(allTable$ID %in% isFFna),]
# # naTable[is.na(naTable)] <- 0
# 
# useNames <- names(ffGenes)
# isSigMatrix <- matrix(NA, nrow=length(isFFna), ncol=length(useNames))
# rownames(isSigMatrix) <- isFFna
# colnames(isSigMatrix) <- useNames
# 
# endStr <- c("Pvalue","Count")
# 
# for (nN in 1:length(useNames)){
# 	searchStr <- paste(useNames[nN],endStr, sep=".")
# 	pLoc <- grep(paste("^",searchStr[1],sep=""), names(naTable))
# 	nLoc <- grep(paste("^",searchStr[2],sep=""), names(naTable))
# 	isSig <- (naTable[,pLoc] <= pvalueCutoff(ffRes$BP)) & (naTable[,nLoc] >= minCount(ffRes$BP))
# 	isSigMatrix[,nN] <- isSig
# }
# isSigMatrix[is.na(isSigMatrix)] <- FALSE
# 
# # now tell us which lists were significant for the various GO terms
# sigGO <- apply(isSigMatrix, 1, function(x){
# 	#browser(expr=TRUE)
# 	names(x)[x]
# })
# 
## there seem to be ~30 that are due to being in opposite directions (musc.DN, musc.UP)
## the rest seem to be due to being in an individual list (only the genes from that list) as well as other lists. The question is, do we roll these into the other lists, or do we consider them truly as not belonging? But we found them somewhere, right? i.e. there are GO terms that were found in skin.T7.UP and musc.UP, but not skin.T7.UP..musc.UP, do we combine them all together?
