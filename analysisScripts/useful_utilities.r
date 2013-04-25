# # set a bioconductor mirror
# makeActiveBinding("biocLite", local({
# 
#     env <- new.env() 
#     function() {
# 
#         if (!exists("biocLite", envir=env, inherits=FALSE)) {
#             evalq(source("http://bioconductor.org/biocLite.R",
#                          local=TRUE),
#                   env)
#         }
#         env[["biocLite"]]
# 
#     } 
# }), .GlobalEnv)

# some core processing functions for outputs from LIMMA.
outLimma <- function(rmaDat,limFit,adjust.method='BH'){
	
	limProbes <- limFit$genes$ID
	limID <- limFit$genes
	rmaProbes <- featureNames(rmaDat)
	
	rmaOrder <- match(rmaProbes,limProbes,nomatch=0)
	rmaDat <- rmaDat[rmaOrder,]
	
	if (length(annotation(rmaDat)) > 0){
		limID$Entrez <- getAnnData(limID[,'ID'],annotation(rmaDat),"ENTREZID")
		limID$Symbol <- getAnnData(limID[,'ID'],annotation(rmaDat),"SYMBOL")
		limID$Name <- getAnnData(limID[,'ID'],annotation(rmaDat),"GENENAME")
	}
	
	confInt <- sqrt(limFit$s2.post) * limFit$stdev.unscaled * qnorm(0.975)
	colnames(confInt) <- paste(colnames(confInt),'Conf975',  sep='.')
	
	colnames(limFit$coefficients) <- paste(colnames(limFit$coefficients),'logFC',sep='.')
	logFC <- as.data.frame(cbind(limFit$coefficients,confInt))
	logFC$ID <- rownames(logFC)
	
	rawExpr <- exprs(rmaDat)
	uniqCols <- colnames(limFit$design)
	nDCols <- ncol(limFit$design)
	
	aveExpr <- matrix(0,nrow(rawExpr),nDCols)
	stdDev <- matrix(0,nrow(rawExpr),nDCols)
	for (iCol in 1:nDCols){
		aveExpr[,iCol] <- apply(rawExpr[,as.logical(limFit$design[,iCol])],1,mean)
		stdDev[,iCol] <- apply(rawExpr[,as.logical(limFit$design[,iCol])],1,sd)
	}
	colnames(aveExpr) <- paste(uniqCols,'AveExpr',sep='.')
	colnames(stdDev) <- paste(uniqCols,'SD',sep='.')
	rownames(aveExpr) <- rownames(stdDev) <- rownames(rawExpr)
	
	
	adjPVal <- matrix(0,nrow(limFit$p.value),ncol(limFit$p.value))
	for (iCol in 1:ncol(limFit$p.value)){
		adjPVal[,iCol] <- p.adjust(limFit$p.value[,iCol], method=adjust.method)
	}
	rownames(adjPVal) <- rownames(limFit$p.value)
	colnames(adjPVal) <- paste(colnames(limFit$p.value),'Adj.P.Value', sep='.')
	colnames(limFit$p.value) <- paste(colnames(limFit$p.value),'P.Value',sep='.')
	
	statDat <- cbind(aveExpr,stdDev,limFit$p.value,adjPVal)
	statDat <- as.data.frame(statDat)
	statDat$ID <- rownames(statDat)
	
	allDat <- merge(limID,logFC)
	allDat <- merge(allDat,statDat)
	names(allDat) <- gsub(' - ','v',names(allDat))
	names(allDat) <- gsub(' ','.',names(allDat))
	return(allDat)
	
}


getAnnData <- function(id,annPackage,type){
	require("annotate")
  annEnv <- getAnnMap(type, annPackage, load=TRUE)
  unlist(mget(id, annEnv, ifnotfound=NA),use.names=FALSE)
}

saveDat <- function(inList, file=stop("'file' must be specified"), useDate=TRUE, .sessionInfo=TRUE, saveFunc=FALSE){
	if (useDate){
		dateStr <- paste(format(Sys.Date(),"%d%m%y"),'.RData',sep='')
	} else {
		dateStr <- ".RData"
	}
	rdataLoc <- regexpr('RData',file)
	if (rdataLoc != -1){
		file <- substr(file,1,rdataLoc-2)
	}
	fileStr <- paste(file,dateStr,sep="_")
	
	if (!saveFunc){
		
		tmpEnv <- environment()
		passFuncs <- as.vector(lsf.str(envir=parent.env(tmpEnv), inList))
		inList <- inList[!(inList %in% passFuncs)]
		
	}
	
	if (.sessionInfo){
		.sessionInfo <- sessionInfo()
		.sessionInfo$date <- Sys.time()
		inList <- c(inList, ".sessionInfo")
	}
	
	save(list=inList, file=fileStr)
}


# filtering out duplicates from DNA microarray data sets
dupCons <- function(datTable,IDCol,RatCol,ValCol,KeepType){

	switch(KeepType,
		min = tmpOrder <- order(datTable[,ValCol],decreasing=FALSE),
		max = tmpOrder <- order(datTable[,ValCol],decreasing=TRUE)
	)
	
	datTable <- datTable[tmpOrder,]
	isDup <- duplicated(datTable[,IDCol],fromLast=FALSE)
	isDup <- unique(datTable[,IDCol][isDup])
	nDup <- length(isDup)
	
	remDup <- vector('logical',nrow(datTable))
	
	# go through and check the consistency of the duplicates
	for (iDup in 1:nDup){
		tmpDup <- datTable[,IDCol] %in% isDup[iDup]
		dupLoc <- which(tmpDup)
		nDupLoc <- sum(tmpDup)
		dupRat <- datTable[,RatCol][dupLoc]
		
		# check if they are in the same direction, then set the first one to keep and get rid of the others
		if (sum(dupRat > 0) == nDupLoc){
			tmpDup[dupLoc[1]] <- FALSE
		} else if (sum(dupRat < 0) == nDupLoc){
			tmpDup[dupLoc[1]] <- FALSE
		}
		remDup <- remDup | tmpDup	
	}

	datTable <- datTable[!remDup,]
}

### generating feature level comparisons of two experiments
# generate a table of only those features that are present in both experiments
compTableGen <- function(compObj1, compObj2){
	tmpID1 <- lapply(compObj1, function(x){
		x$data[, x$useID]
	})
	allID1 <- unique(unlist(tmpID1))
	
	tmpID2 <- lapply(compObj2, function(x){
		x$data[, x$useID]
	})
	allID2 <- unique(unlist(tmpID2))
	
	keepID1 <- lapply(tmpID1, function(x){
		intersect(x, allID2)
	})
	
	keepID2 <- lapply(tmpID2, function(x){
		intersect(x, allID1)
	})
	
	# now grab the data out of the respective tables, and shove it into a new table for all of them
	newTable <- data.frame(ID=NA, stringsAsFactors=F)
	n1 <- length(compObj1)
	for (i1 in 1:n1){
		tmpDat <- compObj1[[i1]]$data
		useIndx <- (tmpDat[,compObj1[[i1]]$useID] %in% keepID1[[i1]])
		tmpTable <- data.frame(ID=tmpDat[useIndx, compObj1[[i1]]$useID], FC=tmpDat[useIndx, compObj1[[i1]]$useFC], orgID=tmpDat[useIndx, compObj1[[i1]]$orgID], stringsAsFactors=F)
		class(tmpTable$FC) <- "numeric"
		names(tmpTable)[2:3] <- paste(names(compObj1)[i1], names(tmpTable[2:3]), sep=".")
		newTable <- merge(newTable, tmpTable, by="ID", all=T)
	}
	
	n2 <- length(compObj2)
	for (i2 in 1:n2){
		tmpDat <- compObj2[[i2]]$data
		useIndx <- (tmpDat[,compObj2[[i2]]$useID] %in% keepID2[[i2]])
		tmpTable <- data.frame(ID=tmpDat[useIndx, compObj2[[i2]]$useID], FC=tmpDat[useIndx, compObj2[[i2]]$useFC], orgID=tmpDat[useIndx, compObj2[[i2]]$orgID], stringsAsFactors=F)
		class(tmpTable$FC) <- "numeric"
		names(tmpTable)[2:3] <- paste(names(compObj2)[i2], names(tmpTable[2:3]), sep=".")
		newTable <- merge(newTable, tmpTable, by="ID", all=T)
	}
	
	return(newTable)
	
}

compTableGenAll <- function(compObj1, compObj2){
	
	# now grab the data out of the respective tables, and shove it into a new table for all of them
	newTable <- data.frame(ID=NA, stringsAsFactors=F)
	n1 <- length(compObj1)
	for (i1 in 1:n1){
		tmpDat <- compObj1[[i1]]$data
		tmpTable <- data.frame(ID=tmpDat[, compObj1[[i1]]$useID], FC=tmpDat[, compObj1[[i1]]$useFC], orgID=tmpDat[, compObj1[[i1]]$orgID], P.Value=tmpDat[, compObj1[[i1]]$pVal], stringsAsFactors=F)
		class(tmpTable$FC) <- "numeric"
		names(tmpTable)[2:4] <- paste(names(compObj1)[i1], names(tmpTable[2:4]), sep=".")
		newTable <- merge(newTable, tmpTable, by="ID", all=T)
	}
	
	n2 <- length(compObj2)
	for (i2 in 1:n2){
		tmpDat <- compObj2[[i2]]$data
		tmpTable <- data.frame(ID=tmpDat[, compObj2[[i2]]$useID], FC=tmpDat[, compObj2[[i2]]$useFC], orgID=tmpDat[, compObj2[[i2]]$orgID], P.Value=tmpDat[, compObj2[[i2]]$pVal], stringsAsFactors=F)
		class(tmpTable$FC) <- "numeric"
		names(tmpTable)[2:4] <- paste(names(compObj2)[i2], names(tmpTable[2:4]), sep=".")
		newTable <- merge(newTable, tmpTable, by="ID", all=T)
	}
	
	return(newTable)
	
}

# and compare those that are down and up-regulated in all the reported 
createVennUpDown <- function(compObj1, useName1, compObj2, useName2, fileName){
	# take each object, and split it into up and down lists
	UPDN1 <- lapply(compObj1, function(x){ 
		tmpDat <- x$data
		useFC <- x$useFC
		useID <- x$useID
		isUp <- tmpDat[,useFC] > 0
		return(list(UP=tmpDat[isUp, useID], DN=tmpDat[!isUp, useID]))
	})
	
	allUp1 <- sapply(UPDN1, function(x){
		x$UP
	})
	allUp1 <- unique(unlist(allUp1))
	
	allDn1 <- sapply(UPDN1, function(x){
		x$DN
	})
	allDn1 <- unique(unlist(allDn1))
	
	UPDN2 <- lapply(compObj2, function(x){ 
		tmpDat <- x$data
		useFC <- x$useFC
		useID <- x$useID
		isUp <- tmpDat[,useFC] > 0
		return(list(UP=tmpDat[isUp, useID], DN=tmpDat[!isUp, useID]))
	})
	
	allUp2 <- sapply(UPDN2, function(x){
		x$UP
	})
	allUp2 <- unique(unlist(allUp2))
	
	allDn2 <- sapply(UPDN2, function(x){
		x$DN
	})
	allDn2 <- unique(unlist(allDn2))
	
	vennClass <- list(allUp1, allDn1, allUp2, allDn2)
	names(vennClass) <- c(paste(useName1, c("UP", "DN"), sep="."), paste(useName2, c("UP", "DN"), sep="."))
	
	vennCount <- Venn(vennClass)
	vennNames <- colnames(vennCount@IndicatorWeight)[1:ncol(vennCount@IndicatorWeight)-1]
	allVenn <- vennCount@IndicatorWeight[,".Weight"]
	vennCounts <- allVenn[!(allVenn == 0)]

	grpNames <- sapply(names(vennCounts), function(x){
		logIndx <- as.logical(as.numeric(strsplit(x, "")[[1]]))
		paste(vennNames[logIndx], collapse="::")
	})

	vennList <- vennCount@IntersectionSets[!(allVenn == 0)]
	names(vennCounts) <- grpNames
	names(vennList) <- grpNames
	outList <- data.frame(compLabel=names(vennCounts), Count=vennCounts, Genes=sapply(vennList, function(x){paste(x, collapse=", ")}), stringsAsFactors=F)
	#write.table(outList, file="upDown_geneOverlap.txt", sep="\t", row.names=F)
	
	return(list(vennClass=vennClass, outList=outList))
}
