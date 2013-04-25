# Did something foolish and didn't save the Cytoscape network to a file, and now want to regenerate all of the images of the selected nodes based on a cytOutNodes object. Take the cytOutNodes object, go through each of the list items, and re-save the images in the correct directory

# use the new SVG image export in CyRPC
saveImage.SVG <- function(cw, fileName){
	result <- xml.rpc (cw@uri, 'Cytoscape.exportViewToSVG', cw@window.id, fileName)
	invisible(result)
}

cytOutNodes.ResaveImages <- function(cytOutNodeObj, cytWindowObject, saveDir){
	nOut <- length(cytOutNodeObj)
	allNodes <- getAllNodes(cytWindowObject)
	
	for (iOut in 1:nOut){
		currNodes <- cytOutNodeObj[[iOut]]$nodes
		nodeIntersect <- intersect(allNodes, currNodes)
		if (length(nodeIntersect) == length(currNodes)){
			selectNodes(cytWindowObject, currNodes, preserve.current.selection=F)
			fitSelectedContent(cytWindowObject)
			clearSelection(cytWindowObject)
			descStr <- cytOutNodeObj[[iOut]]$descStr
			fileName <- file.path(saveDir, paste(descStr, "svg", sep="."))
			saveImage.SVG(cytWindowObject, fileName)
		} else { warning(paste("Missing some nodes on ", as.character(iOut), ", image not written!", collapse="", sep=""))}
	}
}

# also set up a function to modify the SVG to remove the ID text from the SVG images
deleteSVGText <- function(svgFile){
	baseName <- substr(svgFile, 1, nchar(svgFile)-4)
	outName <- paste(baseName, "noText.svg", sep="_")
	fileCon <- file(svgFile, open="r+")
	svgText <- readLines(fileCon)
	close(fileCon)
	
	hasText <- grep("<text", svgText, fixed=TRUE)
	keepSVG <- rep(TRUE, length(svgText))
	keepSVG[hasText] <- FALSE
	svgText <- svgText[keepSVG]
	fileOut <- file(outName, open="w+")
	writeLines(svgText, fileOut)
	close(fileOut)
}

deleteSVGText.notThis <- function(svgFile, notThis){
	baseName <- substr(svgFile, 1, nchar(svgFile)-4)
	outName <- paste(baseName, "notThis.svg", sep="_")
	fileCon <- file(svgFile, open="r+")
	svgText <- readLines(fileCon)
	close(fileCon)
	
	hasNot <- grep(notThis, svgText)
	
	hasText <- grep("<text", svgText, fixed=TRUE)
	keepSVG <- rep(TRUE, length(svgText))
	keepSVG[hasText] <- FALSE
	keepSVG[hasNot] <- TRUE
	svgText <- svgText[keepSVG]
	fileOut <- file(outName, open="w+")
	writeLines(svgText, fileOut)
	close(fileOut)
}

# create a color scheme for a given ccOptions object, and then spit it out to Cytoscape
createLegend <- function(inOptions, legendName, ...){
	listMem <- compareNames(inOptions)
	listCol <- compareColors(inOptions)
	#names(listCol) <- NULL
		
	allNodes <- listMem
	gLegend <- new("graphNEL", nodes=allNodes)
	cw <- CytoscapeWindow(legendName, graph=gLegend, ...)
	displayGraph(cw)
	layout(cw, "grid")
	setVal <- sapply(listMem, function(x){
		setNodeColorDirect(cw, x, listCol[x])
	})
	redraw(cw)
	return(cw)
}