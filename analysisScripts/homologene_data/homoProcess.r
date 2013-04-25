# homoProcess
# this file takes the homologene.data.txt file and reads it into a data frame in R, gives right titles to the columns, and then saves the data as an RData object for use by homoLookup 

currDir <- getwd()
workDir <- 'H:/datasets/homologene_data'
setwd(workDir)
homoData <- read.table("homologene.data.txt", quote="", header=FALSE, sep="\t", stringsAsFactors=FALSE)
names(homoData) <- c('HomoGene','Taxon','Entrez','Gene.Symbol','ProteinID','ProteinAcc')

dateStr <- as.character(Sys.Date(),'%d%m%y')
save('homoData',file=paste('homologeneData_',dateStr,'.RData',sep=''))
save('homoData',file='homologeneData.RData')