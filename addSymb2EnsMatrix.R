### 10/10/19 utility script for adding gene symbols to matrix
### from csv file with 1st column as row names

source("geneIDconversion.R")

### add symbols to matrix from CSV file
addSymbolsCSV <- function(fileIn,species){

    countmat <- read.table(fileIn,header=T,row.names=1,sep=",")
    outMat <- addSymbolsMatrix(countmat=countmat,species=species)
    fileOut <- paste0(fileIn,"geneSymbAdded.csv")
    write.csv(outMat,fileOut)

}

### add symbols to just matrix
addSymbolsMatrix <- function(countmat,species){

    ## add gene symbols to the matrix as the 2nd column
    geneSymb <- Ens2Sym_ordered(ensemblVec=row.names(countmat),
				species=species)
    countmat <- cbind(countmat,as.data.frame(geneSymb))
    countmat <- countmat[,c(ncol(countmat),1:(ncol(countmat)-1))]
    colnames(countmat)[1] <- "Gene symbol"    	

    return(countmat)
}