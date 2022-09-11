### S. Tkachenko 05/21/21 routine given filtered 10x directory, makes raw count matrix
###              assumes usual 10x file names
###              Follows get10xdata.R script
makeRawMat <- function(dataDir,geneID = c("ensembl","symbol"),
                       save2file=T, fileOut="rawCountMatrix.csv"){

    geneID <- match.arg(geneID)

    if(!require("Matrix")){
      install.packages("Matrix")
      require("Matrix")
    }

    geneFile <- paste(dataDir,"features.tsv.gz",sep="/")
    cellFile <- paste(dataDir,"barcodes.tsv.gz",sep="/")
    mtxFile <- paste(dataDir,"matrix.mtx.gz",sep="/")
    #
    if(geneID=="ensembl"){
       geneVec <- read.table(gzfile(geneFile),stringsAsFactors=FALSE)[,1]
    } else {
       geneVec <- read.table(gzfile(geneFile),stringsAsFactors=FALSE)[,2]
    }
    cellVec <- read.table(gzfile(cellFile),stringsAsFactors=FALSE)[,1]
    exprMat <- Matrix::readMM(gzfile(mtxFile))
    ##
    rownames(exprMat) <- geneVec
    colnames(exprMat) <- cellVec
    #
    if(save2file){
      write.csv(as.matrix(exprMat),fileOut)
    }

    return(exprMat)
}