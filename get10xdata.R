### S.Tkachenko 11/22/19

### routine getting data from filtered 10x directory
### singled out from qc10x.R and slightly modified

### 07/01/19 function reading filtered_feature_bc_matrix directory; generic
### format is assumed.
### gets data, converts to single cell experiment object
### input arguments: path to 10X data directory OR rdsInFile to read from
###                  saveRDS: if read from directory, to save time in future
###                  rdsOutFile: out file name, not given, = dataDir+".rds"
### output:         SingleCellExperiment object

### 11/23/2020 Added arg (geneID) to use gene names or ensembl IDs as row names 

get10Xdata <- function(inputMatrix=NULL,dataDir=NULL,rdsInFile=NULL,
                       geneID = c("symbol","ensembl"),
                       runName=NULL,saveRDSfile=FALSE,rdsOutFile=NULL){

    ## checking inputs
    if(is.null(dataDir) && is.null(rdsInFile) && is.null(inputMatrix)){
        stop("Need to provide either data directory or RDS file")
    }
    if(!is.null(dataDir) && !is.null(rdsInFile) && !is.null(inputMatrix)){
        stop("Need to provide only one of: data directory, RDS file, inputMatrix")
    }
    if(!is.null(rdsInFile) && saveRDSfile){
        message("Already reading from rds file, ignoring saveRDSfile argument")
    }
    if(!is.null(rdsOutFile) && !saveRDSfile){
        message("saveRDSfile is false, but rdsOutFile given, ignoring and not saving")
    }
#    if(!is.null(rdsInFile) && !is.null(runName)){
#        message("reading from rds file, runName argument ignored")
#    }
    if(!is.null(inputMatrix) && !is.numeric(as.matrix(inputMatrix))){
        stop("Expecting numeric input matrix")
    }
  
    geneID <- match.arg(geneID)
  
    if(!require("Matrix")){
        install.packages("Matrix")
        require("Matrix")
    }

    if(!require(SingleCellExperiment)){
        if (!requireNamespace("BiocManager", quietly = TRUE))
            install.packages("BiocManager")
        BiocManager::install("SingleCellExperiment")
        require(SingleCellExperiment)
    }

    if(!is.null(dataDir)){
        geneFile <- paste(dataDir,"features.tsv.gz",sep="/")
        cellFile <- paste(dataDir,"barcodes.tsv.gz",sep="/")
        mtxFile <- paste(dataDir,"matrix.mtx.gz",sep="/")
        ##
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
        ##
        sceset <- SingleCellExperiment(assays = list(counts = as.matrix(exprMat)))
        ## int_metadata() "is not recommended for common users", but I do not see
        ## where else to add run name       
        int_metadata(sceset)$runName <- runName
        ##
        if(saveRDSfile){
            if(is.null(rdsOutFile)) {
                rdsOutFile <- paste0(dataDir,".rds")
                message("rdsOutFile argument not provided, saving to ",rdsOutFile)
            }
            saveRDS(sceset,rdsOutFile)
        }
    } else if(!is.null(rdsInFile)){
        sceset <- readRDS(rdsInFile)
	int_metadata(sceset)$runName <- runName
    } else if(!is.null(inputMatrix)){
        sceset <- SingleCellExperiment(assays = list(counts = as.matrix(inputMatrix)))
        int_metadata(sceset)$runName <- runName
    } else { #better redundant than sorry
        stop("Need to provide either data directory or RDS file")
    }

    return(sceset)
}
