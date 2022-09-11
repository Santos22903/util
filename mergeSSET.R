### 11/26/19 Separated merging scesets to this file. Starts with simonDE
###          routine. Info except counts is discarded (problems with 
###          rowData if doing cbind).

### Input: sets: list of scesets
###        labels for the sets
###        doGeneFiltering - to filter on genes after merging or not
###        signCells - if filtering, in at least many cells more tha 1 gene count

### 11/22/19 merging and gene-filtering data sets
### due to problems with merging scesets, combine count matrices
### then convert to singlecellexperiment object
mergeSSET <- function(sets,labels,doGeneFiltering=TRUE,signCells=2){

    stopifnot(length(sets)==length(labels))

    if(!require('Matrix')){
        install.packages('Matrix', repos='http://cran.us.r-project.org')
        require('Matrix')
    }

    ## renaming columns for distinguishing
    for(ctr in 1:length(sets)){
        colnames(sets[[ctr]]) <- paste(colnames(sets[[ctr]]),labels[ctr],sep="_")
    }
    ## merging count matrices
    mergeTwo<-function(x,y){
        stopifnot(identical(rownames(x),rownames(y)))
        ## saving space?
        mat1 <- Matrix(counts(x),sparse=T)
        mat2 <- Matrix(counts(y),sparse=T)
        return(cbind(mat1,mat2))
    }
    combMat <- Reduce(f=mergeTwo,x=sets)

    ## make SingleCellExperiment object out of the combined matrix
    combSet <- SingleCellExperiment(assays = list(counts = combMat))

    ## since to merge usually QC filtering leaves genes alone, gene-filter here
    if(doGeneFiltering){
        ## no expression
        keep_feature <- rowSums(counts(combSet) > 0) > 0
        combSet <- combSet[keep_feature,]
        ## low expression
        keep_feature <- apply(counts(combSet),1,function(x) length(x[x > 1]) >= signCells)
        combSet <- combSet[keep_feature,]
    }

    return(combSet)
}