### S.Tkachenko 1/10/20. Newer version, using column number (rownames=0)
### for adding gene symbols to table with ensembl ids (just table or from csv file)

source("~/UtilScripts/geneIDconversion.R") # core ID conversion

### 1/10/20 adding symbols to table
### Input: table with ensembl IDs
###        column number to convert (can have multiple, simpler this way anyway)
###        species ("human" or "mouse" supported for now, 01/10/20)
### Output: table with gene symbol column added

addSym2EnsTable <- function(dfIn,ensCol,species="human"){

    if(ensCol==0){
        ensVec <- rownames(dfIn)
    } else {
        ensVec <- dfIn[,ensCol]
    }
    ## conversion
    symVec <- Ens2Sym_ordered(ensemblVec=ensVec,species=species,leaveNC.asEnsembl=TRUE,
			      rmDuplSym=FALSE, # not using as row names, not needed
			      localRead=FALSE,localSave=FALSE,localFilename="EnsSymTable.csv")
    ## adding to original table
    dfIn$symbol <- symVec			      

    return(dfIn)
}

### 1/10/20 Uitlizing addSym2EnsTable to add symbols to tables from csv files
### Input: csv file with ensembl IDs
###        column number to convert (can have multiple, simpler this way anyway)
###        species ("human" or "mouse" supported for now, 01/10/20)
### Output: csv file with gene symbol column added

addSym2Ens_csvFile <- function(fileIn,ensCol,species="human"){

    ## get table		   
    ensDF <- read.csv(fileIn)
    print(head(ensDF))		      
    ## convert
    symDF <- addSym2EnsTable(ensDF,ensCol,species)
    print(head(symDF))
    ## write to file
    fileOut <- paste0(fileIn,"_geneSymAdded.csv")
    write.csv(symDF,fileOut)

    return(symDF)
}