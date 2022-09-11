### 09/10/19 Moved to a separate file in UtilScripts directory
### 05/28/19 Added reading from/saving to local file

##########################################
### 03/12/19 routines for converting ensembl IDs to symbols

possSpecies <- c("human","mouse")

### 06/28/19 core part of affy to symbol conversion
internalAffy2Symbol <- function(affyID_vec,species,platform="affy_mouse430a_2"){
    
    stopifnot(species %in% possSpecies)
    
    if(!require(biomaRt)){
       if (!requireNamespace("BiocManager", quietly = TRUE))
           install.packages("BiocManager")
       BiocManager::install("biomaRt")
       require(biomaRt)
    }	

    if(species=="human"){
        dataSet <- "hsapiens_gene_ensembl"
        geneSymbol <- "hgnc_symbol"      
    } else if(species=="mouse"){
        dataSet <- "mmusculus_gene_ensembl"
        geneSymbol <- "mgi_symbol"
    } else {
        stop("I don't know this species")
    }

    if(platform=="mouse430a2"){
      filter <- "affy_mouse430a_2"
    } else {
      stop("check platform")
    }

    mart <- useMart("ensembl",dataset=dataSet,host = "www.ensembl.org",
                    ensemblRedirect = FALSE)
    geneSymb <- getBM(filters=filter,attributes=c(filter,geneSymbol),
                      values=affyID_vec,mart=mart)

    return(geneSymb)
}

### 10/18/18 Function from scSim.R
### leaveNC.asEnsemble - 09/17/19 addition mainly for using output as row names: instead of
### blanks, leave ensembl IDs where not converted
Ens2Sym_ordered <- function(ensemblVec,species,leaveNC.asEnsembl=TRUE,rmDuplSym=TRUE,
                            localRead=FALSE,localSave=FALSE,localFilename="EnsSymTable.csv"){
    if(!require(biomaRt)){
       if (!requireNamespace("BiocManager", quietly = TRUE))
           install.packages("BiocManager")
       BiocManager::install("biomaRt")
       require(biomaRt)
    }
    
    stopifnot(species %in% possSpecies)

    if(localRead){
        if(file.exists(localFilename)){
            if(localSave){
                message("local file ",localFilename," already present, localSave ignored")
            }
            message("Reading from file ",localFilename)
            returnVec <- SymFromFile(ensemblVec,localFilename)
            if(leaveNC.asEnsembl) returnVec <- ifelse(returnVec=="",ensemblVec,returnVec)
            return(returnVec)
        } else {
            message("local file ",localFilename,"  not found, querying biomaRt") 
        }
    }
    
    symbolTable <- internalEnsemble2Symbol(ensembleID_vec=ensemblVec,
                                     species=species)
    ## possible to return duplicated ensemble IDs w diff symbols
    ## remove second occurences
    symbolTable <- symbolTable[!duplicated(symbolTable$ensembl_gene_id,fromLast=T),]
###    print(head(symbolTable))
    ## remove duplicated symbols
    if(rmDuplSym){
        if(species=="human"){
	    geneSymbol <- "hgnc_symbol"
        } else if(species=="mouse"){
            dataSet <- "mmusculus_gene_ensembl"
            geneSymbol <- "mgi_symbol"
        } else {
            stop("I don't know this species")
        }
	symbolTable <- symbolTable[!duplicated(symbolTable[,geneSymbol],fromLast=T),]	
    }
    ## if some entries skipped, not converted, add dummies
    ids2add <- setdiff(ensemblVec,symbolTable$ensembl_gene_id)
    if(species=="human"){
        df2add <- data.frame(ensembl_gene_id=ids2add,hgnc_symbol=rep("",length(ids2add)))
    } else if(species=="mouse"){
        df2add <- data.frame(ensembl_gene_id=ids2add,mgi_symbol=rep("",length(ids2add)))
    } else {
        stop("I don't know this species")
    }
    symbolTableMod <- rbind(symbolTable,df2add)
    ##
    symbolTableMod <- merge(as.data.frame(ensemblVec),symbolTableMod,
                            by.x="ensemblVec",by.y="ensembl_gene_id",
                            all.x=TRUE,all.y=FALSE, sort=FALSE)

    if(species=="human"){
        returnVec <- symbolTableMod$hgnc_symbol
    } else if(species=="mouse"){
        returnVec <- symbolTableMod$mgi_symbol
    } else {
        stop("I don't know this species")
    }

    ## if requested leaving unconverted as ensemble, do
    if(leaveNC.asEnsembl){
        returnVec <- ifelse(returnVec=="",ensemblVec,returnVec)
    }

    if(localSave){
        df2save <- cbind(ensemblVec,returnVec)
        colnames(df2save) <- c("Ensembl","Symbol")
        write.table(x=df2save,file=localFilename,quote=FALSE,
                    sep=",",row.names=FALSE,col.names=TRUE)
    }

    return(returnVec)
}

### 05/28/19 Routine getting return vector from local file
SymFromFile <- function(ensemblVec,localFilename){

    stopifnot(file.exists(localFilename))
    ## get table from file
    geneTable <- read.csv(localFilename)
    ids2add <- setdiff(ensemblVec,geneTable$Ensembl)
    df2add <- data.frame(Ensembl=ids2add,Symbol=rep("",length(ids2add)))
    geneTableMod <- rbind(geneTable,df2add)

    ## match with input vector
    geneTableMod <- merge(as.data.frame(ensemblVec),geneTableMod,
                            by.x="ensemblVec",by.y="Ensembl",
                            all.x=TRUE,all.y=FALSE, sort=FALSE)
    
    returnVec <- as.character(geneTableMod$Symbol)

    ## return the match
    return(returnVec)
    
}

### 10/18/18 Function from scSim.R
internalEnsemble2Symbol <- function(ensembleID_vec,species){

    stopifnot(species %in% possSpecies)
    
    if(!require(biomaRt)){
       if (!requireNamespace("BiocManager", quietly = TRUE))
           install.packages("BiocManager")
       BiocManager::install("biomaRt")
       require(biomaRt)
    }	

    if(species=="human"){
        dataSet <- "hsapiens_gene_ensembl"
        geneSymbol <- "hgnc_symbol"      
    } else if(species=="mouse"){
        dataSet <- "mmusculus_gene_ensembl"
        geneSymbol <- "mgi_symbol"
    } else {
        stop("I don't know this species")
    }

    mart <- useMart("ensembl",dataset=dataSet,host = "www.ensembl.org",
                    ensemblRedirect = FALSE)
    geneSymb <- getBM(filters="ensembl_gene_id",attributes=c("ensembl_gene_id",geneSymbol),
                      values=ensembleID_vec,mart=mart)
    return(geneSymb)
}

### 08/07/20 mouse to human ortholog conversion
mouse2human <- function(mouseSymbVec){

   if(!require(biomaRt)){
       if (!requireNamespace("BiocManager", quietly = TRUE))
           install.packages("BiocManager")
       BiocManager::install("biomaRt")
       require(biomaRt)
   }
   human = useMart("ensembl", dataset = "hsapiens_gene_ensembl",host = "www.ensembl.org",
                    ensemblRedirect = FALSE)
   mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl",host = "www.ensembl.org",
                    ensemblRedirect = FALSE)
   genesV2 = getLDS(attributes = c("mgi_symbol"),filters = "mgi_symbol", 
                    values = mouseSymbVec, mart = mouse, attributesL = c("hgnc_symbol"), 
		    martL = human, uniqueRows=T)

   return(genesV2)
}

########################################################################
