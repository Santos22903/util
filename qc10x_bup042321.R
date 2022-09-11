### S.Tkachenko 07/01/19

### "generic" routine for performing basic QC on 10X data
### using filtered_feature_bc_matrix as input. Written generalizing
### routines from olgaSeq.R and erzGR.R. Manual tinkering will still be needed
### for applying manual cuts

### 03/03/20 add outDir (change to outPlotDir) and dirOut (change to outSetDir) as input args

### 02/12/20 Commented out get10X routine in lieu of calling one from 
###          UtilScripts

### 11/05/19 Added option to use matrix provided as argument

### 09/24/19 Added doGeneFiltering for cases when will combine
###          with other data sets later (DE)

### 09/19/19 Beautified histograms by changing main/axis labels;
###          Changed "automatic" to "PCA" on the Venn diagram 

### 07/01/19 main routine calling
### 1) function getting the data
### 2) qcFilter function

### findManCuts means: "just check, do not do cutting in  sceset"
###           produces pics for distributions for manual cuts and
###           compares with "auto" approach
###           checkManCuts and doCuts (poss redundant) check - drawcuts
###           on figures, no "auto"; do - cut sceset

source("~/UtilScripts/get10xdata.R")

possSpecies <- c("mouse","human")
###outDir <- "filterOutput"
###tablesOut <- FALSE

### signCells - # of cells with given gene count >= 1
runQC <- function(inputMatrix=NULL,                   # args for getData
                  dataDir=NULL,rdsInFile=NULL,        # args for getData
                  saveRDSfile=FALSE,rdsOutFile=NULL,  # args for getData
                  runName="10Xrun",                   # args for getData
                  geneID=c("ensembl","symbol"),       # args for qcFilter
                  doGeneFiltering=TRUE,               # args for qcFilter
                  signCells=2,                        # args for qcFilter
                  species="human",addSpikeGenes="",   # args for qcFilter
                  manFiltFlag=TRUE,                   # args for qcFilter
                  autoFiltFlag=!manFiltFlag,          # args for qcFilter
                  compareFiltFlag=FALSE,              # args for qcFilter
                  findManCuts=FALSE, doCuts=FALSE,    # args for qcFilter
                  checkManCuts=FALSE,countCut=NULL,   # args for qcFilter
                  geneCut=NULL, mitoCut=NULL,         # args for qcFilter
		  hiGeneCut=Inf,                      # args for qcFilter
                  tablesOut=TRUE,                     # send cut table to file
                  outPlotDir=NULL,                    # for saving plots
                  outSetDir=NULL,                     # for saving cut sets
		  saveSSET=TRUE){                     # if doCuts, save result
    

    ## checking input
    stopifnot(species %in% possSpecies)
    if(findManCuts && (checkManCuts || doCuts)){
        stop("findManCuts is not compatible with checkManCuts or doCuts")
    }
    if((checkManCuts || doCuts) && (manFiltFlag || compareFiltFlag)
       && is.null(countCut) && is.null(countCut) && is.null(countCut)){
        stop("Please, provide cut values")
    }
    if(manFiltFlag && autoFiltFlag)
        stop("Either manual filtering or automatic can be requested")
    if(is.null(dataDir) && is.null(rdsInFile) && is.null(inputMatrix)){
        stop("Need to provide either data directory or RDS file or input matrix")
    }
    if(!is.null(dataDir) && !is.null(rdsInFile) && !is.null(inputMatrix)){
        stop("Need to provide only one of: data directory, RDS file, inputMatrix")
    }
    if(!is.null(inputMatrix) && !is.numeric(as.matrix(inputMatrix))){
        stop("Expecting numeric input matrix")
    }

    geneID <- match.arg(geneID)
    if(tablesOut || findManCuts || checkManCuts){
        if(!dir.exists(outPlotDir)){
            dir.create(outPlotDir,recursive=T)
        }
    }
    
    ## read in data (data dir or rds file), create a SCE object
    sceset <- get10Xdata(inputMatrix=inputMatrix,
                         dataDir=dataDir,rdsInFile=rdsInFile,
			 geneID=geneID, runName=runName,
                         saveRDSfile=saveRDSfile,rdsOutFile=rdsOutFile)
    message("dim(sceset)")
    print(dim(sceset))
    gc()
    ## use sceset for filtering
    sceset <- qcFilter(sceset=sceset, species=species, addSpikeGenes=addSpikeGenes,
                       geneID=geneID, doGeneFiltering=doGeneFiltering, signCells=signCells,
                       manFiltFlag=manFiltFlag,autoFiltFlag=autoFiltFlag,
                       compareFiltFlag=compareFiltFlag,
                       findManCuts=findManCuts, doCuts=doCuts, checkManCuts=checkManCuts,
                       countCut=countCut, geneCut=geneCut,mitoCut=mitoCut,hiGeneCut=hiGeneCut,
                       tablesOut=tablesOut,outPlotDir=outPlotDir)

    ## puny attempt to help R shitty memory handling
    gc()

    ## if saveSSET, save to file
    if(doCuts && saveSSET){
        saveSet(sceset,runName,outSetDir=outSetDir)
    }
    
    return(sceset)
}

### 11/25/19 Saving sceset to rds file
saveSet <- function(sceset,runName,outSetDir=NULL){

    if(is.null(outSetDir)){		
        outSetDir <- "Output/QCfiltered"
    }
    if(!dir.exists(outSetDir)){dir.create(outSetDir,recursive=T)}
    fileOut <- paste0(outSetDir,"/",runName,"_filtered.rds")
    saveRDS(object=sceset,file=fileOut)

    return(0)
}

### 07/01/19 function performing filtering
###          Input: raw sceset
###                 species (only human or mouse for now, 07/01/19)
###                 addSpikeGenes-genes to add to mito to use as "spikes" (e.g. eyfp)
###          Output: filtered sceset
qcFilter <- function(sceset,species="human", addSpikeGenes="",
                     geneID=c("ensembl","symbol"), doGeneFiltering=TRUE, signCells=2,
                     manFiltFlag=TRUE,autoFiltFlag=FALSE,compareFiltFlag=FALSE,
                     findManCuts=FALSE, doCuts=FALSE, checkManCuts=FALSE,
                     countCut=NULL, geneCut=NULL, mitoCut=NULL,hiGeneCut=Inf,
		     tablesOut=TRUE,outPlotDir){

    stopifnot(species %in% possSpecies)
    if(findManCuts && (checkManCuts || doCuts)){
        stop("findManCuts is not compatible with checkManCuts || doCuts")
    }
    if(checkManCuts && is.null(countCut) && is.null(geneCut)#hiGeneCut can be NULL
       && is.null(mitoCut)){
        stop("Please, provide cut values for checking")
    }
    if((checkManCuts || doCuts) && (manFiltFlag || compareFiltFlag)
       && is.null(countCut) && is.null(geneCut) && is.null(mitoCut)){
        stop("Please, provide cut values")
    }
    geneID <- match.arg(geneID)
    
    if(manFiltFlag && autoFiltFlag)
        stop("Either manual filtering or automatic can be requested")

    if(!require(SingleCellExperiment)){
        if (!requireNamespace("BiocManager", quietly = TRUE))
            install.packages("BiocManager")
        BiocManager::install("SingleCellExperiment")
        require(SingleCellExperiment)
    }
    if(!require(scater)){
        if (!requireNamespace("BiocManager", quietly = TRUE))
            install.packages("BiocManager")
        BiocManager::install("scater",dependencies=TRUE)
        require(scater)
    }
    run <- int_metadata(sceset)$runName
    if(tablesOut){
        runDir <- paste(outPlotDir,run,sep="/")
        if(!dir.exists(runDir)) dir.create(runDir,recursive=T)
        today <- Sys.Date()
        filterFile <- paste0(runDir,"/geneFilters",run,"_",today,".txt")
    }
### following "the tutorial" this cut is done first even if doCuts=F
    ## remove genes not expressed anywhere

    if(doGeneFiltering){
        keep_feature <- rowSums(counts(sceset) > 0) > 0 
        rows1 <- nrow(assay(sceset))
        sceset <- sceset[keep_feature,]
        rows2 <- nrow(assay(sceset))
        filterMessage <- paste0(">>>Removing genes not expressed in any cell: filtered out ",rows1-rows2," genes out of ",rows1,"\n")
        message(filterMessage)
        if(tablesOut){
            cat(filterMessage,file=filterFile)
        }
    }
###

    ## mitochondrial genes (spike-in substitutes)
    if(geneID=="ensembl"){
        if(species=="mouse"){
            mtList <- c("ENSMUSG00000064341","ENSMUSG00000064345","ENSMUSG00000064351",
                        "ENSMUSG00000064354","ENSMUSG00000064356","ENSMUSG00000064357",
                        "ENSMUSG00000064358","ENSMUSG00000064360","ENSMUSG00000065947",
                        "ENSMUSG00000064363","ENSMUSG00000064367","ENSMUSG00000064368",
                        "ENSMUSG00000064370")
        } else if(species=="human"){
            mtList <- c("ENSG00000198888","ENSG00000198763","ENSG00000198804",
                        "ENSG00000198712","ENSG00000228253","ENSG00000198899",
                        "ENSG00000198938","ENSG00000198840","ENSG00000212907",
                        "ENSG00000198886","ENSG00000198786","ENSG00000198695",
                        "ENSG00000198727")
        } else {stop("Unknown species")}
    } else if(geneID=="symbol"){
        mtList <- grep("^MT-",rownames(sceset),value=T,ignore.case=T)
    } else {stop("Unknown geneID")}
    mtList <- c(mtList,addSpikeGenes)
    ##
    isSpike(sceset,"MT") <- rownames(sceset) %in% mtList
    message(">>>Found ",length(which(isSpike(sceset,"MT")))," mitochondrial genes")

### Calculating QC metrics using scater package (lots of stuff to use later)

    sceset <- calculateQCMetrics(
        sceset,
        feature_controls = list(
            MT = isSpike(sceset,"MT")
        )
    )
    gc()
############# TEMPORARY
### checking if removing + tails affects PCA-manual concordance
### for 194698_1 and 194698_2
### cutoffs: _1: 25000/3000; _2: counts=35000, genes=3000
    ## count_cut_1 <- 25000; gene_cut_1=3000
    ## count_cut_2 <- 35000; gene_cut_2=3000
    ## keep_cell <- colData(sceset)[,"total_counts"]<count_cut_1 & colData(sceset)[,"total_features_by_counts"]<gene_cut_1
    ## sceset <- sceset[,keep_cell]

### pictures of all 3 measures (2 pics for mito: 1d&2d),
### if checkManCuts TRUE, add lines at cuts
    if(findManCuts || checkManCuts){
        plotMeasures(sset=sceset,plotCuts=checkManCuts,
                     countCut=countCut,geneCut=geneCut, mitoCut=mitoCut,hiGeneCut=hiGeneCut,
                     outPlotDir=outPlotDir)
    }
    ## if cutting is needed, do it

    if(doCuts || compareFiltFlag){
        if(manFiltFlag || compareFiltFlag){
            if(is.null(countCut) || is.null(geneCut) || is.null(mitoCut)){
                stop("manual filtering is requested, but no cuts specified")
            }
            sceset <- manualFilter(sset=sceset,tablesOut=tablesOut,
                                   countCut=countCut,geneCut=geneCut, mitoCut=mitoCut,hiGeneCut=hiGeneCut,
				   outPlotDir=outPlotDir)
        }

        if(autoFiltFlag || compareFiltFlag){
            sceset <- autoFilter(sset=sceset,tablesOut=tablesOut,outPlotDir=outPlotDir)
        }

        if(compareFiltFlag){
            compareFilters(sceset,outPlotDir=outPlotDir)
        }

        if(doCuts){
            if(manFiltFlag){
                if(doGeneFiltering){
                    ## gene level filtering, done on cells deemed good to avoid
                    ## erroneous expression from "bad" cells
                    filter_genes <- apply(counts(sceset[,colData(sceset)$use]),1,function(x) length(x[x > 1]) >= signCells)
                    rowData(sceset)$use <- filter_genes
                    if(tablesOut){
                        cat(">>>Manual filtering genes with expression>1 in less than 2 cells\n",
                            file=filterFile,append=TRUE)
                        filterMessage <- paste0("Total genes ",length(rowData(sceset)$use),"\n")
                        cat(filterMessage,file=filterFile,append=TRUE)
                        filterMessage <- paste0("Using ",sum(rowData(sceset)$use)," genes\n")
                        cat(filterMessage,file=filterFile,append=TRUE)
                        filterMessage <- paste0("Removing ",sum(!rowData(sceset)$use)," genes\n")
                        cat(filterMessage,file=filterFile,append=TRUE)
                    }
                    sceset <- sceset[rowData(sceset)$use, colData(sceset)$use]
                } else {
                    sceset <- sceset[, colData(sceset)$use]
                }
            } else if(autoFiltFlag){
                if(doGeneFiltering){
                    filter_genes <- apply(counts(sceset[,!colData(sceset)$outlier]),1,function(x) length(x[x > 1]) >= signCells)
                    rowData(sceset)$use <- filter_genes
                    cat(">>>Automatic filtering genes with expression>1 in less than 2 cells\n",
                        file=filterFile,append=TRUE)
                    filterMessage <- paste0("Total genes ",nrow(sset),"\n")
                    cat(filterMessage,file=filterFile,append=TRUE)
                    filterMessage <- paste0("Using ",sum(rowData(sceset)$use)," genes\n")
                    cat(filterMessage,file=filterFile,append=TRUE)
                    filterMessage <- paste0("Removing ",sum(!rowData(sceset)$use)," genes\n")
                    cat(filterMessage,file=filterFile,append=TRUE)
                    sceset <- sceset[rowData(sceset)$use, !colData(sceset)$outlier]
                } else {
                    sceset <- sceset[, !colData(sceset)$outlier]
                }
            } else stop("Need either manualFilter or autoFilter")
        }
    }

    ## final strokes b4 returning (following tutorial)
    assay(sceset, "logcounts_raw") <- log2(counts(sceset) + 1) # adding log-count slot
    
    return(sceset)
}

### 07/17/19 "Comparing" manual and PCA filtering
compareFilters <- function(sset,outPlotDir){

    if(!require(SingleCellExperiment)){
        if (!requireNamespace("BiocManager", quietly = TRUE))
            install.packages("BiocManager")
        BiocManager::install("SingleCellExperiment")
        require(SingleCellExperiment)
    }
    if(!require(scater)){
        if (!requireNamespace("BiocManager", quietly = TRUE))
            install.packages("BiocManager")
        BiocManager::install("scater",dependencies=TRUE)
        require(scater)
    }
    
    ## some preparation
    run <- int_metadata(sset)$runName
    runDir <- paste(outPlotDir,run,sep="/")
    if(!dir.exists(runDir)) dir.create(runDir,recursive=T)
    today <- Sys.Date()
    
    ## plotting reduced dimensions like in autoFilter, but also
    ## overlaying manual filter "use"
    fileOut <- paste0(runDir,"/man_auto_PCA_plot_run_",run,"_",today,".png")
    png(fileOut)
    pcaplot <- plotReducedDim(sset,
                              use_dimred = "PCA_coldata",
                              size_by = "total_features_by_counts",
                              shape_by = "use", 
                              colour_by = "outlier"
                              )
    print(pcaplot)
    dev.off()

    ## venn diagram of manual vs PCA filters
    if(!require(limma)){
        if (!requireNamespace("BiocManager", quietly = TRUE))
            install.packages("BiocManager")
        BiocManager::install("limma",dependencies=TRUE)
        require(limma)
    }
    auto <- colnames(sset)[sset$outlier]
    man <- colnames(sset)[!sset$use]
    venn.diag <- vennCounts(
        cbind(colnames(sset) %in% auto,
              colnames(sset) %in% man)
    )
    fileOut <- paste0(runDir,"/use_PCA_venn_run_",run,"_",today,".png")
    png(fileOut)
    vennDiagram(
        venn.diag,
        names = c("PCA", "Manual"),
        circle.col = c("blue", "green")
    )
    dev.off()
}

### 07/16/19 "automatic" filtering using scater PCA approach
autoFilter <- function(sset,tablesOut=TRUE,outPlotDir){

    if(!require(scater)){
        if (!requireNamespace("BiocManager", quietly = TRUE))
            install.packages("BiocManager")
        BiocManager::install("scater",dependencies=TRUE)
        require(scater)
    }
     if(!require(mvoutlier)){
        if (!requireNamespace("BiocManager", quietly = TRUE))
            install.packages("BiocManager")
        BiocManager::install("mvoutlier",dependencies=TRUE)
        require(mvoutlier)
    }

###    sset <- runPCA(sset,scale=TRUE,use_coldata=TRUE,detect_outliers=TRUE)
###    sset <- runPCA(sset,use_coldata=TRUE,detect_outliers=TRUE,
###                   selected_variables=c("pct_counts_MT","total_features_by####_counts","total_counts"))
    sset <- runColDataPCA(x = sset,
			  variables=c("pct_counts_MT","total_features_by_counts","total_counts"),
			  outliers = TRUE)

    ## plot of 2 dimensions (small subroutine, keep all here)
    run <- int_metadata(sset)$runName
    runDir <- paste(outPlotDir,run,sep="/")
    if(!dir.exists(runDir)) dir.create(runDir,recursive=T)
    today <- Sys.Date()
    fileOut <- paste0(runDir,"/scaterPCA_run_",run,"_",today,".png")
    png(fileOut)
    pcaplot <- plotReducedDim(sset,
                              dimred = "PCA_coldata",
                              size_by = "total_features_by_counts",
                              colour_by = "outlier"
                              )
    print(pcaplot)
    dev.off()
    ## table of auto cuts
    if(tablesOut){
        fileOut <- paste0(runDir,"/autoCutTable_run",run,"_",Sys.Date(),".txt")
        sink(fileOut)
        cat("table(sset$outlier)\n")
        print(table(sset$outlier))
        sink()
    }

    return(sset)
}

### 07/05/19 manual filter using the 3 measures
### 07/16/19 don't do cuts here, just define $use
###          to accomodate manual/auto choice
manualFilter <- function(sset,tablesOut=TRUE,
                         countCut=NULL,geneCut=NULL,mitoCut=NULL,hiGeneCut=Inf,
			 outPlotDir){

    if(is.null(countCut) || is.null(geneCut) || is.null(mitoCut)){
        stop("filtering is requested, but no cuts specified")
    }
    
    filter_by_total_counts <- (sset$total_counts > countCut)
    filter_by_gene_counts <- (sset$total_features_by_counts > geneCut) & (sset$total_features_by_counts < hiGeneCut)
    filter_by_MT <- (sset$pct_counts_MT < mitoCut)

    sset$use <- filter_by_total_counts &
        filter_by_gene_counts &
        filter_by_MT


    ## tables of manual cuts
    run <- int_metadata(sset)$runName
    if(tablesOut){
        runDir <- paste(outPlotDir,run,sep="/")
        if(!dir.exists(runDir)) dir.create(runDir,recursive=T)
        fileOut <- paste0(runDir,"/manualCutTables_run",run,"_",Sys.Date(),".txt")
        sink(fileOut)        
        cat("table(filter_by_total_counts)\n")
        print(table(filter_by_total_counts))
        cat("table(filter_by_gene_counts)\n")
        print(table(filter_by_gene_counts))
        cat("table(filter_by_MT)\n")
        print(table(filter_by_MT))
        cat("table(sset$use)\n")
        print(table(sset$use))
        sink()
    }
    return(sset)
}

### 07/05/19 plotting 3 measures with or wout cuts (dep on checkManCuts flag
###       of the main routing, called here as plotCuts)
plotMeasures <- function(sset,plotCuts=TRUE,countCut=NULL,
                         geneCut=NULL, mitoCut=NULL, hiGeneCut=Inf,
                         outPlotDir){

    if(plotCuts && is.null(countCut) && is.null(geneCut) && is.null(mitoCut)){
        message("plotting cuts requested, but no cuts given, ignoring")
        plotCuts <- FALSE
    }


    run <- int_metadata(sset)$runName
    runDir <- paste(outPlotDir,run,sep="/")
    if(!dir.exists(runDir)) dir.create(runDir,recursive=T)
    today <- Sys.Date()
### library size
    if(plotCuts && !is.null(countCut)){
        fileOut <- paste0(runDir,"/librarySizeHist_run_",run,"_cut",countCut,"_",today,".png")
    } else {
        fileOut <- paste0(runDir,"/librarySizeHist_run_",run,"_",today,".png")
    }
    png(fileOut)
    hist(sset$total_counts,breaks=1000,main="Library size",
         xlab="Number of reads",cex.main=2.0,cex.lab=1.5,cex.axis=1.2)
    if(plotCuts && !is.null(countCut)) abline(v = countCut, col = "red")
    dev.off()
###
### detected genes
    if(plotCuts && !is.null(geneCut)){
        fileOut <- paste0(runDir,"/detectedGenesHist_run_",run,"_cut",geneCut,
	                  "hiCut",hiGeneCut,"_",today,".png")
    } else {
        fileOut <- paste0(runDir,"/detectedGenesHist_run_",run,"_",today,".png")
    }
    png(fileOut)
    hist(sset$total_features_by_counts,breaks=1000,main="Detected genes",
         xlab="Number of genes",cex.main=2.0,cex.lab=1.5,cex.axis=1.2)
    if(plotCuts && !is.null(geneCut)) abline(v = geneCut, col = "red")
    if(plotCuts && !is.infinite(hiGeneCut)) abline(v = hiGeneCut, col = "red")
    dev.off()
###
### mitochondrial genes
    if(plotCuts && !is.null(mitoCut)){
        fileOut <- paste0(runDir,"/mitoGenesHist_run_",run,"_cut",mitoCut,"_",today,".png")
    } else {
        fileOut <- paste0(runDir,"/mitoGenesHist_run_",run,"_",today,".png")
    }
    png(fileOut)
    par(mfrow=c(1,2))
    ## 1d part
    hist(sset$pct_counts_MT,breaks=1000,main='',xlab="Mitochondrial read %",
         cex.lab=1.5,cex.axis=1.2)
    if(plotCuts && !is.null(mitoCut)) abline(v = mitoCut, col = "red")
    ## 2d part
    if(plotCuts && !is.null(mitoCut)){
        hcols<-ifelse(sset$pct_counts_MT>mitoCut,"red","black")
    } else {
        hcols <- "black"
    }
    plot(sset$total_features_by_counts,sset$pct_counts_MT,col=hcols,
         xlab="Detected genes",ylab="Mitochondrial read %",
         cex.axis=1.2,cex.lab=1.5)
    mtext("Mitochondrial genes", outer = TRUE, cex = 2.2,
          side = 3, line = -3)
    dev.off()
}

### 07/01/19 function reading filtered_feature_bc_matrix directory; generic
### format is assumed.
### gets data, converts to single cell experiment object
### input arguments: path to 10X data directory OR rdsInFile to read from
###                  saveRDS: if read from directory, to save time in future
###                  rdsOutFile: out file name, not given, = dataDir+".rds"
### output:         SingleCellExperiment object

#get10Xdata <- function(inputMatrix=NULL,dataDir=NULL,rdsInFile=NULL,
#                       runName=NULL,saveRDSfile=FALSE,rdsOutFile=NULL){
#
#    ## checking inputs
#    if(is.null(dataDir) && is.null(rdsInFile) && is.null(inputMatrix)){
#        stop("Need to provide either data directory or RDS file")
#    }
#    if(!is.null(dataDir) && !is.null(rdsInFile) && !is.null(inputMatrix)){
#        stop("Need to provide only one of: data directory, RDS file, inputMatrix")
#    }
#    if(!is.null(rdsInFile) && saveRDSfile){
#        message("Already reading from rds file, ignoring saveRDSfile argument")
#    }
#    if(!is.null(rdsOutFile) && !saveRDSfile){
#        message("saveRDSfile is false, but rdsOutFile given, ignoring and not saving")
#    }
#    if(!is.null(rdsInFile) && !is.null(runName)){
#        message("reading from rds file, runName argument ignored")
#    }
#    if(!is.null(inputMatrix) && !is.numeric(as.matrix(inputMatrix))){
#        stop("Expecting numeric input matrix")
#    }
#    
#    if(!require("Matrix")){
#        install.packages("Matrix")
#        require("Matrix")
#    }
#
#    if(!require(SingleCellExperiment)){
#        if (!requireNamespace("BiocManager", quietly = TRUE))
#            install.packages("BiocManager")
#        BiocManager::install("SingleCellExperiment")
#        require(SingleCellExperiment)
#    }
#
#    if(!is.null(dataDir)){
#        geneFile <- paste(dataDir,"features.tsv.gz",sep="/")
#        cellFile <- paste(dataDir,"barcodes.tsv.gz",sep="/")
#        mtxFile <- paste(dataDir,"matrix.mtx.gz",sep="/")
#        ##
#        geneVec <- read.table(gzfile(geneFile),stringsAsFactors=FALSE)[,1]
#        cellVec <- read.table(gzfile(cellFile),stringsAsFactors=FALSE)[,1]
#        exprMat <- Matrix::readMM(gzfile(mtxFile))
#        ##
#        rownames(exprMat) <- geneVec
#        colnames(exprMat) <- cellVec
#        ##
#        sceset <- SingleCellExperiment(assays = list(counts = as.matrix(exprMat)))
#        ## int_metadata() "is not recommended for common users", but I do not see
#        ## where else to add run name       
#        int_metadata(sceset)$runName <- runName
#        ##
#        if(saveRDSfile){
#            if(is.null(rdsOutFile)) {
#                rdsOutFile <- paste0(dataDir,".rds")
#                message("rdsOutFile argument not provided, saving to ",rdsOutFile)
#            }
#            saveRDS(sceset,rdsOutFile)
#        }
#    } else if(!is.null(rdsInFile)){
#        sceset <- readRDS(rdsInFile)
#    } else if(!is.null(inputMatrix)){
#        sceset <- SingleCellExperiment(assays = list(counts = as.matrix(inputMatrix)))
#        int_metadata(sceset)$runName <- runName
#    } else { #better redundant than sorry
#        stop("Need to provide either data directory or RDS file")
#    }
#
#    return(sceset)
####    return(counts(sceset))
#}
