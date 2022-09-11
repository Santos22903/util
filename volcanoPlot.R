### 05/20/19 Function making volcano plot
###          Adopted/modified from Byzova addit_analysis.R

### 02/01/20 Added new arguments: low/high FC cuts using
###          values of logs, parameter to label lowest in q and
###          highest in logFC genes; fileOut name as argument

###          So, added adding hiq/lowFC to genes 2 label 
###	     Also, added showing in blue/red # of genes above/below
###            logFC cuts
### 	     Changed label coloring: blue/red for above/below FC
###          Changed output to pdf - better quality

### Arguments:
###           geneDF - data frame with genes as row names, containing
###                    columns log2FoldChange and padj
###           interactive - to print to screen as opposed to to file
###           localFCcutoff - draw "cutoff line" at this value of FC
###           qcutoff - draw "cutoff line" at this value of padj
###           genes2highlight - show this genes in red and print names
###	      fileOut - optional parameter, name for volcano plot file
###           low_q_out - to send/not list of marked low q genes to file
###	      high_fc_out - to send/not list of marked high FC genes to file

makeVolcanoPlot <- function(geneDF,interactive=TRUE,qcut=0.05,
			    localFCcut_hi=1.5,localFCcut_lo=2^(-log2(localFCcut_hi)),
			    Nlo_q_tolabel=NULL, Nhi_FC_tolabel=NULL,
			    genes2highlight=c(""),fileOut=NULL,
			    low_q_out=T,high_fc_out=T){

    if(is.null(fileOut)){fileOut <- paste0("volcanoPlot",Sys.Date(),".pdf")}

    ## removing NAs
    geneDF <- geneDF[complete.cases(geneDF),]

    ## if not null, add low q/high FC genes to genes to label
    ## low q
    if(!is.null(Nlo_q_tolabel)){
        geneDF <- geneDF[order(geneDF$padj),]
    	addedGenes <- rownames(geneDF)[1:Nlo_q_tolabel]
	if(low_q_out){
	   sink(paste0(fileOut,"_low_q_genes.txt"))
	   cat("Showing lowest padj genes:\n")
	   print(addedGenes)
	   sink()
	}
    	genes2highlight <- c(genes2highlight,addedGenes)
    }
    ## hi logFC
    if(!is.null(Nhi_FC_tolabel)){
        geneDF <- geneDF[order(abs(geneDF$log2FoldChange),decreasing=T),]
    	addedGenes <- rownames(geneDF)[1:Nhi_FC_tolabel]
	if(high_fc_out){
	   sink(paste0(fileOut,"_high_fc_genes.txt"))
	   cat("Showing highest FC genes:\n")	
	   print(addedGenes)
	   sink()
	}
    	genes2highlight <- c(genes2highlight,addedGenes)
    }
    ## numbers of high/low (beyond +/-localFCcut)
    message("localFCcut_hi ",localFCcut_hi)
    nHigh <- length(which(geneDF$log2FoldChange>log2(localFCcut_hi)))
    message("localFCcut_lo ",localFCcut_lo)
    nLow <- length(which(geneDF$log2FoldChange<log2(localFCcut_lo)))
    
    ## add column with label to show on the plot
    geneDF$label2show <- ifelse(rownames(geneDF) %in% genes2highlight,
                                   rownames(geneDF),"")
    labelColor <- ifelse(geneDF$log2FoldChange>log2(localFCcut_hi),"red",
                  ifelse(geneDF$log2FoldChange<log2(localFCcut_lo),"blue","black"))

    if(!require("ggplot2")){
        install.packages("ggplot2")
        require(ggplot2)
    }
    if(!require("ggrepel")){
        install.packages("ggrepel")
        require(ggrepel)
    }

    absXmax <- max(max(geneDF$log2FoldChange),  # X position for hi/lo labels,
                   abs(min(geneDF$log2FoldChange))) # will also make plot symmetric
    txtYpos <- max(-log10(geneDF$padj[geneDF$padj>0]))

    volc = ggplot(geneDF, aes(log2FoldChange, -log10(padj))) +
        xlim(-1.25*absXmax,1.25*absXmax) +
    	geom_vline(aes(xintercept=log2(localFCcut_lo)),color="black") +
        annotate("text", x=log2(localFCcut_lo)-0.1, y=txtYpos, hjust=1,vjust=0,
                 label=paste0("FC = ",round(localFCcut_lo,2)),colour="black", angle=90) +
        geom_vline(aes(xintercept=log2(localFCcut_hi)),color="black")  +
        annotate("text", x=log2(localFCcut_hi)-0.1, y=txtYpos, hjust=1,vjust=0,
                 label=paste0("FC = ",round(localFCcut_hi,2)),colour="black", angle=90) +
        geom_hline(aes(yintercept=-log10(qcut)),color="black") +
        annotate("text", x=absXmax, y=-log10(qcut), hjust=0.5, vjust=1,
                 label=paste0("q value of ",qcut),colour="black") +
        geom_point(col=labelColor) +
	theme(plot.margin = unit(c(1,3,1,1), "lines")) +  # Make room for the grob
    	annotate("text", x=-1.1*absXmax,y=txtYpos,hjust=1,label=nLow,color="blue",size=5) + 
        annotate("text", x=1.1*absXmax,y=txtYpos,hjust=0,label=nHigh,color="red",size=5)  +
        geom_text_repel(data=geneDF, aes(label=label2show), color=labelColor)

    if(!interactive){
        pdf(fileOut)
    }
    print(volc)
    if(!interactive){
        dev.off()
    }

    return(0)
}
