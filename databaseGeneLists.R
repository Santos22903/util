### Moved to separate file and modified on 11/11/20

### 10/08/20 get gene lists for select data bases
###          return each database as list with pathways as list names,
###          gene sets in that pathway as list entries
dbGeneLists <- function(databases=c("KEGG","Reactome","GOBP")){

   if(!require(msigdbr)){
     install.packages("msigdbr",repos='http://cran.us.r-project.org')
     require(msigdbr)
   }
   if(!require(dplyr)){
     install.packages("dplyr",repos='http://cran.us.r-project.org')
   }
   
   stopifnot(length(databases)>0)
   nameList <- c()
   dfList <- list()

   ### get data frames pathway name - gene name
   if("KEGG" %in% databases){
     sub_kegg <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG") %>%
                 dplyr::select(gs_name, gene_symbol) %>%
                 as.data.frame()
     dfList[["KEGG"]] <- sub_kegg
   }
   if("Reactome" %in% databases){
     sub_react <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") %>%
                  dplyr::select(gs_name, gene_symbol) %>%
                  as.data.frame()
     dfList[["Reactome"]] <- sub_react      
   }
   if("GOBP" %in% databases){
     sub_gobp <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP") %>%
                 dplyr::select(gs_name, gene_symbol) %>%
                 as.data.frame()
     dfList[["GOBP"]] <- sub_gobp
   }
   ### make lists of gene sets for each pathway to run genes2cg with lapply
   makeList <- function(sub_df){
     paths <- unique(sub_df[,1])
     databaseList <- list()
     for(ctr in 1:length(paths)){
        databaseList[[paths[ctr]]] <- sub_df[sub_df[,1]==paths[ctr],2]
     }
     return(databaseList)
   }
   dataLists <- lapply(dfList,makeList)

#   names(dataLists) <- nameList

   return(dataLists)
}
