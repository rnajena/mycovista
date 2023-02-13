
######################################################################################################################
## R source code to read GWAS results from PPanGGOLin and analyse gene-wise associations of genes with the clusters
## Requires PPanGGOLin gene presence absence file (gene_presence_absence.Rtab) and 
## a csv file assigning cluster number to sample IDs (cluster_anno_numeric.csv)

library(dplyr)
library(tidyverse)
library(writexl)
## Read the results from PPanGGOLiN analysis
presData <- read.delim("gene_presence_absence.Rtab")
colnames(presData) <- gsub("^X","",colnames(presData))

## Read a cluster table file (columns ID and cluster: specifies integer cluster number)
clusterAnno <- read.delim("cluster_anno_numeric.csv")
clustAnnoList <- clusterAnno %>% group_split(cluster)

resList <- list()
for(i in 1:nrow(presData)){
  
  idxA <- clustAnnoList[[1]]$ID
  idxA <- idxA[idxA %in% colnames(presData)]
  idxB <- clustAnnoList[[2]]$ID
  idxB <- idxB[idxB %in% colnames(presData)]
  idxC <- clustAnnoList[[3]]$ID
  idxC <- idxC[idxC %in% colnames(presData)]
  clustA <- round(sum(presData[i,idxA]) / length(idxA),2)
  clustB <- round(sum(presData[i,idxB]) / length(idxB),2)
  clustC <- round(sum(presData[i,idxC]) / length(idxC),2)
  
  resList[[i]] <- data.frame(clustA,clustB,clustC)
}
resData <- do.call("rbind",resList)
resData$Gene <- as.character(presData$Gene)


## Statistical significance with Fisher exact test ##########################################################################################
#############################################################################################################################################

## Function runFisherTest
## Run Fisher Exact Test for two set of IDs
## Input:
## presData : data.frame from gene_presence_absence.Rtab
## setOne : set of sample IDs cluster 1
## setTwo : set of sample IDs cluster 2

runFisherTest <- function(pData,setOne,setTwo){
  
  resultList <- lapply(pData$Gene,function(x){
    
    geneRow <- pData %>% filter(Gene == x)
    ## Get the cluster A and B genes
    geneRowA <- geneRow[, colnames(geneRow) %in% setOne]
    geneRowB <- geneRow[, colnames(geneRow) %in% setTwo]
    
    contMatrix <- matrix(0,ncol=2,nrow=2,dimnames = list(c("Present", "Absent"),c("ClustA", "ClustB")))
    contMatrix["Present","ClustA"] <- sum(geneRowA == 1)
    contMatrix["Present","ClustB"] <- sum(geneRowB == 1)
    contMatrix["Absent","ClustA"] <- sum(geneRowA == 0)
    contMatrix["Absent","ClustB"] <- sum(geneRowB == 0)
    
    resDF <- data.frame(Gene=x,Present_SetOne=contMatrix["Present","ClustA"],Absent_SetOne=contMatrix["Absent","ClustA"],Present_SetTwo=contMatrix["Present","ClustB"],
                        Absent_SetTwo=contMatrix["Absent","ClustB"],Pval=fisher.test(contMatrix,alternative="two.sided")$p.value,adj.Pval=1)
    return(resDF)
  })
  
  resultTable <- do.call("rbind",resultList)
  resultTable[,"Pval"] <- signif(resultTable[,"Pval"],3)
  
  adjpval <- p.adjust(resultTable[,"Pval"],method="BH")
  resultTable[,"adj.Pval"] <- adjpval
  
  resultTable <- resultTable %>% filter(adj.Pval < 0.01)
  
  resultTable <- resultTable %>% arrange(adj.Pval)
  
  resultTable$Pval <- formatC(resultTable$Pval,format="e",digits=2)
  resultTable$adj.Pval <- formatC(resultTable$adj.Pval,format="e",digits=2)
  
  ## Addtional filter max 2 columns larger 2
  resultTable <- resultTable %>% mutate(PresentCalls=(Present_SetOne > 4) + (Absent_SetOne > 4) + (Present_SetTwo > 4) + (Absent_SetTwo > 4))
  
  ## Additional filter
  resultTable <- resultTable %>% filter(PresentCalls <= 2)
  
  return(resultTable)
}



setA <- clustAnnoList[[1]]$ID
setB <- clustAnnoList[[2]]$ID
setC <- clustAnnoList[[3]]$ID


setOne <- c(setA)
setTwo <- c(setB,setC)
resA <- runFisherTest(setOne,setTwo)

setOne <- c(setB)
setTwo <- c(setA,setC)  
resB <- runFisherTest(setOne,setTwo)

setOne <- c(setC)
setTwo <- c(setA,setB)
resC <- runFisherTest(setOne,setTwo)

write_xlsx(list(ClusterA=resA,ClusterB=resB,ClusterC=resC),"GWAS_Fisher_Table.xlsx")


