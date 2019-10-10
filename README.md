# Gene-Cluster-Shinny-App
Gene Cluster Shiny Application using RnaSeq files and Go.Db Kegg.Db and Biomart characteristics

This is a shiny Application that gives you the opportunity to use an **RnaSeq file** , enrich this dataset with specific characteristcs of those genes using **Go.Db**, **Kegg.Db** and **Biomart** , in order to **cluster** this with three different algorithms (**kmodes,kproto,kmeans** for every different)

The basic steps that this application has are:
1. The user inserts a RnaSeq file with row names the ensembl ids of genes and as data integer and normalized values choosing the right separation
2. After that he has to choose one or more variables he wants to add in his data. These attributes areQ: 
   * Terms from Go.db
   * Orthology,Pathway and Motif from Kegg.Db
   * Band,Gene biotype and Transript count from Biomart
3. Wait until the Upload and the execution of the algorithms
   After that there will be three different data sets, one with arithmetical data,one with categorical and one with mixed type data 
   These will be used for kmeans,kmodes and kproto respectively
4. Then choose the algorithm and the parameters he wants in order to show the graphs and the visualization of the results
5. Finally the visualization must be refered, and they are:
   * Arithmetical analysis of every cluster with boxplots
   * Categorical analysis of every cluster
   * Upset diagram for genes with common gcharacteristics
   * Sankey diagram for comon genes between the clusters
   * Graphs showing the evaluation metrics of every clustering
   
The packages you will need are in a command( you need the latest Rstudio version):
```
load.lib<-c("BiocManager", "shiny", "networkD3", "shinybusy",  "UpSetR", "dplyr",  "DT", "klaR" ,"clustMixType" ,"cluster" ,"readr" ,"Rtsne")

install.lib<-load.lib[!load.lib %in% installed.packages()]

for(lib in install.lib) install.packages(lib,dependencies=TRUE)




if (!requireNamespace("BiocManager", quietly = TRUE))

install.packages("BiocManager")
    
BiocManager::install("Organism.dplyr")

BiocManager::install("GenomicRanges")

BiocManager::install("KEGGREST")

BiocManager::install("AnnotationDbi")

BiocManager::install("GO.db")

BiocManager::install("hugene20sttranscriptcluster.db")

BiocManager::install("biomaRt")

BiocManager::install("org.Hs.eg.db")

BiocManager::install("DESeq2")

BiocManager::install("BiocGenerics")

BiocManager::install("Biobase")


library(shiny)

library(networkD3)

library(shinybusy)

library(AnnotationDbi)

library(org.Hs.eg.db)

library(hugene20sttranscriptcluster.db)
library(GO.db)

library(KEGGREST)

library(UpSetR)

library(cluster)

library(dplyr)

library(readr)

library(Rtsne)

library(DESeq2)


library(BiocGenerics)
library(biomaRt)

library(Biobase)

library(DT)
```
