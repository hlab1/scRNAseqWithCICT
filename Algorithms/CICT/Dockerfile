FROM r-base:4.1.1

LABEL maintainer = "Abbas Shojaee <abbas.shojaee@gmail.com>"

RUN apt-get update && apt-get install time #Installl time command to compute time taken

USER root

WORKDIR /

RUN R -e "install.packages(c('BiocManager','GenNet' ,'myTAI','DREAM4','Seurat','RCytoscape','phateR',''mpmi','BiocManager','caret','Spectrum','Rtsne','umap','devtools','PerformanceAnalytics','sn','tidyr','infotheo','mpmi','stringr','RCy3' ,'tidyr','infotheo','mpmi','stringr','igraph','DREAM4','data.table','plyr','tidyverse','dplyr','zeallot')"

RUN R -e "BiocManager::install('lattice');  BiocManager::install('GenomicRanges');  BiocManager::install('SummarizedExperiment');  install_github('cran/PCIT');  BiocManager::install('netbenchmark');  BiocManager::install('grndata');  BiocManager::install('DREAM4');  BiocManager::install('RCytoscape');  BiocManager::install('RCy3');  BiocManager::install('graph');"
  
COPY runCICT.R /
COPY CICT.R /

RUN mkdir data/

RUN apt-get update && apt-get install time