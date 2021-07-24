#library(mlbench)
library(caret)
#library(Spectrum)
library(Rtsne)
library(umap)
library(devtools)
library(PerformanceAnalytics)
library(sn)

library(tidyr)
library(infotheo)
library(mpmi)
library(stringr)
library(RCy3) #https://github.com/cytoscape/cytoscape-automation/wiki
library(tidyr)
library(infotheo)
library(mpmi)
library(stringr)
library(igraph)
library(DREAM4)
library(data.table)
library(plyr)

library(tidyverse)
library(dplyr)
select<-dplyr::select
mutate<-dplyr::mutate
arrange<-dplyr::arrange
filter<-dplyr::filter
group_by<-dplyr::group_by
require(zeallot)

# data preparation------------
#Benchmarking HESC
if(TRUE){
  # In this challenge we explore the use of Systems Genetics data for elucidating causal network models among genes, 
  # i.e. Gene Networks (DREAM5 SYSGEN A) and predicting complex disease phenotypes (DREAM5 SYSGEN B)
  #Challenge https://www.synapse.org/#!Synapse:syn2787209/wiki/70349
  #DATA: https://static-content.springer.com/esm/art%3A10.1038%2Fnmeth.2016/MediaObjects/41592_2012_BFnmeth2016_MOESM584_ESM.zip
  studyDataset = 'hesc'
  hesc.baseurl = "/Data"
  
  hseq.ptime = fread(paste0(hesc.baseurl,"/BEELINE-data/inputs/scRNA-Seq/hESC/PseudoTime.csv"),header = T);hist(hseq.ptime$PseudoTime)
  recipes::discretize(hseq.ptime$PseudoTime,cuts=5)
  hseq.ptime = hseq.ptime %>%  dplyr::mutate(cellset = infotheo::discretize(hseq.ptime$PseudoTime, "equalfreq", 5)$X);
  hist(hseq.ptime$cellset);table(hseq.ptime$cellset)
  #hseq.cellset1 = hseq.ptime[PseudoTime<=0.02,]
  
  hsec.gtchip = read.csv(paste0(hesc.baseurl,"/BEELINE-Networks/Networks/human/hESC-ChIP-seq-network.csv"),header = T)
  hsec.gtchip$edgetyptruth = 'chipseq'
  colnames(hsec.gtchip) <-c('src','trgt','edgetyptruth') 
  tbl.goldStandard = hsec.gtchip
  
  hesc = read.csv(paste0(hesc.baseurl,"/BEELINE-data/inputs/scRNA-Seq/hESC/ExpressionData.csv"),header = T)
  
  all.dt = hesc %>% rename(gene = X) # %>% select(any_of('vID',hseq.cellset1$V1))
  rownames(all.dt)<- all.dt$gene 
  #hsec1 =  cbind(X = hesc.sbst1$X, exp = rowMeans(hesc.sbst1[,2:ncol(hesc.sbst1)], dims = 1))
  
  length(unique(c(tbl.goldStandard$src,tbl.goldStandard$trgt)))
  
  Nobservations=ncol(all.dt)-1
  ngens = nrow(all.dt)
  
  a=sapply(similarityMatrices[[i]], function(x) sum(is.na(x)));a[a>0]
  
  #Filtering gene set according to the benchmark paper
  geneorder = read.csv(paste0(hesc.baseurl,"/BEELINE-data/inputs/scRNA-Seq/hESC/GeneOrdering.csv"),header = T)
  geneorder.1e3 = geneorder %>% filter(VGAMpValue<0.01) %>% arrange(desc(Variance)) %>% slice_head(n=1000)
  geneorder.5e2 = geneorder %>% filter(VGAMpValue<0.01) %>% arrange(desc(Variance)) %>% slice_head(n=500)
  
  tf.all = unique(tbl.goldStandard$src) # when Chip_seq , gives number of TFs
  #all TFs whose variance had p-value at most 0.01
  tf.set1= geneorder %>% anti_join(geneorder.1e3,by="X") %>% 
    filter(X %in% tf.all &  VGAMpValue>=0.01) %>% arrange(desc(Variance)) 
  
  
  geneorder.1e3 = geneorder %>% filter(X %in% tf.set1$X) %>% rbind(geneorder.1e3  ) %>% unique()
  
  randomgenes = sample_n(setDT(geneorder)[! X %in% geneorder.1e3$X ],size=1000)
  
  geneorder.2e3 = rbind(geneorder.1e3,randomgenes)
  #n.itm.e = n.itm.e %>% filter(src %in% geneorder.1e3$X | trgt %in% geneorder.1e3$X)
  all.dt = all.dt %>% filter(gene %in% geneorder.2e3$X | gene %in% geneorder.2e3$X) 
}


#Mutual information steady state Multiple measures parallel
{
  #Parallel partition for each measure, multiple measures
  if(TRUE){ 
    print(" - Processing in parallel - ")
    #s0m3
    library(doParallel);
    
    # getSimilarityMatrix_MI = function(actualDataset,actualDatasetNNodes,estimators, discretization = TRUE, discretizator){
    #   build.mim(actualDataset, disc = "equalwidth", estimator = "mi.mm")
    # }
    
    # scoreMatrices[[sc]] <- getNormalizedMatrix(scoreMatrices[[sc]],normalization="minimax");
    
    source("import/rnR_Framework.R")
    
    actualDataset <- all.dt
    actualDatasetNNodes <- nrow(actualDataset) + 1;
    actualDatasetNObservations <- ncol(actualDataset);
    actualDatasetName <- "hseq";
    actualDatasetSymbolicPatterns=0;
    actualDatasetPatterns=0
    
    simsGroup=c("MI","CORR","DIST")
    availableGroups <- c("MI","CORR","DIST","SYM");
    availableSimilarities <- vector("list", length(availableGroups));
    names(availableSimilarities) <- availableGroups;
    
    availableSimilarities$MI <- c("efMImm","ewMImm","efMIempirical","ewMIempirical", "efMIshrink","ewMIshrink");
    availableSimilarities$CORR <- c("Pearson","Spearman","Kendall");
    availableSimilarities$DIST <- c("Manhattan","Euclidean","L10Norm");
    #availableSimilarities$SYM <- c("efSym","efSymMI","ewSymMI","efAvgSymMI","ewAvgSymMI","Qual");
    sims <- unlist(lapply(simsGroup, function(group){availableSimilarities[[group]]}), recursive=TRUE, use.names=FALSE);
    
    
    library(doSNOW)
    #processingCluster <-makeCluster(no_cores =5,  outfile = paste("./log/global_prediction_",actualDatasetName,"_.log",sep=""));
    processingCluster<-makeCluster(13,outfile = paste("./log/global_prediction_",actualDatasetName,"_.log",sep=""))
    registerDoSNOW(processingCluster)
    
    clusterEvalQ(processingCluster, library(infotheo));
    clusterEvalQ(processingCluster, library(igraph));
    clusterEvalQ(processingCluster, source("import/rnR_Framework.R"));
    
    similarityMatrices <- parSapply(processingCluster, sims, simplify = FALSE, USE.NAMES = TRUE,
                                    function(sim, actualDataset, actualDatasetNNodes, actualDatasetNObservations, actualDatasetName, actualDatasetSymbolicPatterns, patterns, numCores){
                                      print(paste("[",Sys.time(),"]","Processing",sim,"over",actualDatasetName,"...",sep=" "));
                                      browser()
                                      ##Perform similarity/distance step
                                      
                                      firstStepMatrix <- tryCatch({
                                        
                                        locMatrix <- switch(sim,
                                                            efMImm = getSimilarityMatrix_MI(actualDataset,actualDatasetNNodes,estimators="mi.mm", discretization = TRUE, discretizator = "equalfreq"),
                                                            ewMImm = getSimilarityMatrix_MI(actualDataset,actualDatasetNNodes,estimators="mi.mm", discretization = TRUE, discretizator = "equalwidth"),
                                                            efMIempirical = getSimilarityMatrix_MI(actualDataset,actualDatasetNNodes,estimators="mi.empirical", discretization = TRUE, discretizator = "equalfreq"),
                                                            ewMIempirical = getSimilarityMatrix_MI(actualDataset,actualDatasetNNodes,estimators="mi.empirical", discretization = TRUE, discretizator = "equalwidth"),
                                                            efMIshrink = getSimilarityMatrix_MI(actualDataset,actualDatasetNNodes,estimators="mi.shrink", discretization = TRUE, discretizator = "equalfreq"),
                                                            ewMIshrink = getSimilarityMatrix_MI(actualDataset,actualDatasetNNodes,estimators="mi.shrink", discretization = TRUE, discretizator = "equalwidth"),
                                                            Pearson = getSimilarityMatrix_MI(actualDataset,actualDatasetNNodes,estimators="pearson"),
                                                            Spearman = getSimilarityMatrix_MI(actualDataset,actualDatasetNNodes,estimators="spearman"),
                                                            Kendall = getSimilarityMatrix_MI(actualDataset,actualDatasetNNodes,estimators="kendall"),
                                                            Manhattan = getSimilarityMatrix_DISTANCES(actualDataset,actualDatasetNNodes,norms=1),
                                                            Euclidean = getSimilarityMatrix_DISTANCES(actualDataset,actualDatasetNNodes,norms=2),
                                                            L10Norm = getSimilarityMatrix_DISTANCES(actualDataset,actualDatasetNNodes,norms=10),
                                                            DTWasym = getSimilarityMatrix_DTW(actualDataset,steppatern=asymmetric),
                                                            DTWsym1 = getSimilarityMatrix_DTW(actualDataset,steppatern=symmetric1),
                                                            DTWsym2 = getSimilarityMatrix_DTW(actualDataset,steppatern=symmetric2),
                                                            efSym = getSimilarityMatrix_SYMBOLIC(actualDataset,actualDatasetNNodes,actualDatasetNObservations, simmethod="sym", npatterns=actualDatasetSymbolicPatterns, discretizator = "equalfreq", patterns=patterns, numCores),
                                                            ewSym = getSimilarityMatrix_SYMBOLIC(actualDataset,actualDatasetNNodes,actualDatasetNObservations, simmethod="sym", npatterns=actualDatasetSymbolicPatterns, discretizator = "equalwidth", patterns=patterns, numCores),
                                                            efSymMI = getSimilarityMatrix_SYMBOLIC(actualDataset,actualDatasetNNodes,actualDatasetNObservations,  simmethod="sym.mi", npatterns=actualDatasetSymbolicPatterns, discretizator = "equalfreq", patterns=patterns, numCores),
                                                            ewSymMI = getSimilarityMatrix_SYMBOLIC(actualDataset,actualDatasetNNodes,actualDatasetNObservations,  simmethod="sym.mi", npatterns=actualDatasetSymbolicPatterns, discretizator = "equalwidth", patterns=patterns, numCores),
                                                            efAvgSymMI = getSimilarityMatrix_SYMBOLIC(actualDataset,actualDatasetNNodes,actualDatasetNObservations,  simmethod="avg.sym.mi", npatterns=actualDatasetSymbolicPatterns, discretizator = "equalfreq", patterns=patterns, numCores),
                                                            ewAvgSymMI = getSimilarityMatrix_SYMBOLIC(actualDataset,actualDatasetNNodes,actualDatasetNObservations,  simmethod="avg.sym.mi", npatterns=actualDatasetSymbolicPatterns, discretizator = "equalwidth", patterns=patterns, numCores),
                                                            Qual = getSimilarityMatrix_QUAL(actualDataset,actualDatasetNNodes,actualDatasetNObservations)
                                                            #       efMIml_gc = getSimilarityMatrix_MI(actualDataset,actualDatasetNNodes,estimators="coarse.grained", subestimators="ML",discretizator = "equalfreq"),
                                                            #       ewMIml_gc = getSimilarityMatrix_MI(actualDataset,actualDatasetNNodes,estimators="coarse.grained", subestimators="ML",discretizator = "equalwidth"),
                                                            #       efMImm_gc = getSimilarityMatrix_MI(actualDataset,actualDatasetNNodes,estimators="coarse.grained", subestimators="mm",discretizator = "equalfreq"),
                                                            #       ewMImm_gc = getSimilarityMatrix_MI(actualDataset,actualDatasetNNodes,estimators="coarse.grained", subestimators="mm",discretizator = "equalwidth"),
                                                            #       efMIshrink_gc = getSimilarityMatrix_MI(actualDataset,actualDatasetNNodes,estimators="coarse.grained", subestimators="shrink",discretizator = "equalfreq"),
                                                            #       ewMIshrink_gc = getSimilarityMatrix_MI(actualDataset,actualDatasetNNodes,estimators="coarse.grained", subestimators="shrink",discretizator = "equalwidth"),
                                                            #       Granger = getSimilarityMatrix_MI(actualDataset,actualDatasetNNodes,estimators="granger")
                                        )
                                        
                                        row.names(locMatrix) <- row.names(actualDataset);
                                        colnames(locMatrix) <- row.names(actualDataset);
                                        print(paste("[",Sys.time(),"]","DONE processing",sim,"over",actualDatasetName,"...",sep=" "));
                                        return(locMatrix);
                                        
                                      },error=function(err){
                                        print("Error thrown in thread!")
                                        print(err);
                                        return(NULL);
                                      },warning=function(warn){
                                        print("Warning thrown in thread!")
                                        print(warn);
                                        return(NULL);
                                      });
                                      
                                      
                                      firstStepMatrix
                                      #getNormalizedMatrix(firstStepMatrix,normalization="minimax");
                                      
                                    },
                                    actualDataset, actualDatasetNNodes, actualDatasetNObservations, actualDatasetName, actualDatasetSymbolicPatterns, actualDatasetPatterns, 4); #max(as.numeric((detectCores()-no_cores)/6)-1,1)
    stopCluster(processingCluster)
    
    similarityMatrices=list(res)
    
    n.itm.e = data.frame()
    for(i in 1:length(similarityMatrices)){
      tmp=tmp.1=NULL
      tmp = as.data.frame(similarityMatrices[[i]]) %>% tibble::rownames_to_column()
      tmp.1  = pivot_longer(tmp,cols=colnames(tmp)[2:ncol(tmp)]) 
      names(tmp.1)=c('src','trgt',names(similarityMatrices[i]) )
      if(nrow(n.itm.e)==0) n.itm.e <- tmp.1 else n.itm.e <- inner_join(n.itm.e,tmp.1,on=c("src","trgt"))
    }
    
    #new cols paste0(colnames(n.itm.e),collapse="','")
    #c('efMImm','ewMImm','efMIempirical','ewMIempirical','efMIshrink','ewMIshrink','Pearson','Spearman','Kendall','Manhattan','Euclidean','L10Norm')
    
  }
