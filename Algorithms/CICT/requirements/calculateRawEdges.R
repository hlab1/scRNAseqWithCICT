here::i_am("Algorithms/CICT/requirements/calculateRawEdges.R")

nParallelThreads = 12

#Mutual information steady state Multiple measures parallel
{
  #Parallel partition for each measure, multiple measures
  calculateRawEdges<-function(edgeTypes,url.data,processingCluster,n.workers  =5){ 
    print(" - Processing in parallel - ")
    
    
    source(here::here('Algorithms','CICT','requirements','rnR_Framework.R'))
    #source('/scratch/as15096/eric/Algorithms/CICT/requirements/rnR_Framework.R')
    #s0m3
    library(doParallel);
    #print(url.data)
    outFolder = dirname (url.data)
    actualDataset <-read.csv(url.data) 
    genecol = str_subset(colnames(actualDataset),'X|^(G|g)ene$')
    if(length(genecol)>0) actualDataset =actualDataset %>% column_to_rownames(genecol) 
    actualDataset = actualDataset %>% select_if(is.numeric) #genes in rows and cells in columns  #  stop('Set correct data source') #  all.tdt
      
    
    actualDatasetNNodes <- nrow(actualDataset) + 1;
    actualDatasetNObservations <- ncol(actualDataset);
                                        #actualDatasetName <- basename(url.input)
    actualDatasetName <- basename(url.data)
    actualDatasetSymbolicPatterns=0;
    actualDatasetPatterns=0
    
    simsGroup=c("MI","CORR","DIST")
    availableGroups <- c("MI","CORR","DIST","SYM");
    availableSimilarities <- vector("list", length(availableGroups));
    names(availableSimilarities) <- availableGroups;
    
    availableSimilarities$MI <- c("ewMImm","ewMIempirical","ewMIshrink","efMIempirical", "efMIshrink") #,"efMImm");
    availableSimilarities$CORR <- c("Pearson","Kendall","Spearman") #);
    availableSimilarities$DIST <- c("Euclidean","Granger"); # ,"L10Norm""Manhattan",
    #availableSimilarities$SYM <- c("efSym","efSymMI","ewSymMI","efAvgSymMI","ewAvgSymMI","Qual");
    #sims <- unlist(lapply(simsGroup, function(group){availableSimilarities[[group]]}), recursive=TRUE, use.names=FALSE);
    sims = intersect(edgeTypes,unlist(availableSimilarities))
    
    library(doSNOW)
    #library(parallelly)
    library(doFuture)
    # registerDoFuture()
    # plan(multiprocess)
    
    
    #try({ processingCluster <-getMPIcluster()}) #Uses parallelly loaded by doFuture
    #if(is.null(processingCluster)) processingCluster <-parallelly::makeClusterMPI(n.workers, autoStop = TRUE)
    processingCluster <-parallel::makeCluster(n.workers) #, autoStop = TRUE)
    try({
      clusterEvalQ(processingCluster, library(infotheo));
      clusterEvalQ(processingCluster, library(igraph));
      clusterEvalQ(processingCluster, library(minet));
      clusterEvalQ(processingCluster, here::i_am("Algorithms/CICT/requirements/calculateRawEdges.R"))
      clusterEvalQ(processingCluster, source(here::here('Algorithms','CICT','requirements','rnR_Framework.R')));
      #clusterEvalQ(processingCluster, source('Algorithms/CICT/requirements/rnR_Framework.R'));
    })

    #If has more cores available then the number of simulations, let them be used on parallelizing subprocesses
    n.workers.subprocess = min(n.workers,length(sims))
    
    #similarityMatrices <- sapply(
    similarityMatrices <- parSapply(processingCluster,
      sims, simplify = FALSE, USE.NAMES = TRUE,
      FUN= function(sim, actualDataset, actualDatasetNNodes, actualDatasetNObservations,
                    actualDatasetName, actualDatasetSymbolicPatterns, patterns, numCores)
        {
        
          print(paste("[",Sys.time(),"]","Processing",sim,"over",actualDatasetName,"...",sep=" "));
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
                                #Qual = getSimilarityMatrix_QUAL(actualDataset,actualDatasetNNodes,actualDatasetNObservations)
                                #       efMIml_gc = getSimilarityMatrix_MI(actualDataset,actualDatasetNNodes,estimators="coarse.grained", subestimators="ML",discretizator = "equalfreq"),
                                #       ewMIml_gc = getSimilarityMatrix_MI(actualDataset,actualDatasetNNodes,estimators="coarse.grained", subestimators="ML",discretizator = "equalwidth"),
                                #       efMImm_gc = getSimilarityMatrix_MI(actualDataset,actualDatasetNNodes,estimators="coarse.grained", subestimators="mm",discretizator = "equalfreq"),
                                #       ewMImm_gc = getSimilarityMatrix_MI(actualDataset,actualDatasetNNodes,estimators="coarse.grained", subestimators="mm",discretizator = "equalwidth"),
                                #       efMIshrink_gc = getSimilarityMatrix_MI(actualDataset,actualDatasetNNodes,estimators="coarse.grained", subestimators="shrink",discretizator = "equalfreq"),
                                #       ewMIshrink_gc = getSimilarityMatrix_MI(actualDataset,actualDatasetNNodes,estimators="coarse.grained", subestimators="shrink",discretizator = "equalwidth"),
                                Granger = getSimilarityMatrix_MI(actualDataset,actualDatasetNNodes,estimators="granger")
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
      actualDataset, actualDatasetNNodes, actualDatasetNObservations, actualDatasetName, 
      actualDatasetSymbolicPatterns, actualDatasetPatterns, numCores = n.workers.subprocess
    ); #max(as.numeric((detectCores()-no_cores)/6)-1,1)
       
    #similarityMatrices=list(res)
    n.itm.e = data.frame()
    
    #Merge similarityMatrices into a data.frame and return it 
    for(i in 1:length(similarityMatrices)){
      tmp=tmp.1=NULL
      
      sm = similarityMatrices[i]
      if(is.null(sm)) next
      sm.name = names(sm)
      
      
      tmp = as.data.frame(sm[[1]])
      
      if(!is.data.frame(tmp)) next
      if(nrow(tmp)==0) next
      
      ncol(tmp);nrow(tmp)
      if(all(!  rownames(tmp) %>%as.numeric() %>% is.na())){
          colnames(tmp) = rownames(actualDataset)
          rownames(tmp) = rownames(actualDataset)
          tmp =tmp %>% mutate(src =rownames(actualDataset)  ) %>% select(src,everything())
      } else tmp = tmp %>% tibble::rownames_to_column()
      tmp.1  = pivot_longer(tmp,cols=colnames(tmp)[2:ncol(tmp)],names_to = 'trgt')  
      
      names(tmp.1)=c('src','trgt',sm.name)
      tmp.1= tmp.1 %>% dplyr::filter(!is.na(src) & !is.na(trgt))
      
      
      if(nrow(n.itm.e)==0) 
        n.itm.e <- tmp.1 else 
          n.itm.e <- merge(tmp.1,n.itm.e,all.y=T, by=c("src"='src',"trgt"='trgt'))
    }

    parallel::stopCluster(processingCluster)
    n.itm.e
    #new cols paste0(colnames(n.itm.e),collapse="','")
    #c('efMImm','ewMImm','efMIempirical','ewMIempirical','efMIshrink','ewMIshrink','Pearson','Spearman','Kendall','Manhattan','Euclidean','L10Norm')
    
  }
  
}
