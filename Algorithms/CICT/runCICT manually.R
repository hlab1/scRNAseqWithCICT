#Set working directory and paths
url.base = "/scratch/as15096/eric"
setwd(url.base)

#TODO source all libraries and functions and set up working director
library(doSNOW)
library(doFuture)
# registerDoFuture()
# plan(multiprocess)
#library(parallelly)

n.workers=12
try({ processingCluster <-getMPIcluster()}) #Uses parallelly loaded by doFuture
if(is.null(processingCluster)) processingCluster <-parallelly::makeClusterMPI(n.workers, autoStop = TRUE)

try({
  clusterEvalQ(processingCluster, library(infotheo));
  clusterEvalQ(processingCluster, library(igraph));
  clusterEvalQ(processingCluster, library(minet));
  clusterEvalQ(processingCluster, source('Algorithms/CICT/requirements/rnR_Framework.R'));
})

source('Algorithms/CICT/requirements/0 CICT LibsFunctions.R')
args <- commandArgs(trailingOnly = T)

#Read Config alternative
if(T){

  library(yaml)
  args = read_yaml(paste0(url.base,'/config-files/config_L0.yaml'))
  databases = c('hESC','hHEP','mESC','mHSC-E') # 'dearm5_1','dearm5_3','dearm5_4','DREAM5_1','DREAM5_3','DREAM5_4','hESC',
  for(arg.dname in databases){
    try({
      recalcSimMatrices=F
      #arg.dname <- args$input_settings$datasets[[1]][['name']]
      arg.inFile <- args$input_settings$datasets[[1]][['exprData']]
      arg.pseudoTime <-  args$input_settings$datasets[[1]][['cellData']]
      arg.gtFile <-  args$input_settings$datasets[[1]][['trueEdges']] #ground truth
      arg.edgeTypes <- args$input_settings$algorithms[[1]]$params[['edgeTypes']] # comma separated  'Pearson, ewMImm'
    
      
      # supervised.positiveClass <- args[5] # c:causal edges, 'c,rc': causal and reversecausal
      # supervised.negativeClass<- args[6]  # r:random edges, 'rc': reversecausal  
      # supervised.gtProp<- args[7] #proportion of GT used for classification
      # supervised.train<- args[8] #proportion of classification records used for training
      # fp.fn.cost<- args[9]  #False positive to False Negative cost ratio for threshold identification
    
    
      # Calculate raw edges
      if(T){
        #data/inputs/L1/DREAM5_1/CICT/ExpressionData.csv data/outputs/L1/DREAM5_1/CICT/outFile.txt
        #TODO check if raw edges are not ready prepare it
        url.input = paste0("/scratch/as15096/eric",'/inputs','/L0/',arg.dname,'/' , arg.inFile)
        url.rawedgefile = paste0("/scratch/as15096/eric",'/inputs','/L0/',arg.dname,'/rawEdges.csv')
        
        if(file.exists(url.rawedgefile) & !recalcSimMatrices)
        {
          similarityMatrices = fread(file= url.rawedgefile)
          arg.edgeTypes.existing = colnames(similarityMatrices)
        
        } else {rm(similarityMatrices);arg.edgeTypes.existing = c()}
            
        
        arg.edgeTypes.lst = str_trim(unlist(str_split(arg.edgeTypes,',')))
        arg.edgeTypes.abs = arg.edgeTypes.lst[which(! arg.edgeTypes.lst %in% arg.edgeTypes.existing)]
      
        #produce raw edges
        source('Algorithms/CICT/requirements/calculateRawEdges.R')
        if(length(arg.edgeTypes.abs) > 0 ){
          similarityMatrices.new = calculateRawEdges(arg.edgeTypes.abs,url.input,processingCluster,n.workers  =6)
          if(exists('similarityMatrices') ) 
            {  
            similarityMatrices = 
              similarityMatrices %>% select(-colnames(similarityMatrices.new),src ,trgt) %>%
              merge(similarityMatrices.new,all=T,by=c('src','trgt')) 
            }else similarityMatrices = similarityMatrices.new
        
          write.table(similarityMatrices, url.rawedgefile, sep = ",", quote = FALSE, row.names = FALSE)
          print(sprintf('%s raw edge calculations (%s) finished. %s', 
                        arg.dname, paste0(arg.edgeTypes.lst,collapse = ', ') ,Sys.time()))
        }
       
      }
    })
  }
}
