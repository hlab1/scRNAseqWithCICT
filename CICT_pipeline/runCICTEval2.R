#Set working directory and paths



print('step 000')

url.base = "/scratch/as15096/eric"

setwd(url.base)

#source all libraries and functions and set up working director
source('/scratch/as15096/eric/Algorithms/CICT/requirements/CICT_LibsFunctions.R')

args.cmnd <- commandArgs(trailingOnly = T)
#args.cmnd = c('runCICT_par','/scratch/as15096/eric/outputs/cict_par/sens_sparsity/parConf_7.yaml','FALSE') #sens_edgeType runCICT_par  #runCICT_par
#args.cmnd = c('runCICT','config_L2_ns.yaml','TRUE') #runCICT_par  #runCICT_par  config_SERGIO_DS4.yaml
#args.cmnd = c('calcEdges','config_L2.yaml','TRUE') 
#args.cmnd = c('runCICT','config_SERGIO_DS7.yaml','TRUE') #runCICT_par  #runCICT_par  config_SERGIO_DS4.yaml
#args.cmnd = c('calcEdges','mESC-scalingtmp.yaml','TRUE') #runCICT_par  #runCICT_par  config_SERGIO_DS4.yaml
#args.cmnd = c('config_par','mESC-scalingtmp.yaml','TRUE') #runCICT_par  #runCICT_par  config_SERGIO_DS4.yaml
#args.cmnd = c('runCICT_par','/scratch/as15096/eric/outputs/cict_par/cict_scaling/parConf_6.yaml','TRUE') #runCICT_par  #runCICT_par  config_SERGIO_DS4.yaml
#args.cmnd = c('runSupervised','config_L2_lofgof.yaml','TRUE') #runCICT_par  #runCICT_par  config_SERGIO_DS4.yaml
#args.cmnd = c('runCICT','config_L2_lofgof.yaml','TRUE') # For multisets mESC_lofgof


print('step 001')

if(length(args.cmnd)<2 ){
  print('Please provide [calcEdges|runCICT|runCICT_par|calcEdges_par|runSupervised|] 
                          [config file(s) path or config file name]  
                          [ignore and overwrite existing outputs FALSE|TRUE]') 
  return(-1)
} 

operation<-args.cmnd[1]
configFilePath <- args.cmnd[2]
forceOutput<- as.logical(args.cmnd[3])

print('step 002')
if(is.na(forceOutput)) forceOutput = FALSE
#arg.experiment <- args.cmnd[4]

if(!operation %in% c('calcEdges','runCICT','runSupervised','runCICT_par','install','calcEdges_par')){
  print(paste0('Invalid operation requested: ', operation))
  return(-1)
} else print(sprintf("============================== Operation: %s, Config: %s, forceOutputReplace: %s", 
                     operation, configFilePath, forceOutput))


#CICT caller definition
{
  CICT<-function(theJobID,
                 url.input,
                 url.rawedgefile,
                 url.name.map,
                 url.gt,
                 url.output,
                 url.logfile='',
                 cictRawEdgeCol,
                 earlyThresholdForGraphAnalysis,
                 minGroundTruth.ratio.learning,
                 maxunseenTest.ratio,
                 maxGroundTruth,
                 randomEdgesFoldCausal,
                 sampling.c_rc.ratio,
                 trainingTarget='class2',
                 tstPrecent = .3,
                 forceOutput=F,
                 arg.experiment = NA,
                 FLAG_runOnAllEdges=T,
                 FLAG_exportRankedEdges = T,
                 FLAG_exportTrainAndTest=F,
                 Debug=F,
                 RF_max_depth =20,RF_ntrees=50,
                 preset.train=NA, preset.test=NA,
                 maxNetSize=NA,
                 ...)
  {
    
    name.map=n.itm.e=all.dt=tbl.goldStandard=NULL
    arguments <- list(...)
    #Define urls/paths
    rcrd=list()
    rcrd$outputAlreadyPresent=F
    rcrd$trainingTrgt = trainingTarget
    
    if(is.null(trainingTarget)) stop("Training Target is null")
    
    if(is.null(earlyThresholdForGraphAnalysis))earlyThresholdForGraphAnalysis=0
    #Adds a subfolder to CICT for outputs of an experiemntal run, e.g. different params
    if(!is.na(arg.experiment)) url.output=paste0(url.output,'/',arg.experiment) 
    
    if(!dir.exists(url.output)) dir.create(url.output,recursive = T)
    if(url.logfile=='') url.logfile = paste0(url.output , '/',arg.dname,'_cict log.txt')
    write('Start',  file = url.logfile, append = F)
    
    
    #reporting
    {      
      msg = c("Operation: Running CICT ============================ \n",
              sprintf("Data= %s | operation= %s | experiment= %s | ground truth= %s \n",
                      arg.dname,operation,arg.experiment,arg.gtFile ))
      cat(msg)
      write(msg,file = url.logfile, append = TRUE)
      
      print('step 01')
    } 
    
    #Read rawedge file and goldstandard
    if(file.exists(url.rawedgefile) )
    {
      n.itm.e = fread(url.rawedgefile) #%>% mutate(src=toupper(src), trgt = toupper(trgt) ) 
      if(!any(colnames(n.itm.e)== cictRawEdgeCol) )
        stop("Requested raw edge measure does not exist in rawEdgeFile")
      
      if(file.exists(url.name.map)) {
        try({name.map = fread(url.name.map) })
        if(!is.null(name.map))
          suppressMessages({
            n.itm.e=n.itm.e %>% mutate(
              src = mapvalues(src,name.map$`#ID`,name.map$Name),
              trgt = mapvalues(trgt,name.map$`#ID`,name.map$Name))
          })
      }
      all.dt = fread(url.input)  
      colnames(all.dt)[1] <-'gene'
      Nobservations=ncol(all.dt)-1
      all.dt = all.dt %>%  dplyr::mutate(gene = toupper(gene))
      
      
      tbl.goldStandard = fread(url.gt,fill=TRUE )  
      print('step 03')
      
      
      try({tbl.goldStandard=tbl.goldStandard %>%  
        mutate(edgetruth=ifelse('Type' %in% names(.),Type,1)) })
      
      try({    
        tbl.goldStandard=tbl.goldStandard%>% 
          rename(src = Gene1,trgt=Gene2) %>% select(-any_of('V3')) %>%
          dplyr::filter(edgetruth==1 | edgetruth== '+') %>%
          mutate(edgetyptruth='chipseq')
      })
      intersect(n.itm.e$src , tbl.goldStandard$src);intersect(n.itm.e$trgt , tbl.goldStandard$trgt)
      
      print(paste0('Number of edges from n.itm.e in goldstandard = ' ,
                   nrow(n.itm.e %>% inner_join(tbl.goldStandard, by = c('src'='src','trgt'='trgt')))
      )
      )
      
      
    } else {
      rm(similarityMatrices);
      arg.edgeTypes.existing = c()
      return(list())
    }
    
    if(Debug) browser()
    
    if(!is.na(maxNetSize)) 
      if(maxNetSize>0){
        
        #Making sure that the final sample has some ground truth
        genes.gtsmpl = c(tbl.goldStandard$src) %>% head(size = maxNetSize/10) %>% na.omit()  %>% unique() 
        genes.sampled= c(n.itm.e$src,n.itm.e$trgt) %>% unique() %>% 
          setdiff(genes.gtsmpl) %>% sample(size = maxNetSize - length(genes.gtsmpl) ) 
        genes.included = c(genes.gtsmpl,genes.sampled)
        n.itm.e = n.itm.e %>% filter(src %in% genes.included & trgt %in% genes.included)
      }
    
    
    if( (!file.exists(url.rankedEdges) | !file.exists(url.rankedEdgesGated)) | forceOutput) {
      #try({   
      
      source('Algorithms/CICT/CICT.R',local=TRUE)
      
      c(rcrd,edges,vertices) %<-z% prepareEdgeFeatures(Debug=Debug)
      rcrd = cictTrainTestReport(edges,rcrd,Debug=Debug,preset.train=preset.train, preset.test=preset.test)
      
      return(rcrd)
      #})
    } else {
      print(paste0('Output files already present. Passing on ',arg.dname))
      rcrd$outputAlreadyPresent=T
      return(rcrd)
    }
    
  }
  
} 
#Parallel run


if(operation=='config_par'){
  
  #https://future.batchtools.futureverse.org/
  #https://mllg.github.io/batchtools/articles/batchtools.html
  #https://github.com/mllg/batchtools/blob/master/inst/templates/slurm-simple.tmpl
  
  #“apply algorithm A to problem instance P and store results”
  library(future.batchtools)
  library(batchtools)
  library("future")
  library("listenv")
  library(data.table)
  options(knitr.table.format = 'pipe')
  
  
  runs=1:1
  size.groundTruth=500
  trainingTarget ='class2' # c('class2','class3')
  #Add problems and algorithms 
  {
    
    records = list()
    
    ########### GENERAL settings for config generation
    
    #configFiles = list.files(paste0(url.base,'/config-files/'), pattern = 'SERGIO')
    configFiles=list.files(paste0(url.base,'/config-files/'), pattern = 'mESC-scalingtmp' ) #'L\\d{1}[.]') #"config_L2.yaml"
    #configFiles = 'config_SERGIO_DS4.yaml' ;cnfg=configFiles
    
    #datasetnames = NULL, # c('mDC','mHSC-E','mESC','hESC','hHEP'),
    #edgeTypes = c('Pearson','ewMImm' ),#'Euclidean','ewMIempirical',

    
    theExperimentName = "cict_scaling" # "sens_sparsity" #"sens_multipleRuns"  #"sens_edgeType" #"sens_modeling_choices"   #  'calcEdgesConfs'  # 
    jobConfigBase.url = paste0(url.base,'/outputs/cict_par/',theExperimentName,'/')
    if(!dir.exists(jobConfigBase.url)) dir.create(jobConfigBase.url,recursive = T)
    addToExperiment=F  #Add new config files to experiment directory or start anew
    
    if(addToExperiment == F){
      if(readline(prompt = paste0("Do you want to delete existing files in ", theExperimentName, " folder? (y/n)")) == 'y')
        file.remove(paste0(jobConfigBase.url,'/', list.files(jobConfigBase.url,pattern = '*.yaml',all.files=F)))
      jobConfigBase.idx=0
    }else{
      jobConfigBase.idx = 
        paste0(jobConfigBase.url) %>% 
        list.files( pattern = '*.yaml') %>%
        str_extract_all("(?<=_)\\d{1,4}(?=[.]yaml)") %>%
        unlist() %>% as.numeric() %>% max()
    }
    
    
    #Create combinations of settings ***************************
    rm(prd1,prd2,prd3)
    prd2 = data.frame()
    for(cnfg in configFiles ){ #cnfg = configFiles[1]
      try({
        rm(prd1)
        print(paste0("=====> ", cnfg))  
        cnfg.url = paste0(url.base,'/config-files/',cnfg)
        args = read_yaml(cnfg.url) #
        
        
        expLevel = args$input_settings$dataset_dir
        datasetnames = rbindlist(args$input_settings$datasets) %>% pull('name')
        print(paste0(datasetnames,collapse=','))
        
        #**********************************************************
        #Temporary, removing dream data from datasets
        datasetnames = str_subset(datasetnames,'dream',negate = T) 
        #***********************************************************
        
        #****** produce parCalcEdges configs *************************************
        if(F){
          runs = 1
          
          cictRawEdgeCol ='Pearson'
          
          prd1 = purrr::cross(list(runs,cictRawEdgeCol,datasetnames))
          prd1 = rbindlist(prd1)
          colnames(prd1)<- c('runs','cictRawEdgeCol','dataset')
          
          prd1$size.groundTruth=250
          prd1$training.targets = 'class2'
          
        }
        
        #****** Sensitivity Analysis Sparsity *************************************
        if(F){
          runs = 1
          RF_max_depthS = 5 #c(5,10)
          size.groundTruthS= 500 #c(250,500)
          #RF_ntrees = c(20,50)
          
          cictRawEdgeCol =args$input_settings$algorithms[[1]]$params[['edgeTypes']] # comma separated  'Pearson, ewMImm'
          cictRawEdgeCol  = str_split(cictRawEdgeCol,',',simplify = T)  %>% trim() %>% str_subset(pattern=".{2}")
          cictRawEdgeCol='ewMIshrink' 
          
          prd1 = purrr::cross(list(runs,cictRawEdgeCol,datasetnames,RF_max_depthS,size.groundTruthS))
          prd1 = rbindlist(prd1)
          colnames(prd1)<- c('runs','cictRawEdgeCol','dataset','RF_max_depth','size.groundTruth')
          
          prd1$RF_ntrees = 50
          prd1$trainingTarget = 'class2'
          
        }
        
        #****** Sensitivity Analysis Multiple runs   *************************************
        if(F){
          runs=1:100
          
          prd1 = purrr::cross(list(runs,expLevel,datasetnames))
          prd1 = rbindlist(prd1)
          colnames(prd1)<- c('runs','expLevel','dataset')
          prd1$cictRawEdgeCol ='Pearson' # 'ewMIshrink' # 'ewMIshrink'
          prd1$size.groundTruth=500
          prd1$RF_ntrees = 50
          prd1$RF_max_depth = 10
          prd1$trainingTarget = 'class2'
          
        }
        
        #****** Sensitivity Analysis Modeling Choices  *************************************
        if(F){
          runs = 1:5
          size.groundTruth=c(250,500,1000,2000,4000,8000)
          training.targets = c('class2','class3')
          RF.max_depth = c(5,10,15,20)
          RF.ntrees = c(5,10,20,50)
          
          
          prd1 = purrr::cross(list(runs,expLevel,datasetnames,size.groundTruth,training.targets,RF.max_depth,RF.ntrees))
          prd1 = rbindlist(prd1)
          colnames(prd1)<- c('runs','expLevel','dataset','size.groundTruth','trainingTarget','RF_max_depth','RF_ntrees')
          
          prd1$cictRawEdgeCol = 'Pearson' #'ewMIshrink' # 'ewMIshrink'
        }
        
        #****** CICT Scaling  *************************************
        if(T){
          runs = 1
          maxNetSize=c(125,250,500,1000,2000,4000,8000)
          datasetnames = 'mESC'
          training.targets = c('class2')
          RF.max_depth = c(20)
          RF.ntrees = c(10)
          
          
          prd1 = purrr::cross(list(runs,expLevel,datasetnames,size.groundTruth,training.targets,RF.max_depth,RF.ntrees,maxNetSize))
          prd1 = rbindlist(prd1)
          colnames(prd1)<- c('runs','expLevel','dataset','size.groundTruth','trainingTarget','RF_max_depth','RF_ntrees','maxNetSize')
          
          prd1$cictRawEdgeCol = 'Pearson' #'ewMIshrink' # 'ewMIshrink'
        }
        
        prd1$datasetdir = expLevel
        
        prd1 = prd1 %>% mutate(
          dataFolder= args$input_settings$dataset_dir[[1]],
          configFile = cnfg.url
        ) %>% mutate(
          jobDesc = paste0('gt_',size.groundTruth,'_', datasetdir,'_', dataset,'_',cictRawEdgeCol,'_',runs)
        )
        
        
        prd2 = rbind(prd2,prd1)
      })
      #addProblem(name = 'mDC', data = url.inputFolder , fun = subsample, seed = 42) #Feed config file to data
      #problemDeisgns = append(problemDeisgns,cnfg)
    }

    
    #Scanning config files *******************************
    library("parallel")
    library("foreach")
    library("doParallel")
    
    parcnfTemplate = paste0(url.base,'/CICT_parallel_config_template.yaml') %>% read_yaml()
    if(!dir.exists(jobConfigBase.url)) dir.create(jobConfigBase.url,recursive = T)
    
    numCores=14
    cl<-makeCluster(numCores, outfile="")
    registerDoParallel(cl)
    print("Cluster registered...")
    
    setDT(prd2)
    prd3=foreach(idx = jobConfigBase.idx:(jobConfigBase.idx+nrow(prd2)),.combine = rbind
    ) %dopar%{
      print(idx)
      p=prd2[idx,]
      tmplt = parcnfTemplate
      tmplt$experiment = theExperimentName
      tmplt$given_job_id = idx 
      tmplt$given_job_desc = p$jobDesc 
      
      tmplt$expLevel <- p$expLevel
      tmplt$dataset_dir = p$dataFolder        
      tmplt$datasets[[1]]$name = p$dataset
      tmplt$CICT$maxNetSize = p$maxNetSize
      
      tmplt$CICT$edgeType= p$cictRawEdgeCol
      tmplt$CICT$earlyThresholdForGraphAnalysis =0
      tmplt$CICT$supervised.positiveClass = c('c','rc')# c:causal edges, 'c,rc': causal and reversecausal
      tmplt$CICT$supervised.negativeClass= c('r')  # r:random edges, 'rc': reversecausal
      tmplt$CICT$minGroundTruth.ratio.learning = .2
      tmplt$CICT$maxGroundTruth  = p$size.groundTruth  #controls the total size of learning set
      tmplt$CICT$randomEdgesFoldCausal = 5
      tmplt$CICT$sampling.c_rc.ratio = .3 
      tmplt$CICT$tstPrecent = .3
      tmplt$CICT$trainingTarget = p$trainingTarget
      tmplt$CICT$tstPrecent = .3
      tmplt$CICT$maxunseenTest.ratio = .1
      
      tmplt$CICT$RF_max_depth = ifelse(is.null(p$RF_max_depth),10,p$RF_max_depth) #20 is the default
      tmplt$CICT$RF_ntrees = ifelse(is.null(p$RF_ntrees),50,p$RF_ntrees) #50 is the default
      
      tmplt$originalconfigFile = p$configFile
      
      outputUrl = paste0('parConf','_',idx)
      tmplt$parConf = outputUrl
      #outputUrl = paste0(idx,'_',tmplt$given_job_desc)
      tmplt$output_record = paste0(outputUrl,'_output.txt')
      
      library(data.table)
      prd2=prd2[idx,`:=`(given_job_id = idx,parConf = outputUrl, 
                         rslt = paste0('rslt_',idx,'.rds'))]
      
      
      #write(toString(tmplt),file =paste0(jobConfigBase.url,outputUrl,'.yaml')) 
      yaml::write_yaml(tmplt,file =paste0(jobConfigBase.url,outputUrl,'.yaml') )
      prd2[idx,]
    }
    
    stopCluster(cl)
    
    
    #saveRDS(results.all,file =url.cluster.results )
    #a=results.all %>% kable(format='html')
    #save(a,file=url.cluster.results)
    #fwrite(results.all,url.cluster.results)
    
    
    
    # library(tableHTML)
    # tableHTML(results.all)
    # #browseURL(rmarkdown::render(input = "view_template.Rmd", params = list(myinput = iris)))
    # 
    # 
    # tab = ijoin(pars, results)
    # head(tab)
    # 
    # unwrap(getJobResources())
    # 
    # tab[ratio == 0.67, list(mmce = mean(mce)),
    #     by = c("algorithm", "kernel", "epsilon", "ntree")]
    
    
  }
}else
  if(operation=='calcEdges_par'){
    
    #Only takes two arguments one for config file and one for operation
    library(yaml)
    #configFilePath = "/scratch/as15096/eric/outputs/cict_par/calcEdgesConfs/parConf_3.yaml"
    rm(cnf)
    cnf = read_yaml(configFilePath) 
    
    arg.dname <- cnf$datasets[[1]]$name  #
    #cictRawEdgeCol<-cnf$CICT$edgeType
    arg.edgeTypes= cnf$calcEdges$edgeTypes
    #Reading other information that does not affect scaling, from the config file, 
    # allows less complex code and customizability without code modification
    arg.inputsdir <-  cnf$input_dir
    arg.dataFolder <- cnf$dataset_dir
    
    
    arg.inFile <- cnf$datasets[[1]]$exprData
    arg.pseudoTime <-  cnf$datasets[[1]]$cellData
    arg.gtFile <-  cnf$datasets[[1]]$trueEdges
    
    url.inputFolder = paste0("/scratch/as15096/eric",'/',cnf$input_dir,'/',cnf$dataset_dir,'/',cnf$datasets[[1]]$name ,'/')
    #url.output - url.inputFolder
    #url.output = paste0("/scratch/as15096/eric",'/outputs/',arg.dataFolder,'/',arg.dname,'/CICT' )
    
    url.input = paste0(url.inputFolder, arg.inFile)
    url.rawedgefile = paste0(url.inputFolder,'rawEdges.csv')
    
    url.name.map = paste0(url.inputFolder,'name_map.csv')
    
    
    #**************************************************************
    
    print('Operation: Calculating raw edges ===========')
    
    try({file.remove(url.rawedgefile)},silent=T)
    rm(similarityMatrices);arg.edgeTypes.existing = c()
    
    arg.edgeTypes.lst = str_split(arg.edgeTypes,',') %>% unlist() %>% str_subset('.{2}') %>% trim()
    arg.edgeTypes.abs = arg.edgeTypes.lst[which(! arg.edgeTypes.lst %in% arg.edgeTypes.existing)]
    
    
    #produce raw edges
    rm(similarityMatrices.new)
    source('Algorithms/CICT/requirements/calculateRawEdges.R')
    if(length(arg.edgeTypes.abs) > 0 ){
      similarityMatrices.new = calculateRawEdges(arg.edgeTypes.abs,url.input,n.workers  =10)
      
      if(exists('similarityMatrices') & nrow(similarityMatrices.new)>0 ) 
      {  
        
        similarityMatrices = 
          similarityMatrices %>% select(-colnames(similarityMatrices.new),src ,trgt) %>%
          merge(similarityMatrices.new,all=T,by=c('src','trgt')) 
      }else similarityMatrices = similarityMatrices.new
      
      try({write.table(similarityMatrices, url.rawedgefile, sep = ",", quote = FALSE, row.names = FALSE)})
      print(sprintf('%s, %s raw edge calculations (%s) finished. %s', 
                    arg.dataFolder, arg.dname, paste0(arg.edgeTypes.lst,collapse = ', ') ,Sys.time()))
    }
    
    
    
  }else 
    if(operation=='runCICT_par'){
      print('step 003')
      #Only takes two arguments one for config file and one for operation
      library(yaml)
      #configFile='parConf_1.yaml'
      #cnf = read_yaml('/scratch/as15096/eric/outputs/cict_par/sens_modeling_choices/parConf_1.yaml')
      #configFilePath = paste0("/scratch/as15096/eric",'/outputs/cict_par/',  'sens_sparsity','/parConf_211.yaml')
      
      
      cnf = read_yaml(configFilePath) 
      
      print(toString(cnf))
      
      arg.dname <- cnf$datasets[[1]]$name  #
      cictRawEdgeCol<-cnf$CICT$edgeType
      
      #Reading other information that does not affect scaling, from the config file, 
      # allows less complex code and customizability without code modification
      arg.inputsdir <-  cnf$input_dir
      arg.dataFolder <- cnf$dataset_dir
      arg.expLevel <- cnf$expLevel
      
      arg.inFile <- cnf$datasets[[1]]$exprData
      arg.pseudoTime <-  cnf$datasets[[1]]$cellData
      arg.gtFile <-  cnf$datasets[[1]]$trueEdges
      
      url.inputFolder = paste0("/scratch/as15096/eric",'/',cnf$input_dir,'/',cnf$dataset_dir,'/',cnf$datasets[[1]]$name ,'/')
      url.input = paste0(url.inputFolder, arg.inFile)
      url.rawedgefile = paste0(url.inputFolder,'rawEdges.csv')
      
      url.gt = paste0(url.inputFolder,arg.gtFile)
      url.name.map = paste0(url.inputFolder,'name_map.csv')
      
      url.output = paste0("/scratch/as15096/eric",'/outputs/cict_par/',  cnf$experiment) # arg.dataFolder,'/',arg.dname,'/CICT' )
      dir.create(url.output,recursive = T)
      
      url.rankedEdges = paste0(url.output , '/rankedEdges',cnf$given_job_id,'.csv')
      url.rankedEdgesGated = paste0(url.output , '/rankedEdgesGated',cnf$given_job_id,'.csv')
      url.randomPreds = paste0(url.output,'/randomRankedEdges.csv')
      FLAG_exportRankedEdges=T
      
      url.results= paste0(url.output,"/rslt_",cnf$given_job_id,".rds")
      if(file.exists(url.results) & forceOutput==T)
        if(file.size(url.results) >600 ) {
          print("===========> Results already exist (size >500 bytes) and forceOutput = FALSE")
          return(-1)
        }
      
      options(debug=NULL)
      
      rm(rcrd)
      rcrd =CICT(cnf$given_job_id,
                 url.input,
                 url.rawedgefile,
                 url.name.map,
                 url.gt,
                 url.output,
                 url.logfile="noLog"  ,
                 cictRawEdgeCol,
                 earlyThresholdForGraphAnalysis = cnf$CICT$earlyThresholdForGraphAnalysis,
                 minGroundTruth.ratio.learning =cnf$CICT$minGroundTruth.ratio.learning,
                 maxunseenTest.ratio = cnf$CICT$maxunseenTest.ratio,
                 maxGroundTruth = cnf$CICT$maxGroundTruth,
                 randomEdgesFoldCausal =cnf$CICT$randomEdgesFoldCausal,
                 sampling.c_rc.ratio = cnf$CICT$sampling.c_rc.ratio,
                 trainingTarget=cnf$CICT$trainingTarget,
                 tstPrecent =cnf$CICT$tstPrecent,
                 forceOutput=forceOutput,
                 FLAG_runOnAllEdges = T,
                 FLAG_exportRankedEdges = FLAG_exportRankedEdges,
                 FLAG_exportTrainAndTest=F,
                 RF_max_depth =cnf$CICT$RF_max_depth,
                 RF_ntrees=cnf$CICT$RF_ntrees,
                 Debug=F,
                 maxNetSize = cnf$CICT$maxNetSize)
      #rcrd$job  = job
      
      
      
      try({saveRDS(rcrd,file=url.results)})
      print(paste0("===========> Results saved successfully at: ",url.results ))
      #rcrd$jobDesc  = data$jobDesc
    }else 
      if(operation=='runCICT' | operation=='calcEdges' | operation=='runSupervised')
    {
      
      library(yaml)
      #cnfyaml = read_yaml('/scratch/as15096/eric/outputs/cict_par/sens_modeling_choices/parConf_1.yaml')
      
      #STEP 1  reads the mother config file *******************
      cnfyaml = read_yaml(paste0(url.base,'/config-files/',configFilePath)) #configFilePath='config_L2.yaml'
      inputslocation ="inputs_beeline2" # "inputs"
      ds = cnfyaml$input_settings$datasets %>% rbindlist() 
      databases = str_subset(ds$name, "dream",negate= T)  # removing dreams data
      #databases='hHep'
      #databases = c("mHSC-GM" ,"mHSC-L")
      for(idx in 1:length(databases))  #idx =1
      { #arg.dname in databases){ #
        #arg.dname <- args$input_settings$datasets[[1]][['name']]
        #databases = c('dream5_1','dream5_3','dream5_4','mDC','mHSC-E','mESC','hESC','hHEP')
        arg.dname = databases[idx]#[['name']]
        arg.dataFolder <- cnfyaml$input_settings$dataset_dir[[1]]
        arg.inFile <- cnfyaml$input_settings$datasets[[idx]][['exprData']]
        arg.pseudoTime <-  cnfyaml$input_settings$datasets[[idx]][['cellData']]
        arg.gtFile <-  cnfyaml$input_settings$datasets[[idx]][['trueEdges']] #ground truth
        arg.edgeTypes <- cnfyaml$input_settings$algorithms[[1]]$params[['edgeTypes']] # comma separated  'Pearson, ewMImm'
        
        url.inputbase = paste0("/scratch/as15096/eric",'/',inputslocation,'/',arg.dataFolder)
        url.inputFolder = paste0(url.inputbase ,'/',arg.dname) #,'/CICT')
        
        url.input = paste0(url.inputFolder,'/' , arg.inFile) #arg.dataFolder,'_',
        url.rawedgefile = paste0(url.inputbase,'/',arg.dname,'/rawEdges.csv')
        url.gt = paste0(url.inputbase,'/',arg.dname,'/',arg.gtFile)
        url.name.map = paste0(url.inputbase,'/',arg.dname,'/','name_map.csv')
        
        
        url.output = paste0("/scratch/as15096/eric",'/outputs/',  arg.dataFolder,'/',arg.dname) # arg.dataFolder,'/',arg.dname,'/CICT' )
        dir.create(url.output,recursive = T)
        
        # supervised.positiveClass <- args[5] # c:causal edges, 'c,rc': causal and reversecausal
        # supervised.negativeClass<- args[6]  # r:random edges, 'rc': reversecausal  
        # supervised.gtProp<- args[7] #proportion of GT used for classification
        # supervised.train<- args[8] #proportion of classification records used for training
        # fp.fn.cost<- args[9]  #False positive to False Negative cost ratio for threshold identification
        
        try({
          print('===================================================================')
          print(paste0(idx, ' - PROCESSING:   ', arg.dname, '----------------------'))
          
          
          # Calculate raw edges
          if(operation=='calcEdges'){
            #data/inputs/L1/DREAM5_1/CICT/ExpressionData.csv data/outputs/L1/DREAM5_1/CICT/outFile.txt
            #TODO check if raw edges are not ready prepare it
            recalcSimMatrices=forceOutput
            
            url.output = paste0("/scratch/as15096/eric",'/outputs/',arg.dataFolder,'/',arg.dname,'/CICT' )
            #Adds a subfolder to CICT for outputs of an experiemntal run, e.g. different params
            if(!is.na(arg.experiment)) url.output=paste0(url.output,'/',arg.experiment) 
            
            if(!dir.exists(url.output)) dir.create(url.output,recursive = T)
            url.logfile = paste0(url.output , '/',arg.dname,'_cict log.txt')
            write('Start',  file = url.logfile, append = F)
            
            msg = c("Calculating raw edges =====" )
            cat(msg)
            write(msg,file = url.logfile, append = TRUE)
            
            print('Operation: Calculating raw edges ===========')
            

            if(file.exists(url.rawedgefile) & !recalcSimMatrices)
            {
              print(paste0("RawEdge file exist and overwriting not enforced. Exiting. ", url.rawedgefile))
              
              similarityMatrices = fread(file= url.rawedgefile)
              arg.edgeTypes.existing = colnames(similarityMatrices)
              next
            } else {
              try({file.remove(url.rawedgefile)},silent=T)
              rm(similarityMatrices);arg.edgeTypes.existing = c()
            }
            
            
            
            arg.edgeTypes.lst = str_split(arg.edgeTypes,',') %>% unlist() %>% str_subset('.{2}') %>% trim()
            arg.edgeTypes.abs = arg.edgeTypes.lst[which(! arg.edgeTypes.lst %in% arg.edgeTypes.existing)]
            
            
            #produce raw edges
            rm(similarityMatrices.new)
            source('Algorithms/CICT/requirements/calculateRawEdges.R')
            if(length(arg.edgeTypes.abs) > 0 ){
              similarityMatrices.new = calculateRawEdges(arg.edgeTypes.abs,url.input,n.workers  =10)
              
              if(exists('similarityMatrices') & nrow(similarityMatrices.new)>0 ) 
              {  
                
                similarityMatrices = 
                  similarityMatrices %>% select(-colnames(similarityMatrices.new),src ,trgt) %>%
                  merge(similarityMatrices.new,all=T,by=c('src','trgt')) 
              }else similarityMatrices = similarityMatrices.new
              
              try({write.table(similarityMatrices, url.rawedgefile, sep = ",", quote = FALSE, row.names = FALSE)})
              try({write.table(similarityMatrices, paste0(url.output,'/rawEdges.csv'), sep = ",", quote = FALSE, row.names = FALSE)})
              print(sprintf('%s raw edge calculations (%s) finished. %s', 
                            arg.dname, paste0(arg.edgeTypes.lst,collapse = ', ') ,Sys.time()))
            }
            
          }
          
          # Run CICT Predict all edges
          if(operation=='runCICT'){
            #data/inputs/L1/DREAM5_1/CICT/ExpressionData.csv data/outputs/L1/DREAM5_1/CICT/outFile.txt
            #TODO check if raw edges are not ready prepare it
            #
            options(knitr.table.format = "pipe")
            #CICT parameter setting
            {
              cictRawEdgeCol ='Pearson' # 'efMIempirical' #  'ewMIshrink' # 
              earlyThresholdForGraphAnalysis=0#An ealy threshold lets the graph to produce more useful information
              minGroundTruth.ratio.learning = .2
              maxunseenTest.ratio = .1
              maxGroundTruth =500
              randomEdgesFoldCausal = 5
              sampling.c_rc.ratio = .3   #Percentage of casual& reversecausal sampled 
              
              trainingTarget = "class2"
              tstPrecent = .3
              
            }
            
            theJobID=''
            
            if(theJobID !='') theJobID.tmp = paste0(theJobID,'_') else theJobID.tmp=''
            url.outputFolder =  paste0(url.output ,'/CICT/',theJobID.tmp) 
            #url.outputFolder = "/scratch/as15096/eric/outputs/L2_lofgof/mESC/mESC_lofgof_training_sets/2/"
            url.rankedEdges = paste0(url.output ,'/CICT/',theJobID.tmp, 'rankedEdges.csv')
            url.rankedEdgesGated = paste0(url.output ,'/CICT/',theJobID.tmp, 'rankedEdgesGated.csv')
           
            if(!dir.exists(url.output)) dir.create(url.output,recursive = T)
            url.logfile = paste0(url.output , '/',arg.dname,'_cict log.txt')
            write('Start',  file = url.logfile, append = F)
            
            
            print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
            if(file.exists(url.rankedEdges) & !forceOutput) {
              print(paste0("Results exist and overwriting not enforced. Exiting. ", url.rankedEdges))
              next
            } else print (url.rankedEdges)
              
            rcrd=CICT(theJobID,
                      url.input,
                      url.rawedgefile,
                      url.name.map,
                      url.gt,
                      url.output,
                      url.logfile=url.logfile, #"noLog"  ,
                      cictRawEdgeCol,
                      earlyThresholdForGraphAnalysis,
                      minGroundTruth.ratio.learning,
                      maxunseenTest.ratio,
                      maxGroundTruth,
                      randomEdgesFoldCausal,
                      sampling.c_rc.ratio,
                      trainingTarget='class2',
                      tstPrecent = .3,
                      forceOutput=forceOutput,
                      FLAG_runOnAllEdges =T,
                      FLAG_exportRankedEdges = T,
                      FLAG_exportTrainAndTest = T,
                      Debug=F,
                      RF_max_depth =10,RF_ntrees=50)
            


            try({saveRDS(rcrd,file=paste0(url.outputFolder,"_rslt_",".rds"))})
              
              #When learning set is provided e.g. for mESC_lofGOF
            if(F){  
              url.output.back = url.output
              for(i in 1:10){
                  print("=============================================================")
                  print(i)
                  preset.train=NA; preset.test=NA
                  path.learningset = paste0(url.output.back,'/mESC_lofgof_training_sets/',i)
                  url.outputFolder = path.learningset
                  url.rankedEdges = paste0(path.learningset , '/rankedEdges.csv')
                  url.rankedEdgesGated = paste0(path.learningset , '/rankedEdgesGated.csv')
                  url.randomPreds = paste0(path.learningset,'/randomRankedEdges.csv')
                  
                  
                  preset.train= read.csv(paste0(path.learningset,'/training.csv'))
                  preset.test= read.csv(paste0(path.learningset,'/test.csv'))
                  
                  
                  
                  # preset.test =  preset.test%>%
                  #   mutate(Gene1 = tolower(Gene1), Gene2 = tolower(Gene2)) %>%
                  #   mutate(
                  #     outcomes =1,
                  #     Gene1 = str_replace_all(Gene1, "^\\w{1}", toupper) %>% unlist(),
                  #     Gene2 = str_replace_all(Gene2, "^\\w{1}", toupper) %>% unlist()
                  #   )     
                  
                  url.output=path.learningset
                  rcrd=CICT(theJobID,
                           url.input,
                           url.rawedgefile,
                           url.name.map,
                           url.gt,
                           url.output,
                           url.logfile=url.logfile, #"noLog"  ,
                           cictRawEdgeCol,
                           earlyThresholdForGraphAnalysis,
                           minGroundTruth.ratio.learning,
                           maxunseenTest.ratio,
                           maxGroundTruth,
                           randomEdgesFoldCausal,
                           sampling.c_rc.ratio,
                           trainingTarget='class2',
                           tstPrecent = .3,
                           forceOutput=forceOutput,
                           FLAG_runOnAllEdges =T,
                           FLAG_exportRankedEdges = T,
                           FLAG_exportTrainAndTest = T,
                           Debug=F,
                           RF_max_depth =10,RF_ntrees=20,
                           preset.train=preset.train, preset.test=preset.test)
                
                  try({saveRDS(rcrd,file=paste0(path.learningset,"/_rslt_",i,".rds"))})
              }
  
              }
          }
          
          # Run CICT Predict all edges
          if(operation=='runSupervised'){
            
            arg.experiment = 'TF-aware'
            additional_provided_cols ='edgeGoingFromTF' # '' # 
            
            RF_max_depth =5;RF_ntrees=20
            
            cictRawEdgeCol = 'Pearson'             
            theJobID=''
            
            if(theJobID !='') theJobID.tmp = paste0(theJobID,'_') else theJobID.tmp=''
            url.outputFolder =  paste0(url.output ,'/RandomForest/',theJobID.tmp) 
            #url.outputFolder = "/scratch/as15096/eric/outputs/L2_lofgof/mESC/mESC_lofgof_training_sets/2/"
            url.rankedEdges = paste0(url.output ,'/RandomForest/',theJobID.tmp, 'rankedEdges.csv')
            url.rankedEdgesGated = paste0(url.output ,'/RandomForest/',theJobID.tmp, 'rankedEdgesGated.csv')
            
            #Adds a subfolder to CICT for outputs of an experiemntal run, e.g. different params
            if(!is.na(arg.experiment)) url.output=paste0(url.outputFolder,arg.experiment) 
            if(!dir.exists(url.output)) dir.create(url.output,recursive = T)
            
            url.logfile = paste0(url.output , '/',arg.dname,'_supervised log.txt')
            write('Start',  file = url.logfile, append = F)
            
            
            msg = c("Operation: Running Supervised methods ============================ \n",
                    sprintf("Data= %s | Config= %s | operation= %s | experiment= %s | ground truth= %s \n",
                            arg.dname,configFilePath ,operation,arg.experiment,arg.gtFile ))
            cat(msg)
            write(msg,file = url.logfile, append = TRUE)
            
            
            print('step 01')
            
            earlyThresholdForGraphAnalysis=0#An ealy threshold lets the graph to produce more useful information
            minGroundTruth.ratio.learning = .2
            maxunseenTest.ratio = .1
            maxGroundTruth = 5000
            randomEdgesFoldCausal = 5
            sampling.c_rc.ratio = .3   #Percentage of casual& reversecausal sampled 
            
            trainingTarget = "class2"
            tstPrecent = .3
            
            tbl.goldStandard=NULL
            tbl.goldStandard = fread(url.gt,fill=TRUE )  
            if(ncol(tbl.goldStandard)==2)tbl.goldStandard =tbl.goldStandard%>%  mutate(edgetruth=1)
            try({tbl.goldStandard=tbl.goldStandard%>%  rename(edgetruth=Type) },silent = T)
            try({tbl.goldStandard=tbl.goldStandard%>%  rename(src=Gene1, trgt = Gene2) },silent = T)
            
            try({    
              tbl.goldStandard=tbl.goldStandard%>% 
                select(-any_of('V3')) %>%
                dplyr::filter(edgetruth==1 | edgetruth== '+') %>%
                mutate(edgetyptruth='chipseq')
            },silent = T)
            
            name.map=n.itm.e=all.dt=NULL
            print('step 02')
            if(file.exists(url.rawedgefile) )
            {
              n.itm.e = fread(url.rawedgefile) #%>% mutate(src=toupper(src), trgt = toupper(trgt) ) 
              if(file.exists(url.name.map)) {
                try({name.map = fread(url.name.map) })
                if(!is.null(name.map))
                  suppressMessages({
                    n.itm.e=n.itm.e %>% mutate(
                      src = mapvalues(src,name.map$`#ID`,name.map$Name),
                      trgt = mapvalues(trgt,name.map$`#ID`,name.map$Name))
                    
                   })
              
              } 
              
              if(additional_provided_cols == 'edgeGoingFromTF' ) 
                n.itm.e = n.itm.e %>% mutate((!!additional_provided_cols) := ifelse(src %in% tbl.goldStandard$src,1,0))
              #table(n.itm.e$edgeGoingFromTF)
              
              if(! cictRawEdgeCol %in% colnames(n.itm.e))
                stop("Requested cictRawEdgeCol does not exist in rawEdges.csv and n.itm.e ")
              
              all.dt = fread(url.input)  
              try({all.dt = all.dt %>% rename(vID = V1)})
              Nobservations=ncol(all.dt)-1
  
              print('step 03')
              

              #intersect(n.itm.e$src , tbl.goldStandard$src);intersect(n.itm.e$trgt , tbl.goldStandard$trgt)
              
              # print(paste0('Number of edges from n.itm.e in goldstandard = ' ,
              #              nrow(n.itm.e %>% inner_join(tbl.goldStandard, by = c('src'='src','trgt'='trgt')))
              # )
              # )
             
              
              
              
              try({saveRDS(rcrd,file=paste0(url.outputFolder,"_rslt_",".rds"))})
              
            } else {
              print(sprintf("RawEdge file %s does not exist",url.rawedgefile))
              rm(similarityMatrices);arg.edgeTypes.existing = c()
            }
            
            print(paste0('Calculating supervised methods outputs '))
            
            source('Algorithms/CICT/Supervised.R')
            
            
          }
          

        })
        # 
      }
    }





#
if(F)
{
  datasets = c("hESC","hHep","mESC" ,"mDC","mHSC-E" ,"mHSC-GM","mHSC-L"  )
  Levels = c('L0','L1','L2','L0_ns','L1_ns','L2_ns')
  
  i=0
  
  for(lvl in Levels)
    for(dt in datasets){
      url.input= paste0("/scratch/as15096/eric",'/inputs_beeline2/',lvl,'/',dt)
      url.output = paste0("/scratch/as15096/eric",'/outputs/',  lvl,'/',dt,'/CICT') # arg.dataFolder,'/',arg.dname,'/CICT' )
      url.rawedges.i = paste0(url.input,'/rawEdges.csv')
      url.rawedges.o = paste0(url.output,'/rawEdges.csv')
      
      i=i+1   
      #Report raw edges
      if(T)
      {
        if(file.exists(url.rawedges.o))
        {
          if(file.exists(url.rawedges.i) & 
             file.info(url.rawedges.o)$ctime > file.info(url.rawedges.i)$ctime ) {
                file.remove(paste0(url.rawedges.i,'.back')) 
                file.rename(url.rawedges.i, paste0(url.rawedges.i,'.back')) 
                if(!dir.exists(url.input)) dir.create(url.input,recursive = T)
                file.copy(url.rawedges.o,url.rawedges.i)
             }
         }
  
        if(file.exists( url.rawedges.i))
          print(paste0(lvl, ' ', dt, ' ',  file.info(url.rawedges.i)$ctime)) else
            print(paste0(lvl, ' ', dt, ' is missing <<<<-----------'))
      }

      #report learningset export
      if(F){
        url.train =  paste0(url.output,'/train.csv')
        if(!file.exists(url.train))
        {
          print(paste0('Train ', lvl, ' ', dt, ' is missing <<<<-----------'))
        }
      }
    }
  
}

#Calcualte correlation in parallel
if(F)
{
  
  library(doSNOW)
  #library(parallelly)
  library(doFuture)
  library(parallel)
  data <- read.csv(url.input) 
  data= data[!duplicated(data$gene),] %>% column_to_rownames('gene')
  data[is.na(data)] <- 0
  data1=data %>% mutate_all(as.numeric) %>% t()
  which(!sapply(data1,is.numeric,simplify = T))
  
  nCores <- detectCores() %>% min(10)
  cl <- parallel::makeCluster(nCores)
  
  cor_matrix <- parApply(cl, data1 , 2, function(x,d) cor(x, d),data1)
  row.names(cor_matrix) <-colnames(cor_matrix)
  cor_matrix.lng = cor_matrix %>% as.data.frame() %>% 
    rownames_to_column('Gene1') %>% pivot_longer(everything() & -Gene1,names_to = 'Gene2')
  
  write.csv(cor_matrix.lng,paste0(url.output,'/','rawEdges_', arg.inFile))
  stopCluster(cl)
  
  
  #NEEDS debugging
  # require(WGCNA)
  # allowWGCNAThreads()
  # system.time({ mim <- WGCNA::corFast(ExpressionMatrix)})
  
  
}

