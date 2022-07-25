#Set working directory and paths




recalcSimMatrices=T
url.base = "/scratch/as15096/eric"
setwd(url.base)

#TODO source all libraries and functions and set up working director
  source('/scratch/as15096/eric/Algorithms/CICT/requirements/CICT_LibsFunctions.R')
  
  args.cmnd = c('config_SERGIO_DS4.yaml','runCICT_par','FALSE') #runCICT_par  #runCICT_par
  
  args.cmnd <- commandArgs(trailingOnly = T)
  
  if(length(args.cmnd)<2 ){
    print('Please provide [config file name e.g. config_L0.yaml] [calcEdges|runCICT|runCICT_par] [FALSE|TRUE]') 
    return(-1)
  } 
  
  configFile <- args.cmnd[1]
  operation<-args.cmnd[2]
  forceOutput<- as.logical(args.cmnd[3])
  if(is.na(forceOutput)) forceOutput = FALSE
  arg.experiment <- args.cmnd[4]
  
if(!operation %in% c('calcEdges','runCICT','runSupervised','runCICT_par','install')){
  print('Invalid operation requested')
  return(-1)
} else print(sprintf("Operation: %s, Config: %s, forceOutputReplace: %s", 
                     configFile,operation , forceOutput))


library(yaml)
args = read_yaml(paste0(url.base,'/config-files/',configFile)) #config_L0.yaml
inputslocation ="inputs_beeline2" # "inputs"


arg.dname='hHep'

ds = args$input_settings$datasets %>% rbindlist() 
databases = str_subset(ds$name, "dream",negate= T)  # removing dreams data
#databases = c('dream5_1','dream5_3','dream5_4','mDC','mHSC-E','mESC','hESC','hHEP')
#databases = c('hHep')

#CICT caller definition
{
  CICT<-function(url.input,
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
                 ...)
  {
    
    name.map=n.itm.e=all.dt=tbl.goldStandard=NULL
    arguments <- list(...)
    #Define urls/paths
    rcrd=list()
    rcrd$outputAlreadyPresent=F
    
    #Adds a subfolder to CICT for outputs of an experiemntal run, e.g. different params
    if(!is.na(arg.experiment)) url.output=paste0(url.output,'/',arg.experiment) 
    
    if(!dir.exists(url.output)) dir.create(url.output,recursive = T)
    if(url.logfile=='') url.logfile = paste0(url.output , '/',arg.dname,'_cict log.txt')
    write('Start',  file = url.logfile, append = F)
    url.rankedEdges = paste0(url.output , '/rankedEdges.csv')
    url.rankedEdgesGated = paste0(url.output , '/rankedEdgesGated.csv')
    url.randomPreds = paste0(url.output,'/randomRankedEdges.csv')
    
    
    #reporting
    {      
      msg = c("Operation: Running CICT ============================ \n",
              sprintf("Data= %s | Config= %s | operation= %s | experiment= %s | ground truth= %s \n",
                      arg.dname,configFile ,operation,arg.experiment,arg.gtFile ))
      cat(msg)
      write(msg,file = url.logfile, append = TRUE)
      
      print('step 01')
    } 
    
    #Read rawedge file and goldstandard
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
      all.dt = fread(url.input) 
      colnames(all.dt)[1] <-'gene'
      Nobservations=ncol(all.dt)-1
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
    
    if( (!file.exists(url.rankedEdges) | !file.exists(url.rankedEdgesGated)) | forceOutput) {
      #try({   
      
      source('Algorithms/CICT/CICT.R',local=TRUE)
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
if(operation=='install'){
  install.packages('future.batchtools')
  install.packages('batchtools')
  install.packages("future")
  install.packages("listenv")
  return()
}

if(operation=='runCICT_par'){
  
  cict.wrapper = function(data, job, instance,outputSubFolder='',
                          earlyThresholdForGraphAnalysis,
                          minGroundTruth.ratio.learning,
                          maxunseenTest.ratio,
                          maxGroundTruth,
                          randomEdgesFoldCausal,
                          sampling.c_rc.ratio,
                          trainingTarget='class2',
                          tstPrecent = .3,
                          forceOutput=T,...) 
  {
    
    
    #Feeding datasetname, edgeType and run allows flexible adjusting of tasks 
    # and scaling up through problem assignment mechanism
    arg.dname <- data$dataset  #
    cictRawEdgeCol<-data$cictRawEdgeCol
    
    #Reading other information that does not affect scaling, from the config file, 
    # allows less complex code and customizability without code modification
    cnfg.url = data$configFile
    args = read_yaml(cnfg.url) #
    
    
    arg.inputsdir <-   args$input_settings$input_dir[[1]]
    arg.dataFolder <- args$input_settings$dataset_dir[[1]]
    
    arg.inFile <- args$input_settings$datasets[[1]][['exprData']]
    arg.pseudoTime <-  args$input_settings$datasets[[1]][['cellData']]
    arg.gtFile <-  args$input_settings$datasets[[1]][['trueEdges']] #ground truth
    
    url.inputFolder = paste0("/scratch/as15096/eric",'/',arg.inputsdir,'/',arg.dataFolder,'/',arg.dname,'/')
    url.input = paste0(url.inputFolder, arg.inFile)
    url.rawedgefile = paste0(url.inputFolder,'rawEdges.csv')
    
    url.gt = paste0(url.inputFolder,arg.gtFile)
    url.name.map = paste0(url.inputFolder,'name_map.csv')
    
    url.output = paste0("/scratch/as15096/eric",'/outputs/',  outputSubFolder) # arg.dataFolder,'/',arg.dname,'/CICT' )
    dir.create(url.output,recursive = T)
    
    #url.logfile="noLog"  
     
    rcrd =CICT(url.input,
               url.rawedgefile,
               url.name.map,
               url.gt,
               url.output,
               url.logfile="noLog"  ,
               cictRawEdgeCol,
               earlyThresholdForGraphAnalysis,
               minGroundTruth.ratio.learning,
               maxunseenTest.ratio,
               maxGroundTruth,
               randomEdgesFoldCausal,
               sampling.c_rc.ratio,
               trainingTarget='class2',
               tstPrecent = .3,
               forceOutput=T,...)
    #rcrd$job  = job
    
    rcrd$jobDesc  = data$jobDesc
    #records=append(records,rcrd)
    
    
    
    rcrd
  }
  
  #Test run cict.wrapper
  if(F){
    tmp.datasetnames = "hHep"
    tmp.arg.edgeTypes = "ewMImm"
    tmp=cict.wrapper(prd2[32,],job=NULL, instance=NULL,
                     outputSubFolder = 'cict_par',
                     earlyThresholdForGraphAnalysis,
                     minGroundTruth.ratio.learning,
                     maxunseenTest.ratio,
                     maxGroundTruth,
                     randomEdgesFoldCausal,
                     sampling.c_rc.ratio,
                     trainingTarget='class2',
                     tstPrecent = .3,
                     forceOutput=T)
  }
    
    
    
    #https://future.batchtools.futureverse.org/
    #https://mllg.github.io/batchtools/articles/batchtools.html
    #https://github.com/mllg/batchtools/blob/master/inst/templates/slurm-simple.tmpl
    
    #“apply algorithm A to problem instance P and store results”
    library(future.batchtools)
    library(batchtools)
    library("future")
    library("listenv")
    options(knitr.table.format = 'pipe')
    
    reg = makeExperimentRegistry(file.dir = NA, seed = 1)
    
    #Add problems and algorithms 
    {
      
      records = list()
      
      configFiles = list.files(paste0(url.base,'/config-files/'), pattern = 'L0.*yaml')
      configFiles=configFiles[1]
      #configFiles = 'config_SERGIO_DS4.yaml' ;cnfg=configFiles
      
      #datasetnames = NULL, # c('mDC','mHSC-E','mESC','hESC','hHEP'),
      #edgeTypes = c('Pearson','ewMImm' ),#'Euclidean','ewMIempirical',
      
      #problemDeisgns = list()
      prd2 = data.frame()
      
      runs=1 
      for(cnfg in configFiles ){
                       
        cnfg.url = paste0(url.base,'/config-files/',cnfg)
        args = read_yaml(cnfg.url) #
        
        cictRawEdgeCol =args$input_settings$algorithms[[1]]$params[['edgeTypes']] # comma separated  'Pearson, ewMImm'
        cictRawEdgeCol  = str_split(cictRawEdgeCol,',',simplify = T)  %>% trim() %>% str_subset(pattern=".{2}")
        
        datasetnames = rbindlist(args$input_settings$datasets) %>% pull('name')
        
        #Temporary, removing dream data from datasets
        datasetnames = str_subset(datasetnames,'dream',negate = T) 
        
        prd1 = purrr::cross(list(runs,cictRawEdgeCol,datasetnames))
        prd1 = rbindlist(prd1)
        colnames(prd1)<- c('runs','cictRawEdgeCol','dataset')
        
        prd1 = prd1 %>% mutate(
          dataFolder= args$input_settings$dataset_dir[[1]],
          configFile = cnfg.url
        ) %>% mutate(
          jobDesc = paste0(dataset,'_',cictRawEdgeCol,'_',runs)
        )
        
      
        prd2 = rbind(prd2,prd1)
        #addProblem(name = 'mDC', data = url.inputFolder , fun = subsample, seed = 42) #Feed config file to data
        #problemDeisgns = append(problemDeisgns,cnfg)
      }
      
      
      
      for(idx in 1:nrow(prd2)) 
        addProblem(name = paste0("_", idx), data = prd2[idx,] , fun = NULL, seed = 42) #Feed config file to data
      
      prd2$jobDesc
      #addAlgorithm(name = 'calcEdges', fun = NULL)
      addAlgorithm(name = 'runCICT', fun = cict.wrapper)
      #addAlgorithm(name = 'runOtherAlgos', fun = NULL) #AWE, ARACNE ...
      
      #Adding range of parameters of each algorithm
      algoDesign = list(
        runCICT = data.table( outputSubFolder = 'cict_par',
                              earlyThresholdForGraphAnalysis=0, #An ealy threshold lets the graph to produce more useful information
                              minGroundTruth.ratio.learning = .2,
                              maxunseenTest.ratio = .1,
                              maxGroundTruth =c(500) , #,1000), # c(500,1500,3000,5000),
                              randomEdgesFoldCausal = 5,
                              sampling.c_rc.ratio = .3 ,  #Percentage of casual& reversecausal sampled 
                              trainingTarget = "class2", #class3
                              tstPrecent = .3)
      )  
      
    }
    
    addExperiments(algo.designs = algoDesign,
                   combine = "crossprod",repls = 1)
    summarizeExperiments()
    
    id1 = head(findExperiments(algo.name = "runCICT"), 1)
    # testJob(id = id1)
    
    plan(list(
      tweak(batchtools_slurm, 
            resources = list( ntasks = 20,
                              cpu_per_task=6,
                              memory_per_cpu = "150gb"#vmem = "5gb"
            )),
      multisession
    ))
    
    system.time({ submitJobs() })
    ## Submitting 90 jobs in 90 chunks using cluster functions 'Interactive' ...
    waitForJobs()
    
    # reduce = function(res) list(mce = (sum(res) - sum(diag(res))) / sum(res))
    # results = unwrap(reduceResultsDataTable(fun = reduce))
    
    a=reduceResultsDataTable()
    results = unnest_wider(reduceResultsDataTable(),'result')
    
    
    job.inf = unwrap(getJobPars())
    job.res=(getJobResources(results$job.id, reg))
    job.cost = batchtools::getJobStatus()
    a=batchtools::getLog(id=1)
    
    results.all = cbind(job.inf,job.res,results)
    
    url.output = paste0("/scratch/as15096/eric",'/outputs/','cict_par')
    dir.create(url.output, recursive = T)
    url.cluster.results = paste0(url.output,'/clusterResults2', '.rds')
    
    saveRDS(results.all,file =url.cluster.results )
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
    
  
  }else for(idx in 1:length(databases)) 
    { #arg.dname in databases){
  #arg.dname <- args$input_settings$datasets[[1]][['name']]
  #JUST TESTING, REMOVETHIS
  #idx=2
  arg.dname = args$input_settings$datasets[[idx]][['name']]
  
  
  arg.dataFolder <- args$input_settings$dataset_dir[[1]]
  arg.inFile <- args$input_settings$datasets[[idx]][['exprData']]
  arg.pseudoTime <-  args$input_settings$datasets[[idx]][['cellData']]
  arg.gtFile <-  args$input_settings$datasets[[idx]][['trueEdges']] #ground truth
  arg.edgeTypes <- args$input_settings$algorithms[[1]]$params[['edgeTypes']] # comma separated  'Pearson, ewMImm'
  
  url.inputbase = paste0("/scratch/as15096/eric",'/',inputslocation,'/',arg.dataFolder)
  url.inputFolder = paste0(url.inputbase ,'/',arg.dname)
  
  url.input = paste0(url.inputFolder,'/' , arg.inFile)
  url.rawedgefile = paste0(url.inputbase,'/',arg.dname,'/rawEdges.csv')
  url.gt = paste0(url.inputbase,'/',arg.dname,'/',arg.gtFile)
  url.name.map = paste0(url.inputbase,'/',arg.dname,'/','name_map.csv')
  
  
  # supervised.positiveClass <- args[5] # c:causal edges, 'c,rc': causal and reversecausal
  # supervised.negativeClass<- args[6]  # r:random edges, 'rc': reversecausal  
  # supervised.gtProp<- args[7] #proportion of GT used for classification
  # supervised.train<- args[8] #proportion of classification records used for training
  # fp.fn.cost<- args[9]  #False positive to False Negative cost ratio for threshold identification
  
  try({
    print('===================================================================')
    print(paste0('PROCESSING:   ', arg.dname))
    
    
    # Calculate raw edges
    if(operation=='calcEdges'){
      #data/inputs/L1/DREAM5_1/CICT/ExpressionData.csv data/outputs/L1/DREAM5_1/CICT/outFile.txt
      #TODO check if raw edges are not ready prepare it
      
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
        similarityMatrices = fread(file= url.rawedgefile)
        arg.edgeTypes.existing = colnames(similarityMatrices)
        
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
        similarityMatrices.new = calculateRawEdges(arg.edgeTypes.abs,url.input,n.workers  =8)
        
        if(exists('similarityMatrices') & nrow(similarityMatrices.new)>0 ) 
        {  
          
          similarityMatrices = 
            similarityMatrices %>% select(-colnames(similarityMatrices.new),src ,trgt) %>%
            merge(similarityMatrices.new,all=T,by=c('src','trgt')) 
        }else similarityMatrices = similarityMatrices.new
        
        try({write.table(similarityMatrices, url.rawedgefile, sep = ",", quote = FALSE, row.names = FALSE)})
        try({write.table(similarityMatrices, paste0(url.output,'rawEdges.csv'), sep = ",", quote = FALSE, row.names = FALSE)})
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
        cictRawEdgeCol = 'ewMImm'
        earlyThresholdForGraphAnalysis=0#An ealy threshold lets the graph to produce more useful information
        minGroundTruth.ratio.learning = .2
        maxunseenTest.ratio = .1
        maxGroundTruth = 5000
        randomEdgesFoldCausal = 5
        sampling.c_rc.ratio = .3   #Percentage of casual& reversecausal sampled 
        
        trainingTarget = "class2"
        tstPrecent = .3
      }
      
       tmp=CICT(url.input,
               url.rawedgefile,
               url.name.map,
               url.gt,
               url.output,
               cictRawEdgeCol,
               earlyThresholdForGraphAnalysis,
               minGroundTruth.ratio.learning,
               maxunseenTest.ratio,
               maxGroundTruth,
               randomEdgesFoldCausal,
               sampling.c_rc.ratio,
               trainingTarget='class2',
               tstPrecent = .3,
               forceOutput=T)
      
    }
      
    
    
    # Run CICT Predict all edges
    if(operation=='runSupervised'){
      
      
      cictRawEdgeCol = 'Pearson'      
      url.output = paste0("/scratch/as15096/eric",'/outputs/',arg.dataFolder,'/',arg.dname,'' )
      #Adds a subfolder to CICT for outputs of an experiemntal run, e.g. different params
      if(!is.na(arg.experiment)) url.output=paste0(url.output,'/',arg.experiment) 
      
      if(!dir.exists(url.output)) dir.create(url.output,recursive = T)
      url.logfile = paste0(url.output , '/',arg.dname,'_supervised log.txt')
      write('Start',  file = url.logfile, append = F)
      
      
      msg = c("Operation: Running Supervised methods ============================ \n",
              sprintf("Data= %s | Config= %s | operation= %s | experiment= %s | ground truth= %s \n",
                      arg.dname,configFile ,operation,arg.experiment,arg.gtFile ))
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
      
      name.map=n.itm.e=all.dt=tbl.goldStandard=NULL
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
        } else print(sprintf("Expression file %s does not exist",url.rawedgefile))
        all.dt = fread(url.input) 
        Nobservations=ncol(all.dt)-1
        tbl.goldStandard = fread(url.gt,fill=TRUE )  
        print('step 03')
        
        try({tbl.goldStandard=tbl.goldStandard%>%  rename(edgetruth=Type) })
        
        try({    
          tbl.goldStandard=tbl.goldStandard%>% 
            rename(src = Gene1,trgt=Gene2) %>% select(-any_of('V3')) %>%
            dplyr::filter(edgetruth==1 | edgetruth== '+') %>%
            mutate(edgetyptruth='chipseq')
        })
        #intersect(n.itm.e$src , tbl.goldStandard$src);intersect(n.itm.e$trgt , tbl.goldStandard$trgt)
        
        print(paste0('Number of edges from n.itm.e in goldstandard = ' ,
                     nrow(n.itm.e %>% inner_join(tbl.goldStandard, by = c('src'='src','trgt'='trgt')))
                    )
              )
        
        
      } else {rm(similarityMatrices);arg.edgeTypes.existing = c()}
      
      print(paste0('Calculating supervised methods outputs '))
      source('Algorithms/CICT/Supervised.R')
      
      
    }
    
    
  })
}





