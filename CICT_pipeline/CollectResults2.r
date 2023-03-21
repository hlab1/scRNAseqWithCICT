#Libraries ----
{
  #library(mlbench)
  library(caret)
  #library(Spectrum)
  library(Rtsne)
  library(umap)
  library(devtools)
  library(hutils)
  #library(PerformanceAnalytics)
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
  
  library(fitdistrplus) #fit distributions
  library(gendist)
  library(sn)
  library(ggplot2)
  library(SimilarityMeasures)
  library(PerformanceAnalytics)
  library(Lmoments)
  #library(TSdist) #timeseries distances
  #library(frechet)
  
  
  library(sna)
  # for building network and visualization 
  library(tidygraph)
  library(graphlayouts)
  library(statnet)
  library(ggplot2)
  library(igraph)
  library(yaml)
  PSUDEO_ZERO = 0.000001
  
}

f.all = list.files(path = '/scratch/as15096/eric/outputs',pattern = 'rankedEdges.csv', recursive = T)

if(F)
  for(f.r in f.all){
    # f=paste0('/scratch/as15096/eric/outputs/',f.r)
    # if(str_detect(dirname(f),'/CICT')) {
    #   tbl = fread(f)
    #   fwrite(tbl,f,sep = '\t',row.names = F)
    # }
    
    f=paste0('/scratch/as15096/eric/inputs/',f.r)
    path.splt = str_split(str_replace(f,'/rankedEdges.csv',''), '/') %>% unlist() 
    db=path.splt[length(path.splt)-1]
    algo = path.splt %>% last()
    
    tbl = fread(f)
    tbl$Type=1
    fwrite(tbl,f,sep = '\t',row.names = F)
  }





#Performance check
#Testing performance
{
  require(verification)
  require(pROC)
  require(ROCR)
  require( OptimalCutpoints)
  require(precrec )
  reportAUC<-function(x,prcsn=7)
  {
    a= attr(x,'aucs');
    b=attr(x,'paucs') %>% rename(aucs=paucs,standardized = spaucs) %>% 
      mutate(curvetypes = paste0('p',curvetypes))
    
    rbindlist(list(a,b),fill=T,use.names = T) %>%select(-modnames, -dsids) %>%
      mutate(aucs = round(aucs,prcsn), standardized = round(standardized,prcsn) ) %>%
      knitr::kable()
  }
  
  
  reportPerf<-function(dataFolder,db,algo, url.allresults )
  {
    
    evals = data.frame()
    
    url.input =paste0("/scratch/as15096/eric",'/inputs_beeline2/')
    url.refnet = paste0(url.input,dataFolder,'/',db,'/refNetwork.csv' )
    outcome = fread(url.refnet) %>% rename(outcomes = Type)
    table(outcome$outcomes)
    
    url.output =paste0("/scratch/as15096/eric",'/outputs/')
    url.pred = paste0(url.output,dataFolder,'/',db, '/',algo,'/rankedEdges.csv' )
    prd = fread(url.pred) %>% rename(predictions = EdgeWeight)
    
    assespreds = prd %>% full_join(outcome, by = c('Gene1', 'Gene2')) %>%
      mutate(outcomes = ifelse(is.na(outcomes),0,outcomes),
             predictions = ifelse(is.na(predictions),0,predictions))
    assespreds=assespreds %>% mutate(rndPred=ifelse(runif(nrow(assespreds)) >=.5,1,0))
    
    #table(assespreds$predictions,assespreds$outcomes)
    
    url.logfile  = paste0(url.output,dataFolder,'/',db, '/',algo,'_evaluation.txt')
    
    
    msg = c('================================================',
            "Reporting final results",
            capture.output(h2o.performance(tst1.mdl,valid=T)))
    cat(paste0(msg,collapse='\n') )
    write( msg,  file = url.logfile, append = TRUE, sep = '\n')
    
    
    #TODO remove the reverse edges before assessment, just keep the one with higher prediction score
    theROC <- roc(assespreds$outcomes, assespreds$predictions, percent = TRUE);theROC
    evals$roc1 = theROC
    
    sscurves <- evalmod(scores = assespreds$predictions, labels = assespreds$outcomes)
    pr.crv=sscurves$prcs[1][[1]];  pr.auc =  attr(pr.crv,'auc')
    evals$pr.crv=pr.crv;evals$pr.auc=pr.auc
    
    theauc = precrec::auc(sscurves); 
    msg = paste0(theauc[[3]],"=", round(theauc[[4]],4));cat(msg)
    write(msg,  file = url.logfile, append = TRUE)
    ##autoplot(sscurves, mode = 'rocpr',ret_grob = F)
    
    #partial precision-recall
    sscurves.part <- part(sscurves, xlim = c(0, 0.2))
    ##plot(sscurves.part,title = "Partial top 20% predictions")
    write('Early precision recall ======',  file = url.logfile, append = TRUE)
    write(reportAUC(sscurves.part),  file = url.logfile, append = TRUE)
    evals$sscurves.part =reportAUC(sscurves.part)  
    
    pr.prtcrv=sscurves.part$prcs[1][[1]];  pr.prtauc =  attr(pr.prtcrv,'pauc')
    evals$pr.prtcrv =reportAUC(pr.prtcrv) 
    
    #sink(file = url.logfile, append = T, type = c("output"),split = T) #write(  file = url.logfile, append = TRUE)
    
    #random classifier
    print("Random Classifier comparison ============")
    write("Random Classifier comparison ============",  file = url.logfile, append = TRUE)
    
    randomClassifierCurves <- evalmod(scores = assespreds$rndPred, labels = assespreds$outcomes)
    rnd.crv=randomClassifierCurves$prcs[1][[1]]; rnd.auc =  attr(rnd.crv,'auc')
    
    rndmClscurves.part <- part(randomClassifierCurves, xlim = c(0, 0.2))
    reportAUC(rndmClscurves.part)
    write('Random classifier precision recall ======',  file = url.logfile, append = TRUE)
    write(reportAUC(rndmClscurves.part),  file = url.logfile, append = TRUE)
    
    rnd.prtcrv=rndmClscurves.part$prcs[1][[1]]; rnd.prtauc =  attr(rnd.prtcrv,'pauc')
    
    msg = sprintf("%s: AUCPR Ratio to Random= %s,  Partial AUCPR Ratio to Random= %s ", 
                  algo,  round(pr.auc/rnd.auc,3), round(pr.prtauc/rnd.prtauc,3))
    print(msg)
    write(msg,  file = url.logfile, append = TRUE)
    
    mmpoins <- evalmod(scores = assespreds$predictions, labels = assespreds$outcomes, mode = "basic")
    #autoplot(mmpoins, c("error", "accuracy"))
    
    # Show normalized ranks vs. specificity, sensitivity, and precision
    ##autoplot(mmpoins, c("specificity", "sensitivity", "precision"))
    
    # Show normalized ranks vs. Matthews correlation coefficient and F-score
    ##autoplot(mmpoins, c("mcc", "fscore"))
    
    #bestcutoff = coords(roc, "best", ret="threshold", transpose = FALSE)
    
    # the relative cost of of a false negative classification (as compared with a false positive classification)
    # the prevalence, or the proportion of cases in the population (n.cases/(n.controls+n.cases)).
    relativeCostfn_fp = 1/5
    prv = table(setDT(tst1.totalset)[,get(trainingTarget)])
    best.weights=c(relativeCostfn_fp, prv[2]/prv[1])  
    bestcutoff =as.double(coords(theROC, "best", best.method = "closest.topleft", best.weights=best.weights,
                                 ret="threshold", transpose = FALSE));bestcutoff
    print(sprintf("Best threshold: %s for fn/fp cost ratio: %s " ,bestcutoff,relativeCostfn_fp ))
    write(sprintf("Best threshold: %s for fn/fp cost ratio: %s " ,bestcutoff,relativeCostfn_fp ),  file = url.logfile, append = TRUE)
    
    #bestcutoff = 0.5
    assespreds =assespreds %>% mutate(thresholdpreds= ifelse(assespreds$predictions>bestcutoff,assespreds$predictions,0)) 
    theROC <- roc(assespreds$outcomes, assespreds$thresholdpreds, percent = TRUE);theROC
    
    sscurves <- evalmod(scores = assespreds$thresholdpreds, labels = assespreds$outcomes)
    pr.crv=sscurves$prcs[1][[1]];  pr.auc =  attr(pr.crv,'auc')
    msg = sprintf("With FP-FN ratio 0.1 => AUCPR Ratio CICT to Random= %s,  AUCPR Ratio CICT to Random= %s ", 
                  round(pr.auc/rnd.auc,3), round(pr.prtauc/rnd.prtauc,3))
    print(msg)
    write(msg,  file = url.logfile, append = TRUE)
    
    #sink(NULL)
    # 
    # pander::pander(ftable(factor(assespreds$thresholdpreds, c("TRUE","FALSE")), 
    #                       factor(assespreds$outcomes, c("TRUE","FALSE")), 
    #                       dnn=c('pred','actual'))) #%>%  knitr::kable()
    
    #table(ifelse(predictions>.99,'c','other'))
    #bestcutoff = coords(roc, "best", best.method = "closest.topleft", ret="threshold", transpose = FALSE);bestcutoff
    
  }
  
  dataFolder= 'L0'
  url.allresults  = paste0(url.output,dataFolder,'/', 'evaluation.csv')
  f.all = list.files(path = '/scratch/as15096/eric/outputs',pattern = 'rankedEdges.csv', recursive = T)
  
  f.allexcept.par = f.all[str_detect(f.all,'cict_par',negate = T)]
  for(f.r in f.allexcept.par){
    f=paste0('/scratch/as15096/eric/inputs/',f.r)
    path.splt = str_split(str_replace(f,'/rankedEdges.csv',''), '/') %>% unlist() 
    
    db=path.splt[length(path.splt)-1]
    algo = path.splt %>% last()
    
    
    reportPerf(dataFolder,db,algo,d.gs, url.allresults )
    
    
    tbl = fread(f)
    tbl$Type=1
    fwrite(tbl,f,sep = '\t',row.names = F)
    
  }
  
}


#Extract and integrate results from batch opration
{
  #Function definitions
  {
    extractResourcesInf <- function(fldr,ptrn = 'time_.*txt'){
      out=data.frame()
      f.all = list.files(path = fldr,pattern =ptrn , recursive = F)
      f.r=f.all[1]
      for(f.r in f.all){
        
        
        job_idx = str_extract(f.r,pattern = '(?<=[-]).*(?=[.]txt)') %>% as.integer()
        
        f=paste0(fldr,f.r)
        a=read_file(f) %>% str_split('\n\t') %>% unlist() 
        b= llply(a,function(x) {
          
          parts = str_split(x,":",simplify = T);
          v=trim(parts[2])
          names(v)<-parts[1]
          v
        }) %>% unlist()
        names(b) = str_replace_all(names(b),'[ ()]','_')
        d=as.data.frame(t(b))
        d$job_idx = job_idx
        
        out=bind_rows(out,d)
        
      }
      
      out
    }
    
    #/scratch/ch153/packages/BEELINE/hlab1/Beeline/outputs/L2/mESC/GENIE3/
    extractResourcesInf.beeline <- function(fldr,ptrn = 'time.*txt'){
      out=data.frame()
      f.all = list.files(path = fldr,pattern =ptrn , recursive = F)
      f.r=f.all[1]
      for(f.r in f.all){
        
        
        #job_idx = str_extract(f.r,pattern = '(?<=[-]).*(?=[.]txt)') %>% as.integer()
        
        f=paste0(fldr,f.r)
        a=read_file(f) %>% str_split('\n\t') %>% unlist() 
        b= llply(a,function(x) {
          
          parts = str_split(x,":",simplify = T);
          v=trim(parts[2])
          names(v)<-parts[1]
          v
        }) %>% unlist()
        names(b) = str_replace_all(names(b),'[ ()]','_')
        d=as.data.frame(t(b))
        d$job_idx = job_idx
        
        out=bind_rows(out,d)
        
      }
      
      out
    }
    
    
    extractNestedInf <- function(p,prefix = 'unsn',type=''){
      
      defaultname = paste0(type,'_empty')
      out =   data.frame(t1=1) %>%rename(!!defaultname:=t1)
      if(is.null(p)) return(out)
      if(!is.list(p) & is.vector(p) & length(p)<=1) return(out)
      try({
        if(type=='stats'){
          #p = r$unseensmpl_stats
          p.names = p %>% str_extract_all(pattern='(?<=[|]|^).*?(?=[=])') %>% unlist ()%>% trim() 
          p.vals= p %>% str_extract_all(pattern='(?<=[=]).*?(?=[|]|$)') %>% unlist ()%>% trim()
          names(p.vals) = p.names
          out=rbind(p.vals)
        }
        if(type=='auc_table'){
          #p = r$unseensmpl_part
          p.auc = p$aucs
          names(p.auc) <- paste0(prefix, "_",p$curvetypes)
          p.auc_std = p$standardized
          names(p.auc_std) =  paste0(prefix, "_std_",p$curvetypes)
          out=rbind(c(p.auc,p.auc_std))
          
        }
        if(type=='max_criteria'){
          #p = r$max_criteria
          p.threshold = p$threshold
          names(p.threshold) <- paste0(prefix, "_",p$metric,"_threshold")
          p.val = p$value
          names(p.val) =  paste0(prefix, "_",p$metric,"")
          out=rbind(c(p.threshold,p.val))
          
        }
        
        
        if(type=='cm'){ #confusion matric
          if(!is.null(p$cm)) p =  p$cm[1:2,1:2] else p = p[1:2,1:2]
          p.cm = unlist(p)
          
          names(p.cm) = c('tn','fp','fn','tp')
          out=rbind(p.cm)
        }
      })
      out
    }
    
    flattenResults<-function(r){
      r.classes = sapply(r,class)
      b=r[ r.classes %in% c('character','integer',"numeric","logical","double") ] 
      d=sapply(b,function(x) {ifelse( is.na(as.numeric(x)), x, as.numeric(x))},simplify = T)%>%  cbind() %>% t()
      
      d = d %>% as.data.frame() %>%
        cbind(extractNestedInf(r$cm,'cm','cm' )) %>%
        cbind(extractNestedInf(r$unseensmpl_stats,'unsn','stats' )) %>%
        cbind(extractNestedInf(r$unseensmpl_part,'unsn','auc_table' )) %>%
        cbind(extractNestedInf(r$unseensmpl_rndm,'rndm','auc_table' )) %>%
        cbind(extractNestedInf(r$max_criteria,'crt','max_criteria' ))
      
      d
    }
  }
  
  #Extract results 
  {
    analysisFolder='cict_par' # 'SERGIO.*'
    path.base = "/scratch/as15096/eric/outputs/"
    r.path = paste0(path.base,analysisFolder,"/")
    #  "sens_edgeType/"
    r.files = list.files(r.path,pattern = '.*[.]rds',recursive = T)
    r.files = str_subset(r.files, '.*_par',negate = T)
    r.all =  data.frame()
    for(r.file in r.files){
      #print(r.file)
      
      r=data.frame()
      try({
        r = readRDS(paste0(r.path, r.file))
        if(length(r)>1){
          print(r.file)
          r.flat =  flattenResults(r)  #%>%as.data.frame()
          r.flat =r.flat %>% as.data.frame() %>% 
            mutate(job_idx = str_extract(r.file,pattern = '(?<=[_]).*(?=[.]rds)') %>%as.integer())
          r.all = rbindlist(list(r.all,r.flat),fill = T,use.names = T)
          
        }
      })
    }
  }
  #Extract and integrate results from batch opration
  
  chck<-function(val,fallback = ''){
    if(is.null(val)) return(fallback)
    if(purrr::is_empty(val)) return(fallback)
    if(length(val)==0) return(fallback)
    val
  }
  
  #Batch opration
  {
    library(doRNG)
    library(foreach)
    library(doFuture)
    registerDoFuture();
    plan(multisession, split = TRUE)
    
    library(progressr)
    handlers(global = TRUE)
    handlers("progress", "beepr")
    
    rm(r.resources,r.evals,r.all,r.final)
    analysisFolder = 'sens_sparsity'   # "sens_multipleRuns" # 'sens_modeling_choices'  #  'sens_edgeType'  #  'sens_learningTrgtAndSize' #
    
    path.base = "/scratch/as15096/eric/outputs/cict_par/"
    r.path = paste0(path.base,analysisFolder,"/")
    
    cnf.files = list.files(r.path,pattern = '.*[.]yaml')
    r.evals = data.frame()
    for(r.file in cnf.files){
      tryCatch({
        cnf.d = data.frame()
        #r.file= 'parConf_416.yaml'
        cnf = read_yaml(paste0(r.path, r.file))
        
        if(length(cnf)>1){
          cnf.d=data.frame(  evaluation = chck(cnf$experiment),
                             parConfFile =chck(cnf$parConf),
                             job_idx = as.integer(cnf$given_job_id),
                             expLevel = ifelse(is.null(cnf$expLevel), cnf$output_prefix , cnf$expLevel),
                             dataset.dir = chck(cnf$dataset_dir),
                             dataset = chck(cnf$datasets[[1]]$name),
                             edgeType = chck(cnf$CICT$edgeType),
                             size.groundTruth = chck(cnf$CICT$maxGroundTruth),
                             size.randomEdgeFoldCausal = cnf$CICT$randomEdgesFoldCausal,
                             trainingTarget = chck(cnf$trainingTrgt , 'class2') ,
                             RF_ntree = chck(cnf$CICT$RF_ntrees,NA),
                             RF_max_depth = chck(cnf$CICT$RF_max_depth,NA)
          )
        }
        
      },error =function(e) {print(r.file);print(e)})
      r.evals =rbindlist(list(r.evals,cnf.d),fill = T,use.names = T)  #rbind(r.evals,cnf.d)
    }
    print("Reading conf files finished")
    
    
    #  "sens_edgeType/"
    r.files = list.files(r.path,pattern = '.*[.]rds')
    r.all =  data.frame()

    
    r.all =foreach(r.file = r.files,.combine = rbind)%dopar%
    {
      r=data.frame()
      try({
        r = readRDS(paste0(r.path, r.file))
        if(length(r)>1){
          #i=i+1;if(i%%50 ==0) print(paste0(i , " result files processed."))
          r.flat =  flattenResults(r)  #%>%as.data.frame()
          r.flat =r.flat %>% as.data.frame() %>% 
            mutate(job_idx = str_extract(r.file,pattern = '(?<=[_]).*(?=[.]rds)') %>%as.integer())
          r.flat
        }
      })
    }
    
    i=0
    if(F)for(r.file in r.files){
      #print(r.file)
      
      r=data.frame()
      try({
        r = readRDS(paste0(r.path, r.file))
        if(length(r)>1){
          r.flat =  flattenResults(r)  #%>%as.data.frame()
          r.flat =r.flat %>% as.data.frame() %>% 
            mutate(job_idx = str_extract(r.file,pattern = '(?<=[_]).*(?=[.]rds)') %>%as.integer())
          r.all = rbind(r.all,r.flat)
          
        i=i+1;if(i%%50 ==0) print(paste0(i , " result files processed."))
        
        }
      })
    }
    print("Reading result files finished")
    
    #Convert data.frame of list columns to standard data.frame
    lst <- lapply(r.all, unlist)
    r.all =data.frame(lapply(lst, `length<-`, median(lengths(lst))))
    
    
    # r.all = sapply(r.all, function(x) unlist(x) ) %>% as.data.frame() %>% 
    #   group_by(rslt) %>%slice_head(n=1) %>%ungroup()
    
    r.resources = extractResourcesInf(r.path)
    
    r.final = r.evals %>% left_join(r.all, by='job_idx') %>% left_join(r.resources,by = 'job_idx')
    glimpse(r.final)
    
    r.path.final= paste0(r.path,analysisFolder, ".csv")
    write.csv(r.final,file = r.path.final)
    
  }
  
  #r.exp = prd3  %>% left_join(r.all,by= c('rslt'='rslt' )
}

#Visualize results
if(F){
  library(metafor)
  
  ### copy BCG vaccine meta-analysis data into 'dat'
  dat = r.final %>% select()
  dat <- read.csv(r.path.final)
  
  ### calculate log risk ratios and corresponding sampling variances (and use
  ### the 'slab' argument to store study labels as part of the data frame)
  dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat,
                slab=paste(author, year, sep=", "))
  
  ### fit random-effects model
  res <- rma(yi, vi, data=dat)
  
  ### a little helper function to add Q-test, I^2, and tau^2 estimate info
  mlabfun <- function(text, res) {
    list(bquote(paste(.(text),
                      " (Q = ", .(formatC(res$QE, digits=2, format="f")),
                      ", df = ", .(res$k - res$p),
                      ", p ", .(metafor:::.pval(res$QEp, digits=2, showeq=TRUE, sep=" ")), "; ",
                      I^2, " = ", .(formatC(res$I2, digits=1, format="f")), "%, ",
                      tau^2, " = ", .(formatC(res$tau2, digits=2, format="f")), ")")))}
  
  ### set up forest plot (with 2x2 table counts added; the 'rows' argument is
  ### used to specify in which rows the outcomes will be plotted)
  forest(res, xlim=c(-16, 4.6), at=log(c(0.05, 0.25, 1, 4)), atransf=exp,
         ilab=cbind(tpos, tneg, cpos, cneg), ilab.xpos=c(-9.5,-8,-6,-4.5),
         cex=0.75, ylim=c(-1, 27), order=alloc, rows=c(3:4,9:15,20:23),
         mlab=mlabfun("RE Model for All Studies", res),
         psize=1, header="Author(s) and Year")
  
  ### set font expansion factor (as in forest() above) and use a bold font
  op <- par(cex=0.75, font=2)
  
  ### add additional column headings to the plot
  text(c(-9.5,-8,-6,-4.5), 26, c("TB+", "TB-", "TB+", "TB-"))
  text(c(-8.75,-5.25),     27, c("Vaccinated", "Control"))
  
  ### switch to bold italic font
  par(font=4)
  
  ### add text for the subgroups
  text(-16, c(24,16,5), pos=4, c("Systematic Allocation",
                                 "Random Allocation",
                                 "Alternate Allocation"))
  
  ### set par back to the original settings
  par(op)
  
  ### fit random-effects model in the three subgroups
  res.s <- rma(yi, vi, subset=(alloc=="systematic"), data=dat)
  res.r <- rma(yi, vi, subset=(alloc=="random"),     data=dat)
  res.a <- rma(yi, vi, subset=(alloc=="alternate"),  data=dat)
  
  ### add summary polygons for the three subgroups
  addpoly(res.s, row=18.5, mlab=mlabfun("RE Model for Subgroup", res.s))
  addpoly(res.r, row= 7.5, mlab=mlabfun("RE Model for Subgroup", res.r))
  addpoly(res.a, row= 1.5, mlab=mlabfun("RE Model for Subgroup", res.a))
  
  ### fit meta-regression model to test for subgroup differences
  res <- rma(yi, vi, mods = ~ alloc, data=dat)
  
  ### add text for the test of subgroup differences
  text(-16, -1.8, pos=4, cex=0.75, bquote(paste("Test for Subgroup Differences: ",
                                                Q[M], " = ", .(formatC(res$QM, digits=2, format="f")), ", df = ", .(res$p - 1),
                                                ", p = ", .(formatC(res$QMp, digits=2, format="f")))))
}


#extract results from OTHER methods
{
  #setDF(tst1.totalset)
  library(doRNG)
  library(foreach)
  library(doFuture)
  registerDoFuture();
  plan(multisession) #, split = TRUE
  #options(future.globals.maxSize= 30891289600)
  
  
  require(verification)
  require(pROC)
  require(ROCR)
  require( OptimalCutpoints)
  require(precrec )
  
  reportAUC<-function(x,prcsn=6)
  {
    a= attr(x,'aucs');
    b=attr(x,'paucs') %>% rename(aucs=paucs,standardized = spaucs) %>% 
      mutate(curvetypes = paste0('p',curvetypes))
    
    rbindlist(list(a,b),fill=T,use.names = T) %>%select(-modnames, -dsids) %>%
      mutate(aucs = round(aucs,prcsn), standardized = round(standardized,prcsn) ) 
  }
  
  url.outputs = '/scratch/ch153/packages/BEELINE/hlab1/Beeline/outputs/'
  url.inputs = '/scratch/as15096/eric/inputs_beeline2/'
  benchmarks = dir(url.outputs)[-1]
  
  
  r.errors = r.all=data.frame()
  for(bnch in benchmarks){ #bnch='L2_lofgof'
  #r.all.1=foreach(bnch = benchmarks)%dopar%{ r.all=data.frame() #bnch='L2'
          
    datafolders = dir(paste0(url.outputs,bnch))
    for(dfolder in datafolders) { #dfolder="mESC"
      algorithms =dir(paste0(url.outputs,bnch,'/',dfolder)) 
      algorithms = setdiff(algorithms,'CICT')#algorithms='SINGE'
      
      for(algo in algorithms){ #algo = "DEEPDRIM8"
        runfolders=''
        runfolders=paste0('run_',0:9)
        
        for(runf in runfolders){
          if(runf !='') run.phrase = paste0('/',runf) else run.phrase=''  #For multiruns mESC under the algo folder
          
          r.err = data.frame( bnch, dataset = dfolder, algorithm = algo)
          rcrd = list()
          r.flat = data.frame()
          
          tryCatch(
            { #/scratch/ch153/packages/BEELINE/hlab1/Beeline/outputs/L2_lofgof/mESC/DEEPDRIM8
            url = paste0(url.outputs,bnch,'/',dfolder,'/',algo,run.phrase,'/rankedEdges.csv')
            if(!file.exists(url.outputs)) {
              r.flat$error =paste0('RankedEdges not found: ', url.refnetwork)
            }
            rawFile = read.csv(url,sep='\t')
            rankedEdges= rawFile %>% select('Gene1','Gene2','EdgeWeight') # str_split(unlist(rawFile[,1]),'\t',simplify = T) %>% as.data.frame() #   ldply(.fun = matrix,nrow=1)
            colnames(rankedEdges) <- c('src','trgt','predictions')
            
            #url.refnetwork ="/scratch/ch153/packages/BEELINE/hlab1/Beeline/inputs_beeline2/L2/mHSC-E/refNetwork.csv"
            url.refnetwork = paste0(url.inputs,bnch,'/',dfolder,'/','refNetwork.csv')
            if(!file.exists(url.refnetwork)) {
              r.flat$error =paste0('Refnetwork not found: ', url.refnetwork)
              
            }
            refnetwork= read.csv(url.refnetwork)
            refnetwork = refnetwork %>% mutate(outcomes = 1) %>%
              rename(src=Gene1,trgt=Gene2)
            # if(ncol(refnetwork)==2) refnetwork$outcomes = 1
            # colnames(refnetwork) <- c('src','trgt','outcomes')
            #
            
            pred_outcome = rankedEdges %>% full_join(refnetwork, by=c("src"="src", "trgt"="trgt") )#left_join
            assespreds = pred_outcome %>% mutate(
              predictions = as.numeric(predictions),
              outcomes = ifelse(is.na(outcomes),0,outcomes))
            
            #random classifier
            rndFile=data.frame()
            if(F){ #Empirical 
              assespreds = assespreds %>% mutate(rndPred = ifelse(runif(nrow(pred_outcome)) >=.5,1,0))
              
              #True edges only from outgoing edges of TFs
              if(F) 
                assespreds = assespreds %>% 
                      mutate(rndPred =  ifelse(runif(nrow(pred_outcome)) >=.5,
                                            ifelse(src %in% refnetwork$src,1,0),
                                            0)
                                 )
              #r.all$ratio = r.all$unsn_pPRC/r.all$rndm_pPRC;paste0(sort(round(r.all$ratio,2)),collapse=' | ')
              
            }else
              {
              url.rnd = paste0(url.outputs,bnch,'/',dfolder,'/','RANDOM/rankedEdges.csv')
              if(!file.exists(url.rnd)) {
                r.flat$error =paste0('Random RankedEdges not found: ', url.refnetwork)
              }
              
              rndFile = read.csv(url.rnd,sep = '\t')
              assespreds=assespreds %>% 
                full_join(rndFile, by = c('src'='Gene1','trgt'='Gene2')) %>% 
                rename(rndPred = EdgeWeight) %>%
                mutate(rndPred = ifelse(is.na(rndPred),0,rndPred),
                       predictions = ifelse(is.na(predictions),0,predictions),
                       outcomes = ifelse(is.na(outcomes),0,outcomes)) 
              #%>%                 mutate(rndPred = ifelse(rndPred<0.5,0,1))
                
            }          
            
            assespreds=assespreds %>% mutate(outcomes= as.factor(outcomes))
            {
              sscurves <- evalmod(scores = assespreds$predictions, labels = assespreds$outcomes)
              pr.crv=sscurves$prcs[1][[1]];  pr.auc =  attr(pr.crv,'auc')
              
              theauc = precrec::auc(sscurves); 
              msg = paste0(theauc[[3]],"=", round(theauc[[4]],3))
              rcrd$unseensmpl_roc_pr = msg
              
              
              #partial precision-recall
              sscurves.part <- part(sscurves, xlim = c(0, 0.2))
              rcrd$unseensmpl_part = as.data.frame(reportAUC(sscurves.part) )
              
  
              randomClassifierCurves <- evalmod(scores = assespreds$rndPred, labels = assespreds$outcomes)
              rnd.crv=randomClassifierCurves$prcs[1][[1]]; rnd.auc =  attr(rnd.crv,'auc')
              
              rndmClscurves.part <- part(randomClassifierCurves, xlim = c(0, 0.2))
              reportAUC(rndmClscurves.part) 
              rnd.prtcrv=rndmClscurves.part$prcs[1][[1]]; rnd.prtauc =  attr(rnd.prtcrv,'pauc')
              rcrd$unseensmpl_rndm = as.data.frame(reportAUC(rndmClscurves.part) )
              
              # mmpoints <- evalmod(scores = assespreds$predictions, labels = assespreds$outcomes, mode = "basic")
              # rcrd$mcc = mmpoints$mcc
            }
            
            
            rcrd$cm = ''
            rcrd$max_criteria = ''
            rcrd$unseensmpl_stats = ''
            r.flat =  flattenResults(rcrd)  #%>%as.data.frame()
            r.flat =r.flat %>% as.data.frame() 
            
          }, 
            error=function(e){
              r.err = r.err %>% mutate(error = toString(e))
            }
          )
          
          tryCatch(
            {
              r.flat = r.flat %>% as.data.frame() %>%
                dplyr::mutate(benchmark = bnch,
                              dataset = dfolder,
                              algorithm = algo,
                              error = '')
              
              r.all <- rbind(r.all,r.flat)
            }, 
          error=function(e){
            r.err = r.err %>% mutate(error = toString(e))
          }
          )
          
          r.errors = rbind(r.errors,r.err)
        }
        try({print(sprintf("%s - %s - %s  Done.", bnch,dfolder,algo))})
        
      }
    }
    #r.all
  } #%dorng% 
  #r.all.2 = rbindlist(r.all.1)
  saveRDS(r.all,'/scratch/as15096/eric/outputs/cict_par/otherMethodsPerformance_5.rds' )
  #saveRDS(r.all,'/scratch/as15096/eric/outputs/cict_par/DEEPDRIM8_multirun1.rds' )
  #REsources used by other methods
  r.errors =r.resources= data.frame()
  idx=0
  for(bnch in benchmarks){ #bnch='L2'
    #r.all.1=foreach(bnch = benchmarks)%dopar%{ r.all=data.frame() #bnch='L2'
    
    datafolders = dir(paste0(url.outputs,bnch))
    for(dfolder in datafolders) { #dfolder="mESC"
      algorithms =dir(paste0(url.outputs,bnch,'/',dfolder)) 
      algorithms = setdiff(algorithms,'CICT')#algorithms='SINGE'
      
      for(algo in algorithms){ #algo ="GENIE3" # "DEEPDRIM6"
        
        r.err = data.frame( bnch, dataset = dfolder, algorithm = algo)

        tryCatch(
          {
            idx=idx+1
            url = paste0(url.outputs,bnch,'/',dfolder,'/',algo,'/')
            r.res = extractResourcesInf(url,ptrn='time.*txt')
            
            r.res$idx  = idx
            r.res$algorithm = algo
            
            r.res$benchmark = bnch
            r.res$dataset = dfolder
            r.resources = rbind(r.resources,r.res)
            
          }, 
          error=function(e){
            r.err = r.err %>% mutate(error = toString(e))
          }
        )
        
        tryCatch(
          {
            r.flat = r.flat %>% as.data.frame() %>%
              dplyr::mutate(benchmark = bnch,
                            dataset = dfolder,
                            algorithm = algo,
                            error = '')
            
            r.all <- rbind(r.all,r.flat)
          }, 
          error=function(e){
            r.err = r.err %>% mutate(error = toString(e))
          }
        )
        
        r.errors = rbind(r.errors,r.err)
        
        try({print(sprintf("%s - %s - %s  Done.", bnch,dfolder,algo))})
        
      }
    }
    #r.all
  } #%dorng% 
  
  #r.all.2 = rbindlist(r.all.1)
  saveRDS(r.resources,'/scratch/as15096/eric/outputs/cict_par/otherMethodsResources_5.rds' )
}


  #extract results from multiple train/test set evaluations mESC
  {
    analysisFolder='mESC_lofgof_training_sets' # 'SERGIO.*'
    path.base = paste0("/scratch/as15096/eric/outputs/L2_lofgof/mESC/")
    r.path = paste0(path.base,analysisFolder,"/")
    #  "sens_edgeType/"
    r.files = list.files(r.path,pattern = '.*[.]rds',recursive = T)
    r.files = str_subset(r.files, '.*_par',negate = T)
    r.all =  data.frame()
    for(r.file in r.files){
      #print(r.file)
      
      r=data.frame()
      try({
        r = readRDS(paste0(r.path, r.file))
        if(length(r)>1){
          print(r.file)
          r.flat =  flattenResults(r)  #%>%as.data.frame()
          r.flat =r.flat %>% as.data.frame() %>% 
            mutate(job_idx = str_extract(r.file,pattern = '(?<=[_]).*(?=[.]rds)') %>%as.integer())
          r.all = rbindlist(list(r.all,r.flat),fill = T,use.names = T)
          
        }
      })
    }
    
    saveRDS(r.all,paste0(path.base,analysisFolder,"/all_eval.rds") )
    
  }

  #extract results from Random forest evaluations
  {
    analysisFolder='mESC_lofgof_training_sets' # 'SERGIO.*'
    path.base = paste0("/scratch/as15096/eric/outputs/L2_lofgof/mESC/")
    r.path = paste0(path.base,analysisFolder,"/")
    #  "sens_edgeType/"
    r.files = list.files(r.path,pattern = '.*[.]rds',recursive = T)
    r.files = str_subset(r.files, '.*_par',negate = T)
    r.all =  data.frame()
    for(r.file in r.files){
      #print(r.file)
      
      r=data.frame()
      try({
        r = readRDS(paste0(r.path, r.file))
        if(length(r)>1){
          print(r.file)
          r.flat =  flattenResults(r)  #%>%as.data.frame()
          r.flat =r.flat %>% as.data.frame() %>% 
            mutate(job_idx = str_extract(r.file,pattern = '(?<=[_]).*(?=[.]rds)') %>%as.integer())
          r.all = rbindlist(list(r.all,r.flat),fill = T,use.names = T)
          
        }
      })
    }
    
    saveRDS(r.all,paste0(path.base,analysisFolder,"/all_eval.rds") )
    
  }

  #extract results from Random forest evaluations
  {
    analysisFolder='mESC_lofgof_training_sets' # 'SERGIO.*'
    path.base = paste0("/scratch/as15096/eric/outputs/L2_lofgof/mESC/")
    r.path = paste0(path.base,analysisFolder,"/")
    #  "sens_edgeType/"
    r.files = list.files(r.path,pattern = '.*[.]rds',recursive = T)
    r.files = str_subset(r.files, '.*_par',negate = T)
    r.all =  data.frame()
    for(r.file in r.files){
      #print(r.file)
      
      r=data.frame()
      try({
        r = readRDS(paste0(r.path, r.file))
        if(length(r)>1){
          print(r.file)
          r.flat =  flattenResults(r)  #%>%as.data.frame()
          r.flat =r.flat %>% as.data.frame() %>% 
            mutate(job_idx = str_extract(r.file,pattern = '(?<=[_]).*(?=[.]rds)') %>%as.integer())
          r.all = rbindlist(list(r.all,r.flat),fill = T,use.names = T)
          
        }
      })
    }
    
    saveRDS(r.all,paste0(path.base,analysisFolder,"/all_eval.rds") )
    
  }

  #extract results from Random forest text files
  {
    path.base = paste0("/scratch/as15096/eric/outputs/")
    r.files = list.files(path.base,pattern = 'supervised log',recursive = T)
    r.files.filtered = str_subset(r.files, 'RandomForest/Not-TF-aware|RandomForest/TF-aware')
    r.all = data.frame()
    for(r.file in r.files.filtered){
      try({
        print(r.file)
        path.full = paste0(path.base, r.file)
        r = readLines(path.full)
        r1 = str_subset(r, '^AUCPR Ratio')
        if(length(r1)>=1){
          network.density = str_subset(r,"Network density") %>%
            str_extract("(?<=Network density\\= )[0-9.e-]{1,8}") %>% as.numeric() 
          
          rAUPR = str_extract(r1,"(?<=AUCPR Ratio to Random\\= )[0-9.e-]{1,8}")%>% as.numeric() 
          rpAUPR = str_extract(r1,"(?<=Partial AUCPR Ratio to Random\\= )[0-9.e-]{1,8}")%>% as.numeric() 
          groundTruth = str_extract(path.full,"(?<=/scratch/as15096/eric/outputs/).+?(?=/)")
          data = str_extract(r.file,"(?<=(L2|L2_ns|L2_lofgof)/).+?(?=/)")
          informed = str_extract(r.file,"(?<=RandomForest/).+?(?=/)")
          
          try({r.flat = data.frame(groundTruth,data,informed,rAUPR,rpAUPR, network.density,r.file)})
          r.all = rbindlist(list(r.all,r.flat),fill = T,use.names = T)
        }
      })
    }
    
    saveRDS(r.all,paste0(path.base,"L2_supervised_eval.rds") )
    
  }

  #/scratch/ch153/packages/BEELINE/hlab1/Beeline/outputs/L2/mESC/DEEPDRIM7/
  #extract results from  DEEPDRIM's multiple train/test set evaluations mESC
  {
    analysisFolder='mESC_lofgof_training_sets' # 'SERGIO.*'
    path.base = paste0("/scratch/ch153/packages/BEELINE/hlab1/Beeline/outputs/L2/mESC/DEEPDRIM7/")
    r.path = paste0(path.base,analysisFolder,"/")
    #  "sens_edgeType/"
    r.files = list.files(r.path,pattern = '.*[.]rds',recursive = T)
    r.files = str_subset(r.files, '.*_par',negate = T)
    r.all =  data.frame()
    for(r.file in r.files){
      #print(r.file)
      
      r=data.frame()
      try({
        r = readRDS(paste0(r.path, r.file))
        if(length(r)>1){
          print(r.file)
          r.flat =  flattenResults(r)  #%>%as.data.frame()
          r.flat =r.flat %>% as.data.frame() %>% 
            mutate(job_idx = str_extract(r.file,pattern = '(?<=[_]).*(?=[.]rds)') %>%as.integer())
          r.all = rbindlist(list(r.all,r.flat),fill = T,use.names = T)
          
        }
      })
    }
    
    saveRDS(r.all,paste0(path.base,analysisFolder,"/all_eval.rds") )
    
  }
  
  #extract RESOURCE usage from  BEELINE's
  {
    
    
    url.outputs = '/scratch/ch153/packages/BEELINE/hlab1/Beeline/outputs/'
    url.inputs = '/scratch/as15096/eric/inputs_beeline2/'
    benchmarks = dir(url.outputs)[-1]
    
    
    r.errors = resources=data.frame()
    for(bnch in benchmarks){ #bnch='L2'
      #r.all.1=foreach(bnch = benchmarks)%dopar%{ r.all=data.frame() #bnch='L2'
      
      datafolders = dir(paste0(url.outputs,bnch))
      for(dfolder in datafolders) { #dfolder="mESC"
        algorithms =dir(paste0(url.outputs,bnch,'/',dfolder)) 
        algorithms = setdiff(algorithms,'CICT')#algorithms='SINGE'
        
        for(algo in algorithms){ #algo = "DEEPDRIM6"
          
          r.err = data.frame( bnch, dataset = dfolder, algorithm = algo)
          rcrd = list()
          r.flat = data.frame()
          tryCatch(
            {
              url = paste0(url.outputs,bnch,'/',dfolder,'/',algo,'/')
              resource = extractResourcesInf.beeline(url,ptrn = 'time.*txt')
              resource$algo = algo
              resource$dataset =dfolder
              resource$expLevel=bnch
              
              resources = rbindlist(list(resources,resource),use.names = T,fill=T)
              
            }, 
            error=function(e){
              r.err = r.err %>% mutate(error = toString(e))
              #browser()
            }
          )
          r.errors = rbind(r.errors,r.err)
          
#          try({print(sprintf("%s - %s - %s  Done.", bnch,dfolder,algo))})
          
        }
      }
      #r.all
    } #%dorng% 

    saveRDS(resources,paste0("/scratch/as15096/eric/outputs/beeline_resources.rds") )
    
  }

  #extract CICT_Scaling RESOURCE usage
{
  url = '/scratch/as15096/eric/outputs/cict_par/cict_scaling/'
  url.inputs = '/scratch/as15096/eric/inputs_beeline2/'
  bnch = 'cict_scaling'   ;dfolder='mESC'; algo='CICT'
  resource = extractResourcesInf(url,ptrn = 'time.*txt')
  saveRDS(resource,paste0("/scratch/as15096/eric/outputs/cict_par/cict_scaling/resources.rds") )
  
}