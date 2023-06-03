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
  
  PSUDEO_ZERO = 0.000001
  
}

f.all = list.files(path = '/scratch/as15096/eric/outputs',pattern = 'rankedEdges.csv', recursive = T)

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
  reportAUC<-function(x)
  {
    a= attr(x,'aucs');
    b=attr(x,'paucs') %>% rename(aucs=paucs,standardized = spaucs) %>% 
      mutate(curvetypes = paste0('p',curvetypes))
    
    rbindlist(list(a,b),fill=T,use.names = T) %>%select(-modnames, -dsids) %>%
      mutate(aucs = round(aucs,3), standardized = round(standardized,3) ) %>%
      knitr::kable()
  }
  
  
  reportPerf<-function(dataFolder,db,algo, url.allresults )
  {
    
      evals = data.frame()
      
      url.input =paste0("/scratch/as15096/eric",'/inputs/')
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
       
       
       #plot important features distributions
       if(F){
         library(plotly)
         library(GGally)
         
         df=setDF(d.new) %>% mutate(eID = paste0('src','_', 'trgt')) %>% select(c('class1',trainingTarget, prd.varimp$variable[1:10]))  #tst1.totalset #newDataEnv$tst1.totalset # 
         dfl = df %>% pivot_longer(prd.varimp$variable[1:10]) %>%  
           dplyr::mutate(logval = sign(value) * log2(abs(value))) 
         
         #df[,prd.varimp$variable[1:6]]%>%select_if(is.double)%>%ggpairs(  title = "important vars")
         
         #Parallelogram  Decision boundaries
         {
           theTitle =paste0(studyDataset, " - Dist. top predictors of causal and random edges")
           
           ggp.vars = prd.varimp$variable[1:15]
           ggp.vars = prd.varimp$variable[1:30]
           ggp.vars = setdiff(ggp.vars, str_subset(ggp.vars,"z"))[1:15]
           dfs = setDF(tst1.totalset) %>% select(unique(c('class1',trainingTarget, ggp.vars))) %>%
             group_by(class1) %>% summarise_if(is.numeric, list(myMedian) , na.rm = TRUE)
           
           ggp =ggparcoord(dfs, columns = 2:16, groupColumn = 'class1', scale = "std", splineFactor=0, 
                           title = "Top 20 CICT predictors of connection class")
           ggp + geom_line(size = 1.2)  + geom_point(size = 3)+ 
             theme_bw()+ theme(axis.text = element_text(size = 12,angle=-45,hjust=.2)) #, 
           }
         #Plotting 3D plot_ly(x=Sepal.Length,y=Sepal.Width,z=Petal.Length,type="scatter3d",mode='markers',size=Petal.Width,color=Species)
         {
           library(plotly)
           fig <- dfs %>%
             plot_ly(width = 1000, height = 600) 
           fig <- fig %>% add_trace(type = 'parcoords',
                                    line = list(color = ~class1,  #TODO set to prediction?
                                                colorscale = 'Jet',
                                                showscale = TRUE,
                                                reversescale = TRUE,
                                                cmin = FALSE,
                                                cmax = TRUE),
                                    dimensions = list(
                                      list(range = c(~min(confdisc),~max(confdisc)),
                                           #constraintrange = c(100000,150000),
                                           label = 'confdisc', values = ~confdisc),
                                      list(range = c(~min(ocbtau3.x),~max(ocbtau3.x)),
                                           label = 'ocbtau3.x', values = ~ocbtau3.x),
                                      list(#tickvals = c(0,0.5,1,2,3),
                                        #ticktext = c('A','AB','B','Y','Z'),
                                        label = 'scbNSkew.x', values = ~scbNSkew.x)
                                      # list(range = c(-1,4),
                                      #      tickvals = c(0,1,2,3),
                                      #      label = 'Block Material', values = ~blockMaterial),
                                      # list(range = c(~min(scbNSkew.x),~max(scbNSkew.x)),
                                      #      visible = TRUE,
                                      #      label = 'scbNSkew.x', values = ~scbNSkew.x),
                                      # list(range = c(~min(assemblyPW),~max(assemblyPW)),
                                      #      label = 'Assembly Penalty Weight', values = ~assemblyPW),
                                      # list(range = c(~min(HstW),~max(HstW)),
                                      #      label = 'Height st Width', values = ~HstW),
                                      # list(range = c(~min(minHW),~max(minHW)),
                                      #      label = 'Min Height Width', values = ~minHW),
                                      # list(range = c(~min(minWD),~max(minWD)),
                                      #      label = 'Min Width Diameter', values = ~minWD),
                                      # list(range = c(~min(rfBlock),~max(rfBlock)),
                                      #      label = 'RF Block', values = ~rfBlock)
                                    )
           )
           
           
           fig
         }
         
         
         #Violin
         #%%%%%%%%%%%%%%%%%%  Distribution of edges top predictors
         studyDataset = "HSEC timeset edges"
         
         theTitle =paste0(studyDataset, " - Dist. top predictors of causal and random edges")
         df=setDF(d.new) %>% mutate(eID = paste0('src','_', 'trgt')) %>% 
           mutate(class10 = class1 =='c' , class11= class1 =='u') %>% 
           select(c('class10','class11',trainingTarget, ggp.vars))  #tst1.totalset #newDataEnv$tst1.totalset # 
         
         dfl = df %>% pivot_longer(ggp.vars) %>%  
           dplyr::mutate(logval = sign(value) * log2(abs(value))) 
         
         myplotViolin(theTitle =theTitle,
                      dfl,ggp.vars[1:10], 
                      group1 = 'class10', group2 = 'class11' ,
                      grp1title= 'causal', grp2title = 'random',
                      grpcolor1 = "red",grpcolor2 = "darkgray",
                      yval = 'logval')
         
         
         
         
         
         #Sankey
         {
           p <- plot_ly(
             type = "sankey",
             domain = c(
               x =  c(0,1),
               y =  c(0,1)
             ),
             orientation = "h",
             valueformat = ".0f",
             valuesuffix = "TWh",
             
             node = list(
               label = dfl[class2 == 'c',]$name,
               color = dfl$class2,
               pad = 15,
               thickness = 15,
               line = list(
                 color = "black",
                 width = 0.5
               )
             ),
             
             link = list(
               source = json_data$data[[1]]$link$source,
               target = json_data$data[[1]]$link$target,
               value =  json_data$data[[1]]$link$value,
               color =  json_data$data[[1]]$link$color,
               label =  json_data$data[[1]]$link$label
             )
           ) %>% 
             layout(
               title = "Energy forecast for 2050, UK - Department of Energy & Climate Changes",
               font = list(
                 size = 10
               ),
               xaxis = list(showgrid = F, zeroline = F),
               yaxis = list(showgrid = F, zeroline = F)
             )
         }
         
         #Heatmap
         {
           
           myZ=function(x){(x-myMean(x))/mySD(x) }
           meanZ = function(x){ myMean(myZ(x))}
           
           setDT(dfl)
           
           df=setDF(d.new) %>% mutate(eID = paste0('src','_', 'trgt')) %>% select(c('class1',trainingTarget, prd.varimp$variable[1:30]))  #tst1.totalset #newDataEnv$tst1.totalset # 
           
           df.cols = str_subset(prd.varimp$variable,"^((?!^tz).)*$")[1:10]
           df.cols = prd.varimp$variable[1:20]
           
           dfs = d.new %>% group_by(class1) %>%  dplyr::summarise_at (df.cols,list(mn= meanZ))
           row.names(dfs)<-mapvalues(dfs$class1, c('c','rc','ir'), c('Causal','Reverse causal','Random'))
           tmp = scale(dfs[,2:ncol(dfs)], center = TRUE, scale = TRUE)
           
           
           
           # library(biclust)
           # tmp1=biclust(tmp,method=BCCC())
           # drawHeatmap(tmp,bicResult=tmp1,number=10) #,plotAll=FALSE)
           
           
           #df.cr = cor(setDT(tst1.tst)[,df.cols,with=FALSE])
           p <- plot_ly(x=colnames(tmp), y=rownames(dfs), 
                        z = tmp, 
                        type = "heatmap", 
                        colorscale= "RdYlBu",
                        showscale = T) %>%
             layout(margin = list(l=120))
           p
           
         }
         
         
         # p <- plot_ly(x=colnames(df.cr), y=rownames(df.cr), 
         #              z = df.cr, 
         #              type = "heatmap", 
         #              #colorscale= "Earth",
         #              showscale = F) %>%
         #   layout(margin = list(l=120))
         # p
         
         library(GGally)
         
         # Create data 
         df.cols = c(trg,str_subset(prd.varimp$variable,"^((?!^tz).)*$")[1:5])
         # Check correlations (as scatterplots), distribution and print corrleation coefficient 
         g = ggpairs(setDT(tst1.tst)[,df.cols,with=FALSE], title="correlogram with ggpairs()") 
         
         ggpairs(setDT(tst1.tst)[,df.cols,with=FALSE], columns = 2:4, ggplot2::aes(colour=species)) 
         
       }
       
       
       #Predicts in smaller chunks also adds random prediction for comparition

       #pred_outcome$predictions %>% table() 
       #write.csv(pred_outcome,file = 'CICT mHSC predictions for clustering.csv')
       #Late pruning using ground truth
       if(FALSE) pred_outcome = pred_outcome %>% filter(src %in% tf.all)
       
       #View(pred_outcome[predictions > .7 | outcomes ==TRUE,])
       # Pruning
       if(F){
         # Simple Removing one of the bidirectionals
         if(FALSE){
           pred_outcome = setDT(pred_outcome)[is.na(rvpred),rvpred:=0]
           table(pred_outcome[rvpred> predictions*1.5,]$outcome)
           pred_outcome[rvpred> predictions*1.5,predictions := 0]#[rvpred<=predictions,predictions := predictions+rvpred]
         }
         rm(d.new1,tmp1,d.new1.rv,tmp.rvpred)
         
         #removing based on inequality ?? TOOD correct
         if(FALSE){
           setDT(pred_outcome)
           setkeyv(pred_outcome,c('src','trgt'))
           data.table::setorderv(pred_outcome,cols='predictions',order = -1)
           inp = pred_outcome[,.(src,trgt,predictions,rvpred)]
           
           vrtcs = unique(n.itm.v$vID)
           cnt=0
           for(i in vrtcs)
           {
             for(j in vrtcs)
             {
               if(i==j) next
               for( k in vrtcs)
               {
                 if(i==k | j==k)next
                 shrtcut = inp[src ==i & trgt == k,]$predictions
                 a =inp[src ==i & trgt == j,]$predictions> shrtcut ; if(length(a)==0) a=FALSE
                 b= inp[src ==j & trgt == k,]$predictions> shrtcut; if(length(b)==0) b=FALSE
                 if( a & b ) inp[src ==i & trgt == k,drop:=T]
                 cnt=cnt+1
                 if(cnt %% 1000 ==0) print(cnt)
               }
             }
           }
           i;j;k
           
           
           
           pred_outcome.top2c = pred_outcome[, head(.SD, 2), by = "trgt"]
           pred_outcome.top2c=pred_outcome.top2c[,`:=`(predictions=NULL,outcomes=NULL,rvpred=NULL , is.causal1 =1)]
           d.new2 = merge(d.new1,pred_outcome.top2c,all.x=TRUE, by=c('src','trgt'))
           table(d.new2$outcomes,d.new2$is.causal1)
           
           
           pred_outcome.top1c = pred_outcome[, head(.SD, 1), by = "trgt"]
           pred_outcome.top1c=pred_outcome.top1c[,`:=`(predictions=NULL,outcomes=NULL,rvpred=NULL , is.causal1 =1)]
           
           d.new2 = merge(d.new1,pred_outcome.top1c,all.x=TRUE, by=c('src','trgt'))
           table(d.new2$outcomes,d.new2$is.causal1)
         }
         
         #removing based on inequality ?? TOOD correct
         if(FALSE){
           setDT(pred_outcome)
           data.table::setorderv(pred_outcome,cols='predictions',order = -1)
           pred_outcome.top2c = pred_outcome[, head(.SD, 2), by = "trgt"]
           pred_outcome.top2c=pred_outcome.top2c[,`:=`(predictions=NULL,outcomes=NULL,rvpred=NULL , is.causal1 =1)]
           d.new2 = merge(d.new1,pred_outcome.top2c,all.x=TRUE, by=c('src','trgt'))
           table(d.new2$outcomes,d.new2$is.causal1)
           
           
           pred_outcome.top1c = pred_outcome[, head(.SD, 1), by = "trgt"]
           pred_outcome.top1c=pred_outcome.top1c[,`:=`(predictions=NULL,outcomes=NULL,rvpred=NULL , is.causal1 =1)]
           
           d.new2 = merge(d.new1,pred_outcome.top1c,all.x=TRUE, by=c('src','trgt'))
           table(d.new2$outcomes,d.new2$is.causal1)
         }
         #==>removing based on inequality ?? TOOD correct
         if(FALSE){
           setDT(pred_outcome)
           tmp.inps = pred_outcome[,  sum(.SD$predictions), by = "trgt"]%>% dplyr::rename(prdInp = V1)
           tmp.outps = pred_outcome[,  sum(.SD$predictions), by = "src"]%>% dplyr::rename(prdOutp = V1)
           pred_outcome.nodetypes = full_join(tmp.inps,tmp.outps,by=c('trgt'='src'))
           pred_outcome.nodetypes = pred_outcome.nodetypes %>% dplyr::rename(nodename = trgt)
           mutate(prdtype = ifelse(prdInp <prdoutps,"effect","cause"  ))
           
           pred.outcome1= pred_outcome %>% left_join(pred_outcome.nodetypes,c('src'='nodename')) %>% left_join(pred_outcome.nodetypes,c('trgt'='nodename'))   
           pred.outcome2 = pred.outcome1 %>% dplyr::filter(is.na(rvpred) | (!is.na(rvpred) & prdOutp.x>prdInp.x & prdInp.y>prdOutp.y))
           
           pred.outcome2 = pred.outcome1 %>% dplyr::mutate(predictions=ifelse(!is.na(rvpred) & prdOutp.x<prdInp.x ,0,predictions))
           
           
           d.new2 = merge(d.new1,pred_outcome.top1c,all.x=TRUE, by=c('src','trgt'))
           table(d.new2$outcomes,d.new2$is.causal1)
         }
         
       }
       
       #introducing error to results for further testing
       if(FALSE){
         errprcnt = 0.1
         rplcIndices = as.integer(runif(length(predictions)*errprcnt, min = 0, max = length(predictions)))
         predictions[rplcIndices] <- runif(length(rplcIndices),min=0,max=1)
       }
       
       
       #TODO remove the reverse edges before assessment, just keep the one with higher prediction score
       theROC <- roc(assespreds$outcomes, assespreds$predictions, percent = TRUE);theROC
       evals$roc1 = theROC
       
       sscurves <- evalmod(scores = assespreds$predictions, labels = assespreds$outcomes)
       pr.crv=sscurves$prcs[1][[1]];  pr.auc =  attr(pr.crv,'auc')
       evals$pr.crv=pr.crv;evals$pr.auc=pr.auc
       
       theauc = precrec::auc(sscurves); 
       msg = paste0(theauc[[3]],"=", round(theauc[[4]],3));cat(msg)
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
                     algo,  round(pr.auc/rnd.auc,2), round(pr.prtauc/rnd.prtauc,2))
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
                     round(pr.auc/rnd.auc,2), round(pr.prtauc/rnd.prtauc,2))
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
  url.allresults  = paste0(url.output,dataFolder,'/',db, '/', 'evaluation.csv')
  f.all = list.files(path = '/scratch/as15096/eric/outputs',pattern = 'rankedEdges.csv', recursive = T)
  
  for(f.r in f.all){
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
    extractResourcesInf <- function(fldr){
      out=data.frame()
      f.all = list.files(path = fldr,pattern = 'time_.*txt', recursive = F)
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
    
    
    extractNestedInf <- function(p,prefix = 'unsn',type=''){
      
      out = data.frame()
      if(is.null(p)) return(out)
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
          p = r$cm[1:2,1:2]
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
      
      d = d %>%
        cbind(extractNestedInf(r$cm,'cm','cm' )) %>%
        cbind(extractNestedInf(r$unseensmpl_stats,'unsn','stats' )) %>%
        cbind(extractNestedInf(r$unseensmpl_part,'unsn','auc_table' )) %>%
        cbind(extractNestedInf(r$unseensmpl_rndm,'rndm','auc_table' )) %>%
        cbind(extractNestedInf(r$max_criteria,'crt','max_criteria' ))
      
      d
    }
  }
  
  analysisFolder = "sens_multipleRuns"
  path.base = "/scratch/as15096/eric/outputs/cict_par/"
  r.path = paste0(path.base,analysisFolder,"/")
  
  cnf.files = list.files(r.path,pattern = '.*[.]yaml')
  r.evals = data.frame()
  for(r.file in cnf.files){
    cnf = read_yaml(paste0(r.path, r.file))
    if(length(cnfs)>1){
      cnf.d=data.frame(  evaluation = cnf$experiment,
                         parConfFile = cnf$parConf,
                         job_idx = as.integer(cnf$given_job_id),
                         expLevel = cnf$expLevel,
                         dataset = cnf$datasets[[1]]$name,
                         edgeType = cnf$CICT$edgeType,
                         size.groundTruth = cnf$CICT$maxGroundTruth,
                         size.randomEdgeFoldCausal = cnf$CICT$randomEdgesFoldCausal,
                         trainingTarget = 'class2')
                         
      r.evals = rbind(r.evals,cnf.d)
      
    }
  }
  
  
   #  "sens_edgeType/"
  r.files = list.files(r.path,pattern = '.*[.]rds')
  r.all =  data.frame()
  for(r.file in r.files){
    print(r.file)
    r=data.frame()
    r = readRDS(paste0(r.path, r.file))
    if(length(r)>1){
      r.flat =  flattenResults(r)  #%>%as.data.frame()
      r.flat =r.flat %>% as.data.frame() %>% 
        mutate(job_idx = str_extract(r.file,pattern = '(?<=[_]).*(?=[.]rds)') %>%as.integer())
      r.all = rbind(r.all,r.flat)
      
    }
  }
  
  #Convert data.frame of list columns to standard data.frame
  lst <- lapply(r.all, unlist)
  r.all =data.frame(lapply(lst, `length<-`, median(lengths(lst))))
  
  glimpse(r.all)
  # r.all = sapply(r.all, function(x) unlist(x) ) %>% as.data.frame() %>% 
  #   group_by(rslt) %>%slice_head(n=1) %>%ungroup()

  r.resources = extractResourcesInf(r.path)
  
  r.final = r.evals %>% inner_join(r.all, by='job_idx') %>% left_join(r.resources,by = 'job_idx')
  
  write.csv(r.final, file = paste0(r.path,analysisFolder, ".csv"))

}

#r.exp = prd3  %>% left_join(r.all,by= c('rslt'='rslt' )
