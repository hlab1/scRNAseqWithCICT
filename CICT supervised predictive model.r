########################################################################
# Abbas Shojaee, development: 04-07-2019 , last modification  07-20-2021
########################################################################

#   TRAINING H2O  For Prediction and feature selection ----

require(h2o)
h2o.shutdown(prompt=FALSE)
gc()
H2OCnn = h2o.init(nthreads = parallel::detectCores()-4, enable_assertions = TRUE,max_mem_size = "16g")
#memory.limit(size=25000)
#Create binary outcome to calculate probaibility
#Pretraining with tst2
{
  tst2.train.h2o = as.h2o(tst4.totalset)
  tst2.rf = h2o.randomForest(mdlColNames,"class2",tst2.train.h2o,#ntrees = 1)
                             ntrees=8,nfolds=3) #,max_depth=6)
  tst2.mdl = tst2.rf
  
  
  
  tst2.gbm = h2o.gbm(mdlColNames,"class1",tst2.train.h2o,ntrees=20,nfolds = 10 )
  tst2.mdl = tst2.gbm
  
  plot(tst2.mdl)
  print(tst2.mdl) #save(tst1,tst1.train,tst1.tst, file="temp training set.robj")
  View(h2o.varimp(tst2.mdl))
  
  tst3.train.h2o = as.h2o(tst3.totalset)
  tst3.rf = h2o.randomForest(mdlColNames,"class1",tst3.train.h2o,ntrees=5) #,max_depth=6)
  tst3.mdl = tst3.rf
  
  tst3.gbm = h2o.gbm(mdlColNames,"class13",tst3.train.h2o,nfolds = 10 )
  tst3.mdl = tst3.gbm
  
  plot(tst3.mdl)
  print(tst3.mdl) #save(tst1,tst1.train,tst1.tst, file="temp training set.robj")
  View(h2o.varimp(tst3.mdl))
}

{ 
  #names(tst2.totalset) [!names(tst2.totalset) %in% names(tst1.train)]
# names(tst1.train) [!names(tst1.train) %in% names(tst2.totalset)]
#tst1.train= tst1.train
#tst1.train$class11 = ifelse(tst1.train$class1== '1'  , 1,0)
#tst1.train$class12 = ifelse(tst1.train$class1== '1' | tst1.train$class1== '2'  , 1,tst1.train$class1)
}
trainingTarget;table(tst1.totalset[,get(trainingTarget)])

NFolds = 5
table(setDT(tst1.totalset)[,get(trainingTarget)])
filterOutFeatures =  c('srctrgtSum','srctrgtProduct') #Weight
#existingmdlColNames = intersect(colnames(t2),mdlColNames)

mdlColNames=setdiff(mdlColNames,filterOutFeatures)
d.new = sample_frac(tst1.tst,size = .2)
tst1.tst.h2o=as.h2o(d.new) # as.h2o(tst1.tst) 
tst1.train.h2o = as.h2o(tst1.train)#totalset)#train)#
#tst1.totalset.h2o = as.h2o(tst1.totalset)
# multiclass categorial
rm(tst1.rf,tst1.mdl)

if(FALSE){
  #sort(names(tst1.totalset));sort(mdlColNames)
  
  tst1.rf = h2o.randomForest(mdlColNames,trainingTarget,tst1.train.h2o,nfolds = NFolds,
                             keep_cross_validation_predictions= TRUE,  validation_frame = tst1.tst.h2o)
  #,ntrees=30,max_depth=6
  #selectedFeatures #mdlColNames
  tst1.mdl=tst1.rf
  h2o::h2o.saveModel(tst1.mdl,path = "hesc_all_high_variations.robj")
}
#optimize for causal-versus with ntrees=10,max_depth=10, training=75%
#OPTIMZE FOR dream4ds2 causal versus random  max_depth =7,ntrees=8,
  tst1.rf = h2o.randomForest(mdlColNames,trainingTarget,
                             tst1.train.h2o,nfolds = NFolds,#checkpoint=tst2.mdl@model_id,
                             max_depth =9,ntrees=12, 
                             keep_cross_validation_predictions= FALSE,  validation_frame = tst1.tst.h2o)
  #,ntrees=30,max_depth=5#selectedFeatures #mdlColNames
  tst1.mdl=tst1.rf

# tst1.gbm = h2o.gbm(mdlColNames,trainingTarget,tst1.train.h2o,nfolds = NFolds,ntree=30,max_depth=10,  validation_frame = tst1.tst.h2o )
# #,ntrees=30,max_depth=5
# #selectedFeatures #mdlColNames
# #,validation_frame = tst1.tst.h2o)
# #,checkpoint="DRF_model_R_1467413037220_19")
# tst1.mdl=tst1.gbm

save(tst1.mdl,t1,t2,tst1.totalset,tst1.train,tst1.tst,file="mdl1.robj")
#Testing performance
{
  plot(tst1.mdl);setDF(tst1.totalset)
  print(tst1.mdl) #save(tst1,tst1.train,tst1.tst, file="temp training set.robj")
  prd.varimp=h2o.varimp(tst1.mdl);View(prd.varimp)
  paste0(prd.varimp$variable[1:30],collapse= "','")
  prd.varimp[1:15,] %>% knitr::kable()
  #h2o.saveModel(tst1.mdl,"mdlD100-1") #h2o.loadModel("E:/HDDs/Work/0 Plant biology/R code/mdlD100-1")
  #library(broom); tidy(tst1.mdl); glance(tst1.mdl)
  
  #plot important features distributions
  {
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
                               line = list(color = ~class,  #TODO set to prediction?
                                           colorscale = 'Jet',
                                           showscale = TRUE,
                                           reversescale = TRUE,
                                           cmin = -4000,
                                           cmax = -100),
                               dimensions = list(
                                 list(range = c(~min(blockHeight),~max(blockHeight)),
                                      constraintrange = c(100000,150000),
                                      label = 'Block Height', values = ~blockHeight),
                                 list(range = c(~min(blockWidth),~max(blockWidth)),
                                      label = 'Block Width', values = ~blockWidth),
                                 list(tickvals = c(0,0.5,1,2,3),
                                      ticktext = c('A','AB','B','Y','Z'),
                                      label = 'Cyclinder Material', values = ~cycMaterial),
                                 list(range = c(-1,4),
                                      tickvals = c(0,1,2,3),
                                      label = 'Block Material', values = ~blockMaterial),
                                 list(range = c(~min(totalWeight),~max(totalWeight)),
                                      visible = TRUE,
                                      label = 'Total Weight', values = ~totalWeight),
                                 list(range = c(~min(assemblyPW),~max(assemblyPW)),
                                      label = 'Assembly Penalty Weight', values = ~assemblyPW),
                                 list(range = c(~min(HstW),~max(HstW)),
                                      label = 'Height st Width', values = ~HstW),
                                 list(range = c(~min(minHW),~max(minHW)),
                                      label = 'Min Height Width', values = ~minHW),
                                 list(range = c(~min(minWD),~max(minWD)),
                                      label = 'Min Width Diameter', values = ~minWD),
                                 list(range = c(~min(rfBlock),~max(rfBlock)),
                                      label = 'RF Block', values = ~rfBlock)
                               )
      )
      
      
      fig
    }
    
    
    #Violin
    #%%%%%%%%%%%%%%%%%%  Distribution of edges top predictors
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
          label = dfl[class1 == 'c',]$name,
          color = dfl$class1,
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
      myMean()
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
  
  
  setDF(tst1.totalset)
  
  require(verification)
  require(pROC)
  require(ROCR)
  require( OptimalCutpoints)
  require(precrec )
  
  
  #newDataEnv = new.env()
  #R.utils::loadToEnv(file="CICT dream4-100-2.robj",envir=newDataEnv)
  
  #Predictalldata
  if(FALSE){
    setkey(t2,src,trgt)
    setkey(n.itm.e,src,trgt)
    d.new = n.itm.e[!t2]
    d.new[,(trainingTarget) := ifelse(is.na(get(trainingTarget)),0,get(trainingTarget))]
  }
  
  #d.new =  sample_frac(tst1.tst,size =.5)  #tst1.totalset #newDataEnv$tst1.totalset # 
  predTest.d.h2o = tst1.tst.h2o #as.h2o(d.new) # 
  h2o.pred = as.data.frame(h2o.predict(tst1.mdl,predTest.d.h2o,keep_cross_validation_predictions=TRUE))
  h20.prediction=as.numeric(as.character(h2o.pred[,3]))
  predictions =h20.prediction #pred #ens.predictions #
  outcomes =unlist(setDT(d.new)[,trainingTarget,with=FALSE]) #outcome # ens.outcome# 
  rndPred = ifelse(runif(nrow(d.new)) >=.5,1,0)
  
  d.new1 = as.data.frame( cbind(d.new[,.(src,trgt,Weight)],predictions,outcomes,rndPred))
  #d.new1 = d.new1[!duplicated(d.new1),]
  d.new1.rv =d.new1 %>% dplyr::rename(rvpred=predictions,src1=src,trgt1=trgt) %>% dplyr::select(-outcomes,-rndPred)
  pred_outcome = d.new1 %>% left_join(d.new1.rv, by=c("src"="trgt1", "trgt"="src1") )
  
  #post processing filtering
  pred_outcome = pred_outcome %>% filter(src %in% tf.all)
  
  
  
  
  View(pred_outcome[predictions > .7 | outcomes ==TRUE,])
  # Pruning
  {
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
    #removing based on inequality ?? TOOD correct
    if(TRUE){
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
  assespreds = pred_outcome
  roc <- roc(assespreds$outcomes, assespreds$predictions, percent = TRUE);roc
  
  sscurves <- evalmod(scores = assespreds$predictions, labels = assespreds$outcomes)
  pr.crv=sscurves$prcs[1][[1]];  pr.auc =  attr(pr.crv,'auc')
  
  theauc = precrec::auc(sscurves); paste0(theauc[[3]],"=", round(theauc[[4]],3))
  theauc = precrec::au(sscurves); paste0(theauc[[3]],"=", round(theauc[[4]],3))
  autoplot(sscurves, mode = 'rocpr',ret_grob = F) 
  
  #partial precision-recall
  sscurves.part <- part(sscurves, xlim = c(0, 0.2))
  plot(sscurves.part)
  sscurves.part
  pr.prtcrv=sscurves.part$prcs[1][[1]];  pr.prtauc =  attr(pr.prtcrv,'pauc')
  
  
  #random classifier
  print("Random Classifier comparison ============")
  randomClassifierCurves <- evalmod(scores = assespreds$rndPred, labels = assespreds$outcomes)
  rnd.crv=randomClassifierCurves$prcs[1][[1]]; rnd.auc =  attr(rnd.crv,'auc')
  
  rndmClscurves.part <- part(randomClassifierCurves, xlim = c(0, 0.2))
  rndmClscurves.part
  rnd.prtcrv=rndmClscurves.part$prcs[1][[1]]; rnd.prtauc =  attr(rnd.prtcrv,'pauc')
  
  print(sprintf("AUCPR Ratio CICT to Random: %s ", pr.auc/rnd.auc))
  
  print(sprintf("Partial AUCPR Ratio CICT to Random: %s ", pr.prtauc/rnd.prtauc))
  
  mmpoins <- evalmod(scores = assespreds$predictions, labels = assespreds$outcomes, mode = "basic")
  #autoplot(mmpoins, c("error", "accuracy"))
  
  # Show normalized ranks vs. specificity, sensitivity, and precision
  autoplot(mmpoins, c("specificity", "sensitivity", "precision"))
  
  # Show normalized ranks vs. Matthews correlation coefficient and F-score
  autoplot(mmpoins, c("mcc", "fscore"))
  
  #bestcutoff = coords(roc, "best", ret="threshold", transpose = FALSE)
  
  # the relative cost of of a false negative classification (as compared with a false positive classification)
  # the prevalence, or the proportion of cases in the population (n.cases/(n.controls+n.cases)).
  relativeCostfn_fp = 1/10
  prv = table(setDT(tst1.totalset)[,get(trainingTarget)])
  best.weights=c(relativeCostfn_fp, prv[2]/prv[1])  
  bestcutoff =as.double(coords(roc, "best", best.method = "closest.topleft", best.weights=best.weights,
                               ret="threshold", transpose = FALSE));bestcutoff
  #bestcutoff = 0.5
  assespreds$hresholdpreds= ifelse(assespreds$predictions>bestcutoff,TRUE,FALSE)
  pander::pander(ftable(factor(assespreds$hresholdpreds, c("TRUE","FALSE")), 
                        factor(assespreds$outcomes, c("TRUE","FALSE")), 
                        dnn=c('pred','actual'))) #%>%  knitr::kable()
  #table(ifelse(predictions>.99,'c','other'))
  #bestcutoff = coords(roc, "best", best.method = "closest.topleft", ret="threshold", transpose = FALSE);bestcutoff
  print(paste0("Best threshold: " ,bestcutoff ))
  
}
