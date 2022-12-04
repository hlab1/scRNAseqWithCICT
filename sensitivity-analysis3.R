#Functions and libs
{
  library('forestploter')
  library(ggplot2)
  library(ggparty)
  library(tidyverse)
  library(data.table)
  library(cowplot)
  
  myplotViolin<-function(theTitle,dt,variables, group1, group2,
                         grp1title,grp2title,grpcolor1 = "green",grpcolor2 = "blue", 
                         yval = 'value',thebox=F,...)
  {
    #dataframe in long format
    dt=setDT(dt)[name %in% variables,]
    
    fig <- setDT(dt) %>%
      plot_ly(type = 'violin',...) 
    
    fig <- fig %>%
      add_trace(
        x = ~dt[class10,get('name')],
        y = ~dt[class10,get(yval)],
        legendgroup = grp1title,
        scalegroup = grp1title,
        name = grp1title,
        side = 'negative',
        box = list(
          visible = thebox
        ),
        meanline = list(
          visible = thebox
        ),
        color = I(grpcolor1)
      ) 
    fig <- fig %>%
      add_trace(
        x = ~dt[class11,get('name')],
        y = ~dt[class11,get(yval)],
        legendgroup = grp2title,
        scalegroup = grp2title,
        name = grp2title,
        side = 'positive',
        box = list(
          visible = thebox
        ),
        meanline = list(
          visible = thebox
        ),
        color = I(grpcolor2)
      ) 
    
    fig <- fig %>%
      layout(
        title=theTitle,
        xaxis = list(
          title = ""  
        ),
        yaxis = list(
          title = "",
          zeroline = F
        ),
        violingap = 0,
        violingroupgap = 0,
        violinmode = 'overlay'
      )
    
    
    fig = fig %>% config(displayModeBar = F, showTips = F)
    fig
  }
  
  
  myhist = function(v,nbins,breaks=NULL,plot=F,prob = T)
  {
    h=""
    if(!is.null(breaks)) h = hist( v,breaks = breaks,plot=F,probability = prob)  else h = hist( v,nbins,plot=F,probability = prob)
    h=c(h,rep(0,nbins-length(h)))
    if(prob ==T) h$density else h$counts
  }
  
  myMean <- function (x,default = 0, na.rm = TRUE,...)
  {
    result = mean(x,na.rm=na.rm,...)
    if(is.na(result)) default else result
  }
  
  myMedian <- function (x,default = 0, na.rm = TRUE,...)
  {
    result = median(x,na.rm=na.rm,...)
    if(is.na(result)) default else result
  }
  
  mySD <- function (x,default = 0, na.rm = TRUE)
  {
    result = sd(x,na.rm)
    if(is.na(result)) default else result
  }

}

#Forest plot
{
  rm(r.others.dt,r.multirun.dt,r.edgeType.dt, r.sparsity.dt)
  desiredOutcome = 'unsn_pPRC' #   'unsn_std_pPRC' # 'unsn_PRC' # 
  desiredRandomOutcome = 'rndm_pPRC' # 'rndm_std_pPRC'
  PSEUDO_Zero.std_pPRC=0.00001 #Maxmimum observed value for the network
  basePath = 'E:/HDDs/Work/NYU Research/0 Plant biology/CICT paper/SensitivityRound2/'
  #r.multirun
  {
    r.multirun = read.csv(paste0(basePath,'sens_multipleRuns.csv'))
    glimpse(r.multirun); table(r.multirun$edgeType)
    
    r.multirun.dt.raw = r.multirun %>% 
      mutate(
        edges.total=as.numeric(str_replace_all(edges.total,',','')),
        outcome = .data[[desiredOutcome]] ,
        subgroup =paste0("    " ,  dataset)) %>%
      filter(!is.na(outcome) & !is.na(density.net)) %>%
      select(outcome, edges.total,everything())
    
      # filter(edgeType == 'Pearson') %>%
    r.multirun.dt=r.multirun.dt.raw %>% 
      group_by(subgroup)  %>%
      #mutate(outcome = unsn_pROC/rndm_pPRC) %>%
      mutate(rRatio= .data[[desiredOutcome]]/
                ifelse(.data[[desiredRandomOutcome]]==0,
                       max(PSEUDO_Zero.std_pPRC,min(.data[[desiredRandomOutcome]],na.rm=T)),
                       .data[[desiredRandomOutcome]])) %>%
      
      summarize(outcome.se = sd(outcome,na.rm=T) / sqrt(n()),
                est = mean(outcome,na.rm = T),
                low = est - 1.96*outcome.se, 
                hi =  est + 1.96*outcome.se,
                `Standardized.partial.AUPR` = paste0(rep(' ',60),collapse=''),
                outcome.ci = sprintf("%.2f (%.2f to %.2f)", est, low, hi),
                
                rRatio.se =  sd(rRatio,na.rm=T) / sqrt(n()),
                rRatio.est = mean(rRatio,na.rm = T),
                rRatio.low =  rRatio.est - 1.96*rRatio.se, 
                rRatio.hi = rRatio.est + 1.96*rRatio.se,
                `rRatio.Standardized.partial.AUPR` = paste0(rep(' ',60),collapse=''),
                rRatio.outcome.ci = sprintf("%.2f (%.2f to %.2f)",rRatio.est, rRatio.low,rRatio.hi),
                
                evalCount = n(),
                description= sprintf("N edges = %0.3g | density= %.3f" ,
                                     mean(edges.total,na.rm=T), mean(density.net,na.rm=T))) %>%
      arrange(desc(est))
  
    
    r.titlerow = data.frame(subgroup = 'L2 dataset, 100 runs, ewMIshrink', est='',low='',hi='',se = '',
                            outcome.ci=' ',`Standardized.partial.AUPR`= ' ',`Samples(N)`='', description=' ')
    r.multirun.dt = rbindlist(list(r.titlerow,r.multirun.dt),fill=T)
    
    rm(r.multirun.dt.raw)
  }
  
  #r.edgeType
  {
    r.edgeType = read.csv(paste0(basePath,'sens_edgeType.csv'))
    glimpse(r.edgeType)
    
    r.edgeType.dt.raw = r.edgeType %>% 
      mutate(#outcome = unsn_pROC/rndm_pPRC,
             outcome = .data[[desiredOutcome]],
             subgroup =paste0("    " , edgeType )) %>% #size.groundTruth  cictRawEdgeCol
      filter( !is.na(outcome) &
              dataset != 'mDC' &  #Very low density, outlier
              str_detect(dataset.dir , 'L2$')) 
      #filter(str_detect(dataset.dir,'lofgof|_ns',negate = T))%>%
      
    r.edgeType.dt = r.edgeType.dt.raw %>%
      group_by(subgroup)  %>%
      mutate(rRatio= .data[[desiredOutcome]]/
               ifelse(.data[[desiredRandomOutcome]]==0,
                      max(PSEUDO_Zero.std_pPRC,min(.data[[desiredRandomOutcome]],na.rm=T)),
                      .data[[desiredRandomOutcome]])) %>%
      
      summarize(outcome.se = sd(outcome,na.rm=T) / sqrt(n()),
                est = mean(outcome,na.rm = T),
                low = est - 1.96*outcome.se, 
                hi =  est + 1.96*outcome.se,
                `Standardized.partial.AUPR` = paste0(rep(' ',60),collapse=''),
                outcome.ci = sprintf("%.2f (%.2f to %.2f)", est, low, hi),
                
                rRatio.se =  sd(rRatio,na.rm=T) / sqrt(n()),
                rRatio.est = mean(rRatio,na.rm = T),
                rRatio.low =  rRatio.est - 1.96*rRatio.se, 
                rRatio.hi = rRatio.est + 1.96*rRatio.se,
                `rRatio.Standardized.partial.AUPR` = paste0(rep(' ',60),collapse=''),
                rRatio.outcome.ci = sprintf("%.2f (%.2f to %.2f)",rRatio.est, rRatio.low,rRatio.hi),
                
                evalCount = n(),
                description='') %>%
      filter(!is.na(est)) %>% arrange(desc(est))
    
    
    r.titlerow = data.frame(subgroup = 'L2 cell specific ChIPseq, directed edges', est=NA,low=NA,hi=NA,se = NA,
                            outcome.ci=' ',`Standardized.partial.AUPR`= ' ',`Samples(N)`='', description=' ')
    r.edgeType.dt = rbindlist(list(r.titlerow,r.edgeType.dt),fill=T)
    
    
    
    
    table(r.edgeType$dataset.dir,r.edgeType$size.groundTruth,r.edgeType$trainingTarget)
    
    if(F) #Testing
    {
      r.edgeType.dt = r.edgeType %>% 
        mutate(subgroup =paste0("    " , dataset )) %>% #size.groundTruth  cictRawEdgeCol
        filter(str_detect(dataset.dir , 'L2')) %>%
        filter(str_detect(edgeType , 'ewMImm')) %>%
        #filter(str_detect(dataset.dir,'lofgof|_ns',negate = T))%>%
        group_by(subgroup)  %>%
        #mutate(outcome = unsn_pROC/rndm_pPRC) %>%
        mutate(outcome = .data[[desiredOutcome]]) %>%
        mutate(rRatio= .data[[desiredOutcome]]/
                 ifelse(.data[[desiredRandomOutcome]]==0,
                        max(PSEUDO_Zero.std_pPRC,min(.data[[desiredRandomOutcome]],na.rm=T)),
                        .data[[desiredRandomOutcome]])) %>%
        
        summarize(outcome.se = sd(outcome,na.rm=T) / sqrt(n()),
                  est = mean(outcome,na.rm = T),
                  low = est - 1.96*outcome.se, 
                  hi =  est + 1.96*outcome.se,
                  `Standardized.partial.AUPR` = paste0(rep(' ',60),collapse=''),
                  outcome.ci = sprintf("%.2f (%.2f to %.2f)", est, low, hi),
                  
                  rRatio.se =  sd(rRatio,na.rm=T) / sqrt(n()),
                  rRatio.est = mean(rRatio,na.rm = T),
                  rRatio.low =  rRatio.est - 1.96*rRatio.se, 
                  rRatio.hi = rRatio.est + 1.96*rRatio.se,
                  `rRatio.Standardized.partial.AUPR` = paste0(rep(' ',60),collapse=''),
                  rRatio.outcome.ci = sprintf("%.2f (%.2f to %.2f)",rRatio.est, rRatio.low,rRatio.hi),
                  
                  evalCount = n(),
                  description= 'L2') %>%
        arrange(desc(est)) %>% filter(!is.na(est))
    }
      
  
    }
  
  #r.gtType ChIPseq specifici, ns and goflof
  {
    r.gtType = read.csv(paste0(basePath,'sens_edgeType.csv'))
    glimpse(r.gtType)
    
    r.gtType.dt.raw = r.gtType %>% 
      mutate(#outcome = unsn_pROC/rndm_pPRC,
        edges.total=as.numeric(str_replace_all(edges.total,',','')),
        outcome = .data[[desiredOutcome]],
        subgroup =paste0("    " , dataset.dir )) %>% #size.groundTruth  cictRawEdgeCol
      filter( !is.na(outcome) &
                dataset != 'mDC' &  #Very low density, outlier
                str_detect(edgeType , 'Pearson|Kendall|ewMImm')) 
    #filter(str_detect(dataset.dir,'lofgof|_ns',negate = T))%>%
    
    r.gtType.dt = r.gtType.dt.raw %>%
      group_by(subgroup)  %>%
      mutate(rRatio= .data[[desiredOutcome]]/
               ifelse(.data[[desiredRandomOutcome]]==0,
                      max(PSEUDO_Zero.std_pPRC,min(.data[[desiredRandomOutcome]],na.rm=T)),
                      .data[[desiredRandomOutcome]])) %>%
      
      summarize(outcome.se = sd(outcome,na.rm=T) / sqrt(n()),
                est = mean(outcome,na.rm = T),
                low = est - 1.96*outcome.se, 
                hi =  est + 1.96*outcome.se,
                `Standardized.partial.AUPR` = paste0(rep(' ',60),collapse=''),
                outcome.ci = sprintf("%.2f (%.2f to %.2f)", est, low, hi),
                
                rRatio.se =  sd(rRatio,na.rm=T) / sqrt(n()),
                rRatio.est = mean(rRatio,na.rm = T),
                rRatio.low =  rRatio.est - 1.96*rRatio.se, 
                rRatio.hi = rRatio.est + 1.96*rRatio.se,
                `rRatio.Standardized.partial.AUPR` = paste0(rep(' ',60),collapse=''),
                rRatio.outcome.ci = sprintf("%.2f (%.2f to %.2f)",rRatio.est, rRatio.low,rRatio.hi),
                
                evalCount = n(),
                description = sprintf("N edges = %0.3g | density= %.3f" ,mean(edges.total,na.rm=T), mean(density.net,na.rm=T))) %>%
                filter(!is.na(est)) %>% arrange(desc(est))
    
    
    r.titlerow = data.frame(subgroup = 'Ground truth types, L2 data', est=NA,low=NA,hi=NA,se = NA,
                            outcome.ci=' ',`Standardized.partial.AUPR`= ' ', `Samples(N)`='', description=' ')
    r.gtType.dt = rbindlist(list(r.titlerow,r.gtType.dt),fill=T)
    
    
    
    
    table(r.gtType$dataset.dir,r.gtType$size.groundTruth,r.gtType$trainingTarget)
    
    if(F) #Testing
    {
      r.gtType.dt = r.gtType %>% 
        mutate(subgroup =paste0("    " , dataset )) %>% #size.groundTruth  cictRawEdgeCol
        filter(str_detect(dataset.dir , 'L2')) %>%
        filter(str_detect(gtType , 'ewMImm')) %>%
        #filter(str_detect(dataset.dir,'lofgof|_ns',negate = T))%>%
        group_by(subgroup)  %>%
        #mutate(outcome = unsn_pROC/rndm_pPRC) %>%
        mutate(outcome = .data[[desiredOutcome]]) %>%
        mutate(rRatio= .data[[desiredOutcome]]/
                 ifelse(.data[[desiredRandomOutcome]]==0,
                        max(PSEUDO_Zero.std_pPRC,min(.data[[desiredRandomOutcome]],na.rm=T)),
                        .data[[desiredRandomOutcome]])) %>%
        
        summarize(outcome.se = sd(outcome,na.rm=T) / sqrt(n()),
                  est = mean(outcome,na.rm = T),
                  low = est - 1.96*outcome.se, 
                  hi =  est + 1.96*outcome.se,
                  `Standardized.partial.AUPR` = paste0(rep(' ',60),collapse=''),
                  outcome.ci = sprintf("%.2f (%.2f to %.2f)", est, low, hi),
                  
                  rRatio.se =  sd(rRatio,na.rm=T) / sqrt(n()),
                  rRatio.est = mean(rRatio,na.rm = T),
                  rRatio.low =  rRatio.est - 1.96*rRatio.se, 
                  rRatio.hi = rRatio.est + 1.96*rRatio.se,
                  `rRatio.Standardized.partial.AUPR` = paste0(rep(' ',60),collapse=''),
                  rRatio.outcome.ci = sprintf("%.2f (%.2f to %.2f)",rRatio.est, rRatio.low,rRatio.hi),
                  
                  evalCount = n(),
                  description = 'L2') %>%
        filter(!is.na(est)) %>% arrange(desc(est))
    }
    
    
  }
  
  #r.modeling_choices
  {
    r.modeling_choices = read.csv(paste0(basePath,'sens_modeling_choices.csv'))
    r.modeling_choices.bestEdges = r.modeling_choices %>% 
      group_by(dataset.dir,edgeType)  %>%
      summarise(edgePerf = mean(unsn_pPRC,na.rm=T) ) %>% 
      group_by(dataset.dir) %>%
      arrange(desc(edgePerf)) %>% slice_head(n=3) #Pearson is the best overall
    
    #r.modeling_choices.gtsize
    {
        r.modeling_choices.gtsize.raw = r.modeling_choices %>% 
        mutate(outcome = .data[[desiredOutcome]],
               subgroup =paste0("    " ,  size.groundTruth)) %>% #size.groundTruth
        filter(
          !is.na(outcome) &
          dataset != 'mDC' &  #Very low density, outlier
          trainingTrgt == 'class2' &
          RF_ntree >=20,
          RF_max_depth >=10,
          #size.groundTruth == 500 &
          edgeType == 'Pearson') %>%
        select(outcome,everything())
      
      r.modeling_choices.gtsize = r.modeling_choices.gtsize.raw %>%
        group_by(subgroup)  %>%
        #mutate(outcome = unsn_pROC/rndm_pPRC) %>%
        mutate(rRatio= .data[[desiredOutcome]]/
                 ifelse(.data[[desiredRandomOutcome]]==0,
                        max(PSEUDO_Zero.std_pPRC,min(.data[[desiredRandomOutcome]],na.rm=T)),
                        .data[[desiredRandomOutcome]])) %>%
        
        summarize(outcome.se = sd(outcome,na.rm=T) / sqrt(n()),
                  est = mean(outcome,na.rm = T),
                  low = est - 1.96*outcome.se, 
                  hi =  est + 1.96*outcome.se,
                  `Standardized.partial.AUPR` = paste0(rep(' ',60),collapse=''),
                  outcome.ci = sprintf("%.2f (%.2f to %.2f)", est, low, hi),
                  
                  rRatio.se =  sd(rRatio,na.rm=T) / sqrt(n()),
                  rRatio.est = mean(rRatio,na.rm = T),
                  rRatio.low =  rRatio.est - 1.96*rRatio.se, 
                  rRatio.hi = rRatio.est + 1.96*rRatio.se,
                  `rRatio.Standardized.partial.AUPR` = paste0(rep(' ',60),collapse=''),
                  rRatio.outcome.ci = sprintf("%.2f (%.2f to %.2f)",rRatio.est, rRatio.low,rRatio.hi),
                  
                  evalCount = n(),
                  description = sprintf("N edges = %0.3g | density= %.3f" ,mean(edges.total,na.rm=T), mean(density.net,na.rm=T))) %>%
        filter(!is.na(est)) %>% arrange(desc(est))
      
      
      r.titlerow = data.frame(subgroup = 'N. TP edges in learning set', est=NA,low=NA,hi=NA,se = NA,
                              outcome.ci=' ',`Standardized.partial.AUPR`= ' ',`Samples(N)`='', description=' ')
      r.modeling_choices.gtsize = rbindlist(list(r.titlerow,r.modeling_choices.gtsize),fill=T)
      
      r.modeling_choices.gtsize.diag = rbind(
        r.modeling_choices.gtsize.raw %>% arrange(outcome) %>% slice_head(n=2),
        r.modeling_choices.gtsize.raw %>% arrange(outcome) %>% slice_tail(n=2)
      )
      rm(r.modeling_choices.gtsize.raw,r.modeling_choices.gtsize.diag)
      }
    
    #r.modeling_choices.trainingTrgt.raw
    {
      r.modeling_choices.trainingTrgt.raw = r.modeling_choices %>% 
        mutate(outcome = .data[[desiredOutcome]], #outcome = unsn_pROC/rndm_pPRC
               subgroup =paste0("    " ,  trainingTrgt)) %>% #size.groundTruth
        filter(
          !is.na(outcome) &
          dataset != 'mDC' &  #Very low density, outlier
            #trainingTrgt == 'class2' &
            RF_ntree >=20,
            RF_max_depth >=10,
            size.groundTruth == 500 &
            edgeType == 'Pearson') %>%
        select(outcome,everything())
      
      
      r.modeling_choices.trainingTrgt = r.modeling_choices.trainingTrgt.raw %>%
        group_by(subgroup)  %>%
        mutate(rRatio= .data[[desiredOutcome]]/
                 ifelse(.data[[desiredRandomOutcome]]==0,
                        max(PSEUDO_Zero.std_pPRC,min(.data[[desiredRandomOutcome]],na.rm=T)),
                        .data[[desiredRandomOutcome]])) %>%
        
        summarize(outcome.se = sd(outcome,na.rm=T) / sqrt(n()),
                  est = mean(outcome,na.rm = T),
                  low = est - 1.96*outcome.se, 
                  hi =  est + 1.96*outcome.se,
                  `Standardized.partial.AUPR` = paste0(rep(' ',60),collapse=''),
                  outcome.ci = sprintf("%.2f (%.2f to %.2f)", est, low, hi),
                  
                  rRatio.se =  sd(rRatio,na.rm=T) / sqrt(n()),
                  rRatio.est = mean(rRatio,na.rm = T),
                  rRatio.low =  rRatio.est - 1.96*rRatio.se, 
                  rRatio.hi = rRatio.est + 1.96*rRatio.se,
                  `rRatio.Standardized.partial.AUPR` = paste0(rep(' ',60),collapse=''),
                  rRatio.outcome.ci = sprintf("%.2f (%.2f to %.2f)",rRatio.est, rRatio.low,rRatio.hi),
                  
                  evalCount = n(),
                  description = '') %>%
        filter(!is.na(est)) %>% arrange(desc(est))
      
      r.modeling_choices.trainingTrgt$subgroup =   plyr::mapvalues(r.modeling_choices.trainingTrgt$subgroup,
                                                                   c('    class2','    class3'),c('    Directed','    Undirected'))
      
      r.titlerow = data.frame(subgroup = 'Training with directed vs. undirected edges', est=NA,low=NA,hi=NA,se = NA,
                              outcome.ci=' ',`Standardized.partial.AUPR`= ' ',`Samples(N)`='', description=' ')
      r.modeling_choices.trainingTrgt = rbindlist(list(r.titlerow,r.modeling_choices.trainingTrgt),fill=T)
      
      r.modeling_choices.trainingTrgt.diag = rbind(
        r.modeling_choices.trainingTrgt.raw %>% arrange(outcome) %>% slice_head(n=2),
        r.modeling_choices.trainingTrgt.raw %>% arrange(outcome) %>% slice_tail(n=2)
      )
      
      rm(r.modeling_choices.trainingTrgt.raw,r.modeling_choices.trainingTrgt.diag)
    }
  }
  
  #r.sparsity
  {
    r.sparsity = read.csv(paste0(basePath,'sens_sparsity.csv'))
    r.sparsity.bestEdges = r.sparsity %>% 
      group_by(dataset.dir,edgeType)  %>%
      summarise(edgePerf = mean(unsn_pPRC,na.rm=T) ) %>% 
      group_by(dataset.dir) %>%
      arrange(desc(edgePerf)) %>% slice_head(n=3) #Pearson is the best overall
    
    r.sparsity.dt = r.sparsity %>% 
      mutate(outcome = .data[[desiredOutcome]],
             subgroup =paste0("    " ,  dataset.dir)) %>% #size.groundTruth
      filter(! is.na(outcome) &
               edgeType %in% c( 'Kendall', 'Pearson')  &
               # RF_max_depth == 5 & 
               size.groundTruth==500) %>%
      group_by(subgroup)  %>%
      #mutate(outcome = unsn_pROC/rndm_pPRC) %>%
      mutate(rRatio= .data[[desiredOutcome]]/
               ifelse(.data[[desiredRandomOutcome]]==0,
                      max(PSEUDO_Zero.std_pPRC,min(.data[[desiredRandomOutcome]],na.rm=T)),
                      .data[[desiredRandomOutcome]])) %>%
      
      summarize(outcome.se = sd(outcome,na.rm=T) / sqrt(n()),
                est = mean(outcome,na.rm = T),
                low = est - 1.96*outcome.se, 
                hi =  est + 1.96*outcome.se,
                `Standardized.partial.AUPR` = paste0(rep(' ',60),collapse=''),
                outcome.ci = sprintf("%.2f (%.2f to %.2f)", est, low, hi),
                
                rRatio.se =  sd(rRatio,na.rm=T) / sqrt(n()),
                rRatio.est = mean(rRatio,na.rm = T),
                rRatio.low =  rRatio.est - 1.96*rRatio.se, 
                rRatio.hi = rRatio.est + 1.96*rRatio.se,
                `rRatio.Standardized.partial.AUPR` = paste0(rep(' ',60),collapse=''),
                rRatio.outcome.ci = sprintf("%.2f (%.2f to %.2f)",rRatio.est, rRatio.low,rRatio.hi),
                
                evalCount = n(),
                description = 'Missing rate = K%') %>%
      filter(!is.na(est)) %>% arrange(desc(est))
  
    
    r.titlerow = data.frame(subgroup = 'Sparsity', est=NA,low=NA,hi=NA,se = NA,
                            outcome.ci=' ',`Standardized.partial.AUPR`= ' ',`Samples(N)`='', description=' ')
    r.sparsity.dt = rbindlist(list(r.titlerow,r.sparsity.dt),fill=T)
    
  }
  
  #r.datasets
  if(F){
    r.edgeType = read.csv(paste0(basePath,'sens_edgeType_results.csv'))
    glimpse(r.edgeType)
    
    r.datasets.dt = r.edgeType %>% 
      filter(str_detect(datasetdir , 'L2')) %>%
      mutate(subgroup =paste0("    " ,  dataset)) %>% #size.groundTruth
      group_by(subgroup)  %>%
      mutate(rRatio= .data[[desiredOutcome]]/
               ifelse(.data[[desiredRandomOutcome]]==0,
                      max(PSEUDO_Zero.std_pPRC,min(.data[[desiredRandomOutcome]],na.rm=T)),
                      .data[[desiredRandomOutcome]])) %>%
      
      summarize(outcome.se = sd(outcome,na.rm=T) / sqrt(n()),
                est = mean(outcome,na.rm = T),
                low = est - 1.96*outcome.se, 
                hi =  est + 1.96*outcome.se,
                `Standardized.partial.AUPR` = paste0(rep(' ',60),collapse=''),
                outcome.ci = sprintf("%.2f (%.2f to %.2f)", est, low, hi),
                
                rRatio.se =  sd(rRatio,na.rm=T) / sqrt(n()),
                rRatio.est = mean(rRatio,na.rm = T),
                rRatio.low =  rRatio.est - 1.96*rRatio.se, 
                rRatio.hi = rRatio.est + 1.96*rRatio.se,
                `rRatio.Standardized.partial.AUPR` = paste0(rep(' ',60),collapse=''),
                rRatio.outcome.ci = sprintf("%.2f (%.2f to %.2f)",rRatio.est, rRatio.low,rRatio.hi),
                
                evalCount = n(),
                description = 'L2') %>%
      filter(!is.na(est)) %>% arrange(desc(est))
    
    r.titlerow = data.frame(subgroup = 'Datasets', est=NA,low=NA,hi=NA,se = NA,
                            outcome.ci=' ',`Standardized.partial.AUPR`= ' ',`Samples(N)`='', description=' ')
    r.datasets.dt = rbindlist(list(r.titlerow,r.datasets.dt),fill=T)
    
  }
  
  #r.levels
  if(F){
    r.edgeType = read.csv(paste0(basePath,'sens_edgeType_results.csv'))
    glimpse(r.edgeType)
    
    r.levels.dt = r.edgeType %>% 
      filter(str_detect(datasetdir , 'L2')) %>%
      mutate(subgroup =paste0("    " ,  str_extract(datasetdir , '^.*?(?=_)'))) %>% #size.groundTruth
      group_by(subgroup)  %>%
      mutate(outcome = unsn_std_pROC) %>%
      mutate(outcome.se = sd(outcome,na.rm=T) / sqrt(n())) %>%
      
      summarize(est = mean(outcome,na.rm = T),
                low =  min(outcome,na.rm = T), 
                hi = max(outcome,na.rm = T),
                se = (hi - est)/1.96,
                `Standardized.partial.AUPR` = paste0(rep(' ',60),collapse=''),
                outcome.ci = sprintf("%.2f (%.2f to %.2f)",
                                     mean(outcome,na.rm = T), min(outcome,na.rm = T), max(outcome,na.rm = T)),
                description = 'L2') %>%
      filter(!is.na(est))
    
    r.titlerow = data.frame(subgroup = 'Levels', est=NA,low=NA,hi=NA,se = NA,
                            outcome.ci=' ',`Standardized.partial.AUPR`= ' ',`Samples(N)`='', description=' ')
    r.levels.dt = rbindlist(list(r.titlerow,r.levels.dt),fill=T)
    
  }
  
  #r.others simulation
  if(F){
    
    r.others.dt = data.frame( subgroup = paste0('   ', c('RANDOM','GENIE3','GRISLI','GRNBOOST2','LEAP','PIDC',
                                                         'PPCOR','SCODE','SCRIBE','SINCERITIES','SINGE')),
                              est = rnorm(n=11,mean= .12,sd=.05),
                              low = rnorm(n=11,mean= .05,sd=.04), 
                              hi = rnorm(n=11,mean= .2,sd=.04)) %>%
      mutate(   se = (hi - est)/1.96,
                `Standardized.partial.AUPR` = paste0(rep(' ',60),collapse=''),
                outcome.ci = sprintf("%.2f (%.2f to %.2f)",
                                     mean(est,na.rm = T), min(est,na.rm = T), max(est,na.rm = T)),
                evalCount = n(),
                description = paste0('L2')) 
    
    r.titlerow = data.frame(subgroup = 'Benchmarked methods ', est=NA,low=NA,hi=NA,se = NA,
                            `Standardized.partial.AUPR`= ' ',outcome.ci=' ',`Samples(N)`='', description=' ')
    r.others.dt = rbindlist(list(r.titlerow,r.others.dt),fill=T,use.names = T)
    
  }
  
  
  #r.others real
  {
    
    r.others = readRDS(paste0(basePath,'otherMethodsPerformance_2.rds'))
    glimpse(r.others)
    
    r.others.dt.raw = r.others %>% 
      mutate(#outcome = unsn_pROC/rndm_pPRC,
        outcome = .data[[desiredOutcome]],
        subgroup =paste0("    " , algorithm )) %>% #size.groundTruth  cictRawEdgeCol
      filter( !is.na(outcome) &
                dataset != 'mDC' &  #Very low density, outlier
                #str_detect(algorithm , 'DEEP',negate = T) &
                str_detect(benchmark , '^L2$')) 
    #filter(str_detect(dataset.dir,'lofgof|_ns',negate = T))%>%
    
    r.others.dt = r.others.dt.raw %>%
      group_by(subgroup)  %>%
      #mutate(rRatio= unsn_pPRC/rndm_pPRC
      mutate(rRatio= .data[[desiredOutcome]]/.data[[desiredRandomOutcome]]
               # ifelse(.data[[desiredRandomOutcome]]==0,
               #        max(PSEUDO_Zero.std_pPRC,min(.data[[desiredRandomOutcome]],na.rm=T)),
               #        .data[[desiredRandomOutcome]])
             ) %>%
      
      summarize(outcome.se = sd(outcome,na.rm=T) / sqrt(n()),
                est = mean(outcome,na.rm = T),
                low = est - 1.96*outcome.se, 
                hi =  est + 1.96*outcome.se,
                `Standardized.partial.AUPR` = paste0(rep(' ',60),collapse=''),
                outcome.ci = sprintf("%.2f (%.2f to %.2f)", est, low, hi),
                
                rRatio.se =  sd(rRatio,na.rm=T) / sqrt(n()),
                rRatio.est = mean(rRatio,na.rm = T),
                rRatio.low =  rRatio.est - 1.96*rRatio.se, 
                rRatio.hi = rRatio.est + 1.96*rRatio.se,
                `rRatio.Standardized.partial.AUPR` = paste0(rep(' ',60),collapse=''),
                rRatio.outcome.ci = sprintf("%.2f (%.2f to %.2f)",rRatio.est, rRatio.low,rRatio.hi),
                
                evalCount = n(),
                description='') %>%
      filter(!is.na(est)) %>% arrange(desc(est))
    
    
    r.titlerow = data.frame(subgroup = 'Benchmarked algorithms, L2', est=NA,low=NA,hi=NA,se = NA,
                            outcome.ci=' ',`Standardized.partial.AUPR`= ' ',`Samples(N)`='', description=' ')
    r.others.dt = rbindlist(list(r.titlerow,r.others.dt),fill=T)
    
  }
  
  
  library(RColorBrewer)
  clrs = RColorBrewer::brewer.pal(8,'YlGnBu')[c(1,2,3,4,5,6,7)] #   'Pastel2')[c(5,6,7)]
  clrs = c(clrs,gray.colors(12)[9])
  
  library(paletteer)
  clrs = paletteer_c("ggthemes::Blue-Green Sequential", 7) 
  clrs=clrs[c(2:7,1)]
  
  
  #Desired measure Forest graph
  {
    dt = rbind(r.modeling_choices.gtsize
               ,r.modeling_choices.trainingTrgt
               ,r.gtType.dt           
               ,r.edgeType.dt
    
               ,r.multirun.dt
               ,r.sparsity.dt
               ,r.others.dt
               ) %>% 
      select(evalCount,everything()) %>%
      mutate(low = ifelse(is.na(low),0,low),
             hi = ifelse(is.na(hi),0,hi),
             evalCount = ifelse(is.na(evalCount),'',evalCount),
             outcome.semult100 = as.integer(outcome.se *200) )
    
    # dt= dt %>% mutate(across(everything(), ~ ifelse(is.na(.),'',.) )) | !str_detect(.$subgroup,'\t')
    # dt= dt %>% mutate(across(everything(), ~ ifelse(str_detect(.$subgroup,'\t'),'',.) ))
    
    #Forest graph
    {
      rm(p)
      p <- forest(select(dt,subgroup,outcome.ci ,`Standardized.partial.AUPR`, evalCount),
                  est = dt$est,
                  lower = dt$low, 
                  upper = dt$hi,
                  sizes =  2,#dt$outcome.semult100 ,#
                  ci_column = 3,
                  ref_line = 0,
                  #arrow_lab = c("Placebo Better", "Treatment Better"),
                  xlim = c(0, 1),
                  ticks_at = c(0.1,0.25, .5,.75,1), #    c(0, 0.05, .1, .15, .3),
                  footnote = "footnote",
                  title = paste0(" Outcome: " , desiredOutcome)
                ); 
      
      
      # Bold grouping text
      headingRows = which(str_detect(dt$subgroup,'^   ',negate = T))
      p <- edit_plot(p,
                     row = headingRows,
                     gp = gpar(fontface = "bold"))
      
      # Edit background of row 5
      headingRows = c(headingRows,nrow(dt)+1)
      subgrouprows = lapply(1:(length(headingRows)-1),
                            function(i)  headingRows[i]:(headingRows[i+1]-1)
        )
      
      
      p <- edit_plot(p, gp = gpar(fontfamily = 'Arial'));
      p <- edit_plot(p, gp = gpar(fontsize =12));
      
      p <- edit_plot(p, row = subgrouprows[[1]], which = "background",
                     gp = gpar(fill = clrs[1]));
      
      p <- edit_plot(p, row = subgrouprows[[2]], which = "background",
                     gp = gpar(fill = clrs[2]));
      
      p <- edit_plot(p, row = subgrouprows[[3]], which = "background",
                     gp = gpar(fill = clrs[3]));
      
      p <- edit_plot(p, row = subgrouprows[[4]], which = "background",
                     gp = gpar(fill = clrs[4]));
      
      p <- edit_plot(p, row = subgrouprows[[5]], which = "background",
                     gp = gpar(fill = clrs[5]));
      
      p <- edit_plot(p, row = subgrouprows[[6]], which = "background",
                     gp = gpar(fill = clrs[6]));
      
      # p <- edit_plot(p, row = subgrouprows[[7]], which = "background",
      #                gp = gpar(fill = clrs[7]));
      # 
      
      thepath = "E:/HDDs/Work/NYU Research/0 Plant biology/R code/results/Forest plot stdpAPRC2.tif"
      ggsave(filename = thepath, plot = p, width = 14, height = 14, device='tiff', dpi=300)
      p
    }
  }
  
  #Ratio of desired measure to random classifier Forest graph
  {
    dt = rbind(r.modeling_choices.gtsize
               ,r.modeling_choices.trainingTrgt
               ,r.gtType.dt           
               ,r.edgeType.dt
               
               ,r.multirun.dt
               ,r.sparsity.dt
               ,r.others.dt
    ) %>% 
      select(evalCount,everything()) %>%
      mutate(rRatio.low = ifelse(is.na(rRatio.low),0,rRatio.low),
             rRatio.hi = ifelse(is.na(rRatio.hi),0,rRatio.hi),
             evalCount = ifelse(is.na(evalCount),'',evalCount),
             rRatio.outcome.ci = ifelse(is.na(rRatio.outcome.ci),'',rRatio.outcome.ci),
             outcome.semult100 = as.integer(outcome.se *200) )
    
    
    #Forest graph
    {
      rm(p)
      p <- forest(select(dt,subgroup,rRatio.outcome.ci ,`Standardized.partial.AUPR`, evalCount), #description
                  est = dt$rRatio.est,
                  lower = dt$rRatio.low, 
                  upper = dt$rRatio.hi,
                  sizes = 3, #dt$rRatio.outcome.se ,#
                  ci_column = 3,
                  ref_line = 0,
                  #arrow_lab = c("Placebo Better", "Treatment Better"),
                  xlim = c(0, 180),
                  ticks_at = c(0,  50,100,150,200), #    c(0, 0.05, .1, .15, .3),
                  footnote = "footnote",
                  title = paste0(" Outcome: " , desiredOutcome)
      ); 
      
      
      # Bold grouping text
      headingRows = which(str_detect(dt$subgroup,'^   ',negate = T))
      p <- edit_plot(p,
                     row = headingRows,
                     gp = gpar(fontface = "bold"))
      
      # Edit background of row 5
      headingRows = c(headingRows,nrow(dt)+1)
      subgrouprows = lapply(1:(length(headingRows)-1),
                            function(i)  headingRows[i]:(headingRows[i+1]-1)
      )
      
      
      p <- edit_plot(p, gp = gpar(fontfamily = 'Arial'));
      p <- edit_plot(p, gp = gpar(fontsize =12));
      
      p <- edit_plot(p, row = subgrouprows[[1]], which = "background",
                     gp = gpar(fill = clrs[1]));
      
      p <- edit_plot(p, row = subgrouprows[[2]], which = "background",
                     gp = gpar(fill = clrs[2]));
      
      p <- edit_plot(p, row = subgrouprows[[3]], which = "background",
                     gp = gpar(fill = clrs[3]));
      
      p <- edit_plot(p, row = subgrouprows[[4]], which = "background",
                     gp = gpar(fill = clrs[4]));
      
      p <- edit_plot(p, row = subgrouprows[[5]], which = "background",
                     gp = gpar(fill = clrs[5]));
      
      p <- edit_plot(p, row = subgrouprows[[6]], which = "background",
                     gp = gpar(fill = clrs[6]));
      
      # p <- edit_plot(p, row = subgrouprows[[7]], which = "background",
      #                gp = gpar(fill = clrs[7]));
      # 
    }
    
    thepath = "E:/HDDs/Work/NYU Research/0 Plant biology/R code/results/Forest plot ratios2.tif"
    ggsave(filename = thepath, plot = p, width = 14, height = 13, device='tiff', dpi=300)
    p
  }
}


#Radar charts
{
  
  library(fmsb)
  library(ggradar)
  library(scales)
  
  r.meta =r.modeling_choices %>%
    rowwise() %>%     
    mutate(rRatio= .data[[desiredOutcome]]/
               ifelse(.data[[desiredRandomOutcome]]==0,
                      max(PSEUDO_Zero.std_pPRC,min(.data[[desiredRandomOutcome]],na.rm=T)),
                                  .data[[desiredRandomOutcome]])) %>%
    mutate(grp = paste0(RF_ntree,'_',RF_max_depth)) %>%
    group_by(size.groundTruth,grp) %>% 
    summarize(smrzRatio = mean(rRatio, na.rm=T))  %>%  ungroup() %>%
    pivot_wider(id_cols = size.groundTruth, names_from = grp,values_from = smrzRatio)


  r.meta =r.modeling_choices %>%
    rowwise() %>%     
    mutate(rRatio= .data[[desiredOutcome]]/
             ifelse(.data[[desiredRandomOutcome]]==0,
                    max(PSEUDO_Zero.std_pPRC,min(.data[[desiredRandomOutcome]],na.rm=T)),
                    .data[[desiredRandomOutcome]])) %>%
    mutate(grp = paste0(RF_ntree,'_',RF_max_depth)) %>%
    group_by(dataset,grp) %>% 
    summarize(smrzRatio = mean(rRatio, na.rm=T))  %>%  ungroup() %>%
    pivot_wider(id_cols = dataset, names_from = grp,values_from = smrzRatio)
  
  
  r.meta =r.modeling_choices %>%
    mutate(density.net.cleaned = ifelse(is.na(density.net), 0, density.net),
            density.disct = arules::discretize(density.net.cleaned,
                                               method = "frequency",breaks=5)) %>% 
      # density.disct = infotheo::discretize(density.net.cleaned, 
      #                                           disc="equalfreq", nbins=5)) %>% #n()^(1/3)
    rowwise() %>%     
    mutate( rRatio= .data[[desiredOutcome]]/
             ifelse(.data[[desiredRandomOutcome]]==0,
                    max(PSEUDO_Zero.std_pPRC,min(.data[[desiredRandomOutcome]],na.rm=T)),
                    .data[[desiredRandomOutcome]])) %>%
    mutate(grp = paste0('', RF_ntree,',  ',RF_max_depth,''),
           density.disct = paste0(density.disct)) %>%
    arrange(RF_ntree,RF_max_depth) %>%
    
    group_by(density.disct,grp) %>% 
    summarize(smrzRatio = mean(rRatio, na.rm=T))  %>% ungroup() %>%
    pivot_wider(id_cols = density.disct, names_from = grp,values_from = smrzRatio) %>%
    na.omit() #%>%column_to_rownames('density.disct') 
  
  
  maxval =r.meta %>% ungroup() %>% dplyr::select_if(is.numeric) %>% max(na.rm = T)
  ggradar::ggradar(r.meta,grid.max=maxval,
            group.line.width = 1, 
            group.point.size = 3,
            #group.colours = c("#00AFBB", "#E7B800", "#FC4E07"),
            # Background and grid lines
            background.circle.colour = "white",
            gridline.mid.colour = "grey",
            legend.position = "bottom"
    ) 
    
  ggradar(
    df, 
    values.radar = c("0", "10", "20"),
    grid.min = 0, grid.mid = 10, grid.max = 20,
    # Polygons
    group.line.width = 1, 
    group.point.size = 3,
    group.colours = c("#00AFBB", "#E7B800", "#FC4E07"),
    # Background and grid lines
    background.circle.colour = "white",
    gridline.mid.colour = "grey",
    legend.position = "bottom"
  )
}

##### Sensitivity analysis modeling**
if(F){
sensitivityAnalysis <- function(model,
                                target = NULL,
                                predictors = NULL, #Data Frame
                                predictorsMeans = NULL, #Data Frame
                                samplesN = 101,
                                from = 0.6,
                                to = 1.4,
                                targetPrediction = "ratio", # c("ratio", "absolute")
                                predictionType = "prediction",
                                level = 0.9) {
    
    if (missing(target)) {
        predictors <- model$model[-1]
        target <- model$model[1]
    }
    targetName <- names(target)
    
    target <- target[[1]]
    
    numberPredictors <- ncol(predictors)     
    
    if (missing(predictorsMeans)) {
        initialPredictorsValues <- predictors %>%
            summarise_all(mean)
        
    } else {
        initialPredictorsValues <- predictorsMeans    
    }
    
    initialTargetValue <- mean(predict(model, newdata = initialPredictorsValues))
    
    sensitivityData <- sapply(initialPredictorsValues, function(x) rep(x, samplesN * numberPredictors))
    changeDF <- seq(from, to, length.out = samplesN)
    
    for (i in seq(1, numberPredictors)) {
        sensitivityData[seq((i-1)*samplesN + 1, i * samplesN), i] <- changeDF * initialPredictorsValues[[i]]    
    }
    sensitivityData <- data.frame(sensitivityData)
    
    predictedTarget <- predict(model,
                               newdata = sensitivityData,
                               interval = predictionType,
                               level = level)
    
    predictedTarget <- as.data.frame(predictedTarget)
    
    if (targetPrediction == "ratio")
        predictedTarget <- mutate_all(predictedTarget, function(x) {x / initialTargetValue})
    
    df <- data.frame(normalized.predictor.change = numeric(0),
                     predictor = character(0),
                     normalized.target.change.fit = numeric(0),
                     normalized.target.change.lower= numeric(0),
                     normalized.target.change.upper = numeric(0)) 
    
    for (i in seq(1, numberPredictors)) {
        df <- rbind(df, 
                    data.frame(
                        normalized.predictor.change = changeDF,
                        predictor = names(predictors)[i],
                        normalized.target.change.fit = predictedTarget$fit[seq((i-1)*samplesN + 1, i * samplesN)],
                        normalized.target.change.lower = predictedTarget$lwr[seq((i-1)*samplesN + 1, i * samplesN)],
                        normalized.target.change.upper = predictedTarget$upr[seq((i-1)*samplesN + 1, i * samplesN)]))
    }
    if (targetPrediction == "ratio") {
        gg <- ggplot(df, aes(x = normalized.predictor.change, 
                             group = predictor)) +
            geom_ribbon(aes(ymin = normalized.target.change.lower,
                            ymax = normalized.target.change.upper,
                            fill = predictor), alpha = 0.1) + 
            geom_vline(xintercept = 1,  color = "grey80") + 
            geom_hline(yintercept = 1, color = "grey80") +
            geom_abline(slope = 1, linetype = "dashed", color = "grey80") + 
            geom_abline(slope = -1, intercept = 2, linetype = "dashed", color = "grey80") +
            geom_line(aes(x = normalized.predictor.change,
                          y = normalized.target.change.fit,
                          color = predictor)) +
            ylab(paste("Normalized", targetName, "change")) +
            xlab("Normalized Predictor Change") + 
            theme_few()
    } else {
        gg <- ggplot(df, aes(x = normalized.predictor.change, 
                             group = predictor)) +
            geom_ribbon(aes(ymin = normalized.target.change.lower,
                            ymax = normalized.target.change.upper,
                            fill = predictor), alpha = 0.1) + 
            geom_vline(xintercept = 1,  color = "grey80") + 
            geom_hline(yintercept = initialTargetValue, color = "grey80") +
            geom_line(aes(x = normalized.predictor.change,
                          y = normalized.target.change.fit,
                          color = predictor)) +
            ylab(targetName) + 
            xlab("Normalized Predictor Change") + 
            theme_few()        
    }
    
    variableMeans <- cbind(data.frame(target = initialTargetValue),
                           data.frame(initialPredictorsValues))
    names(variableMeans) <- c(targetName, names(initialPredictorsValues))
    return(list(ggplot = gg,
                predictionData = sensitivityData,
                resultsData = df,
                variableMeans = variableMeans))
}

#####################################
sensitivityAnalysisCaret <- function(model,
                                target = NULL,
                                predictors = NULL, #Data Frame
                                predictorsMeans = NULL, #Data Frame
                                samplesN = 101,
                                from = 0.6,
                                to = 1.4,
                                targetPrediction = "ratio" # c("ratio", "absolute")
                                ) {
    
    if (missing(target)) {
        predictors <- model$model[-1]
        target <- model$model[1]
    }
    targetName <- names(target)
    
    target <- target[[1]]
    
    numberPredictors <- ncol(predictors)     
    
    if (missing(predictorsMeans)) {
        initialPredictorsValues <- predictors %>%
            summarise_all(mean)
        
    } else {
        initialPredictorsValues <- predictorsMeans    
    }
    
    initialTargetValue <- mean(predict(model, newdata = initialPredictorsValues))
    
    sensitivityData <- sapply(initialPredictorsValues, function(x) rep(x, samplesN * numberPredictors))
    changeDF <- seq(from, to, length.out = samplesN)
    
    for (i in seq(1, numberPredictors)) {
        sensitivityData[seq((i-1)*samplesN + 1, i * samplesN), i] <- changeDF * initialPredictorsValues[[i]]    
    }
    sensitivityData <- data.frame(sensitivityData)
    
    predictedTarget <- predict(model,
                               newdata = sensitivityData)
    

    
    if (targetPrediction == "ratio")
        predictedTarget <- predictedTarget / initialTargetValue
    
    df <- data.frame(normalized.predictor.change = numeric(0),
                     predictor = character(0),
                     normalized.target.change.fit = numeric(0)) 
    
    for (i in seq(1, numberPredictors)) {
        df <- rbind(df, 
                    data.frame(
                        normalized.predictor.change = changeDF,
                        predictor = names(predictors)[i],
                        normalized.target.change.fit = predictedTarget[seq((i-1)*samplesN + 1, i * samplesN)]))
    }
    if (targetPrediction == "ratio") {
        gg <- ggplot(df, aes(x = normalized.predictor.change, 
                             group = predictor)) +
            geom_vline(xintercept = 1,  color = "grey80") + 
            geom_hline(yintercept = 1, color = "grey80") +
            geom_abline(slope = 1, linetype = "dashed", color = "grey80") + 
            geom_abline(slope = -1, intercept = 2, linetype = "dashed", color = "grey80") +
            geom_line(aes(x = normalized.predictor.change,
                          y = normalized.target.change.fit,
                          color = predictor)) +
            ylab(paste("Normalized", targetName, "change")) +
            xlab("Normalized Predictor Change") + 
            theme_few()
    } else {
        gg <- ggplot(df, aes(x = normalized.predictor.change, 
                             group = predictor)) +
            geom_vline(xintercept = 1,  color = "grey80") + 
            geom_hline(yintercept = initialTargetValue, color = "grey80") +
            geom_line(aes(x = normalized.predictor.change,
                          y = normalized.target.change.fit,
                          color = predictor)) +
            ylab(targetName) + 
            xlab("Normalized Predictor Change") + 
            theme_few()        
    }
    
    variableMeans <- cbind(data.frame(target = initialTargetValue),
                           data.frame(initialPredictorsValues))
    names(variableMeans) <- c(targetName, names(initialPredictorsValues))
    return(list(ggplot = gg,
                predictionData = sensitivityData,
                resultsData = df,
                variableMeans = variableMeans))
}
library(dplyr)
library(ggplot2)
library(caret)
library(ggthemes)

jumpData <- read.csv("jump-data.csv", header = TRUE)

# Linear model
model <- lm(jumpHeight ~ maxPower + bodyWeight + slope + pushOffDistance, jumpData) # scale & log-log
summary(model)

results <- sensitivityAnalysis(model, level = .90, predictionType = "prediction", targetPrediction = "raw")
plot(results$ggplot)
results$variableMeans

# Interactions Linear model
model <- lm(jumpHeight ~ maxPower * bodyWeight * slope * pushOffDistance, jumpData)
summary(model)
results <- sensitivityAnalysis(model, level = .9, predictionType = "prediction")
plot(results$ggplot)

# Polynomial
target <- jumpData$jumpHeight
predictors <- select(jumpData, maxPower, bodyWeight, slope)
model <- lm(jumpHeight ~ poly(maxPower, 3, raw = TRUE) * poly(slope, 3, raw = TRUE) * poly(bodyWeight, 3, raw = TRUE),
            jumpData)
summary(model)
results <- sensitivityAnalysis(model, target, predictors, level = .9, predictionType = "prediction", targetPrediction = "raw")
plot(results$ggplot)


# Caret
target <- jumpData$jumpHeight
predictors <- select(jumpData, maxPower, bodyWeight, slope, pushOffDistance)
ctrl <- trainControl(method="cv", allowParallel = TRUE, number = 5, repeats = 1, verboseIter = TRUE)
model <- train(y = target, x = predictors, 
               method = "gam",
               preProcess = c("center", "scale"),
               trControl = ctrl)

results <- sensitivityAnalysisCaret(model, target, predictors, targetPrediction = "raw")
plot(results$ggplot + ylab("Jump Height"))

}

#Distribution plots ------------
{

  library(ggplot2)
  library(zeallot)
  library("minet")
  library(h2o)
  #Load data *********************************************************************88
  url.rslts = 'E:/HDDs/Work/NYU Research/0 Plant biology/CICT paper/'
  netObjs = readRDS(file = paste0( url.rslts,  '/mESC network objects.rds')) #'/mESC_rslt_3.rds')) #
  c(rcrd,prd.varimp,n.itm.e,pred_outcome.e,a,t1) %<-% netObjs #truncated_pred_outcome
  rm(netObjs)
  
  varimp = netObjs$varimp
  
  
  d = t1 %>% select(-Weight) %>%
    inner_join(pred_outcome.e, by = c('src'='Gene1','trgt'='Gene2')) %>% 
    arrange(desc(Weight)) %>% dplyr::filter(Weight > .5 | Weight <.1)
  

  # Bin data ********************************************************************888
  nbins= 30 # as.integer(sqrt(nrow(t1)/2)/20)
  breaks =1:nbins
  

  d.binned =data.table(grpvar=c(rep(T,nbins-1),rep(F,nbins-1)))
  d.distDstnc = c()
  # = d %>% select(src,trgt,outcomes, class2, class3,vaimp[1:10]) %>% 
  #   mutate(grpvar =class2) # as.factor(class3)) 
  d.var = 'scfKurt.x' # 'scftau4.x'
  
  #Find variables with best visual distinct distributions
  for(d.var in  prd.varimp$variable[1:80]){
    print(d.var)
    edges.causal = myhist(infotheo::discretize(d[ class2 == T, d.var,with=F] %>% unlist(), "equalwidth", nbins)$X,nbins,breaks,prob=F)
    edges.random =  myhist(infotheo::discretize(d[ class2 == F, d.var,with=F] %>% unlist(), "equalwidth", nbins)$X,nbins,breaks,prob=F)
    
    d.binned[,(d.var):= c(edges.causal,edges.random)]
    
    tmp = data.frame(edges.causal,edges.random)
    mim <- build.mim(tmp, disc = 'equalwidth', estimator = 'kendall')
    dist = mim[1,2]; names(dist)<-d.var
    d.distDstnc = c(d.distDstnc,dist)
  }
  d.distDstnc=sort(desc(d.distDstnc))
  d.topDstnc = names(d.distDstnc)[c(1,3,6,7,8)]# names(d.distDstnc)[c(1,3,4)]
  var =names(d.distDstnc)[3]
  
  # Density plot *****************************************************************888
  library(easyGgplot2)
  library(cowplot)
  {
    #palette https://personal.sron.nl/~pault/
    
    desiredVars = 'scbSD.x'  #'ocfNMAD.y' #c('scftau4.x') #,'ocfNMAD.y') # 'scbKurt.y' # 
    tmp = setDT(d.binned)[,c(desiredVars,'grpvar'),with=F] %>% setDF() %>%
     mutate(!!desiredVars:= log2(abs(.data[[desiredVars]])) )  #%>% pivot_longer(cols = desiredVars, names_to = 'var')
    rm(plt.d);plt.d = ggplot2.density(data=tmp , xName=desiredVars  , groupName = 'grpvar' ,
                    alpha=0.7, fillGroupDensity=T,
                    #colorGroupDensityLine=T,densityFill  =c('#66CCEE','#EE6677'),
                    scale = 'density',
                    addMeanLine=F, meanLineColor=NULL, meanLineSize=1,xlim=c(0,25)
                    )#+facet_wrap(desiredVars);
    
    plt.d= plt.d +theme_half_open(font_size = 9,
                                  font_family = "Arial") +
      background_grid() + 
      theme(axis.text.x = element_text( size = 9,angle=0,hjust=.2),
            legend.position = c(.1,.98),legend.direction="horizontal",
            legend.text = element_text( size = 9,angle=0),
            legend.title = element_text( size = 9,angle=0),
            
            )+ #,    legend.key.size = unit(1, 'cm'))+
      #labs(colour = 'Causal edges',fill ='Causal edges',  )+
      scale_fill_discrete(name = '', labels =c( "Random Edges", "Causal edges"),
                          type  =c('#66CCEE','#EE6677'))+
      guides(fill = guide_legend(reverse=T))+ 
      xlab('Log2 standard deviation of contributions of source')+ #"Median axis deviation of confidences of target")+
      #xlab("L-Kurtosis of source confidence")+
      ylab("Density");plt.d #, 
    
    
    thepath = "E:/HDDs/Work/NYU Research/0 Plant biology/R code/results/scbSD.x.tif"
    ggsave(filename = thepath, plot = plt.d, width = 3.5, height = 2.5, device='tiff', dpi=300)
    
    if(F) ggplot2.histogram(data=tmp , xName='var'  , groupName = 'grpvar' ,
                            legendPosition="top",  #groupName='',
                            alpha=0.8, densityFill ="cyan",addDensityCurve =TRUE,
                            scale = 'density',
                            addMeanLine=TRUE, meanLineColor=NULL, meanLineSize=1,xlim=c(-1,10)
                            
    )
    

    
    #plt.d + scale_fill_brewer(palette="Moonrise3", direction=-1) # scale_colour_manual(values = rev(brewer.pal(3, "BuPu"))) #
  }
  #ggplot(n.itm.e, aes(EE))+ggtitle("EE") + geom_histogram()  + scale_x_log10()
  

  #Violin plot ***********************************************************************
  {
    studyDataset = "HSEC timeset edges"
    
    theTitle = '' #paste0(studyDataset, " - Dist. top predictors of causal and random edges")
    
    plt.violin.cols.1 =c(prd.varimp$variable[c(1:3,5,7)]) # c(names(d.distDstnc)[1:5]) #
    df=setDF(d.binned) %>%  #setDF(d) %>% mutate(grpvar = class2) %>%
      select(plt.violin.cols.1,'grpvar') %>%
      mutate(eID = row_number()) %>%
      mutate(class10 = grpvar , class11= !grpvar)  %>%
      select(-grpvar)
    #select(c(grpvar,trainingTarget, ggp.vars))  #tst1.totalset #newDataEnv$tst1.totalset # 
    
    names(df)[1:5]<-c('Standard dev. of source contrib.',
                      'L-Kurtosis of source contrib.',
                      'L-Skewness of source conf.',
                      'L-Skewness of source contrib.',
                      'Median of source contrib.')

    # names(df)[1:5]<-c('Kurtosis of source conf.',
    #                   'L-Skewness of source conf.',
    #                   'L-Kurtosis of source conf.',
    #                   'L-Skewness of source contrib.',
    #                   'Median of source contrib.')
    
    plt.violin.cols.2 = names(df)[1:length(plt.violin.cols.1)]
    #The most useful of these are tau _{3}called the L-skewness, and tau _{4} the L-kurtosis. 
    
    dfl = df %>%  pivot_longer(plt.violin.cols.2) %>% # setdiff(colnames(df),c('grpvar','eID')))  
      dplyr::mutate(logval = sign(value) * log2(abs(value))) 
    
    library(plotly)
    plt.v=myplotViolin(theTitle =theTitle,
                 dfl,plt.violin.cols.2, 
                 group1 = 'class10', group2 = 'class11',
                 grp1title= 'Causal edges', grp2title = 'Random edges',
                 grpcolor1 = '#EE6677',grpcolor2 = '#66CCEE',
                 yval = 'logval',thebox=T
                 )
    
    t1 <-list(family = "Arial",size = 22,color = "black",bold=T)
    t2 <- list(family = "Arial",size = 22,color = "black",bold=T)
    t3 <- list(family = 'Arial')
    t4= list(family='Arial', color='black', size=22)
    
    plt.v=layout(plt.v
           ,font=t4 
           #,title= list(text = "",font = t1)  
           ,legend=list(font = t1,title=list(text='',font = t2),orientation = "h",   # show entries horizontally
                       # xanchor = "center",  # use center of legend as anchor
                       y=1, 
                       x = 0.2)
           #,annotations =list( font = list(size = 10)),
           ,yaxis = list(title =list(text="Log2(frequency of binned values)", font = t1), 
                         range = c(-3,25))
           ,xaxis = list(title = "", range = c(-1,5),tickangle = 315,
                         tickfont =t2, font =t2), #crimson
           #autosize=T,
           width=800,
           height=600
           #scale = (5*2.54*10) * 17780.0 
           );plt.v
    
    #scale = width_in_mm * 17780.0 # dpi=300
    # Jpeg <- plotly::plotly_IMAGE(plt.v, format = "jpeg", 
    #                      width = 1800, height=1800,
    #                      out_file = thepath)
    
    thepath = "E:/HDDs/Work/NYU Research/0 Plant biology/R code/results/violinTop5.png"
    #ggsave(filename = thepath, plot = plt.d, width = 3, height = 2, format='tiff', dpi=300)
    
    plotly::export(plt.v, file = thepath) #, width = 3, height = 2, format='png', dpi=300)
  }

}

#Sunkey and parallel graph
{
  #Parallelogram  Decision boundaries
  {
    trainingTarget='class2'
    theTitle =paste0(studyDataset, " - Dist. top predictors of causal and random edges")
    
    ggp.vars =c(prd.varimp$variable[c(1:3,5,7)]) # c(names(d.distDstnc)[1:5]) #

    nbins=3
    #ggp.vars = setdiff(ggp.vars, str_subset(ggp.vars,"z"))[1:15]
    dfs = setDF(d) %>% select(unique(c('class1', ggp.vars))) %>%
       mutate(class1=plyr::mapvalues(class1,c('c','rc','u'),c('Causal','Reverse_causal','Random'))) %>%
      # mutate(class1=as.factor(class1) ) %>%
      # mutate_at(.vars = ggp.vars,
      #                ~ log2(infotheo::discretize(., "equalwidth", nbins)$X)) %>%
      # mutate_at(.vars = ggp.vars, ~ paste0(class1,'_',.))
      # #group_by(class1) %>% summarise_if(is.numeric, list(mySD) , na.rm = TRUE)
      # 
      # group_by(class1) %>% summarise_at(ggp.vars, 
      #                                   ~(infotheo::discretize(., "equalwidth", nbins)$X))
    
      group_by(class1) %>% summarise_if(is.numeric, list(myMedian) , na.rm = TRUE)
    
    names(dfs)<- plyr::mapvalues(colnames(dfs),
                      c("scbSD.x" ,"scbL4.x" ,"scftau3.x", "scbtau4.x" , "scbNMADconst.x"),
                      c('Standard dev. of source contrib.',
                        'L-Kurtosis of source contrib.',
                        'L-Skewness of source conf.',
                        'L-Skewness of source contrib.',
                        'Median of source contrib.'))
    ggp.vars.new = c('Standard dev. of source contrib.',
                     'L-Kurtosis of source contrib.',
                     'L-Skewness of source conf.',
                     'L-Skewness of source contrib.',
                     'Median of source contrib.')
    #plt.violin.cols.2 = names(df)[1:length(plt.violin.cols.1)]
    
    ggparallel::ggparallel(vars = ggp.vars.new, data = dfs,width = 0.5,
                           alpha = .5,ratio = .2,
                           method = "adj.angle") # 'hammock') #hammock
    
    ggp =GGally::ggparcoord(dfs, columns = 2:6, groupColumn = 'class1', scale = "std",
                             splineFactor = 0,  boxplot = F, alphaLines=0.6,#shadeBox ='lightgray',
                    title = "");ggp #Top 5 CICT predictors of connection class
    
    
    
    ggp = ggp + geom_line(size = .8)  + geom_point(size =1.5) + #theme_cowplot()+
      theme_half_open() +
      background_grid() + 
       theme(
         text = element_text(size=12,family='Arial'),
         axis.text.x = element_text(size = 12,angle=-315,hjust=1),
         plot.margin = margin(t = 0,  # Top margin
                              #r = 1,  # Right margin
                              b = 0,  # Bottom margin
                              l = 50,  # Left margin
                              unit = "pt"))+
      theme(legend.position="top")  +
      ylab("Standardized median")  + xlab('')+
      guides(colour=guide_legend(title=""))+
      scale_color_manual(values=c('Reverse causal'='lightgray',
                                  'Causal' =  '#EE6677',
                                  'Random'='#66CCEE'
                                  
                                )) ;ggp
   
    thepath = "E:/HDDs/Work/NYU Research/0 Plant biology/R code/results/ParallelGraph.tif"
    ggsave(filename = thepath, plot = ggp, width = 5, height = 4, device='tiff', dpi=300)
    
    
    
    
    snkn.d = data.table()
    d.binned =data.table(grpvar=c(rep(T,nbins-1),rep(F,nbins-1)))
    d.distDstnc = c()
    # = d %>% select(src,trgt,outcomes, class2, class3,vaimp[1:10]) %>% 
    #   mutate(grpvar =class2) # as.factor(class3)) 
    d.var = 'scfKurt.x' # 'scftau4.x'
    
    for(d.var in  plt.violin.cols.1){
      snkn.d = data.table(
        infotheo::discretize(d[ class2 == T, d.var,with=F] %>% unlist(), "equalwidth", nbins)$X
      )
      edges.random =  myhist(infotheo::discretize(d[ class2 == F, d.var,with=F] %>% unlist(), "equalwidth", nbins)$X,nbins,breaks,prob=F)
      
      d.binned[,(d.var):= c(edges.causal,edges.random)]
      
      tmp = data.frame(edges.causal,edges.random)
      mim <- build.mim(tmp, disc = 'equalwidth', estimator = 'kendall')
      dist = mim[1,2]; names(dist)<-d.var
      d.distDstnc = c(d.distDstnc,dist)
    }
    d.distDstnc=sort(desc(d.distDstnc))
    d.topDstnc = names(d.distDstnc)[c(1,3,4)]
    
    
    df=setDF(d.binned) %>%  #setDF(d) %>% mutate(grpvar = class2) %>%
      select(plt.violin.cols.1,'grpvar') %>%
      mutate(eID = row_number()) %>%
      mutate(class10 = grpvar , class11= !grpvar)  %>%
      select(-grpvar)
    
    
  }
  
  {
    
    df=setDF(d) %>% mutate(grpvar = class2) %>%
      select(plt.violin.cols.1,'grpvar') %>%
      mutate(eID = row_number()) %>%
      mutate(class10 = grpvar , class11= !grpvar)  %>%
      select(-grpvar)
    #select(c(grpvar,trainingTarget, ggp.vars))  #tst1.totalset #newDataEnv$tst1.totalset # 
    
    names(df)[1:5]<-c('Kurtosis of source conf.',
                      'L-Skewness of source conf.',
                      'L-Kurtosis of source conf.',
                      'L-Skewness of source contrib.',
                      'Median of source contrib.')
    
    plt.violin.cols.2 = names(df)[1:length(plt.violin.cols1)]
    #The most useful of these are tau _{3}called the L-skewness, and tau _{4} the L-kurtosis. 
    
    dfl = df %>%  pivot_longer(plt.violin.cols.2) %>% # setdiff(colnames(df),c('grpvar','eID')))  
      dplyr::mutate(logval = sign(value) * log2(abs(value))) 
    
    sankeyNetwork(Links = Energy$links, Nodes = Energy$nodes, Source = "source",
                  Target = "target", Value = "value", NodeID = "name",
                  units = "TWh", fontSize = 12, nodeWidth = 30)
  }
}
#########################***
# 3d graph sensitivity analysis
{
  library(plotly)
  plot_ly(data=sns.map,  x= ~theoffset,y= ~strategyLag, z= ~logReturn) %>% 
    group_by(strategyLag) %>% 
    add_lines(color = ~score,width=3) %>% 
    layout( title = paste0(hint,' ',analysisTyp)) #add_surface()
  
}

#Random forest learning
{
  library(ggplot2)
  library(cowplot)
  s.training = data.frame(Number_of_trees=c(5,10,20,30,40,50),
                          Logloss=c(.6,.25,.17,.13,.125,.123))
  s.training$series = 'Training'
  s.validation = data.frame(Number_of_trees=c(5,10,20,30,40,50),
                            Logloss=c(.3,.18,.14,.125,.121,.120))
  s.validation$series = 'Validation'
  
  s.all = rbind(s.training,s.validation)
  
  p <- ggplot(s.all, aes(Number_of_trees, Logloss,colour = series)) +
    geom_point(size = 4) 
  # Add regression line
  p + geom_line(method = lm,size=2)+
    theme_half_open(font_size = 35,font_family = "Arial") +
    background_grid() +
    theme(legend.position="top")+
    guides(colour=guide_legend(title=""))+
    xlab('Number of trees in CICT model')
  
  
  # ggplot(s.all)+
  # geom_point(aes(x = Number_of_trees, y = Logloss, colour = series), size = 3) +
  #   stat_smooth(aes(x = Number_of_trees, y = Logloss, colour = series), method = "lm",
  #               formula = y ~ poly(x, 2), se = FALSE) #+     coord_cartesian(ylim = c(0, 1.5e7))
  
}


#Random forest AUC
{
  url.rslts = 'E:/HDDs/Work/NYU Research/0 Plant biology/CICT paper/mESC_rslt_.rds'
  rcrd = readRDS(url.rslts)
  p = rcrd$learning_curve_plot_aucpr
  p=rcrd$learning_curve_plot_logloss
  p$layers[c(1,2,5)]<-NULL
  
  p = p+  theme_half_open(font_size = 18,font_family = "Arial") +
    background_grid() +
    theme( legend.position="none")+
    guides(colour=guide_legend(title=""))+
    xlab('Number of trees in CICT model')+
    labs(title = '', subtitle = '')+
    xlim(0,25);p
   
  thepath = "E:/HDDs/Work/NYU Research/0 Plant biology/R code/results/CICT RF training.tif"
  ggsave(filename = thepath, plot = p, width = 4.5, height = 4, device='tiff', dpi=300)
  
  
  p= rcrd$sscurves.prc + theme_half_open(font_size = 18,font_family = "Arial") +
    background_grid() +
    theme( legend.position="none")+
    guides(colour=guide_legend(title=""))+
    xlab('Sensitivity')+ylab('Precision')+
    labs(title = '', subtitle = '', );p
  
  thepath = "E:/HDDs/Work/NYU Research/0 Plant biology/R code/results/CICT RF PRC AUC.tif"
  ggsave(filename = thepath, plot = p, width = 4.5, height = 4, device='tiff', dpi=300)

  p= rcrd$sscurves.roc + theme_half_open(font_size = 18,font_family = "Arial") +
    background_grid() +
    theme( legend.position="none")+
    guides(colour=guide_legend(title=""))+
    xlab('1 - Specificity')+ylab('Sensitivity')+
    labs(title = '', subtitle = '', );p
  
  thepath = "E:/HDDs/Work/NYU Research/0 Plant biology/R code/results/CICT RF ROC AUC.tif"
  ggsave(filename = thepath, plot = p, width = 4.5, height = 4, device='tiff', dpi=300)
  
  
  
}


#Network rendering ----
{
  library(zeallot)
  library(SummarizedExperiment)
  require(verification)
  require(pROC)
  require(ROCR)
  require( OptimalCutpoints)
  require(precrec )
  library(dplyr)
  
  library(RJSONIO)
  library(igraph)
  library(httr)  
  library(networkD3)
  library(RCy3)
  library(htmlwidgets)
  library(htmltools)
  if(!"RCy3" %in% installed.packages()){
    install.packages("BiocManager")
    BiocManager::install("RCy3")
    # Basic settings
    # port.number = 1234
    # base.url = paste("http://localhost:", toString(port.number), "/v1", sep="")
    # version.url = paste(base.url, "version", sep="/")
    # cytoscape.version = GET(version.url)
    # cy.version = fromJSON(rawToChar(cytoscape.version$content))
    # print(cy.version)
  }
  
  # url.rslts = 'E:/HDDs/Work/NYU Research/0 Plant biology/CICT paper/SensitivityRound2/'
  # netObjs = readRDS(file = paste0( url.rslts, 'mESC network objects.rds'))
  # c(rcrd,prd.varimp,n.itm.e,pred_outcome.e,truncated_pred_outcome,t1) %<-% netObjs
  # 
  url.rslts = 'E:/HDDs/Work/NYU Research/0 Plant biology/CICT paper/SensitivityRound2/mESC_rslt_.rds'
  rcrd = readRDS(url.rslts)
  # 
  # rcrd = netObjs[[1]]
  # n.itm.e=netObjs[[3]]
  # pred_outcome.e=netObjs[[4]]
  

  
  
  thepath = "E:/HDDs/Work/NYU Research/0 Plant biology/CICT paper/SensitivityRound2/"
  tbl.goldStandard = read.csv(paste0(thepath,'mESC L2_lofgof-refnetwork.csv'))
  tbl.goldStandard = tbl.goldStandard %>%
    mutate(Gene1 = tolower(Gene1), Gene2 = tolower(Gene2)) %>%
    mutate(
      outcomes =1,
      Gene1 = str_replace_all(Gene1, "^\\w{1}", toupper) %>% unlist(),
      Gene2 = str_replace_all(Gene2, "^\\w{1}", toupper) %>% unlist()
    )     
  
  gnames = c(tbl.goldStandard$Gene1,tbl.goldStandard$Gene2) %>% unique() %>% sort() %>% tolower()
  gnames = gnames %>% str_replace_all( "^\\w{1}", toupper) %>% unlist()
  gnames %>% paste0(collapse=',')
  
  #Threshold gated network
  {
    n.itm.e = n.itm.e %>% 
      mutate(src = tolower(src), trgt = tolower(trgt)) %>%
      mutate(
        src = str_replace_all(src, "^\\w{1}", toupper) %>% unlist(),
        trgt = str_replace_all(trgt, "^\\w{1}", toupper) %>% unlist()
      )     
    
    n.itm.e.rawedge = n.itm.e %>% 
      dplyr::select(src, trgt,Euclidean,Spearman,Kendall,Pearson,efMIempirical,ewMIshrink,ewMIempirical,ewMImm)
  
  

    
    
    pred_outcome = pred_outcome.e %>% 
      dplyr::rename(predictions = EdgeWeight) %>% 
      mutate(Gene1 = tolower(Gene1), Gene2 = tolower(Gene2)) %>%
      mutate(Gene1 = str_replace_all(Gene1, "^\\w{1}", toupper) %>% unlist(),
             Gene2 = str_replace_all(Gene2, "^\\w{1}", toupper) %>% unlist()
              ) %>%
      dplyr::select(!any_of('outcomes'))
    
    pred_outcome= pred_outcome %>%
      dplyr::left_join(n.itm.e.rawedge,by=c(Gene1 = 'src' , Gene2 = 'trgt')) %>%
      dplyr::left_join(tbl.goldStandard,by=c('Gene1','Gene2')) 
    pred_outcome = pred_outcome %>%
      mutate(outcomes = replace_na(outcomes,0))
    #   rename(src='Gene1',trgt='Gene2')
    #table(pred_outcome$outcomes)  
    
    relativeCostfn_fp = 1/5
    
    #find best cutoff
    prv = table(pred_outcome$outcomes)
    best.weights=c(relativeCostfn_fp, prv[2]/prv[1]) 
    theROC <- roc(pred_outcome$outcomes, pred_outcome$predictions, percent = TRUE);theROC
    bestcutoff =as.double(coords(theROC, "best", best.method =  "closest.topleft",
                                 best.weights=best.weights,
                                 ret="threshold", transpose = FALSE));bestcutoff
    
    n3d.e = pred_outcome %>%
      dplyr::filter(
        #predictions > rvpred & 
        predictions> bestcutoff)
  
    n3d.e = n3d.e %>% #inner_join(n.itm.e[,.(src,trgt,Weight)])%>% 
      #dplyr::filter( predictions > bestcutoff | seenInTraining ==1.0 | seenInGS==1.0 ) %>%
      mutate(interactions =  case_when(
        outcomes==T & predictions > bestcutoff ~ 'TP'
        ,outcomes==T & predictions < bestcutoff ~ 'FN'
        ,outcomes==F & predictions > bestcutoff ~ 'FP'
        ,outcomes==F & predictions < bestcutoff ~ 'TN'
      )) 
    table(n3d.e$interactions)  
    
    n3d.e = n3d.e %>% 
      mutate(Gene1 = tolower(Gene1), Gene2 = tolower(Gene2)) %>%
      mutate(
        Gene1 = str_replace_all(Gene1, "^\\w{1}", toupper) %>% unlist(),
        Gene2 = str_replace_all(Gene2, "^\\w{1}", toupper) %>% unlist()
      ) 
    
    #n3d.e %>% filter(Gene2 == 'Pou5f1') %>% pull(Gene1) %>% paste0(collapse=',')
    }
  
  
  #n3d.e1 = n3d.e %>% dplyr::filter(predictions > rvpred)

  n3d.ef=n3d.e%>% 
    dplyr::filter(interactions != 'TN' )
  n.itm.v =data.frame(vID = c(n3d.ef$Gene1,n3d.ef$Gene2) %>% unique())
  v.out = n3d.ef %>% group_by(Gene1) %>% 
    summarise(OcrOut=sum(predictions,na.rm=TRUE)) %>% 
    dplyr::rename(vID = Gene1)
  
  v.in = n3d.ef %>% group_by(Gene2) %>% 
    summarise(OcrIn=sum(predictions,na.rm=TRUE)) %>% 
    ungroup() %>% dplyr::rename(vID = Gene2)
  
  n.itm.v = n.itm.v %>% left_join(v.out) %>% left_join(v.in) 
  n3d.v = n.itm.v
  

  #str(n3d.e);str(n3d.v)
  
  

  library(RCy3)
  cytoscapePing ()
  cytoscapeVersionInfo ()
  

  
  
  #Predicted network
  {
    cynodes <- n3d.v  %>% dplyr::rename(id = vID, score  = OcrOut)  #%>% dplyr::select(vID,OcrOut,group)
    cyedges <- n3d.ef %>% #dplyr::select(src,trgt,Weight,predictions ) %>% 
      dplyr::rename(source = Gene1, target  = Gene2) 
    
    
    createNetworkFromDataFrames(cynodes,cyedges, title="mESC_CICT - cost 1/5", collection="mESC-CICT")
  }
  
  
  #Best network suggested by a relevance edge measure: ewMIshrink
  {
    #   MIX the Data + Ground truth  ----
    {
      #NrandomEdges =  min(.3 * nrow(n.itm.e),200000) ; #NrandomEdges=2e6
      
      
      #or 0 means all
      #n.itm.e = n.itm.e %>% filter(src %in% geneorder.2e3$X | trgt %in% geneorder.2e3$X)
      
      
      t1.c = n.itm.e %>% setDF() %>% 
        inner_join(tbl.goldStandard %>% dplyr::select(-outcomes),
                   by=c("src"="Gene1","trgt"="Gene2"))
      if(nrow(t1.c) < nrow(tbl.goldStandard)/2) warning("Problem in ground truth. More than half of ground truth was not found in edges")
      if(nrow(t1.c) <=0) print("!!! No causal edge in the groundturht? check gene names")

      t1.rc = n.itm.e  %>% dplyr::inner_join(t1.c %>% dplyr::select(src,trgt),by=c("src"="trgt","trgt"="src"))
      
      t1.c$predicate = "CAUSES"
      t1.rc$predicate = "REV_CAUSES"
      
      t1.causal_reversecausal=rbind(t1.c,t1.rc)
      
      #Adding 2000 random edges
      t1.rnd = n.itm.e %>% #dplyr::mutate(edgetyptruth = NA) %>%
        anti_join(t1.c,by=c("src"="src","trgt"="trgt")) %>%
        anti_join(t1.rc,by=c("src"="trgt","trgt"="src"))
      t1.rnd$predicate = "IRRELEVANTM"
      
      t1 = rbind(t1.c,t1.rc,t1.rnd) #%>% unique()
      #plyr::arrange(Weight)%>% dplyr::slice_head(n=1100) #s0c2 L- 
      t1$class1=hutils::Switch(t1$predicate, 
                               CAUSES = 'c', PRECEDES= 'el' ,ASSOCIATED_WITH= 'a',
                               REV_CAUSES= 'rc',IRRELEVANTA= 'ir', IF_NA= 'u', DEFAULT = 'u')
      
      t1 = t1 %>% 
        dplyr::mutate(SUID=row_number(),
                      class2=ifelse(class1 %in% c('c'),TRUE,FALSE),
                      class3=ifelse(class1 %in% c('c','rc'),TRUE,FALSE),
                      #class5=ifelse(edgetyptruth %in% c('+'),TRUE,FALSE),
                      shared_name=paste0(src,"-",trgt))   #class1=""
      
      setDT(t1)[,`:=`(shared_name=NULL,DiagnosesAb.x=NULL,DiagnosesAb.y=NULL)]
      #t2=unique(t2)
      table(t1$predicate,t1$class2, useNA="ifany")
      table(t1$edgetyptruth,t1$class5, useNA="ifany")
      
      
      
      maxGroundTruth=500
      minGroundTruth.ratio.learning=.2
      randomEdgesFoldCausal = 5
      sampling.c_rc.ratio = .3 
      
      NcausalEdges = min(maxGroundTruth, nrow(t1.c)*minGroundTruth.ratio.learning) %>% as.integer()
      if(NcausalEdges<300)NcausalEdges=floor(2*nrow(t1.c)*minGroundTruth.ratio.learning)
      
      NrandomEdges = NcausalEdges *  randomEdgesFoldCausal #min(sampling.rnd.ratio * nrow(n.itm.e),sampling.u.max) %>% as.integer() ; #NrandomEdges=2e6
      
      #t2 is the learning set
      t2=rbind(t1 %>% filter(class1=='c') %>% sample_n(size = NcausalEdges),
               t1 %>% filter(class1=='rc') %>% sample_n(size = NcausalEdges),
               t1 %>% filter(class1=='u')%>% sample_n(size=NrandomEdges-2*NcausalEdges)) #to preserve the minGroundTruth.ratio
      
      t2.complement = t1  %>%  anti_join(t2,by=c("src"="src","trgt"="trgt"))
      
      t2= as.data.frame(t2) %>% 
        mutate(predicate = class1) %>%
        dplyr::select(any_of(c('predicate', 'class2','src', 'trgt',
           'Euclidean','Spearman','Kendall','Pearson','efMIempirical','ewMIshrink','ewMIempirical','ewMImm')))
      table(t2$predicate)#,t2$edgeTyp, useNA="ifany")
      a=sapply(t2, function(x) sum(is.na(x)));a[a>0]
      

      
    }
    
    #Train model
    {
    tst1.totalset = t2
    tstPrecent = .3
    
    set.seed(as.integer(runif(1,1,10000)))
    spltIdx =as.vector( caret::createDataPartition(1:nrow(tst1.totalset),p=(1-tstPrecent),list=FALSE,times=1))
    tst1.train = tst1.totalset[spltIdx,] #tst1.dmy[spltIdx,]
    tst1.tst = tst1.totalset[-spltIdx,]
    
    paste0(colnames(tst1.tst),collapse="','")
    
    table(tst1.train$class2);table(tst1.tst$class2)

    require(h2o)
    H2OCnn = h2o.init(nthreads = 8, enable_assertions = TRUE,
                      strict_version_check=FALSE) #max_mem_size = "5gb",

    NFolds = 5
    mdlColNames=c('class2','Euclidean','Spearman','Kendall','Pearson','efMIempirical','ewMIshrink','ewMIempirical','ewMImm')
    
    tst1.tst.h2o=as.h2o(tst1.tst) 
    tst1.train.h2o = as.h2o(tst1.train)
     tst1.rf = h2o.randomForest(mdlColNames,trainingTarget,
                               tst1.train.h2o,nfolds = NFolds,#checkpoint=tst2.mdl@model_id,
                               max_depth =20 , ntrees=50, 
                               keep_cross_validation_predictions= FALSE,  validation_frame = tst1.tst.h2o)
    #,ntrees=30,max_depth=5#selectedFeatures #mdlColNames
    
    tst1.mdl=tst1.rf
    
    tst1.rf@model$validation_metrics@metrics
    
    d.new=t1  %>%  anti_join(t2,by=c("src"="src","trgt"="trgt"))
    maxunseenTest.ratio = .1
     #%>% sample_n(size = min(50000,nrow(t2.complement) * maxunseenTest.ratio))   #tst1.totalset #newDataEnv$tst1.totalset # 
    
    #d.new = n.itm.e
     splitcount = max(floor(nrow(d.new) / 30000),2)
    splts = split(1:nrow(d.new),          
                  cut(seq_along(1:nrow(d.new)),
                      splitcount,labels = FALSE))
    select<-dplyr::select
    d.new.tmp = lapply(splts,
                       function(thesplit){
                         print(last(thesplit))
                         d.new.slice = d.new[thesplit,]
                         predTest.d.h2o = as.h2o(d.new.slice) # 
                         h2o.pred = as.data.frame(h2o.predict(tst1.mdl,predTest.d.h2o,keep_cross_validation_predictions=F))
                         h20.prediction=as.numeric(as.character(h2o.pred[,3]))
                         predictions =h20.prediction #pred #ens.predictions #
                         outcomes =unlist(setDT(d.new.slice)[,trainingTarget,with=FALSE]) #outcome # ens.outcome# 
                         set.seed(runif(1,1,1000))
                         # = ifelse(runif(nrow(d.new.slice)) >=.5,1,0) #Assign a random classifier results
                         prd_outcomes = d.new.slice %>% select(src,trgt,any_of('Weight')) %>%
                           cbind(predictions,outcomes) %>% as.data.frame()
                         prd_outcomes
                       })
    
    pred_outcome.rlv= rbindlist(d.new.tmp)
    }
    
    relativeCostfn_fp = 1/5
    
    #find best cutoff
    prv = table(pred_outcome.rlv$outcomes)
    best.weights=c(relativeCostfn_fp, prv[2]/prv[1])
    
    pred_outcome.rlv = pred_outcome.rlv %>% mutate(predictions.rlv = predictions )#Pearson
    theROC.rlv <- roc(pred_outcome.rlv$outcomes, pred_outcome.rlv$predictions.rlv, percent = TRUE);theROC
    bestcutoff.rlv =as.double(coords(theROC.rlv, "best", best.method =  "closest.topleft",
                                 best.weights=best.weights,
                                 ret="threshold", transpose = FALSE));bestcutoff.rlv
    
    n3d.e.rlv = pred_outcome.rlv %>% 
      arrange(desc(predictions.rlv)) %>%
      dplyr::slice_head(n= nrow(n3d.e))
      #dplyr::filter(predictions.rlv> bestcutoff.rlv)
        #predictions > rvpred & 
        
    
    n3d.e.rlv = n3d.e.rlv %>% #inner_join(n.itm.e[,.(src,trgt,Weight)])%>% 
      #dplyr::filter( predictions > bestcutoff | seenInTraining ==1.0 | seenInGS==1.0 ) %>%
      mutate(interactions =  case_when(
        outcomes==T & predictions.rlv > bestcutoff ~ 'TP'
        ,outcomes==T & predictions.rlv < bestcutoff ~ 'FN'
        ,outcomes==F & predictions.rlv > bestcutoff ~ 'FP'
        ,outcomes==F & predictions.rlv < bestcutoff ~ 'TN'
      )) #%>%       filter(interactions != 'TN' )
    table(n3d.e.rlv$interactions)  
    
  
    n3d.e.rlv =n3d.e.rlv %>% setDF() %>% 
      mutate(Gene1=src,Gene2= trgt) %>% select(-src,-trgt) %>%
      mutate(Gene1 = tolower(Gene1), Gene2 = tolower(Gene2)) %>%
      mutate(
        Gene1 = str_replace_all(Gene1, "^\\w{1}", toupper) %>% unlist(),
        Gene2 = str_replace_all(Gene2, "^\\w{1}", toupper) %>% unlist()
      ) 
    
  }
  
}

#Remove large objects
sort( sapply(ls(),function(x){object.size(get(x))}),decreasing = F)
#https://advaitabio.com/faq-items/understanding-gene-ontology/
#https://girke.bioinformatics.ucr.edu/GEN242/tutorials/rfea/rfea/



#Enrichment analysis on Gene clusters
{
  
  library(clusterProfiler)
  library(org.Mm.eg.db);
  library(org.Hs.eg.db)
  library(GO.db)
  
  
  columns(org.Mm.eg.db)
  keytypes(org.Mm.eg.db)
  org.Mm.eg.db$a
  
  xx <- as.list(org.Mm.eg.db) #org.Hs.eg.db #"org.Rn.eg.db"
  # gc1 = cl.nodes %>% group_by(type) %>% 
  #   summarise( geneset = c(key) )
  # summarise(vector=paste(A, collapse=" "))
  # 
  
  #select(org.Mm.eg.db, keys=genes, columns=c("SYMBOL","ENTREZID","ENSEMBL"), keytype="GENENAME")
  sample_gene <- sample(keys(org.Mm.eg.db), 100)
  
  
  gene.symbols <- n.itm.v$vID # cl.nodes %>% dplyr::filter(keys %in% names(unlist(gc2)))
  
  onto = Ontology(GOTERM)
  
  goterms =data.frame(term=unlist(Term(GOTERM)))
  
  #"pluripotency"  stem cell population maintenance(GO:0019827)
  goterms.pp = str_subset(goterms$term,"embryo") %>% unlist()
  View(goterms.pp)
  
  
  pp.rslr <- AnnotationDbi::select(org.Mm.eg.db, keys=c("GO:0019827"), 
                                   columns = c('SYMBOL'), keytype = "GOALL")
  pp.symbs <- unique(pp.rslr$SYMBOL) #Condition specific genes
  
  #1 Top 20 genes with most targets in CICT Network
  {
    pp.genesTargetingPP =n3d.e %>% group_by(Gene1) %>% 
      summarize(trgtCount = sum(Gene2 %in% pp.symbs)) %>%
      arrange(desc(trgtCount)) %>% 
      mutate(in_grp= Gene1 %in% pp.symbs) %>%
      slice_head(n=100)
    
    
    
    p<-ggplot(pp.genesTargetingPP, aes(x=reorder(Gene1, trgtCount), y=trgtCount, fill=in_grp)) +#
      geom_bar(stat="identity")+theme_minimal() + coord_flip() +
      xlab('')+ylab('Count of targets identified')+
      background_grid() +
      theme( legend.position="none")+
      guides(colour=guide_legend(title=""));  p

    thepath = "E:/HDDs/Work/NYU Research/0 Plant biology/R code/results/GO enriched perf CICT.tif"
    ggsave(filename = thepath, plot = p, width = 3, height = 4, device='tiff', dpi=300)
  }
  
  
  #2 Enrichment for go term "pluripotency"  stem cell population maintenance(GO:0019827)
  {
    enr_go = enrichGO(n.itm.v$vID ,OrgDb=org.Mm.eg.db, keyType = "SYMBOL", ont = "BP",
                      pvalueCutoff = 1, pAdjustMethod = "BH", qvalueCutoff = 1, 
                      minGSSize = 10, maxGSSize = 500, readable = FALSE)
    enr_go_select = enr_go  %>% 
      as.data.frame() %>%
      dplyr::filter(ID=='GO:0019827')
  }
  
  #3 Enrichment for GO terms in CICT network
  {
    rm(enr_go,enr_go_select,p)
    
    #n.itm.v$vID  #pp.genesTargetingPP$Gene1
    enr_go = enrichGO(pp.genesTargetingPP$Gene1,OrgDb=org.Mm.eg.db, keyType = "SYMBOL", ont = "BP",
                       pvalueCutoff = 1, pAdjustMethod = "BH", qvalueCutoff = 1, 
                       minGSSize = 10, maxGSSize = 500, readable = FALSE)
    
    enr_go_select = enr_go  %>% as.data.frame() %>% dplyr::arrange(desc(Count)) %>%
      slice_head(n=20) %>% mutate(pvalNegLog10 = -log10(pvalue)) %>% dplyr::arrange(pvalNegLog10)
    
    p<-ggplot(enr_go_select, aes(x=reorder(Description, pvalNegLog10), y=pvalNegLog10)) +#, fill=in_grp
      geom_bar(stat="identity")+theme_minimal() + coord_flip() +
      labs(title= 'CICT Enriching top 100 TFs with most targets')+ #'CICT Enriching all then top 20 count')+ # 
      xlab('')+ylab('Negative Log 10 of p-value')+
      background_grid() +
      theme( legend.position="none",axis.text.y = element_text(angle = 45))+
      guides(colour=guide_legend(title=""));  p
    
    thepath = "E:/HDDs/Work/NYU Research/0 Plant biology/R code/results/GO terms enriched CICT2.tif"
    ggsave(filename = thepath, plot = p, width = 4, height = 3.5, device='tiff', dpi=300)
    
    
  }  
  
  #4 Top 20 predicted TFs Enrichment with their targets for the term "pluripotency"  stem cell population maintenance(GO:0019827)
  {
    rm(enr_go,enr_go_select,p)
    
    top30TFwithMostTrgts =n3d.e %>% group_by(Gene1) %>% 
      summarize(trgtCount =  sum(Gene2 %in% pp.symbs)) %>%
      arrange(desc(trgtCount)) %>% 
      mutate(in_grp= Gene1 %in% pp.symbs) %>%
      slice_head(n=30)
    
    pp.genesTargetingPP =n3d.e %>% group_by(Gene1) %>% 
      summarize(trgtCount = n()) %>% # sum(Gene2 %in% pp.symbs)) %>%
      arrange(desc(trgtCount)) %>% 
      mutate(in_grp= Gene1 %in% pp.symbs) %>%
      slice_head(n=30)
    
    predTFs.enriched =
      lapply(top30TFwithMostTrgts$Gene1 %>% head(20),
           function(predTF, desiredGOTermID){
            predTF.targets = n3d.e %>% arrange(desc(Weight)) %>% filter(Gene1 ==predTF)  %>% pull(Gene2) %>% unique()
            #n.itm.v$vID  #pp.genesTargetingPP$Gene1
            enr_go = enrichGO(predTF.targets,OrgDb=org.Mm.eg.db, keyType = "SYMBOL", ont = "BP",
                              pvalueCutoff = 1, pAdjustMethod = "BH", qvalueCutoff = 1, 
                              minGSSize = 10, maxGSSize = 500, readable = FALSE)
            
            enr_go_select = enr_go  %>%   as.data.frame() %>%   dplyr::filter(ID==desiredGOTermID)
            enr_go_select = enr_go_select %>% 
                               mutate(predictedTF = predTF, 
                                      predictedTargets = paste0(predTF.targets,collapse = ',')) %>%
                               select(predictedTF,everything() ,predictedTargets )
        },desiredGOTermID='GO:0019827')
    
    predTFs.enriched.d = rbindlist(predTFs.enriched) %>%
      mutate(pvalNegLog10 = -log10(pvalue), in_grp = predictedTF %in% pp.symbs)  %>% 
      dplyr::arrange(pvalNegLog10)
    
    p<-ggplot(predTFs.enriched.d, aes(x=reorder(predictedTF, pvalNegLog10), y=pvalNegLog10, fill=in_grp)) +#
      geom_bar(stat="identity")+theme_minimal() + coord_flip() +
      labs(title= 'CICT Enriching top 20 TFs using their targets for GO:0019827')+ #'CICT Enriching all then top 20 count')+ # 
      xlab('')+
      scale_y_continuous(name='- Log10 (p-value)', limits=c(0, 26),breaks = seq(0,30,5))+
      background_grid() +
      #theme( legend.position="none")+
      #guides(colour=guide_legend(title=""))+
    
      theme(axis.text.x = element_text( size = 8,angle=0,hjust=.2),
            legend.position = "none", # c(.9,.2),#legend.direction="horizontal",
            legend.text = element_text( size = 8,angle=0),
            legend.title = element_text( size = 8,angle=0),
            
      )+    scale_fill_discrete(name = '', labels =c( "False Positive", "True Positive"),
                            type  =c('#66CCEE','#EE6677'))+
        guides(fill = guide_legend(reverse=T)) ;  p
    #grpcolor1 = '#EE6677',grpcolor2 = '#66CCEE'
    thepath = "E:/HDDs/Work/NYU Research/0 Plant biology/R code/results/TFs enriched CICT.tif"
    ggsave(filename = thepath, plot = p, width = 2, height = 4, device='tiff', dpi=300)
    
    
  }   
  

  
  #1 Top 20 Random forest on relevance network Top 20 genes with most targets in 
  {
    
    n3d.v.rlv.genes = c(n3d.e.rlv$Gene1,n3d.e.rlv$Gene2) %>% unique()
    pp.genesTargetingPP.rlv =n3d.e.rlv %>%  group_by(Gene1) %>% 
      summarize(trgtCount = sum(Gene2 %in% pp.symbs)) %>%
      arrange(desc(trgtCount)) %>% 
      mutate(in_grp= Gene1 %in% pp.symbs) %>%
      slice_head(n=27)
    View(pp.genesTargetingPP.rlv)
    
    p.rlv<-ggplot(pp.genesTargetingPP.rlv, aes(x=reorder(Gene1, trgtCount), y=trgtCount, fill=in_grp)) +#
      geom_bar(stat="identity")+theme_minimal() + coord_flip()+
      xlab('')+ylab('Count of targets identified')+
      background_grid() +
      theme( legend.position="none")+
      guides(colour=guide_legend(title=""));  p.rlv
    
    thepath = "E:/HDDs/Work/NYU Research/0 Plant biology/R code/results/GO enriched perf Relevance.tif"
    ggsave(filename = thepath, plot = p.rlv, width = 3, height = 4, device='tiff', dpi=300)
  }
  
  #2 Enrichment for go term "pluripotency"  stem cell population maintenance(GO:0019827)
  {
    enr_go.rlv = enrichGO(n3d.v.rlv.genes,OrgDb=org.Mm.eg.db, keyType = "SYMBOL", ont = "BP",
                          pvalueCutoff = 1, pAdjustMethod = "BH", qvalueCutoff = 1, 
                          minGSSize = 10, maxGSSize = 500, readable = FALSE)
    enr_go_select.rlv = enr_go.rlv  %>% 
      as.data.frame() %>%
      dplyr::filter(ID=='GO:0019827')
  }
  
  #3 Enrichment for GO terms in RF network
  {
    rm(enr_go,enr_go_select,p)
    #pp.genesTargetingPP.rlv$Gene1 #n3d.e.rlv$Gene1
    enr_go = enrichGO( pp.genesTargetingPP.rlv$Gene1 , 
                       n,OrgDb=org.Mm.eg.db, keyType = "SYMBOL", ont = "BP",
                      pvalueCutoff = 1, pAdjustMethod = "BH", qvalueCutoff = 1, 
                      minGSSize = 10, maxGSSize = 500, readable = FALSE)
    
    enr_go_select = enr_go  %>% as.data.frame() %>% dplyr::arrange(desc(Count)) %>%
      slice_head(n=20) %>% mutate(pvalNegLog10 = -log10(pvalue)) %>% dplyr::arrange(pvalNegLog10)
    
    p<-ggplot(enr_go_select, aes(x=reorder(Description, pvalNegLog10), y=pvalNegLog10)) +#, fill=in_grp
      geom_bar(stat="identity")+theme_minimal() + coord_flip() +
      labs(title='RF Enriching top 100 TFs with most targets') + #'RF Enriching all then top 20 count')+ # 
      xlab('')+ylab('Negative Log 10 of p-value')+
      background_grid() +
      theme( legend.position="none",axis.text.y = element_text(angle = 45))+
      guides(colour=guide_legend(title=""));  p
    
    thepath = "E:/HDDs/Work/NYU Research/0 Plant biology/R code/results/GO terms enriched relev2_equalnet.tif"
    ggsave(filename = thepath, plot = p, width = 4, height = 3.5, device='tiff', dpi=300)
    
    
  } 
    
  #4  RF Top 20 predicted TFs Enrichment with their targets for the term "pluripotency"  stem cell population maintenance(GO:0019827)
  {
    rm(enr_go,enr_go_select,p,pp.genesTargetingPP.rlv)
    
    top30TFwithMostTrgts.rlv = n3d.e.rlv %>% group_by(Gene1) %>% 
      summarize(trgtCount = n()) %>% #sum(Gene2 %in% pp.symbs)
      arrange(desc(trgtCount)) %>% 
      mutate(in_grp= Gene1 %in% pp.symbs) %>%
      slice_head(n=30)
    
    pp.genesTargetingPP.rlv =n3d.e.rlv %>% group_by(Gene1) %>% 
      summarize(trgtCount =sum(Gene2 %in% pp.symbs)) %>% #
      arrange(desc(trgtCount)) %>% 
      mutate(in_grp= Gene1 %in% pp.symbs) %>%
      slice_head(n=30)
    
    predTFs.rlv.enriched =
      lapply(top30TFwithMostTrgts.rlv$Gene1%>% head(20),
             function(predTF, desiredGOTermID){
               predTF.targets = n3d.e.rlv %>% arrange(desc(Weight)) %>% filter(Gene1 ==predTF)  %>% pull(Gene2) %>% unique()
               #n.itm.v$vID  #pp.genesTargetingPP.rlv$Gene1
               enr_go = enrichGO(predTF.targets,OrgDb=org.Mm.eg.db, keyType = "SYMBOL", ont = "BP",
                                 pvalueCutoff = 1, pAdjustMethod = "BH", qvalueCutoff = 1, 
                                 minGSSize = 10, maxGSSize = 500, readable = FALSE)
               
               if(is.null(enr_go)) 
                 return (data.frame(predictedTF = predTF, pvalue=1,
                          predictedTargets = paste0(predTF.targets,collapse = ',')))
               
               enr_go_select = enr_go  %>%   as.data.frame() %>%   dplyr::filter(ID==desiredGOTermID)
               enr_go_select = enr_go_select %>% 
                 mutate(predictedTF = predTF, 
                        predictedTargets = paste0(predTF.targets,collapse = ',')) %>%
                 select(predictedTF,everything() ,predictedTargets )
               enr_go_select
               
             },desiredGOTermID='GO:0019827')
    
    predTFs.rlv.enriched.d = rbindlist(predTFs.rlv.enriched,fill=T) %>%
      mutate(pvalNegLog10 = -log10(pvalue), in_grp = predictedTF %in% pp.symbs) %>% 
      dplyr::arrange(pvalNegLog10)
    
    p<-ggplot(predTFs.rlv.enriched.d, aes(x=reorder(predictedTF, pvalNegLog10), y=pvalNegLog10, fill=in_grp)) +#
      geom_bar(stat="identity")+theme_minimal() + coord_flip() +
      labs(title= 'RF relevance Enriching top 20 TFs using their targets for GO:0019827')+ #'CICT Enriching all then top 20 count')+ # 
      xlab('')+
      scale_y_continuous(name='- Log10 (p-value)', limits=c(0, 26),breaks = seq(0,30,5))+
      scale_fill_discrete(name = '', labels =c( "False Positive", "True Positive"),
                               type  =c('#66CCEE','#EE6677'))+
      background_grid() +
      theme( legend.position="none")+
      guides(colour=guide_legend(title=""));  p
    
    thepath = "E:/HDDs/Work/NYU Research/0 Plant biology/R code/results/TFs enriched RF relevance equivalent network.tif"
    ggsave(filename = thepath, plot = p, width = 2, height = 4, device='tiff', dpi=300)
    
    
  }   
    
  
  #Comparision
  {
    library(clusterProfiler)
    library(org.Mm.eg.db)
    tidy_name <- function(name, n_char) {
      ifelse(nchar(name) > (n_char - 2), 
             {substr(name, 1, n_char) %>% paste0(., "..")},
             name)
    }
    
    gene_list = list('CICT'=pp.genesTargetingPP$Gene1, 
                     'RL'=pp.genesTargetingPP.rlv$Gene1,
                     'CICT-GO:0019827'= top30TFwithMostTrgts$Gene1,
                     'RL-GO:0019827'= top30TFwithMostTrgts.rlv$Gene1,
                     top30TFwithMostTrgts.rlv)
    CompareGO_BP=compareCluster(gene_list, fun="enrichGO", OrgDb=org.Mm.eg.db, 
                                keyType = "SYMBOL", ont = "BP",
                                pvalueCutoff = 0.1, pAdjustMethod = "BH", qvalueCutoff = 0.01,
                                minGSSize = 10, maxGSSize = 500, readable = FALSE)
    
    dp = dotplot(CompareGO_BP, font.size = 14,
                 showCategory=8, label_format = 80,
                 title="GO - Biological Process") 
          #coord_trans(y = "log10")
    
    dp$labels$colour = "- Log10(p.adjust)"
    dp$data$p.adjust = -log10(dp$data$p.adjust)
    
    dp = dp+ scale_color_gradient2(
          low = 'darkblue', # '#66CCEE', 
          mid = '#66CCEE',# "#9977EE", 
          high = 'red' , # '#EE6677',
          midpoint =3)
    
    #dp + scico::scale_color_scico(palette = "vik", midpoint=3) 
    #dp+scale_color_viridis()
    
    
    
    # dp+scale_colour_gradientn(
    #   values = c(0,.2,1),
    #   colors = c('darkblue','white','red')
    # )
    
    
    #dp+  scale_y_discrete(labels = stringr::str_wrap(Description, 15))
    
    
    dp +theme(
      binwidth = .7,
      #axis.text = element_text(hjust=0.5, vjust = 1),
      plot.margin = margin(0, 0, 0, 2, 'cm')) +
      coord_cartesian(clip = "off")
      
    
    #scale_y_discrete(labels = abbreviate)
      #stringr::str_wrap(V1, 15) ;dp
    # dp+background_grid() +
    #   theme( font_size = 8)
    
    thepath = "E:/HDDs/Work/NYU Research/0 Plant biology/R code/results/CICT RF cluster comparision RF equinet3.tif"
    ggsave(filename = thepath, plot = dp, width = 16, height =8, device='tiff', dpi=300)
    
  }
  
  gene.symbols <- n.itm.v$vID # cl.nodes %>% dplyr::filter(keys %in% names(unlist(gc2)))
  
  gene.ids.entrz<- AnnotationDbi::mapIds(org.Mm.eg.db, gene.symbols, 'ENTREZID', 'ALIAS')
  gene.ids.entrz=gene.ids.entrz %>%na.omit()
  
  gc1= data.frame(gene.id= n.itm.v$vID , key=names(gene.ids.entrz)) # %>% 
    #merge(cl.nodes, by.x = 'key', by.y='key') 
  
  gc2= gc1%>% #filter(key %in% names(gene.ids.entrz)) %>%
    group_by(type) %>% summarise(vec.id= list(gene.id),
                                 vec.symb= list(key))
  
  names(gc2$vec.id)<- paste0('X',gc2$type)
  names(gc2$vec.symb)<- paste0('X',gc2$type)
  
  
  ck <- compareCluster(geneCluster = gc2$vec.id, fun = enrichGO, # fun = enrichKEGG, # 
                       OrgDb =org.Mm.eg.db)
  
  ck <- compareCluster(geneCluster = gc2$vec.id, fun = enrichKEGG)# OrgDb =org.Mm.eg.db)
  dotplot(ck)
  
  ck <- compareCluster(geneCluster = gc2$vec.id, org,fun = enrichWP, 
                       organism = 'hsp',universe = gc1$gene.id)# OrgDb =org.Mm.eg.db)
  dotplot(ck)
  
  
  ck.ego = lapply(gc2$vec.id, function(x){
    go <- enrichGO(gene = gc2$vec.id$X1, OrgDb =org.Hs.eg.db,# OrgDb =org.Mm.eg.db)
                   keyType = 'ENTREZID' ,ont= 'MF' ,
                   universe = gc1$gene.id,
                   readable=TRUE);go
    go%>% as.data.frame() %>% filter(Count>2)
  })
  ck.ego.d = rbindlist(ck.ego,use.names=T, fill=T, idcol='cluster')
  
  kegg <- enrichKEGG(gene = gc2$vec.id$X8, organism = 'hsa',# OrgDb =org.Mm.eg.db)
                     universe = gc1$gene.id);kegg
  
  go%>% as.data.frame() %>% filter(Count>2)
  
  #ck <- setReadable(ck, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
  

  
  
  
  dotplot(formula_res)
  library(GSEAplot)
  
  
  
  library(ReactomeContentService4R)
  #BiocManager::install("Clusterprofiler")
  ## Load GOstats library
  library(GOstats); library(GO.db)
  ## Print complete GO term information for "GO:0003700"
  GOTERM$"GO:0003700";    GOMFPARENTS$"GO:0003700"; GOMFCHILDREN$"GO:0003700"
  ## Print number of GO terms in each of the 3 ontologies
  zz <- eapply(GOTERM, function(x) x@Ontology); table(unlist(zz))
  ## Gene to GO mappings for an organism (here Arabidopsis)
  library(org.Mm.eg.db) # For human use org.Hs.eg.db
  xx <- as.list(org.At.tairGO2ALLTAIRS)
  
  #Pathway DBs
  #KEGG
  if(F){
    ## Define function to create KEGG pathway list db 
    load_keggList <- function(org="ath") {
      suppressMessages(suppressWarnings(library(KEGG.db))) 
      kegg_gene_list <- as.list(KEGGPATHID2EXTID) # All organisms in kegg
      kegg_gene_list <- kegg_gene_list[grepl(org, names(kegg_gene_list))] # Only human
      kegg_name_list <- unlist(as.list(KEGGPATHID2NAME)) # All organisms in kegg
      kegg_name_list <- kegg_name_list[gsub(paste0("^", org), "", names(kegg_gene_list))]
      names(kegg_gene_list) <- paste0(names(kegg_gene_list), " (", names(kegg_name_list), ") - ", kegg_name_list)
      return(kegg_gene_list)
    }
    ## Usage:
    keggdb <- load_keggList(org="ath") # org can be: hsa, ath, dme, mmu, ... 
  }
  
  #reactome db
  if(F){## Define function to create Reactome pathway list db 
    load_reacList <- function(org="Mus") { #Mus for mouse
      library(reactome.db)
      reac_gene_list <- as.list(reactomePATHID2EXTID) # All organisms in reactome
      reac_gene_list <- reac_gene_list[grepl(org, names(reac_gene_list))] # Only human
      reac_name_list <- unlist(as.list(reactomePATHID2NAME)) # All organisms in reactome
      reac_name_list <- reac_name_list[names(reac_gene_list)]
      names(reac_gene_list) <- paste0(names(reac_gene_list), " (", names(reac_name_list), ") - ", gsub("^.*: ", "", reac_name_list))
      return(reac_gene_list)
    }
    ## Usage:
    reacdb <- load_reacList(org="Mus")
  }
}




#supplementary data
