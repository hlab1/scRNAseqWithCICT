#Functions and libs
{
  library('forestploter')
  library(ggplot2)
  library(ggparty)
  library(tidyverse)
  library(data.table)
  library(cowplot)
  library(dplyr)
  library(knitr)
  
  myplotViolin_old<-function(theTitle,dt,variables, group1, group2,
                         grp1title,grp2title,grpcolor1 = "green",grpcolor2 = "blue", 
                         yval = 'value',thebox=F,...)
  {
    #dataframe in long format
    dt=setDT(dt)[name %in% variables,]
    
    fig <- setDT(dt) %>%
      plot_ly(type = 'violin',...) 
    
    fig <- fig %>%
      add_trace(
        x = ~dt[classReg,get('name')],
        y = ~dt[classReg,get(yval)],
        legendgroup = grp1title,
        scalegroup = grp1title,
        name = grp1title,
        side = 'negative',
        box = list(
          visible = FALSE #thebox
        ),
      color = I(grpcolor1),
        
      meanline = list(
          visible = TRUE,
          line = list(color = 'black', width = 3, dash = 'dash')
        )
      ) 
    fig <- fig %>%
      add_trace(
        x = ~dt[classRnd,get('name')],
        y = ~dt[classRnd,get(yval)],
        legendgroup = grp2title,
        scalegroup = grp2title,
        name = grp2title,
        side = 'positive',
        box = list(
          visible = FALSE
        ),
        meanline = list(
          visible = TRUE,
          line = list(color = 'black', width = 3, dash = 'dash')
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
    
    fig <- fig 
    
    fig = fig %>% config(displayModeBar = F, showTips = F)
    fig
  }
  
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
        x = ~dt[classReg,get('name')],
        y = ~dt[classReg,get(yval)],
        legendgroup = grp1title,
        scalegroup = grp1title,
        name = grp1title,
        side = 'negative',
        box = list(
          visible = FALSE #thebox
        ),
        color = I(grpcolor1),
        scalemode = "count",
        
        meanline = list(
          visible = TRUE,
          line = list(color = 'black', width = 3, dash = 'dash')
        )
      ) 
    fig <- fig %>%
      add_trace(
        x = ~dt[classRnd,get('name')],
        y = ~dt[classRnd,get(yval)],
        legendgroup = grp2title,
        scalegroup = grp2title,
        name = grp2title,
        side = 'positive',
        box = list(
          visible = FALSE
        ),
        meanline = list(
          visible = TRUE,
          line = list(color = 'black', width = 3, dash = 'dash')
        ),
        color = I(grpcolor2),
        scalemode = "count"
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
    
    fig <- fig 
    
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


  library(dplyr)
  mutate <- dplyr::mutate
  select <- dplyr::select
  filter <- dplyr::filter
  rm(r.others.dt,r.multirun.dt,r.edgeType.dt, r.sparsity.dt)
  desiredOutcome = 'unsn_std_pPRC' #'unsn_pPRC' #    'unsn_PRC' # 
  desiredRandomOutcome = 'rndm_std_pPRC' #'rndm_pPRC' # 
  PSEUDO_Zero.std_pPRC=0.00001 #Maxmimum observed value for the network
  basePath = 'E:/HDDs/Work/NYU Research/0 Plant biology/CICT paper/SensitivityRound3/'

#Forest plot
{
  
  #A r.modeling_choices
  {
    r.modeling_choices = read.csv(paste0(basePath,'sens_modeling_choices.csv'))
    r.modeling_choices.bestEdges = r.modeling_choices %>% 
      group_by(dataset.dir,edgeType)  %>%
      dplyr::summarise(edgePerf = mean(unsn_pPRC,na.rm=T) ) %>% 
      group_by(dataset.dir) %>%
      arrange(desc(edgePerf)) %>% slice_head(n=3) #Pearson is the best overall
    
    #r.modeling_choices.gtsize
    {
      rm(r.modeling_choices.gtsize.raw);
      r.modeling_choices.gtsize.raw = r.modeling_choices %>%  setDF() %>%
        dplyr::mutate(outcome = .data[[desiredOutcome]],
               subgroup = paste0("    " ,  size.groundTruth)) %>% #size.groundTruth
        dplyr::filter(
          !is.na(outcome) &
          dataset != 'mDC' &  #Very low density, outlier
          trainingTrgt == 'class2' 
          #RF_ntree >=20,
          #RF_max_depth >=10,
          #size.groundTruth == 250 &
          #edgeType == 'Pearson'
          ) %>%
        dplyr::select(outcome,everything())
      
      r.modeling_choices.gtsize = r.modeling_choices.gtsize.raw %>%
        dplyr::group_by(subgroup)  %>%
        sample_n(size = 400) %>%
        #mutate(outcome = unsn_pROC/rndm_pPRC) %>%
        dplyr::mutate(rRatio= .data[[desiredOutcome]]/
                 ifelse(.data[[desiredRandomOutcome]]==0,density.net,
                        #max(PSEUDO_Zero.std_pPRC,min(.data[[desiredRandomOutcome]],na.rm=T)),
                        .data[[desiredRandomOutcome]])) %>%
        
        dplyr::summarize(outcome.se = sd(outcome,na.rm=T) / sqrt(n()),
                  est = mean(outcome,na.rm = T),
                  low = ifelse(est - 1.96*outcome.se<0, 0,est - 1.96*outcome.se),
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
        dplyr::filter(!is.na(est)) %>% arrange(desc(est))
      
      
      #N. TP edges in learning set
      r.titlerow = data.frame(subgroup = 'A. Size of learning set', est=NA,low=NA,hi=NA,se = NA, 
                              outcome.ci=' ',`Standardized.partial.AUPR`= ' ',`Samples(N)`='', description=' ')
      r.modeling_choices.gtsize = rbindlist(list(r.titlerow,r.modeling_choices.gtsize),fill=T)
      
      r.modeling_choices.gtsize.diag = rbind(
        r.modeling_choices.gtsize.raw %>% arrange(outcome) %>% slice_head(n=2),
        r.modeling_choices.gtsize.raw %>% arrange(outcome) %>% slice_tail(n=2)
      )
      rm(r.modeling_choices.gtsize.raw,r.modeling_choices.gtsize.diag)
      }
    

  }
  
  #B.r.modeling_choices.trainingTrgt.raw  Directed vs. undirected edges
  {
    r.modeling_choices.trainingTrgt.raw = r.modeling_choices %>% 
      dplyr::mutate(outcome = .data[[desiredOutcome]], #outcome = unsn_pROC/rndm_pPRC
             subgroup =paste0("    " ,  trainingTrgt)) %>% #size.groundTruth
      dplyr::filter(
        !is.na(outcome) &
        dataset != 'mDC' ,  #Very low density, outlier
          #trainingTrgt == 'class2' &
          RF_ntree >=20,
          RF_max_depth >=10,
          size.groundTruth == 250 &
          edgeType == 'Spearman') %>%
      dplyr::select(outcome,everything())
    
    
    r.modeling_choices.trainingTrgt = r.modeling_choices.trainingTrgt.raw %>%
      dplyr::group_by(subgroup)  %>%
      sample_n(size = 230) %>%
      dplyr::mutate(rRatio= .data[[desiredOutcome]]/
               ifelse(.data[[desiredRandomOutcome]]==0,density.net,
                      # max(PSEUDO_Zero.std_pPRC,min(.data[[desiredRandomOutcome]],na.rm=T)),
                      .data[[desiredRandomOutcome]])) %>%
      
      dplyr::summarize(outcome.se = sd(outcome,na.rm=T) / sqrt(n()),
                est = mean(outcome,na.rm = T),
                low = ifelse(est - 1.96*outcome.se<0, 0,est - 1.96*outcome.se),
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
      dplyr::filter(!is.na(est)) %>% arrange(desc(est))
    
    r.modeling_choices.trainingTrgt$subgroup =   plyr::mapvalues(r.modeling_choices.trainingTrgt$subgroup,
                                                                 c('    class2','    class3'),c('    Directed','    Undirected'))
    
    #Training with directed vs. undirected edges
    r.titlerow = data.frame(subgroup = 'B. Directed vs. undirected edges', est=NA,low=NA,hi=NA,se = NA,
                            outcome.ci=' ',`Standardized.partial.AUPR`= ' ',`Samples(N)`='', description=' ')
    r.modeling_choices.trainingTrgt = rbindlist(list(r.titlerow,r.modeling_choices.trainingTrgt),fill=T)
    
    r.modeling_choices.trainingTrgt.diag = rbind(
      r.modeling_choices.trainingTrgt.raw %>% arrange(outcome) %>% slice_head(n=2),
      r.modeling_choices.trainingTrgt.raw %>% arrange(outcome) %>% slice_tail(n=2)
    )
    
    rm(r.modeling_choices.trainingTrgt.raw,r.modeling_choices.trainingTrgt.diag)
  } 
  
  #C r.gtType 
  {
    r.gtType = read.csv(paste0(basePath,'sens_edgeType.csv'))
    glimpse(r.gtType)
    
    r.gtType.dt.raw = r.gtType %>% 
      dplyr::mutate(#outcome = unsn_pROC/rndm_pPRC,
        edges.total=as.numeric(str_replace_all(edges.total,',','')),
        outcome = .data[[desiredOutcome]],
        subgroup =paste0("    " , dataset.dir )) %>% #size.groundTruth  cictRawEdgeCol
      dplyr::filter( !is.na(outcome) &
                dataset != 'mDC' &  #Very low density, outlier
                str_detect(edgeType , 'Pearson|Kendall|ewMImm')) 
    #filter(str_detect(dataset.dir,'lofgof|_ns',negate = T))%>%
    
    r.gtType.dt = r.gtType.dt.raw %>%
      dplyr::group_by(subgroup)  %>%
      dplyr::mutate(rRatio= .data[[desiredOutcome]]/
               ifelse(.data[[desiredRandomOutcome]]==0,density.net,
                      #max(PSEUDO_Zero.std_pPRC,min(.data[[desiredRandomOutcome]],na.rm=T)),
                      .data[[desiredRandomOutcome]])) %>%
      
      dplyr::summarize(outcome.se = sd(outcome,na.rm=T) / sqrt(n()),
                est = mean(outcome,na.rm = T),
                low = ifelse(est - 1.96*outcome.se<0, 0,est - 1.96*outcome.se),
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
    
    r.gtType.dt[,"subgroup"] <- c('    L2 cell-specific ChIPseq','    L2 lofgof', '    L2 non-specific ChIPseq')
    r.titlerow = data.frame(subgroup = 'C. Type of gold standard', est=NA,low=NA,hi=NA,se = NA,
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
  
  #D r.edgeType
  {
    r.edgeType = read.csv(paste0(basePath,'sens_edgeType.csv'))
    glimpse(r.edgeType)
    
    r.edgeType.dt.raw = r.edgeType %>% 
      dplyr::mutate(#outcome = unsn_pROC/rndm_pPRC,
        outcome = .data[[desiredOutcome]],
        subgroup =paste0("    " , edgeType )) %>% #size.groundTruth  cictRawEdgeCol
      dplyr::filter( !is.na(outcome) &
                dataset != 'mDC' &  #Very low density, outlier
                str_detect(dataset.dir , 'L2$')) 
    #filter(str_detect(dataset.dir,'lofgof|_ns',negate = T))%>%
    
    r.edgeType.dt = r.edgeType.dt.raw %>%
      dplyr::group_by(subgroup)  %>%
      sample_n(size = 100) %>%
      dplyr::mutate(rRatio= .data[[desiredOutcome]]/
               ifelse(.data[[desiredRandomOutcome]]==0,density.net,
                      # max(PSEUDO_Zero.std_pPRC,min(.data[[desiredRandomOutcome]],na.rm=T)),
                      .data[[desiredRandomOutcome]])) %>%
      
      dplyr::summarize(outcome.se = sd(outcome,na.rm=T) / sqrt(n()),
                est = mean(outcome,na.rm = T),
                low = ifelse(est - 1.96*outcome.se<0, 0,est - 1.96*outcome.se),
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
      dplyr::filter(!is.na(est)) %>% arrange(desc(est))
    
    
    r.titlerow = data.frame(subgroup = 'D. Edge types', est=NA,low=NA,hi=NA,se = NA, #L2 cell specific ChIPseq, directed edges
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
                  low = ifelse(est - 1.96*outcome.se<0, 0,est - 1.96*outcome.se),
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
  
  #E r.multirun
  {
    r.multirun = read.csv(paste0(basePath,'sens_multipleRuns.csv'))
    glimpse(r.multirun); table(r.multirun$edgeType)
    
    r.multirun.dt.raw = r.multirun %>% 
      dplyr::filter(str_detect(expLevel, '^L2.*')) %>%
      dplyr::mutate(
        edges.total=as.numeric(str_replace_all(edges.total,',','')),
        outcome = .data[[desiredOutcome]] ,
        subgroup =paste0("    " ,  dataset)) %>%
      dplyr::filter(!is.na(outcome) & !is.na(density.net)) %>%
      dplyr::select(outcome, edges.total,everything())
    
      # filter(edgeType == 'Pearson') %>%
    r.multirun.dt=r.multirun.dt.raw %>% 

      dplyr::group_by(subgroup)  %>%
      sample_n(size = 150) %>%
      #mutate(outcome = unsn_pROC/rndm_pPRC) %>%
      dplyr::mutate(rRatio= .data[[desiredOutcome]]/
                ifelse(.data[[desiredRandomOutcome]]==0,density.net,
                       #max(PSEUDO_Zero.std_pPRC,min(.data[[desiredRandomOutcome]],na.rm=T)),
                       .data[[desiredRandomOutcome]])) %>%
      
      dplyr::summarize(outcome.se = sd(outcome,na.rm=T) / sqrt(n()),
                est = mean(outcome,na.rm = T),
                low = ifelse(est - 1.96*outcome.se<0, 0,est - 1.96*outcome.se),
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
      dplyr::arrange(desc(est))
  
    
    r.titlerow = data.frame(subgroup = 'E. Learning set sampling', est='',low='',hi='',se = '',
                            outcome.ci=' ',`Standardized.partial.AUPR`= ' ',`Samples(N)`='', description=' ')
    r.multirun.dt = rbindlist(list(r.titlerow,r.multirun.dt),fill=T)
    
    rm(r.multirun.dt.raw)
  }


  #r.others Based on BEELINE results
  {
    
    bloa = read.csv(paste0(basePath,'beeline_eval_summary09_2022-05-31-metrics_alg_out.csv'))
    bloa$met %>% unique()
    bloa %>% glimpse()
    
    algo_list = c('CICT', 'DEEPDRIM','Inferelator.Prior','Inferelator.NoPrior','GENIE3','GRNBOOST2','LEAP','PIDC','PPCOR','SCODE','SCRIBE','SINCERITIES','SINGE','RANDOM')
    #colnames(bloa)[12:39] %>% paste0(collapse = "','")
    
    # bloa1 = bloa
    # bloa1 = bloa1 %>%  dplyr::rename(ds_dir = `ï..ds_dir`) 
    # bloa1 = bloa1 %>%   pivot_wider(id_cols=c('ds_dir', 'ds_name','met'), names_from = 'alg', values_from='eval_value')
    
    bloa2 = bloa1 %>% 
      dplyr::select(-matches('R_.*')) %>%
      dplyr::mutate(Ratio_type = paste0('r',met)) %>%
      dplyr::mutate(across(any_of(algo_list),
                    ~ ./RANDOM,.names = 'r_{.col}'))
    
     #spAUPRC4 
    ratio_list = paste0('r_',algo_list)
    bloa.L2 =bloa2 %>% 
      #select(-any_of(algo_list))%>%
      pivot_longer(any_of(c(algo_list,ratio_list)),names_to = 'algo',values_to = "outcome") %>%
      dplyr::mutate(algoType = ifelse(str_detect(algo,'r_'), 'rRatio', 'outcome')) %>%
      dplyr::filter(met=='spAUPRC4') %>% 
      dplyr::mutate(algo = str_replace_all(algo,"r_","")) %>%
      dplyr::filter(str_detect(ds_dir,'^L2') & !str_detect(ds_dir,'_ns')& !str_detect(ds_dir,'_lofgof') ) %>%
      dplyr::filter(met == 'spAUPRC4') %>%
      #dplyr::filter(algo != 'CICT') %>%
      pivot_wider(names_from = 'algoType' , values_from = 'outcome') %>%
      dplyr::mutate(subgroup = paste0('    ', algo) , algo=NULL) %>%
      group_by(subgroup) %>%
      dplyr::summarize(outcome.se = sd(outcome,na.rm=T) / sqrt(n()),
                       outcome.se = ifelse(is.na(outcome.se),0,outcome.se),
                       est = mean(outcome,na.rm = T),
                       low = ifelse(est - 1.96*outcome.se<0, 0,est - 1.96*outcome.se),
                       hi =  est + 1.96*outcome.se,
                       `Standardized.partial.AUPR` = paste0(rep(' ',60),collapse=''),
                       outcome.ci = sprintf("%.2f (%.2f to %.2f)", est, low, hi),
                       
                       rRatio.se =  sd(rRatio,na.rm=T) / sqrt(n()),
                       rRatio.se = ifelse(is.na(rRatio.se),0,rRatio.se),
                       rRatio.est = mean(rRatio,na.rm = T),
                       rRatio.low =  ifelse( rRatio.est - 1.96*rRatio.se<0, 0, rRatio.est - 1.96*rRatio.se), 
                       rRatio.hi = rRatio.est + 1.96*rRatio.se,
                       `rRatio.Standardized.partial.AUPR` = paste0(rep(' ',60),collapse=''),
                       rRatio.outcome.ci = sprintf("%.2f (%.2f to %.2f)",rRatio.est, rRatio.low,rRatio.hi),
                       
                       evalCount = n(),
                       description='') %>%
      filter(!is.na(est)) %>% arrange(desc(est))
    
    bloa.L2$subgroup = str_replace_all(bloa.L2$subgroup,'[.]','_')
    r.titlerow = data.frame(subgroup = 'F. Benchmarked algorithms', est=NA,low=NA,hi=NA,se = NA,
                            outcome.ci=' ',`Standardized.partial.AUPR`= ' ',`Samples(N)`='', description=' ')
    r.others.dt = rbindlist(list(r.titlerow,bloa.L2),fill=T)
    

    
  }
  
  #Mean and median reports on other datasets
  if(F){
    bloa.L2$outcome %>% hist()
    bloa.L2$algo %>% unique
    
    bloa.L2 %>% 
      dplyr::filter(str_detect(ds_dir,'SERGIO')) %>% 
      dplyr::filter(algo !='r_CICT') %>% 
      filter( Ratio_type=='rpAUPRC4') %>%
      pull(outcome) %>% summary()
    
    bloa.L2 %>% 
      dplyr::filter(str_detect(ds_dir,'SERGIO')) %>% 
      dplyr::filter(str_detect(algo,'CICT|r_Inferelator.Prior')) %>% 
      filter( Ratio_type=='rpAUPRC4') %>%
      pull(outcome) %>% summary()
    
    #BLC are only ratios
    #============================================ @@@@
    blc = read.csv(paste0(basePath,  'beeline_eval_summary09_2022-05-31-metrics_alg_out2.csv'))# 'beeline_eval_summary09_2022-05-31-metrics_alg_fold_wide.csv'))
    algo_list = colnames(blc)[11:25] #%>% paste0(collapse=",")
    names(blc)[1]<-"ds_dir"
    # blc1 = blc %>% 
    #   dplyr::select(-matches('R_.*')) %>%
    #   #dplyr::mutate(RANDOM = ifelse(is.na(RANDOM) | RANDOM==0, .0001, RANDOM)) %>%
    #   dplyr::mutate(across(any_of(algo_list),~ as.numeric(.))) %>%
    #   dplyr::mutate(across(any_of(algo_list),
    #                        ~ ./RANDOM,.names = 'r_{.col}')) %>%
    #   pivot_longer(any_of(c(algo_list)),names_to = 'algo',values_to = "outcome")
    
    summary_custom <- function(data,rndDigits = 2) {
      med = median(data, na.rm = TRUE)%>% round(rndDigits)
      rng = range(data, na.rm = TRUE) %>% round(rndDigits)
      #paste0(med, " (",paste0(rng,collapse=' - '),") " )
      paste0(paste0(rng,collapse=' - '), " with median ", med )
    }
    
    blc.l = blc %>%
      #filter(met =='rpAUPRC4') %>%
      dplyr::mutate(across(any_of(algo_list),~ as.numeric(.))) %>%
      pivot_longer(any_of(c(algo_list)),names_to = 'algo',values_to = "outcome")
    
    
    
    blc.l$met %>% unique()
    
    #CICT with dropout
    blc.l %>% 
      dplyr::filter(str_detect(ds_dir,'SERGIO')) %>% 
      filter( met== 'rpAUPRC4') %>% #'rEPr'
      filter(noise_desc =='With dropout') %>% #No dropout
      filter(str_detect(algo,'CICT')) %>%
      pull(outcome) %>% summary_custom()
    
    blc.l %>% 
      dplyr::filter(str_detect(ds_dir,'SERGIO')) %>% 
      filter( met== 'rpAUPRC4') %>% #'rEPr'
      filter(noise_desc =='With dropout') %>% #No dropout
      filter(!str_detect(algo,'CICT|Inferelator')) %>%
      pull(outcome) %>% summary_custom()
    
    a=blc.l %>% 
      dplyr::filter(str_detect(ds_dir,'SERGIO')) %>% 
      filter( met== 'rEPr-TF') %>% #'rEPr'
      filter(noise_desc =='With dropout') %>% #No dropout
      filter(!str_detect(algo,'CICT|Inferelator')) %>%
      pull(outcome) %>% summary_custom()
    
    bloa.L2 %>% 
      dplyr::filter(str_detect(ds_dir,'SERGIO')) %>% 
      filter(str_detect(algo,'CICT')) %>%
      filter(noise_desc =='With dropout') %>% #No dropout
      filter( Ratio_type=='rpAUPRC4') %>%
      pull(outcome) %>% summary()
    
    
    
  }
  
  library(RColorBrewer)
  clrs = RColorBrewer::brewer.pal(8,'YlGnBu')[c(1,2,3,4,5,6,7)] #   'Pastel2')[c(5,6,7)]
  clrs = c(clrs,gray.colors(12)[9])
  
  library(paletteer)
  #paletteer::paletteer_c()
  clrs = paletteer_c("ggthemes::Blue-Green Sequential", 7) 
  clrs=c(clrs[c(1:5)],"#F0F0F0")
  
  
  #Desired measure Forest graph
  {
    
    dt = rbind(r.modeling_choices.gtsize
               ,r.modeling_choices.trainingTrgt
               ,r.gtType.dt           
               ,r.edgeType.dt
    
               ,r.multirun.dt
               #,r.sparsity.dt
               
               ,r.others.dt
               ) %>% 
      dplyr::select(evalCount,everything()) %>%
      dplyr::mutate(low = ifelse(is.na(low),0,low),
             hi = ifelse(is.na(hi),0,hi),
             evalCount = ifelse(is.na(evalCount),'',evalCount),
             outcome.semult100 = as.integer(outcome.se *200) )%>% 
      #mutate(evalCount = str_pad(evalCount, width = 3, side = "left",pad="\t")) %>%
      dplyr::mutate(across(everything() & !subgroup, ~ ifelse(is.na(est),'',.) ))
    
    # dt= dt %>% mutate(across(everything(), ~ ifelse(is.na(.),'',.) )) #| !str_detect(.data[['subgroup']],'\t')
    # dt= dt %>% mutate(across(everything(), ~ ifelse(str_detect(.$subgroup,'\t'),'',.) ))
    
    #Forest graph
    {
      rm(p)
      p <- forest(dplyr::select(dt,subgroup,outcome.ci ,`Standardized.partial.AUPR`, evalCount),
                  est = dt$est,
                  lower = dt$low, 
                  upper = dt$hi,
                  sizes =  2,#dt$outcome.semult100 ,#
                  ci_column = 3,
                  ref_line = 0,
                  #arrow_lab = c("Placebo Better", "Treatment Better"),
                  xlim = c(0, 1.05),
                  ticks_at = c(0.1,0.25, .5,.75,1), #    c(0, 0.05, .1, .15, .3),
                  footnote = "footnote",
                  title = paste0(" Outcome: " , desiredOutcome)
                ); p
      
      
    
      # Edit background of row 5
      # Bold grouping text
      headingRows = which(str_detect(dt$subgroup,'^   ',negate = T))
      p <- edit_plot(p,
                     row = headingRows,
                     gp = gpar(fontface = "bold"))
      
      # headingRows = c(headingRows,nrow(dt)+1)
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
      
      # p <- edit_plot(p, row = subgrouprows[[6]], which = "background",
      #                gp = gpar(fill = clrs[6]));
      
      # p <- edit_plot(p, row = subgrouprows[[7]], which = "background",
      #                gp = gpar(fill = clrs[7]));
      # 
      
      thepath = paste0(basePath,"Forest plot stdpAPRC2.tif")
      #thepath = "E:/HDDs/Work/NYU Research/0 Plant biology/R code/results/Forest plot stdpAPRC2.tif"
      print(thepath)
      
      library(extrafont)
      #font_import(paths = "C:/Windows/Fonts")
      loadfonts()
      ggsave(filename = thepath, plot = p, width = 12, height = 14, dpi=300,
             device='tiff')#, family = "Arial")
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
               #,r.sparsity.dt
               ,r.others.dt
    ) %>% 
      dplyr::select(evalCount,everything()) %>%
      dplyr::mutate(rRatio.low = ifelse(is.na(rRatio.low),0,rRatio.low),
             rRatio.hi = ifelse(is.na(rRatio.hi),0,rRatio.hi),
             evalCount = ifelse(is.na(evalCount),'',evalCount),
             rRatio.outcome.ci = ifelse(is.na(rRatio.outcome.ci),'',rRatio.outcome.ci),
             outcome.semult100 = as.integer(outcome.se *200) ) 

    sd(dt$rRatio.est,na.rm = T)
    #Forest graph
    {
      rm(p)
      maxtick=max(dt$rRatio.est,na.rm = T)+10
      
      maxtick= max(as.numeric(dt$rRatio.hi),na.rm = T)%>% round(digits=0)+10
      
      p <- forest(dplyr::select(dt,subgroup,rRatio.outcome.ci ,`Standardized.partial.AUPR`, evalCount), #description
                  est = dt$rRatio.est,
                  lower = dt$rRatio.low, 
                  upper = dt$rRatio.hi,
                  sizes = 3, #dt$rRatio.outcome.se ,#
                  ci_column = 3,
                  ref_line = 0,
                  #arrow_lab = c("Placebo Better", "Treatment Better"),
                  xlim = c(0, maxtick),
                  ticks_at = c(0,  50,100,200,maxtick), #    c(0, 0.05, .1, .15, .3),
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
      
      # p <- edit_plot(p, row = subgrouprows[[6]], which = "background",
      #                gp = gpar(fill = clrs[6]));
      
      # p <- edit_plot(p, row = subgrouprows[[7]], which = "background",
      #                gp = gpar(fill = clrs[7]));
      # 
    }
    
    thepath = paste0(basePath,"Forest plot ratios2.tif")
    print(thepath)
    ggsave(filename = thepath, plot = p, width = 12, height = 14, device='tiff', dpi=300)
    p
  }
}

# Ground truth size versus RF_ntree,'_',RF_max_depth
{
  table(r.modeling_choices$RF_ntree)
  table(r.modeling_choices$RF_max_depth)
  table(r.modeling_choices$size.groundTruth)
  table(r.modeling_choices$edgeType)
  table(r.modeling_choices$dataset)
  table(r.modeling_choices$trainingTarget)
  table(r.modeling_choices$dataset)

  
  
  library(fmsb)
  library(ggradar)
  library(scales)
  library(knitr)
  library(dplyr)
  #size of ground truth
  r.meta =r.modeling_choices %>% group_by(dataset) %>%
    mutate(meanRandom = mean(.data[[desiredRandomOutcome]],na.rm=T)) %>% 
    ungroup() %>%
    rowwise() %>%     
    mutate(rRatio= .data[[desiredOutcome]]/
               ifelse(.data[[desiredRandomOutcome]]==0,
                      meanRandom,
                      .data[[desiredRandomOutcome]])) %>%
    #mutate(grp = paste0(RF_ntree,'_',RF_max_depth)) %>% 
    group_by(!!desiredOutcome,RF_max_depth,size.groundTruth,RF_ntree, add=T) %>%
    summarize(smrzRatio = mean(.data[[desiredOutcome]], na.rm=T)) 
  # %>%  ungroup() 
    #group_by(RF_ntree,RF_max_depth,size.groundTruth,grp) %>%
    #summarize(smrzRatio = mean(rRatio, na.rm=T))  %>%  ungroup() %>% select(-grp) %>%
  
  if(F){
    r.meta.wide = r.meta%>% #mutate(grp = paste0("g_",row_number())) %>% #select(-grp) %>%
    arrange(size.groundTruth,RF_ntree)  %>%
    pivot_wider(id_cols =c('size.groundTruth'), names_from = grp,values_from = smrzRatio) %>%
    #mutate_at(vars(-size.groundTruth), rescale) %>%
    mutate(size.groundTruth=as.character(size.groundTruth))
    if(F){
      # ggplot(r.meta, aes(x =RF_ntree , y = smrzRatio, color = factor(RF_max_depth))) +
      #   geom_line() +
      #   facet_wrap(~size.groundTruth) +
      #   labs(x = "Size Ground Truth", y = "smrzRatio") +
      #   theme_bw()
      
      ggplot(r.meta, aes(x = RF_ntree, y = smrzRatio, color = factor(RF_max_depth))) +
        geom_line(size = 1.2) +
        facet_wrap(~ size.groundTruth) +
        labs(x = "Size Ground Truth", y = "smrzRatio") +
        cowplot::theme_cowplot()+
        theme(
          legend.position = "top",
          #panel.grid.major = element_blank(),
          #panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(size = 0.5),
          axis.text = element_text(size = 8),
          axis.title = element_text(size = 10),
          strip.text = element_text(size = 10, face = "bold"),
          plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
        ) +
        
        guides(color = guide_legend(ncol = 2))
    }
  

  #datasets
  r.meta =r.modeling_choices %>%group_by(dataset) %>%
    mutate(meanRandom = mean(.data[[desiredRandomOutcome]],na.rm=T)) %>% 
    ungroup() %>%
    rowwise() %>%     
    mutate(rRatio= .data[[desiredOutcome]]/
             ifelse(.data[[desiredRandomOutcome]]==0,
                     meanRandom,#max(PSEUDO_Zero.std_pPRC,min(.data[[desiredRandomOutcome]],na.rm=T)),
                    .data[[desiredRandomOutcome]])) %>%
    #mutate(grp = paste0(RF_ntree,'_',RF_max_depth)) %>%
    group_by(dataset,RF_ntree,RF_max_depth) %>% 
    summarize(smrzRatio = mean(rRatio, na.rm=T))  %>%  
    ungroup() %>% mutate(grp = paste0('g_',row_number())) %>% 
    mutate_at(vars(-dataset,-grp), rescale) %>%
    pivot_wider(id_cols = dataset, names_from = grp,values_from = smrzRatio)
  }
  
r.meta %>%  kable()
 glimpse(r.meta)
  
 #geomline
 if(T){
   plt.d = ggplot(r.meta, aes(x =RF_ntree , y = smrzRatio, color = factor(RF_max_depth))) +
    geom_line(size = 1) +
    facet_wrap(~size.groundTruth,nrow=1) +
    labs(x = "Number of trees", y = "Standardized partial AUPR") +
     expand_limits(y = 0)+
    cowplot::theme_cowplot()+
     
    theme(
      legend.position = "top",legend.direction = "horizontal", 
      panel.grid.major = element_line(colour = "lightgray", linetype = "solid" ,size = .3),
      #panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(size = 0.5),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 8),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 10),
      strip.text = element_text(size = 8, face = "bold"),
      plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
    ) +
    guides(color = guide_legend(nrow = 1,title="Tree depth"));plt.d
   
   thepath = paste0(basePath, "Learningset size, tree depth, tree count.png")
   ggsave(filename = thepath, plot = plt.d, width = 8, height = 4, device='png', dpi=300)
  }


}

#r.multi.mESC_lofgof
{
  rm(r.others.dt,r.multirun.dt,r.edgeType.dt, r.sparsity.dt)
  desiredOutcome = 'unsn_std_pPRC' #'unsn_pPRC' #    'unsn_PRC' # 
  desiredRandomOutcome = 'rndm_std_pPRC' #'rndm_pPRC' # 
  PSEUDO_Zero.std_pPRC=0.00001 #Maxmimum observed value for the network
  basePath = 'E:/HDDs/Work/NYU Research/0 Plant biology/CICT paper/SensitivityRound2/'
  
  r.multi.mESC_lofgof = readRDS(paste0(basePath,'mESC_multiple_learnigsets.rds'))
  glimpse(r.multi.mESC_lofgof)
  

  r.multi.mESC_lofgof.raw = r.multi.mESC_lofgof %>% 
    mutate( idx = row_number(),
            edges.total=as.numeric(str_replace_all(edges.total,',','')),
            outcome = .data[[desiredOutcome]] ,
            rRatio= .data[[desiredOutcome]]/
              ifelse(.data[[desiredRandomOutcome]]==0,density.net,
                    #max(PSEUDO_Zero.std_pPRC,min(.data[[desiredRandomOutcome]],na.rm=T)),
                    .data[[desiredRandomOutcome]])) %>%
            
    filter(!is.na(outcome) & !is.na(density.net)) %>%
    select(outcome, edges.total,everything()) %>% arrange(desc(rRatio))
  
  r.multi.mESC_lofgof.f = r.multi.mESC_lofgof.raw %>% select(idx,rRatio,unsn_PRC,unsn_pPRC,unsn_std_pPRC)
  
  library(ggplot2)
  ggplot(r.multi.mESC_lofgof.raw, aes(x=reorder(factor(idx),rRatio), y=rRatio)) + 
    geom_bar(stat="identity", position="dodge")
  
  U}



# Density and violin plot ********************************************
{
  # Data preparation for distribution and violin plots 
  {
  
    library(ggplot2)
    library(zeallot)
    library("minet")
    library(h2o)
    #Load data *********************************************************************88
    url.rslts = 'E:/HDDs/Work/NYU Research/0 Plant biology/CICT paper/SensitivityRound3/'
    netObjs = readRDS(file = paste0( url.rslts,  '/mESC network objects.rds')) #'/mESC_rslt_3.rds')) #
    c(rcrd,prd.varimp,n.itm.e,pred_outcome.e,t1) %<-% netObjs #truncated_pred_outcome
    
    varimp = netObjs$varimp
    rm(netObjs)
    
    
    d.unbalanced = t1 %>% dplyr::select(-Weight) %>%
      inner_join(pred_outcome.e, by = c('src'='Gene1','trgt'='Gene2')) %>% 
      arrange(desc(Weight)) %>% as.data.table() 
    
    n.classCausal = sum(d.unbalanced$class2==TRUE)
    n.classRandom = n.classCausal *20
    setDT(d.unbalanced)
    d = d.unbalanced
    
    # Bin data ********************************************************************888
    nbins= 30 # as.integer(sqrt(nrow(t1)/2)/20)
    breaks =1:nbins
    
  
    d.binned =data.table(grpvar=c(rep(T,nbins-1),rep(F,nbins-1)))
    d.distDstnc = c()
  
    d.var = 'scfKurt.x' 
    ggp.vars = c("scbtau4.x" , "scfNKurt.x" ,"scbNMean.x", "scfL3.x" ,   "scbL1.x")
    ggp.vars = prd.varimp$variable[c(1,2,3,5,8)]
    
    #Find variables with best visual distinct distributions
    for(d.var in ggp.vars){ #prd.varimp$variable[1:80]){ # 
      print(d.var)
      edges.causal = myhist(infotheo::discretize(setDT(d)[ class2 == T, d.var,with=F] %>% 
                                                   unlist(), "equalwidth", nbins)$X,nbins,breaks,prob=F)
      edges.random =  myhist(infotheo::discretize(d[ class2 == F, d.var,with=F] %>% 
                                                    unlist(), "equalwidth", nbins)$X,nbins,breaks,prob=F)
      
      setDT(d.binned)[,(d.var):= c(edges.causal,edges.random)]
      
      tmp = data.frame(edges.causal,edges.random)
      mim <- build.mim(tmp, disc = 'equalwidth', estimator = 'kendall')
      dist = mim[1,2]; names(dist)<-d.var
      d.distDstnc = c(d.distDstnc,dist)
    }
    
    d.distDstnc=sort(desc(d.distDstnc))
    d.topDstnc = names(d.distDstnc)#[c(1,3,6,7,8)]# names(d.distDstnc)[c(1,3,4)]
    var =names(d.distDstnc)[3]
  }
    
  # Density plot *****************************************************************888
  {
    library(easyGgplot2)
    library(cowplot)
    {
      for(desiredVars in d.topDstnc)
      {
        
      d.binned.scaled = setDT(d.binned)[,c(desiredVars,'grpvar'),with=F] %>% setDF() %>%
        mutate(!!desiredVars:= log2(abs(.data[[desiredVars]])) )  #%>% pivot_longer(cols = desiredVars, names_to = 'var')
      
      rm(plt.d);plt.d = ggplot2.density(data=d.binned.scaled , xName=desiredVars  , groupName = 'grpvar' ,
                      alpha=0.7, fillGroupDensity=T,
                      scale = 'area',
                      addMeanLine=F, meanLineColor=NULL, meanLineSize=1,xlim=c(0,25)
                      )#+facet_wrap(desiredVars);
      
      plt.d= plt.d +theme_classic() +
        background_grid() + 
        theme(axis.text.x = element_text(size = 9, angle = 0, hjust = .2),
              legend.background = element_rect(fill = "transparent"),
              legend.position = c(.8, 0.98), legend.justification = c(1, 0.5),legend.direction = "horizontal",
              legend.text = element_text(size = 9, angle = 0),
              legend.title = element_text(size = 9, angle = 0),
              legend.key.size = unit(0.5, "cm"),
              plot.margin = margin(5, 5, 5, 5)) + 
        scale_fill_discrete(name = '', labels =c( "Random edges", "Regulatory edges"),
                            type  =c('#66CCEE','#EE6677'))+
        guides(fill = guide_legend(reverse=T))+ 
        xlab('Log2 L-Kurtosis  of source contrib.')+ 
        ylab("Density");plt.d #, 
      
        ggsave(filename = paste0(basePath,desiredVars,".tif"), plot = plt.d, width = 4, height = 2.5, device='tiff', dpi=300)
     }
      
    }
  }

  #Violin plot ***********************************************************************
  #Additional data preparation for violin plot
  {
    
    studyDataset = "mESC-lofgof"
    theTitle = '' #paste0(studyDataset, " - Dist. top predictors of causal and random edges")
    
    
    plt.violin.cols.1 =ggp.vars # c(names(d.distDstnc)[1:5]) #
    df=setDF(d.binned) %>%  #setDF(d) %>% mutate(grpvar = class2) %>%
      dplyr::select(ggp.vars,'grpvar') %>%
      dplyr::mutate(across(ggp.vars, ~ log2(abs(.)) )) %>%
      dplyr::mutate(eID = row_number()) %>%
      dplyr::mutate(classReg = grpvar , classRnd= !grpvar)  %>%
      dplyr::select(-grpvar)
    
    ggp.vars # "scbtau4.x"  "scfNKurt.x" "scbNMean.x" "scfL3.x"    "scbL1.x"

    ggp.vars.new =c('L-Kurtosis  of source contrib.',
                     'Normalized kurtosis of source conf.',
                     'Normalized mean of source contrib.',
                     'L-Skewness of source conf.',
                     'Median of source contrib.'
                    
                     )
    names(df)[1:5]<-ggp.vars.new

    rm(df.lng)
    df.lng = df %>% 
      dplyr::select(-classRnd,ggp.vars.new) %>%
      pivot_longer(ggp.vars.new )
        #mutate(value = log2(value+PSEUDO_Zero.std_pPRC))
      

    
    #The most useful of these are tau _{3}called the L-skewness, and tau _{4} the L-kurtosis. 
    plt.violin.cols.2 = names(df)[1:length(plt.violin.cols.1)]
    PSEUDO_Zero.std_pPRC  = 1e-5
    dfl = df %>%  pivot_longer(plt.violin.cols.2) %>% # setdiff(colnames(df),c('grpvar','eID')))  
      dplyr::mutate(logval = sign(value) * log2(abs(value)+PSEUDO_Zero.std_pPRC)) 


  }

  #GGplot2 Violin plot
  {
    #devtools::install_github("psyteachr/introdataviz")
    ggp.vars 
    ggp.vars.new.1 = ggp.vars.new[c(1,2,4,5,3)]
    df.lng$name <- factor(df.lng$name, levels = ggp.vars.new.1)
    #df.lng$classReg <- factor(df.lng$classReg, levels = c(FALSE, TRUE))

      #geom_boxplot(width=.1, outlier.colour=NA, position = position_dodge()) 
    plt.d=ggplot(df.lng, aes(name, value,fill = as.factor(classReg))) +
      geom_violin(position = position_dodge(width = 0.3), scale = 'area', orientation = "x",trim=F)+
      theme_half_open() +
      background_grid(size.major = 0.25) + 
        theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
              legend.position = "top",
              legend.background = element_rect(fill = "transparent"),
              #legend.position = c(.95, 0.95), legend.justification = c(1, 0.5),
              # legend.text = element_text(size = 12, angle = 0),
              # legend.title = element_text(size = 12, angle = 0),
              text = element_text(size=12,family='Arial'),
              
              plot.margin = margin(t = 0,  # Top margin
                                   #r = 1,  # Right margin
                                   b = 0,  # Bottom margin
                                   l = 50,  # Left margin
                                   unit = "pt")) + 
        scale_fill_discrete(name = '', labels =c(  "Regulatory edges","Random edges"),
                            type  =c('#EE6677','#66CCEE'))+
        guides(fill = guide_legend(reverse=F))+ 
        ylab("Log2 (frequency of binned values)") ;plt.d
    
    
   

    thepath = paste0(url.rslts, "violinTop2.png")
    ggsave(filename = thepath, plot = plt.d, width = 5, height = 5,  dpi=300)
 
    
    
  }


}



#Random forest TF-aware and not TF-Aware
{
  library(khroma)
  url.rslts = "E:/HDDs/Work/NYU Research/0 Plant biology/R code/results/L2_supervised_eval.rds"
  rcrd = readRDS(url.rslts)
  
  rcrd1 = rcrd %>% select(-rAUPR,-r.file) %>% 
    pivot_wider(id_cols = c('groundTruth','data','network.density','informed'),
                values_from = 'rpAUPR',
                names_from = 'informed') #                names_glue = "rpAUPR_{informed}" ) 
  rcrd1 = rcrd1 %>% mutate(raiseRatio = `TF-aware`/`Not-TF-aware`) 
  rcrd1 %>% knitr::kable()
  
  thepath = "E:/HDDs/Work/NYU Research/0 Plant biology/R code/results/L2_supervised TF aware summary.csv"
  write.csv(rcrd1,thepath)
  
  tol_vib7 = khroma::color('vibrant')(7)
  alg = c(tol_vib7[['orange']],tol_vib7[['teal']])
  
  #bar chart
  {
    agg_data=data.frame(
      run = as.character(1:9),
      rnd_nottfaware = c(30.11 , 121.22 , 139.07 , 141.31 , 45.56 , 165.80 , 180.67 , 206.16 , 212.14),
      rnd_tfaware =c(0.63 , 2.64 , 2.85 , 3 , 3.07 , 3.39 , 3.85 , 4.33 , 4.35 )
    )
    
     # reshape the data
    melt_data <- agg_data %>% 
      pivot_longer(cols = c(rnd_nottfaware, rnd_tfaware),
                   names_to = "variable",
                   values_to = "value")
    
    # plot the data using ggplot2
    p=ggplot(melt_data, aes(x = run, y = value, fill = variable)) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
      labs(x="Runs", y="rpAUPR") +
      cowplot::theme_cowplot()+
      scale_fill_manual(name = "", 
                        values = c("rnd_nottfaware" = alg[1], 
                                   "rnd_tfaware" = alg[2]), 
                        labels = c("Non-TF-informed", "TF-informed"))+
    
      theme(legend.position = "top",
            axis.text = element_text(size = 8),
            legend.text = element_text(size = 8, angle = 0),
            legend.title = element_text(size = 8, angle = 0),
            axis.title = element_text(size = 10)) ;p
    
    thepath = "E:/HDDs/Work/NYU Research/0 Plant biology/R code/results/Random Classifier comparision TFaware VS not 1.tif"
    ggsave(filename = thepath, plot = p, width =4, height = 3, device='tiff', dpi=300)
    
  }
  
  #Box plot
  if(F){
    df_long <- rcrd1 %>%
      dplyr::rename(`Non-TF-informed` = `Not-TF-aware`, `TF-informed`=`TF-aware`) %>%
      pivot_longer(cols = c("Non-TF-informed", "TF-informed"), 
                   names_to = "variable", 
                   values_to = "value") %>%
      select(-data, -network.density, -raiseRatio) %>%
      filter(value < median(value) + 5*IQR(value) & value > median(value) -5*IQR(value)) %>%
      mutate(groundTruth = str_replace_all(groundTruth,'^L2$','L2_cs'))
      
      
    
    
    p=ggplot(df_long, aes(x=variable, y=value, color=variable)) +
      geom_boxplot() +
      stat_boxplot(geom ='errorbar',coef=1000) +
      geom_point(position=position_jitterdodge(jitter.width=0.2,dodge.width=0.75),size=.7,alpha=0.8) +
      labs(x="", y="rpAUPR") +
      scale_color_manual(values=alg, labels=c("Non-TF-informed", "TF-informed")) +
      facet_wrap(~groundTruth, scales = "free_y", ncol=1) +
      theme_cowplot() +
      coord_flip() +
      theme(panel.grid.major = element_line(colour = "lightgray", linetype = "solid" ,size = .3)) +
      theme(legend.position = "top", legend.direction = "horizontal",
            axis.text = element_text(size = 8),
            legend.text = element_text(size = 8),
            strip.text = element_text(size = 8),
            axis.title = element_text(size = 10)
            
      ) +
      guides(color = guide_legend(title = "", nrow = 1));p
    
    thepath = "E:/HDDs/Work/NYU Research/0 Plant biology/R code/results/RF comparision TFaware VS not 1.tif"
    ggsave(filename = thepath, plot = p, width =4, height = 3, device='tiff', dpi=300)
    
    
    }

   
  
}

#Network rendering 
{
  library(zeallot)
  library(SummarizedExperiment)
  require(verification)
  require(pROC)
  require(ROCR)
  require( OptimalCutpoints)
  require(precrec )
  library(dplyr)
  library(magrittr)
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
  # netObjs = readRDS(file = paste0( basePath, 'mESC network objects.rds'))
  # c(rcrd,prd.varimp,n.itm.e,pred_outcome.e,truncated_pred_outcome,t1) %<-% netObjs
  # n.itm.e = netObjs[[3]];prd.varimp=netObjs[[2]]
  url.rslts = paste0(basePath,'mESC_rslt_.rds')
  rcrd = readRDS(url.rslts)
  # 
  # rcrd = netObjs[[1]]
  # n.itm.e=netObjs[[3]]
  # pred_outcome.e=netObjs[[4]]
  

  
  
  thepath = basePath
  tbl.goldStandard = read.csv(paste0(basePath,'mESC L2_lofgof-refnetwork.csv'))
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
    

    
    if(T){
      relativeCostfn_fp =1/2 # 1/5
      #find best cutoff
      prv = table(pred_outcome$outcomes)
      best.weights=c(relativeCostfn_fp, prv[2]/prv[1]) 
      theROC <- roc(pred_outcome$outcomes, pred_outcome$predictions, percent = TRUE);theROC
      bestcutoff =as.double(coords(theROC, "best", best.method =  "closest.topleft",
                                   best.weights=best.weights,
                                   ret="threshold", transpose = FALSE));bestcutoff
    
    }
    
    
    
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
        anti_join(t1.causal_reversecausal,by=c("src"="src","trgt"="trgt")) 
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
    
    #Train Random Forest model
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
    mdlColNames=c('Euclidean','Spearman','Kendall','Pearson','efMIempirical','ewMIshrink','ewMIempirical','ewMImm')
    trainingTarget='class2'
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
    
    relativeCostfn_fp = 1/2
    
    #find best cutoff
    prv = table(pred_outcome.rlv$outcomes)
    best.weights=c(relativeCostfn_fp, prv[2]/prv[1])
    
    pred_outcome.rlv = pred_outcome.rlv %>% mutate(predictions.rlv = predictions )#Pearson
    theROC.rlv <- roc(pred_outcome.rlv$outcomes, pred_outcome.rlv$predictions.rlv, percent = TRUE);
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


#Enrichment analysis on Gene clusters 
{
  
  library(clusterProfiler)
  library(org.Mm.eg.db);
  library(org.Hs.eg.db)
  library(GO.db)
  
  
  columns(org.Mm.eg.db)
  keytypes(org.Mm.eg.db)

  xx <- as.list(org.Mm.eg.db) #org.Hs.eg.db #"org.Rn.eg.db"

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
      slice_head(n=20)
    
    
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
      slice_head(n=100) %>% mutate(pvalNegLog10 = -log10(pvalue)) %>% dplyr::arrange(pvalNegLog10)
    
    p<-ggplot(enr_go_select, aes(x=reorder(Description, pvalNegLog10), y=pvalNegLog10)) +#, fill=in_grp
      geom_bar(stat="identity")+theme_minimal() + coord_flip() +
      labs(title= 'CICT Enriching top 100 TFs with most targets')+ #'CICT Enriching all then top 20 count')+ # 
      xlab('')+ylab('Negative Log 10 of p-value')+
      background_grid() +
      theme( legend.position="none",axis.text.y = element_text(angle = 45))+
      guides(colour=guide_legend(title=""));  p
    
    thepath = paste0(url.rslts, "GO terms enriched CICT2.tif")
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
            
            enr_go_select = enr_go  %>%   as.data.frame() %>%   dplyr::filter(ID %in% desiredGOTermID)
            enr_go_select = enr_go_select %>% 
                               dplyr::mutate(predictedTF = predTF, 
                                      predictedTargets = paste0(predTF.targets,collapse = ',')) %>%
                               dplyr::select(predictedTF,everything() ,predictedTargets )
        },
                      desiredGOTermID=c('GO:0019827' ))
    
    predTFs.enriched.d = rbindlist(predTFs.enriched) %>%
      mutate(pvalNegLog10 = -log10(pvalue), in_grp = predictedTF %in% pp.symbs)  %>% 
      dplyr::arrange(pvalNegLog10)
    
    p<-ggplot(predTFs.enriched.d, aes(x=reorder(predictedTF, pvalNegLog10), y=pvalNegLog10, fill=in_grp)) +#
      geom_bar(stat="identity")+theme_minimal() + coord_flip() +
      labs(title= 'CICT Enriching top 20 TFs using their targets for GO:0019827')+ #'CICT Enriching all then top 20 count')+ # 
      xlab('')+
      scale_y_continuous(name='- Log10 (p-value)', limits=c(0, 33),breaks = seq(0,50,5))+
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
    thepath =paste0(basePath, "TFs enriched CICT.tif")
    ggsave(filename = thepath, plot = p, width = 2, height = 4, device='tiff', dpi=300)
    
    
  }   
  

  
  #1 Top 20 Random forest on relevance network Top 20 genes with most targets in 
  {
    
    n3d.v.rlv.genes = c(n3d.e.rlv$Gene1,n3d.e.rlv$Gene2) %>% unique()
    pp.genesTargetingPP.rlv =n3d.e.rlv %>%  group_by(Gene1) %>% 
      summarize(trgtCount = sum(Gene2 %in% pp.symbs)) %>%
      arrange(desc(trgtCount)) %>% 
      mutate(in_grp= Gene1 %in% pp.symbs) %>%
      slice_head(n=26) #27
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
    
    thepath = paste0(url.rslts, "GO terms enriched relev2_equalnet.tif")
    ggsave(filename = thepath, plot = p, width = 4, height = 3.5, device='tiff', dpi=300)
    
    
  } 
    
  #4  RF Top 30 predicted TFs Enrichment with their targets for the term "pluripotency"  stem cell population maintenance(GO:0019827)
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
               
               enr_go_select = enr_go  %>%   as.data.frame() %>%   dplyr::filter(ID %in% desiredGOTermID)
               enr_go_select = enr_go_select %>% 
                 mutate(predictedTF = predTF, 
                        predictedTargets = paste0(predTF.targets,collapse = ',')) %>%
                 select(predictedTF,everything() ,predictedTargets )
               enr_go_select
               
             },desiredGOTermID=c('GO:0019827','GO:0098727' ))
    
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
    

  
  
}

