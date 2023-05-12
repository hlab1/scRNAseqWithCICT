################################################################################@
# © 2016 Abbas Shojaee <abbas.shojaee@gmail.com>
# Causal Inference Using Composition of Transactions CICT
# All rights reserved, please do not use or distribute without written permission
################################################################################@

################################################################################@
# © 2021 Abbas Shojaee <abbas.shojaee@gmail.com> - Carol Huang Lab, NYU
# Causal Inference Using Composition of Transactions CICT
# All rights reserved, please do not use or distribute without written permission
################################################################################@



#Global vars
#Early threshold
MAX_MEM_SIZE= "20g"

#Libraries ----
{
  library(caret)
  library(Rtsne)
  library(umap)
  library(devtools)
  library(hutils)
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
  library(knitr)
  `%<-z%` <- zeallot::`%<-%`
  PSUDEO_ZERO = 0.000001
  
}

# Function_definitions -----

#adds different formats of data as columns of a data.frame r
addToRecord<-function(r,val,name = NA){
  if(is.na(name) & !is.data.frame(val) & length(val)==1 ) name = names(val)
  if(length(val)==1 & !is.data.frame(val) ) r= bind_cols(r,val,!!name:=name,.name_repair='minimal')
  if(length(val)>1 & !is.data.frame(val) ) r= bind_cols(r,t(val),.name_repair='minimal')
  if(is.data.frame(val)) r= bind_cols(r,val,.name_repair='minimal')
  return(r)
}

#adds different types of message to a log or csv file instantly
addToReport<-function(msg,url,append=T,display=F,type = NA, addColnames = NA,addrownames = F){
  
  if(file.exists(url) & append == F) addColnames = F else addColnames = T 
  if(is.character(msg) & (type == 'text' | is.na(type))) write(msg, file = url, append = TRUE)
  
  if((is.data.frame(msg) | is.vector(msg)) & (type == 'data.frame' | is.na(type))){
    write.table(msg, file = url, sep = ",",quote = FALSE,
                append = append, col.names = FALSE, row.names = addrownames)
  }
  
  if(display) cat(msg) #paste0(msg,collapse= ','))
}


#Remove large objects
#sort( sapply(ls(),function(x){object.size(get(x))}),decreasing = F)
#rm(thedist,e,recipe_1,t1.rnd,dt.edge.BeforeAugmenting)
#gc()

#Removes duplicated cols
removeDups<- function(dt,excluded,samplesize=NA){
  if(!is.na(samplesize)) dt = sample_n(dt,samplesize)
  library(digest)
  dupslist=lapply(dt, digest)
  hashs = data.frame(cls = names(dupslist), hash = unlist(dupslist),dups = duplicated(dupslist))
  hashs =  setDT(hashs)[hashs,, on = "hash"][cls != i.cls,]
  dups = hashs[!(cls %in% excluded | i.cls %in% excluded) & dups==T,]$i.cls
  unique(dups)
}

removeDups1<- function(dt,excluded,samplesize=NA){
  if(!is.na(samplesize)) dt = sample_n(dt,samplesize)
  library(digest)
  dupslist=lapply(dt, digest)
  hashs = data.frame(cls = names(dupslist), hash = unlist(dupslist),dups = duplicated(dupslist))
  dups = duplicated(hashs$hash);table(dups)
  
  unique(hashs$cls[dups])
}

changeColTypes<-function (d, type,theCols){
  setDT(d)
  theCols = intersect(theCols,colnames(d))
  if(type== "tonumeric") cnvrt <- function(x) as.numeric(as.character(x))
  if(type== "tofactor") cnvrt <- function(x) as.factor(x)
  d1 = d[, setdiff(names(d) ,theCols), with = FALSE]
  d = d[, lapply(.SD, cnvrt), .SDcols = theCols];
  d = cbind(d, d1)
  return(d)
}


trim <- function (x) gsub("^\\s+|\\s+$", "", x)

asjson<-function(raw,file)
{
  c1<-sapply(raw,rawToChar)
  c<-prettify(c1)
  write(c,file=file)
  #save(c,file="style.json",ascii=TRUE)
}

asobject<-function(raw,thefile)
{
  c1<-sapply(raw,rawToChar)
  c<-prettify(c1)
  write(c,file=thefile)
  c2<-readLines(thefile)
  o<-fromJSON(as.character(c))
  return(o)
}

loadcyobject<-function(thefile)
{
  c<-readLines(thefile)
  o<-fromJSON(c)
  return(o)
}

gettable<-function(url,filename)
{
  tmp<-GET(url)
  tmp1<-asobject(tmp[6],filename)
  return(tmp1$rows)
}


abreviateDesc <- function(description) {
  
  output=gsub("Congestive heart failure","CHF",description, ignore.case=TRUE)
  output=gsub("Acute myocardial infarction","AMI",output, ignore.case=TRUE)
  output=gsub("Chronic obstructive pulmonary disease","COPD",output, ignore.case=TRUE)
  #  output=gsub("[Oo]ther","Othr",output, ignore.case=TRUE)
  output=gsub("[Dd]iagnostic","Dx",output, ignore.case=TRUE)
  output=gsub("[Pp]rocedure","Prc",output, ignore.case=TRUE)
  output=gsub("[Ss]ervice(s)?","Srv",output, ignore.case=TRUE)
  output=gsub("[Ee]lectro","E.",output, ignore.case=TRUE)
  output=gsub("cardiogram","C.gram",output, ignore.case=TRUE)
  output=gsub("[Ii]nterview","Intv",output, ignore.case=TRUE)
  output=gsub("[Ee]valuation","Evl",output, ignore.case=TRUE)
  output=gsub("[Cc]onsultation","Cnslt",output, ignore.case=TRUE)
  output=gsub("[Ll]aboratory","Lab",output, ignore.case=TRUE)
  output=gsub("[Cc]hemistry","Chem",output, ignore.case=TRUE)
  output=gsub("[Hh]ematology","Hem",output, ignore.case=TRUE)
  output=gsub("[Mm]icroscopic","Micr",output, ignore.case=TRUE)
  output=gsub("[Tt]herapeutic","Trap",output, ignore.case=TRUE)
  output=gsub("[Tt]ransfusion","trans.",output, ignore.case=TRUE)
  output=gsub("[Cc]athetrization","cath.",output, ignore.case=TRUE)
  return(output)}

library(arules)
myDiscretize<- function(v,cats,justcuts=TRUE,mthd="interval"){
  v1= na.omit(v)
  v1=as.numeric(v1)
  v2=v1[ is.finite(v1)] #complete.cases(v1) &&
  v3= discretize(v2,method=mthd,categories=cats, onlycuts=justcuts)
  round(v3,3)
}

myNormalize <- function(x) {
  if(is.numeric(x)){
    (x - min(x, na.rm=TRUE))/(max(x,na.rm=TRUE) - min(x, na.rm=TRUE))
  } else x}

# 'method = "mode" [default]: calculates the mode for a unimodal vector, else returns an NA
#     method = "nmodes": calculates the number of modes in the vector
#     method = "modes": lists all the modes for a unimodal or polymodal vector'

modeav <- function(x, method = "mode", na.rm = FALSE){
  x <- unlist(x)
  if (na.rm)
    x <- x[!is.na(x)]
  u <- unique(x)
  n <- length(u)
  #get frequencies of each of the unique values in the vector
  frequencies <- rep(0, n)
  for (i in seq_len(n)) {
    if (is.na(u[i])) {
      frequencies[i] <- sum(is.na(x))
    }
    else {
      frequencies[i] <- sum(x == u[i], na.rm = TRUE)
    }
  }
  #mode if a unimodal vector, else NA
  if (method == "mode" | is.na(method) | method == "")
  {return(ifelse(length(frequencies[frequencies==max(frequencies)])>1,NA,u[which.max(frequencies)]))}
  #number of modes
  if(method == "nmode" | method == "nmodes")
  {return(length(frequencies[frequencies==max(frequencies)]))}
  #list of all modes
  if (method == "modes" | method == "modevalues")
  {return(u[which(frequencies==max(frequencies), arr.ind = FALSE, useNames = FALSE)])}
  #error trap the method
  warning("Warning: method not recognised.  Valid methods are 'mode' [default], 'nmodes' and 'modes'")
  return()
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

mySkewness <- function (x,default = 0, na.rm = TRUE,...)
{
  result = skewness(x,na.rm,...)
  if(is.na(result)) default else result
}

myKurtosis <- function (x,default =3, na.rm = TRUE,...)
{
  result = kurtosis(x,na.rm,...)
  if(is.na(result)) default else result
}

myMedianAbsoluteDeviation <- function (x, center = median(x), constant = 1.4826, na.rm = FALSE,
                                       low = FALSE, high = FALSE,default =0) {
  result = mad(x, center , constant , na.rm,low , high )
  if(is.na(result)) default else result
}

myMin=function(x,lowwerLimit=-Inf) {
  minim= Inf;
  t=lapply(x, function(x) if(is.finite(x) & x < minim) x else minim )
  mn = min(unlist(t),na.rm=TRUE)
  if(lowwerLimit==-Inf) mn else max(mn,lowwerLimit)
}

myMax=function(x,upperLimit=Inf) {
  maxim= -Inf;
  t=lapply(x, function(x) if(is.finite(x) & x > maxim) x else maxim )
  mx=max(unlist(t),na.rm=TRUE)
  if(upperLimit==Inf) mx else min(mx,upperLimit)
}

my_replace_val<-function(df,valToRpl,rplist)
{
  rplptrns=names(rplist)
  rplvals=unlist(rplist)
  
  dfcols=colnames(df)
  rplcols=lapply(rplptrns,function(x) dfcols[str_detect(dfcols,paste0(".*",x))])
  
  rplist2=list()
  for(i in 1:length(rplist)){
    l=unlist(rplcols[i]);v=rplvals[i]
    l1=rep(v,length(l))
    names(l1)<-l
    rplist2=append(rplist2,l1)
  }
  
  replaceFunc = function(i,val,env){
    y=rplist2[i]
    replacement=rplist2[[i]]
    colname=names(y)
    #colName=name
    
    env$df[,colname]=ifelse(env$df[,colname]==val,replacement,env$df[,colname])
  }
  setDF(df)
  lapply(1:length(rplist2), replaceFunc,valToRpl,environment())
  #df = do.call(data.frame,lapply(1:length(rplist2), replaceFunc,val))
  df
}

myAbbreviate<-function(x,ln, collapse="", maxLength=25 ){
  xs=str_split(x,"[, ]")
  xf=lapply(xs,function(strWords) toupper(str_extract_all(strWords,"^.")))
  xr=lapply(xs,function(strWords){
    strWords=unlist(strWords)
    strWords1= ifelse(str_length(strWords)>ln, str_replace(strWords, "([aieou???]{2})",""),strWords)
    str_extract_all(strWords1,paste0("(?<=^.).{1,",ln,"}"))
  })
  xC=lapply(1:length(xs),function(i){ stf=unlist(xf[i]);st=unlist(xr[i]);paste0(stf,st)})
  f=unlist(lapply(xC,function(strAbs) str_sub( paste(strAbs,collapse=collapse),1,maxLength)))
}

myAbbreviate1<-function(x,ln, collapse="", maxLength=25 ){
  x=str_to_title(x)
  xs=str_split(x,"[, ]")
  xC=lapply(xs,function(strWords){ unlist(str_extract_all(strWords,paste0("^.{1,",ln,"}")))})
  f=unlist(lapply(xC,function(strAbs) str_sub( paste(strAbs,collapse=collapse),1,maxLength)))
}


#defults
extractLmoments = function(v)
{
    #browseURL("https://en.wikipedia.org/wiki/L-moment")
  require(Lmoments)
  if(length(v)==1){
    # if just one value reutrns defualst for normal distirbution
    #Normal 	L1=mean, L2=	σ /√π ,tau3=l_skewness=	0 ,tau4=L_kurtosis=	0.1226
    data.frame(L1=v,L2=0,L3=NA,L4=NA,tau3=0,tau4=0.1226)
  } else
  {
    l =Lmoments(v,returnobject = TRUE)
    if(is.null(l$ratios)) l$ratios=c(NA,NA,0,0.1226)  #happens when length(v)==2
    p=data.frame(L1=l$lambdas[1],L2=l$lambdas[2],L3=l$lambdas[3],L4=l$lambdas[4],
                 tau3=l$ratios[3],tau4=l$ratios[4])
    p
  }
}

.myqsnorm<- function(p, xi)
{
  m1 = 2/sqrt(2 * pi)
  mu = m1 * (xi - 1/xi)
  sigma = sqrt((1 - m1^2) * (xi^2 + 1/xi^2) + 2 * m1^2 - 1)
  g = 2/(xi + 1/xi)
  sig = sign(p - 1/2)
  Xi = xi^sig
  p = (Heaviside(p - 1/2) - sig * p)/(g * Xi)
  Quantile = (-sig * qnorm(p = p, sd = Xi) - mu)/sigma
  Quantile
}

myQSnorm = function (p, mean = 0, sd = 1, xi = 1.5)
{
  p=na.omit(p)
  sd=ifelse(is.na(sd),1,sd)
  xi=ifelse(is.na(xi),1.5,xi)
  mean=ifelse(is.na(xi),1.5,mean)
  
  result = .myqsnorm(p = p, xi = xi) * sd + mean
}

read.tsv <- function(..., header=F, sep='\t', quote='', comment='', na.strings='', stringsAsFactors=FALSE) {
  args = list(...)
  args$header = header
  if (!is.null(args$col.names)) {
    # read.delim() is not smart about this.  Yikes.
    args$header = FALSE
  }
  args$sep = sep
  args$quote = quote
  args$comment = comment
  args$stringsAsFactors = stringsAsFactors
  args$na.strings = na.strings
  do.call(read.delim, args)
}

write.tsv <- function(..., header=NA, col.names=F, row.names=F, sep='\t', na='', quote=F) {
  # 'header' to 'col.names' naming consistency with read.table()
  if (is.finite(header)) col.names = header
  write.table(..., col.names=col.names, row.names=row.names, sep=sep, na=na, quote=quote)
}

# Load data extract ground truth -----
print('stepA0')
if(exists('url.logfile'))
  if(url.logfile!="noLog"){
    write(paste0('===== ====== CICT analysis, Abbas Shojaee, Huang Lab, Now: ', lubridate::now(), '======= '),  file = url.logfile, append = TRUE)
    write('Step 0',  file = url.logfile, append = TRUE)
  }


#' prepareEdgeFeatures
#'
#' Constructs CICT features
#'
#' @param Debug If TRUE, function enters debugging mode in critical stops along the exectuion allowing verification of variables and features
#' @param cictRawEdgeCol Name of the column of the desired gene-gene association measure 
#' @param dt.geneexp Gene expression data, expects a data.frame with genes in the rows and evaluations or cells in the columns 
#' @param earlyThresholdForGraphAnalysis A threshold on rawEdgeCol that removes values below threshold before furhter processing. Can be used to reduce memory requirement 
#' @param url.logfile The path to desired location for the log file  
#' @param PSUDEO_ZERO Defualt is 0.000001. This value replace 0 in denominator of various calculations.
#' @param globalDistDiscretization Default is TRUE. Uses the global range of values when discretizing gene Conf. and Contrib. If FALSE uses gene context 
#' @return Returns a list consisted of three objects 
#' # rcrd: is a list object of intermediary objects
#' # edges: a dataframe of edge objects and CICT features for edges
#' # Vertices: a dataframe of vertices objects and CICT features for vertices
#' @examples
#' # Example usage of the function
#' c(rcrd,edges,vertices) %<-z% prepareEdgeFeatures(Debug=Debug)
#' @export
#' 
#' 

prepareEdgeFeatures<-function(Debug = F,
                              cictRawEdgeCol ='Weight',
                              dt.edge, # changed n.itm.e -> dt.edge
                              dt.geneexp, # changed all.dt -> dt.geneexp
                              earlyThresholdForGraphAnalysis,
                              url.logfile,
                              PSUDEO_ZERO = 0.000001){

  dt.edge.back00 = dt.edge
  dt.edge$edgeTyp = '' #'ewMIempirical' #tlag1.mi  'ewMImm' #"tlag1.mi" #ewMIempirical #Pearson #'corP'
  dt.edge=dt.edge %>% dplyr::mutate(Weight=.data[[cictRawEdgeCol]]  ) %>% 
    select(src,trgt,everything()) 
  
  
  dt.vertices = data.table()  #change n.itm.v -> dt.vertices
  dt.vertices=setDF(dt.geneexp) %>%  
    dplyr::mutate(ocr = rowSums(.[which(!colnames(setDF(dt.geneexp)) %in% c('vID','gene'))])) %>% 
    dplyr::mutate(ifelse(ocr ==0 | is.na(ocr), PSUDEO_ZERO,ocr)) %>%
    dplyr::select(gene,ocr) %>% rename(vID=gene) 
  
  
  dt.vertices=setDT(dt.vertices)[,`:=`(OcrInp=ocr)][,ocr:=NULL]
  
  msg = paste0(sprintf("Raw edges= %s, vertices= %s",  nrow(dt.edge),nrow(dt.vertices)),
               paste0(" Early pruning threshold= ",earlyThresholdForGraphAnalysis),
               sprintf(" Remained edges= %s \n",  nrow(dt.edge[dt.edge$Weight>earlyThresholdForGraphAnalysis,])))
  cat(msg)
  rcrd=c(rcrd, c(rawEdges=nrow(dt.edge),
                 vertices=nrow(dt.vertices),
                 earlyThreshold =earlyThresholdForGraphAnalysis))
  
  if(url.logfile!="noLog")   write(msg, file = url.logfile, append = TRUE)
  
  
  
  rm(ig)
  ig=graph.data.frame(dt.edge[dt.edge$Weight>earlyThresholdForGraphAnalysis,],directed=FALSE) 

  ind = igraph::degree(ig,v=V(ig),mode = c("in")); ind=data.frame(Indegree=ind,subcat=names(ind))

  
  dt.vertices=dt.vertices %>% left_join(ind, by=c( "vID" ="subcat"))
  dt.vertices=dt.vertices %>% left_join(outd, by=c( "vID" ="subcat"))
  
  dt.vertices = dt.vertices %>% mutate(Indegree =ifelse(is.na(Indegree),PSUDEO_ZERO,Indegree),
                               Outdegree =ifelse(is.na(Outdegree),PSUDEO_ZERO,Outdegree))
  
  
  vertexOcrSumIn = sum(dt.vertices$OcrInp)
  # vertexOcrSumOut = sum(dt.vertices$OcrOut)
  setDT(dt.vertices);dt.vertices[,c('probInp'):=list(OcrInp/vertexOcrSumIn) ] #probability
  setDF(dt.vertices)  ;setDF(dt.edge)
  vinfcols=c("vID","OcrInp","probInp","Indegree","Outdegree")
  
  dt.edge = dt.edge %>% left_join(dt.vertices[,vinfcols],by=c("src"="vID"))
  dt.edge = dt.edge %>% left_join(dt.vertices[,vinfcols],by=c("trgt"="vID"))
  

  dt.vertices.srcsum = dt.vertices%>% left_join(dt.edge,by=c("vID"="trgt")) %>%
    group_by(vID) %>% summarise(SumOcrInp=sum(OcrInp.x),na.rm=TRUE)
  
  dt.vertices = dt.vertices %>% left_join(dt.vertices.trgtsum,by=c("vID"="vID"),copy=TRUE) %>%
    inner_join(dt.vertices.srcsum,by=c("vID"="vID"),copy=TRUE) #%>%
  
  cols= c("SumOcrInp","OcrInp") 
  factorToNumeric <- function(x) as.numeric(as.character(x))
  
  dt.vertices.num = setDT(dt.vertices)[, lapply(.SD, factorToNumeric), .SDcols = cols];
  dt.vertices = dt.vertices[, !names(dt.vertices) %in% cols, with = FALSE];gc()
  dt.vertices = cbind(dt.vertices.num, dt.vertices);rm(dt.vertices.num);gc()
  
  setDT(dt.vertices);
  dt.vertices[is.na(SumOcrInp),SumOcrInp:=1];

  if(Debug ) browser()
  vinfcols=c("vID","SumOcrInp")
  setDF(dt.edge);setDF(dt.vertices)
  dt.edge = dt.edge %>% inner_join(dt.vertices[,vinfcols],by=c("src"="vID"))
  dt.edge = dt.edge %>% inner_join(dt.vertices[,vinfcols],by=c("trgt"="vID"))
  
  #a=sapply(dt.edge, function(x) sum(is.na(x)));a[a>0]
  
  dt.edge = replace_na(dt.edge,replace=list(SumOcrInp.x=0))
  
  emissionRate = 1
  maxItmstOcr = max(c(dt.vertices$OcrInp))
  
  
  dt.edge$OcrInp.y = as.numeric(dt.edge$OcrInp.x )
  dt.edge$OcrInp.y = as.numeric(dt.edge$OcrInp.y )
  
  dt.edge = setDF(dt.edge) %>% dplyr::mutate(conf=Weight/ OcrInp.x , contrib=Weight/OcrInp.y) 
  dt.edge = dt.edge %>% mutate(conf = ifelse(is.na(conf) | is.infinite(conf),0,conf),
                               contrib=ifelse(is.na(contrib) | is.infinite(contrib),0,contrib))
  
  dt.edge.back0 = dt.edge;dt.vertices.back0 = dt.vertices #dt.vertices=dt.vertices.back0;dt.edge=dt.edge.back0
  rm(dt.vertices.srcsum,dt.vertices.trgtsum,dt.vertices.influx,dt.vertices.outflux,dt.edge.tmp)
  
  # Add CICT zones -----
  {
    
    
    ngens = nrow(dt.geneexp)
    Nobservations=ncol(dt.geneexp)-1
    nbins=as.integer(sqrt(Nobservations)*1)
    breaks =1:nbins;breaks = NULL
    #breaks = seq(0,maxEdgeWeight,maxEdgeWeight/nbins+0.0001)
    
    frchtabscisse = 1:nbins
    maxEdgeWeight = max(dt.edge$conf,dt.edge$contrib)
    
    myhist = function(v,nbins,breaks=NULL,plot=F,prob = T)
    {
      h=""
      if(!is.null(breaks)) h = hist( v,breaks = breaks,plot=F,probability = prob)  else h = hist( v,nbins,plot=F,probability = prob)
      h=c(h,rep(0,nbins-length(h)))
      if(prob ==T) h$density else h$counts
    }
  
    
    
    dt.edge=dt.edge  %>% dplyr::mutate(confdisc = NULL, contribdisc=NULL)
    dt.edge = dt.edge  %>% dplyr::mutate(confdisc = unlist(infotheo::discretize(conf, "equalfreq", nbins)),
                                         contribdisc = unlist(infotheo::discretize(contrib, "equalfreq", nbins))) 
    #a=infotheo::discretize(na.omit(dt.edge$conf), "equalfreq", 10)$X
    e.cnfcnt= dt.edge%>% dplyr::select(src,trgt, confdisc,contribdisc) #conf, contrib,
    
    setDT(dt.vertices);dt.vertices[,`:=`(confhist=NULL,contribhist=NULL)]
    setDT(e.cnfcnt)
    setkeyv(e.cnfcnt, c('src','trgt'))
    
    
    
  }
  
  #1- Add distribution of conf and contribs of nodes edges to each node 
  thisrow = 0;j=1;i=1
  
  errCnt = 0
  for(j in (1):ngens){
    tryCatch({
      e.a=e.b=e.ab=NULL
      firstg = dt.vertices[j,]$vID
      
      if(!globalDistDiscretization){ #Discretization was applied on edges of each particular node
        a.confdisc = myhist(infotheo::discretize(e.cnfcnt[ src == firstg,]$conf, "equalwidth", nbins)$X,nbins,breaks)
        a.contribdisc =  myhist(infotheo::discretize(e.cnfcnt[ src == firstg,]$contrib, "equalwidth", nbins)$X,nbins,breaks)
        
      }else if(globalDistDiscretization){ #Discretization was applied on all edges
        suppressMessages({
          a.confhist = myhist( e.cnfcnt[ src == firstg,]$confdisc,nbins,breaks,plot=F,prob = F)
          a.contribhist =  myhist(e.cnfcnt[ src == firstg,]$contribdisc,nbins,breaks,plot=F,prob = F)
        })
      }
      dt.vertices[vID==firstg, `:=`(confhist = list(a.confhist))]
      dt.vertices[vID==firstg, `:=`(contribhist = list(a.contribhist))]  
    }, error = function(e) {
      errCnt<<-errCnt+1
      # print(sprintf("Error in row: %s, first: %s, " ,j,firstg))
    }
    )
    
  }
  print(sprintf(" %s nodes didn't have edges for confhist or contribhist calculations" ,errCnt)) 

  rm(dt.edge.back3,dt.edge.BeforeAugmenting,dt.vertices.BeforeMoments,dt.edge.noselfedge)
  
  #   Add Harmonized transition rates HTR-----
  
  dt.edge = dt.edge %>%dplyr::mutate(srctrgtSum = OcrInp.x+OcrInp.y,
                                     srctrgtProduct=OcrInp.x * OcrInp.y) %>%
    mutate(srctrgtSum=ifelse(is.na(srctrgtSum)| srctrgtSum==0,PSUDEO_ZERO,srctrgtSum),
           srctrgtProduct=ifelse(is.na(srctrgtProduct)|srctrgtProduct==0 ,PSUDEO_ZERO,srctrgtProduct) ) %>%
    mutate(HTR=(Weight*srctrgtSum)/srctrgtProduct,
           EE=OcrInp.x^2*OcrInp.y/srctrgtSum)#Harmonized transition rate
  
  b = dt.edge %>% filter(is.na(HTR))
  
  self=dt.edge %>% group_by(src) %>% summarise( scfSumHTR=sum(HTR,na.rm= TRUE),scfSumEE=sum(EE,na.rm=TRUE))
  others=dt.edge %>% group_by(trgt) %>%  summarise( ocbSumHTR=sum(HTR,na.rm= TRUE),ocbSumEE=sum(EE,na.rm=TRUE))
  
  dt.vertices.back=dt.vertices
  dt.vertices = dt.vertices %>% left_join(self,by=c("vID"="src"))
  dt.vertices = dt.vertices %>% left_join(others,by=c("vID"="trgt"))
  
  a=sapply(dt.vertices, function(x) sum(is.na(x)));a[a>0]
  #View(dt.vertices[,c("vID","scfSumHTR","ocbSumHTR"),with=FALSE])
  setDF(dt.vertices);dt.vertices = tidyr::replace_na(dt.vertices,replace=list(scfSumHTR=0,ocbSumHTR=0,
                                                                  scfSumEE=0,ocbSumEE=0,DiagnosesAb="",SumOcrInp=0))
  
  dt.vertices1=dt.vertices[,c("vID","scfSumHTR","ocbSumHTR","scfSumEE","ocbSumEE")]
  
  setDF(dt.edge); setDF(dt.vertices1)
  dt.edge=dt.edge %>% left_join(dt.vertices1,by=c("src"="vID")) 
  dt.edge=dt.edge %>% left_join(dt.vertices1,by=c("trgt"="vID"))
  
  
  #   Enhance edges conf contrib -----
  if(!exists("PSUDEO_ZERO")) stop("define PSUDEO_ZERO")
  dt.edge = setDF(dt.edge) %>% #dplyr::mutate(Weight=mf.mi) %>%
    dplyr::mutate(V=abs(OcrInp.x-OcrInp.y), R=V/Weight,P=V*Weight,
                  
                  NHTRfromSource=HTR/scfSumHTR.x,
                  NHTRtoTarget=HTR/ocbSumHTR.y,
                  tNHTR=NHTRfromSource+NHTRtoTarget,
                  
                  NEEfromSource=EE/scfSumEE.x,
                  NEEtoTarget = EE/ocbSumEE.y,
                  
                  
                  NEE_NHTR_Ratio=NHTRfromSource/ifelse(NEEfromSource==0,PSUDEO_ZERO,NEEfromSource),
                  
                  confN = conf*(OcrInp.y/srctrgtSum), #NHTRfromSource
                  contribN = contrib*(OcrInp.x/srctrgtSum), #NHTRtoTarget, #
                  tn=confN+contribN,t=conf+contrib,
                  #+1 is necessary to protect recursive edges from going infinity
                  
                  Pout=(OcrInp.x-OcrInp.y)*Weight/OcrInp.x,
                  Pin=(OcrInp.x-OcrInp.y)*Weight/OcrInp.y,#/maxItmstOcr
                  slope = abs(OcrInp.x-OcrInp.y)/Weight, #tng
                  hypotenuse = sqrt(V^2+Weight^2),
                  EO =emissionRate *(srctrgtProduct)/(SumOcrInp.x+1),
                  EI =emissionRate *(srctrgtProduct)/(SumOcrInp.y+1),
                  OEER=Weight/EO,OERR=Weight/EI,
                  toe=OEER + OERR,
                  AM=(OEER+OERR)/2,GM=(OEER*OERR)^.5, HM=2/((1/OEER)+(1/OERR)),
                  UR=V/(Weight-EO) ,URR= V/(Weight-EI) ) #unexpected resistance
  
  a=sapply(dt.edge, function(x) sum(is.na(x)));a[a>0]
  b=sapply(dt.edge, function(x) sum(x==0 & !is.na(x)) );b[b>0]
  c=sapply(dt.edge, function(x) sum(is.infinite(x)));c[c>0]
  #d=setDT(dt.edge)[is.infinite(contribZ),.(tz,odds,contribZ,OcrInp.x,OcrInp.y,I,IBA,DiagnosesAb.x,DiagnosesAb.y)];View(d)
  #setDT(dt.edge);rpt=dt.edge[Weight==0,.(src,trgt,conf,contrib,Weight)];View(rpt)
  #Outdegree.x, #scaled value of expected weight
  
  if(url.logfile!="noLog")   write('Step 1',  file = url.logfile, append = TRUE)
  
  #   Add Power parameters to vertices and edges ----
  #all the power input from different sources to each trgt
  powerParamsIn=dt.edge %>% group_by(trgt) %>%
    summarise( Pinsum=sum(Pin,na.rm=TRUE),PinAbsSum = sum(abs(Pin),na.rm=TRUE),
               PinSD=sd(Pin,na.rm=TRUE),PinMean=mean(Pin,na.rm=TRUE))
  powerParamsOut=dt.edge %>% group_by(src) %>%
    summarise(Poutsum=sum(Pout,na.rm=TRUE),PoutAbsSum = sum(abs(Pout),na.rm=TRUE),
              PoutSD=sd(Pout,na.rm=TRUE),PoutMean=mean(Pout,na.rm=TRUE))
  
  
  a=sapply(powerParamsOut, function(x) sum(is.na(x)));a[a>0]
  a=sapply(powerParamsIn, function(x) sum(is.na(x)));a[a>0]
  
  dt.vertices = dt.vertices %>% left_join(powerParamsIn,by=c("vID"="trgt"))
  dt.vertices = dt.vertices %>% left_join(powerParamsOut,by=c("vID"="src"))
  
  dt.vertices = replace_na(dt.vertices,replace=list(PinSD=0,PoutSD=0))
  
  vinfcols=c("vID","Pinsum","PinAbsSum","PinSD","PinMean","Poutsum","PoutAbsSum","PoutSD","PoutMean")
  setDF(dt.edge);setDF(dt.vertices)
  dt.edge = dt.edge %>% inner_join(dt.vertices[,vinfcols],by=c("src"="vID"))
  dt.edge = dt.edge %>% inner_join(dt.vertices[,vinfcols],by=c("trgt"="vID"))
  
  dt.edge.back=dt.edge #dt.edge=dt.edge.back
  dt.vertices.back=dt.vertices 
  #dt.vertices=dt.vertices.back;dt.edge=dt.edge.back
  
  
  
  
  #remove duplicated columns 
  if(F){
    library(digest)
    dupcols=c()
    e.set1 = sample_n(dt.edge,20000)
    clnames=colnames(dt.edge)
    dupcols = clnames[duplicated(lapply(e.set1, digest))];dupcols
    dupcols = dupcols[-c(1,3,2,17)] #[-c(1,2,4,5,7)]
    dt.edge = dt.edge %>% select(! any_of(dupcols))
  }
  gc()
  
  #   EARLY attach Ground truth  ----
  if(FALSE){rm(gt1,gt1.c,gt1.rc,gt1.c.r,gt1.abscentGT,gt2.rc,gt2.rnd)
    breaks = seq(0,1,.1) #breaks = seq(0,1,10)
    gt1=tbl.goldStandard %>% dplyr::select(src,trgt)
    
    
    setDF(dt.edge)
    #Test
    tmp =dt.edge.back0 %>% inner_join(gt1,by=c("first"="src","second"="trgt"))
    gt1.abscentGT =gt1 %>% anti_join(dt.edge.back0,by=c("src"="first","trgt"="second"))
    gt1.abscentGT =gt1.abscentGT %>% anti_join(dt.edge.back0,by=c("src"="second","trgt"="first"))
    
    gt1.c = dt.edge %>% inner_join(gt1,by=c("src"="src","trgt"="trgt"))
    gt1.c$predicate = "CAUSES"
    hist(gt1.c$confdisc, 10)
    
    
    gt1.rc = dt.edge %>% inner_join(gt1.c[,c('src','trgt')],by=c("src"="trgt","trgt"="src"))
    gt1.rc$predicate = "REV_CAUSES"
    hist(gt1.rc$confdisc, 10)
    
    gt1.causal_reversecausal = rbind(gt1.c,gt1.rc)
    #tmp = rbind(gt1,gt1.c,gt1.rc)
    
    #Adding 2000 random edges
    gt1.rnd = dt.edge %>% 
      anti_join(gt1.causal_reversecausal,by=c("src"="src","trgt"="trgt")) %>%
      mutate(predicate = "IRRELEVANTA")
    
    #tmp.rdndnt = setDT(gt1.rnd)[duplicated(gt1.rnd[,.(src,trgt)]),]
    
    gt1.rnd$predicate = "IRRELEVANTA"
    gt1.rnd.smpl = gt1.rnd %>% dplyr::slice_sample(n=2000)
    
    gt2=rbind(gt1.causal_reversecausal,gt1.rnd)
    
    gt2= as.data.frame(gt2) %>% dplyr::select(edgeTyp ,src, trgt,Weight,,
                                              OcrInp.x, OcrInp.y,#intvl_avg,intvl_sd,intvl_median,
                                              #scfNMAD.y,scbMedian.x,  ocbNMAD.x,scfL1.x,
                                              EI,EO,EE,everything())
    library(stringr)
    gt2$class1=unlist(lapply(gt2$predicate,function(x) {if(is.na(x) ) 'u'
      else if(x=='CAUSES') 'c'
      else if(x=='PRECEDES') 'el'
      else if(x=='ASSOCIATED_WITH') 'a'
      else if(x=='REV_CAUSES') 'rc'
      else if(str_detect(x,'IRRELEVANTA')) 'ir'
      else  'u'} ))
    table(gt2$class1,useNA="ifany")
    
    table(gt2$predicate,gt2$edgeTyp, useNA="ifany")
    
    a=sapply(gt2, function(x) sum(is.na(x)));a[a>0]
    
    gt2 = gt2 %>% rowwise() %>%
      dplyr::mutate(SUID=row_number(),class3=ifelse(class1 %in% c('c','rc'),'c','other'),class3="0",
                    shared_name=paste0(src,"-",trgt))  #class1=""
    gt2=unique(gt2)
    table(gt2$predicate,gt2$class2, useNA="ifany")
    
    colnames(dt.edge.back0)
    setDT(gt2)
    table(gt2[ mutualInf >.4,]$class1)
    table(gt2[ abs(corP) >.2,]$class1)
    table(gt2[ abs(corS) >.2,]$class1)
    table(gt2[ abs(corK) >.2,]$class1)
    
    table(gt2[ tlag1.mi >.4,]$class1)
    table(gt2[ abs(tlag1.corP) >.2,]$class1)
    table(gt2[ abs(tlag1.corS) >.2,]$class1)
    table(gt2[ abs(tlag1.corK) >.2,]$class1)
  }
  #write.table(gt2,file="Predicates gt2 Data100-1.txt",sep = "\t")
  
  #   check histogram of new variables -----
  if(FALSE)
  {
    
    
    
    tmp.e=dt.edge %>% dplyr::select(src,trgt,OcrInp.x,OcrInp.y,Weight,Weight,V,R,EO,EI,OEER,OERR,AM,GM,HM,UR)
    
    summary(dt.edge$EE);summary(dt.edge$NEEfromSource);summary(dt.edge$NEEtoTarget)
    summary(dt.edge$UR);summary(dt.edge$EO);summary(dt.edge$EI); summary(dt.edge$OERR);summary(dt.edge$OEER)
    summary(dt.edge$AM);summary(dt.edge$HM);summary(dt.edge$GM)
    summary(dt.edge$Pin);summary(dt.edge$Pout)
    summary(dt.edge$EI)
    
    ggplot(dt.edge, aes(EE))+ggtitle("EE") + geom_histogram()  + scale_x_log10()
    ggplot(dt.edge, aes(NEEfromSource))+ggtitle("NEEfromSource") + geom_histogram()  + scale_x_log10()
    ggplot(dt.edge, aes(NHTRfromSource))+ggtitle("NHTRfromSource") + geom_histogram()  + scale_x_log10()
    ggplot(dt.edge, aes(NEE_NHTR_Diff))+ggtitle("NEE_NHTR_Diff") + geom_histogram()  + scale_x_log10()
    ggplot(dt.edge, aes(NEEtoTarget))+ggtitle("NEEtoTarget") + geom_histogram()  + scale_x_log10()
    ggplot(dt.edge, aes(NHTRtoTarget))+ggtitle("NHTRtoTarget") + geom_histogram()  + scale_x_log10()
    ggplot(dt.edge, aes(tNHTR))+ggtitle("tNHTR") + geom_histogram()  + scale_x_log10()
    
    ggplot(dt.edge, aes(cbMean.x))+ggtitle("cbMean.x - log scale") + geom_histogram() +
      scale_x_log10()+geom_vline(xintercept=.006)
    ggplot(dt.edge, aes(OEER))+ggtitle("OEER - log scale")  +   scale_x_log10()  +geom_vline(xintercept=.42)
    ggplot(dt.edge, aes(EO))+ggtitle("EO - log scale")  +   scale_x_log10()
    ggplot(dt.edge, aes(EI))+ggtitle("EI - log scale") +geom_histogram()   +   scale_x_log10()
    ggplot(dt.edge, aes(cfNMad.x))+ggtitle("cfNMad.x - log scale")  +   scale_x_log10()  +geom_vline(xintercept=.00003)
    ggplot(dt.edge, aes(UR))+ggtitle("UR") + geom_histogram() + scale_x_log10()+geom_vline(xintercept=0.0003)
    ggplot(dt.edge, aes(tn))+ggtitle("tn") + geom_histogram() + scale_x_log10()
    ggplot(dt.edge, aes(t))+ggtitle("t") + geom_histogram() + scale_x_log10()
    ggplot(dt.edge, aes(AM))+ggtitle("AM") + geom_histogram() + scale_x_log10()
    ggplot(dt.edge, aes(confN))+ggtitle("confN") + geom_histogram() + scale_x_log10()
    ggplot(dt.edge, aes(Pin))+ggtitle("Pin") + geom_histogram() + scale_x_log10()
    ggplot(dt.edge, aes(Pout))+ggtitle("Pout") + geom_histogram() + scale_x_log10()
    ggplot(dt.edge, aes(RABADiff))+ggtitle("RABADiff") + geom_histogram() + scale_x_log10()
  }
 
   rm(dt.edge1,dt.edge2,dt.vertices1)
  #View(dt.edge %>% filter(abs(Idiff)>.1)%>% arrange(desc(Idiff)) %>% select(src,trgt,Weight,IBA, everything()))
  #dt.edge = replace_na(dt.edge,replace=
  #                          list(confDiff=0,confNDiff=0,directionLogRatio=0))
  
  #   Add transparency -----
 
  dt.edge=as.data.frame(dt.edge)
  n.notinVertices= dt.edge %>% anti_join(dt.vertices,by=c("src"="vID"))
  print(paste0("!!! nodes:" , paste(unique(n.notinVertices$src),collapse=","), "  do not exist in vertices"))

  dt.edge = dt.edge %>%dplyr::mutate(trnsparency= conf/max(conf))
  #tmp.e=dt.edge  %>%  select(SUID,shared_name,Weight,src,trgt,OcrInp.x,OcrInp.y,conf,contrib)#%>% filter(src!=trgt)
  
  
  if(url.logfile!="noLog") write('Step5',  file = url.logfile, append = TRUE)
  
  #   Enhance vertices influx outflux----
  dt.edge=dt.edge[, unique(colnames(dt.edge))]
  dt.vertices=dt.vertices[, unique(colnames(dt.vertices))] #removing duplicated columns
  setDF(dt.edge)
  
  dt.vertices.back=dt.vertices
  dt.edge.tmp=
    #dt.vertices=dt.vertices.back
    dt.vertices.outflux = dt.vertices %>% inner_join(dt.edge[,c("src","Weight")],by=c("vID"="src")) %>%
    group_by(vID) %>% summarise(outflux=sum(Weight))
  
  dt.vertices.influx = dt.vertices %>% inner_join(dt.edge[,c("trgt","Weight")],by=c("vID"="trgt")) %>%
    group_by(vID) %>% summarise(influx=sum(Weight))
  
  dt.vertices= dt.vertices %>% left_join(dt.vertices.outflux,by=c("vID"="vID"))
  dt.vertices= dt.vertices %>% left_join(dt.vertices.influx,by=c("vID"="vID"))
  
  
  a=sapply(dt.vertices, function(x) sum(is.na(x)));a[a>0] # missing outflux means no target for some nodes and
  #New Added
  dt.vertices=as.data.frame(dt.vertices)
  dt.vertices = replace_na(dt.vertices,replace=list(Pinsum=0,
                                            PinAbsSum=0,PinSD=0,PinMean=0,
                                            Poutsum=0,PoutAbsSum=0,PoutSD=0,PoutMean=0,
                                            outflux=0,influx=0) )
  
  dupcols = removeDups(dt.vertices,excluded = c('Indegree','Outdegree','','probOut','probIn'))
  dt.vertices = dt.vertices%>% select(! any_of(dupcols))
  
  rm(dt.vertices.influx,dt.vertices.outflux,dt.vertices1,powerParamsIn,powerParamsOut,dt.edge2,dt.edge1,dt.edge.tmp)
  #   Enhance vertices Add moments of  contribution ----
  #higher kurtosis means more of the variance is the result of infrequent extreme deviations, as opposed to frequent modestly sized deviations.The kurtosiof any univariate normal distribution is 3. platykurtic<3 , leptoKurtic >>3
  
  setDT(dt.edge);dt.edge.noselfedge=dt.edge [src!=trgt,.(src,contrib,conf,trgt)]#removing self edges from the calculation
  setDF(dt.edge);setDF(dt.vertices)
  if(Debug) browser()
  #selfContribNLMoments=dt.edge %>% dplyr::group_by(src) %>% do(extractLmoments(.$contrib))
  a = dt.edge %>% filter(is.na(contribN))
  selfContribLMoments=dt.edge %>% dplyr::group_by(src) %>% do(extractLmoments(.$contribN))
  selfContribs=dt.edge %>% group_by(src) %>% summarise(
    Total=sum(contrib,na.rm= TRUE),
    Mean=myMean(contrib,0,na.rm=TRUE),
    Median=myMedian(contrib,0,na.rm=TRUE),
    SD=mySD(contrib,0, na.rm = TRUE), # SD of less than 2 inputs =
    Skew=mySkewness(contrib,0, na.rm = TRUE),
    Kurt=myKurtosis(contrib,3, na.rm = TRUE),
    MADconst=ifelse(is.na(sn::qsc(.75,Mean,SD,Skew)),1.488,sn::qsc(.75,Mean,SD,Skew)), #myQSnorm
    MAD= myMedianAbsoluteDeviation(contrib,Median, MADconst,na.rm=TRUE,default=0),
    
    NTotal=sum(contribN,na.rm= TRUE),
    NMean=myMean(contribN, 0,na.rm= TRUE),
    NMedian=myMedian(contribN, 0,na.rm= TRUE),
    NSD=mySD(contribN, na.rm = TRUE,0),
    NSkew=mySkewness(contribN, 0,na.rm =TRUE),
    NKurt=myKurtosis(contribN, 3,na.rm = TRUE),
    NMADconst=ifelse(is.na(sn::qsc(.75,NMean,NSD,NSkew)),1.488,sn::qsc(.75,NMean,NSD,NSkew)),
    NMAD= myMedianAbsoluteDeviation(contribN,NMedian,NMADconst,na.rm=TRUE,default=0))
  
  selfContribsFull=selfContribs %>% inner_join(selfContribLMoments,by=c("src"="src"))
  names(selfContribsFull)<-paste0("scb",names(selfContribsFull))
  
  selfConfLMoments=dt.edge %>% dplyr::group_by(src) %>% do(extractLmoments(.$confN))
  selfConf=dt.edge %>% group_by(src) %>%  summarise(
    Total=sum(conf,na.rm= TRUE),
    Mean=myMean(conf,0,na.rm=TRUE),
    Median=myMedian(conf,0,na.rm=TRUE),
    SD=mySD(conf,0, na.rm = TRUE), # SD of less than 2 inputs =
    Skew=mySkewness(conf,0, na.rm = TRUE),
    Kurt=myKurtosis(conf,3, na.rm = TRUE),
    MADconst=ifelse(is.na(sn::qsc(.75,Mean,SD,Skew)),1.488,sn::qsc(.75,Mean,SD,Skew)), #myQSnorm
    MAD= myMedianAbsoluteDeviation(conf,Median, MADconst,na.rm=TRUE,default=0),
    
    NTotal=sum(confN,na.rm= TRUE),
    NMean=myMean(confN, 0,na.rm= TRUE),
    NMedian=myMedian(confN, 0,na.rm= TRUE),
    NSD=mySD(confN, na.rm = TRUE,0),
    NSkew=mySkewness(confN, 0,na.rm =TRUE),
    NKurt=myKurtosis(confN, 3,na.rm = TRUE),
    NMADconst=ifelse(is.na(sn::qsc(.75,NMean,NSD,NSkew)),1.488,sn::qsc(.75,NMean,NSD,NSkew)),
    NMAD= myMedianAbsoluteDeviation(confN,NMedian,NMADconst,na.rm=TRUE,default=0),
    
    tnTotal=sum(tn,na.rm= TRUE),
    tnMean=myMean(tn, 0,na.rm= TRUE),
    tnMedian=myMedian(tn, 0,na.rm= TRUE),
    tnSD=mySD(tn, na.rm = TRUE,0),
    tnSkew=mySkewness(tn, 0,na.rm =TRUE),
    tnKurt=myKurtosis(tn, 3,na.rm = TRUE),
    tnMADconst=ifelse(is.na(sn::qsc(.75,tnMean,tnSD,tnSkew)),1.488,sn::qsc(.75,tnMean,tnSD,tnSkew)),
    tnMAD= myMedianAbsoluteDeviation(tn,tnMedian,tnMADconst,na.rm=TRUE,default=0)
  )
  
  selfConfFull=selfConf %>% inner_join(selfConfLMoments,by=c("src"="src"))
  names(selfConfFull)<-paste0("scf",names(selfConfFull))
  
  dt.vertices.selfparams = selfContribsFull %>% inner_join(selfConfFull,by=c("scbsrc"="scfsrc"))
  
  a=sapply(dt.vertices.selfparams, function(x) sum(is.na(x)));a[a>0]
  
  dt.edge.back1 = dt.edge;dt.vertices.back1 = dt.vertices #dt.vertices=dt.vertices.back1;dt.edge=dt.edge.back1
  
  #   Enhance vertices Add moments of Confidence and contribution ----
  T<-TRUE
  othersConfLMoments=dt.edge %>% dplyr::group_by(trgt) %>% do(extractLmoments(.$confN))
  othersConfs=dt.edge %>% group_by(trgt) %>% summarise(
    Total=sum(conf,na.rm= TRUE),
    Mean=myMean(conf,0,na.rm=TRUE),
    Median=myMedian(conf,0,na.rm=TRUE),
    SD=mySD(conf,0, na.rm = TRUE), # SD of less than 2 inputs =
    Skew=mySkewness(conf,0, na.rm = TRUE),
    Kurt=myKurtosis(conf,3, na.rm = TRUE),
    MADconst=ifelse(is.na(sn::qsc(.75,Mean,SD,Skew)), 1.488,sn::qsc(.75,Mean,SD,Skew)), #myQSnorm
    MAD= myMedianAbsoluteDeviation(conf,Median, MADconst,na.rm=TRUE,default=0),
    
    NTotal=sum(confN,na.rm= TRUE),
    NMean=myMean(confN, 0,na.rm= TRUE),
    NMedian=myMedian(confN, 0,na.rm= TRUE),
    NSD=mySD(confN, na.rm = TRUE,0),
    NSkew=mySkewness(confN, 0,na.rm =TRUE),
    NKurt=myKurtosis(confN, 3,na.rm = TRUE),
    NMADconst=ifelse(is.na(sn::qsc(.75,NMean,NSD,NSkew)),1.488,sn::qsc(.75,Mean,SD,Skew)),
    NMAD= myMedianAbsoluteDeviation(confN,NMedian,NMADconst,na.rm=TRUE,default=0))
  
  othersConfsFull=othersConfs %>% inner_join(othersConfLMoments,by=c("trgt"="trgt"))
  names(othersConfsFull)<-paste0("ocf",names(othersConfsFull))
  
  othersContribLMoments=dt.edge %>% dplyr::group_by(trgt) %>% do(extractLmoments(.$contribN))
  othersContribs=dt.edge %>% group_by(trgt) %>%  summarise(
    Total=sum(contrib,na.rm= TRUE),
    Mean=myMean(contrib,0,na.rm=TRUE),
    Median=myMedian(contrib,0,na.rm=TRUE),
    SD=mySD(contrib,0, na.rm = TRUE), # SD of less than 2 inputs =
    Skew=mySkewness(contrib,0, na.rm = TRUE),
    Kurt=myKurtosis(contrib,3, na.rm = TRUE),
    MADconst=ifelse(is.na(sn::qsc(.75,Mean,SD,Skew)),1.488,sn::qsc(.75,Mean,SD,Skew)), #myQSnorm
    MAD= myMedianAbsoluteDeviation(contrib,Median, MADconst,na.rm=TRUE,default=0),
    
    NTotal=sum(contribN,na.rm= TRUE),
    NMean=myMean(contribN, 0,na.rm= TRUE),
    NMedian=myMedian(contribN, 0,na.rm= TRUE),
    NSD=mySD(contribN, na.rm = TRUE,0),
    NSkew=mySkewness(contribN, 0,na.rm =TRUE),
    NKurt=myKurtosis(contribN, 3,na.rm = TRUE),
    NMADconst=ifelse(is.na(sn::qsc(.75,NMean,NSD,NSkew)),1.488,sn::qsc(.75,NMean,NSD,NSkew)),
    NMAD= myMedianAbsoluteDeviation(contribN,NMedian,NMADconst,na.rm=TRUE,default=0)
  )
  
  othersContribsFull=othersContribs %>% inner_join(othersContribLMoments,by=c("trgt"="trgt"))
  names(othersContribsFull)<-paste0("ocb",names(othersContribsFull))
  
  dt.vertices.othersparams = othersConfsFull %>% inner_join(othersContribsFull,by=c("ocftrgt"="ocbtrgt"))
  
  a=sapply(dt.vertices.othersparams, function(x) sum(is.na(x)));a[a>0]
  
  
  
  # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ -----
  
  
  #   Integrating moments with dt.vertices -----
  a=sapply(dt.vertices, function(x) sum(is.na(x)));a[a>0]
  dt.vertices.BeforeMoments=dt.vertices
  dt.vertices = dt.vertices %>% left_join(dt.vertices.othersparams,by=c("vID"="ocftrgt"))
  dt.vertices = dt.vertices %>% left_join(dt.vertices.selfparams,by=c("vID"="scbsrc"))
  
  a=sapply(dt.vertices, function(x) sum(is.na(x)));a[a>0]
  
  my_replace_na<-function(df,rplist)
  {
    rplptrns=names(rplist)
    rplvals=unlist(rplist)
    
    dfcols=colnames(df)
    rplcols=lapply(rplptrns,function(x) dfcols[str_detect(dfcols,paste0(".*",x))]) #pattern checking
    
    rplist2=list()
    for(i in 1:length(rplist)){
      l=unlist(rplcols[i]);v=rplvals[i]
      l1=rep(v,length(l))
      names(l1)<-l
      rplist2=append(rplist2,l1)
    }
    
    df=replace_na(df,replace=rplist2)
    df
  }
  
  dt.vertices=my_replace_na(dt.vertices,rplist = list(Mean=0,Median=0,Skew=0,Kurt=3,
                                              SD=0,Total=0,NSD=0,flux=0,
                                              sum.x=0,sum.y=0,Pinsum=0,Poutsum=0,
                                              MADcost=1.488,L1=0,L2=0,tau3=0,tau4=0 )) #tau4=0.1226
  
  
  
  
  rm(dt.vertices.othersparams,dt.vertices.selfparams,othersConfs,othersConfsFull,othersContribs,othersContribsFull,othersConfLMoments,othersContribLMoments)
  rm(selfConf,selfConfFull,selfContribs,selfContribsFull,selfConfLMoments,selfContribLMoments)
  
  #   Add vertice level 2 distribution parameters ----  
  #TODO probably wrong should use edges not vertices to find those distribution params
  setDF(dt.vertices)
  L2 = dt.vertices %>% ungroup() %>% summarise(
    MeanOcfSkew=myMean(ocfSkew,0,na.rm=TRUE),
    MedianOcfSkew=myMedian(ocfSkew,0,na.rm=TRUE),
    SDOcfSkew=mySD(ocfSkew,0, na.rm = TRUE), # SD of less than 2 inputs =
    SkewOcfSkew=mySkewness(ocfSkew,0, na.rm = TRUE),
    KurtOcfSkew=myKurtosis(ocfSkew,3, na.rm = TRUE),
    MADconstOcfSkew=sn::qsc(.75,MeanOcfSkew,SDOcfSkew,SkewOcfSkew), #myQSnorm
    MADOcfSkew= myMedianAbsoluteDeviation(ocfSkew,MedianOcfSkew,
                                          MADconstOcfSkew,na.rm=TRUE,default=0),
    
    MeanScbNSkew=myMean(scbNSkew,0,na.rm=TRUE),
    MedianScbNSkew=myMedian(scbNSkew,0,na.rm=TRUE),
    SDScbNSkew=mySD(scbNSkew,0, na.rm = TRUE), # SD of less than 2 inputs =
    SkewScbNSkew=mySkewness(scbNSkew,0, na.rm = TRUE),
    KurtScbNSkew=myKurtosis(scbNSkew,3, na.rm = TRUE),
    MADconstScbNSkew=sn::qsc(.75,MeanScbNSkew,SDScbNSkew,SkewScbNSkew), #myQSnorm
    MADScbNSkew= myMedianAbsoluteDeviation(scbNSkew,MedianScbNSkew,
                                           MADconstScbNSkew,na.rm=TRUE,default=0),
    
    MeanOcfNMedian=myMean(ocfNMedian,0,na.rm=TRUE),
    MedianOcfNMedian=myMedian(ocfNMedian,0,na.rm=TRUE),
    SDOcfNMedian=mySD(ocfNMedian,0, na.rm = TRUE), # SD of less than 2 inputs =
    SkewOcfNMedian=mySkewness(ocfNMedian,0, na.rm = TRUE),
    KurtOcfNMedian=myKurtosis(ocfNMedian,3, na.rm = TRUE),
    MADconstOcfNMedian=sn::qsc(.75,MeanOcfNMedian,SDOcfNMedian,SkewOcfNMedian), #myQSnorm
    MADOcfNMedian= myMedianAbsoluteDeviation(ocfNMedian,MedianOcfNMedian,
                                             MADconstOcfNMedian,na.rm=TRUE,default=0),
    
    MeanPoutsum=myMean(Poutsum,0,na.rm=TRUE),
    MedianPoutsum=myMedian(Poutsum,0,na.rm=TRUE),
    SDPoutsum=mySD(Poutsum,0, na.rm = TRUE), # SD of less than 2 inputs =
    SkewPoutsum=mySkewness(Poutsum,0, na.rm = TRUE),
    KurtPoutsum=myKurtosis(Poutsum,3, na.rm = TRUE),
    MADconstPoutsum=sn::qsc(.75,MeanPoutsum,SDPoutsum,SkewPoutsum), #myQSnorm
    MADPoutsum= myMedianAbsoluteDeviation(Poutsum,MedianPoutsum,
                                          MADconstPoutsum,na.rm=TRUE,default=0),
    
    
    MeanScbL3=myMean(scbL3,0,na.rm=TRUE),
    MedianScbL3=myMedian(scbL3,0,na.rm=TRUE),
    SDScbL3=mySD(scbL3,0, na.rm = TRUE), # SD of less than 2 inputs =
    SkewScbL3=mySkewness(scbL3,0, na.rm = TRUE),
    KurtScbL3=myKurtosis(scbL3,3, na.rm = TRUE),
    MADconstScbL3=sn::qsc(.75,MeanScbL3,SDScbL3,SkewScbL3), #myQSnorm
    MADScbL3= myMedianAbsoluteDeviation(scbL3,MedianScbL3,
                                        MADconstScbL3,na.rm=TRUE,default=0))
  
  a=sapply(L2, function(x) sum(x==0));a[a>0]
  #TODO: add more L2 calculations to dt.vertices
  dt.vertices = dt.vertices %>%dplyr::mutate(
    ZOcfSkew= (ocfSkew-L2$MedianOcfSkew)/L2$MADOcfSkew,
    ZOcfNMedian= (Poutsum-L2$MedianOcfNMedian)/L2$MADOcfNMedian,
    ZScbNSkew= (Poutsum-L2$MedianScbNSkew)/L2$MADScbNSkew,
    ZScbScbL3= (scbL3-L2$MedianScbL3)/L2$MADScbL3,
    ZPoutsum= (Poutsum-L2$MedianPoutsum)/L2$MADPoutsum )
  
  a=sapply(dt.vertices, function(x) sum(is.na(x)));a[a>0]
  #Todo verify that replacement is correct
  # Mean and Median should be replaced with  !?
  dt.vertices=my_replace_na(dt.vertices,rplist = list(Mean=0,Median=0,Skew=0,Kurt=3,MAD=0,
                                              SD=0,Total=0,NSD=0,flux=0,sum.x=0,sum.y=0,
                                              MADcost=1.488,L1=0,L2=0,L3=0,L4=0,tau3=0,tau4=0 ))
  
  #   Add source and target charectristics to edges -----
  dt.edge.BeforeAugmenting=dt.edge
  #dt.vertices = dt.vertices %>% dplyr::mutate(DiagnosesAb.y=NULL) %>% rename(DiagnosesAb = DiagnosesAb.x)
  collist= c('MAD','Median','MADconst','Mean','SD','Skew','Kurt','Total',
             'NMAD','NMedian','NMADconst','NMean','NSD','NSkew','NKurt','NTotal',
             'tnMedian','tnMADconst','tnMAD',
             'L1','L2','L3','L4','tau3','tau4')
  collist=c(paste0("scb",collist),paste0("scf",collist),paste0("ocb",collist),paste0("ocf",collist))
  colList=c('vID','DiagnosesAb','evenTyp','Pinsum','Poutsum',
            'ZOcfSkew','ZPoutsum','ZOcfNMedian','ZScbNSkew','ZScbScbL3',collist)
  
  #print(colList[!colList %in% colnames(dt.vertices)])
  finalCoolList = colList[colList %in% colnames(dt.vertices)]
  dt.vertices1=dt.vertices[,finalCoolList]
  
  
  a=sapply(dt.vertices1, function(x) sum(is.na(x)));a[a>0]
  
  setDF(dt.edge); setDF(dt.vertices1)
  if(url.logfile!="noLog") write('Step5.1',  file = url.logfile, append = TRUE)
  
  #!!!NOTE that parameters with .x point to source parameters and those with .y point to target parameters
  if(FALSE) dt.edge = dt.edge %>% dplyr::filter(abs(Weight)>.2);gc()
  dt.edge.back2 = dt.edge
  dt.edge=as.data.frame(dt.edge)
  dt.edge=dt.edge %>% left_join(dt.vertices1,by=c("src"="vID")) 
  dt.edge = dt.edge %>%   left_join(dt.vertices1,by=c("trgt"="vID"))
  
  dupcols = removeDups1(dt.edge,excluded = c('OcrInp.x','Indegree.x'), samplesize = nrow(dt.edge) * 0.02)
  dt.edge = dt.edge%>% select(! any_of(dupcols), any_of('Weight'))
  
  a=sapply(dt.edge, function(x) sum(is.na(x)));a[a>0]
  #a=sapply(dt.vertices, function(x) sum(x==0));a[a>0]
  
  dt.edge=my_replace_na(dt.edge,rplist =
                          list(L3.x=0,L3.y=0,L4.x=0,L4.y=0 ))
  
  rm(dt.vertices.influx,dt.vertices.outflux,dt.vertices1,dt.edge.tmp,dt.edge.noselfedge)
  #   Adding confz and contribz ----
  if(T){  
  
    dt.edge= dt.edge %>%  
      dplyr::mutate(
        sconfZ.x=(conf-scfMean.x)/ifelse(scfSD.x==0,PSUDEO_ZERO_2,scfSD.x), #Zscore based on Mean and MAD
        scontribZ.x=(contrib-scbMean.x)/ifelse( scbSD.x==0,PSUDEO_ZERO_2 ,scbSD.x ) ,
        
        sconfZ.y=(conf-scfMean.y)/ifelse(scfSD.y==0,PSUDEO_ZERO_2,scfSD.y), #Zscore based on Mean and MAD
        scontribZ.y=(contrib-scbMean.y)/ifelse( scbSD.y==0,PSUDEO_ZERO_2 ,scbSD.y ) 
      )
  }
  
  #   Add inxformation theoretic indices -----
  #pr = P(y|X) , conf
  #e=e %>%dplyr::mutate(LgPy_=log(probInp.y),Lg1=log(p_/probInp.y),)
  #e=dt.edge %>%dplyr::mutate(LgPy_=log(probInp.y),Lg1=log(p_/probInp.y))
  #   Imputing infinites dt.edge NO NA should be in data ----
  setDT(dt.edge)
  #a=sapply(dt.edge, function(x) sum(is.na(x)));a[a>0]
  #c=sapply(dt.edge, function(x) sum(is.infinite(x)));c[c>0]
  #b=dt.edge[src=='0166D',.(src,trgt,Weight,OcrInp.x,OcrInp.y,Pin,Pinsum.x)];#View(b)
  
  #=do.call(data.frame,lapply(t, function(x) replace(x, is.infinite(x),NA)))
  # replaceList=list(confZ=0,contribZ=0,confNZ=0,
  #                 contribNZ=0, oddsN=myMax(dt.edge$oddsN ),odds=myMax(dt.edge$odds ),
   # dt.edge=my_replace_val(dt.edge,Inf,replaceList)
  
  #myMax(dt.edge$confZ)
  #RENV3 saved
  
  
  
  #   Finalize calculations on edges -----
  
  dt.edge=dt.edge %>% mutate(DiagnosesAb.x='',DiagnosesAb.y='')
  rm(dt.vertices.back,dt.edge.back,dt.edge.BeforeAugmenting,dt.edge.BeforeAugmenting,dt.edge.noselfedge,n3d.v,dt.vertices1,dt.vertices.BeforeMoments)
  
  
  #   Correct field names x denotes source and y denotes target ----
  names(dt.edge) <- gsub("_x", ".x", names(dt.edge), fixed = TRUE)
  names(dt.edge) <- gsub("_y", ".y", names(dt.edge), fixed = TRUE)
  # names(causal.v) <- gsub("_x", ".x", names(causal.v), fixed = TRUE)
  # names(causal.v) <- gsub("_y", ".y", names(causal.v), fixed = TRUE)
  #dt.edge=dt.edge[,-grep("x$",colnames(dt.edge))]; dt.edge=dt.edge[,-grep("y$",colnames(dt.edge))]
  
  if(url.logfile!="noLog") write('Step5.5',  file = url.logfile, append = TRUE)
  
  
  return(list(rcrd,dt.edge,dt.vertices))
}

#Method could be CICT, and RF or XGB which will only use the basic non CICT features

#' cictTrainTestReport
#'
#' Implements CICT supervised learning and prediction of regulatory edges. Currently heavily depends on global variables
#'
#' @param dt.edge CICT edges produced by prepareEdgeFeatures
#' @param rcrd A list object that accumulates intermediary objects, results and performance measures and will be stored as an rds file for later use
#' @param method Default value is: "CICT", which uses CICT features in a supervised learning. Also can take "RF" (random forest) or "XGB" (xgboost) which these supervised methods on relevance measures
#' @param evaluateUnseenSample Defualt is: TRUE, calculates regulatory edges on a subsample of up to 50000 edges.
#' @param evaluateAllEdges Defualt is: FALSE. If TRUE, CICT tries predicting all potential edges in the given network which might take longer time.
#' @param Debug If TRUE, function enters debugging mode in critical stops along the exectuion allowing verification of variables and features
#' @param preset.train Defualt is: NA. If provided a path to proper CSV, uses that for training. Useful for sensitivity analysis as well as comparision with other methods on similar set of edges/features
#' @param preset.test Defualt is: NA. If provided a path to proper CSV, uses that for training. Useful for sensitivity analysis as well as comparision with other methods on similar set of edges/features
#' @return Returns a list consisted of three objects 
#' # rcrd: is a list object of intermediary objects
#' # edges: a dataframe of edge objects and CICT features for edges
#' # Vertices: a dataframe of vertices objects and CICT features for vertices
#' @examples
#' # Example usage of the function
#' c(rcrd,edges,vertices) %<-z% prepareEdgeFeatures(Debug=Debug)
#' @export
#' 
cictTrainTestReport <-function(dt.edge,rcrd,method='CICT', 
                               evaluateUnseenSample = T, 
                               evaluateAllEdges = F,
                               Debug = F,
                               preset.train=NA, 
                               preset.test=NA){
  #   MIX the Data + Ground truth  ----
  rm(tst1.totalset,t2,tmp,tst1.train,tst1.tst,tst1.tst.h2o,tst1.totalset.h2o,tst1.train.h2o)
  rm(dt.edge.back3,L2,MisLinks,MisNodes,self,rds1,rds2,rds3,outd,others,t2.tmp1,b,e,e.exmpl.c,e.exmpl.e,e.exmpl.r,n.notinVertices,v,t1.c,t1.rc,t1.causal_reversecausal,rpt,t1.rnd,t1,total.rds,tmp.rdndnt)
  rm(tst1.train,tst1.totalset.h2o,tst1.totalset,tst1.train.h2o,tst1.tst,tst1.tst.h2o,tst1.mdl,tst1.rf,tmp.inps,tmp.outps,tmp1)
  gc()
  
  {
    #NrandomEdges =  min(.3 * nrow(dt.edge),200000) ; #NrandomEdges=2e6
    
    
    #or 0 means all
    #dt.edge = dt.edge %>% filter(src %in% geneorder.2e3$X | trgt %in% geneorder.2e3$X)
    
    
    setDF(dt.edge)
    #Test
    # tmp =dt.edge.back0 %>% inner_join(t1,by=c("first"="src","second"="trgt"))
    # t1.abscentGT =t1 %>% anti_join(dt.edge.back0,by=c("src"="first","trgt"="second"))
    # t1.abscentGT =t1.abscentGT %>% anti_join(dt.edge.back0,by=c("src"="second","trgt"="first"))
    
    #dt.edge= dt.edge %>% dplyr::mutate(src=first,trgt=second)
    
    dt.edge.back3 = dt.edge 
    rm(t1,t1.c,t1.rc,t1.c.r,t1.abscentGT,t2.rc,t2.rnd)
    
    #intersect(dt.edge$src,tbl.goldStandard$src) #'CXXC1') 
    
    t1.c = tbl.goldStandard %>% dplyr::select(src,trgt) %>% 
      inner_join(dt.edge,by=c("src"="src","trgt"="trgt"))
    if(nrow(t1.c) < nrow(tbl.goldStandard)/2) warning("Problem in ground truth. More than half of ground truth was not found in edges")
    # t1.c.r = dt.edge %>% inner_join(t1,by=c("src"="trgt","trgt"="src"))
    # t1.c = rbind(t1.c,t1.c.r)
    if(nrow(t1.c) <=0) print("!!! No causal edge in the groundturht? check gene names")
    
    t1.c$predicate = "CAUSES"
    
    t1.rc = dt.edge %>% inner_join(t1.c %>% select(src,trgt),by=c("src"="trgt","trgt"="src"))
    t1.rc$predicate = "REV_CAUSES"
    
    t1.causal_reversecausal=rbind(t1.c,t1.rc)
    
    #Adding 2000 random edges
    t1.rnd = dt.edge %>% #dplyr::mutate(edgetyptruth = NA) %>%
      anti_join(t1.c,by=c("src"="src","trgt"="trgt")) %>%
      anti_join(t1.rc,by=c("src"="trgt","trgt"="src"))
    t1.rnd$predicate = "IRRELEVANTM"
    
    t1 = rbind(t1.c,t1.rc,t1.rnd) #%>% unique()
    #%>%   dplyr::slice_sample(n=NrandomEdges)
    #t1.causal_reversecausal.sub = sample_frac(t1.causal_reversecausal, size = sampling.c_rc.ratio)  
    #if(maxGroundTruth>0) t1.causal_reversecausal = slice_sample(t1.causal_reversecausal,n=maxGroundTruth)
    #tmp = rbind(t1,t1.c,t1.rc)
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
    
    
    
    
    NcausalEdges = min(maxGroundTruth, nrow(t1.c)*minGroundTruth.ratio.learning) %>% as.integer()
    if(NcausalEdges<300)NcausalEdges=floor(2*nrow(t1.c)*minGroundTruth.ratio.learning)
    
    NrandomEdges = NcausalEdges *  randomEdgesFoldCausal #min(sampling.rnd.ratio * nrow(dt.edge),sampling.u.max) %>% as.integer() ; #NrandomEdges=2e6
    
    #If a predifined train and test set is provided
    
    
    if(!is.na(preset.train)){
      t2 = rbind( preset.train %>% select(Gene1,Gene2) %>% inner_join(t1,by=c("Gene1"="src","Gene2"="trgt")),
                  preset.train %>% select(Gene1,Gene2) %>% inner_join(t1,by=c("Gene2"="src","Gene1"="trgt")),
                  preset.test %>% select(Gene1,Gene2) %>% inner_join(t1,by=c("Gene1"="src","Gene2"="trgt"))
      ) %>% rename(src = Gene1, trgt = Gene2)
      
    }else{
      #t2 is the learning set
      t2=rbind(t1 %>% filter(class1=='c') %>% sample_n(size = NcausalEdges),
               t1 %>% filter(class1=='rc') %>% sample_n(size = NcausalEdges),
               t1 %>% filter(class1=='u')%>% sample_n(size=NrandomEdges-2*NcausalEdges)) #to preserve the minGroundTruth.ratio
    }
    t2.complement = t1  %>%  anti_join(t2,by=c("src"="src","trgt"="trgt"))
    
    msg=sprintf(' Network density= %s, learning set density= %s, \n total edges= %s, learning set=%s, test fraction = %s , \n learning set random= %s, learning set causal edges= %s',  
                round(nrow(t1.c)/nrow(dt.edge),4),round(sum(t2$predicate=='CAUSES')/nrow(t2),4),
                prettyNum(nrow(dt.edge),','),prettyNum(nrow(t2),','),tstPrecent,
                round(NrandomEdges,2),round(NcausalEdges,2)
    ) #tstPrecent
    
    cat(msg)
    if(url.logfile!="noLog") write(msg, file = url.logfile, append = TRUE)
    rcrd=c(rcrd, c(density.net=round(nrow(t1.c)/nrow(dt.edge),4),
                   density.learningset=round(sum(t2$predicate=='CAUSES')/nrow(t2),4),
                   edges.total=prettyNum(nrow(dt.edge),','),
                   set.learning=prettyNum(nrow(t2),','),
                   set.learning.testPcnt=tstPrecent,
                   set.learning.rndm=round(NrandomEdges,2),
                   set.learning.causal =round(NcausalEdges,2)))
    
    
    #t2 = rbind(t1.c,t1.rnd)
    #t2$edgeTyp = ''
    t2= as.data.frame(t2) %>% dplyr::select(any_of(c('edgeTyp' ,'edgetyptruth','src', 'trgt','Weight',
                                                     'OcrInp.x', 'OcrInp.y',#intvl_avg,intvl_sd,intvl_median,
                                                     'scfNMAD.y','scbMedian.x',  'ocbNMAD.x','scfL1.x',
                                                     'EI','EO','EE')),everything())
    table(t2$predicate)#,t2$edgeTyp, useNA="ifany")
    
    
    write(table(t2$predicate) %>% knitr::kable(),  file = url.logfile, append = TRUE)
    
    
    a=sapply(t2, function(x) sum(is.na(x)));a[a>0]
    
    library(stringr)
    
    
  }
  
  
  #Early pruning
  if(F){
    dt.edge.back3 = dt.edge  #dt.edge = dt.edge.back3
    #dt.edge = dt.edge %>% filter(src %in% geneorder.2e3$X | trgt %in% geneorder.2e3$X)
    
    rm(t1,t1.c,t1.rc,t1.c.r,t1.abscentGT,t2.rc,t2.rnd,e,tst1.totalset,tst1.tst,tst1.train,t1,t2)
    
    NrandomEdges =  min(.9 * nrow(dt.edge),500000) ; 
    t1.gtfrac = 1
    
    # t1.gt=sc.gtchip %>% dplyr::select(src,trgt,edgetyptruth)
    t1.c = dt.edge %>% inner_join(t1.gt,by=c("src"="src","trgt"="trgt"))%>%
      mutate(predicate = "CAUSES")
    
    t1.rc = dt.edge %>% inner_join(t1.gt,by=c("src"="trgt","trgt"="src")) %>%
      mutate(predicate = "REV_CAUSES")
    
    t1.causal_reversecausal=rbind(t1.c,t1.rc)
    t1.causal_reversecausal.sub = if(t1.gtfrac==1) t1.causal_reversecausal else sample_frac(t1.causal_reversecausal, size = t1.gtfrac)  
    
    t1.rnd = dt.edge %>% dplyr::mutate(edgetyptruth = NA,predicate = "IRRELEVANTM") %>%
      anti_join(t1.causal_reversecausal,by=c("src"="src","trgt"="trgt")) %>% 
      dplyr::slice_sample(n=NrandomEdges) 
    
    
    t1=rbind(t1.causal_reversecausal.sub,t1.rnd)
    rm(t1.rnd,t1.c,t1.rc,t1.causal_reversecausal)
    
    library(stringr)
    a=sapply(t1, function(x) sum(is.na(x)));a[a>0]
    t1$class1=hutils::Switch(t1$predicate, 
                             CAUSES = 'c', PRECEDES= 'el' ,ASSOCIATED_WITH= 'a',
                             REV_CAUSES= 'rc',IRRELEVANTA= 'ir', IF_NA= 'u', DEFAULT = 'u')
    
    t1 = t1 %>% dplyr::mutate(SUID=row_number(),
                              class2=ifelse(class1 %in% c('c'),TRUE,FALSE),
                              class3=ifelse(class1 %in% c('c','rc'),TRUE,FALSE),
                              class5=ifelse(edgetyptruth %in% c('+'),TRUE,FALSE),
                              shared_name=paste0(src,"-",trgt))   #class1=""
    
    setDT(t1)[,`:=`(shared_name=NULL,DiagnosesAb.x=NULL,DiagnosesAb.y=NULL)]
    #t2=unique(t2)
    table(t1$class1,useNA="ifany")
    table(t1$predicate,t1$class2, useNA="ifany")
    table(t1$edgetyptruth,t1$class5, useNA="ifany")
    
    maxGroundTruth = 6000; t2.frac = .1 #or 0 means all
    t2.orig  = t1 #%>% sample_frac(size = t2.frac);table(t2$predicate)
    
  }
  #Test
  # tmp =dt.edge.back0 %>% inner_join(t1,by=c("first"="src","second"="trgt"))
  # t1.abscentGT =t1 %>% anti_join(dt.edge.back0,by=c("src"="first","trgt"="second"))
  # t1.abscentGT =t1.abscentGT %>% anti_join(dt.edge.back0,by=c("src"="second","trgt"="first"))
  
  maxGroundTruth = 6000; t2.frac = .2 #or 0 means all
  
  t2= as.data.frame(t2) %>% dplyr::select(any_of(c('edgeTyp' ,'edgetyptruth','src', 'trgt','Weight',
                                                   'OcrInp.x', 'OcrInp.y',#intvl_avg,intvl_sd,intvl_median,
                                                   'scfNMAD.y','scbMedian.x',  'ocbNMAD.x','scfL1.x',
                                                   'EI','EO','EE')),everything())
  try({
    table(t2$class1,useNA="ifany")
    table(t2$predicate,t2$class2, useNA="ifany")
    table(t2$edgetyptruth,t2$class5, useNA="ifany")
  })
  # write.table(t2,file="Predicates T2 Data100-2.txt",sep = "\t")
  
  # Model col names -----
  require(caret)
  
  t2$trnsparency<-NULL
  
  #print(colList[!colList %in% colnames(dt.vertices)])
  #mdlColNames = colnames(t2)
  mdlColNames=c('corP','corS','corK','mf.mi','tlag1.corP','tlag1.corS','tlag1.corK','tlag1.mi',
                "efMImm","ewMImm","efMIempirical","ewMIempirical", "efMIshrink","ewMIshrink",
                'Indegree.x','Indegree.y','Outdegree.x','Outdegree.y',
                'intvl_median','intvl_kurtosis','intvl_avg','intvl_skew','intvl_sd',
                'ZOcfNMedian.x','ZOcfNMedian.y','ZScbNSkew.x','ZScbNSkew.y',
                'ZOcfSkew.x','ZPoutsum.x','ZOcfSkew.y','ZPoutsum.y','ZScbScbL3.x','ZScbScbL3.y',
                'HTR','tNHTR','NHTRfromSource','NHTRtoTarget',
                'NHTRRatio_fromSource',
                'NHTRRatio_toTarget',
                'EE','NEEfromSource','NEEtoTarget','NEE_NHTR_Ratio',
                'Weight','R','V','P','UR','URR','Tendency',
                'EI','EO','OEER','OERR','AM','HM','GM',
                
                'PAdd', 'PdirRatio',
                'directionLogRatio',
                'hypotenuse','slope',
                'deltaAttitude', 'IRatio',
                't','tn','toe','tRatio', 'tnRatio',
                'toeRatio','tz','tnz',
                'conf','contrib',
                'OcrInp.x','OcrInp.y', 'OcrOut.x2','trgtOcr2',
                'srctrgtSum','srctrgtProduct',
                'confdisc', 'contribdisc', 'cnfxcnfy.emi','cnfxcnty.emi','cntxcnfy.emi','cntxcnty.emi','cnfxcnfy.frcht','cnfxcnty.frcht','cntxcnfy.frcht','cntxcnty.frcht',
                
                'efMImm','ewMImm','efMIempirical','ewMIempirical','efMIshrink','ewMIshrink','Pearson','Spearman','Kendall','Manhattan','Euclidean','L10Norm',
                'conf.fitdist.x','conf.beta1STshape.x','conf.beta2ndshape.x','conf.loglinmean.x','conf.loglinsd.x','cntrb.fitdist.x','cntrb.beta1STshape.x','cntrb.beta2ndshape.x','cntrb.loglinmean.x','cntrb.loglinsd.x',
                'conf.fitdist.y','conf.beta1STshape.y','conf.beta2ndshape.y','conf.loglinmean.y','conf.loglinsd.y','cntrb.fitdist.y','cntrb.beta1STshape.y','cntrb.beta2ndshape.y','cntrb.loglinmean.y','cntrb.loglinsd.y',
                
                'odds','confZ','contribZ','causal1','causal11','causal2','causal3',
                'confN','contribN','oddsN',
                'confNZ','contribNZ',    #'causalN1','causalN11','causalN2','causalN3' ,'causalM','causalM1' ,
                grep("ocf.*",names(t2),value=T),grep("ocb.*",names(t2),value=T),
                grep("scb.*",names(t2),value=T),grep("scf.*",names(t2),value=T),
                grep("Pout.*",names(t2),value=T),grep("Pin.*",names(t2),value=T),
                grep(".*Ocr.*",names(t2),value=T),grep("prob.*",names(t2),value=T))
  
  
  mdlColNames=unique(mdlColNames)
  selectedFeatures = mdlColNames
  mdlColNames=mdlColNames[mdlColNames %in% colnames(t2)]
  
  
  objectiveColNames=c("class1","class2","class3","class5","class13","class14","class16");
  objectiveColNames=objectiveColNames[objectiveColNames %in% colnames(t2)]
  
  evaluationColNames=c('edgeTyp',"DiagnosesAb.x","predicate","DiagnosesAb.y","src","trgt","SUID","wGroup")
  evaluationColNames=evaluationColNames[evaluationColNames %in% colnames(t2)]
  
  colnamesToExport= unique(c(mdlColNames,objectiveColNames, c('srctrgtSum','srctrgtProduct')))
  mdlChkColNames = c(evaluationColNames,objectiveColNames,mdlColNames)
  
  
  mdlColNames[!mdlColNames %in% colnames(t2)]
  colnames(t2)[!colnames(t2) %in% mdlColNames]
  
  mdlChkColNames=mdlChkColNames[mdlChkColNames %in% colnames(t2)]
  
  #sort(mdlChkColNames)
  gpAn = c('ZOcfSkew.x','ZPoutsum.x','ZOcfSkew.y','ZPoutsum.y','HTR','tNHTR','NHTRfromSource','NHTRtoTarget',
           
           'NHTRRatio_fromSource',
           'NHTRRatio_toTarget', 'ocfSkew.x', 'Poutsum.x'  )
  
  
  #Remove large objects
  #sort( sapply(ls(),function(x){object.size(get(x))}),decreasing = F)
  rm(dfl,df,d.new,mtx.lng.prtrb,tzcolvec,tst1.totalset);gc()
  # tst1 ,tst2 Building EXPERIMENT data:  - training and test sets --------  
  require(caret)
  setDT(t2)
  #existingCols = intersect(colnames(t2),c(mdlChkColNames,objectiveColNames))
  tst1.totalset = t2[,intersect(colnames(t2),unique(c(mdlChkColNames,objectiveColNames))),with=FALSE]
  #tst1.totalset=unique(tst1.totalset)
  a=sapply(tst1.totalset, function(x) sum(is.na(x)));a[a>0]
  
  rp = rep(0,length(a[a>0]));names(rp)<-names(a[a>0]);as.list(rp)
  tst1.totalset =replace_na(tst1.totalset, as.list(rp))
  
  try({
    tst1.totalset = na.omit(tst1.totalset)  # cols does not work ,cols=mdlColNames)
    addmargins(table(tst1.totalset$class1,tst1.totalset$edgeTyp));table(tst1.totalset$class1)
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  })
  ntrgtClass =  nrow(tst1.totalset[get(trainingTarget)== TRUE,])
  
  #@@@@@@@@@@@@@@@@@@@@@@@
  if(!is.na(preset.train)){
    tst1.tst =  preset.test %>% select(Gene1,Gene2) %>% inner_join(tst1.totalset,by=c("Gene1"="src","Gene2"="trgt"))
    tst1.train = tst1.totalset %>% anti_join(tst1.tst,by=c("src"="Gene1","trgt"="Gene2"))
                                                                   
  } else{ #Selects a random subset for train and test set
    while(TRUE){
      set.seed(as.integer(runif(1,1,10000)))
      spltIdx =as.vector( caret::createDataPartition(1:nrow(tst1.totalset),p=(1-tstPrecent),list=FALSE,times=1))
      tst1.train = tst1.totalset[spltIdx,] #tst1.dmy[spltIdx,]
      tst1.tst = tst1.totalset[-spltIdx,]
      
      print(nrow(tst1.tst[get(trainingTarget)== TRUE,]))
      if(nrow(tst1.tst[get(trainingTarget)== TRUE,]) >= (tstPrecent-0.01)* ntrgtClass) break
      # 
      # spltIdx =as.vector( caret::createDataPartition(1:nrow(dt.vertices),p=(1-tstPrecent),list=FALSE,times=1))
      # trainV = paste0('G',spltIdx)
      # tst1.train = setDT(tst1.totalset)[src %in% trainV | trgt %in% trainV,] #tst1.dmy[spltIdx,]
      # tst1.tst = tst1.totalset[!(src %in% trainV | trgt %in% trainV),]
      
    }
  }
  
  if(url.logfile!="noLog") {
    write(paste0('tst1.train$class2 : ',table(tst1.train$class2)),  file = url.logfile, append = TRUE)
    write(paste0('tst1.tst$class2 : ',table(tst1.tst$class2)),  file = url.logfile, append = TRUE)
  }
  
  
  tmp = table(tst1.train$class2); names(tmp)<-paste0('train.',names(tmp))
  rcrd=c(rcrd, tmp)
  tmp = table(tst1.tst$class2); names(tmp)<-paste0('test.',names(tmp))
  rcrd=c(rcrd, tmp)
  
  if(FLAG_exportTrainAndTest){
    try({

      tst1.train %>% select(src , trgt, class2, class3) %>% 
        rename(Gene1 = src, Gene2 = trgt, Type=class2,Association=class3) %>% 
        fwrite(file=paste0(url.outputFolder,'train.csv'),row.names = F, sep='\t')
      tst1.tst %>% select(src , trgt, class2, class3) %>% 
        rename(Gene1 = src, Gene2 = trgt, Type=class2,Association=class3) %>% 
        fwrite(file=paste0(url.outputFolder,'test.csv'),row.names = F, sep='\t')
    })
  }
  
  table(tst1.train$class2);table(tst1.tst$class2)
  table(tst1.train$class3);table(tst1.tst$class3)
  table(tst1.train$class5);table(tst1.tst$class5)
  
  
  #save(dt.edge,dt.vertices,t2,tst1.totalset,mdlChkColNames,mdlColNames,file="CICT dream4-100-1-1 MF two way.robj")
  #load(file="CICT dream4-100-2.robj")
  # TRAINING H2O  For Prediction and feature selection ----
  #Remove large objects
  if(url.logfile!="noLog") write('Step5.6',  file = url.logfile, append = TRUE)
  #sort( sapply(ls(),function(x){object.size(get(x))}),decreasing = F)
  
  #require(retry)
  require(h2o)
  #try({h2o.shutdown(prompt=FALSE)})
  gc()
  
  #retry::retry({
  H2OCnn = h2o.init(nthreads = parallel::detectCores(), enable_assertions = TRUE,
                    max_mem_size = MAX_MEM_SIZE,strict_version_check=FALSE) 
  #},when='OverflowError|54321: Connection refused',interval=15,max_tries=3) 
  
  #memory.limit(size=25000)
  #Create binary outcome to calculate probaibility
  Sys.sleep(5)
  #trainingTarget;table(tst1.totalset[,get(trainingTarget)])
  
  NFolds = 5
  table(setDT(tst1.totalset)[,get(trainingTarget)])
  filterOutFeatures =  c('srctrgtSum','srctrgtProduct') #Weight
  #mdlColNames = intersect(mdlColNames, prd.varimp$variable[1:150])
  mdlColNames=setdiff(mdlColNames,filterOutFeatures)
  d.new = tst1.tst # sample_frac(tst1.tst,size = .3)
  tst1.tst.h2o=as.h2o(d.new) 
  tst1.train.h2o = as.h2o(tst1.train)
  #tst1.totalset.h2o = as.h2o(tst1.totalset)
  rm(tst1.rf,tst1.mdl)
  
  #optimize for causal-versus with ntrees=10,max_depth=10, training=75%
  #OPTIMZE FOR dream4ds2 causal versus random  max_depth =7,ntrees=8,
  tst1.rf = h2o.randomForest(mdlColNames,trainingTarget,
                             tst1.train.h2o,nfolds = NFolds,#checkpoint=tst2.mdl@model_id,
                             max_depth =ifelse(is.null(RF_max_depth) ,20,RF_max_depth ), 
                             ntrees=ifelse(is.null(RF_ntrees) ,50,RF_ntrees ), 
                             #max_depth =9,ntrees=12, 
                             keep_cross_validation_predictions= FALSE,  validation_frame = tst1.tst.h2o)
  #,ntrees=30,max_depth=5#selectedFeatures #mdlColNames
  
  tst1.mdl=tst1.rf
  
  mtrc= tst1.rf@model$validation_metrics@metrics
  
  rcrd$AUC = mtrc$AUC
  rcrd$pr_auc = mtrc$pr_auc
  rcrd$MSE = mtrc$MSE
  rcrd$RMSE = mtrc$RMSE
  
  mtrc.idx.f1=1
  mtrc.idx.accuracy=3
  mtrc.idx.mcc=8
  mtrc$max_criteria_and_metric_scores[c(mtrc.idx.f1,mtrc.idx.accuracy,mtrc.idx.mcc),]
  
  rcrd$max_criteria = mtrc$max_criteria_and_metric_scores[c(mtrc.idx.f1,mtrc.idx.accuracy,mtrc.idx.mcc),]
  rcrd$cmdesc = 'Validation data, confusion matrix. Row labels: Actual class; Column labels: Predicted class'
  rcrd$cm = mtrc$cm$table %>% data.frame()
  
  
  # tst1.gbm = h2o.gbm(mdlColNames,trainingTarget,
  #                    tst1.train.h2o,nfolds = NFolds,#checkpoint=tst2.mdl@model_id,
  #                    max_depth =3, 
  #                    #ntrees=12, 
  #                    keep_cross_validation_predictions= FALSE,  validation_frame = tst1.tst.h2o ) #ntrees=20
  # tst1.mdl = tst1.gbm
  
  mdl.class2 = tst1.mdl
  
  if(url.logfile!="noLog")  write('Step5.7 Model trained successfuly',  file = url.logfile, append = TRUE)
  
  #save(tst1.mdl,t1,t2,tst1.totalset,tst1.train,tst1.tst,file="mdl1.robj")
  
  # 25000 random sample & bestCutoff calculations -----------------
  {
    
    #Predictalldata
    if(FALSE){
      setkey(t2,src,trgt)
      setkey(dt.edge,src,trgt)
      d.new = dt.edge[!t2]
      d.new[,(trainingTarget) := ifelse(is.na(get(trainingTarget)),0,get(trainingTarget))]
    }
    
    #Predicts in smaller chunks also adds random prediction for comparition
    msg = c('================================================',
            "Reporting results on unseen sample")
    cat(paste0(msg,collapse='\n') )
    write( msg,  file = url.logfile, append = TRUE, sep = '\n')
    
    d.new = t2.complement %>% sample_n(size = min(50000,nrow(t2.complement) * maxunseenTest.ratio))   #tst1.totalset #newDataEnv$tst1.totalset # 
    table(d.new$predicate)
    
    msg = sprintf("Learning set= %s | training= %s | validation= %s | unseen sample = %s | t2 + comp = %s | all Edges= %s",  
                  nrow(tst1.totalset),nrow(tst1.train),nrow(tst1.tst),nrow(d.new), nrow(t2)+nrow(t2.complement),nrow(dt.edge))
    
    cat(msg)
    if(url.logfile!="noLog") write(msg, file = url.logfile, append = TRUE)
    rcrd$unseensmpl_stats = msg
    
    
    #d.new = dt.edge
    splitcount = 2
    splts = split(1:nrow(d.new),          
                  cut(seq_along(1:nrow(d.new)),
                      splitcount,labels = FALSE))
    
    #Creates and saves random predictions for all edges
    randomPredictions = dt.edge %>% select(src,trgt) %>% mutate(rndPred = ifelse(runif(nrow(dt.edge)) >=.5,1,0))
    
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
    
    d.new1= rbindlist(d.new.tmp) %>% left_join(randomPredictions,by=c('src'='src','trgt'='trgt'))
    
    
    d.new1.rv =d.new1 %>% dplyr::rename(revWeight = Weight,rvpred=predictions,src1=src,trgt1=trgt) %>% 
      dplyr::select(-outcomes,-rndPred)
    pred_outcome = d.new1 %>% left_join(d.new1.rv, by=c("src"="trgt1", "trgt"="src1") )
    pred_outcome.back = pred_outcome
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
        
        vrtcs = unique(dt.vertices$vID)
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
    
    # Testing performance -----
    {
      prd.varimp=h2o.varimp(tst1.mdl)#;View(prd.varimp)
      
      rcrd$varimp = prd.varimp[1:20,] %>% as.data.frame()
      if(url.logfile!="noLog")  write( prd.varimp[1:20,]%>% knitr::kable(),  file = url.logfile, append = TRUE, sep = '\n')
      #plot(tst1.mdl);setDF(tst1.totalset)
      #sink(NULL)
      #sink(file = paste0(url.logfile,'_t.txt'), append = T, type = c("output"),split = T) #write(  file = url.logfile, append = TRUE)
      h2o.performance(tst1.mdl,valid=T) #save(tst1,tst1.train,tst1.tst, file="temp training set.robj")
      
      msg = c('================================================',
              "Reporting model performance on validation set",
              capture.output(h2o.performance(tst1.mdl,valid=T)))
      cat(paste0(msg,collapse='\n') )
      if(url.logfile!="noLog") write( msg,  file = url.logfile, append = TRUE, sep = '\n')
      
      if(F){
        tstexplain = tst1.totalset # t1 %>% anti_join(t2,by=c("src"="src", "trgt"="trgt")) %>% sample_frac(size =.02)
        h2o.explain(tst1.mdl,as.h2o(tstexplain))
        h2o.ice_plot(tst1.mdl,as.h2o(tstexplain),'scftau4.x')
      }
      #paste0(prd.varimp$variable[1:30],collapse= "','")
      #prd.varimp[1:15,] %>% knitr::kable()
      
      #h2o.saveModel(tst1.mdl,"mdlD100-1") #h2o.loadModel("E:/HDDs/Work/0 Plant biology/R code/mdlD100-1")
      #library(broom); tidy(tst1.mdl); glance(tst1.mdl)
      
      
      
      setDF(tst1.totalset)
      
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
          mutate(aucs = round(aucs,3), standardized = round(standardized,3) ) 
      }
      #newDataEnv = new.env()
      #R.utils::loadToEnv(file="CICT dream4-100-2.robj",envir=newDataEnv)
      
      
      #TODO remove the reverse edges before assessment, just keep the one with higher prediction score
      assespreds = pred_outcome #%>% dplyr::filter(predictions>=rvpred)
      theROC <- roc(assespreds$outcomes, assespreds$predictions, percent = TRUE);theROC
      
      sscurves <- evalmod(scores = assespreds$predictions, labels = assespreds$outcomes)
      pr.crv=sscurves$prcs[1][[1]];  pr.auc =  attr(pr.crv,'auc')
      
      theauc = precrec::auc(sscurves); 
      msg = paste0(theauc[[3]],"=", round(theauc[[4]],3))
      if(url.logfile!="noLog") write(msg,  file = url.logfile, append = TRUE)
      ##autoplot(sscurves, mode = 'rocpr',ret_grob = F)
      rcrd$unseensmpl_roc_pr = msg
      
      
      #partial precision-recall
      sscurves.part <- part(sscurves, xlim = c(0, 0.2))
      ##plot(sscurves.part,title = "Partial top 20% predictions")
      if(url.logfile!="noLog"){ 
        write('Early precision recall ======',  file = url.logfile, append = TRUE)
        write(reportAUC(sscurves.part) %>% kable(),  file = url.logfile, append = TRUE)
      }
      rcrd$unseensmpl_part = as.data.frame(reportAUC(sscurves.part) )
      
      pr.prtcrv=sscurves.part$prcs[1][[1]];  pr.prtauc =  attr(pr.prtcrv,'pauc')
      
      #sink(file = url.logfile, append = T, type = c("output"),split = T) #write(  file = url.logfile, append = TRUE)
      
      #random classifier
      randomClassifierCurves <- evalmod(scores = assespreds$rndPred, labels = assespreds$outcomes)
      rnd.crv=randomClassifierCurves$prcs[1][[1]]; rnd.auc =  attr(rnd.crv,'auc')
      
      rndmClscurves.part <- part(randomClassifierCurves, xlim = c(0, 0.2))
      reportAUC(rndmClscurves.part) 
      rnd.prtcrv=rndmClscurves.part$prcs[1][[1]]; rnd.prtauc =  attr(rnd.prtcrv,'pauc')
      rcrd$unseensmpl_rndm = as.data.frame(reportAUC(rndmClscurves.part) )
      
      msg=""
      if(url.logfile!="noLog") {
        print("Random Classifier comparison ============")
        write("Random Classifier comparison ============",  file = url.logfile, append = TRUE)
        
        
        write('Random classifier precision recall ======',  file = url.logfile, append = TRUE)
        write(reportAUC(rndmClscurves.part) %>% kable(),  file = url.logfile, append = TRUE)
        
        msg = sprintf("%s , AUCPR Ratio CICT to Random= %s,  top 20 percent AUCPR Ratio CICT to Random= %s ", 
                      arg.dname, round(pr.auc/rnd.auc,2), round(pr.prtauc/rnd.prtauc,2))
        write(msg,  file = url.logfile, append = TRUE)
      };  print(msg)
      
      mmpoins <- evalmod(scores = assespreds$predictions, labels = assespreds$outcomes, mode = "basic")
      #autoplot(mmpoins, c("error", "accuracy"))
      
      # Show normalized ranks vs. specificity, sensitivity, and precision
      ##autoplot(mmpoins, c("specificity", "sensitivity", "precision"))
      
      # Show normalized ranks vs. Matthews correlation coefficient and F-score
      ##autoplot(mmpoins, c("mcc", "fscore"))
      
      #bestcutoff = coords(roc, "best", ret="threshold", transpose = FALSE)
      
      # the relative cost of of a false negative classification (as compared with a false positive classification)
      # the prevalence, or the proportion of cases in the population (n.cases/(n.controls+n.cases)).
      relativeCostfn_fp = 1/2
      rcrd$relativeCostfn_fp=relativeCostfn_fp
      
      prv = table(setDT(tst1.totalset)[,get(trainingTarget)])
      best.weights=c(relativeCostfn_fp, prv[2]/prv[1])  
      bestcutoff =as.double(coords(theROC, "best", best.method = "closest.topleft", best.weights=best.weights,
                                   ret="threshold", transpose = FALSE));bestcutoff
      
      
      #bestcutoff = 0.5
      assespreds =assespreds %>% mutate(thresholdpreds= ifelse(assespreds$predictions>bestcutoff,assespreds$predictions,0)) 
      theROC <- roc(assespreds$outcomes, assespreds$thresholdpreds, percent = TRUE);theROC
      
      
      
      sscurves <- evalmod(scores = assespreds$thresholdpreds, labels = assespreds$outcomes)
      pr.crv=sscurves$prcs[1][[1]];  pr.auc =  attr(pr.crv,'auc')
      
      msg=""
      if(url.logfile!="noLog") {
        print(sprintf("Best threshold: %s for fn/fp cost ratio: %s " ,bestcutoff,relativeCostfn_fp ))
        write(sprintf("Best threshold: %s for fn/fp cost ratio: %s " ,bestcutoff,relativeCostfn_fp ),  file = url.logfile, append = TRUE)
        
        msg = sprintf("With FP-FN ratio 0.5 => AUCPR Ratio CICT to Random= %s,  top 20 percent AUCPR Ratio CICT to Random= %s ", 
                      round(pr.auc/rnd.auc,2), round(pr.prtauc/rnd.prtauc,2))
        write(msg,  file = url.logfile, append = TRUE)
        
      }
      
      print(msg)
      rcrd$costFnFp_cutoff_roc = msg
      #sink(NULL)
      # 
      # pander::pander(ftable(factor(assespreds$thresholdpreds, c("TRUE","FALSE")), 
      #                       factor(assespreds$outcomes, c("TRUE","FALSE")), 
      #                       dnn=c('pred','actual'))) #%>%  knitr::kable()
      
      #table(ifelse(predictions>.99,'c','other'))
      #bestcutoff = coords(roc, "best", best.method = "closest.topleft", ret="threshold", transpose = FALSE);bestcutoff
      
      
      
    }
  }
  
  # Final run on all edges ------
  
  rcrd$learning_curve_plot_auto = h2o.learning_curve_plot(tst1.mdl,metric = c(  "Auto"))
  rcrd$learning_curve_plot_aucpr = h2o.learning_curve_plot(tst1.mdl,metric = c(  "aucpr"))
  rcrd$learning_curve_plot_auc = h2o.learning_curve_plot(tst1.mdl,metric = c(  "auc"))
  rcrd$learning_curve_plot_logloss = h2o.learning_curve_plot(tst1.mdl,metric = c(  "logloss"))
  rcrd$learning_curve_plot_amisclassification = h2o.learning_curve_plot(tst1.mdl,metric = c(  "misclassification"))
  rcrd$learning_curve_plot_rmse = h2o.learning_curve_plot(tst1.mdl,metric = c(  "rmse"))
  
  #FLAG_runOnAllEdges=T;FLAG_exportRankedEdges=T
  if(FLAG_runOnAllEdges) {
    write('Step 6, edge ranking',  file = url.logfile, append = TRUE)
    d.new = dt.edge # %>% sample_frac(size =.01)  # %>% anti_join(t2,by=c("src"="src", "trgt"="trgt"))  #tst1.totalset #newDataEnv$tst1.totalset # 
    #d.new = dt.edge
    splitcount = max(floor(nrow(dt.edge) / 30000),2)
    splts = split(1:nrow(d.new),          
                  cut(seq_along(1:nrow(d.new)),
                      splitcount,labels = FALSE))
    
    
    # library("parallel")
    # cl <- makeCluster(getOption("cl.cores", parallel::detectCores()))
    
    d.new.tmp = lapply(#cl = cl,
      splts,
      function(thesplit){
        
        print(last(thesplit))
        d.new.slice = setDT(d.new)[thesplit,]
        predTest.d.h2o = as.h2o(d.new.slice) # 
        h2o.pred = as.data.frame(h2o.predict(tst1.mdl,predTest.d.h2o,keep_cross_validation_predictions=F))
        h20.prediction=as.numeric(as.character(h2o.pred[,3]))
        predictions =h20.prediction #pred #ens.predictions #
        #outcomes =unlist(setDT(d.new.slice)[,trainingTarget,with=FALSE]) #outcome # ens.outcome# 
        set.seed(runif(1,1,1000))
        rndPred = ifelse(runif(nrow(d.new.slice)) >=.5,1,0) #Assign a random classifier results
        prd_outcomes = d.new.slice %>% select(src,trgt,any_of('Weight')) %>%
          cbind(predictions,rndPred) %>% as.data.frame()
        prd_outcomes
      }) #,.parallel = cl
    d.new1= rbindlist(d.new.tmp)
    
    write('Step 6.1, edges ranking finished',  file = url.logfile, append = TRUE)
    
    d.new1.rv =d.new1 %>% dplyr::rename(revWeight = Weight,rvpred=predictions,src1=src,trgt1=trgt) %>% 
      dplyr::select(-rndPred)
    pred_outcome = d.new1 %>% left_join(d.new1.rv, by=c("src"="trgt1", "trgt"="src1") )
    
    if(FLAG_exportRankedEdges)
    {
      pred_outcome.e=pred_outcome %>% rename(Gene1=src,	Gene2=trgt,	EdgeWeight = predictions) %>%
        select(Gene1,Gene2,EdgeWeight,everything()) %>% arrange(desc(EdgeWeight))

      if(Debug) browser()  
      
      
    if(!file.exists(url.rankedEdges) | forceOutput){
      file.remove(url.rankedEdges)
      try({fwrite(pred_outcome.e,url.rankedEdges,row.names = F, sep='\t')},silent=T)
      
      rcrd$rankededges_count = nrow(pred_outcome.e)
      rcrd$rankededges_gated_count = nrow(pred_outcome.e)
      
      
      truncated_pred_outcome = pred_outcome.e %>% 
        mutate(EdgeWeight=ifelse(EdgeWeight>bestcutoff,EdgeWeight,0))
      rcrd$rankededges_gated_count = nrow(truncated_pred_outcome)
      rcrd$rankededges_gated_desc = paste0('EdgeWeight > ',bestcutoff )
      rcrdrankededges_url =url.rankedEdges 
      
      try({fwrite(truncated_pred_outcome,url.rankedEdgesGated,row.names = F, sep='\t')})
      #try({write.csv(randomPredictions ,file= url.randomPreds,row.names=F)})
      
      }
      
    }
    
   
    
    # All results evaluation ==========================
    {
      
      assespreds = t1 %>% select(src,trgt,matches('class'))  %>% 
        right_join(pred_outcome, by=c("src"="src", "trgt"="trgt")  ) %>% 
        mutate(prdEdgeWeight=ifelse(predictions>bestcutoff,predictions,0),
               outcomes = class2) %>%
        select(src,trgt,outcomes,predictions,prdEdgeWeight,rndPred,revWeight,rvpred,matches('class'))
      
      theROC <- roc(assespreds$outcomes, assespreds$predictions, percent = TRUE);theROC
      
      sscurves <- evalmod(scores = assespreds$predictions, labels = assespreds$outcomes)
      pr.crv=sscurves$prcs[1][[1]];  pr.auc =  attr(pr.crv,'auc')
      # rcrd$sscurves.roc = autoplot(sscurves,curvetype = 'ROC')
      rcrd$sscurves.prc = autoplot(sscurves,curvetype = 'PRC')
      
      theauc = precrec::auc(sscurves); 
      msg = paste0(theauc[[3]],"=", round(theauc[[4]],3))
      if(url.logfile!="noLog") write(msg,  file = url.logfile, append = TRUE)
      ##autoplot(sscurves, mode = 'rocpr',ret_grob = F)
      rcrd$allEdges_roc_pr = msg
      
      
      #partial precision-recall
      sscurves.part <- part(sscurves, xlim = c(0, 0.2))
      # rcrd$sscurves.part.roc = autoplot(sscurves.part,curvetype = 'ROC')
      rcrd$sscurves.part.prc = autoplot(sscurves.part,curvetype = 'PRC')
      
      ##plot(sscurves.part,title = "Partial top 20% predictions")
      if(url.logfile!="noLog"){ 
        write('Early precision recall ======',  file = url.logfile, append = TRUE)
        write(reportAUC(sscurves.part) %>% kable(),  file = url.logfile, append = TRUE)
      }
      rcrd$allEdges_part = as.data.frame(reportAUC(sscurves.part) )
      
      pr.prtcrv=sscurves.part$prcs[1][[1]];  pr.prtauc =  attr(pr.prtcrv,'pauc')
      
      #sink(file = url.logfile, append = T, type = c("output"),split = T) #write(  file = url.logfile, append = TRUE)
      
      #random classifier
      randomClassifierCurves <- evalmod(scores = assespreds$rndPred, labels = assespreds$outcomes)
      rnd.crv=randomClassifierCurves$prcs[1][[1]]; rnd.auc =  attr(rnd.crv,'auc')
      
      rndmClscurves.part <- part(randomClassifierCurves, xlim = c(0, 0.2))
      reportAUC(rndmClscurves.part) 
      rnd.prtcrv=rndmClscurves.part$prcs[1][[1]]; rnd.prtauc =  attr(rnd.prtcrv,'pauc')
      rcrd$allEdges_rndm = as.data.frame(reportAUC(rndmClscurves.part) )
      
      msg=""
      if(url.logfile!="noLog") {
        print("Random Classifier comparison ============")
        write("Random Classifier comparison ============",  file = url.logfile, append = TRUE)
        
        
        write('Random classifier precision recall ======',  file = url.logfile, append = TRUE)
        write(reportAUC(rndmClscurves.part) %>% kable(),  file = url.logfile, append = TRUE)
        
        msg = sprintf("AUCPR Ratio CICT to Random= %s,  top 20 percent AUCPR Ratio CICT to Random= %s ", 
                      round(pr.auc/rnd.auc,2), round(pr.prtauc/rnd.prtauc,2))
        write(msg,  file = url.logfile, append = TRUE)
      };  print(msg)
      
      mmpoins <- evalmod(scores = assespreds$predictions, labels = assespreds$outcomes, mode = "basic")
      #autoplot(mmpoins, c("error", "accuracy"))
      
      # Show normalized ranks vs. specificity, sensitivity, and precision
      ##autoplot(mmpoins, c("specificity", "sensitivity", "precision"))
      
      # Show normalized ranks vs. Matthews correlation coefficient and F-score
      ##autoplot(mmpoins, c("mcc", "fscore"))
      
      #bestcutoff = coords(roc, "best", ret="threshold", transpose = FALSE)
      
      # the relative cost of of a false negative classification (as compared with a false positive classification)
      # the prevalence, or the proportion of cases in the population (n.cases/(n.controls+n.cases)).
      relativeCostfn_fp = 1/2
      rcrd$relativeCostfn_fp=relativeCostfn_fp
      
      prv = table(setDT(tst1.totalset)[,get(trainingTarget)])
      best.weights=c(relativeCostfn_fp, prv[2]/prv[1])  
      bestcutoff =as.double(coords(theROC, "best", best.method = "closest.topleft", best.weights=best.weights,
                                   ret="threshold", transpose = FALSE));bestcutoff
      
      #bestcutoff = 0.5
      assespreds =assespreds %>% mutate(thresholdpreds= ifelse(assespreds$predictions>bestcutoff,assespreds$predictions,0)) 
      theROC <- roc(assespreds$outcomes, assespreds$thresholdpreds, percent = TRUE);theROC
      
      
      
      sscurves <- evalmod(scores = assespreds$thresholdpreds, labels = assespreds$outcomes)
      pr.crv=sscurves$prcs[1][[1]];  pr.auc =  attr(pr.crv,'auc')
      # rcrd$sscurves.gated.roc = autoplot(sscurves,curvetype = 'ROC')
      # rcrd$sscurves.gated.prc = autoplot(sscurves,curvetype = 'PRC')
      
      msg=""
      if(url.logfile!="noLog") {
        print(sprintf("Best threshold: %s for fn/fp cost ratio: %s " ,bestcutoff,relativeCostfn_fp ))
        write(sprintf("Best threshold: %s for fn/fp cost ratio: %s " ,bestcutoff,relativeCostfn_fp ),  file = url.logfile, append = TRUE)
        
        msg = sprintf("With FP-FN ratio 0.5 => AUCPR Ratio CICT to Random= %s,  AUCPR Ratio CICT to Random= %s ", 
                      round(pr.auc/rnd.auc,2), round(pr.prtauc/rnd.prtauc,2))
        write(msg,  file = url.logfile, append = TRUE)
        
      }
      
      print(msg)
      rcrd$bestcutoff = bestcutoff
      rcrd$allEdges_costFnFp_cutoff_roc = msg
      #sink(NULL)
      # 
      # pander::pander(ftable(factor(assespreds$thresholdpreds, c("TRUE","FALSE")), 
      #                       factor(assespreds$outcomes, c("TRUE","FALSE")), 
      #                       dnn=c('pred','actual'))) #%>%  knitr::kable()
      
      #table(ifelse(predictions>.99,'c','other'))
      #bestcutoff = coords(roc, "best", best.method = "closest.topleft", ret="threshold", transpose = FALSE);bestcutoff
      
      
    }  

  } 
  
  if(Debug) browser() 
  #exporting network objects for further evaluations
  if(Debug) saveRDS(list(rcrd,prd.varimp,dt.edge,pred_outcome.e,,t1),file = paste0(url.output , '/CICT network objects.rds'))
  
  print('Data produced successfuly ==================================')
  return(rcrd)
  
  #Visualize
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
  
  
  ################################################################################@
  # © 2016 Abbas Shojaee <abbas.shojaee@gmail.com>
  # Causal Inference Using Composition of Transactions CICT
  # All rights reserved, please do not use or distribute without written permission
  ################################################################################@
  
  ################################################################################@
  # © 2021 Abbas Shojaee <abbas.shojaee@gmail.com> - Carol Huang Lab, NYU
  # Causal Inference Using Composition of Transactions CICT
  # All rights reserved, please do not use or distribute without written permission
  ################################################################################@
  
  
  
}
