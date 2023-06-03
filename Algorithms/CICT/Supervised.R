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
# Function_definitions -----

#Remove large objects
sort( sapply(ls(),function(x){object.size(get(x))}),decreasing = F)
#rm(thedist,e,recipe_1,t1.rnd,n.itm.e.BeforeAugmenting)
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
  # hashs =  setDT(hashs)[hashs,, on = "hash"][cls != i.cls,]
  # dups = hashs[!(cls %in% excluded | i.cls %in% excluded) & dups==T,]$i.cls
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

getComorbCategoryfromLabel=function(dxc){
  dxc = gsub("DXCCS","", dxc)
  ctgry=as.numeric(dxc)
  return(ctgry)
}

getDxCategory=function(dxccs){
  ctgry=as.numeric(substr(dxccs,0,regexpr("[.]",dxccs)-1),"")
  return(ctgry)
}

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
  #browser()
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

#Finds modes of a multimodal distribution or you can use discretize with expectation maxima
'method = "mode" [default]: calculates the mode for a unimodal vector, else returns an NA
    method = "nmodes": calculates the number of modes in the vector
    method = "modes": lists all the modes for a unimodal or polymodal vector'

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

#if result of sd is NA replace it with default like when the number of elements in x are less than 2
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

#df:data frame, rplist, value to replace # tmp2=my_replace_val(tmp1,list(confZ=myMax(tmp1$confZ)),Inf,)
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
  # xs=str_split(x,"[, ]")
  # xf=lapply(xs,function(strWords) toupper(str_extract_all(strWords,"^.")))
  # xf=Filter(Negate(is.na), xf)
  # xr=lapply(xs,function(strWords){
  #   strWords=unlist(strWords)  #strWords=str_split(x[14],"[, ]")
  #   strWords1= ifelse(str_length(strWords)>ln, str_replace(strWords, "([aieou???]{2})",""),strWords)
  #   strWords2=unlist(str_extract_all(strWords1,paste0("(?<=^.).{1,",ln,"}")))
  #   #strWords3 =Filter(Negate(is.na), strWords2)
  #   strWords3=strWords2[!(is.na(strWords2)| identical(strWords2, character(0)))]
  #   strWords3
  #   })
  # xC=lapply(1:length(xs),function(i){ stf=unlist(xf[i]);st=unlist(xr[i]);paste0(stf,st)})
  f=unlist(lapply(xC,function(strAbs) str_sub( paste(strAbs,collapse=collapse),1,maxLength)))
}


#defults
extractLmoments = function(v)
{
  require(Lmoments)
  if(length(v)==1){
    # if just one value reutrns defualst for normal distirbution
    #Normal 	L1=mean, L2=	σ /√π ,tau3=l_skewness=	0 ,tau4=L_kurtosis=	0.1226
    #browseURL("https://en.wikipedia.org/wiki/L-moment")
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

# definition from package fGrach Has error in calculation like with p=.75 and xi= 0.007588785 returns NA
#use sn::qsc instead qsc(p, xi = 0, omega = 1, alpha = 0, dp = NULL)
#xi vector of location parameters.  = mean
#omega vector of (positive) scale parameters. = SD
#alpha vector of slant parameters = Skewness

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



printNiceTable<- function(xtablePrint,output,css){
  'The idea is to :
  
      Create a css where you format "stylize" your table using some css features
      Create a html table using print.xtable
      Create a file including a link to the css file and the created html table
      So here the code creating the "res.html" file:'
  
  require(xtable)
  ## here I am using a link to mystyle.css
  html.head <- paste("<head>" ,
                     paste0('<link rel="stylesheet" type="text/css" href=',css,'/>'),
                     "</head>",sep='\n')
  ## the html body
  html.table <- paste(xtablePrint, collapse = "\n")
  html.body <- paste('<body><div class="CSS_Table_Example"> ', html.table,"</div></body>")
  ## the html file
  write(paste(html.head,html.body,sep='\n'),output)
  
}


silhouette.plot=function(clstmethod,clst.avg.silhouette,title){
  #sil=silhouette(clstmethod,distm)
  #plot(sil,cex.names=0.6,nmax=98,main="Silhouette Plot", col=c("red","green", "blue", "purple", "yellow"))
  bndry=(max(clst.avg.silhouette)-min(clst.avg.silhouette))*0.5
  barplot(clst.avg.silhouette, col = rainbow(20),xlab="cluster Number", ylab="Value of Silhouette",
          main=paste("Average of Silhouette Plot",title),ylim=c((min(clst.avg.silhouette)-bndry),(max(clst.avg.silhouette)+bndry)))
}

plot.parallelcord=function(clust.method,clust.size,title){
  clst.df=as.data.frame(p.node) %>%dplyr::mutate(clst=clust.method)
  m=as.data.frame(ddply(clst.df,.(clst),summarize,NDXMed=mean(NDXMed),NDXHigh=mean(NDXHigh),NDXLow=mean(NDXLow)))
  ggparallel(
    list("NDXHigh","NDXLow","NDXMed"),
    text.offset = c(0.03, 0,-0.03),
    data = m,
    width = 0.1,
    order = c(1, 0),
    angle = 0,
    color = "white",
    ylab="Cluster Number",
    main=paste("Parallel Coordinate: Attributes in Clusters",title)
  )
  
  ggparcoord(data = m,columns = c(2:4),groupColumn = 1,showPoints=TRUE,shadeBox = 7,alphaLines=1,
             scale = "uniminmax",boxplot = TRUE,title = paste("Parallel Coordinate: Attributes in Clusters",title))
  #,mapping = ggplot2::ylab(label = "Cluster Number"))
}

Scatter.Plot=function(clust.method,clust.size,start,end,title){
  super.sym <- trellis.par.get("superpose.symbol")
  clst.df=as.data.frame(p.node) %>%dplyr::mutate(clst=clust.method)
  splom(~clst.df[start:end], groups = clst, data = clst.df,panel = panel.superpose,
        key = list(title = paste("The Varieties of Clustering ",title), columns = clust.size,
                   points = list(pch = super.sym$pch[1:clust.size],col = super.sym$col[1:clust.size])),
        text = list(c(1:clust.size)))
  #parallelplot(~clst.df[start:end] | clst, clst.df)
}

myplotViolin<-function(theTitle,dt,variables, group1, group2,
                       grp1title,grp2title,grpcolor1 = "green",grpcolor2 = "blue", 
                       yval = 'value')
{
  #dataframe in long format
  dt=setDT(dt)[name %in% variables,]
  
  fig <- setDT(dt) %>%
    plot_ly(type = 'violin') 
  
  fig <- fig %>%
    add_trace(
      x = ~dt[class10,get('name')],
      y = ~dt[class10,get(yval)],
      legendgroup = grp1title,
      scalegroup = grp1title,
      name = grp1title,
      side = 'negative',
      box = list(
        visible = F
      ),
      meanline = list(
        visible = F
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
        visible = F
      ),
      meanline = list(
        visible = F
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

read.tsv <- function(..., header=F, sep='\t', quote='', comment='', na.strings='', stringsAsFactors=FALSE) {
  # read.table() wrapper with default settings for no-nonsense, pure TSV
  # Typical use case is output from another program.
  # (R's defaults are more geared for human-readable datafiles, which is less
  # feasible for large-scale data anyway.)
  # These options are substantially faster than read.table() defaults.
  #   (see e.g. LINK)
  # stringsAsFactors is the devil.
  
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

translateVars<-function(d,vars,setLargestGroupAsReference=TRUE,replaceUnknownLevelsWithNA=FALSE,warnMissing=FALSE){
  #browseURL("https://it.unt.edu/interpreting-glm-coefficients") #interpreting reference group
  library(plyr)
  if(nrow(d)<=1) return(NULL)
  #Assumes all missing levels got a negative value or alphabetic
  addMissingLevels<-function(d,varName,replaceUnknownLevelsWithNA){
    setDF(d)
    d[,varName] = factor(as.character(d[,varName]),exclude = NULL)
    l=levels(d[,varName])
    if(length(l)>0){
      #missingLevels = which(as.numeric(l)<0 )
      #missingLevels = c(missingLevels,which(unlist(lapply(as.list(l),function(x) {!is.numeric(as.numeric(x)) | is.na(as.numeric(x))}))))
      missingLevels = c(unlist(
        lapply(as.list(l),function(x) {x == "." | x =="A" | x =="B" | x =="-9"  | is.na(x)})))
      if(all(missingLevels==FALSE) ) list() else  {
        if(replaceUnknownLevelsWithNA)  rplcmnts = rep(NA,length(missingLevels[missingLevels])) else rplcmnts= paste0("_",l[missingLevels])
        list(l[missingLevels],rplcmnts )
      }
    }else(list())
  }
  
  largestCategory<-function(d,varName){
    setDF(d)
    categories =table(d[,varName])
    names(categories[categories==(max(categories))])
  }
  
  vars = vars[vars %in% names(d)]
  
  if(any(match("RACE",toupper(vars)),na.rm = TRUE)){
    d$RACE = factor(as.character(d$RACE),exclude=NULL)
    ml = addMissingLevels(d,"RACE",replaceUnknownLevelsWithNA)
    levels(d$RACE)<-mapvalues(levels(d$RACE),c(unlist(ml[1]),"1", "2" ,"3" ,"4" ,"5", "6"),
                              c(unlist(ml[2]),'White','Black','Hispanic',
                                'Asian or Pacific Islander', 'Native American','Other'
                              ), warn_missing = warnMissing)
    
    d$RACE=relevel(d$RACE,ref=largestCategory(d,'RACE'))
  }
  
  #urban-Rural
  if(any(match("PL_UR_CAT5",toupper(vars)),na.rm = TRUE)){
    if(length(levels(d$PL_UR_CAT5))>0){
      d$PL_UR_CAT5 = factor(as.character(d$PL_UR_CAT5),exclude=NULL)
      ml = addMissingLevels(d,"PL_UR_CAT5",replaceUnknownLevelsWithNA)
      levels(d$PL_UR_CAT5)<-mapvalues(levels(d$PL_UR_CAT5),c(unlist(ml[1]),1,2,3,4,5),c(unlist(ml[2]),'Large Metro','Small Metro','Micropolitan'
                                                                                        ,'Non-Urban Adjacent Metro', 'Non-Urban Not Adjacent'),
                                      warn_missing = warnMissing)
      d$PL_UR_CAT5=relevel(d$PL_UR_CAT5,ref=largestCategory(d,'PL_UR_CAT5'))
    }
  }
  
  if(any(match("PAY1",toupper(vars)),na.rm = TRUE)){
    d$PAY1 = factor(as.character(d$PAY1),exclude=NULL)
    ml = addMissingLevels(d,"PAY1",replaceUnknownLevelsWithNA)
    
    levels(d$PAY1)<-mapvalues(levels(d$PAY1),c(unlist(ml[1]),1,2,3,4,5,6),c(unlist(ml[2]),'Medicare','Medicaid','Private Ins'
                                                                            ,'Self-pay', 'No charge','Other'),
                              warn_missing = warnMissing)
    d$PAY1=relevel(d$PAY1,ref=largestCategory(d,'PAY1'))
  }
  
  if(any(match("DISPUNIFORM",toupper(vars)),na.rm = TRUE)){
    d$DISPUNIFORM = factor(as.character(d$DISPUNIFORM),exclude=NULL)
    ml = addMissingLevels(d,"DISPUNIFORM",replaceUnknownLevelsWithNA)
    levels(d$DISPUNIFORM)<-mapvalues(levels(d$DISPUNIFORM),c(unlist(ml[1]),1,2,5,6,7,20,99),
                                     c(unlist(ml[2]),'Routine','Short-term Hospital',
                                       'Transfer Other','Home Health Care', 'Against Medical Advice',
                                       'Died', 'destination unknown'),
                                     warn_missing = warnMissing)
    d$DISPUNIFORM=relevel(d$DISPUNIFORM,ref=largestCategory(d,'DISPUNIFORM'))
  }
  d
}

MyTableOne <- function( .descData,covarTableOne,headNote, footNote,
                        .strata,UseDetailedCategories = TRUE,
                        strataNames=c("G One","G Two"),nonnormal =c('TOTCHG','lengthOfObsrv'),
                        SMDQuantile = .7)
{
  setDT(.descData)
  a=unlist(sapply(.descData,function(x) {if(!is.numeric(x) ) length(unique(x))},simplify = TRUE))
  if(length(a[a>10])>1) {
    warning(paste0("Multiple nominal variables with more than 10 unique values:", paste0(names(a[a>10]),collapse=", ") ))
    stop()
  }
  
  
  incs = apply(.descData[,covarTableOne,with=FALSE],2,function(x) length(unique(x)))>=2 #Removing confounders who have less than 2 levels and cause TOne to throw an error
  covarTableOne=covarTableOne[incs]
  if(!any(.strata %in% covarTableOne)) {
    warning("Exposure just has one level or .strata is anot among the given covariates, Table One is skipped")
    return(NULL)
  }
  
  toneData=.descData[,covarTableOne,with=FALSE];setDF(toneData)
  
  # Change variable names ~~~~~~~~~~~~~~~
  # require(plyr)
  # names(toneData)= mapvalues(names(toneData), from=  c("DISPUNIFORM","PAY1","PL_UR_CAT5"),
  #                            to= c('Discharged_To','Insurance', 'Urban_Rural'))
  # Add factor labels ~~~~~~~~~~~~~~~
  
  addLevels <- function(x,lvls, labels){
    if(is.factor(x)) return(factor(x, levels=c(levels(x), lvls[which(!lvls %in% levels(x))])))
    return(x)
  }
  
  if(UseDetailedCategories){
    
    d=translateVars(toneData,covarTableOne[which(covarTableOne %in% c('RACE',"PL_UR_CAT5","PAY1","DISPUNIFORM"))])
    
  }else #MAPPING FACTORS
  {
    #toneData=unique(head(p,10000))
    require(gdata)
    setDT(toneData)
    if(match("ethnicity",toupper(covarTableOne))){
      #unique(toneData$RACE)
      dic=data.frame(oldcode=c('1','2','3','4','5'),
                     newcode=c('2','1','1','1','1'))
      toneData$ethnicity =dic[match(toneData$ethnicity ,dic$oldcode),2]
      toneData$ethnicity  = factor(as.character(toneData$ethnicity ),c(1,2),
                                   labels = c('White','Other'))
      mapLevels(toneData$ethnicity );#table(toneData$RACE,useNA = "ifany")
    }
    
    #urban-Rural
    if(match("PL_UR_CAT5",toupper(covarTableOne))){
      #unique(toneData$PL_UR_CAT5)
      toneData$PL_UR_CAT5 = factor(as.character(toneData$PL_UR_CAT5),c(1,2,3,4,5),
                                   labels = c('Large Metropolitan','Small Metropolitan',
                                              'Miccropolitan','Non-Urban Adjacent Metro', 'Non-Urban Not Adjacent'))
      dic=data.frame(oldcode=c(1,2,3,4,5),
                     newcode=c('Large Metropolitan','Other','Other','Other','Other'))
      toneData$PL_UR_CAT5=dic[match(as.numeric(toneData$PL_UR_CAT5),dic$oldcode),2]
      mapLevels(toneData$PL_UR_CAT5)
    }
    
    if(match("PAY1",toupper(covarTableOne))){
      #unique(toneData$PAY1)
      toneData$PAY1 = factor(as.character(toneData$PAY1),c(1,2,3,4,5,6),
                             labels = c('Medicare','Medicaid','Private insurance','Self-pay', 'No charge','Other'))
      dic=data.frame(oldcode=c(1,2,3,4,5,6),
                     newcode=c('Medicare','Other','Private ins.','Other','Other','Other'))
      toneData$PAY1=dic[match(as.numeric(toneData$PAY1),dic$oldcode),2]
      
      mapLevels(toneData$PAY1)
    }
    
    if(match("DISPUNIFORM",toupper(covarTableOne))){
      #unique(toneData$DISPUNIFORM)
      
      #table(toneData$DISPUNIFORM)
      dic=data.frame(oldcode=c('1','2','5','6','7','20','99'), newcode=c(1,1,1,1,1,2,1)) #'Died','Other')
      toneData$DISPUNIFORM=dic[match(toneData$DISPUNIFORM,dic$oldcode),2]
      toneData$DISPUNIFORM = factor(as.character(toneData$DISPUNIFORM),
                                    c(2,1),
                                    labels = c('Died','Other'))
      #table(toneData$DISPUNIFORM);levels(toneData$DISPUNIFORM1)
      mapLevels(toneData$DISPUNIFORM)
    }
  }
  
  # adjustColTypes<-function(x) {
  #   if( length(unique(x))<10 & is.numeric(x)) x=as.character(x) else x
  #   if(typeof(x)=="logical"  ) x=ifelse(is.na(x),FALSE,x) else x
  #   if(typeof(x)=="character"  ) x=ifelse(is.na(x),"",x) else x
  # }
  #
  # setDT(d);toneData = d[,lapply(.SD,adjustColTypes),.SDcols = names(d) ]
  
  #Adds a prefix to names to distinguish variable names from their sub categories in table1
  nonnormal <-paste0(".",str_replace_all(nonnormal,"cnf_",""))
  names(toneData)<-paste0(".",str_replace_all(names(toneData),"cnf_",""))
  .strata = paste0(".",.strata)
  
  # expPath = paste0(.studyName,"/statTables/")
  # if(!dir.exists(expPath)) dir.create(expPath)
  
  numericCols = names(toneData)[unlist(lapply(toneData,is.numeric))]
  nominalCols = names(toneData)[unlist(lapply(toneData,function(x) !is.numeric(x)))]
  # cldata1=changeColTypes(cldata1,"tonumeric", c('spo2_max.lst','o2flow.lst','rrt.fst','peep.lst'))
  # table(toneData$.peep.lst,useNA = "ifany")
  
  #CreateTableOne Tutorial
  #browseURL("http://rstudio-pubs-static.s3.amazonaws.com/13321_da314633db924dc78986a850813a50d5.html")
  #browseURL("https://cran.r-project.org/web/packages/tableone/vignettes/introduction.html")
  if(debugging) browser()
  require(tableone)
  #tOne = tableone::CreateTableOne(covarDescriptive,strata = "exposed",data=.descData)
  tOne = tableone::CreateTableOne(colnames(toneData),strata = .strata,data=toneData,
                                  factorVars = nominalCols, includeNA=FALSE, test=TRUE,smd = TRUE)
  
  ## If you want to center-align values in Word, use noSpaces option.
  ## quote = TRUE to copy paste to excel
  tOne1<-print(tOne,nonnormal=nonnormal,   quote = FALSE, noSpaces = TRUE, printToggle = FALSE,smd = TRUE)
  if(nrow(tOne1) == nrow(toneData)) stop("You are passing an ID filed to CreateTableOne")
  # Adjust row names ~~~~~~~~~~~~~~~
  require(flextable);require(officer)
  rnames= row.names(tOne1)
  rnames1=vector(mode = "character", length = length(rnames))
  curName=""
  
  for(i in 1:length(rnames)) {
    if (str_detect(rnames[i],"^[.]") )
      curName=rnames[i] #else curName=""
    rnames1[i]=curName
    #if(trimws(rnames[i],which="both")=="")  rnames1[i] = "" else rnames1[i]=curName
  }
  spn = data.frame(spnName =rnames1)%>% dplyr::group_by(spnName) %>% dplyr::summarise(spn=n())
  newNames = dplyr::left_join(data.frame(spnName =rnames),
                              spn,by =c("spnName"="spnName") ) %>%
    dplyr::select(spn,spnName) %>% dplyr::rename(Name = spnName) %>%
    dplyr::mutate(Name=str_replace_all(Name,"^[.]","")) %>%
    dplyr::mutate(Name=str_replace_all(Name," = TRUE","")) %>%
    dplyr::mutate(Name=str_replace_all(Name," = 1","")) %>%
    dplyr::mutate(Name=str_replace_all(Name,"DISPUNIFORM",'Discharged_To')) %>%
    dplyr::mutate(Name=str_replace_all(Name,"PAY1",'Insurance')) %>%
    dplyr::mutate(Name=str_replace_all(Name,"PL_UR_CAT5",'Urban_Rural'))
  
  
  tOne2 = cbind(grps = as.character(rnames1),Variable = newNames$Name,tOne1,vspan= newNames$spn)
  
  tOne2= as.data.frame(tOne2)
  tOne2=setDT(tOne2)[Variable=="n",`:=`(grps="n",SMD="inc")]
  grporder = tOne2 %>% filter(!is.na(SMD) & !is.nan(SMD) & SMD !="") %>% arrange(desc(SMD)) %>%
    dplyr::mutate(grpOrder =1:n())  %>% dplyr::select(grps,grpOrder)
  
  smdThreshold = quantile(na.omit( as.double(as.character( tOne2$SMD))),SMDQuantile)
  grpsToRemove = tOne2 %>%  dplyr::filter(as.double(as.character(SMD)) <smdThreshold & !(SMD =="NaN"| SMD ==""| SMD =="inc"))
  
  tOne2 = tOne2 %>% dplyr::filter(!grps %in% grpsToRemove$grps) %>% left_join(grporder, by= "grps") %>%  dplyr::group_by(grps) %>%
    dplyr::arrange(grpOrder) %>%   ungroup() %>% dplyr::select(-c(grps,grpOrder,test))
  #%>%  dplyr::filter(as.double(as.character(SMD)) >=smdThreshold | SMD =="NaN"| SMD ==""| SMD =="inc") %>%
  
  setDT(toneData)
  strataNames=c()
  strataNames = na.omit(as.vector( unique(toneData[,get(.strata)])))
  strataNames = paste0(.strata," = " ,strataNames)
  
  if(.strata == ".exposed") strataNames = c("Unexposed","Exposed")
  if(.strata == ".event") strataNames = c("No Event","Event Pos.")
  if(.strata == ".Matched") strataNames = c("Unmatched","Matched")
  #if(length(.strata) >0)  colnames(tOne2)<-c("Variable",strataNames,"P Value","Test", "SMD","VSpan")
  
  #tOne2 = tOne2[,lapply(.SD, as.character),.SDcols = names(tOne2)]
  # More detailed view of categorical variables can be obtained via summary.CatTable methods
  #summary(tableOne$CatTable)
  # See the continuous part (ContTable) only using $ operator
  #tableOne$ContTable
  
  require(huxtable); #browseURL("https://hughjonesd.github.io/huxtable/huxtable.html")
  hx2 = huxtable(tOne2)  #(table1)
  # for(i in 1:nrow(hx2)) rowspan(hx2)[i,2]<-
  #   ifelse(!is.na(hx2[i,"VSpan"]),as.integer(hx2[i,"VSpan"]),1)
  
  colnames(hx2)<-c("Variable",strataNames,"P Value", "SMD","VSpan")
  
  hx2 =hx2[,-7] %>% set_width(1) %>%
    set_number_format(list(function(x) prettyNum(x, big.mark = ','))) %>%
    set_col_width(c(1.5,1.2,1.2,1,1,.5)) %>%
    set_align(everywhere,everywhere,"center")%>%
    set_align(everywhere,1,"right") %>%
    huxtable::set_caption(value=headNote) %>%
    theme_article( header_row = FALSE, header_col = FALSE)%>%
    add_footnote(footNote)
  
  hx2= map_bold(ungroup(hx2),by_rows(c(TRUE,rep(FALSE,nrow(hx2)-1))))#setting bold the caption row
  for(i in 1:nrow(tOne2)){
    if(!is.na(tOne2[i,"vspan"])){
      bottom_border(hx2)[i,] <-1
      hx2 = hx2 %>% set_bold(i,1:ncol(hx2),value = TRUE)
    }
  }
  
  fx2=huxtable::as_flextable(hx2,colnames_to_header = TRUE);#inherits(fx, "FlexTable")
  fx2
}

# Load data extract ground truth -----


rcrd = list()


write(paste0('===== ====== CICT analysis, Abbas Shojaee, Huang Lab, Now: ', lubridate::now(), '======= '),  file = url.logfile, append = TRUE)
write('Step 0',  file = url.logfile, append = TRUE)
# CICT features calculations -------

n.itm.e$edgeTyp = '' #'ewMIempirical' #tlag1.mi  'ewMImm' #"tlag1.mi" #ewMIempirical #Pearson #'corP'
n.itm.e=n.itm.e %>% dplyr::mutate(Weight=.data[[cictRawEdgeCol]]  ) %>% 
  select(src,trgt,everything()) # mf.mi )#)tlag1.mi  #tlag1.corP

originalEdgeMeasures = setdiff(colnames(n.itm.e), c('src','trgt'))

n.itm.v = data.table()
n.itm.v=setDF(n.itm.e) %>%  
  group_by(src) %>% summarize(ocr = sum(.data[[cictRawEdgeCol]])) %>%
  rename(gene =src) %>% mutate(vID=gene) %>% ungroup()


# n.itm.v=setDF(all.dt) %>%  
#   dplyr::mutate(ocr = rowSums(.[which(!colnames(setDF(all.dt)) %in% c('vID','gene'))])) %>% 
#   dplyr::select(gene,ocr) %>% rename(vID=gene)

n.itm.v=setDT(n.itm.v)[,`:=`(OcrInp=ocr,OcrOut=ocr)][,ocr:=NULL]

msg = paste0(sprintf("Raw edges= %s, vertices= %s",  nrow(n.itm.e),nrow(n.itm.v)),
             paste0(" Early pruning threshold= ",earlyThresholdForGraphAnalysis),
             sprintf(" Remained edges= %s \n",  nrow(n.itm.e[n.itm.e$Weight>earlyThresholdForGraphAnalysis,])))
cat(msg)
write(msg, file = url.logfile, append = TRUE)

rm(ig)
ig=graph.data.frame(n.itm.e[n.itm.e$Weight>earlyThresholdForGraphAnalysis,],directed=FALSE) 
#bibliometrix::networkStat(ig)
ind = igraph::degree(ig,v=V(ig),mode = c("in")); ind=data.frame(Indegree=ind,subcat=names(ind))
outd = igraph::degree(ig,v=V(ig),mode = c("out")); outd=data.frame(Outdegree=outd,subcat=names(outd))
n.itm.v=n.itm.v %>% left_join(ind, by=c( "vID" ="subcat"))
n.itm.v=n.itm.v %>% left_join(outd, by=c( "vID" ="subcat"))

vertexOcrSumIn = sum(n.itm.v$OcrInp)
vertexOcrSumOut = sum(n.itm.v$OcrOut)
setDT(n.itm.v);n.itm.v[,c('probInp','probOut'):=list(OcrInp/vertexOcrSumIn,OcrOut/vertexOcrSumOut) ] #probability
setDF(n.itm.v)  ;setDF(n.itm.e)
vinfcols=c("vID","OcrInp","OcrOut","probInp","probOut","Indegree","Outdegree")

n.itm.e = n.itm.e %>% left_join(n.itm.v[,vinfcols],by=c("src"="vID"))
n.itm.e = n.itm.e %>% left_join(n.itm.v[,vinfcols],by=c("trgt"="vID"))
n.itm.v.trgtsum = n.itm.v%>% left_join(n.itm.e,by=c("vID"="src")) %>%
  group_by(vID) %>% summarise(SumOcrInp.y=sum(OcrInp.y,na.rm=TRUE),
                              SumOcrOut.y=sum(OcrOut.y,na.rm=TRUE))
# sumTrgt == SumOcrInp.y, sumSrc=SumOcrOut.x , OcrOut.x=OcrOut.x


n.itm.v.srcsum = n.itm.v%>% left_join(n.itm.e,by=c("vID"="trgt")) %>%
  group_by(vID) %>%summarise(SumOcrInp.x=sum(OcrInp.x),na.rm=TRUE,
                             SumOcrOut.x=sum(OcrOut.x),na.rm=TRUE)

n.itm.v = n.itm.v %>% left_join(n.itm.v.trgtsum,by=c("vID"="vID"),copy=TRUE) %>%
  inner_join(n.itm.v.srcsum,by=c("vID"="vID"),copy=TRUE) #%>%

cols= c("SumOcrInp.x","SumOcrOut.x","SumOcrInp.y","SumOcrOut.y",
        "OcrOut","OcrInp") #"OcrOut.x","OcrInp.x","OcrOut.y","OcrInp.y"

factorToNumeric <- function(x) as.numeric(as.character(x))

n.itm.v.num = setDT(n.itm.v)[, lapply(.SD, factorToNumeric), .SDcols = cols];
n.itm.v = n.itm.v[, !names(n.itm.v) %in% cols, with = FALSE];gc()
n.itm.v = cbind(n.itm.v.num, n.itm.v);rm(n.itm.v.num);gc()

setDT(n.itm.v);n.itm.v[is.na(SumOcrOut.x),SumOcrOut.x:=1];n.itm.v[is.na(SumOcrInp.y),SumOcrInp.y:=1]


vinfcols=c("vID","SumOcrInp.x","SumOcrOut.x","SumOcrInp.y","SumOcrOut.y")
setDF(n.itm.e);setDF(n.itm.v)
n.itm.e = n.itm.e %>% inner_join(n.itm.v[,vinfcols],by=c("src"="vID"))
n.itm.e = n.itm.e %>% inner_join(n.itm.v[,vinfcols],by=c("trgt"="vID"))

#a=sapply(n.itm.e, function(x) sum(is.na(x)));a[a>0]

n.itm.e = replace_na(n.itm.e,replace=list(SumOcrOut.x.x=0,SumOcrOut.x.y=0,
                                          SumOcrInp.x.x=0,SumOcrInp.x.y=0,
                                          SumOcrInp.y.x=0,SumOcrInp.y.y=0,
                                          SumOcrOut.y.x=0,SumOcrOut.y.y=0
))

emissionRate = 1
maxItmstOcr = max(c(n.itm.v$OcrInp, n.itm.v$OcrOut))


n.itm.e$OcrOut.x = as.numeric(n.itm.e$OcrOut.x )
n.itm.e$OcrInp.y = as.numeric(n.itm.e$OcrInp.y )

n.itm.e = setDF(n.itm.e) %>% dplyr::mutate(conf=Weight/ OcrOut.x , contrib=Weight/OcrInp.y) 

n.itm.e.back0 = n.itm.e;n.itm.v.back0 = n.itm.v #n.itm.v=n.itm.v.back0;n.itm.e=n.itm.e.back0
rm(n.itm.v.srcsum,n.itm.v.trgtsum,n.itm.v.influx,n.itm.v.outflux,n.itm.e.tmp)

#   MIX the Data + Ground truth  ----
rm(tst1.totalset,t2,tmp,tst1.train,tst1.tst,tst1.tst.h2o,tst1.totalset.h2o,tst1.train.h2o)
rm(n.itm.e.back3,L2,MisLinks,MisNodes,self,rds1,rds2,rds3,outd,others,t2.tmp1,b,e,e.exmpl.c,e.exmpl.e,e.exmpl.r,n.notinVertices,v,t1.c,t1.rc,t1.causal_reversecausal,rpt,t1.rnd,t1,total.rds,tmp.rdndnt)
rm(tst1.train,tst1.totalset.h2o,tst1.totalset,tst1.train.h2o,tst1.tst,tst1.tst.h2o,tst1.mdl,tst1.rf,tmp.inps,tmp.outps,tmp1)
gc()
{
  #NrandomEdges =  min(.3 * nrow(n.itm.e),200000) ; #NrandomEdges=2e6
  
  
  #or 0 means all
  #n.itm.e = n.itm.e %>% filter(src %in% geneorder.2e3$X | trgt %in% geneorder.2e3$X)
  
  
  setDF(n.itm.e)
  #Test
  # tmp =n.itm.e.back0 %>% inner_join(t1,by=c("first"="src","second"="trgt"))
  # t1.abscentGT =t1 %>% anti_join(n.itm.e.back0,by=c("src"="first","trgt"="second"))
  # t1.abscentGT =t1.abscentGT %>% anti_join(n.itm.e.back0,by=c("src"="second","trgt"="first"))
  
  #n.itm.e= n.itm.e %>% dplyr::mutate(src=first,trgt=second)
  
  n.itm.e.back3 = n.itm.e 
  rm(t1,t1.c,t1.rc,t1.c.r,t1.abscentGT,t2.rc,t2.rnd)
  
  #intersect(n.itm.e$src,tbl.goldStandard$src) #'CXXC1') 
  
  t1.c = tbl.goldStandard %>% dplyr::select(src,trgt) %>% 
    inner_join(n.itm.e,by=c("src"="src","trgt"="trgt"))
  if(nrow(t1.c) < nrow(tbl.goldStandard)/2) warning("Problem in ground truth. More than half of ground truth was not found in edges")
  # t1.c.r = n.itm.e %>% inner_join(t1,by=c("src"="trgt","trgt"="src"))
  # t1.c = rbind(t1.c,t1.c.r)
  if(nrow(t1.c) <=0) print("!!! No causal edge in the groundturht? check gene names")
  
  t1.c$predicate = "CAUSES"
  
  t1.rc = n.itm.e %>% inner_join(t1.c %>% select(src,trgt),by=c("src"="trgt","trgt"="src"))
  t1.rc$predicate = "REV_CAUSES"
  
  t1.causal_reversecausal=rbind(t1.c,t1.rc)
  
  #Adding 2000 random edges
  t1.rnd = n.itm.e %>% #dplyr::mutate(edgetyptruth = NA) %>%
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
  NrandomEdges = NcausalEdges *  randomEdgesFoldCausal #min(sampling.rnd.ratio * nrow(n.itm.e),sampling.u.max) %>% as.integer() ; #NrandomEdges=2e6
  
  t2=rbind(t1 %>% filter(class1=='c') %>% sample_n(size = NcausalEdges),
           t1 %>% filter(class1=='rc') %>% sample_n(size = NcausalEdges),
           t1 %>% filter(class1=='u')%>% sample_n(size=NrandomEdges-2*NcausalEdges)) #to preserve the minGroundTruth.ratio
  
  t2.complement = t1  %>%  anti_join(t2,by=c("src"="src","trgt"="trgt"))
  
  msg=sprintf(' Network density= %s, learning set density= %s, \n total edges= %s, learning set=%s, test fraction = %s , \n learning set random= %s, learning set causal edges= %s',  
              round(nrow(t1.c)/nrow(n.itm.e),4),round(sum(t2$predicate=='CAUSES')/nrow(t2),4),
              prettyNum(nrow(n.itm.e),','),prettyNum(nrow(t2),','),tstPrecent,
              round(NrandomEdges,2),round(NcausalEdges,2)
  ) #tstPrecent
  
  cat(msg)
  write(msg, file = url.logfile, append = TRUE)
  #t2 = rbind(t1.c,t1.rnd)
  #t2$edgeTyp = ''
  t2= as.data.frame(t2) %>% dplyr::select(any_of(c('edgeTyp' ,'edgetyptruth','src', 'trgt','Weight',
                                                   'OcrOut.x', 'OcrInp.y')),everything())
  table(t2$predicate)#,t2$edgeTyp, useNA="ifany")
  
  
  write(table(t2$predicate) %>% knitr::kable(),  file = url.logfile, append = TRUE)
  
  
  a=sapply(t2, function(x) sum(is.na(x)));a[a>0]
  
  library(stringr)
  
  
}

#Early pruning
if(F){
  
  
  n.itm.e.back3 = n.itm.e  #n.itm.e = n.itm.e.back3
  #n.itm.e = n.itm.e %>% filter(src %in% geneorder.2e3$X | trgt %in% geneorder.2e3$X)
  
  rm(t1,t1.c,t1.rc,t1.c.r,t1.abscentGT,t2.rc,t2.rnd,e,tst1.totalset,tst1.tst,tst1.train,t1,t2)
  
  NrandomEdges =  min(.9 * nrow(n.itm.e),500000) ; 
  t1.gtfrac = 1
  
  # t1.gt=sc.gtchip %>% dplyr::select(src,trgt,edgetyptruth)
  t1.c = n.itm.e %>% inner_join(t1.gt,by=c("src"="src","trgt"="trgt"))%>%
    mutate(predicate = "CAUSES")
  
  t1.rc = n.itm.e %>% inner_join(t1.gt,by=c("src"="trgt","trgt"="src")) %>%
    mutate(predicate = "REV_CAUSES")
  
  t1.causal_reversecausal=rbind(t1.c,t1.rc)
  t1.causal_reversecausal.sub = if(t1.gtfrac==1) t1.causal_reversecausal else sample_frac(t1.causal_reversecausal, size = t1.gtfrac)  
  
  t1.rnd = n.itm.e %>% dplyr::mutate(edgetyptruth = NA,predicate = "IRRELEVANTM") %>%
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
# tmp =n.itm.e.back0 %>% inner_join(t1,by=c("first"="src","second"="trgt"))
# t1.abscentGT =t1 %>% anti_join(n.itm.e.back0,by=c("src"="first","trgt"="second"))
# t1.abscentGT =t1.abscentGT %>% anti_join(n.itm.e.back0,by=c("src"="second","trgt"="first"))

maxGroundTruth = 6000; t2.frac = .2 #or 0 means all

t2= as.data.frame(t2) %>% dplyr::select(any_of(c('edgeTyp' ,'edgetyptruth','src', 'trgt','Weight',
                                                 'OcrOut.x', 'OcrInp.y',#intvl_avg,intvl_sd,intvl_median,
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

#print(colList[!colList %in% colnames(n.itm.v)])
#mdlColNames = colnames(t2)
mdlColNames=c(additional_provided_cols, 'corP','corS','corK','mf.mi','tlag1.corP','tlag1.corS','tlag1.corK','tlag1.mi',
              "efMImm","ewMImm","efMIempirical","ewMIempirical", "efMIshrink","ewMIshrink",
              'Indegree.x','Indegree.y','Outdegree.x','Outdegree.y',
              'intvl_median','intvl_kurtosis','intvl_avg','intvl_skew','intvl_sd',
              'ZOcfNMedian.x','ZOcfNMedian.y','ZScbNSkew.x','ZScbNSkew.y',
              'ZOcfSkew.x','ZPoutsum.x','ZOcfSkew.y','ZPoutsum.y','ZScbScbL3.x','ZScbScbL3.y',
              'HTR','tNHTR','NHTRfromSource','NHTRtoTarget',
              'NHTRfromSourceBA','NHTRtoTargetBA','HTRBA','tNHTRBA',
              'NHTRDiff_fromSource','NHTRRatio_fromSource',
              'NHTRDiff_toTarget','NHTRRatio_toTarget',
              'EEBA','NEEfromSourceBA','NEEtoTargetBA',
              'EE','NEEfromSource','NEEtoTarget','NEE_NHTR_Diff','NEE_NHTR_Ratio',
              'Weight','R','V','P','UR','URR','URdiff','Tendency',
              'EI','EO','OEER','OERR','AM','HM','GM',
              'URBA','AMBA','IBA','OERRBA','OEERBA','tBA','tnBA','toeBA',
              'confDiff','confNDiff','contribNDiff','contribDiff',
              'PDiff','PAdd', 'PdirDiff','PdirRatio',
              'directionLogRatio',
              'hypotenuse','slope',
              'URdiff' ,'AMdiff','deltaAttitude', 'Idiff', 'IRatio',
              't','tn','toe','tdiff','tndiff','toediff','tRatio', 'tnRatio',
              'toeRatio','tz','tnz',
              'conf','contrib','confBA','contribBA' , #'vertexOcrSum',
              'OcrOut.x','OcrInp.y', 'OcrOut.x2','trgtOcr2',
              'srctrgtSum','srctrgtProduct',
              'confdisc', 'contribdisc', 'cnfxcnfy.emi','cnfxcnty.emi','cntxcnfy.emi','cntxcnty.emi','cnfxcnfy.frcht','cnfxcnty.frcht','cntxcnfy.frcht','cntxcnty.frcht',
              'cnfxcnfy.emiBA','cnfxcnty.emiBA','cntxcnfy.emiBA','cntxcnty.emiBA','cnfxcnfy.frchtBA','cnfxcnty.frchtBA','cntxcnfy.frchtBA','cntxcnty.frchtBA',
              'efMImm','ewMImm','efMIempirical','ewMIempirical','efMIshrink','ewMIshrink','Pearson','Spearman','Kendall','Manhattan','Euclidean','L10Norm',
              'conf.fitdist.x','conf.beta1STshape.x','conf.beta2ndshape.x','conf.loglinmean.x','conf.loglinsd.x','cntrb.fitdist.x','cntrb.beta1STshape.x','cntrb.beta2ndshape.x','cntrb.loglinmean.x','cntrb.loglinsd.x',
              'conf.fitdist.y','conf.beta1STshape.y','conf.beta2ndshape.y','conf.loglinmean.y','conf.loglinsd.y','cntrb.fitdist.y','cntrb.beta1STshape.y','cntrb.beta2ndshape.y','cntrb.loglinmean.y','cntrb.loglinsd.y',
              
              'odds','confZ','contribZ','causal1','causal11','causal2','causal3',
              'confN','contribN','confNBA','contribNBA','oddsN',
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
         'NHTRfromSourceBA','NHTRtoTargetBA','HTRBA','tNHTRBA',
         'NHTRDiff_fromSource','NHTRRatio_fromSource',
         'NHTRDiff_toTarget','NHTRRatio_toTarget', 'ocfSkew.x', 'Poutsum.x'  )


#Remove large objects
#sort( sapply(ls(),function(x){object.size(get(x))}),decreasing = F)
rm(dfl,df,d.new,mtx.lng.prtrb,tzcolvec,tst1.totalset);gc()
# tst1 ,tst2 Building EXPERIMENT data:  - training and test sets --------  
require(caret)
setDT(t2)
#existingCols = intersect(colnames(t2),c(mdlChkColNames,objectiveColNames))
tst1.totalset = t2[,intersect(colnames(t2),unique(c(mdlChkColNames,objectiveColNames))),with=FALSE]
tst1.totalset = tst1.totalset %>% select('src','trgt', originalEdgeMeasures,matches('class'))
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
while(TRUE){
  set.seed(runif(1,1,10000))
  spltIdx =as.vector( caret::createDataPartition(1:nrow(tst1.totalset),p=(1-tstPrecent),list=FALSE,times=1))
  tst1.train = tst1.totalset[spltIdx,] #tst1.dmy[spltIdx,]
  tst1.tst = tst1.totalset[-spltIdx,]
  # 
  # spltIdx =as.vector( caret::createDataPartition(1:nrow(n.itm.v),p=(1-tstPrecent),list=FALSE,times=1))
  # trainV = paste0('G',spltIdx)
  # tst1.train = setDT(tst1.totalset)[src %in% trainV | trgt %in% trainV,] #tst1.dmy[spltIdx,]
  # tst1.tst = tst1.totalset[!(src %in% trainV | trgt %in% trainV),]
  
  print(nrow(tst1.tst[get(trainingTarget)== TRUE,]))
  if(nrow(tst1.tst[get(trainingTarget)== TRUE,]) >= (tstPrecent-0.01)* ntrgtClass) break
}

write(paste0('tst1.train$class2 : ',table(tst1.train$class2)),  file = url.logfile, append = TRUE)
write(paste0('tst1.tst$class2 : ',table(tst1.tst$class2)),  file = url.logfile, append = TRUE)

table(tst1.train$class2);table(tst1.tst$class2)
table(tst1.train$class3);table(tst1.tst$class3)
table(tst1.train$class5);table(tst1.tst$class5)


#save(n.itm.e,n.itm.v,t2,tst1.totalset,mdlChkColNames,mdlColNames,file="CICT dream4-100-1-1 MF two way.robj")
#load(file="CICT dream4-100-2.robj")
# Performance and output functions -----
#Testing performance
testPerformance<-function(modelName,tst1.mdl,url.logfile,url.output){
  url.output= paste0(url.output,'/',modelName)
  if(!dir.exists(url.output)) dir.create(url.output,recursive = T)
  
  prd.varimp=h2o.varimp(tst1.mdl)#;View(prd.varimp)
  write( prd.varimp[1:20,]%>% knitr::kable(),  file = url.logfile, append = TRUE, sep = '\n')
  #plot(tst1.mdl);setDF(tst1.totalset)
  #sink(NULL)
  #sink(file = paste0(url.logfile,'_t.txt'), append = T, type = c("output"),split = T) #write(  file = url.logfile, append = TRUE)
  h2o.performance(tst1.mdl,valid=T) #save(tst1,tst1.train,tst1.tst, file="temp training set.robj")
  
  msg = c('================================================',
          modelName,
          '================================================',
          "Reporting model performance on validation set",
          capture.output(h2o.performance(tst1.mdl,valid=T)))
  cat(paste0(msg,collapse='\n') )
  write( msg,  file = url.logfile, append = TRUE, sep = '\n')
  
  
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
      mutate(aucs = round(aucs,3), standardized = round(standardized,3) ) %>%
      knitr::kable()
  }
  #newDataEnv = new.env()
  #R.utils::loadToEnv(file="CICT dream4-100-2.robj",envir=newDataEnv)
  
  #Predictalldata
  if(FALSE){
    setkey(t2,src,trgt)
    setkey(n.itm.e,src,trgt)
    d.new = n.itm.e[!t2]
    d.new[,(trainingTarget) := ifelse(is.na(get(trainingTarget)),0,get(trainingTarget))]
  }
  
  #Predicts in smaller chunks also adds random prediction for comparition
  msg = c('================================================',
          "Reporting results on unseen sample")
  cat(paste0(msg,collapse='\n') )
  write( msg,  file = url.logfile, append = TRUE, sep = '\n')
  
  d.new = t2.complement %>% sample_n(size = min(25000,nrow(t2.complement) * maxunseenTest.ratio))   #tst1.totalset #newDataEnv$tst1.totalset # 
  table(d.new$predicate)
  
  msg = sprintf("Learning set= %s | training= %s | validation= %s | unseen sample = %s | t2 + comp = %s | all Edges= %s",  
                nrow(tst1.totalset),nrow(tst1.train),nrow(tst1.tst),nrow(d.new), nrow(t2)+nrow(t2.complement),nrow(n.itm.e))
  
  cat(msg)
  write(msg, file = url.logfile, append = TRUE)
  
  #d.new = n.itm.e
  splitcount = 2
  splts = split(1:nrow(d.new),          
                cut(seq_along(1:nrow(d.new)),
                    splitcount,labels = FALSE))
  
  
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
                       rndPred = ifelse(runif(nrow(d.new.slice)) >=.5,1,0) #Assign a random classifier results
                       prd_outcomes = d.new.slice %>% select(src,trgt,any_of('Weight')) %>%
                         cbind(predictions,outcomes,rndPred) %>% as.data.frame()
                       prd_outcomes
                     })
  d.new1= rbindlist(d.new.tmp)
  
  
  d.new1.rv =d.new1 %>% dplyr::rename(revWeight = Weight,rvpred=predictions,src1=src,trgt1=trgt) %>% 
    dplyr::select(-outcomes,-rndPred)
  pred_outcome = d.new1 %>% left_join(d.new1.rv, by=c("src"="trgt1", "trgt"="src1") )
  pred_outcome.back = pred_outcome
  #pred_outcome$predictions %>% table() 
  
  
  #TODO remove the reverse edges before assessment, just keep the one with higher prediction score
  assespreds = pred_outcome #%>% dplyr::filter(predictions>=rvpred)
  theROC <- roc(assespreds$outcomes, assespreds$predictions, percent = TRUE);theROC
  
  sscurves <- evalmod(scores = assespreds$predictions, labels = assespreds$outcomes)
  pr.crv=sscurves$prcs[1][[1]];  pr.auc =  attr(pr.crv,'auc')
  
  theauc = precrec::auc(sscurves); 
  msg = paste0(theauc[[3]],"=", round(theauc[[4]],3))
  write(msg,  file = url.logfile, append = TRUE)
  ##autoplot(sscurves, mode = 'rocpr',ret_grob = F)
  
  #partial precision-recall
  sscurves.part <- part(sscurves, xlim = c(0, 0.2))
  ##plot(sscurves.part,title = "Partial top 20% predictions")
  write('Early precision recall ======',  file = url.logfile, append = TRUE)
  write(reportAUC(sscurves.part),  file = url.logfile, append = TRUE)
  
  
  pr.prtcrv=sscurves.part$prcs[1][[1]]
  if(is.null(pr.prtcrv)) {i=0;while(i <10000 & is.null(pr.prtcrv)){i=i+1 ;pr.prtcrv<-sscurves.part$prcs[1][[1]]}}
  if(is.null(pr.prtcrv)) browser()
  pr.prtauc =  attr(pr.prtcrv,'pauc')
  
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
  
  msg = sprintf("AUCPR Ratio to Random= %s,  Partial AUCPR Ratio to Random= %s ", 
                round(pr.auc/rnd.auc,2), round(pr.prtauc/rnd.prtauc,2))
  print(msg)
  write(msg,  file = url.logfile, append = TRUE)
  
  mmpoins <- evalmod(scores = assespreds$predictions, labels = assespreds$outcomes, mode = "basic")
  
  # the relative cost of of a false negative classification (as compared with a false positive classification)
  # the prevalence, or the proportion of cases in the population (n.cases/(n.controls+n.cases)).
  relativeCostfn_fp = 1/2
  prv = table(setDT(tst1.totalset)[,get(trainingTarget)])
  best.weights=c(relativeCostfn_fp, prv[2]/prv[1])  
  bestcutoff =as.double(coords(theROC, "best", best.method = "closest.topleft", best.weights=best.weights,
                               ret="threshold", transpose = FALSE));bestcutoff
  print(sprintf("Best threshold: %s for fn/fp cost ratio: %s " ,bestcutoff,relativeCostfn_fp ))
  write(sprintf("Best threshold: %s for fn/fp cost ratio: %s " ,bestcutoff,relativeCostfn_fp ),  file = url.logfile, append = TRUE)
  
  #bestcutoff = 0.5
  tryCatch({  
    assespreds =assespreds %>% mutate(thresholdpreds= ifelse(assespreds$predictions>bestcutoff,assespreds$predictions,0)) 
    theROC <- roc(assespreds$outcomes, assespreds$thresholdpreds, percent = TRUE);theROC
    
    sscurves <- evalmod(scores = assespreds$thresholdpreds, labels = assespreds$outcomes)
    pr.crv=sscurves$prcs[1][[1]];  pr.auc =  attr(pr.crv,'auc')

    msg = sprintf("With FP-FN ratio 0.5 => AUCPR Ratio to Random= %s,  AUCPR Ratio to Random= %s ", 
                round(pr.auc/rnd.auc,2), round(pr.prtauc/rnd.prtauc,2))
    print(msg)
    write(msg,  file = url.logfile, append = TRUE) 
  },error=function(e) browser())
  write('Step5.7 Model trained successfuly',  file = url.logfile, append = TRUE)
  bestcutoff
  
}

#final run on all edges
produceFinalEdges<-function(modelName,tst1.mdl,url.logfile,url.output,bestcutoff=0){
  url.output= paste0(url.output,'/',modelName)
  if(!dir.exists(url.output)) dir.create(url.output,recursive = T)
  url.rankedEdges = paste0(url.output , '/rankedEdges.csv')
  url.rankedEdgesGated = paste0(url.output , '/rankedEdgesGated.csv')
  
  
  write('****************************',  file = url.logfile, append = TRUE)
  write(modelName,  file = url.logfile, append = TRUE)
  write('Step 6, edge ranking',  file = url.logfile, append = TRUE)
  d.new = n.itm.e # %>% sample_frac(size =.01)  # %>% anti_join(t2,by=c("src"="src", "trgt"="trgt"))  #tst1.totalset #newDataEnv$tst1.totalset # 
  #d.new = n.itm.e
  splitcount = max(floor(nrow(n.itm.e) / 30000),2)
  splts = split(1:nrow(d.new),          
                cut(seq_along(1:nrow(d.new)),
                    splitcount,labels = FALSE))
  
  
  d.new.tmp = lapply(splts,
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
                     })
  d.new1= rbindlist(d.new.tmp)
  
  write('Step 6.1, edges ranking finished',  file = url.logfile, append = TRUE)
  
  d.new1.rv =d.new1 %>% dplyr::rename(revWeight = Weight,rvpred=predictions,src1=src,trgt1=trgt) %>% 
    dplyr::select(-rndPred)
  pred_outcome = d.new1 %>% left_join(d.new1.rv, by=c("src"="trgt1", "trgt"="src1") )
  pred_outcome=pred_outcome %>% rename(Gene1=src,	Gene2=trgt,	EdgeWeight = predictions) %>%
    select(Gene1,Gene2,EdgeWeight,everything()) %>% arrange(desc(EdgeWeight))
  
  try({file.remove(url.rankedEdges)},silent=T)
  fwrite(pred_outcome,url.rankedEdges,row.names = F, sep='\t')
  
  truncated_pred_outcome = pred_outcome %>% 
    mutate(EdgeWeight=ifelse(EdgeWeight>bestcutoff,EdgeWeight,0))
  fwrite(truncated_pred_outcome,url.rankedEdgesGated,row.names = F, sep='\t')
  
  
  write('================================================',  file = url.logfile, append = TRUE)
  write('Data produced successfuly',  file = url.logfile, append = TRUE)
  print('Data produced successfuly ==================================')
  
}


#   TRAINING H2O  For Prediction and feature selection ----
#Remove large objects
write('Step5.6',  file = url.logfile, append = TRUE)
#sort( sapply(ls(),function(x){object.size(get(x))}),decreasing = F)

#require(retry)
require(h2o)
#try({h2o.shutdown(prompt=FALSE)})
gc()

#retry::retry({
H2OCnn = h2o.init(nthreads = parallel::detectCores()-2, enable_assertions = TRUE,
                  max_mem_size = "20g",strict_version_check=FALSE) 
#},when='OverflowError|54321: Connection refused',interval=15,max_tries=3) 

#memory.limit(size=25000)
#Create binary outcome to calculate probaibility
Sys.sleep(5)
#trainingTarget;table(tst1.totalset[,get(trainingTarget)])

NFolds = 5
table(setDT(tst1.totalset)[,get(trainingTarget)])
filterOutFeatures =  c('srctrgtSum','srctrgtProduct') #Weight
#mdlColNames = intersect(mdlColNames, prd.varimp$variable[1:150])
mdlColNames= intersect(mdlColNames,colnames(tst1.totalset))
mdlColNames=setdiff(mdlColNames, filterOutFeatures)
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
                           keep_cross_validation_predictions= FALSE,  validation_frame = tst1.tst.h2o)
#,ntrees=30,max_depth=5#selectedFeatures #mdlColNames
tst1.mdl=tst1.rf
h2o::h2o.varimp(tst1.mdl)

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


  
  rcrd$learning_curve_plot_auto = h2o.learning_curve_plot(tst1.mdl,metric = c(  "Auto"))
  rcrd$learning_curve_plot_aucpr = h2o.learning_curve_plot(tst1.mdl,metric = c(  "aucpr"))
  rcrd$learning_curve_plot_auc = h2o.learning_curve_plot(tst1.mdl,metric = c(  "auc"))
  rcrd$learning_curve_plot_logloss = h2o.learning_curve_plot(tst1.mdl,metric = c(  "logloss"))
  rcrd$learning_curve_plot_amisclassification = h2o.learning_curve_plot(tst1.mdl,metric = c(  "misclassification"))
  rcrd$learning_curve_plot_rmse = h2o.learning_curve_plot(tst1.mdl,metric = c(  "rmse"))
  

bestcutoff=testPerformance( modelName='randomForest',tst1.rf,url.logfile,url.output)
produceFinalEdges( 'randomForest',tst1.rf,url.logfile,url.output,bestcutoff)



# tst1.gbm = h2o.gbm(mdlColNames,trainingTarget,
#                    tst1.train.h2o,nfolds = NFolds,#checkpoint=tst2.mdl@model_id,
#                    max_depth =3,
#                    #ntrees=12,
#                    keep_cross_validation_predictions= FALSE,  validation_frame = tst1.tst.h2o ) #ntrees=20
# tst1.mdl = tst1.gbm

# bestcutoff=testPerformance( 'gbm',tst1.rf,url.logfile,url.output)
# produceFinalEdges( 'gbm',tst1.rf,url.logfile,url.output,bestcutoff)



