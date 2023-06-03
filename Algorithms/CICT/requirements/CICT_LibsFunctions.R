#Libraries -----

#library(mlbench)
library(caret)
#library(Spectrum)
library(Rtsne)
library(umap)
library(devtools)
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
#library(DREAM4)
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
require(yaml)
require(Lmoments)
require(GGally)
library(PerformanceAnalytics)
library(hutils)
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


