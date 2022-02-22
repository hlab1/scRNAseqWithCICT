#Edu
# browseURL('https://www.youtube.com/watch?v=l4BAfRekohk')
# browseURL('https://scrnaseq-course.cog.sanger.ac.uk/website/index.html')
# browseURL('(http://bar.utoronto.ca/eplant') #, a visual analytic tool for exploring multiple levels of Arabidopsis thaliana data ')
# Environment ----- 
if(FALSE){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install("arabidopsis.db0")
  
  install.packages("GenNet") #http://www.strimmerlab.org/software/genenet/
  install.packages("myTAI")
  install.packages("DREAM4")
  devtools::install_github("paubellot/netbenchmark")
  install.packages("Seurat") #browseURL("https://satijalab.org/seurat/vignettes.html")
  install.packages("RCytoscape")
  install.packages("phateR")
  install.packages("'mpmi")
  
  data("Arabidopsis")
  
  
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install("lattice")
  BiocManager::install("GenomicRanges")
  
  BiocManager::install("SummarizedExperiment")
  install_github("cran/PCIT")
  BiocManager::install("netbenchmark")
  BiocManager::install("grndata")
  BiocManager::install("DREAM4")
  BiocManager::install("RCytoscape")
  BiocManager::install("RCy3")
  BiocManager::install("graph")
  
  source("https://bioconductor.org/biocLite.R")
  biocLite("RCytoscape")
  #install.packages(c("h2o",'XVector','lattice','zlibbioc'))
  
  #browseURL("https://github.com/tschaffter/gnw")
  #browseURL("https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-015-0728-4)
  # retrieving the taxonomic hierarchy of "Arabidopsis thaliana"
  # from the Integrated Taxonomic Information System
  myTAI::taxonomy( organism = "Arabidopsis thaliana", 
                   db       = "itis",
                   output   = "classification" )
}


#Libraries -----
setwd("C:/E/HDDs/Work/NYU Research/0 Plant biology/R code")
#library(officer)

#library(mlbench)
library(caret)
#library(Spectrum)
library(Rtsne)
library(umap)
library(devtools)
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
library(DREAM4)
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

#IDEAS:
#1- directed edge between each gene and all genes in previous time frame
#2- Use a second round RF on CICT predictions, keep most important influencers as causes
#3- replace mutual information with a better including a conditional on all others form
#4- MI (g1,gx.LAGk ) | MI(g1,gn.LAG1)  #MI(g1,gx ge and particular gx using all consequent time frames, conditioned on g1 and all gj in the previous time frame

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

# EXPERIMENTING Curve distance measure Distribution fitting  ----
#https://daviddalpiaz.github.io/stat3202-au18/notes/fitting.html  
#https://cran.r-project.org/web/packages/fitdistrplus/vignettes/FAQ.html
#http://zoonek2.free.fr/UNIX/48_R/07.html
#https://machinelearningmastery.com/probability-density-estimation/
#https://www.scienceforums.net/topic/56448-need-help-with-code-to-calculate-mutual-information/

# EXPERIMENTING Curve distance measure Distribution fitting  ----
#https://daviddalpiaz.github.io/stat3202-au18/notes/fitting.html  
#https://cran.r-project.org/web/packages/fitdistrplus/vignettes/FAQ.html
#http://zoonek2.free.fr/UNIX/48_R/07.html
#https://machinelearningmastery.com/probability-density-estimation/
#https://www.scienceforums.net/topic/56448-need-help-with-code-to-calculate-mutual-information/

# EXPERIMENTING Curve distance measure Distribution fitting  ----
#https://daviddalpiaz.github.io/stat3202-au18/notes/fitting.html  
#https://cran.r-project.org/web/packages/fitdistrplus/vignettes/FAQ.html
#http://zoonek2.free.fr/UNIX/48_R/07.html
#https://machinelearningmastery.com/probability-density-estimation/
#https://www.scienceforums.net/topic/56448-need-help-with-code-to-calculate-mutual-information/

# USE https://en.wikipedia.org/wiki/Fr%C3%A9chet_distance
# fast frechet, check kmlShape package

library(fitdistrplus) #fit distributions
library(gendist)
library(sn)
library(ggplot2)
library(SimilarityMeasures)
library(TSdist) #timeseries distances
library(frechet)

fitg<-fitdist(ussample,"gamma")


x= seq(0,1,.001)# runif(100, min=0, max=1)
y=dskew(x, spec1="norm", arg1=list(mean=0,sd=1), spec2="logis",
        arg2=list(location=0,scale=2) )

y=dskew(x, spec1="norm", arg1=list(mean=.3,sd=.1), spec2="logis",
        arg2=list(location=-3,scale=2) )

y=dskew(x, spec1="norm", arg1=list(mean=.3,sd=.1), spec2="norm",
        arg2=list(mean=.2,sd=.3) )

integrate(dskew,0,1,spec1="norm", arg1=list(mean=.3,sd=.1), spec2="norm",
          arg2=list(mean=.2,sd=.3))

z.l=y/sum(y)

#Testing if area under curve is 1
integrand=function(x){
  y=dskew(x, spec1="norm", arg1=list(mean=.3,sd=.1), spec2="norm",
          arg2=list(mean=.2,sd=.3) )
}
integrate(integrand,0,1)

ggplot(data.frame(x,z.l),aes(x,z.l))+geom_line()

z.r=rev(z.l)
ggplot(data.frame(x,z.r),aes(x,z.r))+geom_line()


z.nd =  dnorm(x,.5,.1)
z.nd = z.nd/sum(z.nd)
ggplot(data.frame(x,z.nd),aes(x,z.nd))+geom_line()

#Regarding upper bounds for I(X;Y), you should write
#I(X;Y)=H(X)−H(X|Y)⟹I(X;Y)≤H(X)≤log(|X|) 
#where |X| is the cardinality of the alphabet for X (n above)

mpmi::cmi.pw(z.nd,z.nd)$mi
mpmi::cmi.pw(z.nd,z.l)$mi
mpmi::cmi.pw(z.r,z.l)$mi

Frechet(matrix(c(x,z.nd),nrow=100),matrix(c(x,z.nd),nrow=100))
Frechet(matrix(c(x,z.nd),nrow=2),matrix(c(x,z.l),nrow=2))
Frechet(matrix(c(x,z.l),nrow=2),matrix(c(x,z.nd),nrow=2))
Frechet(matrix(c(x,z.l),nrow=100),matrix(c(x,z.r),nrow=100))

#TSdist slower calculation
FrechetDistance(z.nd,z.nd,x,x)
FrechetDistance(z.nd,z.nd)
FrechetDistance(z.nd,z.l,x,x)

#frechet package
library(frechet)
library(kmlShape) #Collection of shape respecting distance functions
l1=list(x=x,y=z.nd)
emd = frechet::dist4den(l1,l1,fctn_type='density')




#x1=(100, min=0, max=1)
y1=dskew(x, spec1="norm", arg1=list(mean=.3,sd=.1), spec2="norm",
         arg2=list(mean=.2,sd=.3) )[1:100]
l2=(list(x=x1,y=y1))
emd =frechet::dist4den(l2,l2,fctn_type='density') #Calculates earth moving distance

a=c(1,2,3,4,5,1.5)
b=1:6
ggplot(data.frame(b,a),aes(b,a))+geom_line()

a1= rev(a)
ggplot(data.frame(b,a1),aes(b,a1))+geom_line()

# par <- nlm(function(p){nloglik(p,spec1 = 'logis', arg1 = list(scale = p[1]),
#                                 spec2 = 'logis', arg2 = list(scale = p[2]))}, p = c(1,1))$estimate
# par2<- nlm(function(p){nloglik(p, spec1 = 'norm', arg1 = list(sd = p [1]),
#                                    spec2 = 'norm', arg2 = list(sd = p [2]))}, p = c(1,1))$estimate
qqplot(x,qskew(x, spec1 = "norm", arg1 = list(sd = par2 [1]), spec2 = "norm",
               arg2 = list(sd = par2 [2]), interval = c(-10,10)), xlim = c(0,4),
       ylim = c(0,4), ylab = 'Theoretical', xlab = 'Empirical',
       main = 'Skew Normal-Normal distribution')

qqplot(x,qskew(x, spec1 = "norm", arg1 = list(sd = .2), spec2 = "norm",
               arg2 = list(sd = .3), interval = c(-10,10)), xlim = c(0,1),
       ylim = c(0,1), ylab = 'Theoretical', xlab = 'Empirical',
       main = 'Skew Normal-Normal distribution')

library(sn)
dsn()

curve( dbeta(x,4,2), add=T, lty=2, lwd=2, col='blue' )
curve( dbeta(x,2,4), add=T, lty=3, lwd=3, col='red' )
title(main="Beta distribution")
legend(par('usr')[1], par('usr')[4], xjust=0,
       c('(1,1)', '(2,1)', '(3,1)', '(4,1)', 
         '(2,2)', '(3,2)', '(4,2)',
         '(2,3)', '(3,3)', '(4,3)' ),
       lwd=1, #c(1,1,1,1, 2,2,2, 3,3,3),
       lty=c(1,1,1,1, 2,2,2, 3,3,3),
       col=c(par('fg'), 'red', 'green', 'blue', 
             'red', 'green', 'blue', 
             'red', 'green', 'blue' ))

y=rskew(10, spec1="norm", arg1=list(mean=0,sd=0.1), spec2="logis", 
        arg2=list(location=0,scale=0.2))

# Load data extract ground truth -----
#https://dreamchallenges.org/closed-challenges/


#library(netbenchmark); Availabledata
PSUDEO_ZERO = 0.000001


#Dream4
if(FALSE){
  studyDataset = 'dream4_100_02'
  data(dream4_100_01)
  data(dream4_100_02)
  
  show(dream4_100_02)
  
  extractDataElements<-function(data){
    show(data)
    mtx.all <- assays (data)[[1]]
    
    print(names(metadata(data)))
    mtx.goldStandard  <-  metadata(data)[[1]]
    
    idx  <-  which(mtx.goldStandard  ==  1)
    idx.m1 <- idx -1
    rows  <-  idx.m1  %%  nrow  (mtx.goldStandard)  +  1
    cols  <-  idx.m1  %/%  nrow  (mtx.goldStandard)  +  1
    tbl.goldStandard <- data.frame(src=rownames(mtx.goldStandard)[rows],
                                   trgt=colnames(mtx.goldStandard)[cols],
                                   Source=rep('goldStandard', length(rows)),
                                   stringsAsFactors=FALSE)
    list(mtx.all,tbl.goldStandard)
    #Transform the gold standard matrix into an explict table of regulators and targets.
  }
  
  mtx.all=NULL;tbl.goldStandard=NULL
  
  c(mtx.all,tbl.goldStandard) %<-% extractDataElements(get(studyDataset)) #First of the five networks
  tbl.goldStandard$edgetyptruth = '?'
  dim(mtx.all)
  colnames(mtx.all)
  mtx.all = mtx.all %>%  as.data.frame() %>% rownames_to_column('gene')
  
  #Collapsing mtx.all timeseries
  
  mtx.all=as.data.frame(mtx.all)
  ts.cols = colnames(mtx.all)[(str_detect(colnames(mtx.all),"perturbation"))]
  mf.cols = colnames(mtx.all)[(str_detect(colnames(mtx.all),"MF[.]"))]
  
  #Previously used all columns of all data
  cict.ts = as.data.frame(mtx.all) %>%  dplyr::select(ts.cols) 
  cict.mf = as.data.frame(mtx.all) %>%  dplyr::select(mf.cols) 
  all.dt = cict.mf#All perturbation data of the 20 time series
  
  Nobservations = ncol(all.dt)
}

#Dream5
if(FALSE){
  # In this challenge we explore the use of Systems Genetics data for elucidating causal network models among genes, 
  # i.e. Gene Networks (DREAM5 SYSGEN A) and predicting complex disease phenotypes (DREAM5 SYSGEN B)
  #Challenge https://www.synapse.org/#!Synapse:syn2787209/wiki/70349
  #DATA: https://static-content.springer.com/esm/art%3A10.1038%2Fnmeth.2016/MediaObjects/41592_2012_BFnmeth2016_MOESM584_ESM.zip
  dream5.baseurl = "C:/E/HDDs/Work/NYU Research/0 Plant biology/R code/Data/DREAM5_network_inference_challenge/Network3/"
  
  dr5n1 = readr::read_tsv(paste0(dream5.baseurl,"input data/net3_expression_data.tsv"),col_names =   T)
  all.dt = t(dr5n1) %>% as.data.frame() %>% tibble::rownames_to_column('vID') 
  tbl.goldStandard =  read_tsv(paste0(dream5.baseurl,"gold standard/net3_gold_standard.tsv"),col_names = F) #net1_gold_standard_signed.tsv
  length(unique(c(tbl.goldStandard$src,tbl.goldStandard$trgt)))
  colnames(tbl.goldStandard) <-c('src','trgt','edgetyptruth') 
  
  Nobservations=ncol(all.dt)-1
  ngens = nrow(all.dt)
  
  n.itm.e = data.frame()
  
  dr5n1.sim = similarityMatrices
  a=sapply(similarityMatrices[[i]], function(x) sum(is.na(x)));a[a>0]
  for(i in 1:length(similarityMatrices)){
    tmp=tmp.1=NULL
    tmp = as.data.frame(similarityMatrices[[i]]) %>% tibble::rownames_to_column()
    tmp.1  = pivot_longer(tmp,cols=colnames(tmp)[2:ncol(tmp)]) 
    names(tmp.1)=c('src','trgt',names(similarityMatrices[i]) )
    if(nrow(n.itm.e)==0) n.itm.e <- tmp.1 else n.itm.e <- inner_join(n.itm.e,tmp.1,on=c("src","trgt"))
  }
  
}

#Benchmarking HESC
if(FALSE){
  # In this challenge we explore the use of Systems Genetics data for elucidating causal network models among genes, 
  # i.e. Gene Networks (DREAM5 SYSGEN A) and predicting complex disease phenotypes (DREAM5 SYSGEN B)
  #Challenge https://www.synapse.org/#!Synapse:syn2787209/wiki/70349
  #DATA: https://static-content.springer.com/esm/art%3A10.1038%2Fnmeth.2016/MediaObjects/41592_2012_BFnmeth2016_MOESM584_ESM.zip
  studyDataset = 'hesc'
  hesc.baseurl = "C:/E/HDDs/Work/NYU Research/0 Plant biology/R code/Data"
  
  hseq.ptime = fread(paste0(hesc.baseurl,"/BEELINE-data/inputs/scRNA-Seq/hESC/PseudoTime.csv"),header = T);hist(hseq.ptime$PseudoTime)
  recipes::discretize(hseq.ptime$PseudoTime,cuts=5)
  hseq.ptime = hseq.ptime %>%  dplyr::mutate(cellset = infotheo::discretize(hseq.ptime$PseudoTime, "equalfreq", 5)$X);
  hist(hseq.ptime$cellset);table(hseq.ptime$cellset)
  #hseq.cellset1 = hseq.ptime[PseudoTime<=0.02,]
  
  hsec.gtchip = read.csv(paste0(hesc.baseurl,"/BEELINE-Networks/Networks/human/hESC-ChIP-seq-network.csv"),header = T)
  hsec.gtchip$edgetyptruth = 'chipseq'
  colnames(hsec.gtchip) <-c('src','trgt','edgetyptruth') 
  tbl.goldStandard = hsec.gtchip
  
  hesc = read.csv(paste0(hesc.baseurl,"/BEELINE-data/inputs/scRNA-Seq/hESC/ExpressionData.csv"),header = T)
  
  all.dt = hesc %>% rename(gene = X) # %>% select(any_of('vID',hseq.cellset1$V1))
  rownames(all.dt)<- all.dt$gene 
  #hsec1 =  cbind(X = hesc.sbst1$X, exp = rowMeans(hesc.sbst1[,2:ncol(hesc.sbst1)], dims = 1))
  
  length(unique(c(tbl.goldStandard$src,tbl.goldStandard$trgt)))
  
  Nobservations=ncol(all.dt)-1
  ngens = nrow(all.dt)
  
  a=sapply(similarityMatrices[[i]], function(x) sum(is.na(x)));a[a>0]
  
  #Filtering gene set according to the benchmark paper
  geneorder = read.csv(paste0(hesc.baseurl,"/BEELINE-data/inputs/scRNA-Seq/hESC/GeneOrdering.csv"),header = T)
  geneorder.1e3 = geneorder %>% filter(VGAMpValue<0.01) %>% arrange(desc(Variance)) %>% slice_head(n=1000)
  geneorder.5e2 = geneorder %>% filter(VGAMpValue<0.01) %>% arrange(desc(Variance)) %>% slice_head(n=500)
  
  tf.all = unique(tbl.goldStandard$src) # when Chip_seq , gives number of TFs
  #all TFs whose variance had p-value at most 0.01
  tf.set1= geneorder %>% anti_join(geneorder.1e3,by="X") %>% 
    filter(X %in% tf.all &  VGAMpValue>=0.01) %>% arrange(desc(Variance)) 
  
  
  geneorder.1e3 = geneorder %>% filter(X %in% tf.set1$X) %>% rbind(geneorder.1e3  ) %>% unique()
  
  randomgenes = sample_n(setDT(geneorder)[! X %in% geneorder.1e3$X ],size=1000)
  
  geneorder.2e3 = rbind(geneorder.1e3,randomgenes)
  #n.itm.e = n.itm.e %>% filter(src %in% geneorder.1e3$X | trgt %in% geneorder.1e3$X)
  #all.dt = all.dt %>% filter(gene %in% geneorder.1e3$X | gene %in% geneorder.1e3$X) 
  all.dt = all.dt %>% filter(gene %in% geneorder.2e3$X | gene %in% geneorder.2e3$X) 
  
  write.csv(all.dt,paste0(hesc.baseurl,"/BEELINE-data/inputs/scRNA-Seq/hESC/myFilteredTop1000Varying.csv"))
}

#Benchmarking mHSC Mouse Hematopoietic Stem and Progenitor Cells (mHSC
if(FALSE){
  studyDataset = 'mHSC'
  sc.baseurl = "C:/E/HDDs/Work/NYU Research/0 Plant biology/R code/Data"
  
  sc.ptime = fread(paste0(hesc.baseurl,"/BEELINE-data/inputs/scRNA-Seq/mHSC-E/PseudoTime.csv"),header = T);hist(hseq.ptime$PseudoTime)
  recipes::discretize(sc.ptime$PseudoTime,cuts=5)
  sc.ptime = sc.ptime %>%  dplyr::mutate(cellset = infotheo::discretize(sc.ptime$PseudoTime, "equalfreq", 5)$X);
  hist(sc.ptime$cellset);table(sc.ptime$cellset)
  #sc.cellset1 = sc.ptime[PseudoTime<=0.02,]
  
  sc.gtchip = read.csv(paste0(hesc.baseurl,"/BEELINE-Networks/Networks/mouse/mHSC-ChIP-seq-network.csv"),header = T)
  sc.gtchip$edgetyptruth = 'chipseq'
  colnames(sc.gtchip) <-c('src','trgt', 'score','edgetyptruth') 
  tbl.goldStandard = sc.gtchip
  
  
  sc.gtchip.ns = read.csv(paste0(hesc.baseurl,"/BEELINE-Networks/Networks/mouse/Non-Specific-ChIP-seq-network.csv"),header = T)
  sc.gtchip.ns$edgetyptruth = 'chipseq.ns'
  colnames(sc.gtchip.ns) <-c('src','trgt', 'edgetyptruth') 
  tbl.goldStandard.ns = sc.gtchip.ns
 
  sc.gtstring = read.csv(paste0(hesc.baseurl,"/BEELINE-Networks/Networks/mouse/STRING-network.csv"),header = T)
  sc.gtstring$edgetyptruth = 'string'
  colnames(sc.gtstring) <-c('src','trgt', 'edgetyptruth') 
  tbl.goldStandard.string = sc.gtstring
  
  
  sc = read.csv(paste0(sc.baseurl,"/BEELINE-data/inputs/scRNA-Seq/mHSC-E/ExpressionData.csv"),header = T)
  
  all.dt = sc %>% rename(gene = X) # %>% select(any_of('vID',hseq.cellset1$V1))
  rownames(all.dt)<- all.dt$gene 
  #hsec1 =  cbind(X = sc.sbst1$X, exp = rowMeans(sc.sbst1[,2:ncol(sc.sbst1)], dims = 1))
  
  length(unique(c(tbl.goldStandard$src,tbl.goldStandard$trgt)))
  
  Nobservations=ncol(all.dt)-1
  ngens = nrow(all.dt)
  
  a=sapply(similarityMatrices[[i]], function(x) sum(is.na(x)));a[a>0]
  
  #Filtering gene set according to the benchmark paper
  geneorder = read.csv(paste0(sc.baseurl,"/BEELINE-data/inputs/scRNA-Seq/mHSC-E/GeneOrdering.csv"),header = T)
  geneorder.1e3 = geneorder %>% filter(VGAMpValue<0.01) %>% arrange(desc(Variance)) %>% slice_head(n=1000)
  geneorder.5e2 = geneorder %>% filter(VGAMpValue<0.01) %>% arrange(desc(Variance)) %>% slice_head(n=500)
  
  tf.all = unique(tbl.goldStandard$src) # when Chip_seq , gives number of TFs
  #all TFs whose variance had p-value at most 0.01
  tf.set1= geneorder %>% anti_join(geneorder.1e3,by="X") %>% 
    filter(X %in% tf.all ) %>% arrange(desc(Variance))  #&  VGAMpValue>=0.01
  
  
  geneorder.1e3 = geneorder %>% filter(X %in% tf.set1$X) %>% rbind(geneorder.1e3  ) %>% unique()
  
  randomgenes = sample_n(setDT(geneorder)[! X %in% geneorder.1e3$X ],size=1000)
  
  geneorder.2e3 = rbind(geneorder.1e3,randomgenes)
  #n.itm.e = n.itm.e %>% filter(src %in% geneorder.1e3$X | trgt %in% geneorder.1e3$X)
  all.dt = all.dt %>% filter(gene %in% geneorder.2e3$X | gene %in% geneorder.2e3$X) 
}

#Benchmarking Mouse Embryonic Stem Cells (mESC
if(FALSE){
  studyDataset = 'mESC'
  sc.baseurl = "C:/E/HDDs/Work/NYU Research/0 Plant biology/R code/Data"
  
  sc.ptime = fread(paste0(hesc.baseurl,"/BEELINE-data/inputs/scRNA-Seq/mESC/PseudoTime.csv"),header = T);hist(sc.ptime$PseudoTime)
  recipes::discretize(sc.ptime$PseudoTime,cuts=5)
  sc.ptime = sc.ptime %>%  dplyr::mutate(cellset = infotheo::discretize(sc.ptime$PseudoTime, "equalfreq", 5)$X);
  hist(sc.ptime$cellset);table(sc.ptime$cellset)
  #sc.cellset1 = sc.ptime[PseudoTime<=0.02,]
  
  sc.gtchip = read.csv(paste0(hesc.baseurl,"/BEELINE-Networks/Networks/mouse/mESC-ChIP-seq-network.csv"),header = T)
  sc.gtchip$edgetyptruth = 'chipseq'
  colnames(sc.gtchip) <-c('src','trgt', 'edgetyptruth') 
  tbl.goldStandard = sc.gtchip
  
  
  sc.gtchip.ns = read.csv(paste0(hesc.baseurl,"/BEELINE-Networks/Networks/mouse/Non-Specific-ChIP-seq-network.csv"),header = T)
  sc.gtchip.ns$edgetyptruth = 'chipseq.ns'
  colnames(sc.gtchip.ns) <-c('src','trgt', 'edgetyptruth') 
  tbl.goldStandard.ns = sc.gtchip.ns
  
  sc.gtstring = read.csv(paste0(hesc.baseurl,"/BEELINE-Networks/Networks/mouse/STRING-network.csv"),header = T)
  sc.gtstring$edgetyptruth = 'string'
  colnames(sc.gtstring) <-c('src','trgt', 'edgetyptruth') 
  tbl.goldStandard.string = sc.gtstring
  
  
  sc = read.csv(paste0(sc.baseurl,"/BEELINE-data/inputs/scRNA-Seq/mESC/ExpressionData.csv"),header = T)
  
  all.dt = sc %>% rename(gene = X) # %>% select(any_of('vID',hseq.cellset1$V1))
  rownames(all.dt)<- all.dt$gene 
  #hsec1 =  cbind(X = sc.sbst1$X, exp = rowMeans(sc.sbst1[,2:ncol(sc.sbst1)], dims = 1))
  
  length(unique(c(tbl.goldStandard$src,tbl.goldStandard$trgt)))
  
  Nobservations=ncol(all.dt)-1
  ngens = nrow(all.dt)
  
  a=sapply(similarityMatrices[[i]], function(x) sum(is.na(x)));a[a>0]
  
  #Filtering gene set according to the benchmark paper
  geneorder = read.csv(paste0(sc.baseurl,"/BEELINE-data/inputs/scRNA-Seq/mHSC-E/GeneOrdering.csv"),header = T)
  geneorder.1e3 = geneorder %>% filter(VGAMpValue<0.01) %>% arrange(desc(Variance)) %>% slice_head(n=1000)
  geneorder.5e2 = geneorder %>% filter(VGAMpValue<0.01) %>% arrange(desc(Variance)) %>% slice_head(n=500)
  
  tf.all = unique(tbl.goldStandard$src) # when Chip_seq , gives number of TFs
  #all TFs whose variance had p-value at most 0.01
  tf.set1= geneorder %>% anti_join(geneorder.1e3,by="X") %>% 
    filter(X %in% tf.all ) %>% arrange(desc(Variance))  #&  VGAMpValue>=0.01
  
  
  geneorder.1e3 = geneorder %>% filter(X %in% tf.set1$X) %>% rbind(geneorder.1e3  ) %>% unique()
  
  randomgenes = sample_n(setDT(geneorder)[! X %in% geneorder.1e3$X ],size=1000)
  
  geneorder.2e3 = rbind(geneorder.1e3,randomgenes)
  #n.itm.e = n.itm.e %>% filter(src %in% geneorder.1e3$X | trgt %in% geneorder.1e3$X)
  all.dt = all.dt %>% filter(gene %in% geneorder.2e3$X | gene %in% geneorder.2e3$X) 
}


#Groundtruth Network statistics
{
  library(sna)
  # for building network and visualization 
  library(tidygraph)
  library(graphlayouts)
  library(statnet)
  library(ggplot2)
  #library(ggnetwork)
  
  help(package = sna)
  
  gs.e <- setDT(tbl.goldStandard )
  gs.e[,weight := ifelse(edgetyptruth == '+', 1 , -1)]
  gs.sn=network(gs.e,matrix.type="edgelist",directed=TRUE) 
  
  gs.vnames = network.vertex.names(gs.sn) 
  gs.v = data.frame(vID = unique(c(gs.e$src,gs.e$trg)))
  gs.v = setDT(gs.v)[match(gs.vnames, gs.v$vID),] #Reorders gs.v accorging to network vertices
  sum(gs.v$vID == gs.vnames)
  head(gs.v$vID,10);head(gs.vnames,10)
  
  
  
  # outinhi, outexci, ininhi, inexci
  
  gs.v.outps = tbl.goldStandard %>% group_by(src) %>% summarise(outinhi = 0, outexci=0 )
  gs.v.inps=tbl.goldStandard %>%  group_by(trgt) %>% summarise(outinhi = 0, outexci=0 )#for version 4
  
  gs.v.outps = tbl.goldStandard %>% group_by(src) %>% summarise(outinhi = sum(edgetyptruth=='-'), outexci = sum(edgetyptruth=='+'))
  gs.v.inps = tbl.goldStandard %>% group_by(trgt) %>% summarise(ininhi = sum(edgetyptruth=='-'), inexci = sum(edgetyptruth=='+'))
  
  gs.g=graph.data.frame(gs.e,directed=TRUE)
  gs.v.ind = igraph::degree(gs.g,v=V(gs.g),mode = c("in")); gs.v.ind=data.frame(Indegree=gs.v.ind,subcat=names(gs.v.ind))
  gs.v.outd = igraph::degree(gs.g,v=V(gs.g),mode = c("out")); gs.v.outd=data.frame(Outdegree=gs.v.outd,subcat=names(gs.v.outd))
  
  gs.v=gs.v %>% left_join(gs.v.ind, by=c( "vID" ="subcat")) %>% left_join(gs.v.inps, by = c('vID'='trgt'))
  gs.v=gs.v %>% left_join(gs.v.outd, by=c( "vID" ="subcat"))%>% left_join(gs.v.outps, by = c('vID'='src'))
  
  rm(gs.v.ind,gs.v.inps,gs.v.outd,gs.v.outps)
  
  setDT(gs.v)[ vID %in% tbl.goldStandard$src & !vID %in% tbl.goldStandard$trgt, grp1:=  'cause']
  gs.v[ !vID %in% tbl.goldStandard$src & vID %in% tbl.goldStandard$trgt, grp1:='effect']
  gs.v[ vID %in% tbl.goldStandard$src & vID %in% tbl.goldStandard$trgt, grp1:='both']
  gs.v[ !vID %in% tbl.goldStandard$src & !vID %in% tbl.goldStandard$trgt, grp1:='none']
  gs.v[,color:=mapvalues(gs.v$grp1, from = c('cause','effect','both','none'), to = c('red','lightblue','green','gray' ))]
  
  a=sapply(gs.v, function(x) sum(is.na(x)));a[a>0]
  gs.v = replace_na(gs.v,replace=list(ininhi=0,inexci=0,outinhi=0,outexci=0))
  
  gs.v$eign <- evcent(gs.sn)                          #Computing the eigenvector centrality of each node
  InDegree <- degree(gs.sn, cmode="indegree")     #Computing the in-degree of each node
  InDegree <- InDegree * .25                            #Scaling in-degree to avoid high in-degree nodes from crowding out the rest of the nodes
  #Assigning Attributes to Vertices
  #gs.vnames = get.vertex.attribute(gs.sn,"vertex.names")
  set.vertex.attribute(gs.sn,"grp1",gs.v$grp1)
  set.vertex.attribute(gs.sn,"color",gs.v$color)#,1:nrow(gs.v)
  
  #tst if vertices in the dataset and network match
  if(FALSE){
    gs.v$netname = get.vertex.attribute(gs.sn,"vertex.names")
    gs.v$color = get.vertex.attribute(gs.sn,"color")
  }
  
  tmp = merge(gs.v[,.(vID,color)])
  
  #Step 5: Visualizing the Network
  gs.sn
  summary(gs.sn)                                        #Get numerical summaries of the network
  
  set.seed(12345)
  ggnetwork::ggnetwork(gs.sn) %>%
    ggplot(aes(x = x, y = y, xend = xend, yend = yend)) + 
    geom_edges(color = "lightgray") +
    geom_nodes(color = gs.v$color, size = InDegree) +       
    #geom_nodelabel_repel (color = Race, label = Race) +#   For networks with fewer nodes, we might want to label
    theme_blank() + 
    geom_density_2d()
  
  if(FALSE){ 
    ################################################
    
    
    #Step 3: Calculating Network Measures to Create Network Attributes for Visualization Purposes, More on the Measures Soon
    
    #Step 4: Creating Network Attributes
    #Specifying Colors for Gender and Race
    gs.v <- gs.v %>% 
      mutate (Color_Female = ifelse(sex == 2, 'red', ifelse(sex != 2, 'black', 'black')))
    
    gs.v <- gs.v %>% 
      mutate (Color_Race = ifelse(race5 == 0, 'gold', ifelse(race5 == 1, 'chartreuse4', 
                                                             ifelse(race5 == 2, 'blue1', ifelse(race5 == 3, 'brown', ifelse(race5 == 4, 'purple', 'gray0'))))))
    
    ############################
    #Di-Graph: A set of nodes and arcs. The relationships are directional and can either be dichotomous or weighted. 
    
    #Network: A graph or di-graph where the nodes have attributes assigned to them such as names, genders, or sizes. 
    
    network.size(gs.sn)
    network.edgecount(gs.sn) 
    network.dyadcount(gs.sn) 
    
    #Density: The ratio of Observed Ties/All Possible Ties
    gden(gs.sn, mode = 'digraph')
    
    #Degree Distribution
    #Calculating In-Degree and Out-Degree to Visualize the Total Degree Distribution: What is the distribution of Connectiveness?
    InDegree <- degree(gs.sn, cmode="indegree")     #Computing the in-degree of each node
    OutDegree <- degree(gs.sn, cmode="outdegree")   #Computing the out-degree of each node
    
    par(mar = rep(2, 4))
    par(mfrow=c(2,2)) # Set up a 2x2 display
    hist(InDegree, xlab="Indegree", main="In-Degree Distribution", prob=FALSE)
    hist(OutDegree, xlab="Outdegree", main="Out-Degree Distribution", prob=FALSE)
    hist(InDegree+OutDegree, xlab="Total Degree", main="Total Degree Distribution", prob=FALSE)
    par(mfrow=c(1,1)) # Restore display
    
    #Average Path Length 
    #Walks: A walk is a sequence of nodes and ties, starting and ending with nodes, in which each node is incident with the edges
    #...following and preceding it in the sequence (Wasserman and Faust 1994, p. 105).
    # The beginning and ending node of a walk may be differeent, some nodes may be included more than once, and some ties may be included more than once.
    #Paths: A path is a walk where all the nodes and all the ties are distinct.
    #A shortest path between two nodes is refrred to as a geodesic (Wasserman and Faust 1994, p. 110)
    #Average path length or the geodesic distance is the average number of steps along the shortest paths for all possible pairs of nodes.
    
    # By default, nodes that cannot reach each other have a geodesic distance of infinity. 
    # Because, Inf is the constant for infinity, we need to replace INF values to calculate the shortest path length.
    # Here we replace infinity values with 0 for visualization purposes.
    
    AHS_Geo <- geodist(gs.sn, inf.replace=0)
    #AHS_Geo <- geodist(gs.sn)                #Matrix with Infinity
    
    #AHS_Geo$gdist #The length of the shortest path for all pairs of nodes.
    #AHS_Geo$counts  #The number of shortest path for all pairs of nodes.
    
    Geo_Dist = AHS_Geo$gdist #Shortest Path Matrix
    hist(Geo_Dist)
    
    #Global Clustering Coefficient: Transitivity
    #Transitivity: A triad involving actors i, j, and k is transitive if whenever i --> j and j --> k then i --> k (Wasserman and Faust 1994, p. 243)
    gtrans(gs.sn)
    #Weak and Weak Census
    #Weak transitivity is the most common understanding, the one reflected in Wasserman's and Faust's definition.
    #When 'weak' is specified as the measure, R returns the fraction of potentially intransitive triads obeying the weak condition
    #Transitive Triads/Transtive and Intransitive Triads.
    #In contrast, when 'weak census' is specfified, R returns the count of transitive triads.
    gtrans(gs.sn, mode='digraph', measure='weak')
    gtrans(gs.sn, mode='digraph', measure='weakcensus')
    
    #CUG (Conditional Uniform Graph) Tests:  IS this Graph More Clustered than We Would Expect by Chance
    #See Wasserman and Faust 1994, p. 543-545 for more information.
    #Note: These tests are somewhat computationally intensive.
    #Conducting these tests, we find that athough the transitivity is higher than would be expect by chance given the network's size;
    #...it is not greater than would be expected given either the number of edges or dyads.
    
    #Take long time
    #Test transitivity against size
    if(FALSE)  Cug_Size <- cug.test(gs.sn,gtrans,cmode="size");plot(Cug_Size)
    
    
    #Test transitivity against density
    Cug_Edges <- cug.test(gs.sn,gtrans,cmode="edges");plot(Cug_Edges)
    
    #Test Transitivity against the Dyad Census
    if(FALSE) Cug_Dyad <- cug.test(gs.sn,gtrans,cmode="dyad.census");plot(Cug_Dyad)
    
    ###########################
    #   MESO-LEVEL MEASURES   #
    ###########################
    
    #Dyads
    #Null-Dyads: Pairs of nodes with no arcs between them
    #Asymmetric dyads: Pairs of nodes that have an arc between the two nodes going in one direction or the other, but not both
    #Mutual/Symmetric Dyad: Pairs of nodes that have arcs going to and from both nodes  <--> 
    
    #Number of Symmetric Dyads
    mutuality(gs.sn)
    
    #Dyadic Ratio: Ratio of Dyads where (i,j)==(j,i) to all Dyads
    grecip(gs.sn, measure="dyadic")
    
    #Edgwise Ratio: Ratio of Reciprocated Edges to All Edges
    grecip(gs.sn, measure="edgewise")
    
    #Directed Triad Census
    #Triads can be in Four States
    #Empty: A, B, C
    #An Edge: A -> B, C
    #A Star (2 Edges): A->B->C
    #Closed: A->B->C->A
    
    #Triad types (per Davis & Leinhardt):
    #003  A, B, C, empty triad.
    #012  A->B, C 
    #102  A<->B, C  
    #021D A<-B->C 
    #021U A->B<-C 
    #021C A->B->C
    #111D A<->B<-C
    #111U A<->B->C
    #030T A->B<-C, A->C
    #030C A<-B<-C, A->C.
    #201  A<->B<->C.
    #120D A<-B->C, A<->C.
    #120U A->B<-C, A<->C.
    #120C A->B->C, A<->C.
    #210  A->B<->C, A<->C.
    #300  A<->B<->C, A<->C, completely connected.
    
    triad.census(gs.sn) %>% kable()
    
    #Hierarchy Measures: Components,Cut Points, K-Cores, and Cliques
    #Components: Components are maximally connected subgraphs (Wasserman and Faust 1994, p. 109). 
    #Recall that community 7 has two large components and several small dyads and triads.
    #There are two types of components: strong and weak.
    #Strong components are components connected through directed paths (i --> j, j --> i)
    #Weak components are components connected through semi-paths (--> i <-- j --> k)
    components(gs.sn, connected="strong")
    components(gs.sn, connected="weak")
    
    #Which node belongs to which component?
    AHS_Comp <- component.dist(gs.sn, connected="strong")
    AHS_Comp
    
    AHS_Comp$membership # The component each node belongs to
    AHS_Comp$csize      # The size of each component
    AHS_Comp$cdist      # The distribution of component sizes
    
    #Cut-Sets and Cut-Points: Cut-sets describe the connectivity of the graph based on the removal of nodes, while cut-points describe
    #...the connectivity of the graph based on the removal of lines (Harary 1969)
    #k refers to the number of nodes or lines that would need to be removed to reduce the graph to a disconnected state.
    cutpoints(gs.sn, connected="weak")
    gplot(gs.sn,vertex.col=2+cutpoints(gs.sn,mode="graph",return.indicator=T))
    #The plot only shows subgraphs consisting of nodes with a degree of 2 or more.
    #The green nodes indicate cut-ponts where the removal of the node would separate one subgraph from another.
    
    #Let's remove one of the cutpoints and count components again.
    AHS_Cut <- gs.sn[-11,-11]
    #"-11" selects all the elments in the first row/column.
    #So, AHS_Cut will be gs.sn with node 1 removed.
    
    components(AHS_Cut, connected="strong")  #There are 74 strong components in AHS_Cut compared to 73 in gs.sn
    
    #Bi-Components: Bi-Components refer to subgraphs that require at least the removal of two nodes or two lines to transform it into a 
    #...disconnected set of nodes. 
    #In large highly connected networks, we frequently analyze the properties of the largest bi-component to get a better understanding
    #...of the social system represented by the network.
    bicomponent.dist(gs.sn) 
    
    #Identify Cohesive Subgroups
    #K-Cores: A k-core is a subgraph in which each node is adjacent to at least a minimum number of, k, to the other nodes in the subgraph.
    #..., while a k-plex specifies the acceptable number of lines that can be absent from each node (Wasserman and Faust 1994, p. 266). 
    kcores(gs.sn) 
    #Show the nesting of cores
    AHS_kc<-kcores(gs.sn,cmode="indegree")
    gplot(gs.sn,vertex.col=rainbow(max(AHS_kc)+1)[AHS_kc+1])
    
    #Now, showing members of the 4-core only (All Nodes Have to Have a Degree of 4)
    gplot(gs.sn[AHS_kc>3,AHS_kc>3],vertex.col=rainbow(max(AHS_kc)+1)[AHS_kc[AHS_kc>3]+1])
    
    #Cliques:  A clique is a maximally complete subgraph of three or more nodes.
    #In other words, a clique consists of a subset of nodes, all of which are adjacent to each other, and where there are no other 
    #...nodes that are also adjacent to all of the members of the clique (Luce and Perry 1949)
    
    #We need to symmetrize recover all ties between i and j.
    set.network.attribute(gs.sn, "directed", FALSE) 
    
    #The clique census returns a list with several important elements 
    #Let's assign that list to an object we'll call AHS_Cliques.
    #The clique.comembership parameter takes values "none" (no co-membership is computed),
    #"sum" (the total number of shared cliques for each pair of nodes is computed),
    #bysize" (separate clique co-membership is computed for each clique size)
    
    AHS_Cliques <- clique.census(gs.sn, mode = "graph", clique.comembership="sum")
    AHS_Cliques # an object that now contains the results of the clique census
    
    #The first element of the result list is clique.count: a matrix containing the number of cliques of different 
    #...sizes (size = number of nodes in the clique).
    #The first column (named Agg) gives you the total  number of cliqies of each size,
    #The rest of the columns show the number of cliques each node participates in.
    
    #Note that this includes cliques of sizes 1 & 2. We have those when the largest fully connected structure includes just 1 or 2 nodes.
    AHS_Cliques$clique.count
    
    #The second element is the clique co-membership matrix:
    View(AHS_Cliques$clique.comemb)
    
    # The third element of the clique census result is a list of all found cliques:
    # (Remember that a list can have another list as its element)
    AHS_Cliques$cliques # a full list of cliques, all sizes
    
    AHS_Cliques$cliques[[1]] # cliques size 1
    AHS_Cliques$cliques[[2]] # cliques of size 2
    AHS_Cliques$cliques[[3]] # cliques of size 3
    AHS_Cliques$cliques[[4]] # cliques of size 4
    
  }
  ###########################
  #   NODE LEVEL MEASURES   #
  ###########################
  
  #Restoring Our Directed Network
  set.network.attribute(gs.sn, "directed", TRUE) 
  
  #Reachability
  #An actor is "reachable" by another if there exists any set of connections by which we can trace from the source to the target actor, 
  #regardless of how many other nodes fall between them (Wasserman and Faust 1994, p. 132).
  #If the network is a directed network, then it possible for actor i to be able to reach actor j, but for j not to be able to reach i.
  #We can classify how connected one node is to another by considering the types of paths connecting them.
  #Weakly Connected: The nodes are connected by a semi-path (--> i <--- j ---> k)
  #Unilaterally Connected: The nodes are connected by a path (i --> j --> k)
  #Strongly Connected: The nodes are connected by a path from i to k and a path from k to i.
  #Recursively Connected: The nodes are strongly connected, and the nodes along the path from i to k and from k to i are the same in reverse order.
  #e.g., i <--> j <--> k 
  
  #Interpreting the reachability matrix, the first column indicates a specific node, the second an alter (alters can occur multiple times),
  #and the third column indicates the number of paths connecting the two (total is a cumulative count of the number of paths in the network).
  #For example, interpreting row 2, node 2 can reach node 235 through 235 paths (470-235), whereas in the middle of the list node 343 can reach node 1 through only 1 path.
  reachability(gs.sn) 
  ??reachablity #For more information on this measure
  
  #Degree Centraltiy: Total, In-Degree, Out-Degree
  
  #In-Degree Centrality: The number of nodes adjacent to node i (Wasserman and Faust 1994, p. 126). i <--
  InDegree <- degree(gs.sn, cmode="indegree")
  InDegree <- InDegree * .15                #Scaling in-degree to avoid high in-degree nodes from crowding out the rest of the nodes
  
  set.vertex.attribute(gs.sn, "InDegree", InDegree)
  
  #Out-Degree Centrality: The number of nodes adjacent from node i (Wasserman and Faust, p. 126). i -->
  OutDegree <- degree(gs.sn, cmode="outdegree")
  OutDegree <- OutDegree * .5                 #Scaling in-degree to avoid high in-degree nodes from crowding out the rest of the nodes
  
  set.vertex.attribute(gs.sn, "OutDegree", OutDegree)
  
  #Total Degree Centrality: The Total Number of Adjacent Nodes (In-Degree + Out-Degree)
  gs.v$TotalDegree <- degree(gs.sn, cmode="outdegree") + degree(gs.sn, cmode="indegree")
  TotalDegree <- gs.v$TotalDegree * .4
  
  set.vertex.attribute(gs.sn, "TotalDegree", TotalDegree)
  
  #Try Sizing by the Different Degrees
  set.seed(12345)
  ggnetwork::ggnetwork(gs.sn) %>%
    ggplot(aes(x = x, y = y, xend = xend, yend = yend)) + 
    geom_edges(color = "lightgray") +
    geom_nodes(color = gs.v$color, size = InDegree*3) +       
    #geom_nodelabel_repel (color = Race, label = Race) +#   For networks with fewer nodes, we might want to label
    theme_blank() + 
    geom_density_2d()
  
  #Path Centralities: Closeness Centrality, Information Centrality, Betweenness Centrality
  
  #Closeness Centrality: Closeness centrality measures the geodesic distances of node i to all other nodes.
  #Functionally, this measures range from 0 to 1, and is the inverse average distance between actor i and all other actors (Wasserman and Faust 1994, p. 185)
  #This measure does not work well when there are disconnected components because the distances between components cannot be summed as
  #...they are technically infinite. There are several work arounds, see Acton and Jasny's alternative below.
  
  AHS_Closeness <- closeness(gs.sn, gmode="digraph", cmode="directed")
  
  hist(AHS_Closeness , xlab="Closness", prob=TRUE) 
  
  #Alternative Approach to Measuring Closesness from the Geodesic Distances Matrix from Acton's and Jasny's Statnet Tutorial
  Closeness <- function(x){ # Create an alternate closeness function!
    geo <- 1/geodist(x)$gdist # Get the matrix of 1/geodesic distance
    diag(geo) <- 0 # Define self-ties as 0
    apply(geo, 1, sum) # Return sum(1/geodist) for each vertex
  }
  
  gs.v$closeness <-  Closeness(gs.sn) ;setDT(gs.v)                       #Applying the function
  #View(gs.v[closeness!=0,])
  hist( gs.v$closeness , xlab="Alt. Closeness", prob=TRUE)         #Better behaved!
  
  #Information Centrality: Information Centrality measures the information flowing from node i.
  #In general, actors with higher information centrality are predicted to have greater control over the flow of information within a network.
  #Highly information-central individuals tend to have a large number of short paths to many others within the social structure.
  ?infocent  #For more information
  
  gs.v$infocent <- infocent(gs.sn, rescale=TRUE)
  hist(gs.v$infocent , xlab="Information Centrality", prob=TRUE) 
  
  gplot(gs.sn, vertex.cex=(AHS_Info)*250, gmode="graph") # Use w/gplot
  #As suggested by the histogram there is relatively little variation in information centrality in this graph.
  
  #Betweenness Centrality: The basic intuition behind Betweenness Centrality is that the actor between all the other actors in the 
  #...has some control over the paths in the network. 
  #Functionally, Betweenness Centrality is the ratio of the sum of all shortest paths linking j and k that includes node i over 
  #...all the shortest paths linking j and k (Wasserman and Faust 1994, p. 191)
  
  gs.v$betweenness <- betweenness(gs.sn, gmode="digraph")  
  hist(gs.v$betweenness , xlab="Betweenness Centrality", prob=TRUE) 
  
  gplot(gs.sn, vertex.cex=sqrt(gs.v$betweenness)/25, gmode="digraph") 
  
  #Comparing Closeness and Betweenness Centralities
  cor(gs.v$closeness, gs.v$betweenness)                  #Correlate our adjusted measure of closeness with betweenness
  plot(gs.v$closeness, gs.v$betweenness)                            #Plot the bivariate relationship
  
  #Measures of Power in Influence Networks: Bonachich and Eigenvector Centrality
  
  #Bonachich Centrality: The intuition behind Bonachich Power Centrality is that the power of node i is recursively defined 
  #...by the sum of the power of its alters. 
  #The nature of the recursion involved is then controlled by the power exponent: positive values imply that vertices become 
  #...more powerful as their alters become more powerful (as occurs in cooperative relations), while negative values imply 
  #...that vertices become more powerful only as their alters become weaker (as occurs in competitive or antagonistic relations).
  ?bonpow   #For more information about the measure
  
  #Eigenvector Centrality: Conceptually, the logic behind eigenvectory centrality is that node i's influence is proportional to the 
  #...to the centraltities' of the nodes adjacent to node i. In other words, we are important because we know highly connected people.
  #Mathematically, we capture this concept by calculating the values of the first eigenvector of the graph's adjacency matrix.
  ?evcent   #For more information.
  
  gs.v$eigen <- evcent(gs.sn)
  hist(gs.v$eigen , xlab="Eigenvector Centrality", prob=TRUE) 
  
  gplot(gs.sn, vertex.cex=gs.v$eigen*10, gmode="digraph") 
  
  ###########################
  #   POSITIONAL ANALYSIS   #
  ###########################
  
  #Burt's (1992) measures of structural holes are supported by iGraph and ego network variants of these measures are supported by egonet
  #...the egonet package is compatable with the sna package.
  
  #You can find descriptions and code to run Burt's measures in igraph at: http://igraph.org/r/doc/constraint.html
  
  #Brokerage: The brokerage measure included in the SNA package builds on past work on borkerage (Marsden 1982), but is a more 
  #...explicitly group oriented measure. Unlike Burt's (1992) measure, the Gould-Fernandez measure requires specifying a group variable
  #...based on an attribute. I use race in the example below.
  
  #Brokerage Roles: Group-Based Concept
  #w_I: Coordinator Role (Mediates Within Group Contact)
  #w_O: Itinerant Broker Role (Mediates Contact between Individuals in a group to which the actor does not belong)
  #b_{IO}: Representative: (Mediates incoming contact from out-group members)
  #b_{OI}: Gatekeeper: (Mediates outgoing contact from in-group members)
  #b_O: Liason Role: (Mediates contact between individuals of two differnt groups, neither of which the actor belongs)
  #t: Total or Cumulative Brokerage (Any of the above paths)
  ?brokerage   #for more information
  
  AHS_Brokerage <- brokerage(gs.sn, gs.v$grp1)
  gs.v = cbind(gs.v, AHS_Brokerage$z.nli)
  AHS_Brokerage$raw.gli
  hist(AHS_Brokerage$cl, xlab="Cumulative Brokerage", prob=TRUE) 
  
  AHS_CBrokerage <- (AHS_Brokerage$cl)
  gplot(gs.sn, vertex.cex=AHS_CBrokerage*.5, gmode="digraph") 
  
  #Structural Equivalence
  #Structural equivalence: Similarity/Distance Measures Include:
  #Correlation
  #Euclidean Distance
  #Hamming Distance
  #Gamma Correlation
  sedist(gs.sn, mode="digraph", method="correlation")
  
  #Cluster based on structural equivalence:
  AHS_Clustering <- equiv.clust(gs.sn, mode="digraph",plabels=network.vertex.names(gs.sn))
  AHS_Clustering                        #Specification of the equivalence method used
  plot(AHS_Clustering)                  #Plot the dendrogram
  rect.hclust(AHS_Clustering$cluster, h=30)
  
  #Generating a Block Model based on the Structural Equivalence Clustering
  AHS_BM <- blockmodel(gs.sn, AHS_Clustering, h=200)
  gs.v$clst_block = AHS_BM$block.membership
  unique(gs.v$clst_block)
  
  
  #Extract the block image for Visualization
  bimage <- AHS_BM$block.model
  bimage
  bimage[is.nan(bimage)] <- 1
  
  #Visualizing the block image (with self-reflexive ties)
  gplot(bimage, diag=TRUE, edge.lwd=bimage*5, vertex.cex=sqrt(table(AHS_BM$block.membership))/2,
        gmode="graph", vertex.sides=50, vertex.col=gray(1-diag(bimage)))
  
  rm(AHS_BM,AHS_Brokerage,AHS_kc,AHS_Cut,AHS_Comp,AHS_Geo,AHS_BM,AHS_CBrokerage,AHS_Clustering,h2o.pred)
}

#Mutual information steady state Multiple measures parallel
{
  #Parallel partition for each measure, multiple measures
  if(TRUE){ 
    print(" - Processing in parallel - ")
    #s0m3
    library(doParallel);
    
    # getSimilarityMatrix_MI = function(actualDataset,actualDatasetNNodes,estimators, discretization = TRUE, discretizator){
    #   build.mim(actualDataset, disc = "equalwidth", estimator = "mi.mm")
    # }
    
    # scoreMatrices[[sc]] <- getNormalizedMatrix(scoreMatrices[[sc]],normalization="minimax");
    
    source("c:/E/HDDs/Work/NYU Research/0 Plant biology/R code/Comparison/vkuzmanovski-rn-approach/rnR_Framework.R")
    
    actualDataset <-all.dt #genes in rows and cells in columns  #  stop('Set correct data source') #  all.tdt
    actualDatasetNNodes <- nrow(actualDataset) + 1;
    actualDatasetNObservations <- ncol(actualDataset);
    actualDatasetName <- "hseq";
    actualDatasetSymbolicPatterns=0;
    actualDatasetPatterns=0
    
    simsGroup=c("MI","CORR","DIST")
    availableGroups <- c("MI","CORR","DIST","SYM");
    availableSimilarities <- vector("list", length(availableGroups));
    names(availableSimilarities) <- availableGroups;
    
    availableSimilarities$MI <- c("efMImm","ewMImm","efMIempirical","ewMIempirical", "efMIshrink","ewMIshrink");
    availableSimilarities$CORR <- c("Pearson","Spearman") #,"Kendall");
    availableSimilarities$DIST <- c("Manhattan","Euclidean","L10Norm");
    #availableSimilarities$SYM <- c("efSym","efSymMI","ewSymMI","efAvgSymMI","ewAvgSymMI","Qual");
    sims <- unlist(lapply(simsGroup, function(group){availableSimilarities[[group]]}), recursive=TRUE, use.names=FALSE);
    
    
    library(doSNOW)
    #processingCluster <-makeCluster(no_cores =5,  outfile = paste("./log/global_prediction_",actualDatasetName,"_.log",sep=""));
    processingCluster<-makeCluster(14,outfile = paste("./log/global_prediction_",actualDatasetName,"_.log",sep=""))
    registerDoSNOW(processingCluster)
    
    clusterEvalQ(processingCluster, library(infotheo));
    clusterEvalQ(processingCluster, library(igraph));
    clusterEvalQ(processingCluster, source("c:/E/HDDs/Work/NYU Research/0 Plant biology/R code/Comparison/vkuzmanovski-rn-approach/rnR_Framework.R"));
    
    similarityMatrices <- parSapply(processingCluster, sims, simplify = FALSE, USE.NAMES = TRUE,
                                    function(sim, actualDataset, actualDatasetNNodes, actualDatasetNObservations, actualDatasetName, actualDatasetSymbolicPatterns, patterns, numCores){
                                      print(paste("[",Sys.time(),"]","Processing",sim,"over",actualDatasetName,"...",sep=" "));
                                      browser()
                                      ##Perform similarity/distance step
                                      
                                      firstStepMatrix <- tryCatch({
                                        
                                        locMatrix <- switch(sim,
                                                            efMImm = getSimilarityMatrix_MI(actualDataset,actualDatasetNNodes,estimators="mi.mm", discretization = TRUE, discretizator = "equalfreq"),
                                                            ewMImm = getSimilarityMatrix_MI(actualDataset,actualDatasetNNodes,estimators="mi.mm", discretization = TRUE, discretizator = "equalwidth"),
                                                            efMIempirical = getSimilarityMatrix_MI(actualDataset,actualDatasetNNodes,estimators="mi.empirical", discretization = TRUE, discretizator = "equalfreq"),
                                                            ewMIempirical = getSimilarityMatrix_MI(actualDataset,actualDatasetNNodes,estimators="mi.empirical", discretization = TRUE, discretizator = "equalwidth"),
                                                            efMIshrink = getSimilarityMatrix_MI(actualDataset,actualDatasetNNodes,estimators="mi.shrink", discretization = TRUE, discretizator = "equalfreq"),
                                                            ewMIshrink = getSimilarityMatrix_MI(actualDataset,actualDatasetNNodes,estimators="mi.shrink", discretization = TRUE, discretizator = "equalwidth"),
                                                            Pearson = getSimilarityMatrix_MI(actualDataset,actualDatasetNNodes,estimators="pearson"),
                                                            Spearman = getSimilarityMatrix_MI(actualDataset,actualDatasetNNodes,estimators="spearman"),
                                                            Kendall = getSimilarityMatrix_MI(actualDataset,actualDatasetNNodes,estimators="kendall"),
                                                            Manhattan = getSimilarityMatrix_DISTANCES(actualDataset,actualDatasetNNodes,norms=1),
                                                            Euclidean = getSimilarityMatrix_DISTANCES(actualDataset,actualDatasetNNodes,norms=2),
                                                            L10Norm = getSimilarityMatrix_DISTANCES(actualDataset,actualDatasetNNodes,norms=10),
                                                            DTWasym = getSimilarityMatrix_DTW(actualDataset,steppatern=asymmetric),
                                                            DTWsym1 = getSimilarityMatrix_DTW(actualDataset,steppatern=symmetric1),
                                                            DTWsym2 = getSimilarityMatrix_DTW(actualDataset,steppatern=symmetric2),
                                                            efSym = getSimilarityMatrix_SYMBOLIC(actualDataset,actualDatasetNNodes,actualDatasetNObservations, simmethod="sym", npatterns=actualDatasetSymbolicPatterns, discretizator = "equalfreq", patterns=patterns, numCores),
                                                            ewSym = getSimilarityMatrix_SYMBOLIC(actualDataset,actualDatasetNNodes,actualDatasetNObservations, simmethod="sym", npatterns=actualDatasetSymbolicPatterns, discretizator = "equalwidth", patterns=patterns, numCores),
                                                            efSymMI = getSimilarityMatrix_SYMBOLIC(actualDataset,actualDatasetNNodes,actualDatasetNObservations,  simmethod="sym.mi", npatterns=actualDatasetSymbolicPatterns, discretizator = "equalfreq", patterns=patterns, numCores),
                                                            ewSymMI = getSimilarityMatrix_SYMBOLIC(actualDataset,actualDatasetNNodes,actualDatasetNObservations,  simmethod="sym.mi", npatterns=actualDatasetSymbolicPatterns, discretizator = "equalwidth", patterns=patterns, numCores),
                                                            efAvgSymMI = getSimilarityMatrix_SYMBOLIC(actualDataset,actualDatasetNNodes,actualDatasetNObservations,  simmethod="avg.sym.mi", npatterns=actualDatasetSymbolicPatterns, discretizator = "equalfreq", patterns=patterns, numCores),
                                                            ewAvgSymMI = getSimilarityMatrix_SYMBOLIC(actualDataset,actualDatasetNNodes,actualDatasetNObservations,  simmethod="avg.sym.mi", npatterns=actualDatasetSymbolicPatterns, discretizator = "equalwidth", patterns=patterns, numCores),
                                                            Qual = getSimilarityMatrix_QUAL(actualDataset,actualDatasetNNodes,actualDatasetNObservations)
                                                            #       efMIml_gc = getSimilarityMatrix_MI(actualDataset,actualDatasetNNodes,estimators="coarse.grained", subestimators="ML",discretizator = "equalfreq"),
                                                            #       ewMIml_gc = getSimilarityMatrix_MI(actualDataset,actualDatasetNNodes,estimators="coarse.grained", subestimators="ML",discretizator = "equalwidth"),
                                                            #       efMImm_gc = getSimilarityMatrix_MI(actualDataset,actualDatasetNNodes,estimators="coarse.grained", subestimators="mm",discretizator = "equalfreq"),
                                                            #       ewMImm_gc = getSimilarityMatrix_MI(actualDataset,actualDatasetNNodes,estimators="coarse.grained", subestimators="mm",discretizator = "equalwidth"),
                                                            #       efMIshrink_gc = getSimilarityMatrix_MI(actualDataset,actualDatasetNNodes,estimators="coarse.grained", subestimators="shrink",discretizator = "equalfreq"),
                                                            #       ewMIshrink_gc = getSimilarityMatrix_MI(actualDataset,actualDatasetNNodes,estimators="coarse.grained", subestimators="shrink",discretizator = "equalwidth"),
                                                            #       Granger = getSimilarityMatrix_MI(actualDataset,actualDatasetNNodes,estimators="granger")
                                        )
                                        
                                        row.names(locMatrix) <- row.names(actualDataset);
                                        colnames(locMatrix) <- row.names(actualDataset);
                                        print(paste("[",Sys.time(),"]","DONE processing",sim,"over",actualDatasetName,"...",sep=" "));
                                        return(locMatrix);
                                        
                                      },error=function(err){
                                        print("Error thrown in thread!")
                                        print(err);
                                        return(NULL);
                                      },warning=function(warn){
                                        print("Warning thrown in thread!")
                                        print(warn);
                                        return(NULL);
                                      });
                                      
                                      
                                      firstStepMatrix
                                      #getNormalizedMatrix(firstStepMatrix,normalization="minimax");
                                      
                                    },
                                    actualDataset, actualDatasetNNodes, actualDatasetNObservations, actualDatasetName, actualDatasetSymbolicPatterns, actualDatasetPatterns, 4); #max(as.numeric((detectCores()-no_cores)/6)-1,1)
    stopCluster(processingCluster)
    
    similarityMatrices=list(res)
    
    n.itm.e = data.frame()
    for(i in 1:length(similarityMatrices)){
      tmp=tmp.1=NULL
      tmp = as.data.frame(similarityMatrices[[i]]) %>% tibble::rownames_to_column()
      tmp.1  = pivot_longer(tmp,cols=colnames(tmp)[2:ncol(tmp)]) 
      names(tmp.1)=c('src','trgt',names(similarityMatrices[i]) )
      if(nrow(n.itm.e)==0) n.itm.e <- tmp.1 else n.itm.e <- inner_join(n.itm.e,tmp.1,on=c("src","trgt"))
    }
    
    #new cols paste0(colnames(n.itm.e),collapse="','")
    #c('efMImm','ewMImm','efMIempirical','ewMIempirical','efMIshrink','ewMIshrink','Pearson','Spearman','Kendall','Manhattan','Euclidean','L10Norm')
    
  }
  
  #Old implementation
  if(FALSE){
      mf.mi = mpmi::cmi.pw(v1,v2)$mi,
      corP = cor(v1,v2,method='pearson'),
      corK = cor(v1,v2,method='kendall'),
      corS = cor(v1,v2,method='spearman')   
  
  
}

  #Parallel partition for gene set, one or more measures
  if(TRUE){
    library(doSNOW)
    library(foreach)
    library(doFuture)
    library(progressr)
    library(parallel)
   
    myMI <- function(X,Y,methode="emp",discretizers="equalwidth")
    {
      library("infotheo")
      
      Xd <- unlist(discretize(X,disc=discretizers))
      Yd <- unlist(discretize(Y,disc=discretizers))
      XYd <- array(0,c(length(X),2))
      XYd[,1] <- Xd
      XYd[,2] <- Yd
      
      I <- entropy(Xd,method=methode) + entropy(Yd,method=methode) - entropy(XYd,method=methode)
      return(I)
    }
    
    
    myOuterJoin=function(a,b){
      setDT(a);setDT(b);setkey(a,VISITLINK);setkey(b,VISITLINK)
      merge(a,b, all=TRUE)
      
    }
    
    
    handlers(global = TRUE)
    #handlers("progress", "beepr")
    #handlers(handler_progress(complete = "#"))
    
    all.tdt = all.dt %>% data.frame() %>% tibble::column_to_rownames('vID') %>% mutate(across(!any_of('vID'),~as.numeric(.x)))%>%  t()
    ngene = ncol(all.tdt)
    #genesets = split(all.dt$gene[1:26],1:ncores)
    
    ncores = detectCores()-3
    cl=makeCluster(ncores) # 
    registerDoFuture()  #registerDoSNOW(cl)
    plan(multisession)
    
    library(progressr)
    handlers('progress')
    with_progress({    
      
      prgrs<-progressor(steps = ngene)
      res <- foreach(i = 1:ngene,
                     .combine = rbind,
                     .multicombine = TRUE,
                     .inorder = FALSE,
                     .verbose = F,
                     .packages = c('data.table', 'doParallel','minet','mpmi','progressr')) %:%
        foreach(j=1:ngene, .combine='cbind') %dopar% 
        {
          prgrs()        
          
          a=data.frame(corK = cor(all.tdt,all.tdt,method='kendall'), corP = cor(all.tdt[,i],all.tdt[,j],method=c("pearson")) ,
                       mf.mi = myMI(all.tdt[,i],all.tdt[,j]) ) #
          mpmi::cmi.pw(all.tdt[,i],all.tdt[,j])$mi
          a
          # 
          # corS = cor(all.tdt,all.tdt,method='spearman')
          #mim <- build.mim((all.tdt), estimator='pearson')
          #mim
        }
      row.names(res) <- colnames(all.tdt[,1:ngene]);
      colnames(res) <- tidyr::crossing(x=colnames(all.tdt[,1:ngene]),y=colnames(a)) %>% unite(z,x:y) %>% pull(z); #View(head(res))
    })
    
    calcor <- function(expm,ngene=NA) {
      if(is.na(ngene)) ngene = ncol(expm)
      p <- progressr::progressor(along = ngene)
      res <- foreach(i = 1:ngene,
                     .combine = rbind,
                     .multicombine = TRUE,
                     .inorder = FALSE,
                     .packages = c('data.table', 'doParallel','minet','mpmi','progressr')) %dopar% {
                       # mf.mi = mpmi::cmi.pw(v1,v2)$mi,
                       corP = cor(expm[,i],expm,method='pearson')
                       
                       # corK = cor(all.tdt,all.tdt,method='kendall'),
                       # corS = cor(all.tdt,all.tdt,method='spearman')
                       
                       
                       #mim <- build.mim((all.tdt), estimator='pearson')
                       #mim
                       
                       p(message = sprintf("Added %g", 1))
                     }
      res
    }
    #system.time({res<-calcor(all.tdt,ngene = 1000) })
    
    
    stopCluster(cl)
    
    similarityMatrices = list(corP=res)
    n.itm.e = data.frame()
    for(i in 1:length(similarityMatrices)){
      tmp=tmp.1=NULL
      tmp = as.data.frame(similarityMatrices[[i]]) %>% tibble::rownames_to_column() #View(head(tmp))
      tmp.1  = pivot_longer(tmp,cols=colnames(tmp)[2:ncol(tmp)]) 
      names(tmp.1)=c('src','trgt',names(similarityMatrices[i]) )#View(head(tmp.1))
      if(nrow(n.itm.e)==0) n.itm.e <- tmp.1 else n.itm.e <- inner_join(n.itm.e,tmp.1,on=c("src","trgt"))
    }    
    #n.itm.e = readRDS("baseline n.itm.e.rds")
    rm(tmp,tmp.1,similarityMatrices)
    rm(hesc,hsec.gt,hsec1,hseq.cellset1,ind,outd, all.tdt)
  }
}

#Mutual information with previous timeframe
{
  #Hesc
  maxpsteps = max(hseq.ptime$cellset)
  maxCellRank = min(table(hseq.ptime$cellset)) # makes all cellsets having equal number of cells
  hseq.ptime = hseq.ptime %>% group_by(cellset) %>% mutate(cellrank = row_number()) %>% ungroup()
  all.dt=as.data.frame(all.dt)  %>% dplyr::mutate(geneId = row_number()) %>% select(geneId,everything())
  
  mtx.lng.prtrb = all.dt %>%
    dplyr::select(-geneId,-vID) %>% select(gene,everything()) %>% #tibble::rownames_to_column(var = "vID") %>%
    pivot_longer(cols = !matches('gene'),names_to = "celltags") %>%
    left_join(hseq.ptime, by = c('celltags'='V1')) %>% #Assigns each cell to a psuedotime tag
    filter(cellrank<=maxCellRank) %>% #make sure all celsets have the same number of cells
    dplyr::mutate(pstep= cellset,prtrb = cellrank)  %>%  #prtrb just used to ensure a one by one cell correspondance in each cell set to prevent cascade relations in the e.ab = inner_join
    dplyr::select(celltags,prtrb,pstep,gene,value) %>% ungroup()
  
  table(mtx.lng.prtrb$prtrb)
  
  mtx.lng.prtrb = mtx.lng.prtrb %>% group_by(pstep) %>%
    mutate(pstep.lead = ifelse(pstep==maxpsteps,maxpsteps+1,pstep+1))  %>%ungroup()
  
  
  
  #Dream4, in each prturbation, each step with the next step
  {
    mtx.lng =as.data.frame(mtx.all) %>% tibble::rownames_to_column(var = "gene") %>%
      pivot_longer(cols = setdiff(colnames(mtx.all),'gene'),names_to = "set") %>%
      dplyr::select(set,gene,value) %>% 
      dplyr::mutate(prtrb = as.double(str_extract(set,'(?<=perturbation[.])[0-9]+')),
                    ptime =as.double( str_extract(set,'[0-9]+$'))) %>%
      dplyr::mutate(pstep= ptime/50) 
    
    mtx.lng.prtrb =mtx.lng %>%dplyr::filter(!is.na(prtrb)) %>% #arrange(prtrb,pstep) %>%
      group_by(prtrb) %>%
      mutate(pstep.lead = dplyr::lead(pstep,1)) %>% 
      mutate(pstep.lead = ifelse(pstep==20,NA,pstep.lead)) %>% ungroup()
    
    
    }
  setDT(mtx.lng.prtrb)
  setkeyv(mtx.lng.prtrb, c('gene','prtrb','pstep','pstep.lead'))
  
  
  library(doSNOW)
  library(foreach)
  library(doFuture)
  
  ncores = detectCores()-3
  cl=makeCluster(ncores) # 
  registerDoSNOW(cl)
  #registerDoFuture();plan(multisession) #Didnt work for nested foreach
  

  #n.itm.e = data.table(idx=1:(ngens*(ngens-1)))
  
  #Parallel
  #New implementations
  {
    sort( sapply(ls(),function(x){object.size(get(x))}),decreasing = F)
    #rm(n.itm.e.back3,n.itm.e.BeforeAugmenting,n.itm.e.back,t2,tst1.totalset,d.new,tst1.rnd);gc()
    library(progressr)
    library(mpmi)
    my_cmi<-function(x,y){
      tryCatch({
        #minerva::mine(x,y)$MIC
        mpmi::cmi(unlist(x),unlist(y))$mi
      },error = function(e) 0 ) # list(mi=0,bcmi=0,zvalues=0))
    }
    ngens = nrow(all.dt)
    thisrow = 0;j=1;i=1
    e.ts = data.frame()
    setDT(n.itm.e)
    
    progressr::handler_progress()
    pb <- txtProgressBar(max = ngens, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)

    system.time({
      e.ts<- foreach(firstg = all.dt$gene,.combine=rbind,.packages =c( "data.table",'mpmi','dplyr'),
                   .options.snow = list(progress = progress)) %dopar%{
        e.a=e.b=e.ab=e.sum=NULL
        if(firstg!=secondg) {
          e.a = mtx.lng.prtrb[ gene == firstg ,.(gene,prtrb,value,pstep.lead)]
          #e.b = mtx.lng.prtrb[ gene == secondg,.(gene,prtrb,value,pstep)]
          e.ab = merge(e.a,mtx.lng.prtrb[gene != firstg ,] ,
                       by.x=c('prtrb','pstep.lead'), by.y=c('prtrb','pstep'))[,.(gene.y,value.x,value.y)]

          e.sum = setDF(e.ab) %>% dplyr::group_by(gene.y) %>% 
            summarise(tlag1.mi = my_cmi(value.x,value.y),
                       tlag1.corP = cor(value.x,value.y,method='pearson'))
          

          e.sum =e.sum %>% mutate(first = firstg) %>% rename(second = gene.y)
          #actualDatasetNNodes <- nrow(actualDataset) + 1;
          #getSimilarityMatrix_MI(e.ab,actualDatasetNNodes,estimators="mi.empirical", discretization = TRUE, discretizator = "equalwidth")
          
        } else(e.sum=data.frame(first=NA,second=NA,tlag1.mi=NA,tlag1.corP=NA))
        e.sum

                   }
    })
    #firstg = all.dt$gene[1];headgenes = all.dt$gene[1:500]
    #for foreach strategies Check https://github.com/HenrikBengtsson/future/issues/95
    # e.ts<- foreach(firstg = all.dt$gene,.combine=function(...) rbindlist(list(...)),
    #                .multicombine=TRUE) %:%
    #   foreach(secondg =  all.dt$gene,.combine=rbind,.packages =c( "data.table",'mpmi')) 
    
    close(pb)
    stopCluster(cl)
      
    #Todo Remove
    tryCatch({
    #   
    #   if(!(myMean(e.ab$value.y)==0 | myMean(e.ab$value.x)==0 ) & 
    #      length(e.ab$value.x)>0 & length(e.ab$value.y)>0) {
    #     
    #     therow=data.frame(first=firstg,second=secondg,
    #                       tlag1.mi = mpmi::cmi.pw(e.ab$value.x,e.ab$value.y)$mi,
    #                       tlag1.corP = cor(e.ab$value.x,e.ab$value.y,method='pearson')
    #                       #tlag1.corK = cor(e.ab$value.x,e.ab$value.y,method='kendall'),
    #                       #tlag1.corS = cor(e.ab$value.x,e.ab$value.y,method='spearman')
    #     )
    #   }
    # }, error = function(e) {
    #   print(sprintf("Error in row: %s, first: %s, %f second %s %f : %s" ,
    #                 thisrow,firstg,myMean(e.ab$value.x), secondg,myMean(e.ab$value.y),e))
    #   browser()
    })  
    
  }
  
 
  n.itm.e=na.omit(e.ts)
  
  
  stopCluster(cl)
  
  #Serial
  if(FALSE){
    with_progress({
      prgrs<-progressor(ngens*(ngens-1))
      for(i in 1:nrow(all.dt))
        #for(j in (i):nrow(all.dt)){
        for(j in (1):nrow(all.dt)){
          tryCatch({
            e.a=e.b=e.ab=NULL
            firstg = all.dt[i,'gene']
            secondg = all.dt[j,'gene']
            if(firstg==secondg) next
            e.a = mtx.lng.prtrb[ gene == firstg & pstep.lead==2,.(gene,prtrb,value,pstep.lead)] 
            e.b = mtx.lng.prtrb[ gene == secondg& pstep==2,.(gene,prtrb,value,pstep)] 
            #TODO
            e.ab = inner_join(e.a,e.b,by=c('prtrb'='prtrb','pstep.lead'='pstep'),all=FALSE)
            #todo calculate MI per perturbation(time series) and average them all
            thisrow = thisrow+1
            
            if(all(e.ab$value.y==0) | all(e.ab$value.x==0) ) next
            n.itm.e[thisrow, `:=`(first=firstg,second=secondg,  
                                  tlag1.mi = mpmi::cmi.pw(e.ab$value.x,e.ab$value.y)$mi,
                                  tlag1.corP = cor(e.ab$value.x,e.ab$value.y,method='pearson'),
                                  tlag1.corK = cor(e.ab$value.x,e.ab$value.y,method='kendall'),
                                  tlag1.corS = cor(e.ab$value.x,e.ab$value.y,method='spearman')
            )] #Calculate edge between two genes
            #n.itm.e=rbind(n.itm.e,tmp)
            prgrs()
            #if(nrow(n.itm.e[first==firstg & second==secondg,])==0 ) print(sprintf("Error in adding: %s, first: %s, second %s" ,thisrow,firstg,secondg))
            if(thisrow %% 100 == 0) print(sprintf(" %s 100 calculations done", as.integer(thisrow/100)))
          }, error = function(e) {
            print(sprintf("Error in row: %s, first: %s, %f second %s %f : %s" ,
                          thisrow,firstg,myMean(e.ab$value.x), secondg,myMean(e.ab$value.y),e))
            browser()
          })
        }
      
    })} 
  
  
}

rm(tmp.1,tmp,similarityMatrices,n3d.e,actualDataset,mtx.lng,total.rds,rds3,rds2,t1)
gc()
#save(n.itm.e,file = paste0('n.itm.e two way 0 -',studyDataset,'.robj'))
#load(file = paste0('n.itm.e two way 0 -',studyDataset,'.robj'))
n.itm.e = n.itm.e %>% dplyr::mutate(src = first , trgt=second)
n.itm.e = n.itm.e %>% dplyr::mutate(DiagnosesAb.x = "" , DiagnosesAb.y="")

all.dt=as.data.frame(all.dt)
if(!any(colnames(all.dt)=='vID')) all.dt = all.dt %>% tibble::rownames_to_column('vID')
setDF(all.dt) #lapply(all.dt,class);

n.itm.v = data.table()
n.itm.v=all.dt %>% dplyr::mutate(ocr = rowSums(.[which(!colnames(setDF(all.dt)) %in% c('vID','gene'))])) %>% 
  dplyr::select(vID,ocr)


#the sum of all genes are the same 100.9978 , that is weired and not a natural distribution

# for(i in 1:length(unique(n.itm.e$src))){
#   n.itm.v=rbind(n.itm.v,data.frame(vID =all.dt[i,'rowname'],ocr = sum(all.dt[i,2:ncol(all.dt)])))
#   
# }
n.itm.v=setDT(n.itm.v)[,`:=`(OcrInp=ocr,OcrOut=ocr)][,ocr:=NULL]

n.itm.e.back0=n.itm.e 
n.itm.v.back0=n.itm.v #tt

#n.itm.e=n.itm.e.back0;n.itm.v=n.itm.v.back0
# Preliminary data assessement ---------

#plot important features distributions
{
  library(plotly)
  library(GGally)
  
  
  # pheatmap to create & plot clusters of similarly expressed genes. So instead of plotting 30000 genes, you will be plotting x number (can be 25, 50, 100 or more) of clusters of similarly expressed genes by providing a value to k_means parameter in the pheatmap function. If you want to cluster rows, use cluster_rows=T and to cluster columns, use cluster_cols=T (you may want to do both because of the large dataset).
   {
    all.dt.dg <- DGEList(all.dt.n);all.dt.dg$samples
    
    # we can also adjust the labelling if we want
    barplot(all.dt.dg$samples$lib.size/1e5, names=colnames(all.dt.dg), las=2, ann=FALSE, cex.names=0.75)  
    
    boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2)
    # Let's add a blue horizontal line that corresponds to the median logCPM
    abline(h=median(logcounts),col="blue")
    
    
    # Get log2 counts per million
    logcounts <- cpm(all.dt.dg,log=TRUE)
    # Check distributions of samples using boxplots
    boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2)
  }
  
  # df=setDF(tst1.totalset) %>% mutate(eID = paste0('src','_', 'trgt')) %>% select(c('class1',trainingTarget, prd.varimp$variable[1:10]))  #tst1.totalset #newDataEnv$tst1.totalset # 
  # dfl = df %>% pivot_longer(prd.varimp$variable[1:10]) %>%  
  #   dplyr::mutate(logval = sign(value) * log2(abs(value))) 
  
  #https://gdevailly.netlify.app/post/plotting-big-matrices-in-r/
  # reduce matrix size, using a summarising function (default, mean)
  {
    redim_matrix <- function(
      mat,
      target_height = 100,
      target_width = 100,
      summary_func = function(x) mean(x, na.rm = TRUE),
      output_type = 0.0, #vapply style
      n_core =detectCores()-3 # parallel processing
    )
    {
      library(doFuture)
      if(target_height > nrow(mat) | target_width > ncol(mat)) {
        stop("Input matrix must be bigger than target width and height.")
      }
      
      try({
        cl=makeCluster(ncores) # 
        registerDoFuture()  #registerDoSNOW(cl)
        plan(multisession)
        
        
        seq_height <- round(seq(1, nrow(mat), length.out = target_height + 1))
        seq_width  <- round(seq(1, ncol(mat), length.out = target_width  + 1))
        
        # complicated way to write a double for loop
        redimed = do.call(rbind, parLapply(cl,seq_len(target_height), function(i,mat,seq_height,seq_width) { # i is row
          vapply(seq_len(target_width), function(j) { # j is column
            summary_func(
              as.numeric(mat[
                seq(seq_height[i], seq_height[i + 1]),
                seq(seq_width[j] , seq_width[j + 1] )
              ])
            )
          }, output_type)
        },mat,seq_height,seq_width))
      })
      stopCluster(cl)
      redimed
    }
    
    # left matrix, we take the mean of the -log10 of the p-values
    all.dt.heat1 = redim_matrix(
      t(as.matrix(all.dt))[-1,],
      target_height = 59, target_width = 600,
      summary_func = function(x) mean(x, na.rm = TRUE),
      n_core = 10
    )
    
    
    
    
    # layout(matrix(c(1, 2), nrow = 1))
    # 
    # # left plot, original matrix
    # image(
    #   t(genmat),
    #   axes = FALSE,
    #   col = colorRampPalette(c("white", "darkorange", "black"))(30),
    #   breaks = c(seq(0, 3, length.out = 30), 100),
    #   main = "Original matrix"
    # )
    # box()
    
    missedGenes = all.dt.n %>% filter(rownames(.) %in% names(frequentlyMissedGenes ))
    # right plot, summarised matrix
    image(
      t(missedGenes),
      axes = FALSE,
      col = colorRampPalette(c("white", "darkorange", "black"))(30),
      breaks = c(seq(0, 3, length.out = 30), 100),
      main = "Reduced matrix"
    )
    box()
    
  }
  
  #Heatmap
  {
    
    
    library(edgeR)
    library(limma)
    library(Glimma)
    library(org.Mm.eg.db)
    library(gplots)
    library(RColorBrewer)
    library(NMF)
    ## Get some nicer colours
    mypalette <- brewer.pal(11,"RdYlBu")
    morecols <- colorRampPalette(mypalette)
    # Set up colour vector for celltype variable
    #col.cell <- c("purple","orange")[sampleinfo$CellType]
    
    # Plot the heatmap
    heatmap.2(as.matrix(all.dt.n[6e3:1e4,]),trace="none", main="Histogram of genes 6000:10000",
              scale="row")
    
    heatmap(as.matrix(all.dt.n[6e3:1e4,]))
    # heatmap.2(all.dt.n,col=rev(morecols(50)),trace="none", main="Top 500 most variable genes across samples",
    #           ColSideColors=col.cell,scale="row")
    
    
    
    
    
    all.dt.n = all.dt %>% select(-vID) %>% mutate_all(as.numeric) 
    threshold = 1e-5
    
    all.dt.t = as.data.frame(t(all.dt.n))
    cellsMissingGenes=sapply(all.dt.t, function(x) sum(x<threshold));cellsMissingGenes[cellsMissingGenes>0]
    hist(cellsMissingGenes,main = "Historgram of cells missing genes", ylab = "Missing Gene Freq.")#   geneMIssing[geneMIssing>0]
    
    cellsMissingGenes=as.data.frame(cellsMissingGenes) %>% tibble::rownames_to_column('gene')
    plot_ly(cellsMissingGenes,x='gene', type='histogram')
    
    
    genesAbscentInCells=sapply(all.dt, function(x) sum(x<threshold));genesAbscentInCells[genesAbscentInCells>0]
    hist(genesAbscentInCells,main = "Historgram of genes abscent in cells", ylab = "Freq cells with missed genes")#   geneMIssing[geneMIssing>0]
    
    frequentlyMissedGenes = genesMissed[genesMissed>50]
    
    histogram(n.itm.e$corP)
    fltrd.srcCnt =  n.itm.e %>% filter(abs(corP) > .2 ) %>% group_by(src ) %>% summarise(srcCnt = count(.))
    fltrd.trgtCnt =  n.itm.e %>% filter(abs(corP) < .2 ) %>% group_by(trgt ) %>% summarise(trgtCnt = count(.))
    
    
    
    
    
    
    
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
  #Plotting the Iris dataset in 3D plot_ly(x=Sepal.Length,y=Sepal.Width,z=Petal.Length,type="scatter3d",mode='markers',size=Petal.Width,color=Species)
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
  df=setDF(tst1.totalset) %>% mutate(eID = paste0('src','_', 'trgt')) %>% 
    mutate(class10 = class1 =='c' , class11= class1 =='ir') %>% 
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
  ))))
}



# CICT features calculations -------

#n.itm.e=n.itm.e.back0;n.itm.v=n.itm.v.back0
n.itm.e$edgeTyp = '' #'ewMIempirical' #tlag1.mi  'ewMImm' #"tlag1.mi" #ewMIempirical #Pearson #'corP'
n.itm.e=n.itm.e %>% dplyr::mutate(Weight=ewMImm  ) %>% mutate(corP=NULL) %>% 
  #select(-first,-second) %>% 
  select(src,trgt,everything()) # mf.mi )#)tlag1.mi  #tlag1.corP

# gt1=tbl.goldStandard %>% mutate(gtweight=1) %>% dplyr::select(src,trgt,gtweight)
# n.itm.e = merge(n.itm.e,gt1,all.x= T,by=c('src','trgt'))
# hist(n.itm.e$Weight)
# n.itm.e.rmvd=setDT(n.itm.e)[I<=.15 | gtweight !=1,]
# setDT(n.itm.e)[I>.15  | gtweight ==1,Weight :=0]
# n.itm.e=n.itm.e[,gtweight:=NULL] #s0c0 L-
# hist(n.itm.e$Weight)

rm(ig)
ig=graph.data.frame(n.itm.e[n.itm.e$Weight>0.5,],directed=FALSE) #An ealy threshold lets the graph to produce more useful information

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

try({n.itm.e = dplyr::rename(n.itm.e,src = first, trgt = second)})
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

a=sapply(n.itm.v, function(x) sum(is.na(x)));a[a>0]
a=sapply(n.itm.e, function(x) sum(is.na(x)));a[a>0]

a = n.itm.v[!vID %in% n.itm.e$trgt]


vinfcols=c("vID","SumOcrInp.x","SumOcrOut.x","SumOcrInp.y","SumOcrOut.y")
setDF(n.itm.e);setDF(n.itm.v)
n.itm.e = n.itm.e %>% inner_join(n.itm.v[,vinfcols],by=c("src"="vID"))
n.itm.e = n.itm.e %>% inner_join(n.itm.v[,vinfcols],by=c("trgt"="vID"))

a=sapply(n.itm.e, function(x) sum(is.na(x)));a[a>0]

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


rm(n.itm.v.srcsum,n.itm.v.trgtsum,n.itm.v.influx,n.itm.v.outflux,n.itm.e.tmp)
# Add CICT zone interactions -----
#For each pair of nodes calculate interactions(e.g. mutual information) of their four zones
#Preparation for the 3 following stpes
{
  globalDistDiscretization = TRUE
  
  ngens = nrow(all.dt)
  nbins=as.integer(sqrt(Nobservations)*1)
  breaks =1:nbins;breaks = NULL
  #breaks = seq(0,maxEdgeWeight,maxEdgeWeight/nbins+0.0001)
  
  frchtabscisse = 1:nbins
  maxEdgeWeight = max(n.itm.e$conf,n.itm.e$contrib)
  
  myhist = function(v,nbins,breaks=NULL,plot=F,prob = T)
  {
    h=""
    if(!is.null(breaks)) h = hist( v,breaks = breaks,plot=F,probability = prob)  else h = hist( v,nbins,plot=F,probability = prob)
    h=c(h,rep(0,nbins-length(h)))
    if(prob ==T) h$density else h$counts
  }
  
  n.itm.e=n.itm.e  %>% dplyr::mutate(confdisc = NULL, contribdisc=NULL)
  n.itm.e = n.itm.e  %>% dplyr::mutate(confdisc = unlist(infotheo::discretize(na.omit(conf), "equalfreq", nbins)),
                                       contribdisc = unlist(infotheo::discretize(na.omit(contrib), "equalfreq", nbins))) 
  #a=infotheo::discretize(na.omit(n.itm.e$conf), "equalfreq", 10)$X
  e.cnfcnt= n.itm.e%>% dplyr::select(src,trgt, confdisc,contribdisc) #conf, contrib,
  
  setDT(n.itm.v);n.itm.v[,`:=`(confhist=NULL,contribhist=NULL)]
  setDT(e.cnfcnt)
  setkeyv(e.cnfcnt, c('src','trgt'))
  
  
  
}

#1- Add distribution of conf and contribs of nodes edges to each node 
thisrow = 0;j=1;i=1

for(j in (1):ngens){
  tryCatch({
    e.a=e.b=e.ab=NULL
    firstg = n.itm.v[j,]$vID
    
    #a.confdisc and a.contribdisc have all confidence and contributions that  originates from source, zone 1 and 2
    #b.confdisc and b.contribdisc have all confidence and contributions that  ends to target, zone 4 and 3
    
    if(!globalDistDiscretization){ #Discretization was applied on edges of each particular node
      a.confdisc = myhist(infotheo::discretize(e.cnfcnt[ src == firstg,]$conf, "equalwidth", nbins)$X,nbins,breaks)
      a.contribdisc =  myhist(infotheo::discretize(e.cnfcnt[ src == firstg,]$contrib, "equalwidth", nbins)$X,nbins,breaks)
      
    }else if(globalDistDiscretization){ #Discretization was applied on all edges
      suppressMessages({
        a.confhist = myhist( e.cnfcnt[ src == firstg,]$confdisc,nbins,breaks,plot=F,prob = F)
        a.contribhist =  myhist(e.cnfcnt[ src == firstg,]$contribdisc,nbins,breaks,plot=F,prob = F)
      })
      
      
    }
    
    #Calculating interactions between zone 1 and 2 of source with target's zones 3 and 4.
    #Attention, interactions between zone 1,2 of target and sources's zones 3 & 4 for this edge will be received in BA calculations.
    n.itm.v[vID==firstg, `:=`(confhist = list(a.confhist))]
    n.itm.v[vID==firstg, `:=`(contribhist = list(a.contribhist))]  
  }, error = function(e) print(sprintf("Error in row: %s, first: %s, " ,thisrow,firstg)))
}

#2- ADD src and trgt zone interactions for each edge
if(FALSE){
  e.cnfcnt = merge(e.cnfcnt,n.itm.v[,.(vID,confhist,contribhist)],by.x="src", by.y="vID")
  e.cnfcnt = merge(e.cnfcnt,n.itm.v[,.(vID,confhist,contribhist)],by.x="trgt", by.y="vID")
  
  mydistFrechet<-function(x,y,frchtabscisse,timeScale=0.1, FrechetSumOrMax = "sum"){
    lapply(1:length(x),function(i,x,y,frchtabscisse,timeScale,FrechetSumOrMax){
      kmlShape::distFrechet(frchtabscisse,x[[i]],frchtabscisse, y[[i]], timeScale=timeScale, FrechetSumOrMax = FrechetSumOrMax)
    },x,y,frchtabscisse,timeScale,FrechetSumOrMax)
  }
  
  myMutinformation<-function(x,y,method){
    lapply(1:length(x),function(i,x,y,method){infotheo::mutinformation(x[[i]],y[[i]],method = method)},
           x,y,method)
    
  }
  
  e.cnfcnt = e.cnfcnt %>% 
    dplyr::mutate(    
      cnfxcnfy.emi = myMutinformation(confhist.x,confhist.y,method = "emp"),
      cnfxcnty.emi = myMutinformation(confhist.x,contribhist.y,method = "emp"),
      cntxcnfy.emi = myMutinformation(contribhist.x,confhist.y,method = "emp"),
      cntxcnty.emi = myMutinformation(contribhist.x,contribhist.y,method = "emp"),
      
      cnfxcnfy.frcht = mydistFrechet(confhist.x,confhist.y, frchtabscisse, timeScale=0.1, FrechetSumOrMax = "sum"), #timescale is lambda
      cnfxcnty.frcht = mydistFrechet(confhist.x,contribhist.y,frchtabscisse,  timeScale=0.1, FrechetSumOrMax = "sum"),
      cntxcnfy.frcht = mydistFrechet(contribhist.x,confhist.y,frchtabscisse,  timeScale=0.1, FrechetSumOrMax = "sum"),
      cntxcnty.frcht = mydistFrechet(contribhist.x,contribhist.y,frchtabscisse,  timeScale=0.1, FrechetSumOrMax = "sum")
      # cnfxcnfy.emd = frechet::dist4den(confhist.x,confhist.y,fctn_type='density'), #emd : earth moving distance
      # cnfxcnty.emd = frechet::dist4den(confhist.x,contribhist.y,fctn_type='density'),
      # cntxcnfy.emd = frechet::dist4den(contribhist.x,confhist.y,fctn_type='density'),
      # cntxcnty.emd = frechet::dist4den(contribhist.x,contribhist.y,fctn_type='density'),
    )
  
  e.cnfcnt = e.cnfcnt %>% select(!starts_with('conf') & !starts_with('contrib'))
  n.itm.e = merge(n.itm.e,e.cnfcnt, by =c('src','trgt'))
  #n.itm.e = n.itm.e %>% dplyr::select(-confhist.x,-contribhist.x,-confhist.y,-contribhist.y)
}

#3- TODO Fit distributions and add dist. parameters to  thenodes #Not a strong feature
if(FALSE){
  #Todo fit distributions store params for each subtype
  require(fitdistrplus)
  
  myFitDist = function(x,discrete = F){
    # print(x)
    # x = unlist( n.itm.v[x,distCol,with=FALSE])
    
    x= ifelse(x==0,PSUDEO_ZERO,x)
    x=x/sum(x)
    fg <- fitdist(x, "gamma",discrete = discrete)
    fln <- fitdist(x, "lnorm",discrete = discrete)
    try({fb <- fitdist(x, "beta",discrete =discrete)})
    if(is.null(fb))try({fb <- fitdist(x, "beta",discrete =discrete,method = "mse")})
    gof = NULL
    tryCatch({
      gof = gofstat(list(fln,fg, fb ),fitnames = c("lnorm", "gamma", "beta"))
      bestAICfit = gof$aic[gof$aic==min(gof$aic)]
    }, error = function(e) {})
    if(is.null(gof)){
      try({
        gof = gofstat(list(fln,fg, fb ),fitnames = c("lnorm", "gamma", "beta"),chisqbreaks = seq(0,1,.1) )
        bestAICfit = gof$aic[gof$aic==min(gof$aic)]
      })
    }
    
    if(is.null(gof)) bestAICfit = NA 
    list(fitdist = names(bestAICfit),fitAIC=bestAICfit,
         beta1STshape = fb$estimate[1],beta2ndshape = fb$estimate[2], 
         loglinmean = fln$estimate[1],loglinsd = fln$estimate[2],
         beta=fb,gamma=fg,loglin =fln)
  }
  

}
#2@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@222

#save(n.itm.e,n.itm.v,file="n.itm.e pearson weight, discretized zone interactions and n.itm.v.robj")
# 
rm(n.itm.e.back3,n.itm.e.BeforeAugmenting,n.itm.v.BeforeMoments,n.itm.e.noselfedge)

#n.itm.e=n.itm.e.back;n.itm.v=n.itm.v.back

#   Add Harmonized transition rates HTR-----

n.itm.e = n.itm.e %>%dplyr::mutate(srctrgtSum = OcrOut.x+OcrInp.y,
                                   srctrgtProduct=OcrOut.x * OcrInp.y,
                                   HTR=(Weight*srctrgtSum)/srctrgtProduct,
                                   EE=OcrOut.x^2*OcrInp.y/srctrgtSum)#Harmonized transition rate

self=n.itm.e %>% group_by(src) %>% summarise( scfSumHTR=sum(HTR,na.rm= TRUE),scfSumEE=sum(EE,na.rm=TRUE))
others=n.itm.e %>% group_by(trgt) %>%  summarise( ocbSumHTR=sum(HTR,na.rm= TRUE),ocbSumEE=sum(EE,na.rm=TRUE))

n.itm.v.back=n.itm.v
n.itm.v = n.itm.v %>% left_join(self,by=c("vID"="src"))
n.itm.v = n.itm.v %>% left_join(others,by=c("vID"="trgt"))

a=sapply(n.itm.v, function(x) sum(is.na(x)));a[a>0]
#View(n.itm.v[,c("vID","scfSumHTR","ocbSumHTR"),with=FALSE])
setDF(n.itm.v);n.itm.v = tidyr::replace_na(n.itm.v,replace=list(scfSumHTR=0,ocbSumHTR=0,
                                                                scfSumEE=0,ocbSumEE=0,DiagnosesAb="",SumOcrInp.x=0))

n.itm.v1=n.itm.v[,c("vID","scfSumHTR","ocbSumHTR","scfSumEE","ocbSumEE")]

setDF(n.itm.e); setDF(n.itm.v1)
n.itm.e=n.itm.e %>% left_join(n.itm.v1,by=c("src"="vID")) 
n.itm.e=n.itm.e %>% left_join(n.itm.v1,by=c("trgt"="vID"))


#   Enhance edges conf contrib -----
if(!exists("PSUDEO_ZERO")) stop("define PSUDEO_ZERO")
n.itm.e = setDF(n.itm.e) %>% #dplyr::mutate(Weight=mf.mi) %>%
  dplyr::mutate(V=abs(OcrOut.x-OcrInp.y), R=V/Weight,P=V*Weight,
                
                NHTRfromSource=HTR/scfSumHTR.x,
                NHTRtoTarget=HTR/ocbSumHTR.y,
                tNHTR=NHTRfromSource+NHTRtoTarget,
                
                NEEfromSource=EE/scfSumEE.x,
                NEEtoTarget = EE/ocbSumEE.y,
                
                NEE_NHTR_Diff=NHTRfromSource-NEEfromSource,
                NEE_NHTR_Ratio=NHTRfromSource/ifelse(NEEfromSource==0,PSUDEO_ZERO,NEEfromSource),
                
                confN = conf*(OcrInp.y/srctrgtSum), #NHTRfromSource
                contribN = contrib*(OcrOut.x/srctrgtSum), #NHTRtoTarget, #
                tn=confN+contribN,t=conf+contrib,
                #+1 is necessary to protect recursive edges from going infinity
                
                Pout=(OcrOut.x-OcrInp.y)*Weight/OcrOut.x,
                Pin=(OcrOut.x-OcrInp.y)*Weight/OcrInp.y,#/maxItmstOcr
                slope = abs(OcrOut.x-OcrInp.y)/Weight, #tng
                hypotenuse = sqrt(V^2+Weight^2),
                EO =emissionRate *(srctrgtProduct)/(SumOcrOut.x.x+1),
                EI =emissionRate *(srctrgtProduct)/(SumOcrOut.x.y+1),
                OEER=Weight/EO,OERR=Weight/EI,
                toe=OEER + OERR,
                AM=(OEER+OERR)/2,GM=(OEER*OERR)^.5, HM=2/((1/OEER)+(1/OERR)),
                UR=V/(Weight-EO) ,URR= V/(Weight-EI) ) #unexpected resistance

a=sapply(n.itm.e, function(x) sum(is.na(x)));a[a>0]
b=sapply(n.itm.e, function(x) sum(x==0 & !is.na(x)) );b[b>0]
c=sapply(n.itm.e, function(x) sum(is.infinite(x)));c[c>0]
#d=setDT(n.itm.e)[is.infinite(contribZ),.(tz,odds,contribZ,OcrOut.x,OcrInp.y,I,IBA,DiagnosesAb.x,DiagnosesAb.y)];View(d)
setDT(n.itm.e);rpt=n.itm.e[Weight==0,.(src,trgt,conf,contrib,Weight)];View(rpt)
#Outdegree.x, #scaled value of expected weight
#Unexpected resistance


' %>% #(EWeight-Weight)/(Weight*EWeight)
  dplyr::mutate(conf=round(conf,6) , contrib=round(contrib,6),
  confN=round(confN,6) , contribN=round(contribN,6),
  Weight=round(Weight,6) , V=round(V,6), P=round(P,6), R=round(R,6),
  EI=round(EI,6),ER=round(ER,6),EP=round(EP,6),Xs=round(Xs,6))'
#Resistance Adjusted
#Rate = ( Scaled MAX - Scaled MIN ) / ( Input MAX - Input MIN )
#Offset = Scaled MIN - ( Input MIN * Rate )
#Scaled Value = (Input Value * Rate) + Offset

#   Add Power parameters to vertices and edges ----
#all the power input from different sources to each trgt
powerParamsIn=n.itm.e %>% group_by(trgt) %>%
  summarise( Pinsum=sum(Pin,na.rm=TRUE),PinAbsSum = sum(abs(Pin),na.rm=TRUE),
             PinSD=sd(Pin,na.rm=TRUE),PinMean=mean(Pin,na.rm=TRUE))
powerParamsOut=n.itm.e %>% group_by(src) %>%
  summarise(Poutsum=sum(Pout,na.rm=TRUE),PoutAbsSum = sum(abs(Pout),na.rm=TRUE),
            PoutSD=sd(Pout,na.rm=TRUE),PoutMean=mean(Pout,na.rm=TRUE))


a=sapply(powerParamsOut, function(x) sum(is.na(x)));a[a>0]
a=sapply(powerParamsIn, function(x) sum(is.na(x)));a[a>0]

n.itm.v = n.itm.v %>% left_join(powerParamsIn,by=c("vID"="trgt"))
n.itm.v = n.itm.v %>% left_join(powerParamsOut,by=c("vID"="src"))

n.itm.v = replace_na(n.itm.v,replace=list(PinSD=0,PoutSD=0))

vinfcols=c("vID","Pinsum","PinAbsSum","PinSD","PinMean","Poutsum","PoutAbsSum","PoutSD","PoutMean")
setDF(n.itm.e);setDF(n.itm.v)
n.itm.e = n.itm.e %>% inner_join(n.itm.v[,vinfcols],by=c("src"="vID"))
n.itm.e = n.itm.e %>% inner_join(n.itm.v[,vinfcols],by=c("trgt"="vID"))

n.itm.e.back=n.itm.e #n.itm.e=n.itm.e.back
n.itm.v.back=n.itm.v #n.itm.v=n.itm.v.back;n.itm.e=n.itm.e.back




#remove duplicated columns 
library(digest)
dupcols=c()
e.set1 = sample_n(n.itm.e,20000)
dupcols = clnames[duplicated(lapply(e.set1, digest))];dupcols
dupcols = dupcols[-c(1,3)] #[-c(1,2,4,5,7)]
n.itm.e = n.itm.e %>% select(! any_of(dupcols))

gc()


# <=== 21st December 2021 ------------------

#   EARLY attach Ground truth  ----
  if(FALSE){rm(gt1,gt1.c,gt1.rc,gt1.c.r,gt1.abscentGT,gt2.rc,gt2.rnd)
  breaks = seq(0,1,.1) #breaks = seq(0,1,10)
  gt1=tbl.goldStandard %>% dplyr::select(src,trgt)
  
  
  setDF(n.itm.e)
  #Test
  tmp =n.itm.e.back0 %>% inner_join(gt1,by=c("first"="src","second"="trgt"))
  gt1.abscentGT =gt1 %>% anti_join(n.itm.e.back0,by=c("src"="first","trgt"="second"))
  gt1.abscentGT =gt1.abscentGT %>% anti_join(n.itm.e.back0,by=c("src"="second","trgt"="first"))
  
  gt1.c = n.itm.e %>% inner_join(gt1,by=c("src"="src","trgt"="trgt"))
  gt1.c$predicate = "CAUSES"
  hist(gt1.c$confdisc, 10)
  
  gt1.rc = n.itm.e %>% inner_join(gt1.c[,c('src','trgt')],by=c("src"="trgt","trgt"="src"))
  gt1.rc$predicate = "REV_CAUSES"
  hist(gt1.rc$confdisc, 10)
  
  gt1.causal_reversecausal = rbind(gt1.c,gt1.rc)
  #tmp = rbind(gt1,gt1.c,gt1.rc)
  
  #Adding 2000 random edges
  gt1.rnd = n.itm.e %>% 
    anti_join(gt1.causal_reversecausal,by=c("src"="src","trgt"="trgt")) %>%
    mutate(predicate = "IRRELEVANTA")
  
  #tmp.rdndnt = setDT(gt1.rnd)[duplicated(gt1.rnd[,.(src,trgt)]),]
  
  gt1.rnd$predicate = "IRRELEVANTA"
  gt1.rnd.smpl = gt1.rnd %>% dplyr::slice_sample(n=2000)
  
  gt2=rbind(gt1.causal_reversecausal,gt1.rnd)
  
  gt2= as.data.frame(gt2) %>% dplyr::select(edgeTyp ,src, trgt,Weight,,
                                            OcrOut.x, OcrInp.y,#intvl_avg,intvl_sd,intvl_median,
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
  
  colnames(n.itm.e.back0)
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
  
  
  
  tmp.e=n.itm.e %>% dplyr::select(src,trgt,OcrOut.x,OcrInp.y,Weight,SumOcrOut.x.x,SumOcrInp.y.x,Weight,V,R,EO,EI,OEER,OERR,AM,GM,HM,UR)
  
  summary(n.itm.e$EE);summary(n.itm.e$NEEfromSource);summary(n.itm.e$NEEtoTarget)
  summary(n.itm.e$UR);summary(n.itm.e$EO);summary(n.itm.e$EI); summary(n.itm.e$OERR);summary(n.itm.e$OEER)
  summary(n.itm.e$AM);summary(n.itm.e$HM);summary(n.itm.e$GM)
  summary(n.itm.e$Pin);summary(n.itm.e$Pout)
  summary(n.itm.e$EI)
  
  #summary(n.itm.e$SumOcrInp.y.x);summary(n.itm.e$R);summary(n.itm.e$ER);summary(n.itm.e$RA)
  ggplot(n.itm.e, aes(EE))+ggtitle("EE") + geom_histogram()  + scale_x_log10()
  ggplot(n.itm.e, aes(NEEfromSource))+ggtitle("NEEfromSource") + geom_histogram()  + scale_x_log10()
  ggplot(n.itm.e, aes(NHTRfromSource))+ggtitle("NHTRfromSource") + geom_histogram()  + scale_x_log10()
  ggplot(n.itm.e, aes(NEE_NHTR_Diff))+ggtitle("NEE_NHTR_Diff") + geom_histogram()  + scale_x_log10()
  ggplot(n.itm.e, aes(NEEtoTarget))+ggtitle("NEEtoTarget") + geom_histogram()  + scale_x_log10()
  ggplot(n.itm.e, aes(NHTRtoTarget))+ggtitle("NHTRtoTarget") + geom_histogram()  + scale_x_log10()
  ggplot(n.itm.e, aes(tNHTR))+ggtitle("tNHTR") + geom_histogram()  + scale_x_log10()
  
  ggplot(n.itm.e, aes(cbMean.x))+ggtitle("cbMean.x - log scale") + geom_histogram() +
    scale_x_log10()+geom_vline(xintercept=.006)
  ggplot(n.itm.e, aes(OEER))+ggtitle("OEER - log scale")  +   scale_x_log10()  +geom_vline(xintercept=.42)
  ggplot(n.itm.e, aes(EO))+ggtitle("EO - log scale")  +   scale_x_log10()
  ggplot(n.itm.e, aes(EI))+ggtitle("EI - log scale") +geom_histogram()   +   scale_x_log10()
  ggplot(n.itm.e, aes(cfNMad.x))+ggtitle("cfNMad.x - log scale")  +   scale_x_log10()  +geom_vline(xintercept=.00003)
  ggplot(n.itm.e, aes(UR))+ggtitle("UR") + geom_histogram() + scale_x_log10()+geom_vline(xintercept=0.0003)
  ggplot(n.itm.e, aes(tn))+ggtitle("tn") + geom_histogram() + scale_x_log10()
  ggplot(n.itm.e, aes(t))+ggtitle("t") + geom_histogram() + scale_x_log10()
  ggplot(n.itm.e, aes(AM))+ggtitle("AM") + geom_histogram() + scale_x_log10()
  ggplot(n.itm.e, aes(confN))+ggtitle("confN") + geom_histogram() + scale_x_log10()
  ggplot(n.itm.e, aes(Pin))+ggtitle("Pin") + geom_histogram() + scale_x_log10()
  ggplot(n.itm.e, aes(Pout))+ggtitle("Pout") + geom_histogram() + scale_x_log10()
  ggplot(n.itm.e, aes(RABADiff))+ggtitle("RABADiff") + geom_histogram() + scale_x_log10()
}
#    Add BA factors and combined measures -----
if(TRUE){ 
  setDF(n.itm.e)
  #n.itm.e<-subset(n.itm.e,select=-c(confBA,contribBA,confBA1,contribBA1))
  gc()
  n.itm.e.back=n.itm.e
  #n.itm.e1= n.itm.e %>% dplyr::select(src,trgt,IBA=Weight,URBA=UR,AMBA=AM,OEERBA=OEER,OERRBA=OERR,
  # confBA=conf,contribBA=contrib,confNBA=confN,contribNBA=contribN,
  # tBA=t,tnBA=tn,toeBA=toe,HTRBA=HTR,tNHTRBA=tNHTR,
  # PinBA=Pin,PoutBA=Pout,
  # NHTRfromSourceBA=NHTRfromSource,NHTRtoTargetBA=NHTRtoTarget,
  # EEBA=EE,NEEfromSourceBA=NEEfromSource,NEEtoTargetBA=NEEtoTarget)
  # 
  n.itm.e1 = n.itm.e   
  n.itm.e1=plyr::rename(n.itm.e, replace=c( 'Weight' ='IBA' ,'UR' ='URBA' ,'AM' ='AMBA' ,'OEER' ='OEERBA' ,'OERR' ='OERRBA' ,
                                            'conf' ='confBA' ,'contrib' = 'contribBA','confN' ='confNBA' ,'contribN' ='contribNBA' ,
                                            't' ='tBA' ,'tn' ='tnBA' ,'toe' ='toeBA' ,'HTR' = 'HTRBA','tNHTR' ='tNHTRBA' ,
                                            'Pin' = 'PinBA','Pout' = 'PoutBA',
                                            'NHTRfromSource' = 'NHTRfromSourceBA','NHTRtoTarget' ='NHTRtoTargetBA' ,
                                            'EE' ='EEBA' ,'NEEfromSource' ='NEEfromSourceBA' ,'NEEtoTarget'='NEEtoTargetBA'), warn_missing = F)   %>% 
    
    select(any_of(c ('src','trgt','IBA','URBA','AMBA','OEERBA','OERRBA',
                     'confBA','contribBA','confNBA','contribNBA',
                     'tBA','tnBA','toeBA','HTRBA','tNHTRBA',
                     'PinBA','PoutBA',
                     'NHTRfromSourceBA','NHTRtoTargetBA',
                     'EEBA','NEEfromSourceBA','NEEtoTargetBA')))
  
  # cnfxcnfy.emiBA = cnfxcnfy.emi, cnfxcnty.emiBA= cnfxcnty.emi, 
  # cntxcnfy.emiBA= cntxcnfy.emi  ,cntxcnty.emiBA =cntxcnty.emi,
  # cnfxcnfy.frchtBA = cnfxcnfy.frcht, cnfxcnty.frchtBA= cnfxcnty.frcht, 
  # cntxcnfy.frchtBA= cntxcnfy.frcht  ,cntxcnty.frchtBA =cntxcnty.frcht,
  
  
  #NHTRtoTargetBA=HTRBA/ocbSumHTR.x , HTRBA/scbSumHTR.y
  
  #"Pinsum","PinAbsSum","PinSD","PinMean","Poutsum","PoutAbsSum","PoutSD","PoutMean"
  n.itm.e2 = n.itm.e %>% left_join(n.itm.e1, by=c("src"="trgt","trgt"="src"))
  rm(n.itm.e1)
  #checking 630925 edges nonexisting edges results in 630925 conBA=0
  #tmp=n.itm.e%>% dplyr::filter(confBA==0) %>% group_by(DiagnosesAb.y) %>% summarize(n=n())
  a=sapply(n.itm.e2, function(x) sum(is.na(x)));a[a>0]
  b=sapply(n.itm.e2, function(x) sum(x==0 & !is.na(x)));b[b>0]
  
  #n.itm.e2$IBA = ifelse(is.na(n.itm.e2$IBA),PSUDEO_ZERO,n.itm.e2$IBA)
  #Bounds finding a minimum to a specific boundary
  PSUDEO_ZERO_1=1e-5
  tBAreplacement=min(myMin(n.itm.e2$tBA),PSUDEO_ZERO_1)
  if(tBAreplacement==0)tBAreplacement=PSUDEO_ZERO_1
  tnBAreplacement=min(myMin(n.itm.e2$tnBA),PSUDEO_ZERO_1)
  if(tnBAreplacement==0)tnBAreplacement=PSUDEO_ZERO_1
  toeBAreplacement=min(myMin(n.itm.e2$toeBA),PSUDEO_ZERO_1)
  if(toeBAreplacement==0)toeBAreplacement=PSUDEO_ZERO_1
  confNBAreplacement=min(myMin(n.itm.e2$confNBA),PSUDEO_ZERO_1)
  if(confNBAreplacement==0)confNBAreplacement=PSUDEO_ZERO_1
  contribNBAreplacement=min(myMin(n.itm.e2$contribNBA),PSUDEO_ZERO_1)
  if(contribNBAreplacement==0)contribNBAreplacement=PSUDEO_ZERO_1
  tnBAreplacement=min(myMin(n.itm.e2$tnBA),PSUDEO_ZERO_1)
  if(tnBAreplacement==0)tnBAreplacement=PSUDEO_ZERO_1
  tBAreplacement=min(myMin(n.itm.e2$tBA),PSUDEO_ZERO_1)
  if(tBAreplacement==0)tBAreplacement=PSUDEO_ZERO_1
  toereplacement=min(myMin(n.itm.e2$toeBA),PSUDEO_ZERO_1)
  if(toeBAreplacement==0)toeBAreplacement=PSUDEO_ZERO_1
  
  PSUDEO_ZERO_2 = 0.0001
  
  n.itm.e2 = replace_na(n.itm.e2,replace=
                          list(IBA=PSUDEO_ZERO_2,URBA=PSUDEO_ZERO_2,
                               AMBA=0,OEERBA=0,OERRBA=0,
                               confBA=0,contribBA=0,confNBA=0,contribNBA=0,
                               toeBA=toeBAreplacement,PinBA=0,PoutBA=0,HTRBA=0,
                               tBA=tBAreplacement,tnBA=tnBAreplacement,
                               EEBA=0,NEEfromSourceBA=0,NEEtoTargetBA=0,
                               NHTRfromSourceBA=0,NHTRtoTargetBA=0,tNHTRBA=0 ,HTRBA=0,
                               Pinsum.x=0 ,PinAbsSum.x=0,PinMean.x=0,Poutsum.y=0,
                               PoutAbsSum.y=0,PoutMean.y=0) )
  
  a=sapply(n.itm.e2, function(x) sum(is.na(x)));a[a>0]
  #do.call(data.frame,lapply(df, function(x) if(x==0),replacement,val)
  #TODO check if this is valid or not
  #n.itm.e2$toeBA <-ifelse(n.itm.e2$toeBA==0,toeBAreplacement,n.itm.e2$toeBA)
  rm(n.itm.e);gc()
  n.itm.e= n.itm.e2 %>%dplyr::mutate(confDiff=conf-confBA,contribDiff=contrib-contribBA,
                                     confNDiff=confN-confNBA,contribNDiff=contribN-contribNBA,
                                     directionLogRatio=log((Weight/OcrOut.x)/(IBA/OcrInp.y)),
                                     Idiff = Weight - IBA, IRatio = Weight/IBA,
                                     URdiff = UR-URBA,
                                     AMdiff=AM-AMBA,
                                     tdiff= t-tBA,tRatio = t/tBA,
                                     tndiff= tn-tnBA,tnRatio= tn/tnBA,
                                     toediff=toe-toeBA,toeRatio = toe/toeBA,
                                     #observed to expected emission  minus observed to expected reception for each node
                                     deltaAttitude=OEER-OERRBA,
                                     Tendency=ifelse(is.na(UR)| UR==0,0,1/UR)-ifelse(is.na(URBA)| URBA==0,0,1/URBA))
  
  a=sapply(n.itm.e, function(x) sum(is.na(x)));a[a>0]
  c=sapply(n.itm.e, function(x) sum(is.infinite(x)));c[c>0]
  print("X==0 ------------------------------------>")
  b=sapply(n.itm.e, function(x) sum(x==0 & !is.na(x)));b[b>0]
}
rm(n.itm.e1,n.itm.e2,n.itm.v1)
#View(n.itm.e %>% filter(abs(Idiff)>.1)%>% arrange(desc(Idiff)) %>% select(src,trgt,Weight,IBA, everything()))
#n.itm.e = replace_na(n.itm.e,replace=
#                          list(confDiff=0,confNDiff=0,directionLogRatio=0))


#   check histogram of specific Edges or parameters -----
if(FALSE){ 
  vn=nrow(n.itm.v) ; print(paste0("Total possible edges: ",vn/2*(vn-1)))
  sapply(n.itm.e, function(x) sum(is.na(x)))
  
  ggplot(n.itm.e, aes(confNDiff))+ggtitle("confNDiff") + geom_histogram() + scale_x_log10()
  ggplot(n.itm.e, aes(contribN))+ggtitle("contribN") + geom_histogram() + scale_x_log10()
  ggplot(n.itm.e, aes(confN))+ggtitle("confN") + geom_histogram() + scale_x_log10()
  ggplot(n.itm.e, aes(contribN))+ggtitle("contribN") + geom_histogram() + scale_x_log10()
  ggplot(n.itm.e, aes(PDiff))+ggtitle("PDiff") + geom_histogram() + scale_x_log10()
  ggplot(n.itm.e, aes(Directional))+ggtitle("Directional") + geom_histogram() + scale_x_log10()
  ggplot(n.itm.e, aes(AMdiff))+ggtitle("AMdiff") + geom_histogram() + scale_x_log10()
  ggplot(n.itm.e, aes(URdiff))+ggtitle("URdiff") + geom_histogram() + scale_x_log10()
  ggplot(n.itm.e, aes(tdiff))+ggtitle("tdiff") + geom_histogram() + scale_x_log10()
  ggplot(n.itm.e[!is.na(n.itm.e$tRatio),], aes(tRatio))+ggtitle("tRatio") + geom_histogram() + scale_x_log10()
  summary(n.itm.e$AMdiff);summary(n.itm.e$URdiff);summary(n.itm.e$deltaAttitude)
  
  summary(n.itm.e$conf)
  
  ggplot2.histogram(data=n.itm.e$conf[n.itm.e$conf<.002],
                    legendPosition="top", xName='', #groupName='',
                    alpha=0.5, fill="cyan",addDensity=TRUE,
                    addMeanLine=TRUE, meanLineColor="white", meanLineSize=1.5)
  
  summary(n.itm.e$confN)
  ggplot2.histogram(data=n.itm.e$conf[n.itm.e$confN<.002],
                    legendPosition="top", xName='', #groupName='',
                    alpha=0.5, fill="cyan",addDensity=TRUE,
                    addMeanLine=TRUE, meanLineColor="white", meanLineSize=1.5)
  
  summary(n.itm.e$conf)
  histogram(n.itm.e$conf[n.itm.e$conf<.005])
}
#View(n.itm.e[str_detect(shared_name,"2860") & str_detect(shared_name,"209"),.(shared_name,causal2,causal11)])
#View(n.itm.e[,c("shared_name","OcrOut.x","OcrInp.y","Weight","conf","contribBA","confBA","contrib")])
# View(n.itm.e[n.itm.e$src  %in% c(12676,9355) & n.itm.e$trgt %in% c(12676,9355),c("shared_name","OcrOut.x","OcrInp.y","Weight","conf","contribBA","confBA","contrib")])

#   Add transparency -----
n.itm.e=as.data.frame(n.itm.e)
n.notinVertices= n.itm.e %>% anti_join(n.itm.v,by=c("src"="vID"))
print(paste0("!!! nodes:" , paste(unique(n.notinVertices$src),collapse=","), "  do not exist in vertices"))

n.itm.e = n.itm.e %>%dplyr::mutate(trnsparency= conf/max(conf))
#tmp.e=n.itm.e  %>%  select(SUID,shared_name,Weight,src,trgt,OcrOut.x,OcrInp.y,conf,contrib)#%>% filter(src!=trgt)

#   Enhance vertices influx outflux----
n.itm.e=n.itm.e[, unique(colnames(n.itm.e))]
n.itm.v=n.itm.v[, unique(colnames(n.itm.v))] #removing duplicated columns
setDF(n.itm.e)

n.itm.v.back=n.itm.v
n.itm.e.tmp=
  #n.itm.v=n.itm.v.back
  n.itm.v.outflux = n.itm.v %>% inner_join(n.itm.e[,c("src","Weight")],by=c("vID"="src")) %>%
  group_by(vID) %>% summarise(outflux=sum(Weight))

n.itm.v.influx = n.itm.v %>% inner_join(n.itm.e[,c("trgt","Weight")],by=c("vID"="trgt")) %>%
  group_by(vID) %>% summarise(influx=sum(Weight))

n.itm.v= n.itm.v %>% left_join(n.itm.v.outflux,by=c("vID"="vID"))
n.itm.v= n.itm.v %>% left_join(n.itm.v.influx,by=c("vID"="vID"))


a=sapply(n.itm.v, function(x) sum(is.na(x)));a[a>0] # missing outflux means no target for some nodes and
#No output results in  missing outflux and contribution parameters
#No input results in  missing influx and contribution parameters
#New Added
n.itm.v=as.data.frame(n.itm.v)
n.itm.v = replace_na(n.itm.v,replace=list(Pinsum=0,
                                          PinAbsSum=0,PinSD=0,PinMean=0,
                                          Poutsum=0,PoutAbsSum=0,PoutSD=0,PoutMean=0,
                                          outflux=0,influx=0) )

'n.itm.v = n.itm.v %>% rowwise()%>%
  dplyr::mutate(outflux=sum(n.itm.e[n.itm.e$src==vID,"Weight"],na.rm=TRUE),
  influx =sum(n.itm.e[which(n.itm.e$trgt==vID),"Weight"],na.rm=TRUE))%>%
  dplyr::mutate(deltaFlow=(influx-outflux-DISPUNIFORM20*Ocr)/Ocr )'


dupcols = removeDups(n.itm.v,excluded = c('Indegree','Outdegree','','probOut','probIn'))
n.itm.v = n.itm.v%>% select(! any_of(dupcols))


rm(n.itm.v.influx,n.itm.v.outflux,n.itm.v1,powerParamsIn,powerParamsOut,n.itm.e2,n.itm.e1,n.itm.e.tmp)
#   Enhance vertices Add moments of  contribution ----

# Calculate Node total contrib and node total recieved confidences
#No output edge results in  missing contribution parameters
#Lile SD , the MAD is the median of the absolute values of the residuals (deviations) from the data's median. in Madconst we use cbMADconst= Q(0.75) is the 0.75 quantile of that underlying distribution.(Huber, 1981).#use sn::qsc instead qsc(p, xi = 0, omega = 1, alpha = 0, dp = NULL)
#xi vector of location parameters.  = mean
#omega vector of (positive) scale parameters. = SD
#alpha vector of slant parameters = Skewness

#higher kurtosis means more of the variance is the result of infrequent extreme deviations, as opposed to frequent modestly sized deviations.The kurtosiof any univariate normal distribution is 3. platykurtic<3 , leptoKurtic >>3
library(fitdistrplus)
require(moments)
require(fGarch) #install.packages("fGarch")
require(sn)
require(Lmoments)
setDT(n.itm.e);n.itm.e.noselfedge=n.itm.e [src!=trgt,.(src,contrib,conf,trgt)]#removing self edges from the calculation
setDF(n.itm.e);setDF(n.itm.v)

#selfContribNLMoments=n.itm.e %>% dplyr::group_by(src) %>% do(extractLmoments(.$contrib))
selfContribLMoments=n.itm.e %>% dplyr::group_by(src) %>% do(extractLmoments(.$contribN))
selfContribs=n.itm.e %>% group_by(src) %>% summarise(
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

selfConfLMoments=n.itm.e %>% dplyr::group_by(src) %>% do(extractLmoments(.$confN))
selfConf=n.itm.e %>% group_by(src) %>%  summarise(
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

n.itm.v.selfparams = selfContribsFull %>% inner_join(selfConfFull,by=c("scbsrc"="scfsrc"))

a=sapply(n.itm.v.selfparams, function(x) sum(is.na(x)));a[a>0]


#   Testing  examples of effect,random and cause ----
t2.tmp1=unique(t2) %>% dplyr::select(src,DiagnosesAb.x,trgt,DiagnosesAb.y ,predicate,Weight,Indegree.x,Outdegree.y) %>%
  #dplyr::filter(str_detect(predicate,"CAUSES")) %>%
  group_by(src,predicate) %>% dplyr::summarize(mxInDeg.x=max(Indegree.x),mxOutDeg.y=max(Outdegree.y),count=n(),sum=sum(Weight), averageWeight=(sum(Weight)/n()))

setDT(n.itm.e)
e=n.itm.e
e.exmpl.c=e[e$src=="G5",]; #Rheumatoid Artritis
e.exmpl.c$exmpTyp="CAUSE";e.exmpl.c$exmpDistrib=NA;e.exmpl.c$exmpEvent=e.exmpl.c$src
ggplot(e.exmpl.c, aes(conf))+ggtitle("contribZ") + geom_histogram() + scale_x_log10()

e.exmpl.r=e[e$src=="G2",]; #Pnemonia
e.exmpl.r$exmpTyp="RANDOM" ;e.exmpl.r$exmpDistrib=NA;e.exmpl.r$exmpEvent =e.exmpl.r$src
ggplot(e.exmpl.r, aes(tz))+ggtitle("Random contrib") + geom_histogram() + scale_x_log10()

e.exmpl.e=e[e$src=="G9",]; #Syncope and collapse has 8 causes and  caused nothing
e.exmpl.e$exmpTyp="EFFECT";e.exmpl.e$exmpDistrib=NA;e.exmpl.e$exmpEvent=e.exmpl.e$src
ggplot(e.exmpl.e, aes(tz))+ggtitle("EFFECT contribBA") + geom_histogram() + scale_x_log10()


e.exmpl=rbind(e.exmpl.c,e.exmpl.e,e.exmpl.r)%>% dplyr::select(DiagnosesAb.x,DiagnosesAb.y,Weight,exmpTyp,everything())
setDT(e.exmpl)
e.exmpl.measure.vars=c("conf","contrib","confBA","contribBA")
e.exmpl.histSmpl=e.exmpl[,c("exmpTyp","exmpEvent","DiagnosesAb.x","DiagnosesAb.y", e.exmpl.measure.vars),with=FALSE]

e.exmpl.histSmpl

e.exmpl.histSmpl.mlt=melt(e.exmpl.histSmpl,id.vars=c("exmpTyp","exmpEvent"),
                          measure.vars=e.exmpl.measure.vars)

facetLabeller <- function( string){
  lapply(levels(string$variable), function(l)switch(l,
                                                    conf="CONF kj",
                                                    contrib="CONTRIB kj",
                                                    confBA="CONF jk",
                                                    contribBA="CONTRIB jk"))
  
}

clrs <- c("causalHist"="red","randomHist"="blue","tmp"="yellow")
g=ggplot(data=e.exmpl.histSmpl.mlt,aes(x=value,fill=exmpTyp))+
  #facet_wrap(c("CONF kj","CONTRIB kj","CONF jk","CONTRIB jk"),ncol=4,scale="free_x")+
  facet_wrap(~variable,ncol=4,scale="free_x")+#,labeller=facetLabeller)+ #,y=..density..
  #geom_histogram(binwidth=.1, alpha=.3, position="identity")+
  geom_density(alpha=.4)+
  scale_y_continuous( trans="log1p", expand=c(0,0.1))+#scale_y_log10()+
  scale_x_log10()+
  scale_fill_manual( values  = c("red","yellow","blue"))+
  #labels=c("Cause: Rheumatoid Artritis  ","Effect: Syncope and collapse  ","Random: Pnemonia"  ))
  guides(fill=guide_legend(title=NULL) ,color = guide_legend(title=NULL)) # override.aes = list(shape=25)


gtheme= theme(strip.text.x = element_text(size = 22),
              legend.text=element_text(size=22),
              axis.text=element_text(size=22),
              axis.title=element_text(size=22),#,face="bold")
              legend.position="bottom",legend.box = "horizontal")#,legend.title=element_blank()

gmedianLines=
  geom_vline(data = ddply(e.exmpl.histSmpl.mlt, c("variable","exmpTyp"),plyr::summarize, wavg = median(value)),
             aes(xintercept=wavg,colour=exmpTyp,guides="none"), #colour=factor(class13,labels=c("Causal","Random"))
             linetype="dashed", size=1,show.legend=FALSE)

g+gtheme+gmedianLines #+gpoints






# TODO Enhance vertices Add moments of Confidence and contribution

# <=== 4th Jan 2022
#Add transparency -----
  n.itm.e=as.data.frame(n.itm.e)
n.notinVertices= n.itm.e %>% anti_join(n.itm.v,by=c("src"="vID"))
print(paste0("!!! nodes:" , paste(unique(n.notinVertices$src),collapse=","), "  do not exist in vertices"))

n.itm.e = n.itm.e %>%dplyr::mutate(trnsparency= conf/max(conf))
#tmp.e=n.itm.e  %>%  select(SUID,shared_name,Weight,src,trgt,OcrOut.x,OcrInp.y,conf,contrib)#%>% filter(src!=trgt)

#   Enhance vertices influx outflux----
n.itm.e=n.itm.e[, unique(colnames(n.itm.e))]
n.itm.v=n.itm.v[, unique(colnames(n.itm.v))] #removing duplicated columns
setDF(n.itm.e)

n.itm.v.back=n.itm.v
n.itm.e.tmp=
  #n.itm.v=n.itm.v.back
  n.itm.v.outflux = n.itm.v %>% inner_join(n.itm.e[,c("src","Weight")],by=c("vID"="src")) %>%
  group_by(vID) %>% summarise(outflux=sum(Weight))

n.itm.v.influx = n.itm.v %>% inner_join(n.itm.e[,c("trgt","Weight")],by=c("vID"="trgt")) %>%
  group_by(vID) %>% summarise(influx=sum(Weight))

n.itm.v= n.itm.v %>% left_join(n.itm.v.outflux,by=c("vID"="vID"))
n.itm.v= n.itm.v %>% left_join(n.itm.v.influx,by=c("vID"="vID"))


a=sapply(n.itm.v, function(x) sum(is.na(x)));a[a>0] # missing outflux means no target for some nodes and
#No output results in  missing outflux and contribution parameters
#No input results in  missing influx and contribution parameters
#New Added
n.itm.v=as.data.frame(n.itm.v)
n.itm.v = replace_na(n.itm.v,replace=list(Pinsum=0,
                                          PinAbsSum=0,PinSD=0,PinMean=0,
                                          Poutsum=0,PoutAbsSum=0,PoutSD=0,PoutMean=0,
                                          outflux=0,influx=0) )

'n.itm.v = n.itm.v %>% rowwise()%>%
  dplyr::mutate(outflux=sum(n.itm.e[n.itm.e$src==vID,"Weight"],na.rm=TRUE),
  influx =sum(n.itm.e[which(n.itm.e$trgt==vID),"Weight"],na.rm=TRUE))%>%
  dplyr::mutate(deltaFlow=(influx-outflux-DISPUNIFORM20*Ocr)/Ocr )'


dupcols = removeDups(n.itm.v,excluded = c('Indegree','Outdegree','','probOut','probIn'))
n.itm.v = n.itm.v%>% select(! any_of(dupcols))


rm(n.itm.v.influx,n.itm.v.outflux,n.itm.v1,powerParamsIn,powerParamsOut,n.itm.e2,n.itm.e1,n.itm.e.tmp)
#   Enhance vertices Add moments of  contribution ----

# Calculate Node total contrib and node total recieved confidences
#No output edge results in  missing contribution parameters
#Lile SD , the MAD is the median of the absolute values of the residuals (deviations) from the data's median. in Madconst we use cbMADconst= Q(0.75) is the 0.75 quantile of that underlying distribution.(Huber, 1981).#use sn::qsc instead qsc(p, xi = 0, omega = 1, alpha = 0, dp = NULL)
#xi vector of location parameters.  = mean
#omega vector of (positive) scale parameters. = SD
#alpha vector of slant parameters = Skewness

#higher kurtosis means more of the variance is the result of infrequent extreme deviations, as opposed to frequent modestly sized deviations.The kurtosiof any univariate normal distribution is 3. platykurtic<3 , leptoKurtic >>3
library(fitdistrplus)
require(moments)
require(fGarch) #install.packages("fGarch")
require(sn)
require(Lmoments)
setDT(n.itm.e);n.itm.e.noselfedge=n.itm.e [src!=trgt,.(src,contrib,conf,trgt)]#removing self edges from the calculation
setDF(n.itm.e);setDF(n.itm.v)

#selfContribNLMoments=n.itm.e %>% dplyr::group_by(src) %>% do(extractLmoments(.$contrib))
selfContribLMoments=n.itm.e %>% dplyr::group_by(src) %>% do(extractLmoments(.$contribN))
selfContribs=n.itm.e %>% group_by(src) %>% summarise(
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

selfConfLMoments=n.itm.e %>% dplyr::group_by(src) %>% do(extractLmoments(.$confN))
selfConf=n.itm.e %>% group_by(src) %>%  summarise(
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

n.itm.v.selfparams = selfContribsFull %>% inner_join(selfConfFull,by=c("scbsrc"="scfsrc"))

a=sapply(n.itm.v.selfparams, function(x) sum(is.na(x)));a[a>0]

#   Enhance vertices Add moments of Confidence and contribution ----
T<-TRUE
othersConfLMoments=n.itm.e %>% dplyr::group_by(trgt) %>% do(extractLmoments(.$confN))
othersConfs=n.itm.e %>% group_by(trgt) %>% summarise(
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

othersContribLMoments=n.itm.e %>% dplyr::group_by(trgt) %>% do(extractLmoments(.$contribN))
othersContribs=n.itm.e %>% group_by(trgt) %>%  summarise(
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

n.itm.v.othersparams = othersConfsFull %>% inner_join(othersContribsFull,by=c("ocftrgt"="ocbtrgt"))

a=sapply(n.itm.v.othersparams, function(x) sum(is.na(x)));a[a>0]




 # 11th of Jan 2020================
  #   Integrating moments with n.itm.v -----
a=sapply(n.itm.v, function(x) sum(is.na(x)));a[a>0]
n.itm.v.BeforeMoments=n.itm.v
n.itm.v = n.itm.v %>% left_join(n.itm.v.othersparams,by=c("vID"="ocftrgt"))
n.itm.v = n.itm.v %>% left_join(n.itm.v.selfparams,by=c("vID"="scbsrc"))
#Attention, lots of missings will be introduced due to
#    nodes who does not have a target like death (outflux is missing)
# nodes who does not have a source (influx is missing)

a=sapply(n.itm.v, function(x) sum(is.na(x)));a[a>0]

#replaces all columns that their name ends with given pattern
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

#This operation turns every missing value to its defualt assuming
#that if there is no information it should be normal or 0
#Values for L moments are defined based on the uniform distribution for random events
#save(n.itm.v,n.itm.e,file="vertxEdgeBackup1.robj");load(file="vertxEdgeBackup1.robj")
n.itm.v=my_replace_na(n.itm.v,rplist = list(Mean=0,Median=0,Skew=0,Kurt=3,
                                            SD=0,Total=0,NSD=0,flux=0,
                                            sum.x=0,sum.y=0,Pinsum=0,PinBAsum=0,Poutsum=0,PoutBAsum=0,
                                            MADcost=1.488,L1=0,L2=0,tau3=0,tau4=0 )) #tau4=0.1226




rm(n.itm.v.othersparams,n.itm.v.selfparams,othersConfs,othersConfsFull,othersContribs,othersContribsFull,othersConfLMoments,othersContribLMoments)
rm(selfConf,selfConfFull,selfContribs,selfContribsFull,selfConfLMoments,selfContribLMoments)

#   Add vertice level 2 distribution parameters ----  
#TODO probably wrong should use edges not vertices to find those distribution params
setDF(n.itm.v)
L2 = n.itm.v %>% ungroup() %>% summarise(
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
#TODO: add more L2 calculations to n.itm.v
n.itm.v = n.itm.v %>%dplyr::mutate(
  ZOcfSkew= (ocfSkew-L2$MedianOcfSkew)/L2$MADOcfSkew,
  ZOcfNMedian= (Poutsum-L2$MedianOcfNMedian)/L2$MADOcfNMedian,
  ZScbNSkew= (Poutsum-L2$MedianScbNSkew)/L2$MADScbNSkew,
  ZScbScbL3= (scbL3-L2$MedianScbL3)/L2$MADScbL3,
  ZPoutsum= (Poutsum-L2$MedianPoutsum)/L2$MADPoutsum )

a=sapply(n.itm.v, function(x) sum(is.na(x)));a[a>0]
#Todo verify that replacement is correct
# Mean and Median should be replaced with  !?
n.itm.v=my_replace_na(n.itm.v,rplist = list(Mean=0,Median=0,Skew=0,Kurt=3,MAD=0,
                                            SD=0,Total=0,NSD=0,flux=0,sum.x=0,sum.y=0,
                                            MADcost=1.488,L1=0,L2=0,L3=0,L4=0,tau3=0,tau4=0 ))


#   Descriptive analysis on n.itm.v to identify important variables ----
require(Hmisc)
if(FALSE){
  sink("nItmVReport.html")
  #t=describe(d)
  html(describe(n.itm.v), size=85, tabular=TRUE,
       greek=TRUE, scroll=FALSE, rows=25, cols=20)
  sink()
  browseURL("nItmVReport.html")
}
#   Add source and target charectristics to edges -----
n.itm.e.BeforeAugmenting=n.itm.e
#n.itm.v = n.itm.v %>% dplyr::mutate(DiagnosesAb.y=NULL) %>% rename(DiagnosesAb = DiagnosesAb.x)
collist= c('MAD','Median','MADconst','Mean','SD','Skew','Kurt','Total',
           'NMAD','NMedian','NMADconst','NMean','NSD','NSkew','NKurt','NTotal',
           'tnMedian','tnMADconst','tnMAD',
           'L1','L2','L3','L4','tau3','tau4')
collist=c(paste0("scb",collist),paste0("scf",collist),paste0("ocb",collist),paste0("ocf",collist))
colList=c('vID','DiagnosesAb','evenTyp','Pinsum','PinBAsum','Poutsum','PoutBAsum',
          'ZOcfSkew','ZPoutsum','ZOcfNMedian','ZScbNSkew','ZScbScbL3',collist)

#print(colList[!colList %in% colnames(n.itm.v)])
finalCoolList = colList[colList %in% colnames(n.itm.v)]
n.itm.v1=n.itm.v[,finalCoolList]


a=sapply(n.itm.v1, function(x) sum(is.na(x)));a[a>0]

setDF(n.itm.e); setDF(n.itm.v1)

#!!!NOTE that parameters with .x point to source parameters and those with .y point to target parameters
if(FALSE) n.itm.e = n.itm.e %>% dplyr::filter(abs(Weight)>.2);gc()
n.itm.e=as.data.frame(n.itm.e)
n.itm.e=n.itm.e %>% left_join(n.itm.v1,by=c("src"="vID")) 
n.itm.e = n.itm.e %>%   left_join(n.itm.v1,by=c("trgt"="vID"))

dupcols = removeDups(n.itm.e,excluded = c('OcrInp.x','Indegree.x'), samplesize = 20000)
n.itm.e = n.itm.e%>% select(! any_of(dupcols))

a=sapply(n.itm.e, function(x) sum(is.na(x)));a[a>0]
#a=sapply(n.itm.v, function(x) sum(x==0));a[a>0]

n.itm.e=my_replace_na(n.itm.e,rplist =
                        list(L3.x=0,L3.y=0,L4.x=0,L4.y=0 ))

rm(n.itm.v.influx,n.itm.v.outflux,n.itm.v1,n.itm.e.tmp,n.itm.e.noselfedge)
#  curating missings  ----


#   Imputing missings n.itm.e

a=sapply(n.itm.e, function(x) sum(is.na(x)));a[a>0]
n.itm.e=my_replace_na(n.itm.e,rplist =
                        list(confZ=0,contribZ=0,confNZ=0,contribNZ=0,
                             tz=0,tnz=0,PinN=0,PoutN=0,PDiff=0,PAdd=0,
                             PdirDiff=0,PdirRatio=0, PoutBAN=0,PinBAN=0, Tendency=0 ,
                             Pinsum.x=0,PinAbsSum.x=0,PinMean.x=0, Poutsum.y=0, PoutAbsSum.y=0, PoutMean.y =0 ))
#   Add inxformation theoretic indices -----
#pr = P(y|X) , conf
#e=e %>%dplyr::mutate(LgPy_=log(probInp.y),Lg1=log(p_/probInp.y),)
#e=n.itm.e %>%dplyr::mutate(LgPy_=log(probInp.y),Lg1=log(p_/probInp.y))
#   Imputing infinites n.itm.e NO NA should be in data ----
setDT(n.itm.e)
a=sapply(n.itm.e, function(x) sum(is.na(x)));a[a>0]
c=sapply(n.itm.e, function(x) sum(is.infinite(x)));c[c>0]
b=n.itm.e[src=='0166D',.(src,trgt,Weight,OcrOut.x,OcrInp.y,Pin,Pinsum.x)];View(b)

#=do.call(data.frame,lapply(t, function(x) replace(x, is.infinite(x),NA)))
# replaceList=list(confZ=0,contribZ=0,confNZ=0,
#                 contribNZ=0, oddsN=myMax(n.itm.e$oddsN ),odds=myMax(n.itm.e$odds ),
#                 PinN=myMax(n.itm.e$PinN ),PoutN=myMax(n.itm.e$PoutN ),PDiff=myMax(n.itm.e$PDiff ),
#                 PAdd=myMax(n.itm.e$PAdd ),PdirDiff=myMax(n.itm.e$PdirDiff ),PdirRatio=myMax(n.itm.e$PdirRatio ) )
# n.itm.e=my_replace_val(n.itm.e,Inf,replaceList)

#myMax(n.itm.e$confZ)
#RENV3 saved



#   Finalize calculations on edges -----
# n.itm.e = n.itm.e %>%dplyr::mutate (vertexOcrSum=vertexOcrSumIn,
#                                     OcrOut.x2=OcrOut.x^2,trgtOcr2=OcrInp.y^2)



n.itm.e=n.itm.e %>% mutate(DiagnosesAb.x='',DiagnosesAb.y='')
#   Report Improtant Causal relation ------
setDF(n.itm.e)
n.itm.e$DiagnosesAb.x<-trim(n.itm.e$DiagnosesAb.x);n.itm.e$DiagnosesAb.y<-trim(n.itm.e$DiagnosesAb.y)
e= n.itm.e %>% dplyr::select(src,trgt,OcrOut.x,OcrInp.y,DiagnosesAb.x,DiagnosesAb.y,Weight,everything())
e = e %>%mutate(SUID=row_number(),class1="",class2="0",class3="0",
                class14 ="",shared_name=paste0(src,"-",trgt))
causal.v = n.itm.v


rm(n.itm.v.back,n.itm.e.back,n.itm.e.BeforeAugmenting,n.itm.e.BeforeAugmenting,n.itm.e.noselfedge,n3d.v,n.itm.v1,n.itm.v.BeforeMoments)


#   INACTIVE Use another subgraph  ------
# e<-gettable(paste(url.network,"208871","tables/defaultedge",sep="/"),"edgetable.json")
# causal.v<-gettable(paste(url.network,"208871","tables/defaultnode",sep="/"),"edgetable.json")
#     Add more fields 
#causal.v=causal.v[,c('vID','DiagnosesAb')];e$causal11<-as.numeric(e$causal11)
#e=left_join(e,causal.v,by=c("src"="vID"));e.causal2=left_join(e,causal.v31,by=c("trgt"="vID"))
#   Correct field names x denotes source and y denotes target ----
names(n.itm.e) <- gsub("_x", ".x", names(n.itm.e), fixed = TRUE)
names(n.itm.e) <- gsub("_y", ".y", names(n.itm.e), fixed = TRUE)
names(causal.v) <- gsub("_x", ".x", names(causal.v), fixed = TRUE)
names(causal.v) <- gsub("_y", ".y", names(causal.v), fixed = TRUE)
#n.itm.e=n.itm.e[,-grep("x$",colnames(n.itm.e))]; n.itm.e=n.itm.e[,-grep("y$",colnames(n.itm.e))]
#Environment saved as 2.rdata
#removing duplicated columns
library(stringr)
length(colnames(n.itm.e));length(unique(colnames(n.itm.e)))
n.itm.e <- n.itm.e[ , !duplicated(colnames(n.itm.e))]
PSUDEO_ZERO_2 = 0.000001

# n.itm.e= n.itm.e %>%
#   dplyr::mutate(
#     confZ=(conf-scfMedian.x)/ifelse(scfMAD.x==0,PSUDEO_ZERO_2,scfMAD.x), #Zscore based on median and MAD
#     contribZ=(contrib-ocbMedian.y)/ifelse( ocbMAD.y ==0,PSUDEO_ZERO_2 ,ocbMAD.y ) ,
#     tz=confZ+contribZ)
#TODO   Add vertice level 2 distribution parameters ----  
#TODO probably wrong should use edges not vertices to find those distribution params


#TODO   Descriptive analysis on n.itm.v to identify important variables ----

#TODO   Add source and target charectristics to edges -----

#TODO   Add inxformation theoretic indices -----


#TODO   Finalize calculations on edges -----
#TODO   Report Improtant Causal relation ------


#   TRAINING H2O  For Prediction and feature selection ----
#Remove large objects
sort( sapply(ls(),function(x){object.size(get(x))}),decreasing = F)

require(h2o)
h2o.shutdown(prompt=FALSE)
gc()
H2OCnn = h2o.init(nthreads = parallel::detectCores()-4, enable_assertions = TRUE,max_mem_size = "10g")
#memory.limit(size=25000)
#Create binary outcome to calculate probaibility
#Pretraining with tst2
{
  tst2.train.h2o = as.h2o(tst1.totalset)
  tst2.rf = h2o.randomForest(mdlColNames,"class2",tst2.train.h2o,#ntrees = 1)
                             ntrees=8,nfolds=3) #,max_depth=6)
  tst2.mdl = tst2.rf
  
  
  
  tst2.gbm = h2o.gbm(mdlColNames,"class1",tst2.train.h2o, max_depth=3,nfolds = 2 ) #ntrees=20
  tst2.mdl = tst2.gbm
  
  plot(tst2.mdl)
  h2o.performance(tst2.mdl) #save(tst1,tst1.train,tst1.tst, file="temp training set.robj")
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
#mdlColNames = intersect(mdlColNames, prd.varimp$variable[1:150])
mdlColNames=setdiff(mdlColNames,filterOutFeatures)
d.new = tst1.tst # sample_frac(tst1.tst,size = .3)
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
                           #max_depth =10,ntrees=15, 
                           #max_depth =9,ntrees=12, 
                           keep_cross_validation_predictions= FALSE,  validation_frame = tst1.tst.h2o)
#,ntrees=30,max_depth=5#selectedFeatures #mdlColNames
tst1.mdl=tst1.rf

tst1.gbm = h2o.gbm(mdlColNames,trainingTarget,
                   tst1.train.h2o,nfolds = NFolds,#checkpoint=tst2.mdl@model_id,
                   max_depth =3, 
                   #ntrees=12, 
                   keep_cross_validation_predictions= FALSE,  validation_frame = tst1.tst.h2o ) #ntrees=20
tst1.mdl = tst1.gbm

mdl.class3 = tst1.mdl



#mdlname = h2o.saveModel(mdl.class3, "RF model early prunning mHSC-E",force = T);mdlname 
#mHSC-E\\DRF_model_R_1629820087068_1
#mdl.class3 #DRF_model_R_1629166232553_186"
# mdl.class2 = "DRF_model_R_1629166232553_347"
#mdl.class3= mHSC-E\\DRF_model_R_1629823793412_1"

#mdl.class3=h2o.loadModel( "C:/E/HDDs/Work/NYU Research/0 Plant biology/R code/RF model early prunning HESC/DRF_model_R_1629166232553_186")
# tst1.gbm = h2o.gbm(mdlColNames,trainingTarget,tst1.train.h2o,nfolds = NFolds,ntree=30,max_depth=10,  validation_frame = tst1.tst.h2o )
# #,ntrees=30,max_depth=5
# #selectedFeatures #mdlColNames
# #,validation_frame = tst1.tst.h2o)
# #,checkpoint="DRF_model_R_1467413037220_19")
# tst1.mdl=tst1.gbm

#save(tst1.mdl,t1,t2,tst1.totalset,tst1.train,tst1.tst,file="mdl1.robj")

#Testing performance
{
  plot(tst1.mdl);setDF(tst1.totalset)
  h2o.performance(tst1.mdl) #save(tst1,tst1.train,tst1.tst, file="temp training set.robj")
  prd.varimp=h2o.varimp(tst1.mdl);View(prd.varimp)
  
  tstexplain =  t1 %>% anti_join(t2,by=c("src"="src", "trgt"="trgt")) %>% sample_frac(size =.02)
  h2o.explain(tst1.mdl,as.h2o(tstexplain))
  h2o.ice_plot(tst1.mdl,as.h2o(tstexplain),'scftau4.x')
  
  paste0(prd.varimp$variable[1:30],collapse= "','")
  #prd.varimp[1:15,] %>% knitr::kable()
  
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
  
  d.new = t1 %>% anti_join(t2,by=c("src"="src", "trgt"="trgt")) %>% sample_frac(size =.2)   #tst1.totalset #newDataEnv$tst1.totalset # 
  predTest.d.h2o = as.h2o(d.new) # 
  h2o.pred = as.data.frame(h2o.predict(tst1.mdl,predTest.d.h2o,keep_cross_validation_predictions=TRUE))
  h20.prediction=as.numeric(as.character(h2o.pred[,3]))
  predictions =h20.prediction #pred #ens.predictions #
  outcomes =unlist(setDT(d.new)[,trainingTarget,with=FALSE]) #outcome # ens.outcome# 
  set.seed(runif(1,1,1000))
  rndPred = ifelse(runif(nrow(d.new)) >=.5,1,0) #Assign a random classifier results
  
  d.new1 = as.data.frame( cbind(d.new[,.(src,trgt,Weight)],predictions,outcomes,rndPred))
  #d.new1 = d.new1[!duplicated(d.new1),]
  d.new1.rv =d.new1 %>% dplyr::rename(revWeight = Weight,rvpred=predictions,src1=src,trgt1=trgt) %>% 
    dplyr::select(-outcomes,-rndPred)
  pred_outcome = d.new1 %>% left_join(d.new1.rv, by=c("src"="trgt1", "trgt"="src1") )
  pred_outcome.back = pred_outcome
  
  #Late pruning using ground truth
  if(FALSE) pred_outcome = pred_outcome %>% filter(src %in% tf.all)
  
  
  
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
  theROC <- roc(assespreds$outcomes, assespreds$predictions, percent = TRUE);theROC
  
  sscurves <- evalmod(scores = assespreds$predictions, labels = assespreds$outcomes)
  pr.crv=sscurves$prcs[1][[1]];  pr.auc =  attr(pr.crv,'auc')
  
  theauc = precrec::auc(sscurves); paste0(theauc[[3]],"=", round(theauc[[4]],3))
  autoplot(sscurves, mode = 'rocpr',ret_grob = F)
  
  #partial precision-recall
  sscurves.part <- part(sscurves, xlim = c(0, 0.2))
  plot(sscurves.part,title = "Partial top 20% predictions")
  reportAUC(sscurves.part)  
  pr.prtcrv=sscurves.part$prcs[1][[1]];  pr.prtauc =  attr(pr.prtcrv,'pauc')
  
  
  #random classifier
  print("Random Classifier comparison ============")
  randomClassifierCurves <- evalmod(scores = assespreds$rndPred, labels = assespreds$outcomes)
  rnd.crv=randomClassifierCurves$prcs[1][[1]]; rnd.auc =  attr(rnd.crv,'auc')
  
  rndmClscurves.part <- part(randomClassifierCurves, xlim = c(0, 0.2))
  reportAUC(rndmClscurves.part)
  
  rnd.prtcrv=rndmClscurves.part$prcs[1][[1]]; rnd.prtauc =  attr(rnd.prtcrv,'pauc')
  
  print(sprintf("AUCPR Ratio CICT to Random= %s,  Partial AUCPR Ratio CICT to Random= %s ", 
                round(pr.auc/rnd.auc,2), round(pr.prtauc/rnd.prtauc,2)))

  
  mmpoins <- evalmod(scores = assespreds$predictions, labels = asses`preds$outcomes, mode = "basic")
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
  bestcutoff =as.double(coords(theROC, "best", best.method = "closest.topleft", best.weights=best.weights,
                               ret="threshold", transpose = FALSE));bestcutoff
  print(sprintf("Best threshold: %s for fn/fp cost ratio: %s " ,bestcutoff,relativeCostfn_fp ))
  #bestcutoff = 0.5
  assespreds =assespreds %>% mutate(thresholdpreds= ifelse(assespreds$predictions>bestcutoff,TRUE,FALSE)) 
  
  #post processing pruning%>%
  if(FALSE){
    assespreds =assespreds %>% mutate(thresholdpreds = 
                                        ifelse(src %in% tf.all)
                                      
                                      assespreds =assespreds %>% mutate(thresholdpreds = 
                                                                          ifelse(!(src %in% tbl.goldStandard$src | trgt %in% tbl.goldStandard$trgt), FALSE,thresholdpreds ))
                                      
                                      assespreds =assespreds %>% mutate(thresholdpreds = 
                                                                          ifelse(!((src %in% geneorder.1e3$X &src %in% tbl.goldStandard$src) | 
                                                                                     trgt %in% geneorder.1e3$X & trgt %in% tbl.goldStandard$trgt), FALSE,thresholdpreds ))
  }
  
  pander::pander(ftable(factor(assespreds$thresholdpreds, c("TRUE","FALSE")), 
                        factor(assespreds$outcomes, c("TRUE","FALSE")), 
                        dnn=c('pred','actual'))) #%>%  knitr::kable()
  #table(ifelse(predictions>.99,'c','other'))
  #bestcutoff = coords(roc, "best", best.method = "closest.topleft", ret="threshold", transpose = FALSE);bestcutoff
  
}

