################################################################################@
# This is a modified edition of the code provided by Kuzmanovski et al. to
#  produce raw edge measurements:

# Kuzmanovski V, Todorovski L, Džeroski S. Extensive evaluation of the
#  generalized relevance network approach to inferring gene regulatory networks.
#  Gigascience. 2018 Nov 1;7(11):giy118. doi: 10.1093/gigascience/giy118.
#  PMID: 30239704; PMCID: PMC6420648.
################################################################################@
my.build.mim <- function (dataset, estimator = "spearman", disc = "none", 
                          nbins = sqrt(NROW(dataset))) 
{
  
  if (disc == "equalfreq" || disc == "equalwidth" || 
      disc == "globalequalwidth") 
    dataset <- infotheo::discretize(dataset, disc, nbins)
  if (estimator == "pearson" || estimator == "spearman" || 
      estimator == "kendall") {
    
    mim <- cor(dataset, method = estimator, use = "pairwise.complete.obs")^2  # My change  complete.obs
    diag(mim) <- 0
    maxi <- 0.999999
    mim[which(mim > maxi)] <- maxi
    mim <- -0.5 * log(1 - mim)
  }
  else if (estimator == "mi.mm") 
    estimator = "mm"
  else if (estimator == "mi.empirical") 
    estimator = "emp"
  else if (estimator == "mi.sg") 
    estimator = "sg"
  else if (estimator == "mi.shrink") 
    estimator = "shrink"
  else stop("unknown estimator")
  if (estimator == "mm" || estimator == "emp" || 
      estimator == "sg" || estimator == "shrink") {
    mim <- infotheo::mutinformation(dataset, method = estimator)
    diag(mim) <- 0
  }
  mim[mim < 0] <- 0
  mim
}

################ Get and setup data ###############################
getExpresionMatrix <- function(filePath, nrows, ncolumns, byrow = TRUE,what='character')
{
	print("Reading expresion data")
	if(file.exists(filePath))
	{
		print("File confirmed!");
	}
	E <- scan(file = filePath, what = what, sep = ",")
	EE <- matrix(E,nrows,ncolumns,byrow = byrow)
	return(EE)
}

getDirectedGoldStandard <- function(ExpresionMatrix, nrows, filePath,byrow = TRUE)
{
	print("Reading directed gold standard")

	LNetz <- length(scan(file = filePath, what = 'character'))/3
	Netz <- matrix(scan(file = filePath, what = 'character'),LNetz,3,byrow = byrow) 
	Netznr <- array(0,c(LNetz,3))
	for (i in 1:LNetz)
	{
 		Netznr[i,1] <- which(Netz[i,1] == ExpresionMatrix[,1])-1
 		Netznr[i,3] <- which(Netz[i,3] == ExpresionMatrix[,1])-1
 		if (Netz[i,2] == "ac")
 		{
  			Netznr[i,2] <- 1
 		}
 		if (Netz[i,2] == "re")
 		{
  			Netznr[i,2] <- -1
 		}
 		if (Netz[i,2] == "du")
 		{
  			Netznr[i,2] <- 0
 		}
	}
	NetzMat <- array(0,c((nrows-1),(nrows-1)))
	for(k in 1:LNetz)
	{
 		NetzMat[Netznr[k,1],Netznr[k,3]] <- 1
	}
	grorg <- graph.adjacency(NetzMat,mode=c("DIRECTED"))

	return(NetzMat)
}

getUndirectedGoldStandard <- function(ExpresionMatrix, nrows,filePath,byrow = TRUE)
{
	print("Reading undirected gold standard")

	LNetz <- length(scan(file = filePath, what = 'character'))/3
	Netz <- matrix(scan(file = filePath, what = 'character'),LNetz,3,byrow = byrow) 
	Netznr <- array(0,c(LNetz,3))
	for (i in 1:LNetz)
	{
 		Netznr[i,1] <- which(Netz[i,1] == ExpresionMatrix[,1])-1
 		Netznr[i,3] <- which(Netz[i,3] == ExpresionMatrix[,1])-1
 		if (Netz[i,2] == "ac")
 		{
  			Netznr[i,2] <- 1
 		}
 		if (Netz[i,2] == "re")
 		{
  			Netznr[i,2] <- -1
 		}
 		if (Netz[i,2] == "du")
 		{
  			Netznr[i,2] <- 0
 		}
	}
	NetzMatsym <- array(0,c((nrows-1),(nrows-1)))
	for(k in 1:LNetz)
	{
 		NetzMatsym[Netznr[k,1],Netznr[k,3]] <- 1
 		NetzMatsym[Netznr[k,3],Netznr[k,1]] <- 1
	}
	grorgundir <- graph.adjacency(NetzMatsym,mode=c("UNDIRECTED")) 

	return(NetzMatsym)
}

getCalculatedReplicates<-function(ExpressionMatrix, ngenes, nexpr, nrepls)
{
	print("Calculate mean of technical replicates")
	
	EE_m <- array(0,c((ngenes-1),(nexpr+1)))
	EE_m[,1] <- ExpressionMatrix[2:ngenes,1] 
	b1 <- 2
	j <- 1
	for(k in 1:nexpr)
	{
 		for(i in 2:ngenes)
 		{
  			EE_m[i-1,j+1] <- mean(as.numeric(ExpressionMatrix[i,b1:(b1+nrepls-1)]))
 		}
 		j <- j+1
 		b1 <- b1 + nrepls
	}
	EXRATE <- array(0,c((ngenes-1),nexpr))
	EXRATE[,1:nexpr] <- as.numeric(EE_m[,2:(nexpr+1)])

	return(EXRATE)

}







################ Relevance Networks Inferring ###############################
mutualinformation <- function(X,Y,methode,discretizers="equalwidth")
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

getNormalizedMatrix <- function(ExpressionMatrix,normalization="max")
{
	##.. normalization: max, minimax
  
  if(normalization=="max")
	{
		NormalizedMatrix <- ExpressionMatrix/max(ExpressionMatrix)
	}
	else if(normalization=="minimax")
	{
	  minMat <- min(ExpressionMatrix)
	  maxMat <- max(ExpressionMatrix)
		NormalizedMatrix <- (ExpressionMatrix-minMat)/(maxMat-minMat)
	}

	return(NormalizedMatrix)
}

getRowColNamedMatrix <- function(ExpressionMatrix, prefixs="V")
{
	rownames(ExpressionMatrix) <- rownames(ExpressionMatrix, do.NULL = FALSE, prefix = prefixs)
	colnames(ExpressionMatrix) <- colnames(ExpressionMatrix, do.NULL = FALSE, prefix = prefixs)

	return(ExpressionMatrix)
}

getSimilarityMatrix_MI <- function(ExpressionMatrix, nrows, estimators="pearson", subestimators="mm", discretization = FALSE, discretizator = "equalwidth", diagr=0)
{
  library("minet")
  require(WGCNA)
	##.. estimators[correlation]: pearson, spearman, kendall
	##.. estimators[mutual information]: mi.empirical, mi.mm, mi.shrink, mi.sg 
	##.. estimators[other]: coarse.grained, granger
	##.. subestimators[coarse.grained]: ML, MM, Jeffreys, Laplace, SG, minimax, CS, NSB, shrink
	##.. discretizator: equalwidth, equalfreq
	##.. diagr = replacement value of the main diagonal (default: diagr=0)
	
	if(estimators == "granger")
	{
		source("gc1.R")
		gc_simp <- array(0,c((nrows-1),(nrows-1)))
		for(i in 1:(nrows-1))
		{
 			if(i < (nrows-1))
 			{
  				for(j in (i+1):(nrows-1))
  				{
					T1 <- ExpressionMatrix[i,]
     					T2 <- ExpressionMatrix[j,]

					l1 <- VARselect(T1,lag.max = 8)$selection[[1]]
     					l2 <- VARselect(T2,lag.max = 8)$selection[[1]]
     					if(is.finite(l1) == 0)  l1 <- NA
     					if(is.finite(l2) == 0)  l2 <- NA
     					LAG <- floor(mean(c(l1,l2),na.rm = TRUE))
     					if(is.na(LAG)) LAG <- 1
      				gc_simp[i,j] <- granger(cbind(T2,T1), L=LAG)
      				gc_simp[j,i] <- granger(cbind(T1,T2), L=LAG)
   				}
 			}
		}

		mim <- gc_simp
		diag(mim) <- diagr

	}
	else if(estimators == "coarse.grained")
	{
		DATA <- t(ExpressionMatrix)
		L <- length(DATA[,1])
		Ixy <- array(0,c((nrows-1),(nrows-1)))
		tau_max <- (L-1)
		for(i in 1:(nrows-1))
		{
			
 			for(j in 1:(nrows-1))
 			{ 
  				for(tau in -tau_max: tau_max)
  				{

   					if(tau < 0)
   					{
    						X <- DATA[(-tau+1):L,i]
    						Y <- DATA[1:(L+tau),j]
    						I <- mutualinformation(X,Y,subestimators,discretizers = discretizator) #build.mim(data.frame(X,Y), estimator=subestimators)
   					}
   					if(tau > 0)
   					{
    						X <- DATA[1:(L-tau),i]
    						Y <- DATA[(tau+1):L,j]
    						I <- mutualinformation(X,Y,subestimators,discretizers = discretizator) #build.mim(data.frame(X,Y), estimator=subestimators)
   					}
   					if(tau == 0)
   					{
    						I <- 0
   					}
   					Ixy[i,j] <- Ixy[i,j] + I
  				}
  				Ixy[i,j] <- Ixy[i,j]/(2*tau_max)
 			}
		}
		mim <- Ixy
		diag(mim) <- diagr

	}
  else if(estimators =='pearsonFALSE' ){
    #NEEDS debugging
    require(WGCNA)
    allowWGCNAThreads()
    system.time({ mim <- WGCNA::corFast(ExpressionMatrix)})
  }
	else
	{
	  if(discretization == FALSE)
		{
			mim <- build.mim(t(ExpressionMatrix), estimator=estimators)
		}
		else
		{
		  #mim22 <- build.mim(discretize(t(ExpressionMatrix), discretizator), estimator = estimators)
		  mim <- build.mim(t(ExpressionMatrix), disc = discretizator, estimator = estimators)
		}
		diag(mim) <- diagr
	}
	return(mim)
}

getSimilarityMatrix_DISTANCES <- function(ExpressionMatrix, nrows, norms=10, diagr=0)
{
  ##.. norms: 1-> Manhattan, 2-> Euclidian, m->Lm-Norm (default=10)
	# DistanceMatrix <- array(0,c((nrows-1),(nrows-1)))
	# 
	# for(di in 1:(nrows-1))
	# {
	#  for(dj in 1:(nrows-1))
	#  {
	# 	DistanceMatrix[di,dj] <- (sum(abs(ExpressionMatrix[di,] - ExpressionMatrix[dj,])^norms))^(1/norms)
	#  }
	# }
  
  distMatObject <- dist(ExpressionMatrix,method="minkowski",p=norms)
  
  DistanceMatrix <- as.matrix(distMatObject)
  diag(DistanceMatrix) <- diagr

	return(DistanceMatrix)
}

getSimilarityMatrix_DTW <- function(ExpressionMatrix, distmethod="Euclidean", steppatern=asymmetric)
{
library("dtw")
	##.. distmethod: Manhattan, Euclidian (default),...
	##.. steppattern: asymmetric(default), symmetric1, symmetric2,...
  DTWMatrix <- dtwDist(ExpressionMatrix,method="DTW", keep.internals=TRUE,step.pattern=steppatern)

	return(DTWMatrix)
}

getSimilarityMatrix_SYMBOLIC <- function(ExpressionMatrix, nrows, npoints, simmethod="sym", npatterns=4, patterns = NULL, diagr=0, discretization = TRUE, discretizator = "equalwidth", mitype="mm", numCores=1)
{
source("symbolvector.R")
library("minet")
library("parallel")
library("foreach")
library("doParallel")

	##.. simmethod: sym, sym.mu, avg.sym.mi
	##.. npatterns: 1,2,3... number that maximizes no of combination (=npoints/2)
	##.. discretizator: equalwidth, equalfreq
	##.. diagr = replacement value of the main diagonal (default: diagr=0)


  if(is.null(patterns)){
    SVEC <- symbolvector(ExpressionMatrix,npoints,npatterns)  
  } else {
    SVEC <- patterns;
  }
	
	A <- SVEC$A
	A2 <- SVEC$A2
	l_pattern <- SVEC$l_pattern

	if(simmethod=="sym")
	{
        ##setup parallel backend to use numCores cores
        cl<-makeCluster(numCores, outfile="")
        registerDoParallel(cl)
        print("Cluster registered...")
    
		P1 <- array(0,c((nrows-1),(nrows-1)))
		P2 <- array(0,c((nrows-1),(nrows-1)))
		P3 <- array(0,c((nrows-1),(nrows-1)))
        
		for(i in 1:(nrows-1))
		{
 			if((i+1) <= (nrows-1))
 			{
  				
                P4i = foreach(j = (i+1):(nrows-1), .combine='c') %do%
                {
                    
                    p1 <- 0
   					p2 <- 0
   					for(pl in 1:l_pattern)
   					{
    						p1 <- p1 + ((length(which(A[,i]==pl & A[,j]==pl)))/length(A[,i]))
    						p2 <- p2 + ((length(which(A[,i]==pl & A2[,j]==pl)))/length(A2[,i])) 
   					}
   					
                    P4tmp <- max(c(p1,p2))
                    as.numeric(P4tmp)
                }
                
                
                P4pre <- array(NA,c(1,i))
                P3[i,] <- c(P4pre,P4i)
 			}
			
		}
        
        stopCluster(cl)
        print("Cluster stopped!")
        
        for(i in 1:(nrows-1)){
            for(j in i:(nrows-1)) {
                if(i == j) {
                    P3[i,j] <- 0
                } else {
                    P3[j,i] <- P3[i,j]
                }
            }
        }
        
		SimMilarityMatrix <- P3
	}
	else if(simmethod=="sym.mi")
	{
		if(discretization==TRUE)
		{
			MI_A <- mutinformation(discretize(A,disc=discretizator),method=mitype)
		}
		else
		{
			MI_A <- mutinformation(A,method=mitype)
		}
		SimMilarityMatrix <- MI_A
	}
	else if(simmethod=="avg.sym.mi")
	{
        ##setup parallel backend to use numCores cores
        cl<-makeCluster(numCores, outfile="")
        registerDoParallel(cl)
        print("Cluster registered...")
    
		### Order pattern
		P1 <- array(0,c((nrows-1),(nrows-1)))
		P2 <- array(0,c((nrows-1),(nrows-1)))
		P3 <- array(0,c((nrows-1),(nrows-1)))
		for(i in 1:(nrows-1))
		{
 			if((i+1) <= (nrows-1))
 			{
  				#for(j in (i+1):(nrows-1))
  				#{
   				#	p1 <- 0
   				#	p2 <- 0
   				#	for(pl in 1:l_pattern)
   				#	{
    			#			p1 <- p1 + ((length(which(A[,i]==pl & A[,j]==pl)))/length(A[,i]))
    			#			p2 <- p2 + ((length(which(A[,i]==pl & A2[,j]==pl)))/length(A2[,i])) 
   				#	}
   				#	P1[i,j] <- p1
   				#	P1[j,i] <- p1
   				#	P2[i,j] <- p2
   				#	P2[j,i] <- p2
   				#	P3[i,j] <- max(c(P1[i,j],P2[i,j]))
   				#	P3[j,i] <- max(c(P1[j,i],P2[j,i]))
				#	
  				#}
                
                P4i = foreach(j = (i+1):(nrows-1), .combine='c') %do%
                {
                    
                    p1 <- 0
   					p2 <- 0
   					for(pl in 1:l_pattern)
   					{
    						p1 <- p1 + ((length(which(A[,i]==pl & A[,j]==pl)))/length(A[,i]))
    						p2 <- p2 + ((length(which(A[,i]==pl & A2[,j]==pl)))/length(A2[,i])) 
   					}
   					
                    P4tmp <- max(c(p1,p2))
                    as.numeric(P4tmp)
                }
                
                
                P4pre <- array(NA,c(1,i))
                P3[i,] <- c(P4pre,P4i)
 			}
		}

        stopCluster(cl)
        print("Cluster stopped!")
        
        for(i in 1:(nrows-1)){
            for(j in i:(nrows-1)) {
                if(i == j) {
                    P3[i,j] <- 0
                } else {
                    P3[j,i] <- P3[i,j]
                }
            }
        }
        
        
        
		### Order pattern + mi
		if(discretization==TRUE)
		{
			MI_A <- mutinformation(discretize(A,disc=discretizator),method=mitype)
		}
		else
		{
			MI_A <- mutinformation(A,method=mitype)
		}

		### Finall Avg. Order pattern+mi
		SimMilarityMatrix <- ((P3+MI_A)/2)
	}

	diag(SimMilarityMatrix) <- diagr

	return(SimMilarityMatrix)
}
getSimilarityMatrix_QUAL <- function(ExpressionMatrix, nrows, npoints)
{
source("d_qual.R")
	##.. nrows (a1) is number of genes + 1
	##.. npoints is number of time points within the time series	

	D_Q_MATRIX <- Dq_F_MATRIX(ExpressionMatrix,npoints,(nrows-1))
	D_Q <- Dq(D_Q_MATRIX,(nrows-1))

	return(D_Q)
}






getScorredMatrix <- function(SimilarityMatrix, scorrer="MRNET", aracne_eps=0)
{
library("minet")
	##.. scorrers: mrnet(default), clr, aracne, awe
	##.. aracne_eps is aracne parameter (see minet package manual)
	SimilarityMatrix <- getNormalizedMatrix(SimilarityMatrix,normalization="minimax")
	
	if(scorrer=="MRNET")
	{
	  ScorredMatrix <- tryCatch({
	    return(mrnet(SimilarityMatrix))
	  },error=function(err){
	    print("Error thrown in MRNET!")
	    print(err);
	    return(NULL);
	  },warning=function(warn){
	    print("Warning thrown in MRNET!")
	    print(warn);
	    return(NULL);
	  })
		
    #ScorredMatrix <- mrnet(SimilarityMatrix)
	}
	else if(scorrer=="CLR")
	{
	  ScorredMatrix <- tryCatch({
	    return(clr(SimilarityMatrix))
	  },error=function(err){
	    print(SimilarityMatrix)
	    print("Error thrown in CLR!")
	    print(err);
	    return(NULL);
	  },warning=function(warn){
	    print("Warning thrown in CLR!")
	    print(warn);
	    return(NULL);
	  })
	  #ScorredMatrix <- clr(SimilarityMatrix)
	}
	else if(scorrer=="ARACNE")
	{
	  ScorredMatrix <- tryCatch({
	    return(aracne(SimilarityMatrix,eps=aracne_eps))
	  },error=function(err){
	    print("Error thrown in ARACNE!")
	    print(err);
	    return(NULL);
	  },warning=function(warn){
	    print("Warning thrown in ARACNE!")
	    print(warn);
	    return(NULL);
	  })
		#ScorredMatrix <- aracne(SimilarityMatrix,eps=aracne_eps)
	}
	else if(scorrer=="AWE")
	{
		ScorredMatrix <- SimilarityMatrix
		sm_length <- ncol(SimilarityMatrix)
		for(i in 1:sm_length)
		{
 			ScorredMatrix[,i] <- SimilarityMatrix[,i]/sum(SimilarityMatrix[,i])
		}
	}
	else
	{
		ScorredMatrix <- mrnet(SimilarityMatrix)
	}
	
	return(ScorredMatrix)
}







