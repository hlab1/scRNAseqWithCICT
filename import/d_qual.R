Diff_q <- function(x_i,x_j,y_i,y_j)
{
	#print(paste(x_i,"       ",x_j))
	##print(paste(y_i,"       ",y_j))
	x_i_j <- (x_i-x_j)
	y_i_j <- (y_i-y_j)
	#print(paste(y_i,"-",y_j,"=",y_i_j))
	if(abs(x_i_j) < 0.00001) x_i_j <- 0
	if(abs(y_i_j) < 0.00001) y_i_j <- 0
	
	if(y_i_j==0) {if(x_i_j==0) return(0) else return(0.5)}
	else if((x_i_j/y_i_j)>0) return(0)
	else if((x_i_j/y_i_j)==0) return(0.5)
	else if((x_i_j/y_i_j)<0) return(1)
}
Dq_PAIRWISE <- function(X,Y)
{
	if(length(X)!=length(Y)) return(-1)
	#print(X)
	##print(Y)
	##print(X[2])
	ts_length = length(X)
	d_qual = 0
	for(i in 2:ts_length-1)
	{
		##cat(i,sep="\n")
		j <- (i+1)
		while(j<=ts_length)
		{
			##cat(paste(X[i],X[j],Y[i],Y[j], sep="\t"), sep="\n")
			d_qual = d_qual + (2 - (2*Diff_q(X[i],X[j],Y[i],Y[j])/(ts_length*(ts_length-1))))
			j<-j+1
		}
	}

	return(d_qual)
}
Dq <- function(EG_MATRIX,N_GENES)
{	
	Dqual_MATRIX <- array(-1,c(N_GENES,N_GENES))
	for(i in 1:N_GENES) 
	{
		j <- (i+1)
		while(j<=N_GENES)
		{
			Dqual_MATRIX[i,j] <- Dq_PAIRWISE(EG_MATRIX[,i],EG_MATRIX[,j])
			Dqual_MATRIX[j,i] <- Dq_PAIRWISE(EG_MATRIX[,j],EG_MATRIX[,i])	
			j <- j+1
		}
		Dqual_MATRIX[i,i] <- 0.0
		#cat(paste("Finished with ",i," gene", sep=""), sep="\n")
	}
	return(Dqual_MATRIX)
}
Dq_PAIRWISE_dyn <- function(X,Y)
{
	if(length(X)!=length(Y)) return(-1)
	##print(X)
	##print(Y)
	##print(X[2])
	ts_length = length(X)
	d_qual = 0
	for(i in 2:ts_length-1)
	{
		##cat(i,sep="\n")
		j <- (i+1)
		while(j<=ts_length)
		{
			##cat(paste(X[i],X[j],Y[i],Y[j], sep="\t"), sep="\n")
			d_qual = d_qual + (2 - (2*Diff_q(X[i],X[j],Y[i],Y[j])/((ts_length-i+1))))
			j<-j+1
		}
	}

	return(d_qual)
}
Dq_dyn <- function(EG_MATRIX,N_GENES)
{	
	Dqual_MATRIX <- array(-1,c(N_GENES,N_GENES))
	for(i in 1:N_GENES) 
	{
		j <- (i+1)
		while(j<=N_GENES)
		{
			Dqual_MATRIX[i,j] <- Dq_PAIRWISE_dyn(EG_MATRIX[,i],EG_MATRIX[,j])
			Dqual_MATRIX[j,i] <- Dq_PAIRWISE_dyn(EG_MATRIX[,j],EG_MATRIX[,i])	
			j <- j+1
		}
		Dqual_MATRIX[i,i] <- 1.0
		#cat(paste("Finished with ",i," gene", sep=""), sep="\n")
	}
	return(Dqual_MATRIX)
}
Dq_F_MATRIX <- function(EXRATE,N_EXP,N_GENES)
{
	TMP_G_E_MATRIX <- array(0,c(N_EXP,N_GENES))
	
	for(i in 1:N_EXP)
	{
		for(j in 1:N_GENES)
		{
			TMP_G_E_MATRIX[i,j] = EXRATE[j,i]
		}
	}
	return(TMP_G_E_MATRIX)

}
