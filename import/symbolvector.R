symbolvector <- function(EXRATE,N_EXP,nn)
{
  
 a1 <- (length(EXRATE[,1])+1)
 library(gtools);
 
 m <- ceiling(nn+1)
 l_pattern <- factorial(m)
 le <- ceiling(factorial(N_EXP)/(factorial(m)*factorial(N_EXP-m)))
 tt <- array(NA,c(le,m))
 i <- 0
 z <- seq(1,m)
 
 while(i < le)
 {
  i <- i+1
  for(j in 1:m)
  {
   tt[i,j] <- z[j]
  }
  
  z[m] <- z[m] + 1
  mm <- m
  mmm <- 0
  while(mm > 1)
  {
   if(z[mm] > (N_EXP-m+mm))
   {
    mmm <- mmm + 1
    z[mm-1] <- z[mm-1] + 1
    z[mm] <- 1
   }
   mm <- mm - 1
  }
  while(mmm > 0)
  {
   z[m-mmm+1] <- z[m-mmm] + 1
   mmm <- mmm - 1
  }
 }
 
 ind1 <- array(0,c(le,3))
 A <- array(0,c(le,(a1-1)))
 A2 <- array(0,c(le,(a1-1)))
 muster <- array(NA,c(l_pattern))
 paste_n <- function(ab,n) 
 {
  paste(ab[1:n],collapse="")
 }
 muster <- apply(permutations(n=m,r=m),1, paste_n,m)
 muster2 <- array(NA,c(l_pattern))
 for(mi in 1:l_pattern)
 {
  txt1 <- strsplit(muster[mi],"")
  muster2[mi] <- paste(txt1[[1]][m:1],collapse="")
 }
 for(i in 1:(a1-1))
 {
  for(t_count in 1:le)
  {
   A[t_count,i] <- which(paste(order(EXRATE[i,tt[t_count,]]),collapse="") == muster)
   A2[t_count,i] <- which(paste(order(EXRATE[i,tt[t_count,]]),collapse="") == muster2)
  }
 }
 return(list(A = A, A2 = A2, l_pattern = l_pattern))
}
