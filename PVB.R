if (!("TDA" %in% rownames(installed.packages()))) install.packages("TDA") # install TDA package if it is not installed
library(TDA) # load TDA package

PVB = function(PD,dimension,xSeq,ySeq,alpha){
      # PD - N by 2 matrix (1st and 2nd columns contain birth and persistence values respectively)
      # dimension - homological dimension (0 for H0, 1 for H1, etc.)
      # xSeq - sequence of x (birth) values of the grid vertices
      # ySeq - sequence of y (persistence) values of the grid vertices
      # alpha - parameter (between 0 and 1) controlling block width 
      ##############
      vectorize1D <- function(c,d,PD){
        y <- PD[,2] # vector of persistence values
        lambda <- alpha*y # vector of (block widths)/2
        yMin <- pmax(c,y-lambda)
        yMax <- pmin(d,y+lambda)
        B <- pmax(0,yMax-yMin)
        sum(B*(yMin+yMax)/2) # weight function w(x,y)=y
      }
      ##############
      vectorize2D <- function(a,b,c,d,PD){
        x <- PD[,1] # vector of birth values 
        y <- PD[,2] # vector of persistence values
        lambda <- alpha*y # vector of (block widths)/2
        xMin <- pmax(a,x-lambda)
        xMax <- pmin(b,x+lambda)
        yMin <- pmax(c,y-lambda)
        yMax <- pmin(d,y+lambda)
        A <- pmax(0,xMax-xMin)
        B <- pmax(0,yMax-yMin)
        sum(A*B*(yMin+yMax+xMin+xMax)/2) # weight function w(x,y)=x+y
      }
  # Body of PVB() function
  dy <- diff(ySeq)
  m <- length(dy)
  if (dimension==0){
    pvb <- numeric(length = m)
    for (i in 1:m){
      c <- ySeq[i]
      d <- ySeq[i+1]
      pvb[i] <- vectorize1D(c,d,PD)/dy[i]
    }
  }
  else{
    dx <- diff(xSeq)
    n <- length(dx)
    pvb <- numeric(length = n*m)
    for (i in 1:m){
      c <- ySeq[i]
      d <- ySeq[i+1]
      for (j in 1:n){
        a <- xSeq[j]
        b <- xSeq[j+1]
        pvb[j+(i-1)*n] <- vectorize2D(a,b,c,d,PD)/(dx[j]*dy[i])
      }
    }
  }
  return(pvb)
}

# Example
N <- 100
X <- circleUnif(N,r=1)+rnorm(2*N,mean = 0,sd = 0.1) # sample N points from unit circle and add random noise 
plot(X,asp = 1,pch=19,xlab = 'x',ylab = 'y') 

PD <- ripsDiag(X,maxdimension = 1,maxscale = 2)$diagram # compute PD using Rips filtration
PD[,3] <- PD[,3] - PD[,2] # switch from birth-death to birth-persistence coordinates
colnames(PD)[3] <- "Persistence"

xSeq <- seq(0,0.5,by=0.1) # sequence of x (birth) values of the grid vertices
ySeq <- seq(0,2,by=0.25) # sequence of y (persistence) values of the grid vertices

PVB(PD[PD[,1]==0,2:3],dimension = 0,NULL,ySeq,alpha = 0.5) # compute PVB for homological dimension H0 and alpha=0.5
PVB(PD[PD[,1]==1,2:3],dimension = 1,xSeq,ySeq,alpha = 0.5) # compute PVB for homological dimension H1 and alpha=0.5
