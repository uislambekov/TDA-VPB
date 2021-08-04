if (!("TDAstats" %in% rownames(installed.packages()))) install.packages("TDAstats") # install TDAstats package if it is not installed
library(TDAstats) # load TDAstats package

VPB = function(D,dimension,xSeq=NULL,ySeq,tau){
      # D - N by 3 matrix (columns contain dimension, birth and persistence values respectively)
      # dimension - homological dimension (0 for H0, 1 for H1, etc.)
      # xSeq - sequence of x (birth) values of the grid vertices
      # ySeq - sequence of y (persistence) values of the grid vertices
      # tau - parameter (between 0 and 1) controlling block width 
      ##############
      vectorize1D <- function(c,d,D){
        y <- D[,2] # vector of persistence values
        lambda <- tau*y # vector of (block widths)/2
        yMin <- pmax(c,y-lambda)
        yMax <- pmin(d,y+lambda)
        B <- pmax(0,yMax-yMin)
        sum(B*(yMin+yMax)/2) # weight function w(x,y)=y
      }
      ##############
      vectorize2D <- function(a,b,c,d,D){
        x <- D[,1] # vector of birth values 
        y <- D[,2] # vector of persistence values
        lambda <- tau*y # vector of (block widths)/2
        xMin <- pmax(a,x-lambda)
        xMax <- pmin(b,x+lambda)
        yMin <- pmax(c,y-lambda)
        yMax <- pmin(d,y+lambda)
        A <- pmax(0,xMax-xMin)
        B <- pmax(0,yMax-yMin)
        sum(A*B*(yMin+yMax+xMin+xMax)/2) # weight function w(x,y)=x+y
      }
  # Body of VPB() function
  D = D[D[, 1] == dimension, 2:3, drop = F]
  dy <- diff(ySeq)
  m <- length(dy)
  if (dimension==0){
    vpb <- numeric(length = m)
    for (i in 1:m){
      c <- ySeq[i]
      d <- ySeq[i+1]
      vpb[i] <- vectorize1D(c,d,D)/dy[i]
    }
  }
  else{
    dx <- diff(xSeq)
    n <- length(dx)
    vpb <- numeric(length = n*m)
    for (i in 1:m){
      c <- ySeq[i]
      d <- ySeq[i+1]
      for (j in 1:n){
        a <- xSeq[j]
        b <- xSeq[j+1]
        vpb[j+(i-1)*n] <- vectorize2D(a,b,c,d,D)/(dx[j]*dy[i])
      }
    }
  }
  return(vpb)
}

# Example
X <- circle2d+rnorm(200,mean = 0,sd = 0.1) # sample 100 points from unit circle and add random noise 
plot(X,asp = 1,pch=19,xlab = 'x',ylab = 'y') 

D <- calculate_homology(X,dim = 1,threshold = 2) # compute PD using Rips filtration
D[,3] <- D[,3] - D[,2] # switch from birth-death to birth-persistence coordinates
colnames(D)[3] <- "Persistence"

xSeq <- seq(0,0.5,by=0.1) # sequence of x (birth) values of the grid vertices
ySeq <- seq(0,2,by=0.25) # sequence of y (persistence) values of the grid vertices

VPB(D,dimension = 0,ySeq=ySeq,tau = 0.5) # compute VPB for homological dimension H0 and tau=0.5
VPB(D,dimension = 1,xSeq=xSeq,ySeq=ySeq,tau = 0.5) # compute VPB for homological dimension H1 and tau=0.5
