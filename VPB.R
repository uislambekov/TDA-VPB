if (!("TDAstats" %in% rownames(installed.packages()))) install.packages("TDAstats") # install TDAstats package if it is not installed
library(TDAstats) # load TDAstats package

VPB = function(D,dimension,xSeq,ySeq,tau){
  # D - N by 3 matrix (columns contain dimension, birth and persistence values respectively)
  # dimension - homological dimension (0 for H0, 1 for H1, etc.)
  # xSeq - sequence of x (birth) values of the grid vertices
  # ySeq - sequence of y (persistence) values of the grid vertices
  # tau - parameter (between 0 and 1) controlling block width 
  # weight function w(x,y) = x+y
  
  # Body of VPB() function
  D = D[D[, 1] == dimension, 2:3, drop = F]
  x <- D[,1] 
  y <- D[,2] 
  lambda <- tau*y 
  dy <- diff(ySeq)
  m <- length(dy)
  if (sum(abs(x))==0){
    vpb <- numeric(length = m)
    for (i in 1:m){
      c <- ySeq[i]
      d <- ySeq[i+1]
      yInd <- which((y>c-lambda)&(y<d+lambda))
      if (length(yInd)>0){
        y_cd <- y[yInd]
        lambda_cd <- lambda[yInd]
        yMin <- pmax(c,y_cd-lambda_cd)
        yMax <- pmin(d,y_cd+lambda_cd)
        vpb[i] <- 0.5*sum(yMax^2-yMin^2)/dy[i] 
      } 
      else vpb[i] <- 0
    }
  }
  else{
    dx <- diff(xSeq)
    n <- length(dx)
    vpb <- numeric(length = n*m)
    for (i in 1:m){
      c <- ySeq[i]
      d <- ySeq[i+1]
      yInd <- which((y>c-lambda)&(y<d+lambda))
      if (length(yInd)>0){
        y_cd <- y[yInd]
        x_cd <- x[yInd]
        lambda_cd <- lambda[yInd]
        yMin_cd <- pmax(c,y_cd-lambda_cd)
        yMax_cd <- pmin(d,y_cd+lambda_cd)
        for (j in 1:n){
          a <- xSeq[j]
          b <- xSeq[j+1]
          xInd <- which((x_cd>a-lambda_cd)&(x_cd<b+lambda_cd))
          if (length(xInd)>0){
            x_abcd <- x_cd[xInd]
            lambda_abcd <- lambda_cd[xInd]
            xMin <- pmax(a,x_abcd-lambda_abcd)
            xMax <- pmin(b,x_abcd+lambda_abcd)
            yMin <- yMin_cd[xInd]
            yMax <- yMax_cd[xInd]
            vpb[j+(i-1)*n] <- 0.5*sum((xMax-xMin)*(yMax-yMin)*(xMin+xMax+yMin+yMax))/(dx[j]*dy[i]) 
          }
          else vpb[j+(i-1)*n] <- 0
        }
      }
      else vpb[1:n+(i-1)*n] <- 0
    }
  }
  return(vpb)
}

# Example
# `TDAstats::circle2d` dataset contains the coordinates of 
# 100 points uniformly distributed on the circumference of
# a unit circle
X <- circle2d + rnorm(200,mean = 0,sd = 0.1) # add Gaussian noise  
plot(X,asp = 1,pch=19,xlab = 'x',ylab = 'y') 

# compute PD using Rips filtration
D <- calculate_homology(X,dim = 1,threshold = 2) 
# switch from birth-death to birth-persistence coordinates
D[,3] <- D[,3] - D[,2] 
colnames(D)[3] <- "persistence"

# compute VPB for homological dimension H0 and tau=0.5
ySeq <- unique(quantile(D[D[,1]==0,3],probs = seq(0,1,by=0.1))) 
VPB(D,dimension = 0,xSeq=NULL,ySeq,tau = 0.5) 

# compute VPB for homological dimension H1 and tau=0.5
xSeq <- unique(quantile(D[D[,1]==1,2],probs = seq(0,1,by=0.1)))
ySeq <- unique(quantile(D[D[,1]==1,3],probs = seq(0,1,by=0.1)))
VPB(D,dimension = 1,xSeq,ySeq,tau = 0.5) 
