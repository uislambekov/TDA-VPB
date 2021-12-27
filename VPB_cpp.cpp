#include <Rcpp.h>
#include <iostream>
using namespace Rcpp;

IntegerVector which_C(LogicalVector vec) {
  std::vector<int> outC;
  for(int i=0; i<vec.size(); ++i) {
    if(vec[i]) {
      outC.push_back(i);
    }
  };
  IntegerVector out(outC.cbegin(), outC.cend());
  return out;
}

// [[Rcpp::export]]
NumericVector VPB_C(NumericMatrix D_, int dimension, NumericVector xSeq, NumericVector ySeq,
                     double tau) {
//  D <- D[D[, 1] == dimension, 2:3, drop = F]
//  x <- D[,1] 
//  y <- D[,2] 
  int n_rows = 0; // number of rows with the correct dimension
  for(int i=0; i<D_.nrow(); ++i) {
    if(D_(i,0) == dimension) {
      n_rows = n_rows+1; 
    }
  };
  NumericVector x(n_rows), y(n_rows);
  int n=0;
  for(int i=0; i<D_.nrow(); ++i) {
    if( D_(i,0) == dimension) {
      x[n] = D_(i,1);
      y[n] = D_(i,2);
      ++n;
    }
  }
  NumericVector lambda = tau*y;
  NumericVector dy = diff(ySeq);
  int m = dy.size();
  std::vector<double> vpb;
  if (sum(abs(x))==0) {
    std::vector<double> vpb_(m, 2);
    vpb = vpb_;
    for(int i=0; i<m; ++i) {
      double c = ySeq[i];
      double d = ySeq[i+1];
      vpb[i] = 0;
      for(int ie=0; ie<y.size(); ++ie) {
        if( (y[ie] > c - lambda[ie]) & (y[ie] < d+lambda[ie])) {
          double y_cd_ = y[ie];
          double lambda_cd_ = lambda[ie];
          double yMin_ = std::max(c, y_cd_ - lambda_cd_);
          double yMax_ = std::min(d, y_cd_ + lambda_cd_);
          vpb[i] += 0.5*(yMax_*yMax_-yMin_*yMin_)/dy[i] ;
        }
      }
    }
    return NumericVector(vpb.begin(), vpb.end());
  }
  else{
    NumericVector dx = diff(xSeq);
    int n = dx.size();
    std::vector<double> vpb_(n*m);
    vpb = vpb_;
    for(int i=1; i<=m; ++i) {
      double c = ySeq[i-1];      
      double d =  ySeq[i+1-1];
      IntegerVector yInd =which_C( (y>c-lambda) & (y<d+lambda));
      if (yInd.size()>0){
        NumericVector y_cd = y[yInd];
        NumericVector x_cd = x[yInd];
        NumericVector lambda_cd = lambda[yInd];
        NumericVector yMin_cd = pmax(c,y_cd-lambda_cd);
        NumericVector yMax_cd = pmin(d,y_cd+lambda_cd);
        for (int j=1; j<=n; ++j) {
          double a = xSeq[j-1];
          double b = xSeq[j+1-1];
          IntegerVector xInd = which_C((x_cd>a-lambda_cd) & (x_cd<b+lambda_cd));
          if (xInd.size()>0) {
            NumericVector x_abcd = x_cd[xInd];
            NumericVector lambda_abcd = lambda_cd[xInd];
            NumericVector xMin = pmax(a,x_abcd-lambda_abcd);
            NumericVector xMax = pmin(b,x_abcd+lambda_abcd);
            NumericVector yMin = yMin_cd[xInd];
            NumericVector yMax = yMax_cd[xInd];
            vpb[j+(i-1)*n-1] = 0.5*sum((xMax-xMin)*(yMax-yMin)*(xMin+xMax+yMin+yMax))/(dx[j-1]*dy[i-1]);
          }
          else vpb[j+(i-1)*n-1] = 0;
        }
      }
      //else vpb[1:n+(i-1)*n] <- 0
    }
    return NumericVector(vpb.begin(), vpb.end());
  }
}


NumericVector VPB_C_fast(NumericMatrix D_, int dimension, NumericVector xSeq, NumericVector ySeq,
                    double tau) {
  //  D <- D[D[, 1] == dimension, 2:3, drop = F]
  //  x <- D[,1] 
  //  y <- D[,2] 
  int n_rows = 0; // number of rows with the correct dimension
  for(int i=0; i<D_.nrow(); ++i) {
    if(D_(i,0) == dimension) {
      n_rows = n_rows+1; 
    }
  };
  NumericVector x(n_rows), y(n_rows);
  int n=0;
  for(int i=0; i<D_.nrow(); ++i) {
    if( D_(i,0) == dimension) {
      x[n] = D_(i,1);
      y[n] = D_(i,2);
      ++n;
    }
  }
  NumericVector lambda = tau*y;
  NumericVector dy = diff(ySeq);
  int m = dy.size();
  std::vector<double> vpb;
  if (sum(abs(x))==0) {
    std::vector<double> vpb_(m, 2);
    vpb = vpb_;
    for(int i=0; i<m; ++i) {
      double c = ySeq[i];
      double d = ySeq[i+1];
      vpb[i] = 0;
      for(int ie=0; ie<y.size(); ++ie) {
        if( (y[ie] > c - lambda[ie]) & (y[ie] < d+lambda[ie])) {
          double y_cd_ = y[ie];
          double lambda_cd_ = lambda[ie];
          double yMin_ = std::max(c, y_cd_ - lambda_cd_);
          double yMax_ = std::min(d, y_cd_ + lambda_cd_);
          vpb[i] += 0.5*(yMax_*yMax_-yMin_*yMin_)/dy[i] ;
        }
      }
    }
    return NumericVector(vpb.begin(), vpb.end());
  }
  else{
    NumericVector dx = diff(xSeq);
    int n = dx.size();
    std::vector<double> vpb_(n*m);
    vpb = vpb_;
    for(int i=1; i<=m; ++i) {
      double c = ySeq[i-1];      
      double d =  ySeq[i+1-1];
      IntegerVector yInd =which_C( (y>c-lambda) & (y<d+lambda));
      if (yInd.size()>0){
        NumericVector y_cd = y[yInd];
        NumericVector x_cd = x[yInd];
        NumericVector lambda_cd = lambda[yInd];
        NumericVector yMin_cd = pmax(c,y_cd-lambda_cd);
        NumericVector yMax_cd = pmin(d,y_cd+lambda_cd);
        for (int j=1; j<=n; ++j) {
          double a = xSeq[j-1];
          double b = xSeq[j+1-1];
          IntegerVector xInd = which_C((x_cd>a-lambda_cd) & (x_cd<b+lambda_cd));
          if (xInd.size()>0) {
            NumericVector x_abcd = x_cd[xInd];
            NumericVector lambda_abcd = lambda_cd[xInd];
            NumericVector xMin = pmax(a,x_abcd-lambda_abcd);
            NumericVector xMax = pmin(b,x_abcd+lambda_abcd);
            NumericVector yMin = yMin_cd[xInd];
            NumericVector yMax = yMax_cd[xInd];
            vpb[j+(i-1)*n-1] = 0.5*sum((xMax-xMin)*(yMax-yMin)*(xMin+xMax+yMin+yMax))/(dx[j-1]*dy[i-1]);
          }
          else vpb[j+(i-1)*n-1] = 0;
        }
      }
      //else vpb[1:n+(i-1)*n] <- 0
    }
    return NumericVector(vpb.begin(), vpb.end());
  }
}
