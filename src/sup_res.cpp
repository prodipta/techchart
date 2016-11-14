#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericVector timesTwo(NumericVector x) {
  return x * 2;
}

double max(NumericVector x){
  int n = x.size();
  double xmax = x[0];
  for(int i=0;i<n;i++){
    if(x[i]>xmax){
      xmax = x[i];
    }
  }
  return xmax;
}
double min(NumericVector x){
  int n = x.size();
  double xmin = x[0];
  for(int i=0;i<n;i++){
    if(x[i]<xmin){
      xmin = x[i];
    }
  }
  return xmin;
}
NumericVector subset(NumericVector x, int i, int j){
  NumericVector y(j-i+1);
  for(int k=0; k<y.size(); k++)y[k]=x[i+k];
  return y;
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

// [[Rcpp::export]]
IntegerMatrix findminima(NumericVector xmin, NumericVector xmax,
                         NumericVector threshold){
  int n = xmin.size();
  int imin = 0;
  int k = 0;
  IntegerMatrix y(n,2);

  for(int i=0; i<n; i++){

    if(xmin[i] < xmin[imin])imin = i;
    if(xmax[i]/xmin[imin] > threshold[i]){
      if(max(subset(xmax,(y.row(k))[0],imin))/xmin[imin]> threshold[i]){
        (y.row(k))[0] = imin;
        (y.row(k+1))[0] = imin;
        (y.row(k))[1] = -1;
        k = k+1;
      }
      imin = i+1;
    }
  }
  return(y);
}

// [[Rcpp::export]]
IntegerMatrix findmaxima(NumericVector xmin, NumericVector xmax,
                         NumericVector threshold){
  int n = xmin.size();
  int imax = 0;
  int k = 0;
  IntegerMatrix y(n,2);

  for(int i=0; i<n; i++){

    if(xmax[i] > xmax[imax])imax = i;
    if(xmax[imax]/xmin[i] > threshold[i]){
      if(xmax[imax]/min(subset(xmin,(y.row(k))[0],imax))> threshold[i]){
        (y.row(k))[0] = imax;
        (y.row(k+1))[0] = imax;
        (y.row(k))[1] = 1;
        k = k+1;
      }
      imax = i+1;
    }
  }
  return(y);
}

// [[Rcpp::export]]
IntegerVector sortoptimaposition(IntegerVector pos, IntegerVector sign,
                                 NumericVector value){
  int n = pos.size();
  for(int i=1; i<n; i++){
    if(pos[i] == pos[i-1]){
      if(sign[i]==1){
        if(value[i]>value[i-1]){
          sign[i-1] = 0;
        } else{
          sign[i] = 0;
        }
      }else{
        if(value[i]<value[i-1]){
          sign[i-1] = 0;
        } else{
          sign[i] = 0;
        }
      }
    }
  }
  return sign;
}

// [[Rcpp::export]]
IntegerVector sortoptimasign(IntegerVector pos, IntegerVector sign,
                                 NumericVector value){
  int n = pos.size();
  for(int i=1; i<n; i++){
    if(sign[i] == sign[i-1]){
      if(sign[i]==1){
        if(value[i]>value[i-1]){
          sign[i-1] = 0;
        } else{
          sign[i] = 0;
        }
      }else{
        if(value[i]<value[i-1]){
          sign[i-1] = 0;
        } else{
          sign[i] = 0;
        }
      }
    }
  }
  return sign;
}

// [[Rcpp::export]]
bool checkoptimasign(IntegerVector sign){
  int n = sign.size();
  bool ret = true;
  for(int i=0; i<n; i++){
    if(sign[i]==sign[i-1]){
      ret = false;
      break;
    }
  }
  return ret;
}

// [[Rcpp::export]]
bool checkoptimapos(IntegerVector pos){
  int n = pos.size();
  bool ret = true;
  for(int i=0; i<n; i++){
    if(pos[i]==pos[i-1]){
      ret = false;
      break;
    }
  }
  return ret;
}
