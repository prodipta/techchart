#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;


double square(double x){
  return (x*x);
}
long argmax(NumericVector x){
  long k=0;
  for(int i=1;i<x.size();i++){
    if(x[i]>x[k]){
      k = i;
    }
  }
  return k;
}
long argmin(NumericVector x){
  long k=0;
  for(int i=1;i<x.size();i++){
    if(x[i]<x[k]){
      k = i;
    }
  }
  return k;
}
NumericVector head(NumericVector x, long k){
  NumericVector y(k);
  if(x.size()<k){
    return x;
  }
  for(int i=0;i<k;i++){
    y[i] = x[i];
  }
  return y;
}
NumericVector tail(NumericVector x, long k){
  long n = x.size();
  NumericVector y(n-k);
  if(n<k){
    return x;
  }

  for(int i=k;i<n;i++){
    y[i-k] = x[i];
  }
  return y;
}
NumericVector ols_estimate(NumericVector x, NumericVector y){
  long n;
  double sumx, sumx2, sumxy, sumy, sumy2;
  double denom, a, b;
  NumericVector ret(2);

  ret[0]=ret[1]=-1;

  if(x.size() != y.size()){
    return ret;
  } else{
    n = x.size();
  }

  sumx=0;sumx2=0;sumxy=0;sumy=0;sumy2=0;

  for(int i=0; i<n; i++){
    sumx = sumx + x[i];
    sumx2 = sumx2 + square(x[i]);
    sumxy = sumxy + x[i]*y[i];
    sumy = sumy + y[i];
    sumy2 = sumy2 + square(y[i]);
  }

  denom = sumx2 - square(sumx)/n;

  if(denom==0){
    return ret;
  }

  b = (sumxy - sumx*sumy/n)/denom;
  a = sumy/n - b*sumx/n;

  ret[0] = a; ret[1] = b;

  return ret;

}
double ols_llk (NumericVector x, NumericVector y, long k){
  long n;
  double r, residsq, variance, ret;
  NumericVector ols(2), ols2(2);

  if(x.size() != y.size()){
    return -1;
  } else{
    n = x.size();
  }

  NumericVector x1(k), x2(n-k), y1(k), y2(n-k);

  if(k==0){
    ols = ols_estimate(x,y);
  } else{
    x1 = head(x,k); y1 = head(y,k);
    x2 = tail(x,k); y2 = tail(y,k);
    ols = ols_estimate(x1,y1);
    ols2 = ols_estimate(x2,y2);
  }

  if(ols[0]==-1 && ols[1]==-1){
    return -1;
  }
  if(ols2[0]==-1 && ols2[1]==-1 && k!=0){
    return -1;
  }

  r=residsq=variance=ret=0;

  if(k==0){
    for(int i=0; i<n; i++){
      r = y[i] - (ols[0] + ols[1]*x[i]);
      residsq = residsq + square(r);
    }
  } else{
    for(int i=0; i<n; i++){
      if(i < k){
        r = y[i] - (ols[0] + ols[1]*x[i]);
      } else{
        r = y[i] - (ols2[0] + ols2[1]*x[i]);
      }
      residsq = residsq + square(r);
    }
  }


  variance = residsq/n;
  if(variance == 0){
    return -1;
  }

  //Rcout << "variance: " << y1<< std::endl;

  ret = double(n)/2*log(2*PI) + double(n)/2*log(variance) + double(n)/2;

  return ret;
}
long cpt_trend_AMOC(NumericVector x, NumericVector y, long minseglen,
                    double penalty){
  long n, nn;
  double stat, critical_stat;

  if(x.size() != y.size()){
    return -1;
  } else{
    n = x.size();
  }

  if(n <2*minseglen){
    return -1;
  }

  NumericVector llc(n-2*(minseglen-1));

  for(int i=0; i<llc.size(); i++){
    llc[i] = ols_llk(x,y,minseglen+i);
    if(llc[i]==-1){
      return -1;
    }
  }

  nn = argmin(llc);
  stat = ols_llk(x,y,0) - llc[nn];
  critical_stat = (2 + square(2*nn/n-1))*log(n)*penalty;

  if(stat < critical_stat){
    nn = -1;
  } else{
    nn = nn+minseglen;
  }

  return nn;
}
long cpt_trend_binary(NumericVector x, NumericVector y, IntegerVector cpts,
                      int Q, long minseglen, double penalty,
                      IntegerVector k, long last_pt){
  long pt, pt1, pt2;

  pt=pt1=pt2=0;
  if(k[1]>=cpts.size()){
    return -1;
  }

  if(Q<2){
    pt = cpt_trend_AMOC(x,y,minseglen,penalty);
    if(pt<1)return-1;
    cpts[k[1]] = pt + last_pt;
    k[1] = k[1]+1;
    return 1;
  }
  else{
    pt = cpt_trend_AMOC(x,y,minseglen,penalty);
    if(pt<1)return-1;
    cpts[k[1]] = pt + last_pt;
    k[1] = k[1]+1;
    pt1 = cpt_trend_binary(head(x,pt),head(y,pt),cpts,Q-2,minseglen,penalty,k,last_pt);
    pt2 = cpt_trend_binary(tail(x,pt),tail(y,pt),cpts,Q-2,minseglen,penalty,k,pt+last_pt);
  }
  return 1;
}



// [[Rcpp::export]]
IntegerVector cpt_trend(NumericVector x,NumericVector y, int Q,
                        long minseglen, double penalty){
  int n;
  IntegerVector k(1);

  n = 2*log(Q+1)/log(2) -1;
  if(n%2 !=0)n = n+1;

  k[1]=0;
  IntegerVector cpts(Q);

  cpt_trend_binary(x,y,cpts,n,minseglen,penalty,k,0);

  return(cpts);
}

