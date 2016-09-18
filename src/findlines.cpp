#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;


double pi(){
  return(3.141593);
}
double variance(NumericVector x){
  int n = x.size();
  double total=0, var=0, avg;
  for(int i=0; i<n; i++){
    total = total + x[i];
  }
  avg = total/n;
  for(int i =0; i<n; i++){
    total = total + (x[i]-avg)*(x[i]-avg);
  }
  var = total/n;
  return var;
}
int findmatch(double x, NumericVector r){
  int n = r.size();
  for(int i=1; i<n; i++){
    if(x>=r[i-1] && x<r[i])
      return(i-1);
  }
  return(n-1);
}
IntegerMatrix localmaxima(IntegerMatrix mat, IntegerMatrix minmat,
                          IntegerMatrix maxmat, int nrow, int ncol, int s){
  int k=0;
  IntegerMatrix lines(nrow*ncol,5);

  // check top left
  if(mat.row(0)[0] > mat.row(1)[0] && mat.row(0)[0] > mat.row(0)[1] &&
     mat.row(0)[0] > mat.row(1)[1] && mat.row(0)[0] > s){
    lines.row(k)[0]=0;
    lines.row(k)[1]=0;
    lines.row(k)[2]=mat.row(0)[0];
    lines.row(k)[3]=minmat.row(0)[0];
    lines.row(k)[4]=maxmat.row(0)[0];
    k=k+1;
  }
  // check bottom left
  if(mat.row(nrow-1)[0] > mat.row(nrow-2)[0] && mat.row(nrow-1)[0] > mat.row(nrow-1)[1] &&
     mat.row(nrow-1)[0] > mat.row(nrow-2)[1] && mat.row(nrow-1)[0] > s){
    lines.row(k)[0]=nrow-1;
    lines.row(k)[1]=0;
    lines.row(k)[2]=mat.row(nrow-1)[0];
    lines.row(k)[3]=minmat.row(nrow-1)[0];
    lines.row(k)[4]=maxmat.row(nrow-1)[0];
    k=k+1;
  }
  // check bottom right
  if(mat.row(nrow-1)[ncol-1] > mat.row(nrow-1)[ncol-2] && mat.row(nrow-1)[ncol-1] > mat.row(nrow-2)[ncol-1] &&
     mat.row(nrow-1)[ncol-1] > mat.row(nrow-2)[ncol-2] && mat.row(nrow-1)[ncol-1] > s){
    lines.row(k)[0]=nrow-1;
    lines.row(k)[1]=ncol-1;
    lines.row(k)[2]=mat.row(nrow-1)[ncol-1];
    lines.row(k)[3]=minmat.row(nrow-1)[ncol-1];
    lines.row(k)[4]=maxmat.row(nrow-1)[ncol-1];
    k=k+1;
  }
  // check top left
  if(mat.row(0)[ncol-1] > mat.row(0)[ncol-2] && mat.row(0)[ncol-1] > mat.row(1)[ncol-1] &&
     mat.row(0)[ncol-1] > mat.row(1)[ncol-2] && mat.row(0)[ncol-1] > s){
    lines.row(k)[0]=0;
    lines.row(k)[1]=ncol-1;
    lines.row(k)[2]=mat.row(0)[ncol-1];
    lines.row(k)[3]=minmat.row(0)[ncol-1];
    lines.row(k)[4]=maxmat.row(0)[ncol-1];
    k=k+1;
  }

  for(int i=1; i<nrow-1;i++){
    for(int j=1; j<ncol-1; j++){
      if(mat.row(i)[j] > mat.row(i-1)[j] && mat.row(i)[j] > mat.row(i+1)[j] &&
         mat.row(i)[j] > mat.row(i)[j-1] && mat.row(i)[j] > mat.row(i)[j+1] &&
         mat.row(i)[j] > s){
        lines.row(k)[0]=i;
        lines.row(k)[1]=j;
        lines.row(k)[2]=mat.row(i)[j];
        lines.row(k)[3]=minmat.row(i)[j];
        lines.row(k)[4]=maxmat.row(i)[j];
        k=k+1;
      }
    }
  }
  return(lines);
}
NumericVector envelopescore(IntegerMatrix lines, int n1, NumericVector x,
                            NumericVector y, NumericVector rbucket,
                            NumericVector abucket, int flag){

  int n2, start, end, count=0;
  double dist, r, theta, rquanta, mse=0;
  NumericVector score(n1);

  n2 = x.size();
  rquanta = rbucket[2] - rbucket[1];

  for(int i=0; i<n1; i++){
    start=lines.row(i)[3];
    end=lines.row(i)[4];
    //length = end-start+1;
    for(int j=start;j<(end+1); j++){
      if(j >= n2)break;
      r = rbucket[lines.row(i)[0]];
      theta = abucket[lines.row(i)[1]];
      dist = y[j] - (r/sin(theta*pi()/180) - x[j]/tan(theta*pi()/180));
      if(flag==-1){
        if(dist+rquanta<0){
          mse = mse+dist*dist;
          count++;
        }
      }
      if(flag==1){
        if(dist-rquanta>0){
          mse = mse+dist*dist;
          count++;
        }
      }
    }

    if(count!=0)score[i] = (mse/count);
    mse = 0;
    count =0;
  }
  return(score);
}
NumericVector fitscore(IntegerMatrix lines, int n1, NumericVector x,
                       NumericVector y, NumericVector rbucket,
                       NumericVector abucket, int flag){

  int n2, start, end, length;
  double dist, r, theta, mse=0;
  NumericVector score(n1);

  n2 = x.size();

  for(int i=0; i<n1; i++){
    start=lines.row(i)[3];
    end=lines.row(i)[4];
    length = end-start+1;
    for(int j=start;j<(end+1); j++){
      if(j >= n2)break;
      r = rbucket[lines.row(i)[0]];
      theta = abucket[lines.row(i)[1]];
      dist = y[j] - (r/sin(theta*pi()/180) - x[j]/tan(theta*pi()/180));
      mse = mse + dist*dist;
    }
    score[i] = 1 - (mse/length)/variance(y);
    mse=0;
  }
  return(score);
}

// [[Rcpp::export]]
DataFrame houghtransform(NumericVector x1, NumericVector y1,int flag,
                         NumericVector rbucket, NumericVector abucket){

  int k, n = x1.size();
  int nA = abucket.size(), nR = rbucket.size();
  double r;
  IntegerMatrix accumulator(nR,nA);
  IntegerMatrix accumulatorMin(nR,nA);
  IntegerMatrix accumulatorMax(nR,nA);
  IntegerMatrix lines(nR*nA,5);
  NumericVector scores(nR*nA), fit(nR*nA);

  for(int i=0; i< nR*nA; i++){
    accumulatorMax[i] = -1;
  }

  for(int i=0; i<n; i++){
    for(int j=0; j<nA; j++){
      r = x1[i]*cos(abucket[j]*pi()/180) + y1[i]*sin(abucket[j]*pi()/180);
      k = findmatch(r, rbucket);
      accumulator.row(k)[j] = accumulator.row(k)[j] + 1;
      if(accumulatorMin.row(k)[j]==0){
        accumulatorMin.row(k)[j] = i;
      }
      if(accumulatorMax.row(k)[j]<i){
        accumulatorMax.row(k)[j] = i;
      }
    }
  }
  lines = localmaxima(accumulator,accumulatorMin, accumulatorMax,nR,nA,2);
  scores = envelopescore(lines,nR*nA,x1,y1,rbucket,abucket,flag);
  fit = fitscore(lines,nR*nA,x1,y1,rbucket,abucket,flag);

  return DataFrame::create(_["r"]=lines.column(0), _["theta"]=lines.column(1),
                           _["strength"]=lines.column(2), _["start"]=lines.column(3),
                           _["end"]=lines.column(4), _["score"]=scores,
                           _["fit"]=fit);
}

