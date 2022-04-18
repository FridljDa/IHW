/*
 * 

#include <Rcpp.h>
using namespace Rcpp;


int findInterval(double x, NumericVector vec) {
  // https://www.cplusplus.com/reference/algorithm/upper_bound/ 
  // https://en.cppreference.com/w/cpp/algorithm/upper_bound
  auto upper = std::upper_bound(vec.begin(), vec.end(), x);
  int i = std::distance(vec.begin(), upper);
  return i;
}

// https://stackoverflow.com/questions/321068/returning-multiple-values-from-a-c-function
// https://github.com/RcppCore/Rcpp/issues/974
// [[Rcpp::export]]
NumericVector  invert_subgradient(NumericVector slopeknots, NumericVector xknotsl, double lambda) {
  if (lambda > slopeknots[0]) {
    NumericVector res = { 0,0 };
    return res;
  } else if (lambda <= 0) {
    NumericVector res = { 1,1 };
    return res;
  }else{
    int idx = findInterval(-lambda, -slopeknots);
    double l;
    if(slopeknots[idx-1] == lambda){
       l = xknotsl[idx-1];
    }else{
       l = xknotsl[idx];
    }
    double u = xknotsl[idx];
    NumericVector res = {l,u};
    return res;
  }
}

// [[Rcpp::export]]
DataFrame  invert_subgradient_vec(Rcpp::List slopeknots, Rcpp::List xknotsl, double lambda) {
  int groups = slopeknots.size();
  NumericVector slopeknots_tmp;
  NumericVector xknotsl_tmp;
  //Rcpp::List ts;
  NumericVector ts_L;
  NumericVector ts_R;
  NumericVector tmp;
  for(int i = 0; i < groups; ++i){
    slopeknots_tmp = slopeknots[i];
    xknotsl_tmp = xknotsl[i];
    tmp=invert_subgradient(slopeknots_tmp, xknotsl_tmp, lambda);
    //ts.push_back(tmp);
    ts_L.push_back(tmp[0]);
    ts_R.push_back(tmp[1]);
  }
  DataFrame res = DataFrame::create(Named("ts_L")=ts_L,Named("ts_R")=ts_R);
  return res;
}


//std::vector<double>  grenander_pdf(NumericVector slopeknots, NumericVector xknotsl, double x, double tolerance = 1.0e-10) {
//  if (x > 1 | x < 0) {
//    std::vector<double> res = { 0,0 };
//    return res;
//    }else{
//    //tolerance to handle numerical rounding and discontinuity points
//    int idx1 = findInterval(x, xknotsl - tolerance);//-tolerance
//    int idx2 = findInterval(x, xknotsl + tolerance);//+ tolerance
//    double pdf_min = slopeknots[idx1];
//    double pdf_max = slopeknots[idx2];
//    std::vector<double> res = {pdf_min, pdf_max};
//    return res;
//      }
//  }
 

 //double grenander_cdf(NumericVector xknotsl, NumericVector yknots, double x){
 //  if (x > 1) {
 //  return 1;
 //  }else if(x <0){
 //    return 0;
 //  }else{
 //    int idx = findInterval(x, xknotsl);
 //    double loc_L = xknotsl[idx-1];
 //    double loc_U = xknotsl[idx];
    
    //    double t = (x - loc_L) / (loc_U - loc_L);
    
    //    double f_l = yknots[idx-1];
    //double f_r = yknots[idx];
    
    //double f = (1 - t) * f_l + t * f_r;
    //return f;
    //}
    //}


    //NumericVector  grenander_cdf_vec(Rcpp::List xknotsl,Rcpp::List yknots, NumericVector xs) {
    //int groups = yknots.size();
    //NumericVector xknotsl_tmp;
    //NumericVector yknots_tmp;
    //double tmp;
    //double x;
    //NumericVector fs;
    //for(int i = 0; i < groups; ++i){
    //xknotsl_tmp = xknotsl[i];
    //yknots_tmp = yknots[i];
    //x = xs[i];
    //tmp=grenander_cdf(xknotsl_tmp, yknots_tmp, x);
    //fs.push_back(tmp);
    //  }
    //  return fs;
    //}

    //double Fdr(NumericVector ts, NumericVector Fs, NumericVector m_groups){
    //NumericVector num_vec = ts*m_groups;
    //double num = std::accumulate(num_vec.begin(), num_vec.end(), 0.0);
    //NumericVector denom_vec = Fs*m_groups;
    //double denom = std::accumulate(denom_vec.begin(), denom_vec.end(), 0.0);
    //if(num == 0){
    //return 0;
    //}else{ //if(denom != 0){
    //return num/denom;
    //}//else{
  //  return arma::math::inf();
  //}
  //}

  //double Fdr_linearized(NumericVector ts, NumericVector Fs, NumericVector m_groups, double alpha){
  //NumericVector num_vec = ts*m_groups;
  //double num = std::accumulate(num_vec.begin(), num_vec.end(), 0.0);
  //NumericVector denom_vec = Fs*m_groups;
  //double denom = std::accumulate(denom_vec.begin(), denom_vec.end(), 0.0);
  //double res = num - denom* alpha;
  //return res;
  //}


  //NumericVector lagrange_balance(double lambda,Rcpp::List slopeknotsvec,Rcpp::List xknotslvec,Rcpp::List yknotsvec, NumericVector m_groups, double alpha, bool linearized){
  //DataFrame ts = invert_subgradient_vec(slopeknotsvec, xknotslvec, lambda);
  //NumericVector ts_L = ts["ts_L"];
  //NumericVector ts_R = ts["ts_R"]; //TODO consistency ts_U ts_R
  
  //NumericVector Fs_L = grenander_cdf_vec(xknotslvec, yknotsvec, ts_L);
  //NumericVector  Fs_R = grenander_cdf_vec(xknotslvec, yknotsvec, ts_R);
  //double Fdr_L;
  //double Fdr_R;
  //if (linearized) {
  //   Fdr_L = Fdr_linearized(ts_L, Fs_L, m_groups, alpha);
  //   Fdr_R = Fdr_linearized(ts_R, Fs_R, m_groups, alpha);
  //} else {
  //   Fdr_L = Fdr(ts_L, Fs_L, m_groups);
  //   Fdr_R = Fdr(ts_R, Fs_R, m_groups);
  //}
  
  //NumericVector res =
  //   NumericVector::create( _["Fdr_L"]=Fdr_L, _["Fdr_R"]=Fdr_R);
  // return res; 
  //} 
*/
