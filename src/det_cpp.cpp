// [[Rcpp::depends(RcppArmadillo)]]

# include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::export]]
double det_cpp(NumericMatrix set) {
  arma::mat set_arma(set.begin(), set.nrow(), set.ncol(), false);
  return(arma::det(set_arma));
}