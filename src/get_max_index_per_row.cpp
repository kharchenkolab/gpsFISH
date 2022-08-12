#include <Rcpp.h>
using namespace Rcpp;

//' Find the index of the maximum value of each row in a matrix.
//'
//' @param m Numeric matrix
//'
//' @return A numeric vector containing the index for each row.
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector MaxCol(Rcpp::NumericMatrix m) {
  R_xlen_t nr = m.nrow(), nc = m.ncol(), i = 0;
  Rcpp::NumericVector result(nr);

  for ( ; i < nr; i++) {
    double current = m(i, 0);
    R_xlen_t idx = 0, j = 1;
    for ( ; j < nc; j++) {
      if (m(i, j) > current) {
        current = m(i, j);
        idx = j;
      }
    }
    result[i] = idx + 1;
  }
  return result;
}

