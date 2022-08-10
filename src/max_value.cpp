
#include <Rcpp.h>


// https://stackoverflow.com/questions/35047929/fastest-way-to-find-the-index-of-the-maximum-of-each-row/35048449#35048449

// [[Rcpp::export]]
Rcpp::NumericVector MaxCol(Rcpp::NumericMatrix m) {
    R_xlen_t nr = m.nrow(), nc = m.ncol(), i = 0;
    Rcpp::NumericVector result(nr);

    for (int i = 0; i < nr; i++) {
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

