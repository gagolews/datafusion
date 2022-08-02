#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::plugins("cpp11")]]
#include <string>

/*** R
hamming_dist_max_R <- function(P, X) {
   apply(P, 2, function(p) { # for each p in P
      max(apply(X, 2, function(x, p) sum(x != p), p))
   })
}
*/

// [[Rcpp::export]]
IntegerVector hamming_dist_max(IntegerMatrix Y, IntegerMatrix X) {
   int nx = X.ncol(), ny = Y.ncol(), d = Y.nrow();
   if (X.nrow() != d) stop("X.nrow() != Y.nrow()");

   IntegerVector out(ny);
   for (int i=0; i<ny; ++i) {
      int max_hamming = 0;
      for (int j=0; j<nx; ++j) {
         // Hamming distance between Y[,i] and X[,j]
         int h = 0;
         for (int k=0; k<d; ++k) h += (int)(Y(k,i) != X(k,j));
         if (h > max_hamming) max_hamming = h;
      }
      out[i] = max_hamming;
   }

   return out;
}
