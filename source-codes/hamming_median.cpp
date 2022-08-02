#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::plugins("cpp11")]]
#include <unordered_map>
#include <string>

/*** R
hamming_dist_sum_R <- function(P, X) {
   apply(P, 2, function(p) { # for each p in P
      sum(apply(X, 2, function(x, p) sum(x != p), p))
   })
}
*/

// [[Rcpp::export]]
IntegerVector hamming_dist_sum(IntegerMatrix Y, IntegerMatrix X) {
   int nx = X.ncol();
   int ny = Y.ncol();
   int d = Y.nrow();
   if (X.nrow() != d) stop("X.nrow() != Y.nrow()");

   IntegerVector out(ny);
   for (int i=0; i<ny; ++i) {
      int sum_hamming = 0;
      for (int j=0; j<nx; ++j) {
         // Hamming distance between Y[,i] and X[,j]
         for (int k=0; k<d; ++k) sum_hamming += (int)(Y(k,i) != X(k,j));
      }
      out[i] = sum_hamming;
   }

   return out;
}


// [[Rcpp::export]]
IntegerVector hamming_median(IntegerMatrix X) {
   int n = X.ncol();
   int d = X.nrow();
   IntegerVector out(d);
   for (int i=0; i<d; ++i) {
      std::unordered_map<int, int> hashtable;
      for (int j=0; j<n; ++j) /* ints are default-constructed as 0 */
         hashtable[X(i,j)]++; /* count the number of occurrences of each letter */
      int max = 0, argmax = -1;
      for (auto it=hashtable.cbegin(); it != hashtable.cend(); ++it)
         if (max < (*it).second) { // find a most frequently
            max = (*it).second;    // occurring letter
            argmax = (*it).first;
         }
      out[i] = argmax;
   }
   return out;
}
