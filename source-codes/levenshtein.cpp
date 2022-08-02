#include <Rcpp.h>
using namespace Rcpp;
#include <omp.h>
// [[Rcpp::plugins("openmp")]]
// [[Rcpp::plugins("cpp11")]]
#include <list>

// [[Rcpp::export]]
int levenshtein_bigmem(IntegerVector s1, IntegerVector s2) {
   int n1 = s1.size();
   int n2 = s2.size();
   IntegerMatrix D(n1+1, n2+1);
   D(0,0) = 0;
   for (int i=1; i<=n1; ++i) D(i,0) = D(i-1,0)+1;
   for (int j=1; j<=n2; ++j) D(0,j) = D(0,j-1)+1;
   for (int i=1; i<=n1; ++i) {
      for (int j=1; j<=n2; ++j) {
         int m1 = D(i-1,j-1)+(int)(s1[i-1]!=s2[j-1]);
         int m2 = D(i,j-1)+1;
         int m3 = D(i-1,j)+1;
         D(i,j) = std::min(std::min(m1, m2), m3);
      }
   }
   return D(n1, n2);
}


// #define INF (1.0/0.0)

// [[Rcpp::export]]
IntegerVector levenshtein_centroid2(IntegerVector s1, IntegerVector s2) {
   int n1 = s1.size(), n2 = s2.size();
   NumericMatrix D(n1+1, n2+1);
   IntegerMatrix T(n1+1, n2+1);
   for (int i=1; i<=n1; ++i) { // deletion
      D(i,0) = D(i-1,0)+1; T(i,0) = 4;
   }
   for (int j=1; j<=n2; ++j) { // insertion
      D(0,j) = D(0,j-1)+1; T(0,j) = 2;
   }
   for (int i=1; i<=n1; ++i) {
      for (int j=1; j<=n2; ++j) {
         T(i,j) = 0;
         if (s1[i-1]==s2[j-1])
            D(i,j) = D(i-1,j-1);
         else {
            double m1 = D(i-1,j-1)+1; // sub
            double m2 = D(i,j-1)+1;   // ins
            double m3 = D(i-1,j)+1;   // del
            D(i,j) = std::min(std::min(m1, m2), m3);
            if (D(i,j) == m1) T(i,j) |= 1;
            if (D(i,j) == m2) T(i,j) |= 2;
            if (D(i,j) == m3) T(i,j) |= 4;
         }
      }
   }

   int maxd = (int)(D(n1, n2)*0.5);
   if (maxd <= 0) return s1;

   std::list<int> l1(s1.begin(), s1.end());
   auto it1 = l1.end(); --it1;
   auto it2 = s2.end(); --it2;
   int x = n1, y = n2;
   for (int curd=0; curd < maxd; ) {
      curd += (int)(T(x,y) != 0);
      if (T(x,y) == 0) { // no change needed
         x--; y--;
         --it1; --it2;
      }
      else if (T(x,y) & 1) { // substitution
         x--; y--;
         (*(it1--)) = (*(it2--));
       }
      else if ((T(x, y) & 2) && ((!(T(x,y)&4)) || ((int)l1.size() < std::max(n1,n2)))) {
         // insertion
         y--;
         it1 = l1.insert(++it1, *(it2--)); --it1;
      }
      else { // deletion
         x--;
         it1 = l1.erase(it1); --it1;
      }
   }

   return IntegerVector(l1.begin(), l1.end());
}




int levenshtein_smallmem(int* s1, int* s2, int n1, int n2) {
   if (n1 < n2) {
      std::swap(s1, s2); // pointer swap
      std::swap(n1, n2);
   }

   int* v_cur = new int[n2+1];
   int* v_last = new int[n2+1]; // n2 <= n1
   for (int j=0; j<=n2; ++j) v_cur[j] = j;

   for (int i=1; i<=n1; ++i) {
      std::swap(v_last, v_cur); // pointer swap
      v_cur[0] = i;
      for (int j=1; j<=n2; ++j)
         v_cur[j] = std::min(std::min(
               v_last[j-1]+(int)(s1[i-1]!=s2[j-1]),
               v_cur[j-1]+1),
               v_last[j]+1);
   }

   int ret = v_cur[n2];
   delete [] v_cur;
   delete [] v_last;
   return ret;
}


// [[Rcpp::export]]
int levenshtein_smallmem(IntegerVector s1, IntegerVector s2) {
   return levenshtein_smallmem(INTEGER(s1), INTEGER(s2), LENGTH(s1), LENGTH(s2));
}


// [[Rcpp::export]]
IntegerVector levenshtein_dist_max(List Y, List X) {
   int nx = X.size(), ny = Y.size();
   if (!Rf_isVectorList((SEXP)Y)) stop("!Rf_isInteger((SEXP)Y)");
   if (!Rf_isVectorList((SEXP)X)) stop("!Rf_isInteger((SEXP)X)");
   for (int i=0; i<ny; ++i)
      if (!Rf_isInteger((SEXP)Y[i])) stop("!Rf_isInteger((SEXP)Y[i])");
   for (int j=0; j<nx; ++j)
      if (!Rf_isInteger((SEXP)X[j])) stop("!Rf_isInteger((SEXP)X[j])");

   IntegerVector out(ny);
   for (int i=0; i<ny; ++i) {
      int max_lev = 0;
      // TO DO parallel
      for (int j=0; j<nx; ++j) {
         int lev = levenshtein_smallmem(INTEGER(Y[i]), INTEGER(X[j]), LENGTH(Y[i]), LENGTH(X[j]));
         if (lev > max_lev) max_lev = lev;
      }
      out[i] = max_lev;
   }

   return out;
}



// [[Rcpp::export]]
IntegerVector levenshtein_dist_sum(List Y, List X) {
   int nx = X.size(), ny = Y.size();
   if (!Rf_isVectorList((SEXP)Y)) stop("!Rf_isInteger((SEXP)Y)");
   if (!Rf_isVectorList((SEXP)X)) stop("!Rf_isInteger((SEXP)X)");
   for (int i=0; i<ny; ++i)
      if (!Rf_isInteger((SEXP)Y[i])) stop("!Rf_isInteger((SEXP)Y[i])");
   for (int j=0; j<nx; ++j)
      if (!Rf_isInteger((SEXP)X[j])) stop("!Rf_isInteger((SEXP)X[j])");

   IntegerVector out(ny);
   for (int i=0; i<ny; ++i) {
      int sum_lev = 0;
#pragma omp parallel for schedule(dynamic) reduction(+:sum_lev) num_threads(4)
      for (int j=0; j<nx; ++j) {
         sum_lev += levenshtein_smallmem(INTEGER(Y[i]), INTEGER(X[j]), LENGTH(Y[i]), LENGTH(X[j]));
      }
      out[i] = sum_lev;
   }

   return out;
}

