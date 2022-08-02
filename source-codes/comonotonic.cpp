#include <Rcpp.h>
#include <algorithm>
using namespace Rcpp;
// [[Rcpp::plugins("cpp11")]]


// [[Rcpp::export]]
bool is_comonotonic_old(NumericVector x, NumericVector y) {
   int n = x.size();
   if (y.size() != n) stop("x and y are not of the same length");
   for (int i=0; i<n-1; ++i)
      for (int j=(i+1); j<n; ++j)
         // ^ denotes bitwise XOR
         if (x[i] != x[j] && y[i] != y[j] && bool(x[i] > x[j]) ^ bool(y[i] > y[j]))
            return false;
   return true;
}


bool aux_comonotonic(const double* x, const double* y, int* o, int a, int b) {
   if (a >= b) return true;

   std::swap(o[a+(int)(unif_rand()*(b-a+1))], o[a]); // random pivot element
   for (int i=a+1; i<=b; ++i)
      if (x[o[a]] != x[o[i]] && y[o[a]] != y[o[i]] && bool(x[o[a]] > x[o[i]]) ^ bool(y[o[a]] > y[o[i]]))
         return false;
   // here surely for i=a+1,...,b:
   // x[o[a]] == x[o[i]] AND y[o[a]] ?? y[o[i]] OR
   // x[o[a]] ?? x[o[i]] AND y[o[a]] == y[o[i]] OR
   // x[o[a]] <  x[o[i]] AND y[o[a]] <  y[o[i]] OR
   // x[o[a]] >  x[o[i]] AND y[o[a]] >  y[o[i]]

   // cases excluded:
   // x[o[a]] > x[o[i]] AND x[o[a]] < x[o[i]] OR
   // x[o[a]] < x[o[i]] AND x[o[a]] > x[o[i]]
   int a2 = a+1;
   int b2 = b;
   while (a2 <= b2) {
//      if (x[o[a]] <  x[o[a2]] && y[o[a]] < y[o[a2]]) a2++;
//      else if (x[o[a]] == x[o[a2]] && y[o[a]] < y[o[a2]]) a2++;
//      else if (x[o[a]] <  x[o[a2]] && y[o[a]] == y[o[a2]]) a2++;
      while (a2 <= b2 && x[o[a]] <=  x[o[a2]] && y[o[a]] <= y[o[a2]]) a2++;
      while (a2 <= b2 && x[o[a]] >=  x[o[b2]] && y[o[a]] >= y[o[b2]]) b2--;
      if (a2 < b2) { std::swap(o[a2], o[b2]); a2++; b2--; }
   }
   if (a2 -1 != b2) throw std::exception();
   std::swap(o[a], o[a2-1]);
   return aux_comonotonic(x, y, o, a, a2-2) && aux_comonotonic(x, y, o, a2, b);
}


// [[Rcpp::export]]
bool is_comonotonic3(NumericVector x, NumericVector y) {
   int n = x.size();
   if (y.size() != n) stop("x and y are not of the same length");

   // recall that array elements in C++ are numbered from 0
   std::vector<int> o(n); for (int i=0; i<n; ++i) o[i] = i; // o = (0,1,...,n-1)

   return aux_comonotonic(REAL(x), REAL(y), o.data(), 0, n-1);
}


struct Comparer {
   const double* v;
   Comparer(const double* _v) { v = _v; }
   bool operator()(const int& i, const int& j) {
      // returns true if the first argument is less than
      // (i.e. is ordered before) the second.
      return (v[i] < v[j]);
   }
};


// [[Rcpp::export]]
bool is_comonotonic(NumericVector x, NumericVector y) {
   int n = x.size();
   if (y.size() != n) stop("x and y are not of the same length");

   // recall that array elements in C++ are numbered from 0
   std::vector<int> o(n); for (int i=0; i<n; ++i) o[i] = i; // o = (0,1,...,n-1)

   Comparer lt_x(REAL(x));
   std::sort(o.begin(), o.end(), lt_x);
   // now o is an ordering permutation of x

   Comparer lt_y(REAL(y));
   int i1 = 0;
   while (i1 < n) { // now search for the longest subsequence consisting of equal x's
      int i2 = i1+1;
      while (i2 < n && x[o[i1]] == x[o[i2]]) ++i2;
      // i2 == n or i2 is first such that x[o[i1]] < x[o[i2]]
      if (i2-i1 > 1) std::sort(o.begin()+i1, o.begin()+i2, lt_y);
      if (i1 > 0 && y[o[i1-1]] > y[o[i1]]) return false;
      i1 = i2;
   }
//   for (int i=0; i<n-1; ++i) // we surely have x[o[i]] <= x[o[i+1]]
//      if (y[o[i]] > y[o[i+1]]) return false;
   // as a by-product, (o[0]+1, o[1]+1, ..., o[n-1]+1)
   // is a permutation that orders both x and y
   return true;
}


// Comparer used in std::sort:
// determines whether (x[o[i]], y[o[i]]) <= (x[o[j]], y[o[j]])
struct Comparer2 {
   const double* x;
   const double* y;
   Comparer2(const double* _x, const double* _y) {
      x = _x; y = _y;
   }
   bool operator()(const int& i, const int& j) {
      // ^ denotes bitwise XOR
      if (x[i] != x[j] && y[i] != y[j] && bool(x[i] > x[j]) ^ bool(y[i] > y[j]))
         throw false; // finish sorting immediately
      else
         return (x[i] < x[j] && y[i] <= y[j]) || (x[i] == x[j] && y[i] < y[j]);
         // returns true if i should be ordered before j
   }
};


// [[Rcpp::export]]
bool is_comonotonic2(NumericVector x, NumericVector y) {
   int n = x.size();
   if (y.size() != n) stop("x and y are not of the same length");

   // recall that array elements in C++ are numbered from 0
   std::vector<int> o(n); for (int i=0; i<n; ++i) o[i] = i; // o = (0,1,...,n-1)

   try {
      Comparer2 leq(REAL(x), REAL(y));
      std::sort(o.begin(), o.end(), leq);
      // if no exception thrown => comonotonic;
      // as a by-product, (o[0]+1, o[1]+1, ..., o[n-1]+1)
      // is a permutation that orders both x and y
      return true;
   } catch (bool) { // exception thrown => definitely not comonotonic
      return false;
   }
}

