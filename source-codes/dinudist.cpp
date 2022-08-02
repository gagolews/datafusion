#include <Rcpp.h>
#include <algorithm>
using namespace Rcpp;

struct Comparer {
   const int* v;
   Comparer(const int* _v) { v = _v; }
   bool operator()(const int& i, const int& j) const { return v[i] < v[j]; }
};


// [[Rcpp::export]]
double dinudist(IntegerVector x, IntegerVector y) {
   int nx = x.size();
   std::vector<int> ox(nx); // ordering permutation of x
   for (int i=0; i<nx; ++i) ox[i] = i;
   std::stable_sort(ox.begin(), ox.end(), Comparer(INTEGER(x)));

   int ny = y.size();
   std::vector<int> oy(ny); // ordering permutation of y
   for (int i=0; i<ny; ++i) oy[i] = i;
   std::stable_sort(oy.begin(), oy.end(), Comparer(INTEGER(y)));

   double d = 0.0;
   int ix = 0, iy = 0;
   while (ix < nx && iy < ny) {
      if (x[ox[ix]] == y[oy[iy]])
         d += std::abs((ox[ix++]+1) - (oy[iy++]+1));
      else if (x[ox[ix]] < y[oy[iy]])
         d += std::abs((ox[ix++]+1) - 0);
      else
         d += std::abs(0 - (oy[iy++]+1));
   }
   while (ix < nx) d += std::abs((ox[ix++]+1) - 0);
   while (iy < ny) d += std::abs(0 - (oy[iy++]+1));

   return d;
}
