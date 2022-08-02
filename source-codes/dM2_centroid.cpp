/* Cena, Gagolewski, AGOP 2015 */

#include <iostream>
#include <iomanip>

#ifdef NDEBUG
//#undef NDEBUG
#endif

#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::plugins("cpp11")]]
#include <deque>

// [[Rcpp::export]]
double dM2_dist(List X, NumericVector y, int ny, double p, double r) {
   int l = X.size();
   double dist = 0.0;
   for (int i=0; i<l; ++i) {
      NumericVector x(X[i]);
      int nx = x.size();
      int min_nx_ny = std::min(nx, ny);
      for (int j=0; j<min_nx_ny;  ++j) dist += (x[j]-y[j])*(x[j]-y[j]);
      for (int j=min_nx_ny; j<nx; ++j) dist += x[j]*x[j];
      for (int j=min_nx_ny; j<ny; ++j) dist += y[j]*y[j];
      dist += p*abs(std::pow(nx, r)-std::pow(ny, r));
   }
   return dist;
}


int calculate_max_vector_length(List X) {
   int m = 0;
   for (auto it=X.begin(); it != X.end(); ++it) {
      NumericVector x(*it);
      m = std::max(m, x.size());
   }
   return m;
}


NumericVector calculate_xtilde(List X, int m) {
   NumericVector xtilde(m);  // filled with 0s
   for (auto it=X.begin(); it != X.end(); ++it) {
      NumericVector x(*it);
      int nx = x.size();
      for (int j=0; j<nx; ++j)
         xtilde[j] += x[j];
   }
   return xtilde;
}

// [[Rcpp::export]]
NumericVector dM2_centroid(List X, double p, double r) {
   int l = X.size();
   int m = calculate_max_vector_length(X);
   NumericVector xtilde = calculate_xtilde(X, m);

#ifndef NDEBUG
   Rcout << std::setw(2) <<  " " << "\t" << std::setw(7) << std::setprecision(2)  << std::fixed << "xtilde=" << "\t|";
   for (int j=0; j<m; ++j) Rcout << "\t" << std::setw(7) << std::fixed << std::setprecision(3) << xtilde[j]/l;
   Rcout << std::endl;
   Rcout << "-- ------- | --------------------------------------------\n";
   Rcout << std::setw(2) <<  "n" << "\t" << std::setw(7) << "dist";
   Rcout << "\t|" << "   y1      y2      y3      ..." << std::endl;
#endif

   std::deque< std::pair<int, int> > part; // linked list (a stack)
   NumericVector y(m);
   NumericVector best_y = NumericVector(0);
   double best_dist = INFINITY;

   for (int n=1; n<=m; ++n) {
      // array elements in C++ are numbered from 0
      part.push_front( std::pair<int, int>(n-1, n-1) );
      y[n-1] = xtilde[n-1]/l;

      auto it = part.begin();
      while (it+1 != part.end() && y[(*it).first] > y[(*(it+1)).second]+1e-8) { // merge
         int p1 = (*it).second-(*it).first+1;
         int p2 = (*(it+1)).second-(*(it+1)).first+1;
         y[(*it).second] = (y[(*it).second]*p1 + y[(*(it+1)).second]*p2)/(p1+p2);
         for (int j = (*it).second-1; j>=(*(it+1)).first; --j) y[j] = y[(*it).second];
         (*(it+1)).second = (*it).second;
         it = part.erase(it); // erases current and moves forward (pop stack)
      }

      double cur_dist = dM2_dist(X, y, n, p, r);
      if (cur_dist < best_dist) {
         best_dist = cur_dist;
         best_y = NumericVector(y.begin(), y.begin()+n);
      }

#ifndef NDEBUG
         Rcout << std::setw(2) <<  n << "\t" << std::setw(7) << std::setprecision(2)  << std::fixed << cur_dist << "\t|";
         for (int j=0; j<n; ++j) Rcout << "\t" << std::setw(7) << std::fixed << std::setprecision(3) << y[j];
         Rcout << std::endl;
#endif
         {for (int j=0; j<n-1; ++j) if (y[j] +1e-8< y[j+1]) stop("UNSORTED!");
         auto it = part.begin();
         while (it != part.end()) {
            if ((*it).first != (*it).second) {
               double sumall = 0.0;
               for (int i=(*it).first; i<=(*it).second; ++i)
                  sumall += xtilde[i];
               double sum = 0.0;
               for (int i=(*it).first; i<(*it).second; ++i) {
                  sum += xtilde[i];
                  if ((i-(*it).first+1)*sumall/((*it).second-(*it).first+1) - sum <1e-8) stop("LAMBDAs FAILED!");
               }
            }
            ++it;
         }}
#ifndef NDEBUG
#endif
   }
   return best_y;
}

