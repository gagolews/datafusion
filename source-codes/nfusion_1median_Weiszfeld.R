Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
Rcpp::cppFunction('
NumericVector C_Weiszfeld1median(NumericMatrix X, NumericVector w,
                               NumericVector y0, double eps=1.0e-9) {
   int d = X.nrow();
   int n = X.ncol();
   if (w.length()  != n) stop("w.length()  != n");
   if (y0.length() != d) stop("y0.length() != d");
   NumericVector y_last(d); // a new vector
   NumericVector y_cur(Rcpp::clone(y0)); // a deep copy of y0
   double lasterr;
   do {
      std::swap(y_last, y_cur); // just swaps underlying pointers
      for (int j=0; j<d; ++j) y_cur[j] = 0.0;
      double w_over_d_x_y = 0.0;
      for (int i=0; i<n; ++i) {
         double d_xi_y = 0.0;
         for (int j=0; j<d; ++j)
            d_xi_y += (X(j, i)-y_last[j])*(X(j, i)-y_last[j]);
         d_xi_y = sqrt(d_xi_y);
         if (d_xi_y <= eps) return y_last;
         double w_over_d_xi_y = w[i]/d_xi_y;
         w_over_d_x_y += w_over_d_xi_y;
         for (int j=0; j<d; ++j)
            y_cur[j] += w_over_d_xi_y*X(j, i);
      }
      lasterr = 0.0;
      for (int j=0; j<d; ++j) {
         y_cur[j] /= w_over_d_x_y;
         lasterr += (y_cur[j]-y_last[j])*(y_cur[j]-y_last[j]);
      }
   } while (lasterr > eps*eps);
   return y_cur;
}
', includes=c("#include <algorithm>", "#include <utility>"))



Weiszfeld1median <- compiler::cmpfun(function(X, w=rep(1, ncol(X))/ncol(X), y0=rowMeans(X), eps=1.0e-9)
   C_Weiszfeld1median(X, w, y0, eps))
