#include <CGAL/basic.h>
#include <CGAL/QP_models.h>
#include <CGAL/QP_functions.h>
#include <CGAL/MP_Float.h>
#include <Rcpp.h>
using namespace Rcpp;
// program and solution types
typedef double ET;
typedef CGAL::Quadratic_program_solution<ET> Solution;
typedef CGAL::Quadratic_program_from_iterators<double**, double*,
   CGAL::Const_oneset_iterator<CGAL::Comparison_result>,
   bool*, double*, bool*, double*, double**, double*> Program;


//#include <iostream>

// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::export]]
NumericVector seb_solver(NumericMatrix eD, NumericVector ec) {
   // now construct the quadratic program:
   int n = ec.size();
   double *Ax = new double[n], **A = new double*[n];
   for (int i=0; i<n; ++i) { Ax[i] = 1.0; A[i] = Ax+i; }
   double b[] = {1.0}; // right-hand side
   bool *fl = new bool[n], *fu = new bool[n];
   double *l = new double[n], *u = new double[n];
   for (int i=0; i<n; ++i) { fl[i] = true; fu[i] = false;
                              l[i] = 0.0;   u[i] = 1.0; }
   double* d = REAL((SEXP)eD);  // underlying pointer
   double** D = new double*[n];
   for (int i=0; i<n; ++i) D[i] = d+n*i; // a symmetric matrix
   double* c = REAL((SEXP)ec);  // underlying pointer
   double c0 = 0.0; // constant term
   CGAL::Const_oneset_iterator<CGAL::Comparison_result>
      r(CGAL::EQUAL); // constraints are "="
   Program qp (n, 1, A, b, r, fl, l, fu, u, D, c, c0);

   // solve the program:
   Solution s = CGAL::solve_quadratic_program(qp, ET());
//   std::cout << s << std::endl;

   // generate output solution:
   NumericVector out(n);
   int i=0;
   for (auto it = s.variable_numerators_begin();
         it != s.variable_numerators_end(); ++it)
      out[i++] = (*it)/s.variables_common_denominator();

   delete Ax;   delete A;   delete fl;   delete fu;
   delete l;    delete u;   delete D;
   return out;
}
