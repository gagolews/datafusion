/* *****************************************************************************

CGAL quadratic programming solver interface for R

Compile in R using:
   Sys.setenv(PKG_LIBS="-lCGAL")
   Rcpp::sourceCpp("cgal_qp_solver.cpp")
Note that CGAL-devel libraries must be available in your system.


Copyright (c) 2015 Marek Gagolewski

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

***************************************************************************** */


/*** R
#' @title       CGAL Quadratic Programming Solver
#' @author      Marek Gagolewski
#' @description An R interface to the CGAL library QP solver.
#' Minimizes \eqn{ 0.5 x^T D x + c^T x + c0 }
#' subject to
#'    \eqn{ Ax >=< b },
#'    \eqn{ l <= x <= u }
#' with respect to \eqn{ x = (x_1, ..., x_n) }.
#'
#' @details
#' Note that the CGAL 4.6 User Guide claims that the problem being
#' solved is \eqn{ x^T D x + c^T x + c0 } (without the 0.5 term),
#' but this is not the case.
#'
#' @param D a symmetric positive-semidefinite \eqn{n*n} numeric matrix
#' @param c a numeric vector of length \eqn{n}
#' @param A an \eqn{m*n} numeric matrix
#' @param b a numeric vector of length \eqn{m}
#' @param r a character vector of length \eqn{m}
#'        with elements like \code{<=}, \code{==}, or \code{>=};
#'        specifies types of linear constraints
#' @param l a numeric vector of length \eqn{n}
#'        which gives lower bounds for corresponding \eqn{x} variables,
#'        \code{-Inf} gives no bound
#' @param u a numeric vector of length \eqn{n}
#'        which gives lower bounds for corresponding \eqn{x} variables,
#'        \code{+Inf} gives no bound
#' @param c0 a single numeric value
#'
#' @return
#' A list with the following components:
#' \describe{
#'    \item{\code{par}}{The best set of parameters, \eqn{x}, found}
#'    \item{\code{value}}{The value of the objective function at \code{par}}
#'    \item{\code{counts}}{The number of iterations that it took to solve the program}
#'    \item{\code{status}}{An integer code: 0 indicates successful completion}
#'    \item{\code{message}}{Status message: optimal, infeasible, or unbounded}
#' }
#'
#' @references
#' CGAL 4.6 User Manual, \emph{Linear and Quadratic Programming Solver},
#' \url{http://doc.cgal.org/latest/QP_solver/index.html}
cgal_qp_solver <- function(D, c, A, b, r=rep(">=", length(b)),
      l=rep(-Inf, length(c)), u=rep(Inf, length(c)), c0=0.0)
{
   stopifnot(is.numeric(D), is.finite(D), is.matrix(D))
   stopifnot(is.numeric(A), is.finite(A), is.matrix(A))
   stopifnot(is.numeric(c), is.finite(c))
   stopifnot(is.numeric(b), is.finite(b))
   stopifnot(is.character(r), r %in% c("<=", "==", ">="))
   stopifnot(is.numeric(l), !is.na(l) & !is.nan(l))
   stopifnot(is.numeric(u), !is.na(u) & !is.nan(u))
   stopifnot(is.numeric(c0), is.finite(c0))

   stopifnot(length(b) == nrow(A))
   stopifnot(length(r) == nrow(A))
   stopifnot(isSymmetric(D))
   stopifnot(ncol(D) == ncol(A))
   stopifnot(length(c) == nrow(D))
   stopifnot(length(c0) == 1)
   stopifnot(length(l) == nrow(D))
   stopifnot(length(u) == nrow(D))
   stopifnot(length(c) > 0, length(b) > 0)

   r <- match(r, c("<=", "==", ">="))-2 # gives values in {-1, 0, 1}
   fl <- is.finite(l) # which lower bounds for x are active
   fu <- is.finite(u) # which upper bounds for x are active
   .cgal_qp_solver(length(c), length(b), A, b, r, fl, l, fu, u, D, c, c0)
}
*/




#define ASSERT(x) { \
   if ((!(x))) stop(#x); \
}

#include <CGAL/QP_functions.h>
#include <CGAL/MP_Float.h>
#include <Rcpp.h>
using namespace Rcpp;
// program and solution types
typedef CGAL::MP_Float ET;
typedef CGAL::Quadratic_program_solution<ET> Solution;
typedef CGAL::Quadratic_program_from_iterators<double**, double*,
   CGAL::Comparison_result*,
   int*, double*, int*, double*, double**, double*> Program;

// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::export(".cgal_qp_solver")]]
List cgal_qp_solver(
      int n, int m,
      NumericMatrix A, NumericVector b,
      IntegerVector r,
      LogicalVector fl, NumericVector l,
      LogicalVector fu, NumericVector u,
      NumericMatrix D, NumericVector c, double c0)
{
   double* Aptr = REAL(static_cast<SEXP>(A));
   double* Dptr = REAL(static_cast<SEXP>(D));
   double** _A = new double*[n]; // go columnwise over the constraint matrix A
   double** _D = new double*[n];
   /* _D[i] is a random access iterator for the entries in row i
       below or on the diagonal. The valid range of this iterator is guaranteed
       to have length i+1 but not more;
       note that D is symmetric */
   for (int i=0; i<n; ++i) {
      _A[i] = Aptr+i*m; // ith column
      _D[i] = Dptr+i*n; // ith column
   }

   CGAL::Comparison_result* _r = new CGAL::Comparison_result[m];
   for (int j=0; j<m; ++j) {
      _r[j] = (r[j] < 0 ? CGAL::SMALLER
                        : (r[j] > 0 ? CGAL::LARGER : CGAL::EQUAL));
   }

   // solve the program:
   List retval;
   try {
      Program qp(n, m, _A, REAL(static_cast<SEXP>(b)), _r,
         (int*)LOGICAL(static_cast<SEXP>(fl)), REAL(static_cast<SEXP>(l)),
         (int*)LOGICAL(static_cast<SEXP>(fu)), REAL(static_cast<SEXP>(u)),
         _D, REAL(static_cast<SEXP>(c)), c0);
      Solution s(CGAL::solve_quadratic_program(qp, ET()));
      // generate output solution:
      NumericVector solution(n);
      int i=0;
      for (auto it = s.variable_values_begin();
            it != s.variable_values_end(); ++it)
         solution[i++] = to_double(*it);
      retval = List::create(
         _("par") = solution,
         _("value") = to_double(s.objective_value()),
         _("counts") = s.number_of_iterations(),
         _("status") = (s.status() == CGAL::QP_OPTIMAL ? 0
                        : (s.status() == CGAL::QP_INFEASIBLE ? 1
                        : (s.status() == CGAL::QP_UNBOUNDED ? 2
                        : -1))),
         _("message") = (s.status() == CGAL::QP_OPTIMAL ? "optimal"
                        : (s.status() == CGAL::QP_INFEASIBLE ? "infeasible"
                        : (s.status() == CGAL::QP_UNBOUNDED ? "unbounded"
                        : NULL)))
      );
   } catch(...) {
      stop("Fatal error in .cgal_qp_solver()");
   };

   delete [] _r; delete [] _A; delete [] _D;
   return retval;
}
