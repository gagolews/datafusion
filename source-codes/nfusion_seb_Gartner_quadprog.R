# SMALLEST ENCLOSING BALL using quadprog

Sys.setenv(PKG_LIBS="-lCGAL")
#Rcpp::sourceCpp("cgal_qp_solver.cpp")

seb <- function(X) {
   stopifnot(is.numeric(X), is.matrix(X))
   n <- ncol(X)

   # the CGAL solver determines argmin_v 0.5 v^T D v + c^T v
   XtX <- crossprod(X) # (t(X) %*% X)
   D <- 2.0*XtX
   C <- -diag(XtX)
   A <- matrix(rep(1, n), ncol=n)
   B <- 1

   res <- cgal_qp_solver(D, C, A, B, r="==", l=rep(0, n))
   testthat::expect_true(res$status == 0)
   v <- res$par
   testthat::expect_true(all(v >= 0))
   testthat::expect_equal(sum(v), 1)
   as.numeric(tcrossprod(v, X)) # v %*% t(X)
}


Sys.setenv("PKG_LIBS"="-lCGAL")
#Rcpp::sourceCpp("nfusion_seb_Gartner_quadprog.cpp")

seb_old <- function(X) {
   stopifnot(is.numeric(X), is.matrix(X))

#    ch <- chull(X[1,],X[2,])
#    X <- X[,ch]

   # the CGAL solver determines argmin_v 0.5 v^T D v + c^T v
   D <- 2.0*crossprod(X) # == t(X) %*% X
   c <- -apply(X, 2L, # apply a function on each X's column
      function(xi) crossprod(xi)) # t(xi) %*% xi

   w <- seb_solver(D, c)
   return(colSums(t(X)*w))
}
