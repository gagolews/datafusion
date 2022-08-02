fit_wam_L2_regular <- function(X, Y, lambda) {
   stopifnot(is.matrix(X), is.matrix(Y))
   n <- nrow(X); m <- ncol(X)
   stopifnot(1 == nrow(Y), m == ncol(Y))
   stopifnot(is.numeric(lambda), length(lambda)==1)

   # Linear constraint (sum(w) == 1):
   A <- matrix(1, ncol=n, nrow=1)
   B <- 1

   # Objective function definition:
   D <- tcrossprod(X)+lambda*diag(n)
   C <- -tcrossprod(X, Y)

   res <- cgal_qp_solver(D, C, A, B, r="==", l=rep(0, n))
   testthat::expect_true(res$status == 0)
   w <- res$par
   testthat::expect_true(all(w >= 0))
   testthat::expect_equal(sum(w), 1)
   w
}
