fit_wam_L1_linprog_outrank <- function(X, Y, p) {
   stopifnot(is.matrix(X), is.matrix(Y))
   n <- nrow(X); m <- ncol(X)
   stopifnot(1 == nrow(Y), m == ncol(Y))
   stopifnot(is.numeric(p), length(p) == 1, p > 0)

   sigma <- order(Y)
   X[,] <- X[,sigma]
   Y[,] <- Y[,sigma]

   A <- rbind(
      cbind(t(X), -diag(m), diag(m), matrix(0, ncol=m-1, nrow=m)),
      c(rep(1, n), rep(0, 2*m), rep(0, m-1)),
      cbind(t(X[,-1]-X[,-m]), matrix(0, nrow=m-1, ncol=2*m), diag(m-1))
   )

   B <- c(Y, 1, rep(0, m-1))
   C <- c(rep(0, n), rep(1, 2*m), rep(p, m-1))
   D <- matrix(0, nrow=n+3*m-1, ncol=n+3*m-1) # an LP problem

   res <- cgal_qp_solver(D, C, A, B,
      r=c(rep("==", m+1), rep(">=", m-1)), l=rep(0, n+3*m-1))
   testthat::expect_true(res$status == 0)
   testthat::expect_equal(max(
      apply(matrix(res$par[n+1:(2*m)], nrow=2, byrow=TRUE), 2, min)
   ), 0)
   w <- res$par[1:n]
   testthat::expect_true(all(w >= 0))
   testthat::expect_equal(sum(w), 1)
   w
}


fit_wam_L2_quadprog_outrank <- function(X, Y, p) {
   stopifnot(is.matrix(X), is.matrix(Y))
   n <- nrow(X); m <- ncol(X)
   stopifnot(1 == nrow(Y), m == ncol(Y))
   stopifnot(is.numeric(p), length(p) == 1, p > 0)

   sigma <- order(Y)
   X[,] <- X[,sigma]
   Y[,] <- Y[,sigma]

   B <- c(1, rep(0, m-1)) # c(1, 0, 0, ..., 0)
   A <- rbind(
      c(rep(1, n), rep(0, m-1)),
      cbind(t(X[,-1]-X[,-m]), diag(m-1))
   )

   # Constraints are of the form A^T %*% w {=,>=} B
   # (first meq are equality constraints, the rest are >=)

   # Objective function definition:
#    D <- cbind(rbind(tcrossprod(X), matrix(0, nrow=m-1, ncol=n)), matrix(0, ncol=m-1, nrow=n+m-1))
#    C <- c(-tcrossprod(X, Y), rep(p, m-1))


   D <- cbind(
      rbind(tcrossprod(X), matrix(0, nrow=m-1, ncol=n)),
      rbind(matrix(0, nrow=n, ncol=m-1), p*diag(m-1))
   )
   C <- c(-tcrossprod(X, Y), rep(0, m-1))

   res <- cgal_qp_solver(D, C, A, B,
      r=c(rep("==", 1), rep(">=", m-1)), l=rep(0, n+m-1))

#    print(c(res$message, p))
   testthat::expect_true(res$status == 0)
   w <- res$par[1:n]
   testthat::expect_true(all(w >= 0))
   testthat::expect_equal(sum(w), 1)
   w
}



fit_wam_L1_linprog_outrank_old <- function(X, Y, p) {
   stopifnot(is.matrix(X), is.matrix(Y))
   n <- nrow(X); m <- ncol(X)
   stopifnot(1 == nrow(Y), m == ncol(Y))
   stopifnot(is.numeric(p), length(p) == 1, p >= 0)

   sigma <- order(Y)
   X[,] <- X[,sigma]
   Y[,] <- Y[,sigma]

   A <- rbind(
      cbind(t(X), -diag(m), diag(m), matrix(0, ncol=m-1, nrow=m)),
      c(rep(1, n), rep(0, 2*m), rep(0, m-1)),
      cbind(t(X[,-1]-X[,-m]), matrix(0, nrow=m-1, ncol=2*m), diag(m-1))
   )
   B <- c(Y, 1, rep(0, m-1))
   C <- c(rep(0, n), rep(1, 2*m), rep(p, m-1))

   T <- c(rep("=", m+1), rep(">=", m-1))

   # n weights and 3*m-1 auxiliary variables
   # lpSolveAPI implicitly assumes that all variables are >= 0
   lp <- lpSolveAPI::make.lp(ncol=n+2*m+m-1)
   lpSolveAPI::set.objfn(lp, C)
   for (j in 1:nrow(A))
      lpSolveAPI::add.constraint(lp, A[j,], T[j], B[j])

   testthat::expect_equal(0, lpSolveAPI::solve.lpExtPtr(lp))
   testthat::expect_equal(max(
      apply(matrix(lpSolveAPI::get.variables(lp)[n+1:(2*m)], nrow=2, byrow=TRUE), 2, min)
   ), 0)

   w <- lpSolveAPI::get.variables(lp)[1:n]
   testthat::expect_true(all(w >= 0))
   testthat::expect_equal(sum(w), 1)
   w
}
