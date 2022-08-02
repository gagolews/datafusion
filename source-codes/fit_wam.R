wam <- function(x, w) {
   stopifnot(length(w) == length(x))
   expect_equal(sum(w), 1)
   sum(w*x)
}


wamX <- function(w) {
   expect_equal(sum(w), 1)
   as.numeric(t(X) %*% w)
}


error_wam <- function(w) {
   Yobs <- wamX(w) #apply(X, 2, function(x) wqam(x, w))
   c(L1=sum(abs(  Yobs-Y   )   ),
     L2=sum(   (  Yobs-Y   )^2 )^0.5,
     LInf=max(abs(  Yobs-Y   )   )
   )
}





fit_wam_L1_linprog <- function(X, Y) {
   stopifnot(is.matrix(X), is.matrix(Y))
   n <- nrow(X); m <- ncol(X)
   stopifnot(1 == nrow(Y), m == ncol(Y))

   A <- rbind(
      cbind(t(X), -diag(m), diag(m)),
      c(rep(1, n), rep(0, 2*m))
   )
   B <- c(Y, 1)
   C <- c(rep(0, n), rep(1, 2*m))
   D <- matrix(0, nrow=n+2*m, ncol=n+2*m) # an LP problem

   res <- cgal_qp_solver(D, C, A, B, r=rep("==", nrow(A)),
                                     l=rep(0, n+2*m))
   testthat::expect_true(res$status == 0)
   testthat::expect_equal(max(
      apply(matrix(res$par[-(1:n)], nrow=2, byrow=TRUE), 2, min)
   ), 0)
   w <- res$par[1:n]
   testthat::expect_true(all(w >= 0))
   testthat::expect_equal(sum(w), 1)
   w
}




fit_wam_L2_quadprog <- function(X, Y) {
   stopifnot(is.matrix(X), is.matrix(Y))
   n <- nrow(X); m <- ncol(X)
   stopifnot(1 == nrow(Y), m == ncol(Y))

   # Linear constraint (sum(w) == 1):
   A <- matrix(1, ncol=n, nrow=1)
   B <- 1

   # Objective function definition:
   D <- tcrossprod(X)     #    X %*% t(X)
   C <- -tcrossprod(X, Y) # - (X %*% t(Y))

   res <- cgal_qp_solver(D, C, A, B, r="==", l=rep(0, n))
   testthat::expect_true(res$status == 0)
   w <- res$par
   testthat::expect_true(all(w >= 0))
   testthat::expect_equal(sum(w), 1)
   w
}


fit_wam_LInf_linprog <- function(X, Y) {
   stopifnot(is.matrix(X), is.matrix(Y))
   n <- nrow(X); m <- ncol(X)
   stopifnot(1 == nrow(Y), m == ncol(Y))

   A <- rbind(
      cbind(t(X), -1),
      cbind(t(X), +1),
      c(rep(1, n), 0)
   )
   B <- c(Y, Y, 1)
   C <- c(rep(0, n), 1)
   D <- matrix(0, nrow=n+1, ncol=n+1) # an LP problem

   res <- cgal_qp_solver(D, C, A, B,
      r=c(rep("<=", m), rep(">=", m), "=="),
      l=rep(0, n+1))
   testthat::expect_true(res$status == 0)
   w <- res$par[1:n]
   testthat::expect_true(all(w >= 0))
   testthat::expect_equal(sum(w), 1)
   w
}





fit_wam_L1_linprog_old <- function(X, Y) {
   stopifnot(is.matrix(X), is.matrix(Y))
   n <- nrow(X); m <- ncol(X)
   stopifnot(1 == nrow(Y), m == ncol(Y))

   A <- rbind(
      cbind(t(X), -diag(m), diag(m)),
      c(rep(1, n), rep(0, 2*m))
   )
   B <- c(Y, 1)
   C <- c(rep(0, n), rep(1, 2*m))

   # n weights and 2*m auxiliary variables
   # lpSolveAPI implicitly assumes that all variables are >= 0
   lp <- lpSolveAPI::make.lp(ncol=n+2*m)
   lpSolveAPI::set.objfn(lp, C)
   for (j in 1:nrow(A))
      lpSolveAPI::add.constraint(lp, A[j,], "=", B[j])

   testthat::expect_equal(0, lpSolveAPI::solve.lpExtPtr(lp))
   testthat::expect_equal(max(
      apply(matrix(lpSolveAPI::get.variables(lp)[-(1:n)], nrow=2, byrow=TRUE), 2, min)
   ), 0)

   w <- lpSolveAPI::get.variables(lp)[1:n]
   testthat::expect_true(all(w >= 0))
   testthat::expect_equal(sum(w), 1)
   w
}

fit_wam_L2_quadprog_old <- function(X, Y) {
   stopifnot(is.matrix(X), is.matrix(Y))
   n <- nrow(X); m <- ncol(X)
   stopifnot(1 == nrow(Y), m == ncol(Y))

   meq <- 1 # number of equality constraints
   B <- c(1, rep(0, n)) # c(1, 0, 0, ..., 0)
   A <- cbind(rep(1, n), diag(n)) # 1st column: ones,
                               # then the n*n diagonal matrix
   # Constraints are of the form A^T %*% w {=,>=} B
   # (first meq are equality constraints, the rest are >=)

   # Objective function definition:
   D <- tcrossprod(X)    # X %*% t(X)
   C <- tcrossprod(X, Y) # X %*% t(Y)

   # Solve min(0.5 * w^T %*% D %*% w - C^T %*% w) for w
   w <- quadprog::solve.QP(D, C, A, B, meq)$solution
   testthat::expect_true(all(w >= -1e-16))
   testthat::expect_equal(sum(w), 1)
   w
}



fit_wam_LInf_linprog_old <- function(X, Y) {
   stopifnot(is.matrix(X), is.matrix(Y))
   n <- nrow(X); m <- ncol(X)
   stopifnot(1 == nrow(Y), m == ncol(Y))

   A <- rbind(
      cbind(t(X), -1),
      cbind(t(X), +1),
      c(rep(1, n), 0)
   )
   B <- c(Y, Y, 1)
   C <- c(rep(0, n), 1)
   T <- c(rep("<=", m), rep(">=", m), "=")

   # n weights and 1 auxiliary variable;
   # lpSolveAPI implicitly assumes that all variables are >= 0
   lp <- lpSolveAPI::make.lp(ncol=n+1)
   lpSolveAPI::set.objfn(lp, C)
   for (j in 1:nrow(A))
      lpSolveAPI::add.constraint(lp, A[j,], T[j], B[j])

   testthat::expect_equal(0, lpSolveAPI::solve.lpExtPtr(lp))

   w <- lpSolveAPI::get.variables(lp)[1:n]
   testthat::expect_true(all(w >= 0))
   testthat::expect_equal(sum(w), 1)
   w
}
