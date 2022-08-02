Sys.setenv(PKG_LIBS="-lCGAL")
Rcpp::sourceCpp("cgal_qp_solver.cpp")

# QP

D <- matrix(c(2,0,0,8), nrow=2)
c <- c(0, -32)
c0 <- 64

l <- c(0, 0)
u <- c(Inf, 4)

A <- matrix(c(1, -1, 1, 2), ncol=2)
b <- c(7, 4)
r <- c("<=", "<=")

print(cgal_qp_solver(D, c, A, b, r, l, u, c0))


# LP

D <- matrix(0, nrow=2, ncol=2)
c <- c(0, -32)
c0 <- 64

l <- c(0, 0)
u <- c(Inf, 4)

b <- c(7, 4)
A <- matrix(c(1, -1, 1, 2), nrow=2)
r <- rep("<=", 2)

print(cgal_qp_solver(D, c, A, b, r, l, u, c0))
