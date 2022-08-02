library(testthat)
options(stringsAsFactors=FALSE)
source("fit_wam.R")
source("fit_wam_outrank.R")
Sys.setenv(PKG_LIBS="-lCGAL")
Rcpp::sourceCpp("cgal_qp_solver.cpp")
source("fit_wqam.R")

############# INPUT DATA #########################################

set.seed(123)
n <- 2
m <- 10
X <- t(matrix(runif(n*m, 0, 1), nrow=m))

realw <- runif(n)
realw <- realw/sum(realw)
Y <- matrix(apply(X, 2, function(x) wqam(x, realw)), ncol=m)+rnorm(m, 0, 0.05)

## CHECK:
stopifnot(dim(Y) == c(1, m))
stopifnot(dim(X) == c(n, m))
expect_equal(wqamX(realw, X, phi, phiInv), apply(X, 2, function(x) phiInv(sum(phi(x)*realw))))



##################################################################

results <- matrix(NA_real_, nrow=0, ncol=n+3, dimnames=list(NULL, c(paste0("w", 1:n), "L1", "L2", "LInf")))

############# 1. USE realw & compute WAM fits      ######################

cat("\r1", file=stderr())
w <- realw
results <- rbind(results, exact=c(w, error_wqam(w, X, Y, phi, phiInv)))

cat("\r2", file=stderr())
w <- fit_wam_L1_linprog(phi(X), phi(Y))
results <- rbind(results, WAM_L1_linearize=c(w, error_wqam(w, X, Y, phi, phiInv)))

cat("\r3", file=stderr())
w <- fit_wam_L2_quadprog(phi(X), phi(Y))
results <- rbind(results, WAM_L2_linearize=c(w, error_wqam(w, X, Y, phi, phiInv)))

cat("\r4", file=stderr())
w <- fit_wam_LInf_linprog(phi(X), phi(Y))
results <- rbind(results, WAM_LInf_linearize=c(w, error_wqam(w, X, Y, phi, phiInv)))

#####################################################

cat("\r5", file=stderr())
w <- fit_wqam_L2_optim(X, Y, phi, phiInv, phiInvPrime)
results <- rbind(results, WQAM_L2_optim=c(w, error_wqam(w, X, Y, phi, phiInv)))


cat("\r6", file=stderr())
w <- fit_wqam_L1_optim(X, Y, phi, phiInv)
results <- rbind(results, WQAM_L1_optim=c(w, error_wqam(w, X, Y, phi, phiInv)))

eps <- 1e-8/m
cat("\r7", file=stderr())
w <- fit_wqam_L1_optim_approx(X, Y, phi, phiInv, phiInvPrime, eps)
results <- rbind(results, WQAM_L1_optim_approx=c(w, error_wqam(w, X, Y, phi, phiInv)))

############### OUTPUT #############################

cat("\n", file=stderr())
results <- cbind(results, L1rel=(results[,"L1"]-min(results[,"L1"]))/min(results[,"L1"]))
results <- cbind(results, L2rel=(results[,"L2"]-min(results[,"L2"]))/min(results[,"L2"]))
results <- cbind(results, LInfrel=(results[,"LInf"]-min(results[,"LInf"]))/min(results[,"LInf"]))
old_scipen <- options("scipen")$scipen
options(scipen=16)
print(results, digits=4, scientific = FALSE)
options(scipen=old_scipen)



# ############## 2. LAD - linear programming ####################################
#
#
# A <- rbind(
#    cbind(phi(t(X)), -diag(m), diag(m)),
#    c(rep(1, n), rep(0, 2*m))
# )
# B <- rbind(phi(t(Y)), 1)
# C <- matrix(c(rep(0, n), rep(1, 2*m)), ncol=1)
#
# library("lpSolveAPI")
# lp <- make.lp(ncol=n+2*m)
# set.objfn(lp, C)
# for (j in 1:nrow(A))
#    add.constraint(lp, A[j,], "=", B[j])
#
#
#
# print(solve(lp)) # 0 == optimal solution found
# w <- get.variables(lp)[1:n]
# expect_equal(max(apply(matrix(get.variables(lp)[-(1:n)], nrow=2, byrow=TRUE), 2, min)), 0)
#
# results[[length(results)+1]] <- list("LAD_linprog", w, error_wqam(w))
# str(list(GRAD=error1wp_erf_grad(w[-n])))
# print(error1(w)-error1_erf(w))
#
# ############# 3. LAD -- optim #######################
#
# # res <- optim(runif(n-1), function(w) {
# #    stopifnot(all(w >= 0), sum(w) <= 1)
# #    w = c(w, 1-sum(w))
# #    error(X, Y, w)["L1"]
# # }, method='L-BFGS-B', lower=0, upper=1)
# #
# # print(res$message)
# # w <- c(res$par, 1-sum(res$par))
# # results[[length(results)+1]] <- list("LAS_optim", w, error(X, Y, w))
#
#
# A <- rbind(rep(-1, n-1), diag(n-1))
# B <- c(-1, rep(0, n-1))
# startw <- runif(n-1)
# startw <- startw/sum(startw)*0.9
# # res <- constrOptim(startw,
# #    function(w) {
# #       stopifnot(all(w >= 0), sum(w) <= 1+1e-12)
# #       w <- c(w, 1-sum(w))
# #       Yobs <- wqamX(w) #apply(X, 2, function(x) wqam(x, w))
# # #       sum(abs(Yobs-y)) # L1 - non-smooth
# # #       sum(sqrt( (Yobs-Y)^2 + eps^2 )) # upper bnd
# # #       sum((Yobs-Y)^2/sqrt( (Yobs-Y)^2 + eps^2 )) # lower bnd
# #       sum((Yobs-Y)*erf((Yobs-Y)/eps))
# # }, grad=NULL, A, B, control=list(maxit=1000))
# res <- constrOptim(startw, f=error1wp_erf, grad=error1wp_erf_grad, A, B, control=list(maxit=10000, reltol=1e-16), outer.eps=1e-9)
# str(list(STATUS=res$convergence, MSG=res$message, GRAD=error1wp_erf_grad(res$par), GRAD2=numDeriv::grad(error1wp, res$par)))
# w <- c(res$par, 1-sum(res$par))
# # print(error1(w)-error1_erf(w))
# results[[length(results)+1]] <- list("LAS_optim_erf", w, error(w))
#
# res <- constrOptim(res$par, f=error1wp, grad=NULL, A, B, control=list(maxit=10000, reltol=1e-16), outer.eps=1e-9)
# w <- c(res$par, 1-sum(res$par))
# # print(error1(w)-error1_erf(w))
# results[[length(results)+1]] <- list("LAS_optim_erf+refine", w, error(w))
#
# res <- constrOptim(startw, f=error1wp, grad=NULL, A, B, control=list(maxit=10000, reltol=1e-16), outer.eps=1e-9)
# w <- c(res$par, 1-sum(res$par))
# results[[length(results)+1]] <- list("LAS_optim_abs", w, error(w))
#
#
#
# # wp <- w[-n]
# # repeat {
# #    wpnew <- wp-1e-12*error1wp_erf_grad(wp)
# #    if (error1wp(wp) < error1wp(wpnew)) break
# #    wp <- wpnew
# #    print(error1wp(wp))
# # }
#
# ############## 4. MSE - quad prog #################################
#
# # Constraint definition:
# meq <- 1 # number of equality constraints
# B <- c(1, rep(0, n)) # c(1, 0, 0, ..., 0)
# A <- cbind(rep(1, n), diag(n)) # 1st column: ones,
#                                # then the n*n diagonal matrix
# # Constraints are of the form A^T %*% w {=,>=} B
# # (first meq are equality constraints, the rest is of >= type)
#
# # Objective function definition:
# D <- tcrossprod(phi(X)) # X %*% t(X)
# C <- phi(X) %*% phi(t(Y))
#
# # Solve min(0.5 * w^T %*% D %*% w - C^T %*% w) for w
# w <- quadprog::solve.QP(D, C, A, B, meq)$solution
# results[[length(results)+1]] <- list("MSE_quadprog", w, error(w))
# str(list(GRAD=error2wp_grad(w[-n])))
#
# ############# 5. MSE -- optim #######################
#
# # ui %*% theta >= ci
# # wi >= 0 i=1..n-1
# # sum(wi) <= 1 => sum(-wi) >= -1
#
# A <- rbind(rep(-1, n-1), diag(n-1))
# B <- c(-1, rep(0, n-1))
# startw <- runif(n-1)
# startw <- startw/sum(startw)*0.9
# # res <- constrOptim(startw,
# #    function(w) {
# #       stopifnot(all(w >= 0), sum(w) <= 1+1e-12)
# #       w <- c(w, 1-sum(w))
# #       error(w)["L2"]
# # }, grad=NULL, A, B, control=list(maxit=1000))
#
# res <- constrOptim(startw, error2wp, error2wp_grad, A, B, control=list(maxit=1000, reltol=1e-15))
# str(list(STATUS=res$convergence, MSG=res$message, GRAD=error2wp_grad(res$par)))
#
# w <- c(res$par, 1-sum(res$par))
# results[[length(results)+1]] <- list("MSE_optim", w, error(w))
#
# # res <- optim(runif(n-1), function(w) {
# #    stopifnot(all(w >= 0), sum(w) <= 1)
# #    w = c(w, 1-sum(w))
# #    error(X, Y, w)["L2"]
# # }, method='L-BFGS-B', lower=0, upper=1)
#
#
# ############# 6. Chebyshev -- linprog #######################
#
# w <- fit_wam_LInf_linprog(X, Y)
# results[[length(results)+1]] <- list("Cheb_linprog", w, error(w))
#
#
# ############### OUTPUT #############################
#
# results <- matrix(unlist(results), nrow=length(unlist(results))/(1+3+n), byrow=TRUE)
# results <- as.data.frame(results)
# names(results) <- c("method", paste0("w", 1:n), "L1", "L2", "LInf")
# for (i in 2:ncol(results)) results[[i]] <- as.numeric(results[[i]])
# results$L1rel <- results$L1-min(results$L1)
# results$L2rel <- results$L2-min(results$L2)
# results$LInfrel <- results$LInf-min(results$LInf)
# print(results)
#
# # xtable::xtable(rbind(X, Y), digits=2)
#
# # cat("-----------------------------\n")
# # w <- lm.fit(phi(t(X)), phi(t(Y)), tol=1e-12)$coefficients
# # print(w <- w/sum(w))
# # print(error(X, Y, w))
#
# # library("scatterplot3d")
# # scatterplot3d(cbind(X, y))
# # library("rgl")
# # plot3d(X[,1], X[,2], y)
# # lines3d(c(0,1),c(0,1),c(0,1),col=2)
# # polygon3d(c(0,0,1,1), c(0,1,1,0), c(0,0,1,1), col="#f0f0f020")
# # polygon3d(c(0,0,1,1), c(0,1,1,0), c(0,1,1,0), col="#f0f0f020")
# # polygon3d(c(0,0,1,1), c(0,1,1,0), c(0,0.8,1,0.2), col="#f030f020")
#
#
