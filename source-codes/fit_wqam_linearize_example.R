library(testthat)
options(stringsAsFactors=FALSE)
source("fit_wam.R")
source("fit_wam_outrank.R")
Sys.setenv(PKG_LIBS="-lCGAL")
Rcpp::sourceCpp("cgal_qp_solver.cpp")
source("fit_wqam.R")


###########################################################################

set.seed(19)
n <- 5
m <- 9
X <- t(round(matrix(runif(n*m, 0, 1), nrow=m), 2))

realw <-  runif(n)#c(0.23, 0.77) #
realw <- realw/sum(realw)
realw <- round(realw, 2)
stopifnot(sum(realw) == 1)
Y <- matrix(apply(X, 2, function(x) wqam(x, realw)), ncol=m)+rnorm(m, 0.05, 0.1)
Y <- round(Y, 2)                            # PERTURBATION
Y <- pmin(pmax(Y, apply(X, 2, min)), apply(X, 2, max))

stopifnot(dim(Y) == c(1, m))
stopifnot(dim(X) == c(n, m))
expect_equal(wqamX(realw, X, phi, phiInv), apply(X, 2, function(x) phiInv(sum(phi(x)*realw))))

###########################################################################

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

eps <- 1e-12
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
