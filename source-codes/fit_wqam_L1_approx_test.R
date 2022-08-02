library(parallel)
options("mc.cores"=3L)
options("cl.cores"=3L)

n <- 5
m <- 25
k <- 1
eps1 <- 1e-9
eps2 <- 1e-12

print(environment(phi)$p)
library(testthat)
options(stringsAsFactors=FALSE)
source("utils.R")
source("fit_wam.R")
source("fit_wam_outrank.R")
Sys.setenv(PKG_LIBS="-lCGAL")
Rcpp::sourceCpp("cgal_qp_solver.cpp")
source("fit_wqam.R")




fit_wqam_L1_optim_approx_multi <- compiler::cmpfun(function(X, Y, phi, phiInv, phiInvPrime, eps, k) {
   replicate(k,
      tryCatch({
         w <- fit_wqam_L1_optim_approx(X, Y, phi, phiInv, phiInvPrime, eps)
         c(as.numeric(error_wqam(w, X, Y, phi, phiInv)["L1"]), as.numeric(w))
      }, error=function(e) rep(NA_real_, 1+nrow(X))
   )) -> res

   if (k > 1 && !all(is.na(res[-1,]))) {
      i <- which.min(res[1,])
      res[-1, i]
   }
   else
      as.numeric(res)[2:(nrow(X)+1)]
})


fit_wqam_L1_optim_multi <- compiler::cmpfun(function(X, Y, phi, phiInv, k) {
   replicate(k,
      tryCatch({
         w <- fit_wqam_L1_optim(X, Y, phi, phiInv)
         c(as.numeric(error_wqam(w, X, Y, phi, phiInv)["L1"]), as.numeric(w))
      }, error=function(e) rep(NA_real_, 1+nrow(X))
   )) -> res

   if (k > 1 && !all(is.na(res[-1,]))) {
      i <- which.min(res[1,])
      res[-1, i]
   }
   else
      as.numeric(res)[2:(nrow(X)+1)]
})


print(system.time(replicate2(10000, {
   X <- t(matrix(runif(n*m, 0, 1), nrow=m))
   realw <- runif(n)
   realw <- realw/sum(realw)
   Y <- matrix(apply(X, 2, function(x) wqam(x, realw)), ncol=m)+rnorm(m, 0, 0.05)
   results1 <- tryCatch({w <- fit_wam_L1_linprog(phi(X), phi(Y)); error_wqam(w, X, Y, phi, phiInv)["L1"]}, error=function(e) NA_real_)
   results2 <- tryCatch({w <- fit_wqam_L1_optim_multi(X, Y, phi, phiInv, k); error_wqam(w, X, Y, phi, phiInv)["L1"]}, error=function(e) NA_real_)
   results3 <- tryCatch({w <- fit_wqam_L1_optim_approx_multi(X, Y, phi, phiInv, phiInvPrime, eps1, k); error_wqam(w, X, Y, phi, phiInv)["L1"]}, error=function(e) NA_real_)
   results4 <- tryCatch({w <- fit_wqam_L1_optim_approx_multi(X, Y, phi, phiInv, phiInvPrime, eps2, k); error_wqam(w, X, Y, phi, phiInv)["L1"]}, error=function(e) NA_real_)

   c(lp=results1, abs=results2, erf1=results3, erf2=results4)
}, structure(numeric(4), names=c("LP", "ABS", "APPROX1", "APPROX2"))) -> res))



old_scipen <- options("scipen")$scipen
options(scipen=16)
print(apply(res, 1, quantile, c(0,0.25,0.5,0.75,0.99,0.999,1), na.rm=TRUE), digits=22)
relerr <- (t(res[-1,])-as.numeric(res[1,]))/as.numeric(res[1,])
print(apply(relerr, 2, quantile, c(0,0.25,0.5,0.75,0.99,0.999,1), na.rm=TRUE), digits=22)
print(apply(relerr, 2, function(x) c(NumNA=sum(is.na(x)))))
# add <- min(relerr+, na.rm=TRUE)+1e-16
options(scipen=old_scipen)

boxplot(pmax(relerr, 1e-16), log="y", names=c("abs", sprintf("approx_eps=%g", c(eps1, eps2))))

