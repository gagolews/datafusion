library(testthat)
options(stringsAsFactors=FALSE)
Sys.setenv(PKG_LIBS="-lCGAL")
Rcpp::sourceCpp("cgal_qp_solver.cpp")
source("fit_wqam.R")

powmeanX <- compiler::cmpfun(function(w, X, p) {
   expect_equal(sum(w), 1)
   stopifnot(is.numeric(p), length(p) == 1, p > 0)
#    as.numeric(phiInv(phiXt %*% w))
   as.numeric((t(X^p) %*% w)^(1/p))
})

error_powmeanX <- compiler::cmpfun(function(w, X, Y, p) {
   Yobs <- powmeanX(w, X, p)
   c(L1=sum(abs(  Yobs-Y   )   ),
     L2=sum(   (  Yobs-Y   )^2 )^0.5
     # LInf=max(abs(  Yobs-Y   )   )
   )
})


############# INPUT DATA #########################################



## CHECK:
stopifnot(dim(Y) == c(1, m))
stopifnot(dim(X) == c(n, m))


p <- 2
phi         <- function(x) x^p
phiInv      <- function(x) exp(log(x)/p) # x^(1/p)
phiInvPrime <- function(x) exp((1-p)*log(x)/p)/p # (1/p)*x^(1/p-1)

envir_p <- new.env()
envir_p[["p"]] <- p
environment(phi) <- envir_p
environment(phiInv) <- envir_p
environment(phiInvPrime) <- envir_p


set.seed(132)
n <- 2
m <- 9
X <- t(matrix(runif(n*m, 0, 1), nrow=m))

realw <- runif(n)
realw <- realw/sum(realw)
Y <- matrix(as.numeric((t(X^p) %*% realw)^(1/p))+rt(m, 5)*0.05, ncol=m)

print(environment(phi)$p)

Ps <- seq(0.5, 6, length=101)
res1 <- sapply(Ps, function(p) {
   assign("p", p, environment(phi))
   replicate(10, {tryCatch({
      w <- fit_wqam_L2_optim(X, Y, phi, phiInv, phiInvPrime)
      error_powmeanX(w, X, Y, p)}, error=function(e) c(NA, NA)
   )}) -> tmp
   tmp[,which.min(tmp[2,])]
})

res2 <- sapply(Ps, function(p) {
   assign("p", p, environment(phi))
   replicate(10, {tryCatch({
      w <- fit_wqam_L1_optim_approx(X, Y, phi, phiInv, phiInvPrime, 1e-12)
      error_powmeanX(w, X, Y, p)}, error=function(e) c(NA, NA)
   )}) -> tmp
   tmp[,which.min(tmp[1,])]
})


fit_powmean_L2_optim <- function(X, Y, pmin=0.1, pmax=10) {
   stopifnot(is.matrix(X), is.matrix(Y))
   n <- nrow(X); m <- ncol(X)
   stopifnot(1 == nrow(Y), m == ncol(Y))

   # p <- 1 # this will be a parameter shared by the 3 functions:
   phi         <- function(x) x^p
   phiInv      <- function(x) exp(log(x)/p) # x^(1/p)
   phiInvPrime <- function(x) exp((1-p)*log(x)/p)/p # (1/p)*x^(1/p-1)
   envir_p <- new.env()
   envir_p[["p"]] <- p
   environment(phi) <- envir_p
   environment(phiInv) <- envir_p
   environment(phiInvPrime) <- envir_p

   E <- function(p) {
      assign("p", p, environment(phi)) # affects 3 functions
      w <- fit_wqam_L2_optim(X, Y, phi, phiInv, phiInvPrime)
      sum((as.numeric((t(X^p) %*% w)^(1/p))-Y)^2)
   }

   optimize(E, c(pmin, pmax))$minimum
}



# par(mfrow=c(2,1))
# matplot(Ps, t(res1), type='l', xlab=expression(p), ylab=expression(Error(p)), lty=c(1,2,4), col=c(1,2,4), log="y", ylim=range(c(res1,res2)))
# legend("top", horiz=TRUE, expression(L[1], L[2], L[infinity]), lty=c(1,2,4), col=c(1,2,4))

# matplot(Ps, t(res2), type='l', xlab=expression(p), ylab=expression(Error(p)), lty=c(1,2,4), col=c(1,2,4), log="y", ylim=range(c(res1,res2)))

# pdf(sprintf("~/Publications/Books/Habilitacja/figures/fit_powmean_error_p.pdf"), height=4, width=6)
par(mar=c(4, 4.2, 0.5, 0.5))
# matplot(Ps, cbind(t(res1),t(res2)), type='l', xlab=expression(p), ylab=expression(Error(p)), lty=c(1,2,1,2), col=c(1,1,2,2), log="y", lwd=c(1,2,2,1), las=1)
# legend("bottomright", horiz=FALSE, expression(d[1]~(LSE), d[2]~(LSE), d[1]~(LAD), d[2]~(LAD)), lty=c(1,2,1,2), col=c(1,1,2,2), lwd=c(1,2,2,1))
matplot(Ps, cbind(t(res1)), type='l', xlab=expression(p), ylab=expression(Error(p)), lty=c(1,2,1,2), col=c(1,1,2,2), log="y", lwd=c(1,2,2,1), las=1)
legend("bottomright", horiz=FALSE, expression(d[1]~(LSE), d[2]~(LSE)), lty=c(1,2,1,2), col=c(1,1,2,2), lwd=c(1,2,2,1))
abline(v=fit_powmean_L2_optim(X, Y), col='gray', lty=3)
# dev.off()

##################################################################

