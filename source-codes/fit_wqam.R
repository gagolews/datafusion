################ PHI GENERATING FUNS ###############################

## CASE 1
p <- 2
phi         <- function(x) x^p
phiInv      <- function(x) x^(1/p)
phiInvPrime <- function(x) (1/p)*(x^(1/p-1))

envir_p <- new.env()
envir_p[["p"]] <- p
environment(phi) <- envir_p
environment(phiInv) <- envir_p
environment(phiInvPrime) <- envir_p

# x <- seq(0.1, 10, len=100)
# expect_equal(grad(phiInv, x), phiInvPrime(x))

## CASE 2
# phi <- exp
# phiInv <- log
# phiInvPrime <- function(x) 1/x

## CASE 3
# phi <- function(x) -log(x)
# phiInv <- function(x) exp(-x)
# phiInvPrime <- function(x) -exp(-x)


#########################################################################

wqam <- compiler::cmpfun(function(x, w) {
   stopifnot(length(w) == length(x))
   expect_equal(sum(w), 1)
   expect_equal(phiInv(phi(x)), x)
   phiInv(sum(w*phi(x)))
})


wqamX <- compiler::cmpfun(function(w, X, phi, phiInv) {
   expect_equal(sum(w), 1)
#    as.numeric(phiInv(phiXt %*% w))
   as.numeric(phiInv(t(phi(X)) %*% w))
})


error_wqam <- compiler::cmpfun(function(w, X, Y, phi, phiInv) {
#    stopifnot(identical(all.equal(apply(X, 2, function(x) wqam(x, w)), as.double(phiInv(t(w)%*%phi(X)))), TRUE))
#    print(as.double(phiInv(t(w)%*%phi(X))-y))
#    print(as.double(t(w)%*%phi(X)-phi(y)))
#    stopifnot(identical(all.equal(as.double(phiInv(t(w)%*%phi(X))-y), as.double(t(w)%*%phi(X)-phi(y))), TRUE))
   Yobs <- wqamX(w, X, phi, phiInv) #apply(X, 2, function(x) wqam(x, w))
#    print(X)
#    print(Y)
   c(L1=sum(abs(  Yobs-Y   )   ),
     L2=sum(   (  Yobs-Y   )^2 )^0.5,
     LInf=max(abs(  Yobs-Y   )   )
   )
})

# error2 <- function(w) {
#    error_wqam(w)["L2"]
# }
#
# error1 <- function(w) {
#    error_wqam(w)["L1"]
# }

########################################################################


fit_wqam_L2_optim <- compiler::cmpfun(function(X, Y, phi, phiInv, phiInvPrime) {
   stopifnot(is.matrix(X), is.matrix(Y))
   n <- nrow(X); m <- ncol(X)
   stopifnot(1 == nrow(Y), m == ncol(Y))

   phiX <- phi(X)

   w0 <- runif(n); w0 <- w0/sum(w0)
   lambda0 <- log(w0)

   E <- function(lambda) {
      w <- exp(lambda)/sum(exp(lambda))
      sum((phiInv(t(w) %*% phiX)-Y)^2)
   }

   gradE <- function(lambda) {
      w <- exp(lambda)/sum(exp(lambda))
      Z <- as.numeric(t(w) %*% phiX)
      2*w*(
         ((phiInv(Z)-Y)*phiInvPrime(Z)) %*%
         (t(phiX) - Z)
      )
   }

#    max(replicate(1000, {w0 <- runif(n); w0 <- w0/sum(w0);   lambda0 <- log(w0);    sum((numDeriv::grad(E, lambda0) - gradE(lambda0))^2)}))

#    gradE <- function(lambda) {
#       w <- exp(lambda)/sum(exp(lambda))
#
# #       apply(
# #          t(sapply(1:n, function(k) {
# #          sapply(1:m, function(j) {
# #             2*(phiInv(sum(w*phiX[,j]))-Y[,j])*phiInvPrime(sum(w*phiX[,j]))*w[k]*(phiX[k,j]-sum(w*phiX[,j]))
# #          })
# #       })), 1, sum)
#
#       2*w*(
#          ((phiInv(t(w) %*% phiX)-Y)*phiInvPrime(t(w) %*% phiX)) %*% (t(phiX) -  as.numeric((t(w) %*% phiX)))
#       )
#
# #       2*w*(
# #          ((phiInv(t(w) %*% phiX)-Y)*phiInvPrime(t(w) %*% phiX)) %*% t(phiX)      -
# #          as.numeric(((phiInv(t(w) %*% phiX)-Y)*phiInvPrime(t(w) %*% phiX)) %*% t((t(w) %*% phiX)))
# #       )
#
#    }

   res <- optim(lambda0, E, gradE, method="BFGS", control=list(reltol=1e-16, maxit=10000))
   testthat::expect_equal(res$convergence, 0)
   exp(res$par)/sum(exp(res$par))
})



fit_wqam_L1_optim <- compiler::cmpfun(function(X, Y, phi, phiInv) {
   stopifnot(is.matrix(X), is.matrix(Y))
   n <- nrow(X); m <- ncol(X)
   stopifnot(1 == nrow(Y), m == ncol(Y))

   phiX <- phi(X)
   w0 <- runif(n); w0 <- w0/sum(w0)
   lambda0 <- log(w0)

   E <- function(lambda) {
      w <- exp(lambda)/sum(exp(lambda))
      sum(abs(phiInv(t(w) %*% phiX)-Y))
   }

   res <- optim(lambda0, E, NULL, method="BFGS", control=list(reltol=1e-16, maxit=10000))
   testthat::expect_equal(res$convergence, 0)
   exp(res$par)/sum(exp(res$par))
})

erf <- compiler::cmpfun(function(x) 2*pnorm(x*sqrt(2))-1)
erf_prime <- compiler::cmpfun(function(x) 2*exp(-x*x)/sqrt(pi))

fit_wqam_L1_optim_approx <- compiler::cmpfun(function(X, Y, phi, phiInv, phiInvPrime, eps) {
   stopifnot(is.matrix(X), is.matrix(Y))
   n <- nrow(X); m <- ncol(X)
   stopifnot(1 == nrow(Y), m == ncol(Y))

   phiX <- phi(X)
   w0 <- runif(n); w0 <- w0/sum(w0)
   lambda0 <- log(w0)

   E <- function(lambda) {
      w <- exp(lambda)/sum(exp(lambda))
      e <- phiInv(t(w) %*% phiX)-Y
      sum(sqrt(e^2+eps^2))
      # sum(abs(phiInv(t(w) %*% phiX)-Y))
   }

   gradE <- function(lambda) {
      w <- exp(lambda)/sum(exp(lambda))
      Z <- as.numeric(t(w) %*% phiX)
#       sapply(1:n, function(k) {
#          sum(sapply(1:m, function(j) {
#             (2*(phiInv(Z[j])-Y[,j])*phiInvPrime(Z[j])*w[k]*(phiX[k,j]-Z[j])) / (2*sqrt( (phiInv(Z[j])-Y[j])^2 + eps^2 ))
#          }))
#       })

      w*(
       ((phiInv(Z)-Y)*phiInvPrime(Z)/sqrt((phiInv(Z)-Y)^2 + eps^2)) %*% (t(phiX)-Z)
      )
   }

#    E <- function(lambda) {
#       w <- exp(lambda)/sum(exp(lambda))
#       e <- phiInv(t(w) %*% phiX)-Y
#       sum(e^2/sqrt(e^2+eps^2))
#    }
#
#    gradE <- NULL

#    E <- function(lambda) {
#       w <- exp(lambda)/sum(exp(lambda))
#       e <- phiInv(t(w) %*% phiX)-Y
#       sum(e*erf(e/eps))
#    }
#    gradE <- function(lambda) {
#       w <- exp(lambda)/sum(exp(lambda))
#       Z <- as.numeric(t(w) %*% phiX)
#       sapply(1:n, function(k) {
#          sum(sapply(1:m, function(j) {
#             erf((phiInv(Z[j])-Y[j])/eps)*phiInvPrime(phiInv(Z[j]))*w[k]*(phiX[k,j]-Z[j]) +
#             (phiInv(Z[j])-Y[j])*erf_prime((phiInv(Z[j])-Y[j])/eps)*phiInvPrime(phiInv(Z[j]))*w[k]*(phiX[k,j]-Z[j])/eps
#          }))
#       })
#    }

#    max(replicate(1000, {w0 <- runif(n); w0 <- w0/sum(w0);   lambda0 <- log(w0);    sum((numDeriv::grad(E, lambda0) - gradE(lambda0))^2)}))


   res <- optim(lambda0, E, gradE, method="BFGS", control=list(reltol=1e-16, maxit=10000))
   testthat::expect_equal(res$convergence, 0)
   exp(res$par)/sum(exp(res$par))
})


# error2wp <- function(wp) {
# #    Yobs <- phiInv(sapply(1:m, function(j) sum(wp*phiXmN[,j])+phiX[n,j]) )
#    phiXt   <- t(phi(X))
#    phiXmNt <- as.matrix(phiXt[,-n]-phiXt[,n])
#    phiXmN  <- t(phiXmNt)
#
#    Yobs <- as.numeric(phiInv(phiXmNt %*% wp + phiXt[,n]))
# #    sum( (Yobs-Y)^2 )
#    as.numeric((Yobs-Y) %*% t(Yobs-Y))
# }
#
# error2wp_grad <- function(wp) {
# #    arg <- sapply(1:m, function(j)
# #       sum(wp*phiXmN[,j])
# #    )+phiX[n,]
# #    sapply(1:(n-1), function(k) sum(2*(Y-phiInv(arg))*(-phiInvPrime(arg))*phiXmN[k,]))
#    phiXt   <- t(phi(X))
#    phiXmNt <- as.matrix(phiXt[,-n]-phiXt[,n])
#    phiXmN  <- t(phiXmNt)
#
#    arg <- as.numeric(phiXmNt %*% wp + phiXt[,n])
#    as.numeric(as.matrix((2*(Y-phiInv(arg))*(-phiInvPrime(arg)))) %*% phiXmNt)
# }
#
#
#
#
# # A <- rbind(rep(-1, n-1), diag(n-1))
# # B <- c(-1, rep(0, n-1))
# # startw <- runif(n-1)
# # startw <- startw/sum(startw)*0.9
# # res <- constrOptim(startw, error2wp, error2wp_grad, A, B, control=list(maxit=1000))
# #
# # print(res$convergence)
# # w <- c(res$par, 1-sum(res$par))
# # results[[length(results)+1]] <- list("MSE_optim", w, error(X, Y, w))
#
#
#

#
#
#
# error1_erf <- function(w) {
#    Yobs <- wqamX(w)#apply(X, 2, function(x) wqam(x, w))
# #    sum( (Yobs-Y)*erf((Yobs-Y)/eps) )
#    as.numeric((Yobs-Y) %*% t(erf((Yobs-Y)/eps)))
# }
#
# error1wp <- function(wp) {
#    phiXt   <- t(phi(X))
#    phiXmNt <- as.matrix(phiXt[,-n]-phiXt[,n])
#    phiXmN  <- t(phiXmNt)
#
#    Yobs <- as.numeric(phiInv(phiXmNt %*% wp + phiXt[,n]))
#    sum( abs(Yobs-Y) )
# }
#
# error1wp_erf <- function(wp) {
#    phiXt   <- t(phi(X))
#    phiXmNt <- as.matrix(phiXt[,-n]-phiXt[,n])
#    phiXmN  <- t(phiXmNt)
#
#    Yobs <- as.numeric(phiInv(phiXmNt %*% wp + phiXt[,n]))
# #    sum( (Yobs-Y)*erf((Yobs-Y)/eps) )
#    as.numeric((Yobs-Y) %*% t(erf((Yobs-Y)/eps)))
# }
#
#
# error1wp_erf_grad <- function(wp) {
#    phiXt   <- t(phi(X))
#    phiXmNt <- as.matrix(phiXt[,-n]-phiXt[,n])
#    phiXmN  <- t(phiXmNt)
#
# #    arg <- sapply(1:m, function(j)
# #       sum(wp*phiXmN[,j])
# #    )+phiX[n,]
#    arg <- as.numeric(phiXmNt %*% wp + phiXt[,n])
#
# #    as.numeric(
# #       as.matrix(
# #          -erf( (Y-phiInv(arg))/eps )*phiInvPrime(arg)-
# #          (Y-phiInv(arg))*erf_prime( (Y-phiInv(arg))/eps )*phiInvPrime(arg)/eps
# #       ) %*% phiXmNt
# #    )
#
#    -as.numeric(
#       as.matrix((erf( (Y-phiInv(arg))/eps )+(Y-phiInv(arg))*erf_prime( (Y-phiInv(arg))/eps )/eps)*phiInvPrime(arg))%*%phiXmNt
#    )
#
# #    as.numeric(sapply(1:(n-1), function(k) {
# #       -sum(
# #          (erf( (Y-phiInv(arg))/eps )+(Y-phiInv(arg))*erf_prime( (Y-phiInv(arg))/eps )/eps)*phiInvPrime(arg)*phiXmNt[,k]
# #       )
# #    }))
#
# #    as.numeric(sapply(1:(n-1), function(k) {
# #       sum(
# #          -erf( (Y-phiInv(arg))/eps )*phiInvPrime(arg)*phiXmNt[,k]+
# #          (Y-phiInv(arg))*erf_prime( (Y-phiInv(arg))/eps )*(-1)*phiInvPrime(arg)*phiXmNt[,k]/eps
# #       )
# #    }))
# }
#
#
# set.seed(12342)
# replicate(1000, {
#    w <- runif(n); w <- w/sum(w)
#    wp <- w[-n]
#    error1wp_erf_grad(wp)-numDeriv::grad(error1wp_erf, wp,
#       method.args=list(eps=1e-7, d=0.000001, zero.tol=sqrt(.Machine$double.eps/7e-7), r=4, v=2, show.details=FALSE))
# }) -> e
# print(apply(e, 1, range))
#
# set.seed(12342)
# replicate(1000, {
#    w <- runif(n); w <- w/sum(w)
#    wp <- w[-n]
#    error1(w)-error1_erf(w)
# }) -> e
# print(range(e))
#
# local({
#    set.seed(12345642)
#    replicate(1000, {
#       w <- runif(n); w <- w/sum(w)
#       wp <- w[-n]
#       expect_equal(length(wp), n-1)
#       expect_equal(w, c(wp, 1-sum(wp)))
#       stopifnot(all(wp >= 0))
#       stopifnot(sum(wp) <= 1)
#       expect_equal(error2(w), error2wp(wp))
#       expect_equal(sqrt(error2(w)), unname(error(w)["L2"]))
#       expect_equal(error1(w), unname(error(w)["L1"]))
#       expect_equal(error1_erf(w), error1wp_erf(wp))
#       expect_equal(error2wp_grad(wp), numDeriv::grad(error2wp, wp))
#       expect_equal(error1wp_erf_grad(wp), numDeriv::grad(error1wp_erf, wp,
#          method.args=list(eps=1e-4, d=0.000001, zero.tol=sqrt(.Machine$double.eps/7e-7), r=4, v=2, show.details=FALSE)))
#    })
# })
#
#
#
#
#
# # # # convex ????
# # x <- seq(0.18,0.22,by=0.005)
# # y <- seq(0.45,0.5,by=0.005)
# # z <- matrix(NA_real_, nrow=length(x), ncol=length(y))
# # for (i in seq_along(x))
# #    for (j in seq_along(y))
# #       z[i,j] <- if (x[i] + y[j] <= 1) error1wp(c(x[i], y[j])) else NA_real_
# # image(x, y, z)
# # contour(x, y, z, nlevels = 100, add=TRUE)
#
