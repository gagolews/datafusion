library(testthat)
options(stringsAsFactors=FALSE)
source("fit_wam.R")
source("fit_wam_outrank.R")
source("fit_wam_regularization.R")
Sys.setenv(PKG_LIBS="-lCGAL")
Rcpp::sourceCpp("cgal_qp_solver.cpp")

###########################################################################

set.seed(321)
n <- 10
m <- 100
realw <- runif(n)
realw <- realw/sum(realw)
X <- t(round(matrix(runif(n*m, 0, 1), nrow=m), 2))
Y <- t(realw) %*% X + rnorm(m, 0, 0.05)

train <- sample(1:m, m*0.8)
X_test <- X[,-train,drop=FALSE] # test sample
Y_test <- Y[,-train,drop=FALSE]
X <- X[,train,drop=FALSE]       # training sample
Y <- Y[,train,drop=FALSE]

error_wam_test <- function(w) {
   Yobs <- as.numeric(t(X_test) %*% w)
   c(L1=sum(abs(  Yobs-Y_test   )   ),
     L2=sum(   (  Yobs-Y_test   )^2 )^0.5,
     LInf=max(abs(  Yobs-Y_test   )   )
   )
}

###########################################################################

results <- matrix(NA_real_, nrow=0, ncol=n+3, dimnames=list(NULL, c(paste0("w", 1:n), "L1", "L2", "LInf")))

###### 1. use realw

w <- realw
results <- rbind(results, exact=c(w, error_wam_test(w)))

###### 2. L1 linprog

# w <- fit_wam_L1_linprog_old(X, Y)
# results <- rbind(results, L1_lp_old=c(w, error_wam(w), cor(order(t(w) %*% X), order(Y), method="kendall")))

w <- fit_wam_L1_linprog(X, Y)
results <- rbind(results, L1_lp=c(w, error_wam_test(w)))

# w <- fit_wam_L1_linprog(X_test, Y_test)
# results <- rbind(results, L1_lp_test=c(w, error_wam_test(w)))

###### 4. L2 quadprog

# w <- fit_wam_L2_quadprog_old(X, Y)
# results <- rbind(results, L2_qp_old=c(w, error_wam(w), cor(order(t(w) %*% X), order(Y), method="kendall")))

w <- fit_wam_L2_quadprog(X, Y)
results <- rbind(results, L2_qp=c(w, error_wam_test(w)))

# w <- fit_wam_L2_quadprog(X_test, Y_test)
# results <- rbind(results, L2_qp_test=c(w, error_wam_test(w)))

###### 5. L2 quadprog+regularization

p <- seq(-5, 5, 0.005)
bp <- sapply(p, function(p) {
   w <- fit_wam_L2_regular(X, Y, p)
   error_wam_test(w)
})
# cairo_pdf("figures/wam_regularization_L2.pdf", height=4, width=6)
par(mar=c(4, 4, 0.5, 0.5))
matplot(p, t(bp), type='l', log="y", xlab=expression(lambda), ylab="training sample error", las=1, col=c(1,2,4), lty=c(1,2,4), lwd=c(1,2,1))
legend("top", horiz=TRUE, col=c(1,2,4), lty=c(1,2,4), lwd=c(1,2,1), expression(L[1], L[2], L[infinity]))
# dev.off()
p <- p[which.min(bp["L2",])]
cat(sprintf("fit_wam_L2_regular: p=%g\n", p))
# abline(v=p, lty=2)

w <- fit_wam_L2_regular(X, Y, p)
results <- rbind(results, L2_qp_regular=c(w, error_wam_test(w)))

###### 6. LInf linprog

# w <- fit_wam_LInf_linprog_old(X, Y)
# results <- rbind(results, LInf_lp_old=c(w, error_wam(w), cor(order(t(w) %*% X), order(Y), method="kendall")))

w <- fit_wam_LInf_linprog(X, Y)
results <- rbind(results, LInf_lp=c(w, error_wam_test(w)))

############### OUTPUT #############################

results <- cbind(results, L1rel=results[,"L1"]-min(results[,"L1"]))
results <- cbind(results, L2rel=results[,"L2"]-min(results[,"L2"]))
results <- cbind(results, LInfrel=results[,"LInf"]-min(results[,"LInf"]))
old_scipen <- options("scipen")$scipen
options(scipen=16)
print(results, digits=4, scientific = FALSE)
options(scipen=old_scipen)

# xtable::xtable(rbind(X, Y), digits=2)
