library(testthat)
options(stringsAsFactors=FALSE)
source("fit_wam.R")
source("fit_wam_outrank.R")
Sys.setenv(PKG_LIBS="-lCGAL")
Rcpp::sourceCpp("cgal_qp_solver.cpp")

###########################################################################

set.seed(15464)
n <- 2
m <- 10
X <- t(round(matrix(runif(n*m, 0, 1), nrow=m), 2))

realw <- c(0.23, 0.77)
runif(n)
realw <- realw/sum(realw)
Y <- matrix(apply(X, 2, function(x) wam(x, realw)), ncol=m)   # ERROR = 0
Y <- round(Y+rnorm(m, 0, 0.05), 2)                            # PERTURBATION



set.seed(19)
n <- 5
m <- 9
X <- t(round(matrix(runif(n*m, 0, 1), nrow=m), 2))

realw <-  runif(n)#c(0.23, 0.77) #
realw <- realw/sum(realw)
realw <- round(realw, 2)
stopifnot(sum(realw) == 1)
Y <- matrix(apply(X, 2, function(x) wam(x, realw)), ncol=m)+rnorm(m, 0.05, 0.1)
Y <- round(Y, 2)                            # PERTURBATION
Y <- pmin(pmax(Y, apply(X, 2, min)), apply(X, 2, max))



stopifnot(dim(Y) == c(1, m))
stopifnot(dim(X) == c(n, m))
expect_equal(wamX(realw), apply(X, 2, function(x) wam(x, realw)))

###########################################################################

results <- matrix(NA_real_, nrow=0, ncol=n+3+1, dimnames=list(NULL, c(paste0("w", 1:n), "L1", "L2", "LInf", "Kendall")))

###### 1. use realw

w <- realw
results <- rbind(results, exact=c(w, error_wam(w), cor(order(t(w) %*% X), order(Y), method="kendall")))

###### 2. L1 linprog

# w <- fit_wam_L1_linprog_old(X, Y)
# results <- rbind(results, L1_lp_old=c(w, error_wam(w), cor(order(t(w) %*% X), order(Y), method="kendall")))

w <- fit_wam_L1_linprog(X, Y)
results <- rbind(results, L1_lp=c(w, error_wam(w), cor(order(t(w) %*% X), order(Y), method="kendall")))

###### 3. L1 penalty outrank

par(mfrow=c(2,1))
p <- seq(0, 20, 0.05)[-1]
bp <- sapply(p, function(p) {
   w <- fit_wam_L1_linprog_outrank(X, Y, p)
   c(cor=cor(order(t(w) %*% X), order(Y), method="kendall"), error_wam(w)["L1"])
})
plot(p, bp[1,], pch='.', ylim=c(0, max(bp)))
points(p, bp[2,], pch='.', col=2)
w <- which(bp["cor",] == max(bp["cor",], na.rm=TRUE))
bp <- bp[,w]
p <- p[w]
p <- p[which.min(bp["L1",])]

# p <- p[which.max(bp["cor",])]
cat(sprintf("fit_wam_L1_linprog_outrank: p=%g, #NA=%g\n", p, sum(is.na(bp[1,]))))
abline(v=p, lty=2)
w <- fit_wam_L1_linprog_outrank(X, Y, p)
results <- rbind(results, L1_outrank=c(w, error_wam(w), cor(order(t(w) %*% X), order(Y), method="kendall")))

# w <- fit_wam_L1_linprog_outrank_old(X, Y, p)
# results <- rbind(results, L1_outrank_old=c(w, error_wam(w), cor(order(t(w) %*% X), order(Y), method="kendall")))

###### 4. L2 quadprog

# w <- fit_wam_L2_quadprog_old(X, Y)
# results <- rbind(results, L2_qp_old=c(w, error_wam(w), cor(order(t(w) %*% X), order(Y), method="kendall")))


w <- fit_wam_L2_quadprog(X, Y)
results <- rbind(results, L2_qp=c(w, error_wam(w), cor(order(t(w) %*% X), order(Y), method="kendall")))

###### 5. L2 penalty outrank

p <- seq(0, 10, 0.01)[-1]
bp <- sapply(p, function(p) {
   tryCatch({
      w <- fit_wam_L2_quadprog_outrank(X, Y, p)
      c(cor=cor(order(t(w) %*% X), order(Y), method="kendall"), error_wam(w)["L2"])
   }, error=function(e) c(NA, NA))
})
plot(p, bp[1,], pch='.', ylim=c(0, max(bp, na.rm = TRUE)))
points(p, bp[2,], pch='.', col=2)

w <- which(bp["cor",] == max(bp["cor",], na.rm=TRUE))
bp <- bp[,w]
p <- p[w]
p <- p[which.min(bp["L2",])]

# p <- p[which.max(bp["cor",])]
cat(sprintf("fit_wam_L2_quadprog_outrank: p=%g, #NA=%g\n", p, sum(is.na(bp[1,]))))
abline(v=p, lty=2)
w <- fit_wam_L2_quadprog_outrank(X, Y, p)
results <- rbind(results, L2_outrank=c(w, error_wam(w), cor(order(t(w) %*% X), order(Y), method="kendall")))

###### 6. LInf linprog

# w <- fit_wam_LInf_linprog_old(X, Y)
# results <- rbind(results, LInf_lp_old=c(w, error_wam(w), cor(order(t(w) %*% X), order(Y), method="kendall")))

w <- fit_wam_LInf_linprog(X, Y)
results <- rbind(results, LInf_lp=c(w, error_wam(w), cor(order(t(w) %*% X), order(Y), method="kendall")))

############### OUTPUT #############################

results <- cbind(results, L1rel=(results[,"L1"]-min(results[,"L1"]))/min(results[,"L1"]))
results <- cbind(results, L2rel=(results[,"L2"]-min(results[,"L2"]))/min(results[,"L2"]))
results <- cbind(results, LInfrel=(results[,"LInf"]-min(results[,"LInf"]))/min(results[,"LInf"]))
old_scipen <- options("scipen")$scipen
options(scipen=16)
print(results, digits=4, scientific = FALSE)
options(scipen=old_scipen)

# xtable::xtable(rbind(X, Y), digits=2)
