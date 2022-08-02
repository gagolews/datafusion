# set.seed(1234)
library('copula')
# set.seed(126)
n <- 3
d <- 2

# Weak monotonicity in Simon's sense

w  <- rep(1/n, n)
y0 <- rowMeans(X)


source('nfusion_invariance_test.R')
source('nfusion_1median_Brimberg.R')

X <- matrix(nrow=2,byrow=TRUE,
   c(0,1,20,
     0,-5,1))
plot(X[1,],X[2,],asp=1)


X2 <- matrix(nrow=2,byrow=TRUE,
   c(0,1,2000,
     0,-5,2))
plot(X[1,],X[2,],asp=1)
print(val1 <- Weiszfeld1median(X, w, y0, 1e-16))
print(val2 <- Weiszfeld1median(X2, w, y0, 1e-15))
print(colSums(t(X-val1)/   apply(X, 2, function(x) sqrt(sum((x-val1)^2)))))
print(colSums(t(X2-val2)/   apply(X2, 2, function(x) sqrt(sum((x-val2)^2)))))
stopifnot((any(val1<=val2)))



for (i in 1:10000) {
   # X <- matrix(rcauchy(n*d), nrow=d)
   # copula <- frankCopula(dim=d, 1.2)
# X <- t(rCopula(n, copula))
# X[1,] <- qexp(X[1,], 0.1)
# X[2,] <- qt(X[2,], 1)

   w  <- rep(1/n, n)
   y0 <- rowMeans(X)

   i <- sample(1:ncol(X), 1)
   X2 <- X
   X2[,i] <- X2[,i] + runif(nrow(X2), 0, 1000)
   val2 <- Weiszfeld1median(X2, w, y0)
   val1 <- Weiszfeld1median(X, w, y0)

   # weak monotonicity by Simon James
   if (all(val1>val2))
      m <- sum(val1-val2)
   else
      m <- 0

   stopifnot(m < 1)
}

