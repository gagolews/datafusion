# set.seed(1234)
library('copula')
# set.seed(126)
n <- 25
d <- 4
copula <- claytonCopula(dim=d, 1.2)
X <- t(rCopula(n, copula))
X[1,] <- qexp(X[1,], 10)*20
X[2,] <- qt(X[2,], 3)-100
if (d >= 3) X[3,] <- qt(X[3,], 1)+10
if (d >= 4) X[4,] <- qunif(X[4,],-100, 100)
# X <- matrix(rcauchy(n*d), nrow=d)

# X[1,] <- rnorm(n)
# X[2,] <- rnorm(n, 0, 0.1)

# plot(X[1,], X[2,], asp=1)

# Y <- rortho(2)%*%X

# Y <- X
# Y <- diag(c(2, 0.5)) %*% X

# SVD <- svd(t(Y), n)
# u <- SVD$u
# d <- SVD$d
# v <- SVD$v

# print(u)
# print(d)
# print(v)

# plot(u[,1], u[,2])

# plot(Y[1,], Y[2,])
# z <- as.matrix(c(median(u[,1]), median(u[,2])))
# print(z)
# print(zz <- t(z) %*% diag(d) %*% t(v))
# points(zz[1], zz[2], col=2)

f2 <- function(x) sign(x[1])*max(pmin(1:length(x),sort(abs(x), decreasing=TRUE)))
# f2 <- function(x)            max(pmin(1:length(x),sort(abs(x), decreasing=TRUE)))
# f2 <- median
# f2 <- function(x) sum(x^4)^(1/4)
# f2 <- function(x) quantile(x, 0.74)
# f2 <- function(x) log(sum(exp(x)))
# w <- rexp(n); w <- w/sum(w); f2 <- function(x) sum(w*x)
# w <- abs(rcauchy(n)); f2 <- function(x) sum(w*x)
# w <- rexp(n); w <- w/sum(w); f2 <- function(x) sum(w*sort(x))
# w <- abs(rcauchy(n)); f2 <- function(x) sum(w*sort(x))
# f2 <- function(x) if (max(x)>-min(x)) 100 else -100

comedortho2a <- compiler::cmpfun(function(X, f=f2) {
   centers <- rowMeans(X)
   SVD <- svd(t(X-centers))
   s <- sqrt(sum(SVD$d^2))

   as.numeric(s * (SVD$v %*% as.matrix(apply(t(SVD$v) %*% (X-centers)/s, 1, f)))+centers)
})


comedortho2b <- compiler::cmpfun(function(X, f=f2) {
   centers <- rowMeans(X)
   SVD <- svd(t(X-centers))

   as.numeric(t(as.matrix(apply(SVD$u, 2, f) %*% diag(SVD$d) %*% t(SVD$v))))+centers
   # as.numeric(SVD$v %*% diag(SVD$d) %*% apply(SVD$u, 2, f) + centers)
})


comedortho2 <- comedortho2a

# X <- X-rowMeans(X)
# X <- t(svd(t(X))$v) %*% X
plot(X[1,], X[2,], asp=1)
print(var(X[1,]))
print(var(X[2,]))
print(res <- comedortho2(X))
points(res[1], res[2], pch=3, col='red')

source('nfusion_invariance_test.R')
cat(sprintf("translation %s\n", paste(fivenum(replicate(100, norm(as.matrix(check_translation(X, comedortho2))))), collapse=', ')))
cat(sprintf("scaling     %s\n", paste(fivenum(replicate(100, norm(as.matrix(check_scaling(X, comedortho2))))), collapse=', ')))
cat(sprintf("scalingn    %s\n", paste(fivenum(replicate(100, norm(as.matrix(check_scalingn(X, comedortho2))))), collapse=', ')))
cat(sprintf("monotone    %s\n", paste(fivenum(replicate(100, norm(as.matrix(check_monotone(X, comedortho2))))), collapse=', ')))
cat(sprintf("orhogonal   %s\n", paste(fivenum(replicate(100, norm(as.matrix(check_orthogonal(X, comedortho2))))), collapse=', ')))
cat(sprintf("affine      %s\n", paste(fivenum(replicate(100, norm(as.matrix(check_affine(X, comedortho2))))), collapse=', ')))
cat(sprintf("affine2     %s\n", paste(fivenum(replicate(100, norm(as.matrix(check_affine2(X, comedortho2))))), collapse=', ')))



comed_ddcs <- function(X, f=f2) {
   d <- nrow(X)
   A <- X[1:d, 1:d]
   Y <- solve(A) %*% X
   centers <- rowMeans(Y)
   as.numeric(A%*%(apply(Y, 1, f)))

   centers <- rowMeans(X)
   Xc <- X-centers
   d <- nrow(Xc)
   A <- X[1:d, 1:d]-X[,d+1]
   Y <- solve(A) %*% Xc
   as.numeric(A%*%(apply(Y, 1, f)))+centers

#    centers <- rowMeans(X)
#    Xc <- X-centers
#    d <- nrow(Xc)
#    A <- Xc[1:d, 1:d]
#    Y <- solve(A) %*% Xc
#    as.numeric(A%*%(apply(Y, 1, f)))+centers
}

plot(X[1,], X[2,], asp=1)
print(res <- comed_ddcs(X))
points(res[1], res[2], pch=3, col='red')

source('nfusion_invariance_test.R')
cat(sprintf("translation %s\n", paste(fivenum(replicate(100, norm(as.matrix(check_translation(X, comed_ddcs))))), collapse=', ')))
cat(sprintf("scaling     %s\n", paste(fivenum(replicate(100, norm(as.matrix(check_scaling(X, comed_ddcs))))), collapse=', ')))
cat(sprintf("scalingn    %s\n", paste(fivenum(replicate(100, norm(as.matrix(check_scalingn(X, comed_ddcs))))), collapse=', ')))
cat(sprintf("monotone    %s\n", paste(fivenum(replicate(100, norm(as.matrix(check_monotone(X, comed_ddcs))))), collapse=', ')))
cat(sprintf("orhogonal   %s\n", paste(fivenum(replicate(100, norm(as.matrix(check_orthogonal(X, comed_ddcs))))), collapse=', ')))
cat(sprintf("affine      %s\n", paste(fivenum(replicate(100, norm(as.matrix(check_affine(X, comed_ddcs))))), collapse=', ')))
cat(sprintf("affine2     %s\n", paste(fivenum(replicate(100, norm(as.matrix(check_affine2(X, comed_ddcs))))), collapse=', ')))


#
#
# comedortho <- compiler::cmpfun(function(X) {
#    centers <- rowMeans(X)
#    Xr <- rep(1, nrow(X))#apply(apply(X-centers, 1, range), 2, diff)
#    #Xr <- apply(X-centers, 1, sd)
#    XtXev <- svd(t((X-centers)/Xr))$v
#
#    Y <- t(t((X-centers)/Xr)%*%XtXev)
#
#    y <- apply(Y, 1, function(x) median(x))
#
#    as.numeric(t((y)%*%solve(XtXev))*Xr+centers)
# })
#
#
#
# print(rowMeans(X))
# print(comedortho(X))
#
#
#
# cat(sprintf("translation %s\n", paste(fivenum(replicate(250, norm(as.matrix(check_translation(X, comedortho))))), collapse=', ')))
# cat(sprintf("scaling     %s\n", paste(fivenum(replicate(250, norm(as.matrix(check_scaling(X, comedortho))))), collapse=', ')))
# cat(sprintf("orhogonal    %s\n", paste(fivenum(replicate(250, norm(as.matrix(check_orthogonal(X, comedortho))))), collapse=', ')))


# comean <- function(X) { rowMeans(X) }
#
# cat(sprintf("translation %s\n", paste(fivenum(replicate(250, norm(as.matrix(check_translation(X, comean))))), collapse=', ')))
# cat(sprintf("scaling     %s\n", paste(fivenum(replicate(250, norm(as.matrix(check_scaling(X, comean))))), collapse=', ')))
# cat(sprintf("orhogonal    %s\n", paste(fivenum(replicate(250, norm(as.matrix(check_orthogonal(X, comean))))), collapse=', ')))


# # IS IT AN AGOP?
# set.seed(126)
# n <- 4
# d <- 2
# X <- matrix(byrow=TRUE, nrow=d, c(
#    1, -1,   1, -1,
#    -1, 1, 1,-1
#    ))
# m <- comedortho(X)
#
# replicate(1000, {
#    i <- sample(1:ncol(X), 1)
#    t <- runif(nrow(X), 0, 100)
#    Y <- X
#    Y[,i] <- Y[,i]+t
#    c(comedortho(Y)-m, i, t)
# }) -> res
# summary(t(res[1:2,]))
#
# Y <- X
# Y[,3] <- Y[,3]+c(0,2)
# m2 <- comedortho(Y)
# print(m)
# print(m2)

