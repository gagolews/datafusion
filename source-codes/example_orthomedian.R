# data set
set.seed(12347)
n <- 9
X <- matrix(rnorm(2*n), nrow=2, byrow=TRUE)
par(mar=c(2.5, 2.5, 0.5, 0.5))
X[1,] <- X[1,]-median(X[1,])
X[2,] <- X[2,]-median(X[2,])


source('nfusion_orthomedian.R')

ormed <- orthomedian2(X)

par(mar=c(2.5, 2.5, 0.5, 0.5))
plot(X[1,], X[2,], xlim=range(X[1,]), ylim=range(X[2,]), las=1, xlab=NA, ylab=NA, asp=1)
points(ormed[1], ormed[2], col=2, pch=8, cex=2)

source('nfusion_invariance_test.R')
cat(sprintf("translation %s\n", paste(fivenum(replicate(250, norm(as.matrix(check_translation(X, orthomedian2))))), collapse=', ')))
cat(sprintf("scaling     %s\n", paste(fivenum(replicate(250, norm(as.matrix(check_scaling(X, orthomedian2))))), collapse=', ')))
cat(sprintf("scalingn    %s\n", paste(fivenum(replicate(250, norm(as.matrix(check_scalingn(X, orthomedian2))))), collapse=', ')))
cat(sprintf("monotone    %s\n", paste(fivenum(replicate(250, norm(as.matrix(check_monotone(X, orthomedian2))))), collapse=', ')))
cat(sprintf("monotone2   %s\n", paste(fivenum(replicate(250, norm(as.matrix(check_monotone2(X, orthomedian2))))), collapse=', ')))
cat(sprintf("monotone3   %s\n", paste(fivenum(replicate(250, norm(as.matrix(check_monotone3(X, orthomedian2))))), collapse=', ')))
cat(sprintf("orhogonal   %s\n", paste(fivenum(replicate(250, norm(as.matrix(check_orthogonal(X, orthomedian2))))), collapse=', ')))
cat(sprintf("affine      %s\n", paste(fivenum(replicate(250, norm(as.matrix(check_affine(X, orthomedian2))))), collapse=', ')))


# # how much an orthomedian differs from 1-median?
# d <- 2
# replicate(100, {
#    n <- 10
#    copula <- normalCopula(dim=d, runif(1, 0, 1))
#    X <- t(rCopula(n, copula))
#    X[1,] <- qexp(X[1,])
#    X[2,] <- qt(X[2,], 1)
#
#    m3 <- rowMeans(X)
#    m1 <- Weiszfeld1median(X, rep(1/n, n), m3)
#    m2 <- orthomedian2(X, 1000)
#
#    c(sqrt(sum((m1-m2)^2)), sqrt(sum((m1-m3)^2)), sqrt(sum((m3-m2)^2)))
# }) -> dst
# apply(dst,1,fivenum)
