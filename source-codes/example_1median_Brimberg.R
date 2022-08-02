# set.seed(1234)
library('copula')
# set.seed(126)
n <- 3
d <- 2
copula <- frankCopula(dim=d, 1.2)
X <- t(rCopula(n, copula))
X[1,] <- qexp(X[1,], 0.1)
X[2,] <- qt(X[2,], 1)
# X <- matrix(rcauchy(n*d), nrow=d)

source('nfusion_invariance_test.R')
source('nfusion_1median_Brimberg.R')

w  <- rep(1/n, n)
y0 <- rowMeans(X)
cat(sprintf("translation %s\n", paste(fivenum(replicate(250, norm(as.matrix(check_translation(X, Weiszfeld1median, w,y0))))), collapse=', ')))
cat(sprintf("scaling     %s\n", paste(fivenum(replicate(250, norm(as.matrix(check_scaling(X, Weiszfeld1median, w,y0))))), collapse=', ')))
cat(sprintf("scalingn    %s\n", paste(fivenum(replicate(250, norm(as.matrix(check_scalingn(X, Weiszfeld1median, w,y0))))), collapse=', ')))
cat(sprintf("monotone    %s\n", paste(fivenum(replicate(250, norm(as.matrix(check_monotone(X, Weiszfeld1median, w,y0))))), collapse=', ')))
cat(sprintf("monotone2   %s\n", paste(fivenum(replicate(250, norm(as.matrix(check_monotone2(X, Weiszfeld1median, w,y0))))), collapse=', ')))
cat(sprintf("monotone3   %s\n", paste(fivenum(replicate(250, norm(as.matrix(check_monotone3(X, Weiszfeld1median, w,y0))))), collapse=', ')))
cat(sprintf("orhogonal   %s\n", paste(fivenum(replicate(250, norm(as.matrix(check_orthogonal(X, Weiszfeld1median, w,y0))))), collapse=', ')))
cat(sprintf("affine      %s\n", paste(fivenum(replicate(250, norm(as.matrix(check_affine(X, Weiszfeld1median, w,y0))))), collapse=', ')))


d <- function(x, y) sqrt(sum((x-y)^2))

Weiszfeld1median_old <- function(X, w, y0, eps=1e-9) {
   eps <- 1e-12
   y <- y0
   repeat {
      yp <- y
      dy <- sapply(1:ncol(X), function(i) d(X[,i], y))
      if (any(dy < eps))
         y <- X[which.min(dy),]
      else {
         wdy <- (w/dy)
         y <- sapply(1:nrow(X), function(i) sum(wdy*X[i,]))
         y <- y/sum(w/dy)
      }
      if (d(yp, y) < eps) break
   }
   y
}

y0 <- rowMeans(X)
w  <- rep(1/n, n)


yoptim <- optim(y0, function(y) sum(sapply(1:ncol(X), function(i) d(X[, i], y))), control=list(abstol=1e-12))$par
yrcpp <- Weiszfeld1median(X, w, y0, 1e-10)
# yold <- Weiszfeld1median_old(X, w, y0)

print(d(yoptim,yrcpp))
# print(d(yoptim,yold))
# print(sum(sapply(1:nrow(X), function(i) d(X[,i], y)))-sum(sapply(1:nrow(X), function(i) d(X[,i], yoptim))))

# pdf("~/Publications/Books/habilitacja/figures/FermatWeber.pdf", height=4, width=6)
par(mar=c(2.5, 2.5, 0.5, 0.5))
plot(X[1,], X[2,], las=1, xlab=NA, ylab=NA)
ch <- chull(X[1,],X[2,])
# polygon(X[1,ch], X[2,ch])

points(mean(X[1,]), mean(X[2,]), col=4, pch=2, cex=2)
points((yrcpp[1]), (yrcpp[2]), col=2, pch=8, cex=2)
# points(mean(y0[1]), mean(y0[2]), col=4, pch=5)
# dev.off()

## IS IT AN AGOP???
# set.seed(126)
# n <- 3
# d <- 2
# X <- matrix(byrow=TRUE, nrow=d,
#    c(0, 0, 1,
#      0, 1, 0))
# m <- Weiszfeld1median(X, rep(1/ncol(X), ncol(X)), rowMeans(X), eps=1e-12)
# replicate(100, {
#    i <- sample(1:ncol(X), 1)
#    t <- runif(nrow(X), 0, 100)
#    Y <- X
#    Y[,i] <- Y[,i]+t
#    c(Weiszfeld1median(Y, rep(1/ncol(X), ncol(X)), m, eps=1e-12)-m, i, t)
# }) -> res
# summary(t(res[1:2,]))
# res
#
# Y <- X
# Y[,3] <- Y[,3]+c(2,2)
# plot(c(X[1,], Y[1,]), c(X[2,], Y[2,]), las=1, xlab=NA, ylab=NA, pch=rep(as.character(1:3), times=2), col=rep(1:2, each=3))
# m1 <- Weiszfeld1median(X, rep(1/ncol(X), ncol(X)), rowMeans(X), eps=1e-12)
# m2 <- Weiszfeld1median(Y, rep(1/ncol(Y), ncol(Y)), rowMeans(Y), eps=1e-12)
# points((m1[1]), (m1[2]), col=2, pch=8, cex=2)
# points((m2[1]), (m2[2]), col=4, pch=8, cex=2)
#
# yoptim1 <- optim(y0, function(y) sum(sapply(1:ncol(X), function(i) d(X[, i], y))), control=list(abstol=1e-12))$par
# yoptim2 <- optim(y0, function(y) sum(sapply(1:ncol(Y), function(i) d(Y[, i], y))), control=list(abstol=1e-12))$par
#
#
# x <- seq(0,1,length=100)
# y <- seq(0,2,length=100)
# z <- matrix(NA_real_, nrow=length(x), ncol=length(y))
# for (i in 1:length(x))
#    for (j in 1:length(y))
#       z[i,j] <- sum(sapply(1:ncol(Y), function(k) d(X[, k], c(x[i], y[j]))))
#
# contour(x,y,z,nlevels = 100)
