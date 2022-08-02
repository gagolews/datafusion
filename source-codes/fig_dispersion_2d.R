# x <- matrix(byrow=TRUE, nrow=2, c(
#    1, 0,
#    0, 1,
#    0.5, 0))
#
# plot(x[1,], x[2,])
#
# apply(x, 2, function(y) t(as.matrix(c(0.5,0.5)/sqrt(sum(c(0.5,0.5)^2))))%*%as.matrix(y))

set.seed(123)
n <- 100
d <- 2
copula <- copula::frankCopula(dim=d, 3)
X <- t(copula::rCopula(n, copula))
X[1,] <- qexp(X[1,])*4
X[2,] <- qt(X[2,], 4)-100

X[1,] <-X[1,]-mean(X[1,])
X[2,] <-X[2,]-mean(X[2,])
theta <- pi/6
X <-  matrix(byrow=TRUE, ncol=2, c(
  cos(theta), sin(theta),
   -sin(theta), cos(theta)
)) %*% X


pdf("figures/dispersion_pursuit1.pdf", height=6, width=6)
par(mar=c(2.5, 2.5, 0.5, 0.5))
plot(X[1,], X[2,], asp=1, las=1, xlab=NA, ylab=NA)
dev.off()


theta <- seq(0,2*pi,length.out=999)[-1]
Y <- matrix(NA_real_, nrow=5, ncol=length(theta))
for (i in seq_along(theta)) {
   Y[1:2,i] <- y <- as.matrix(c(cos(theta[i]), sin(theta[i])))
   x <- as.numeric(apply(X, 2, function(x) t(y)%*%x))
   Y[3,i] <- 1.4826*median(abs(x-median(x)))
   Y[4,i] <- sd(x)
   Y[5,i] <- diff(quantile(x, c(0.25, 0.75))) #mean(abs(x-mean(x)))*sqrt(pi/2)
}
theta <- seq(0,2*pi,length.out=31)[-1]
Z <- matrix(NA_real_, nrow=5, ncol=length(theta))
for (i in seq_along(theta)) {
   Z[1:2,i] <- y <- as.matrix(c(cos(theta[i]), sin(theta[i])))
   x <- as.numeric(apply(X, 2, function(x) t(y)%*%x))
   Z[3,i] <- 1.4826*median(abs(x-median(x)))
   Z[4,i] <- sd(x)
   Z[5,i] <- diff(quantile(x, c(0.25, 0.75))) #mean(abs(x-mean(x)))*sqrt(pi/2)
}


for (i in 3:5) {
   pdf(sprintf("figures/dispersion_pursuit%d.pdf", i), height=6, width=6)
   par(mar=c(2.5, 2.5, 0.5, 0.5))
   plot(NA, NA, xlim=range(Y[1,]*Y[i,]), ylim=range(Y[2,]*Y[i,]), asp=1, las=1, xlab=NA, ylab=NA)
   polygon(Y[1,]*Y[i,], Y[2,]*Y[i,], col='gray80', border='gray60')
   arrows(0, 0, Z[1,]*Z[i,], Z[2,]*Z[i,], col=1, asp=1)
   dev.off()
}


# mean(Y[3,])
