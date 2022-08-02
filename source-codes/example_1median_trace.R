library('copula')
set.seed(123)
n <- 25
d <- 2
copula <- claytonCopula(1.2)
X <- t(rCopula(n, copula))
X[1,] <- qt(X[1,], 3)
X[2,] <- qexp(X[2,], 0.06)-20

pdf("figures/1median_trace.pdf", height=4, width=6)
par(mar=c(2.5, 2.5, 0.5, 0.5))
palette(RColorBrewer::brewer.pal(5, "Dark2"))
plot(X[1,], X[2,], las=1, xlab=NA, ylab=NA)

u <- chull(t(X))
polygon(X[1,u], X[2,u], border='gray', lty=3)

one_median_Lp <- function(X, p) {
   optim(X[,13], function(c) {
      sum(colSums(abs(X-c)^p)^(1/p))
   }, method="BFGS", control=list(reltol=1e-16))$par
}

P <- c(seq(1, 25, by=0.1), 100)
sapply(P, function(p) {
   y <- one_median_Lp(X, p)
   points(y[1], y[2], pch=16, col=2)
   y
}) -> res

pi <- which(P == 1)
arrows(res[1,pi]-0.1, -18, res[1,pi], res[2,pi], lty=2, col=2)
text(res[1,pi]-0.1, -18, expression(d[1]), col=2, cex=1.2)

pi <- which(P == 2)
arrows(res[1,pi]+0.5, -18, res[1,pi], res[2,pi], lty=2, col=2)
text(res[1,pi]+0.5, -18, expression(d[2]), col=2, cex=1.2)

pi <- which(P == 100)
arrows(res[1,pi]+0.5, -18, res[1,pi], res[2,pi], lty=2, col=2)
text(res[1,pi]+0.5, -18, expression(d[infinity]), col=2, cex=1.2)


# for (p in c(seq(1, 25, by=0.1), 26:100)) {
#    r <- optim(X[,1], function(c) {
#       sum(colSums(abs(X-c)^p))
#    }, method="BFGS", control=list(reltol=1e-16))
#    stopifnot(r$convergence==0)
#    y <- r$par
#    points(y[1], y[2], pch=16, col=3)
# }
#
# r <- optim(X[,1], function(c) {
#    max(apply(abs(X-c), 2, max))
# }, method="BFGS", control=list(reltol=1e-16))
# stopifnot(r$convergence==0)
# y <- r$par
# points(y[1], y[2], pch=16, col=3)
dev.off()
