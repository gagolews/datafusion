library('geometry')
library(copula)
library('tripack')

pdf("~/Publications/Books/Habilitacja/figures/delaunay.pdf", height=4, width=6)
par(mar=c(2.5, 2.5, 0.5, 0.5))
set.seed(43767218)
n <- 5
copula <- frankCopula(dim=2, 1.2)
# copula <- plackettCopula(5)
X <- t(rCopula(n, copula))
X[1,] <- qt(X[1,], 10)
X[2,] <- qt(X[2,], 14)
plot(X[1,], X[2,], las=1, xlab=NA, ylab=NA, asp=1, xlim=c(-2,4), ylim=c(-4,2))


T <- delaunayn(t(X))
for (i in 1:nrow(T)) {
   polygon(X[1,T[i,]], X[2,T[i,]], lty=3, lwd=2)
   C <- circumcircle(X[1,T[i,]], X[2,T[i,]], plot=FALSE, num=3)
   circles(C$x, C$y, C$radius, col='gray')
}

dev.off()

# T2 <- tri.mesh(X[1,],X[2,])
# plot.tri(T2)
