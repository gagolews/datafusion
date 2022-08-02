library(depth)
library(copula)

# pdf("~/Publications/Books/Habilitacja/figures/tdepth.pdf", height=4, width=6)
par(mar=c(2.5, 2.5, 0.5, 0.5))
set.seed(26659)
n <- 7
copula <- frankCopula(dim=2, 1.2)
# copula <- plackettCopula(5)
X <- t(rCopula(n, copula))
X[1,] <- qexp(X[1,])
X[2,] <- qt(X[2,], 2)
# plot(X[1,], X[2,], )
isodepth(t(X), las=1, xlab=NA, ylab=NA)
d <- isodepth(t(X), las=1, xlab=NA, ylab=NA, output=TRUE)

for (i in 1:(n-1))
   for (j in (i+1):n)
      lines(X[1,c(i,j)], X[2,c(i,j)], lty=3, col='gray')

deep <- d$Contour3
deep <- rbind(deep, deep[1,])
deepn <- nrow(deep)
Cx <- sum( (deep[-1,1]+deep[-deepn,1])*(deep[-deepn,1]*deep[-1,2]-deep[-1,1]*deep[-deepn,2]) )/(3*sum(deep[-deepn,1]*deep[-1,2]-deep[-1,1]*deep[-deepn,2]))
Cy <- sum( (deep[-1,2]+deep[-deepn,2])*(deep[-deepn,1]*deep[-1,2]-deep[-1,1]*deep[-deepn,2]) )/(3*sum(deep[-deepn,1]*deep[-1,2]-deep[-1,1]*deep[-deepn,2]))
points(Cx, Cy, col=2, pch=8)
# med(t(X))
# dev.off()

