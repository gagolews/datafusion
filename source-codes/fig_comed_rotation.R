comed <- function(X) {
   stopifnot(is.numeric(X), is.matrix(X))
   apply(X, 1, median)
}



# pdf("~/Publications/Books/habilitacja/figures/component_median_rotation1.pdf", height=4, width=4)

set.seed(12347)
n <- 9
X <- matrix(rnorm(2*n), nrow=2, byrow=TRUE)
par(mar=c(2.5, 2.5, 0.5, 0.5))
cm <- comed(X)
# X[1,] <- X[1,]-cm[1]
# X[2,] <- X[2,]-cm[2]
# cm <- comed(X)

# check_translation2(X, comed)
# check_scaling2(X, comed)

par(mar=c(2.5, 2.5, 0.5, 0.5))
plot(X[1,], X[2,], xlim=range(X[1,]), ylim=range(X[2,]), las=1, xlab=NA, ylab=NA, asp=1)
abline(h=cm[2], v=cm[1], lty=3, col='gray')
points(cm[1], cm[2], col=2, pch=8, cex=2)

# dev.off()


# pdf("~/Publications/Books/habilitacja/figures/component_median_rotation2.pdf", height=4, width=4)

THETA <- seq(0, 2*pi, length=1000)[-1]
MEDS <- matrix(nrow=2, ncol=length(THETA))
for (j in 1:length(THETA)) {
   theta <- THETA[j]
   A <- matrix(c(
      cos(theta), -sin(theta),
      sin(theta), cos(theta)
   ), ncol=2, byrow=TRUE)

   X2 <- apply(X, 2, function(row) A%*%row)
   MEDS[,j] <- comed(X2)
   MEDS[,j] <- solve(A)%*%MEDS[,j] # A orthogonal => solve(A) == t(A)
}

par(mar=c(2.5, 2.5, 0.5, 0.5))

plot(MEDS[1,], MEDS[2,], col=2, pch=16, xlim=range(X[1,]), ylim=range(X[2,]), las=1, xlab=NA, ylab=NA, cex=0.4, asp=1)
abline(h=cm[2], v=cm[1], lty=3, col='gray')

# dev.off()
