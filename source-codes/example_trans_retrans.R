set.seed(12347)
n <- 50
X <- matrix(rnorm(2*n), nrow=2, byrow=TRUE)
par(mar=c(2.5, 2.5, 0.5, 0.5))

plot(X[1,], X[2,])



# A <- rortho(2)

# X <- A%*%X

# points(Y[1,], Y[2,], col=2)


tret_med <- function(X, IDX=c(1:nrow(X))) {
   S <- rowMeans(X)
   X <- X-S
   T <- X[,IDX]
   Tinv <- solve(T)
   Y <- Tinv %*% X
   M <- apply(Y, 1, median)
   T %*% M + S
}


M <- tret_med(X)
points(M[1], M[2], col=2)
