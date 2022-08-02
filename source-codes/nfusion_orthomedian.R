orthomedian2 <- function(X, k=100) {
   stopifnot(is.numeric(X), is.matrix(X), nrow(X) == 2)
   n <- ncol(X)

   delta <- 0#rowMeans(X)
   s <- 1 #apply(X, 1, sd)
   X <- (X-delta)/s

   theta <- seq(0, 2*pi, length.out=k+1)[-1]
#    theta <- runif(k, 0, 2*pi)
   A <- rbind(cos(theta), sin(theta))

   projmedians <- sapply(1:k, function(i) {
      projX <- apply(X, 2, function(row) t(A[,i])%*%row)
      median(projX) * A[,i]
   })

#    rowMeans(projmedians)*2 # multiplied by d=2
  rowMeans(projmedians)*2*s+delta # multiplied by d=2
}


# orthomedian2 <- function(X, k=10000) {
#    stopifnot(is.numeric(X), is.matrix(X), nrow(X) == 2)
#    n <- ncol(X)
#
#    rowMeans(replicate(k, {
#       A <- rortho(2)
#       projX <- apply(X, 2, function(row) A%*%row)
#       comeds <- apply(projX, 1, median)
#       t(A) %*% comeds
#    }))
# }
