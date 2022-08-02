source("nfusion_seb_Gartner_quadprog.R")

generate_seb_illustration <- function(X, y) {
   R <- apply(X, 2, function(xi) sqrt(sum((xi-y)^2)))
   r <- max(R)
   # print(r)

   theta <- seq(0, 2*pi, length=250)[-1]
   par(mar=c(2.5, 2.5, 0.5, 0.5))
   plot(X[1,], X[2,], las=1, xlim=y[1]+c(-1,1)*r, ylim=y[2]+c(-1,1)*r, asp=1, xlab=NA, ylab=NA)
   # points(mean(X[1,]), mean(X[2,]), col=2, pch=2)

   polygon(y[1]+r*cos(theta), y[2]+r*sin(theta)) # draw the ball
   points(y[1], y[2], col="red", pch=8, cex=2) # draw the SEB center

   whichbound <- which(abs(R-r)<1e-12)
   points(X[1, whichbound], X[2, whichbound], pch=16, col="#00000040")
   for (wh in whichbound)
      arrows(y[1], y[2], X[1,wh], X[2,wh], lty=2, col="gray")

   # polygon(X[1,ch], X[2,ch], lty=2)
}


#################################
set.seed(124632)
n <- 49
d <- 2
X <- matrix(rnorm(n*d), ncol=n)
print(system.time({
   y <- seb(X)
}))

# pdf("~/Publications/Books/habilitacja/figures/seb1.pdf", height=4, width=4)
generate_seb_illustration(X, y)
# dev.off()



X <- t(matrix(ncol=2, byrow=TRUE,
   c(
     0,0,
     1,0,
     0.2, 0.2
   )))

y <- seb(X)
# pdf("~/Publications/Books/habilitacja/figures/seb2.pdf", height=4, width=4)
generate_seb_illustration(X, y)
# dev.off()

# stop("!!!")

##### CHECK - general solver
set.seed(12463)
n <- 60
d <- 5
X <- matrix(rnorm(n*d), ncol=n)

replicate(100, {
   X <- matrix(rnorm(n*d), ncol=n)
   y <- seb(X)
   r1 <- max(apply(X, 2, function(xi) sqrt(sum((xi-y)^2))))

   D <- crossprod(X) # == t(X) %*% X
   c <- apply(X, 2L, # apply a function on each X's column
      function(xi) crossprod(xi)) # t(xi) %*% xi
   res <- optim(rep(1/n, n), function(y) {
      y <- abs(y)
      y <- as.matrix(y/sum(y))
      t(y) %*% D %*% y - sum(c*y)
   }, method="L-BFGS-B", lower=rep(0,n), control=list(pgtol=1e-16))
   w <- abs(res$par) / sum(abs(res$par ))
   y <- colSums(t(X)*w)
   r2 <- max(apply(X, 2, function(xi) sqrt(sum((xi-y)^2))))
   c(r1,r2)
}) -> res
print(summary(res[2,]-res[1,]))

source('nfusion_invariance_test.R')
cat(sprintf("translation %s\n", paste(fivenum(replicate(250, norm(as.matrix(check_translation(X, seb))))), collapse=', ')))
cat(sprintf("scaling     %s\n", paste(fivenum(replicate(250, norm(as.matrix(check_scaling(X, seb))))), collapse=', ')))
cat(sprintf("scalingn    %s\n", paste(fivenum(replicate(250, norm(as.matrix(check_scalingn(X, seb))))), collapse=', ')))
cat(sprintf("monotone    %s\n", paste(fivenum(replicate(250, norm(as.matrix(check_monotone(X, seb))))), collapse=', ')))
cat(sprintf("monotone2   %s\n", paste(fivenum(replicate(250, norm(as.matrix(check_monotone2(X, seb))))), collapse=', ')))
cat(sprintf("monotone3   %s\n", paste(fivenum(replicate(250, norm(as.matrix(check_monotone3(X, seb))))), collapse=', ')))
cat(sprintf("orhogonal   %s\n", paste(fivenum(replicate(250, norm(as.matrix(check_orthogonal(X, seb))))), collapse=', ')))
cat(sprintf("affine      %s\n", paste(fivenum(replicate(250, norm(as.matrix(check_affine(X, seb))))), collapse=', ')))

# A <- cbind(rep(1, n), diag(n))
# b0 <- c(1, rep(0, n))
#
# quadprog::solve.QP(D, d, A, b0, meq=1)
#
# C%*%d




#################################
## Is it an aggregation function? NO
# set.seed(126)
# n <- 3
# d <- 2
# X <- matrix(byrow=TRUE, nrow=d, c(
#    1, -1,   -sqrt(2),
#    -1, 1, 0
#    ))
# m <- seb(X)
# generate_seb_illustration(X, m)
#
# Y <- X
# Y[,1] <- Y[,1]+c(3,0)
# m2 <- seb(Y)
# print(m)
# print(m2)
