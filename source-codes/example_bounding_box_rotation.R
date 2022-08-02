rot <- function(t) {
   matrix(c(cos(t), sin(t), -sin(t), cos(t)), nrow=2)
}

BB <-function(X) {
   matrix(c(
      min(X[1,]),min(X[1,]),max(X[1,]),max(X[1,]),
      min(X[2,]),max(X[2,]),max(X[2,]),min(X[2,])
   ), byrow=TRUE, nrow=2)
}

set.seed(123)
X <- matrix(rnorm(20), nrow=2)+c(3,2)
plot(X[1,], X[2,],asp=1,col=2,pch=16)
b <- BB(X)

for (t in seq(0, 2*pi, len=100)[-1]) {
   A <- rot(t)
   b <- solve(A)%*%BB(A%*%X)
   polygon(b[1,], b[2,], border="#00000030")
}
