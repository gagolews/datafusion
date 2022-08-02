library(splines)
a <- -4
b <- -1

p <- 1 # B-spline degree
k <- 2   # number of internal knots

t <- tail(head(seq(a, b, length=k+2), -1), -1) # uniformly distributed internal knots
stopifnot(length(t) == k)
stopifnot(diff(t) > 0)
knots <- c(rep(a,p+1), t, rep(b,p+1))
stopifnot(length(knots) == 2*p+2+k)            # m+1 == 2*p+2+k
stopifnot(knots[p+1] == a, knots[p+k+2] == b)

x <- seq(a, b, length=1001)
y <- splineDesign(knots, x, ord=p+1)

# there are (n == k+p+1) basis functions

# pdf(sprintf("~/Publications/Books/Habilitacja/figures/bspline_basis_%d_%d.pdf", p, k), height=4, width=4)
par(mar=c(4, 4.2, 0.5, 0.5))
matplot(x, y, type='l', lty=1, xlab=expression(theta), ylab=expression(N["i-p,p"](theta)), las=1)
abline(v=knots, lty=3, col="gray")
# dev.off()

rowSums(y)

n <- k+p+1
stopifnot(ncol(y) == n)

# T <- function(i, t) {
#    if (i < 1) a
#    else if (i > length(t)) b
#    else t[i]
# }
#
# N <- function(i, j, theta, t) {
#    if (j == 0) {
#       as.numeric(theta >= T(i-1, t) && theta <= T(i, t))
#    }
#    else {
#       x1 <- if ((T(i+j-1,t)-T(i-1,t)) == 0) 0 else (theta-T(i-1, t))*N(i,j-1,theta,t)/(T(i+j-1,t)-T(i-1,t))
#       x2 <- if ((T(i+j,t)-T(i,t)) == 0) 0 else (T(i+j,t)-theta)*N(i+1,j-1,theta,t)/(T(i+j,t)-T(i,t))
#       x1+x2
#    }
# }
#
# y2 <- matrix(NA_real_, nrow=length(x), ncol=n)
# for (xi in seq_along(x))
#    for (i in 1:n)
#       y2[xi,i] <- N(i-p, p, x[xi], t)
# matplot(x, y2, type='l', lty=1)
# abline(v=knots, lty=3, col="gray")
#
# print(max(abs(y-y2)))

# pdf(sprintf("~/Publications/Books/Habilitacja/figures/bspline_example_%d_%d.pdf", p, k), height=4, width=4)
par(mar=c(4, 4.2, 0.5, 0.5))
w <- sort(c(a, (runif(n-2, a, b)), b))
# w <- c(0, 0.1, 0.15, 0.2, 0.95, 1)
# w <- c(0, 0.25, 0.8, 1)
stopifnot(length(w) == n)
fun <- as.numeric(t(w) %*% t(y))
plot(x, fun, type='l', xlab=expression(theta), ylab=expression(B[v](theta)), las=1)

abline(h=w, lty=3, col="gray")
abline(v=c(a, t, b), lty=3, col="gray")
# points(c(a, t, b), w)
# dev.off()



# # 2D
# w <- matrix(byrow=TRUE, nrow=2, c(
#    sort(c(a, (runif(n-2, a, b)), b)),
#    sort(c(a, (runif(n-2, a, b)), b))
# ))
#
# f <- interpSpline(w[1,], w[2,], bSpline=TRUE)
# unclass(f)
# par(mfrow=c(2,1))
# plot(f)
# points(w[1,], w[2,])
# plot(backSpline(f))
# f(5)
