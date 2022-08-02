## TEST IF A MULTIDIMENSIONAL FUSION FUNCTION IS
## TRANSLATION, SCALE, AND/OR ORTHOGONAL INVARIANT


check_translation <- compiler::cmpfun(function(X, f, ..., mu=10, sigma=100) {
   stopifnot(is.numeric(X), is.matrix(X))

   delta <- rnorm(nrow(X), mu, sigma)

   val1 <- f(X, ...)+delta
   val2 <- f(X+delta, ...)
   val1-val2
})


check_scaling <- compiler::cmpfun(function(X, f, ..., smin=1, smax=100) {
   stopifnot(is.numeric(X), is.matrix(X))

   scal <- runif(1, smin, smax)

   val1 <- f(X, ...)*scal
   val2 <- f(X*scal, ...)
   val1-val2
})


check_scalingn <- compiler::cmpfun(function(X, f, ..., smin=1, smax=100) {
   stopifnot(is.numeric(X), is.matrix(X))

   scal <- runif(nrow(X), smin, smax)

   val1 <- f(X, ...)*scal
   val2 <- f(diag(scal)%*%X, ...)
   val1-val2
})


# check_rotation2 <- compiler::cmpfun(function(X, f, ...) {
#    stopifnot(is.numeric(X), is.matrix(X), nrow(X) == 2)
#
#    theta <- runif(1, 0, 2*pi)
#
#    A <- matrix(c(
#       cos(theta), -sin(theta),
#       sin(theta), cos(theta)
#    ), ncol=2, byrow=TRUE)
#
#    val1 <- f(X, ...)
#    val2 <- t(A) %*% f(apply(X, 2, function(row) A%*%row), ...)
#    val1-val2
# })


Rcpp::cppFunction('
NumericMatrix rortho(int d) {
   if (d < 1) stop("d < 1");

   NumericMatrix A(d, d);
   double theta = Rf_runif(0.0, 2.0*M_PI);
   double b = double(Rf_runif(0.0, 1.0) < 0.5)*2.0-1.0;
   A(0,0) =    cos(theta);
   A(1,0) = -b*sin(theta);
   A(0,1) =    sin(theta);
   A(1,1) =  b*cos(theta);

   NumericVector x(d);

   for (int i=3; i<=d; ++i) {
      double xnorm = 0.0;
      for (int j=0; j<i; ++j) {
         x[j] = Rf_rnorm(0.0, 1.0);
         xnorm += x[j]*x[j];
      }
      xnorm = sqrt(xnorm);
      x[0] = 1.0-x[0]/xnorm;
      double xnorm2 = x[0]*x[0];
      for (int j=1; j<i; ++j) {
         x[j] = -x[j]/xnorm;
         xnorm2 += x[j]*x[j];
      }
      xnorm2 = sqrt(xnorm2);
      for (int j=0; j<i; ++j) x[j] /= xnorm2;

      for (int k=i-1; k>0; --k)
         for (int j=i-1; j>0; --j)
            A(j,k) = A(j-1, k-1);
      for (int j=1; j<i; ++j) A(0,j) = A(j,0) = 0.0;
      A(0,0) = 1.0;

      for (int k=0; k<i; ++k) {
         double x2 = 0.0;
         for (int j=0; j<i; ++j) x2 += x[j]*A(j,k);
         for (int j=0; j<i; ++j) A(j,k) -=  2*x[j]*x2;
      }
   }

   return A;
}')

# A <- rortho(d<- 300)
# stopifnot(sum((rep(1, d)-apply(A, 2, function(c) sqrt(sum(c*c)))^2)) < 1e-12)


# rortho_old <- function(d) {
#    d <- as.integer(d)
#    stopifnot(d >= 2)
#
#    theta <- runif(1, 0, 2*pi)
#    b <- 2.0*(runif(1) < 0.5)-1.0
#
#    A <- matrix(c(
#       cos(theta), sin(theta),
#       -b*sin(theta), b*cos(theta)
#    ), ncol=2, byrow=TRUE)
#
#    # this is not very intelligent - this may be done in O(d^3)
#    i <- 2
#    while (i < d) {
#       i <- i+1
#       v <- rnorm(i)
#       v <- as.matrix(v/sqrt(sum(v*v))) # random distrib on a unit sphere
#       e1 <- c(1, rep(0, i-1))
#       x <- (e1-v)/sqrt(sum((e1-v)^2))
#       I <- diag(i)
#       A2 <- I
#       A2[2:i, 2:i] <- A
# #       A <- (I-2*x%*%t(x))%*%A2 # d^3
#       A <- A2 - 2*x%*%(t(x)%*%A2) # d^2
#    }
#
#    stopifnot(sum((rep(1, d)-apply(A, 2, function(c) sqrt(sum(c*c)))^2)) < 1e-12)
#    A
# }

check_orthogonal <- compiler::cmpfun(function(X, f, ...) {
   stopifnot(is.numeric(X), is.matrix(X))

   A <- rortho(nrow(X))

   val1 <- A%*%f(X, ...)
   val2 <- f(apply(X, 2, function(row) A%*%row), ...)
   val1-val2
})


check_affine <- compiler::cmpfun(function(X, f, ...) {
   stopifnot(is.numeric(X), is.matrix(X))

   d <- nrow(X)
   A <- switch(sample(1:3, 1),
      '1'=matrix(rnorm(d*d, 0, 2), nrow=d),
      '2'=matrix(runif(d*d, -1, 1), nrow=d),
      '3'=matrix(rexp(d*d, 1), nrow=d)
#       '4'=matrix(c(1,0,0,0), nrow=2)
   )

   val1 <- A%*%f(X, ...)
   val2 <- f(apply(X, 2, function(row) A%*%row), ...)
   val1-val2
})


check_affine2 <- compiler::cmpfun(function(X, f, ..., mu=10, sigma=100) {
   stopifnot(is.numeric(X), is.matrix(X))

   d <- nrow(X)
   A <- switch(sample(1:3, 1),
      '1'=matrix(rnorm(d*d, 0, 2), nrow=d),
      '2'=matrix(runif(d*d, -1, 1), nrow=d),
      '3'=matrix(rexp(d*d, 1), nrow=d)
#       '4'=matrix(c(1,0,0,0), nrow=2)
   )

   delta <- rnorm(nrow(X), mu, sigma)

   val1 <- A%*%f(X, ...)+delta
   val2 <- f(apply(X, 2, function(row) A%*%row+delta), ...)
   val1-val2
})

check_monotone <- compiler::cmpfun(function(X, f, ..., a=0, b=1000) {
   stopifnot(is.numeric(X), is.matrix(X))

   i <- sample(1:ncol(X), 1)

   X2 <- X
   X2[,i] <- X2[,i] + runif(nrow(X2), a, b)

   val2 <- f(X2, ...)
   val1 <- f(X, ...)

   as.numeric(all((val2-val1) >= -1e-14))-1
})

check_monotone3 <- compiler::cmpfun(function(X, f, ..., a=0, b=1000) {
   stopifnot(is.numeric(X), is.matrix(X))

   i <- sample(1:ncol(X), 1)

   X2 <- X
   X2[,i] <- X2[,i] + runif(nrow(X2), a, b)

   val2 <- f(X2, ...)
   val1 <- f(X, ...)

   # weak monotonicity by Simon James
   if (all(val1>val2))
      sum(val1-val2)
   else
      0
})


check_monotone2 <- compiler::cmpfun(function(X, f, ..., a=0, b=1000) {
   stopifnot(is.numeric(X), is.matrix(X))

   X2 <- X + runif(nrow(X), a, b)

   val2 <- f(X2, ...)
   val1 <- f(X, ...)

   as.numeric(all((val2-val1) >= 1e-12))-1

   val2 <- f(X2, ...)
   val1 <- f(X, ...)

   as.numeric(all((val2-val1) >= 0))-1
})

