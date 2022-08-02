Rcpp::sourceCpp("~/Publications/Books/Habilitacja/comonotonic.cpp")


is_comonotonic_new <- is_comonotonic
n <- 1000
replicate(1000, {
   x <- sort(sample(1:n, n, replace=TRUE))
   y <- sort(sample(1:n, n, replace=TRUE))
   o <- sample(1:n)
   stopifnot(is_comonotonic_old(x[o], y[o]))
   stopifnot(is_comonotonic_new(x[o], y[o]))

   s <- sample(2:n, 1)
   x[o[s]] <- x[o[s-1]] - 1
   stopifnot(is_comonotonic_old(x[o], y[o]) == is_comonotonic_new(x[o], y[o]))

   s <- sample(2:n, 1)
   x[o[s]] <- x[o[s-1]] - 1
   stopifnot(is_comonotonic_old(x[o], y[o]) == is_comonotonic_new(x[o], y[o]))

   x <- sample(1:n, n, replace=TRUE)
   y <- sample(1:n, n, replace=TRUE)
   stopifnot(is_comonotonic_old(x[o], y[o]) == is_comonotonic_new(x[o], y[o]))
   stopifnot(is_comonotonic_old(y[o], x[o]) == is_comonotonic_new(y[o], x[o]))

   x <- sample(1:10, n, replace=TRUE)
   y <- sample(1:10, n, replace=TRUE)
   stopifnot(is_comonotonic_old(x[o], y[o]) == is_comonotonic_new(x[o], y[o]))
   stopifnot(is_comonotonic_old(y[o], x[o]) == is_comonotonic_new(y[o], x[o]))

   x <- sample(1:10, n, replace=TRUE)
   y <- sample(1:n, n, replace=FALSE)
   stopifnot(is_comonotonic_old(x[o], y[o]) == is_comonotonic_new(x[o], y[o]))
   stopifnot(is_comonotonic_old(y[o], x[o]) == is_comonotonic_new(y[o], x[o]))

   x <- sort(x)
   y <- sort(y)
   i <- sample(1:n, 2)
   x[i] <- x[rev(i)]
   stopifnot(is_comonotonic_old(x[o], y[o]) == is_comonotonic_new(x[o], y[o]))

   i <- sample(1:n, sample(1:n, 1))
   x[i] <- x[rev(i)]
   stopifnot(is_comonotonic_old(x[o], y[o]) == is_comonotonic_new(x[o], y[o]))
})


x <- sort(sample(1:n, n, replace=TRUE))
y <- sort(sample(1:n, n, replace=TRUE))
o <- sample(1:n)
x2 <- x[o]
y2 <- y[o]
print(microbenchmark::microbenchmark(
   is_comonotonic_old(x, y),
   is_comonotonic(x2, y2),
   is_comonotonic2(x2, y2),
   is_comonotonic3(x2, y2)
))

s <- sample(2:n, 1)
x[o[s]] <- x[o[s-1]] - 1
x2 <- x[o]
y2 <- y[o]
print(microbenchmark::microbenchmark(
   is_comonotonic_old(x, y),
   is_comonotonic(x2, y2),
   is_comonotonic2(x2, y2),
   is_comonotonic3(x2, y2)
))
