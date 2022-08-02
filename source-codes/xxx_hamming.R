###
### Closest and median string w.r.t. the Hamming distance

hamming_dist <- function(x, y)  sum(x!=y)
Rcpp::sourceCpp("hamming_closest.cpp")
Rcpp::sourceCpp("hamming_median.cpp")

#############################################################################

set.seed(36)
S <- 1:4 # alphabet, consecutive integers only
n <- 25
d <- 2
X <- matrix(nrow=d, runif(n*d))
plot(t(X))


cat("----Original---\n")
# print(hamming_dist_max(as.matrix(x), X))

#############################################################################

# x <- sample(S, d, replace=TRUE)
# X <- sapply(1:n, function(...) {
#    e <- d*0.1
#    i <- sample(d, e)
#    y <- x
#    y[i] <- sample(S, e, replace=TRUE)
#    y
# })


cat("---- Median string ----\n")

hamming_median_R <- function(X) {
   d <- nrow(X)
   res <- lapply(1:d, function(i) {
      tab <- table(X[i,])
      as.integer(names(tab)[tab == max(tab)])
   })
   unname(t(do.call(expand.grid, res))) # all solutions
   # Or sapply(res, "[", 1) for just one solution
   # sapply(res, "[", 1)
}

hamming_median_exponential <- function(X) {
   # all possible strings of length d
   X2 <- X
   X2[,] <- as.integer(factor(X2))
   S <- unique(X2)
   A <- t(do.call(expand.grid, rep(list(S), nrow(X2))))
   res_median <- hamming_dist_sum(A, X2)
   print(min(res_median))
   unname(A[, res_median == min(res_median)])
}

# M1 <- hamming_median_R(X)
# M2 <- hamming_median_exponential(X)
M3 <- hamming_median(X)
# print(t(M1))
# print(M2)
# print(M3)

print(M3)
print(hamming_dist_max(as.matrix(M3), X))
