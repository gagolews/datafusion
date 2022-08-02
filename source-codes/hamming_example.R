###
### Closest and median string w.r.t. the Hamming distance

hamming_dist <- function(x, y)  sum(x!=y)
Rcpp::sourceCpp("hamming_closest.cpp")
Rcpp::sourceCpp("hamming_median.cpp")

#############################################################################

set.seed(36)
S <- 1:4 # alphabet, consecutive integers only
n <- 100
d <- 100
# X <- matrix(nrow=d, sample(S, n*d, replace=TRUE))

x <- sample(S, d, replace=TRUE)
X <- sapply(1:n, function(...) {
   e <- min(d, 1+rpois(1, d*0.5))
   i <- sample(1:d, e)
   y <- x
   y[i] <- sample(S, e, replace=TRUE)
   y
})

cat("----Original---\n")
print(hamming_dist_max(as.matrix(x), X))

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
   S <- unique(as.integer(X))
   A <- t(do.call(expand.grid, rep(list(S), nrow(X))))
   res_median <- hamming_dist_sum(A, X)
   print(min(res_median))
   unname(A[, res_median == min(res_median)])
}

# M1 <- hamming_median_R(X)
# M2 <- hamming_median_exponential(X)
M3 <- hamming_median(X)
# print(t(M1))
# print(M2)
# print(M3)

print(hamming_dist_max(as.matrix(M3), X))


##############################################################################

hamming_seboid <- function(X) {
   X[,
      which.min(hamming_dist_max(X, X))
   ]
}

M4 <- hamming_seboid(X)

cat("---- seboid ----\n")
print(hamming_dist_max(as.matrix(M4), X))


##############################################################################


cat("---- Closest string ----\n")
hamming_closest_exponential <- function(X) {
   # all possible strings of length d
   S <- unique(as.integer(X))
   A <- t(do.call(expand.grid, rep(list(S), nrow(X))))
   res_closest <- hamming_dist_max(A, X)
   print(min(res_closest))
   unname(A[,res_closest == min(res_closest)])
}

# if (d <= 7) {
# x   M2 <- hamming_closest_exponential(X)
#    print(M2)
# }


# print(hamming_dist_max(M1, X))
# print(hamming_dist_max(X, X))

hamming_closest_ga <- function(X, k=length(X)*4, niter=1000, lambdaMutMult=0.001) {
   n <- ncol(X)
   d <- nrow(X)
   S <- unique(as.integer(X))
   lambdaMut <- max(1, k*d*lambdaMutMult) # expected value of number of bits to mutate per iteration

   selection <- function(P, f) {
      p <- (d-f+1)^3 # max(f) == d
      p <- p/sum(p)
      P[,sample(k, replace=TRUE, size=2*k, prob=p)]
   }

   crossover <- function(P2) {
      P <- P2[,1:k]
      for (i in 1:k) {
         # b <- sample(d, 1)
         # P[1:b,i] <- P2[1:b,i+k]

         # this works far better:
         b <- sample(d, d/2)
         P[b, i] <- P2[b, i+k]
      }
      P
   }

   mutation <- function(P) {
      m <- sample(length(P), min(k*d, rpois(1, lambdaMut)))
      P[m] <- sample(S, length(m), replace=TRUE)
      P
   }

   # initial population: points in X and random ones (mixed):
   P <- matrix(nrow=d, sample(S, k*d, replace=TRUE))
   P[,sample(k, min(n, k))] <- X[,sample(n, min(n, k))]

   # store the best solution so far:
   f <- hamming_dist_max(P, X)
   bestP <- t(unique(t(P[,f==min(f)])))
   bestF <- min(f)

   for (i in 1:niter) {
      P <- mutation(crossover(selection(P, f)))

      f <- hamming_dist_max(P, X)
      if (bestF > min(f)) { # we got a better solution
         bestP <- t(unique(t(P[,f==min(f)])))
         bestF <- min(f)
      }
      if (i %% 10 == 0) cat(sprintf("#%05d: cur Hamming=[%g, %g], best so far=%g\n", i, min(f), max(f), bestF))
   }

   print(bestF)
   bestP # return value
}

k <- 1000
print(system.time(M4 <- hamming_closest_ga(X, k=3000, niter=2000, lambdaMutMult=0.001)))
# print(M4)
