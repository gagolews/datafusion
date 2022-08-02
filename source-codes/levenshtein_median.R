Rcpp::sourceCpp('levenshtein.cpp')



# a version of sample() which does not generate a permutation
# of 1:x if length(x) == 1
sample3 <- compiler::cmpfun(function(x, size, replace = FALSE, prob = NULL) {
   if (missing(size))
      size <- length(x)
   # print(c(x, NA, length(x), size, replace))
   x[sample.int(length(x), size, replace, prob)]
})



levenshtein_jitter <- compiler::cmpfun(function(x, n, l, S) {
   lapply(1:n, function(...) {
      y <- x
      for (i in 1:l) {
         t <- sample3(1:3, 1)
         if (t == 1) {
            # deletion
            j <- sample3(0:(length(y)), 1)
            y <- c(head(y, max(0,j)), tail(y, max(0,length(y)-j-1)))
         }
         else if (t == 2) {
            # insertion
            j <- sample3(0:(length(y)), 1)
            y <- c(head(y, max(0,j)), sample3(S, 1), tail(y, max(0,length(y)-j)))
         }
         else {
            # substitution
            j <- sample3(seq_along(y), 1)
            y[j] <- sample3(setdiff(S, y[j]), 1)
         }
      }
      y
   })
})


levenshtein_medoid <- function(X) {
   X[[
      which.min(levenshtein_dist_sum(X, X))
   ]]
}




levenshtein_median_ga <- compiler::cmpfun(function(X, k=length(X)*8, niter=2500, lambdaMutMult=0.001, verbose=TRUE) {
   n <- length(X)
   dmin <- min(sapply(X, length))
   dmax <- max(sapply(X, length))
   S <- unique(as.integer(unlist(X)))
   lambdaMut <- max(1, k*mean(c(dmin, dmax))*lambdaMutMult) # expected value of number of bits to mutate per iteration

   selection <- function(P, f, e) {
      p <- (dmax-f/n+1)^(e) # max(f/n) == dmax
      # print(range(sapply(P, length)))
      # print(dmax)
      # print(range(dmax - f/n + 1))
      p <- p/sum(p)
      # print(head(p))
      # print(head(P))
      P[sample(k, replace=TRUE, size=2*k, prob=p)]
   }

   crossover <- function(P2) {
      P <- P2[1:k]
      for (i in 1:k) {
         stopifnot(is.finite(P2[[i]]))
         stopifnot(is.finite(P2[[i+k]]))

# CUT_AND_SPLICE:
#          u1 <- sample3(0:length(P2[[i]]), 1)
#          u2 <- max(0, dmin-u1):min(length(P2[[i+k]]), dmax-u1)
#          # print(dmax-u1)
#          u2 <- if (length(u2) > 1) sample3(u2, 1) else u2
#          # stopifnot(u1+u2 >= dmin, u1+u2 <= dmax)
#          P[[i]] <- c(head(P2[[i]], u1), tail(P2[[i+k]], u2))

# \cite{DinuIonescu2012:efficientrankcloseststring}
#          u1 <- sample3(0:length(P2[[i]]), 1)
#          u2 <- max(0, dmin-u1):min(length(P2[[i+k]]), dmax-u1)
#          # print(dmax-u1)
#          u2 <- if (length(u2) > 1) sample3(u2, 1) else u2
#          # stopifnot(u1+u2 >= dmin, u1+u2 <= dmax)
#          if (runif(1) > 0.5)
#             P[[i]] <- c(head(P2[[i]], u1), sample3(tail(P2[[i+k]], u2)))
#          else
#             P[[i]] <- c(sample3(head(P2[[i]], u1)), tail(P2[[i+k]], u2))

# EXPERIMENT:
         # dput(P2[[i]])
         # dput(P2[[i+k]])
         P[[i]] <- levenshtein_centroid2(P2[[i]], P2[[i+k]])
         # stopifnot(length(P[[i]]) >= dmin, length(P[[i]]) <= dmax)

         stopifnot(is.finite(P2[[i]]))
      }
      P
   }

   mutation <- function(P) {
      # stopifnot(min(sapply(P, length)) >= dmin)
      # stopifnot(max(sapply(P, length)) <= dmax)

      # TO DO: prob. mutacji zaleze od temperatury?

      m <- max(1, rpois(1, lambdaMut)) # number of mutate ops
      U <- sample3(seq_along(P), m, replace=TRUE) # which vector?
      T <- sample3(1:3, m, replace=TRUE) # which op to apply?
      for (i in 1:m) {
         u <- U[i]
         t <- T[i]
         if (t == 1 && length(P[[u]]) > dmin) {
            # deletion
            j <- sample3(0:(length(P[[u]])), 1)
            P[[u]] <- c(head(P[[u]], max(0,j)), tail(P[[u]], max(0,length(P[[u]])-j-1)))
         }
         else if (t == 2 && length(P[[u]]) < dmax) {
            # insertion
            j <- sample3(0:(length(P[[u]])), 1)
            P[[u]] <- c(head(P[[u]], max(0,j)), sample3(S, 1), tail(P[[u]], max(0,length(P[[u]])-j)))
         }
         else {
            # substitution
            j <- sample3(seq_along(P[[u]]), 1)
            P[[u]][j] <- sample3(setdiff(S, P[[u]][j]), 1)
         }
      }
      P
   }

   # initial population: points in X and random ones (mixed):
   P <- lapply(1:k, function(...)
      sample3(S, replace=TRUE, size=sample3(dmin:dmax, 1)) # random vector
   )
   P[1:min(k, length(X))] <- X[1:min(k, length(X))]

   # store the best solution so far:
   f <- levenshtein_dist_sum(P, X)
   bestP <- unique(P[f==min(f)])
   bestF <- min(f)
   bestI <- 0

   for (i in 1:niter) {
      e <- 5-4*((niter-i)/niter)^2
      P <- mutation(crossover(selection(P, f, e)))

      f <- levenshtein_dist_sum(P, X)
      if (bestF > min(f)) { # we got a better solution
         bestP <- unique(P[f==min(f)])
         bestF <- min(f)
         bestI <- i
         # if (verbose) cat(sprintf("#%05d: cur Lev=[%g, %g], best so far=%g; e=%f\n", i, min(f), max(f), bestF, e))
      }
      else if (bestF == min(F)) {
         bestP <- unique(c(bestP, P[f==min(f)]))
         # if (verbose) cat(sprintf("#%05d: cur Lev=[%g, %g], best so far=%g; e=%f\n", i, min(f), max(f), bestF, e))
      }
      else if (verbose && i %% 100 == 1)
         cat(sprintf("#%05d: cur Lev=[%g, %g], best so far=%g; e=%f\n", i, min(f), max(f), bestF, e))
      # stopifnot(min(sapply(P, length)) >= dmin)
      # stopifnot(max(sapply(P, length)) <= dmax)

   }

   list(
      par=bestP, # return value
      value=bestF,
      iter=bestI
   )
})
