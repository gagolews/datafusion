source('levenshtein_median.R')

set.seed(2634)
S <- 1:4 # alphabet
n <- 100
d <- 25
l <- 15
x <- sample3(S, replace=TRUE, size=d)
X <- levenshtein_jitter(x, n, l, S)

cat("---- Original string ----\n")
print(x)
print(levenshtein_dist_sum(list(x), X))


##############################################################################


M <- levenshtein_medoid(X)

cat("---- Medoid ----\n")
print(M)
print(levenshtein_dist_sum(list(M), X))


##############################################################################


print(system.time(resM4 <- levenshtein_median_ga(X)))
cat("genetic median(s):")
M4 <- resM4$par
print(M4)
print(levenshtein_dist_sum(M4, X))

cat("x-based:")
print(x)
print(levenshtein_dist_sum(list(x), X))

