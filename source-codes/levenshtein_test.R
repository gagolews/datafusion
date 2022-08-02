library("stringdist")
library("stringi")


Rcpp::sourceCpp('levenshtein.cpp')


# x <- stri_rand_strings(2, length = c(300, 1000))
# xc <- stri_enc_toutf32(x)
# microbenchmark::microbenchmark(
#    adist(x[1], x[2]),
#    Levenshtein_smallmem(xc[[1]], xc[[2]]),
#    Levenshtein_bigmem(xc[[1]], xc[[2]])
# )


s <- c("kitten", "sitting", "a", "growing", "kite", "konstantynopolitanczykowna")
sc <- stri_enc_toutf32(s)
(a <- adist(s))
(l1 <- outer(sc, sc, Vectorize(levenshtein_bigmem, SIMPLIFY = TRUE)))
(l2 <- outer(sc, sc, Vectorize(levenshtein_smallmem, SIMPLIFY = TRUE)))


sc <- lapply(1:100, function(...) sample(1:10, rpois(1,100), replace=TRUE))
range(outer(sc, sc, Vectorize(function(x, y) {
   m <- levenshtein_mean2(x,y)
   d <- levenshtein_smallmem(x, y)
   abs(levenshtein_smallmem(x, m)+levenshtein_smallmem(m,y)-d)
}, SIMPLIFY=TRUE)))


print(a-l1)
print(a-l2)

x <- list(1:1000)
y <- list(50:1030)
microbenchmark::microbenchmark(
   levenshtein_smallmem(x[[1]], y[[1]]),
   levenshtein_bigmem(x[[1]], y[[1]])
)


levenshtein_dist_max_R <- function(P, X) {
   sapply(P, function(p) { # for each p in P
      max(sapply(X, function(x, p) levenshtein_smallmem(x, p), p))
   })
}

levenshtein_dist_sum_R <- function(P, X) {
   sapply(P, function(p) { # for each p in P
      sum(sapply(X, function(x, p) levenshtein_smallmem(x, p), p))
   })
}

t <- c("sit", "grow", "konstantynopol", "sewastopol", "mike")
tc <- stri_enc_toutf32(s)
levenshtein_dist_max(tc, sc)
levenshtein_dist_max_R(tc, sc)
levenshtein_dist_sum(tc, sc)
levenshtein_dist_sum_R(tc, sc)


sc <- lapply(1:100, function(i) as.integer(rpois(rpois(1, 20), 20)))
st <- lapply(1:1000, function(i) as.integer(rpois(rpois(1, 20), 20)))
levenshtein_dist_sum(sc, st)-levenshtein_dist_sum_R(sc,st)
