S <- 1:4 # alphabet
n <- 100
d <- 20
l <- 10
sx <- 1234567
fname <- "levenshtein_median_bias.csv"

source('levenshtein_median.R')
set.seed(sx)
x <- sample3(S, replace=TRUE, size=d)
M <- 1000
res <- t(sapply(1:M, function(si) {
   print(si)
   set.seed(si)
   X <- levenshtein_jitter(x, n, l, S)

   cat("---- Original string ----\n")
   print(origd <- levenshtein_dist_sum(list(x), X))

   res <- levenshtein_median_ga(X)

   out <- c(n=n, d=d, l=l, S=deparse(S, width.cutoff=500),
      seedx=sx,
      seedi=si,
      gaiter=res$iter, numsol=length(res$par),
      distsolx=levenshtein_dist_sum(res$par[1], list(x)),
      solval=res$value,
      xval=origd,
      sol=deparse(res$par, width.cutoff=500),
      x=deparse(x, width.cutoff=500))

   write.table(matrix(out, nrow=1), fname, append=file.exists(fname), sep=";",
      row.names=FALSE, col.names=if (file.exists(fname)) FALSE else names(out))

   out
}))

print(res)
