require(microbenchmark)
require(grup)
require(copula)
source('devel/MG_generate_data.R')

s <- 174
N <- rep(10000, 5)
D <- rep(200, 5)
METRIC <- c("euclidean", "manhattan", "maximum", "dinu", "levenshtein")
SCENARIO <- c("rnorm", "rnorm", "rnorm",  "actg", "ispell")

stopifnot(length(N) == length(D), length(METRIC) == length(SCENARIO), length(N) == length(METRIC))

i <- 5
res <- replicate(10, {
   d <- D[i]
   n <- N[i]
   metric <- METRIC[i]
   scenario <- SCENARIO[i]

   # set.seed(s)
   x <- generateData(n, d, metric, scenario)

   system2('date')
   cat(sprintf("n = %d, d = %d, s = %d, metric = %s, scenario = %s\n", n, d, s, metric, scenario))

   r1 <- medoid_exact(metric, x)
   r2 <- medoid_approx(metric, x, iters=15, nntry=5)
   c(alg2=r1, alg2=unlist(attributes(r1)), approx=r2, approx=unlist(attributes(r2)))
})

res <- as.data.frame(t(res))
res$alg2.speedup <- res$alg2.dist_theoretical/res$alg2.dist_calls
res$approx.speedup <- res$approx.dist_theoretical/res$approx.dist_calls
res$approx.relerr <- (res$approx.penalty-res$alg2.penalty)/res$alg2.penalty
# print(res)

options(scipen=23)
print(summary(res[,-(1:8)]))
