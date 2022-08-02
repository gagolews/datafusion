set.seed(1234)
library('copula')
set.seed(126)
n <- 25
d <- 2
copula <- frankCopula(dim=d, 1.2)
X <- t(rCopula(n, copula))
X[1,] <- qexp(X[1,])
X[2,] <- qt(X[2,], 3)
# X <- matrix(rcauchy(n*d), nrow=d)

source('nfusion_invariance_test.R')

f <- function(x) {
#    if (all(diff(x) == 0)) x[1]
#    else 0
   sum(-3*x^3)
}

cofusion <- function(X, f) {
   apply(X, 1, f)
}



cat(sprintf("translation %s\n", paste(fivenum(replicate(250, norm(as.matrix(check_translation(X, cofusion, f))))), collapse=', ')))
cat(sprintf("scaling     %s\n", paste(fivenum(replicate(250, norm(as.matrix(check_scaling(X, cofusion, f))))), collapse=', ')))
cat(sprintf("scalingn    %s\n", paste(fivenum(replicate(250, norm(as.matrix(check_scalingn(X, cofusion, f))))), collapse=', ')))
cat(sprintf("monotone    %s\n", paste(fivenum(replicate(250, norm(as.matrix(check_monotone(X, cofusion, f))))), collapse=', ')))
cat(sprintf("monotone2   %s\n", paste(fivenum(replicate(250, norm(as.matrix(check_monotone2(X, cofusion, f))))), collapse=', ')))
cat(sprintf("monotone3   %s\n", paste(fivenum(replicate(250, norm(as.matrix(check_monotone3(X, cofusion, f))))), collapse=', ')))
cat(sprintf("orhogonal   %s\n", paste(fivenum(replicate(250, norm(as.matrix(check_orthogonal(X, cofusion, f))))), collapse=', ')))
cat(sprintf("affine      %s\n", paste(fivenum(replicate(250, norm(as.matrix(check_affine(X, cofusion, f))))), collapse=', ')))
