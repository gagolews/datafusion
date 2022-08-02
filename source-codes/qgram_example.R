###
### median string w.r.t. the qgram distance

Y <- c(
   "abcb",
   "cbac",
   "acab"
   # "abcd",
)

# Y <- c("ab", "ba", "ba", "aba")

qgrams <- function(Y, q) {
   lapply(Y, function(x) {
      stringi::stri_sub(x, 1:(stringi::stri_length(x)-q+1), length=q)
   })
}

g <- qgrams(Y, 2)

ug <- unique(unlist(g))
X <- simplify2array(lapply(g, function(gi)
   table(factor(gi, levels=ug))
))

# source("nfusion_1median_Brimberg.R")
# Weiszfeld1median(g2) # zle -- trzeba robic search w integer domain

d <- nrow(X)
n <- ncol(X)
# (all variables are non-negative by default)
C <- c(rep(0, d), rep(1, 2*n*d))
lp <- lpSolveAPI::make.lp(ncol=d+2*d*n)


Ap <- diag(d); for (i in 1:(n-1)) Ap <- rbind(Ap, diag(d))
A <- cbind(Ap, diag(n*d), -diag(n*d))
B <- as.integer(X)
for (j in 1:nrow(A))
   lpSolveAPI::add.constraint(lp, A[j,], "=", B[j])

set.type(lp, 1:(d+2*d*n), "integer")
lpSolveAPI::set.objfn(lp, C)
testthat::expect_equal(0, lpSolveAPI::solve.lpExtPtr(lp))
w <- lpSolveAPI::get.variables(lp)[1:d]
print(structure(w, names=dimnames(X)[[1]]))
testthat::expect_equal(max(
   apply(matrix(lpSolveAPI::get.variables(lp)[-(1:d)], nrow=2, byrow=TRUE), 2, min)
), 0)
