suppressMessages(suppressWarnings({
   library('copula')
   library('depth')
}))

source('nfusion_invariance_test.R')
source('nfusion_1median_Brimberg.R') # Weiszfeld1median
source('nfusion_seb_Gartner_quadprog.R') # seb
source('nfusion_orthomedian.R') # orthomedian2
source('nfusion_componentwise.R') # comean, comed
source('nfusion_depth.R') # tukeyMedian, spatialMedian, liuMedian, ojaMedian
# Weiszfeld1median == spatialMedian

# set.seed(126777)
n <- 100
d <- 2
copula <- frankCopula(dim=d, 1.2)
# copula <- plackettCopula(5)
X <- t(rCopula(n, copula))
X[1,] <- qexp(X[1,])
X[2,] <- qt(X[2,], 2)
A <- rortho(nrow(X))
X <- apply(X, 2, function(row) A%*%row)
plot(X[1,], X[2,])


# d <- apply(X, 2, function(x) depth(x, t(X), method="Oja"))
# plot(X[1,], X[2,], col=gray(1.0-(0.1+0.9*(d-min(d))/diff(range(d)))), cex=3)
# text(X[1,], X[2,], round(d,2))

# for (i in 2:100) {
# #    X2 <- X
#    X2 <- t(rortho(2))%*%X
#    X2 <- X2-rowMeans(X2)
# # plot(X2[1,], X2[2,])
#    V <- svd(t(X2))$v
#    Y <- t(V) %*% X2
#
#
# #    if (max(Y[1,]) < max(-Y[1,]))
# #       Y[1,] <- -Y[1,]
# #    if (max(Y[2,]) < max(-Y[2,]))
# #       Y[2,] <- -Y[2,]
#    if (i == 2)
#       plot(X2[1,], X2[2,], xlim=c(-15,15), ylim=c(-10,10))
#    else
#       points(X2[1,], X2[2,], col=i)
#
#    y <- as.matrix(apply(Y, 1, median))
#    print(y)
# #    y <- t(solve(V))%*%y
#    points(y[1,1], y[2,1], pch=17, col=2)
# }

funs <- expression(comean, comed, Weiszfeld1median, seb, tukeyMedian, liuMedian, ojaMedian, comedortho, orthomedian2, tret_med)
funs <- funs[c(      F,     F,            F,          F,     F,           F,         F,           F,          F,       T)]
for(i in seq_along(funs))  {
   cat("\n\n", as.character(funs[i]), "\n")
   fe <- get(as.character(funs[i]))
   v <- fe(X)
   print(v)
   points(v[1], v[2], pch=1+i, col=1+i)
#    text(v[1], v[2], labels=as.character(funs[i]), col=1+i)

   cat(sprintf("translation  %s\n", paste(sprintf("%.13f", fivenum(replicate(250, norm(as.matrix(check_translation(X, fe)))))), collapse=', ')))
   cat(sprintf("scaling1     %s\n", paste(sprintf("%.13f", fivenum(replicate(250, norm(as.matrix(check_scaling(X, fe, smin=1,smax=100)))))), collapse=', ')))
   cat(sprintf("scalingn     %s\n", paste(sprintf("%.13f", fivenum(replicate(250, norm(as.matrix(check_scalingn(X, fe, smin=1,smax=100)))))), collapse=', ')))
   cat(sprintf("orhogonal    %s\n", paste(sprintf("%.13f", fivenum(replicate(250, norm(as.matrix(check_orthogonal(X, fe)))))), collapse=', ')))
   cat(sprintf("affine       %s\n", paste(sprintf("%.13f", fivenum(replicate(250, norm(as.matrix(check_affine(X, fe)))))), collapse=', ')))
   cat(sprintf("monotone     %s\n", paste(sprintf("%.13f", apply(replicate(2500, check_monotone(X, fe)), 1, min)), collapse=', ')))
}

