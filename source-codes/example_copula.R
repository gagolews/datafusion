n <- 100
d <- 2
C <- list(
   copula::normalCopula(dim=d, 0.4),
   copula::indepCopula(dim=d),
   copula::claytonCopula(dim=d, param=4) # Clayton copula, theta=2
)
Finv <- list( # marginal c.d.f.s (inverses)
   function(y) qnorm(y, 0, 1),  # F_1 = N(0,1)
   function(y) qexp(y, 0.1) # F_1 = Exp(0.1)
)

for (j in 1:3) {
   set.seed(123)
   X <- t(copula::rCopula(n, C[[j]]))
   for (i in 1:d) X[i,] <- Finv[[i]](X[i,])
   pdf(sprintf("~/Publications/Books/Habilitacja/figures/copula%d.pdf", j), width=4, height=4)
   par(mar=rep(0.1,4))
   plot(X[1,], X[2,], las=1, ann=FALSE, axes=FALSE)
   box()
   dev.off()
}
