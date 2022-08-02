library(circular)

set.seed(1234578)
# x <- circular(rnorm(25, 180, 30)*pi/180)
# rose.diag(x)
# points(x

cairo_pdf("~/Publications/Books/Habilitacja/figures/rose.pdf", height=4, width=6)
par(mar=rep(2,4))
par(xpd=TRUE)
x <- rvonmises(n=50, mu=circular(pi/4), kappa=5)
y <- rose.diag(x, bins=18, shrink=0.8, prop=1.5) # Points fall out of bounds.
points(x)
# points(x, plot.info=y, stack=TRUE, pch=1)
dev.off()
