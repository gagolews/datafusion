library(RColorBrewer)
set.seed(12373)
k <- 5
C <- matrix(rnorm(2*k), ncol=2)

xmin <- min(C[,1])-diff(range(C[,1]))*0.1
xmax <- max(C[,1])+diff(range(C[,1]))*0.1
ymin <- min(C[,2])-diff(range(C[,2]))*0.1
ymax <- max(C[,2])+diff(range(C[,2]))*0.1

# png("figures/voronoi_LInf.png", height=200, width=200)

par(mar=rep(0,4))
plot(C[,1], C[,2], asp=1, xlim=c(xmin, xmax), ylim=c(ymin, ymax), xlab=NA, ylab=NA, las=1, axes=FALSE)

d <- function(x, y) max(abs(x-y))
# d <- function(x, y) sqrt(sum(abs(x-y)^2))
# d <- function(x, y) sum(abs(x-y))
d <- function(x, y) (sum(abs(x-y)^4))^(1/4)
# d <- function(x, y) sum(abs(x-y)/abs(x+y)) # canberra - no sense
# d <- function(x, y) 1-acos(sum(x*y)/sqrt(sum(x^2))/sqrt(sum(y^2)))/pi # angular similarity -- numeric problems....

palette(RColorBrewer::brewer.pal(k+1, "Accent"))

x <- seq(xmin, xmax, length.out=200)
y <- seq(ymin, ymax, length.out=200)
z <- matrix(NA_integer_, nrow=length(y), ncol=length(x))
for (i in 1:(length(y)-1))
   for (j in 1:(length(x)-1)) {
      z[i,j] <- which.min(apply(C, 1, function(Ck) d(Ck, c(x[j], y[i]))))
      rect(x[j], y[i], x[j+1], y[i+1], col=z[i,j]+1, border=z[i,j]+1)
   }

points(C[,1], C[,2], asp=1, pch=16)

# dev.off()
