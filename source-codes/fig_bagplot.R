#install.packages("aplpack")
library("aplpack")

# pdf("~/Publications/Books/habilitacja/figures/bagplot.pdf", height=4, width=6)
par(mar=c(2.5, 2.5, 0.5, 0.5))
set.seed(12356)
x <- rnorm(100)
y <- rnorm(100)
bagplot(x, y, las=1, xlim=range(x), transparency=TRUE, col.baghull="#dddddd", col.loophull="#f5f5f5", cex=1)
box()
# dev.off()

# pdf("~/Publications/Books/habilitacja/figures/boxplot.pdf", height=1, width=6)
par(mar=c(0, 2.5, 0, 0.5))
boxplot(x, horizontal = TRUE, xlim=range(x), axes=FALSE, width=1.5)
# dev.off()



# library(depth)
# x <- c(0,0,1,1,0.5)
# y <- c(0,1,0,1,0.5)
# b <- bagplot(x, y, precision=5)
# med(cbind(x, y), "Tukey", factor=0)
# print(b$center)
# points(x, y)
#
# x <- c(0,1,1,1,0.5)
# y <- c(0,0,1,1,0.5)
# b <- bagplot(x, y, precision=5)
# points(x,y)
# med(cbind(x, y), "Spatial")
# print(b$center)
