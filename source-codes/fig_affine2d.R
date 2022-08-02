#### DRAW EXEMPLES OF AFFINE TRANSFORMATIONS


pdf("~/Publications/Books/habilitacja/figures/affine2d_translation.pdf", height=4, width=4)
A <- matrix(c(
   1, 0,
   0, 1
), ncol=2, byrow=TRUE)

b <- c(0.25, 0.1)


# pdf("~/Publications/Books/habilitacja/figures/affine2d_rotation.pdf", height=4, width=4)
# theta <- pi/4
# A <- matrix(c(
#    cos(theta), -sin(theta),
#    sin(theta), cos(theta)
# ), ncol=2, byrow=TRUE)
#
# b <- c(-0.25, -0.25)


# pdf("~/Publications/Books/habilitacja/figures/affine2d_shearh.pdf", height=4, width=4)
# mx <- 0.5
# A <- matrix(c(
#    1, mx,
#    0, 1
# ), ncol=2, byrow=TRUE)
#
# b <- c(0, 0)


# pdf("~/Publications/Books/habilitacja/figures/affine2d_projectionx.pdf", height=4, width=4)
# A <- matrix(c(
#    1, 0,
#    0, 0
# ), ncol=2, byrow=TRUE)
#
# b <- c(0, 0)

# pdf("~/Publications/Books/habilitacja/figures/affine2d_reflection.pdf", height=4, width=4)
# A <- matrix(c(
#    -1, 0,
#    0, 1
# ), ncol=2, byrow=TRUE)
#
# b <- c(1.5, 0.25)
#
p <- matrix(c(
   0.25,0.25,
   1,0.25,
   0.75,1
), byrow=TRUE, ncol=2)


p2 <- p
for (i in 1:nrow(p))
   p2[i,] <- A%*%p[i,]+as.matrix(b)



par(mar=c(2.5, 2.5, 0.5, 0.5))
plot(p[,1], p[,2], pch=16, las=1, xlim=c(0, 1.25), ylim=c(0,1.25), xlab=NA, ylab=NA)
polygon(p[,1], p[,2], col=rgb(0.5, 0.5, 0.5, 0.05))
points(p2[,1], p2[,2], pch=1, las=1)
polygon(p2[,1], p2[,2], col=rgb(0, 0, 0, 0.05), lty=2)

arrows(p[,1], p[,2], p2[,1], p2[,2], col="gray", length=0.15, angle=15)
dev.off()
