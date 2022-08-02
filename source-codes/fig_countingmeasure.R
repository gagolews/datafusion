library(agop)
x <- c(1, 6, 4.5, 1, 0, 5, 2, 4)

lastfalse <- function(x) { x[length(x)] <- FALSE; x}
t <- seq(0, 7, by=0.5)
St <- sapply(t, function(t) length(which(x >= t)))
pdf('figures/countingmeasure.pdf', width=6, height=6/sqrt(2))
par(mar=c(4,4,2,2))
plot(t[lastfalse(!duplicated(St, fromLast=TRUE))],St[lastfalse(!duplicated(St, fromLast=TRUE))], pch=16,
   xlim=c(0,max(t)), ylim=c(0,max(St)), xlab=expression(t), ylab=expression(S(t)), las=1)
lines(t[], St[], type='S')
dev.off()
