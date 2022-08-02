library("FuzzyNumbers")
A <- TrapezoidalFuzzyNumber(1, 2, 3, 4)
B <- TriangularFuzzyNumber(3, 5.5, 6)
C <- as.PiecewiseLinearFuzzyNumber(A, knot.n=100) * as.PiecewiseLinearFuzzyNumber(B, knot.n=100)
alphacut(C, c(0, 1)) # support and core

core(C)
plot(A, xlim=c(0, 20))
plot(B, add=TRUE)
plot(C, add=TRUE)
