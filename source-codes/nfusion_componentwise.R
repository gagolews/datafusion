comean <- rowMeans

comed <- compiler::cmpfun(function(X)
   apply(X, 1, median)
)
