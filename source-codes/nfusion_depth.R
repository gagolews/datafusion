suppressMessages(suppressWarnings({
   library('depth')
}))


tukeyMedian <- function(X, ...) {
   med(t(X), "Tukey", ...)$median
}

liuMedian <- function(X, ...) {
   med(t(X), "Liu", ...)$median
}

ojaMedian <- function(X, ...) {
   med(t(X), "Oja", ...)$median
}

spatialMedian <- function(X, ...) {
   med(t(X), "Spatial", ...)$median
}
