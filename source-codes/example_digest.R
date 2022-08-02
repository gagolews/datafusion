qmean <- function(x) sqrt(mean(x^2)) # quadratic mean
gmean <- function(x) exp(mean(log(x)))
hmean <- function(x) 1/mean(1/x)
sl <- function(x) min(1, sum(x))

f <- list(mean, qmean, gmean, hmean, median, min, max, prod, sl)[[1]]

dch <- function(y) {
   abs(y-0.5)
}

n <- 100
res <- replicate(100000, {
   x <- runif(n)
   y <- f(x)
   dch(y)
})

print(mean(res))



x <- "fusion functions" # input string
paste(charToRaw(x), collapse="") # hex sequence
digest::digest(x, "crc32")  # CRC-32 checksum
digest::digest(x, "md5")    # MD5 checksum
digest::digest(x, "sha256") # SHA-256

x <- "Fusion functions" # input string
paste(charToRaw(x), collapse="") # hex sequence
digest::digest(x, "crc32")  # CRC-32 checksum
digest::digest(x, "md5")    # MD5 checksum
digest::digest(x, "sha256") # SHA-256
