Rcpp::sourceCpp('dinudist.cpp')

print(dinudist(c(1,2,3,1,1), c(2,1,1,3,3)))
print(dinudist(c(1,2,2,1), c(1,1,2,2)))

print(dinudist(c(1,2,2,1), c(2,1,2,1)))
print(dinudist(c(2,1,2,1), c(2,2,1,1)))
print(dinudist(c(1,2,2,1), c(1,1,2,2,1)))


print(dinudist(c(1,2,3,4), c(1,2,4,3)))
print(dinudist(c(1,2,3,4), c(1,4,3,2)))
print(dinudist(c(1,2,3,4), c(1,2,3,5,4)))
