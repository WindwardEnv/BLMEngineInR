# Define a matrix and an x vector
a = matrix(nrow=4, ncol=4, c(1, 3, 5, 3, 5, 1, 33, 4, 5, 1, 23, 4, 5, 6, 4, 2))
x=c(5, 3, 7, 8)

# create y
y=a %*% x

# now, we have a, x, and y such that a*x = y


# ----------------------------------------------------
# Here's the fun bit -
# pretend we don't know x, find it from a and y
# ----------------------------------------------------

# store the SVD of a
asvd = svd(a)

# store all the decomposition results
d = diag(asvd$d)
dinv = diag(1/asvd$d)
U = asvd$u
V = asvd$v

UT = t(U)
VT = t(V)

# verify that we can get matrix a back
acopy = asvd$u %*% diag(asvd$d) %*% t(asvd$v)
acopy = U %*% d %*% VT

# find the matrix inverse of a by SVD
ainv1 = V %*% dinv %*% UT

# find the matrix inverse of a by solve
# solve is based on Gaussian elimination by partial pivoting, 
# which is also what I use in CHESS
ainv2 = solve(a)

# which method is faster?
time1 = Sys.time()
for (i in 1: 10000) {
  # store the SVD of a
  asvd = svd(a)
  ainv1 = asvd$v %*% diag(1/asvd$d) %*% t(asvd$u)
}
time2 = Sys.time()

for (i in 1: 10000) {
  # inverse of a
  ainv2 = solve(a)
}
time3 = Sys.time()

time2-time1
time3-time2

# which method is faster?
# answer: solve is more than twice as fast. WHAM!  ha ha

# with either method we can solve for x
ainv1 %*% y









