rzad <- function(n, phi, mu, patterns) {

  D <- dim(patterns)[2] - 1
  p <- patterns[, D + 1]
  patterns <- patterns[, 1:D]
  k <- dim(patterns)[1]
  nk <- round( n * p )

  x0 <- NULL
  for ( i in 1:k ) {
    x0 <- rbind(x0, Compositional::rdiri(nk[i], phi * mu * patterns[i, ]) )
  }
  n1 <- n - sum(nk)
  x1 <- Compositional::rdiri(n1, phi * mu)
  rbind(x1, x0)
}
