test_plotClusters <- function() {
  checkEquals(divideBy(4, 2), 2)
  checkTrue(is.na(divideBy(4, 0)))
  checkEqualsNumeric(divideBy(4, 1.2345), 3.24, tolerance=1.0e-4)
}