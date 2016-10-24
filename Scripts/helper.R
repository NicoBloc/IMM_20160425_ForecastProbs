readCBPs <- function() {
  return(read.table('../../Data/cbp.txt', header=TRUE))
}


readSigmaE <- function() {
  return(read.table('../../Data/errorStdev.txt', header=TRUE))
}
