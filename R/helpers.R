# Helper functions for phylonest


## Check hierarchical design of data
.checknested <- function(forstru) {
  n <- ncol(forstru)
  for (i in 1:(n - 1)) {
    tf <- table(forstru[, c(i, i + 1)])
    niv <- apply(tf, 1, function(x) sum(x != 0))
    if (any(niv != 1)) {
      stop(paste("non hierarchical design for structures, column", i, "is not nested in column", i + 1))
    }
  }
}



## from EqRao & EqRSintra
.diversity <- function(x) {
  if (sum(x) == 0) 
    return(0) else return(t(x) %*% d %*% x)
}

