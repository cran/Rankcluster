#' @title Kendall distance between two ranks
#'
#' @description The Kendall distance between two ranks is the number of pairs that are in different order in the two ranks.
#'
#' @param x,y two ranks of size m.
#' @param type type of the rank representation ("ordering" ou "ranking").
#'
#' @return an integer, the Kendall distance between x and y.
#'
#' @references A New Measure of Rank Correlation, M. G. Kendall
#'
#' @examples
#' x <- 1:5
#' y <- c(2, 3, 1, 4, 5)
#' distKendall(x, y, type = "ordering")
#'
#' @author Julien Jacques
#'
#' @family distance
#'
#' @export
distKendall <- function(x, y, type = "ordering") {
  if (type == "ordering") {
    dist <- distKendall_ordering(x, y)
  } else {
    dist <- distKendall_ranking(x, y)
  }

  return(dist)
}

distKendall_ranking <- function(x, y) {
  m <- length(x)
  distKendall_ranking <- 0
  for (i in 1:(m - 1)) {
    for (j in (i + 1):m) {
      if (((x[i] - x[j]) * (y[i] - y[j])) < 0) {
        distKendall_ranking <- distKendall_ranking + 1
      }
    }
  }
  return(distKendall_ranking)
}



distKendall_ordering <- function(x, y) {
  m <- length(x)
  distKendall_ordering <- 0
  x_ordering <- x
  y_ordering <- y
  for (i in 1:m) {
    x_ordering[i] <- which(x == i)
    y_ordering[i] <- which(y == i)
  }
  distKendall_ordering <- distKendall_ranking(x_ordering, y_ordering)
  return(distKendall_ordering)
}

#' @title Spearman distance between two ranks
#'
#' @description The Spearman distance is the square of Euclidean distance between two rank vector.
#'
#' @param x,y two ranks of size m.
#'
#' @return an integer, the Spearman distance between x and y.
#'
#' @examples
#' x <- 1:5
#' y <- c(2, 3, 1, 4, 5)
#' distSpearman(x, y)
#'
#' @author Julien Jacques
#'
#' @family distance
#'
#' @export

distSpearman <- function(x, y) {
  distSpearman <- sum((x - y)^2)
  return(distSpearman)
}

# It calculates the Spearman rank correlation between two ranks.
# @title Spearman rank correlation coefficient
# @author Julien Jacques
# @param x,y two ranks of size m.
# @return a real, the Spearman rank correlation between x and y.
# @examples
# x=1:5
# y=c(2,3,1,4,5)
# CorrelSpearman(x, y)
# @export
correlSpearman <- function(x, y) {
  correlSpearman <- sqrt(sum((x - y)^2))
  return(correlSpearman)
}

#' @title Cayley distance between two ranks
#'
#' @description The Cayley distance between two ranks x and y is the minimum number of transpositions required to
#' transform the ranking x into y.
#'
#'
#' @param x,y two ranks of size m.
#'
#' @return the Cayley distance between x and y.
#'
#' @examples
#' x <- 1:5
#' y <- c(2, 3, 1, 4, 5)
#' distCayley(x, y)
#'
#' @author Julien Jacques
#'
#' @family distance
#'
#' @export
distCayley <- function(x, y) {
  m <- length(x)
  distCayley <- 0
  for (i in 1:(m - 1)) {
    if (!x[i] == y[i]) {
      distCayley <- distCayley + 1
      tmp <- x[i]
      x[i] <- y[i]
      x[which(y[i] == x)] <- tmp
    }
  }
  return(distCayley)
}

#' @title Hamming distance between two ranks
#'
#' @description The Hamming distance between two ranks x and y is the number of difference between the two ranks.
#' For example, the Hamming's distance between x=(1,4,2,5,3) and y=(1,3,4,5,2) is 3 because, only 1 and 5 have the same
#' place in both ranks.
#'
#' @param x,y two ranks of size m.
#'
#' @return an integer, the Hamming distance between x and y.
#'
#' @examples
#' x <- 1:5
#' y <- c(2, 3, 1, 4, 5)
#' distHamming(x, y)
#'
#' @author Julien Jacques
#'
#' @family distance
#'
#' @export
distHamming <- function(x, y) {
  distHamming <- sum(x != y)

  return(distHamming)
}
