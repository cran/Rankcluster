% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{APA}
\alias{APA}
\title{Rank data: APA}
\format{
A list containing:
\describe{
  \item{data}{A matrix of size 5738x5 containing the 5738 observed full ranks in ranking representation.
 The ranking representation r=(r_1,...,r_m) contains the  ranks assigned to the objects, and means that the ith object is
in r_ith position.

For example, if the ranking representation of a rank is (4,3,1,2,5), it means that judge ranks
the first object in 4th position, second object in 3rd position, ...}

  \item{frequency}{matrix of size 120x6. Each row corresponds to one of the different observed rank.
                  The first fifth columns contains the observed ranks (ordering representation) and the sixth column
                   contains the frequency of observation.}

  \item{m}{vector with the size of the ranks (5 here).}
}
}
\source{
"Group representations in probability and statistics", P. Diaconis, 1988.
}
\description{
This dataset contains the 5738 full rankings resulting from the American Psychological Association (APA)
presidential election of 1980. For this election, members of APA had to rank five candidates in order of preference.

For information, a total of 15449 votes have been registered for this election, but only the 5738 full rankings are
reported in the APA dataset. Candidates A and C were research psychologists, candidates D and E were clinical
psychologists and candidate B was a community psychologist.
}
\examples{
data(APA)

}
\seealso{
Other datasets: 
\code{\link{big4}},
\code{\link{eurovision}},
\code{\link{quiz}},
\code{\link{sports}},
\code{\link{words}}
}
\concept{datasets}
\keyword{datasets}
