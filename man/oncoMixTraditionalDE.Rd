% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/oncoMix_examples.R
\name{oncoMixTraditionalDE}
\alias{oncoMixTraditionalDE}
\title{Creating a schematic of a traditional diff. expression experiment}
\usage{
oncoMixTraditionalDE(means = c(3, 7))
}
\arguments{
\item{means}{Set the values for the difference between parameter means}
}
\value{
Returns a ggplot object that shows the traditional method (2 sample
t-test) for mRNA differential expression.
}
\description{
This function allows you to generate a schematic of the assumptions of a
traditional DE expermiment between two known groups.
}
\examples{
oncoMixTraditionalDE(means=c(3,7))
oncoMixTraditionalDE(means=c(3,10))
}
\keyword{differential}
\keyword{expression}
\keyword{idealized,}
\keyword{oncoMix,}
\keyword{theoretical,}
