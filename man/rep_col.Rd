% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{rep_col}
\alias{rep_col}
\title{Replicate a vector column-wise to generate a matrix}
\usage{
rep_col(x, n)
}
\arguments{
\item{x}{A numeric vector.}

\item{n}{Number of times to repeat \code{x}.}
}
\value{
A numeric matrix.
}
\description{
Replicate a vector column-wise to generate a matrix
}
\examples{
x = c(1:10)
rep_col(x, 3)
}
