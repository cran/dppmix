% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distribs.R
\name{dgammapois}
\alias{dgammapois}
\title{Density function for Gamma-Poisson distribution.}
\usage{
dgammapois(x, a, b = 1, log = FALSE)
}
\arguments{
\item{x}{vector of x values}

\item{a}{shape parameter for gamma distribution on mean parameter}

\item{b}{rate parameter for gamma distribution on mean parameter}

\item{log}{whether to return the density in log scale}
}
\value{
density values
}
\description{
Data follow the Poisson distribution parameterized by a mean parameter
that follows a gamma distribution.
}
