% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filter_correlations.R
\name{filter_correlations}
\alias{filter_correlations}
\title{filter_correlations(corrmat = NULL, mincorrelcoeff = NULL)
Given a pairwise correlation matrix, eliminate features which do not present an absolute correlation coefficient smaller than mincorrelcoeff with all other features other than itself.}
\usage{
filter_correlations(corrmat = NULL, mincorrelcoeff = NULL)
}
\description{
filter_correlations(corrmat = NULL, mincorrelcoeff = NULL)
Given a pairwise correlation matrix, eliminate features which do not present an absolute correlation coefficient smaller than mincorrelcoeff with all other features other than itself.
}
