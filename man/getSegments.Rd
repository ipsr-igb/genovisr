% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getSegments.R
\name{getSegments}
\alias{getSegments}
\title{Get Segments Data}
\usage{
getSegments(object, data = "haplotype", sample = NULL)
}
\arguments{
\item{object}{A \code{genovis} class object.}

\item{data}{A character string specifying which data to use, either "haplotype" or "dosage". Default is "haplotype".}

\item{sample}{A character or numeric vector specifying samples to include in the retrieval. Default is \code{NULL}.}
}
\value{
A dataframe containing the segments data.
}
\description{
This function retrieves segments data for either haplotype or dosage from a \code{genovis} object.
The data can then be used for plotting graphical genotypes.
}
