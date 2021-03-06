% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spikein_median_normalization.R
\name{spikein_median_normalization}
\alias{spikein_median_normalization}
\title{Spike-in normalization function}
\usage{
spikein_median_normalization(gr, spikein_factor)
}
\arguments{
\item{gr}{Input signal track data as a \code{GRanges} object (see
\code{?"GRanges-class"} for more details). To load bedGraph data
run \code{\link{import_bedGraph}}. No default.}
}
\value{

}
\description{
Computes median signal (in the \code{spikescore} \code{GRanges} metadata column)
genome-wide.
}
\examples{
\dontrun{
spikein_median_normalization(GRanges_object)
}
}
