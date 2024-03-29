% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Functions.R
\name{plotTotalNASCByPSU}
\alias{plotTotalNASCByPSU}
\title{Plot summed NASC by PSU}
\usage{
plotTotalNASCByPSU(
  x,
  maxDepth = NULL,
  xlab = NULL,
  ylab = NULL,
  scaleNASC = NULL
)
}
\arguments{
\item{x}{A data.table holding at least the the columns "Beam", "PSU", "MaxLayerDepth" and "NASC".}

\item{maxDepth}{If given only the layers down to the maxDepth is included.}

\item{xlab}{The x label.}

\item{ylab}{The y label.}
}
\description{
Plot summed NASC by PSU
}
