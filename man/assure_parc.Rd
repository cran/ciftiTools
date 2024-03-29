% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parc_utils.R
\name{assure_parc}
\alias{assure_parc}
\title{Assure this is a parcellation}
\usage{
assure_parc(parc)
}
\arguments{
\item{parc}{The putative parcellation.}
}
\value{
\code{parc}, if it's a parcellation.
}
\description{
Assure an input \code{"xifti"} object represents a parcellation. Keep only
the first column if multiple columns are present.
}
\keyword{internal}
