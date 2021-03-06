% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/add_surf.R
\name{add_surf}
\alias{add_surf}
\title{Add surface(s) to a \code{"xifti"} object}
\usage{
add_surf(xifti, surfL = NULL, surfR = NULL)
}
\arguments{
\item{xifti}{Object of class "xifti".
See \code{\link{is.xifti}} and \code{\link{make_xifti}}.}

\item{surfL}{(Optional) Left brain surface model. Can be a file
path to a GIFTI surface geometry file (ends in "*.surf.gii"), a
\code{"gifti"} object representing surface geometry, or a \code{"surf"}
object.}

\item{surfR}{(Optional) Right brain surface model. Can be a file
path to a GIFTI surface geometry file (ends in "*.surf.gii"), a
\code{"gifti"} object representing surface geometry, or a \code{"surf"}
object.}
}
\value{
the \code{"xifti"} object with added surface geometry components.
}
\description{
Add left or right cortical surface geometry to a \code{"xifti"} object.
}
\details{
\code{surfL} will be added to \code{xifti$surf$cortex_left} and \code{surfR}
will be added to \code{xifti$surf$cortex_right}. Any existing surfaces will
be overwritten.
}
