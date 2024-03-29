% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is.xifti.R
\name{is.subcort_labs}
\alias{is.subcort_labs}
\title{Validate a factor vector of subcortical labels}
\usage{
is.subcort_labs(x)
}
\arguments{
\item{x}{The putative factor vector of brain substructures.}
}
\value{
Logical. Is \code{x} a factor vector of subcortical
structures?
}
\description{
Check if object is a factor vector of subcortical structures.This is a helper
function for \code{\link{is.xifti}}.
}
\details{
Requirements: see the "Label Levels" section for the expected factor levels.
}
\section{Label Levels}{

\code{xifti$meta$subcort$labels} is a factor with the following levels:

\enumerate{
\item{Cortex-L}
\item{Cortex-R}
\item{Accumbens-L}
\item{Accumbens-R}
\item{Amygdala-L}
\item{Amygdala-R}
\item{Brain Stem}
\item{Caudate-L}
\item{Caudate-R}
\item{Cerebellum-L}
\item{Cerebellum-R}
\item{Diencephalon-L}
\item{Diencephalon-R}
\item{Hippocampus-L}
\item{Hippocampus-R}
\item{Pallidum-L}
\item{Pallidum-R}
\item{Putamen-L}
\item{Putamen-R}
\item{Thalamus-L}
\item{Thalamus-R}
}

These correspond to the same structures as given by
\code{ft_read_cifti} in the \code{cifti-matlab} MATLAB toolbox. Note that
the first two levels (left and right cortex) are not used.
}

\keyword{internal}
