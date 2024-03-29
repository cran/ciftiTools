% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/write_gifti.R
\name{write_metric_gifti}
\alias{write_metric_gifti}
\title{Write a data matrix to a GIFTI metric file}
\usage{
write_metric_gifti(
  x,
  gifti_fname,
  hemisphere = c("left", "right"),
  intent = NULL,
  data_type = NULL,
  encoding = NULL,
  endian = c("LittleEndian", "BigEndian"),
  col_names = NULL,
  label_table = NULL
)
}
\arguments{
\item{x}{A \eqn{V \times T} data matrix (V vertices, T measurements). This can also
be an object from \code{gifti::readgii}, or a length \eqn{T} list of length
\eqn{V} vectors.}

\item{gifti_fname}{Where to write the GIFTI file.}

\item{hemisphere}{\code{"left"} (default) or \code{"right"}. Ignored if
\code{data} is already a \code{"gifti"} object.}

\item{intent}{"NIFTI_INTENT_*". \code{NULL} (default) will use
metadata if \code{data} is a \code{"gifti"} object, or "NONE" if it cannot be
inferred. If not \code{NULL} and \code{data} is a \code{"gifti"} object, it will
overwrite the existing intent. See
https://nifti.nimh.nih.gov/nifti-1/documentation/nifti1fields/nifti1fields_pages/group__NIFTI1__INTENT__CODES.html/document_view .}

\item{data_type}{the type of \code{data}:
"NIFTI_TYPE_*" where * is "INT32" or "FLOAT32". If \code{NULL} (default), the
data type will be inferred. If not \code{NULL} and \code{data} is a
\code{"gifti"} object, it will overwrite the existing data type.}

\item{encoding}{One of "ASCII", "Base64Binary", or "GZipBase64Binary". If
\code{NULL} (default), will use the metadata if \code{data} is a GIFTI object,
or "ASCII" if the \code{data_type} is "NIFTI_TYPE_INT32" and
"GZipBase64Binary" if the \code{data_type} is "NIFTI_TYPE_FLOAT32". If not
\code{NULL} and \code{data} is a \code{"gifti"} object, it will overwrite the
existing data type.}

\item{endian}{"LittleEndian" (default) or "BigEndian". If \code{data} is a
\code{"gifti"} object, it will overwrite the existing endian.}

\item{col_names}{The names of each data column in \code{gii} (or entries in
\code{gii$data}).}

\item{label_table}{A data.frame with labels along rows. The row names should
be the label names. The column names should be among: "Key", "Red", "Green",
"Blue", and "Alpha". The "Key" column is required whereas the others are
optional (but very often included). Values in the "Key" column should be
non-negative integers, typically beginning with 0. The other columns should
be floating-point numbers between 0 and 1.

Although CIFTI files support a different label table for each data column,
GIFTI files only support a single label table. So this label table should be
applicable to each data column.}
}
\value{
Whether the GIFTI was successfully written
}
\description{
Write the data for the left or right cortex to a metric GIFTI file.
}
\seealso{
Other writing: 
\code{\link{separate_cifti}()},
\code{\link{write_cifti}()},
\code{\link{write_subcort_nifti}()},
\code{\link{write_surf_gifti}()},
\code{\link{write_xifti2}()}
}
\concept{writing}
