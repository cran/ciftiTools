% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/edit_mask_surf.R
\name{edit_mask_surf}
\alias{edit_mask_surf}
\alias{erode_mask_surf}
\alias{dilate_mask_surf}
\title{Edit mask on surface}
\usage{
edit_mask_surf(
  x,
  mwall = NULL,
  surf = NULL,
  hemisphere = c("left", "right"),
  do = c("erode", "dilate", "borders"),
  depth = 1
)

erode_mask_surf(
  x,
  mwall = NULL,
  surf = NULL,
  hemisphere = c("left", "right"),
  depth = 1
)

dilate_mask_surf(
  x,
  mwall = NULL,
  surf = NULL,
  hemisphere = c("left", "right"),
  depth = 1
)
}
\arguments{
\item{x, mwall}{Vector of the data mask to edit, and the medial wall mask.
These can be specified in two ways. First, \code{mwall} can be a logical vector
with each entry corresponding to a vertex as the cortical surface, and
using \code{FALSE} values to indicate medial wall vertices. In this first case,
\code{x} should then be a logical vector with each entry corresponding to a
\code{TRUE} value in \code{mwall}. \code{TRUE} values in \code{x} should
indicate the mask to be edited.

Second, \code{mwall} can be \code{NULL} (default) in which case \code{x}
should then be a logical vector with each entry corresponding to a vertex
on the cortical surface. \code{TRUE} values in \code{x} should indicate the
mask to be edited.

In either case, \code{xii$data$cortex_left[,1]} and
\code{xii$meta$cortex$medial_wall_mask$left} should work.}

\item{surf, hemisphere}{Provide one: the surface in the same resolution as the
data, or the name of the hemisphere of the surface to resample and use
(default: resample the left surface).}

\item{do}{\code{"erode"} (default), \code{"dilate"}, or \code{"borders"}.
\code{"erode"} removes faces with at least one vertex not inside the mask.
\code{"dilate"} adds faces with at least one vertex inside the mask.
\code{"borders"} obtains the vertices inside the mask which share a face
with at least one vertex not inside the mask.}

\item{depth}{How many iterations of the edit? Default: \code{1}. Does not
apply to \code{"borders"}.}
}
\value{
\code{x} after erosion or dilation.
}
\description{
Erode, dilate, or get the borders of a mask along the cortical surface
}
\details{
The depth of the edit is determined by the number of edges between the
vertices. To erode or dilate based on spatial distance (mm), see
\code{-cifti-dilate} and \code{-cifti-erode}.
}
\seealso{
Other surface-related: 
\code{\link{add_surf}()},
\code{\link{boundary_mask_surf}()},
\code{\link{even_vert_samp}()},
\code{\link{is.surf}()},
\code{\link{load_surf}()},
\code{\link{mask_surf}()},
\code{\link{read_surf}()},
\code{\link{resample_surf}()},
\code{\link{rotate_surf}()},
\code{\link{surf_area}()},
\code{\link{view_surf}()},
\code{\link{write_surf_gifti}()}
}
\concept{surface-related}
