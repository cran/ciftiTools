% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parc_utils.R
\name{parc_borders}
\alias{parc_borders}
\title{Parcellation borders}
\usage{
parc_borders(parc, surf = NULL, hemisphere = c("left", "right"))
}
\arguments{
\item{parc}{Integer vector the same length as the number of vertices. Each
entry indicates the parcel that vertex belongs to.}

\item{surf}{The surface which the vertices belong to, or just the \code{"faces"}
component (\eqn{F \times 3} matrix where each row indicates the vertices which
comprise a face). If not provided, the (resampled) default \code{hemisphere}
surface included with \code{ciftiTools} will be used.}

\item{hemisphere}{Only used to choose which default surface to use if
\code{is.null(surf)}. Should be \code{"left"} (default) or \code{"right"}.}
}
\value{
Logical vector the same length as \code{parc} indicating if the
vertex lies on a border.
}
\description{
Identify vertices which lie on the border of different parcels.
}
\seealso{
Other parcellation-related: 
\code{\link{apply_parc}()},
\code{\link{load_parc}()},
\code{\link{load_sub_parc}()},
\code{\link{parc_add_subcortex}()},
\code{\link{parc_vals_to_xifti}()}
}
\concept{parcellation-related}
