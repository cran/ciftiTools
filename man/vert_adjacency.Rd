% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_surf.R
\name{vert_adjacency}
\alias{vert_adjacency}
\title{Vertex Adjacency Matrix}
\usage{
vert_adjacency(faces, v1, v2 = NULL)
}
\arguments{
\item{faces}{An \eqn{F \times 3} matrix, where each row contains the vertex
indices for a given triangular face in the mesh. \eqn{F} is the number of
faces in the mesh.}

\item{v1, v2}{The first and second set of vertices. These are logical vectors
the same length as \code{vertices} indicating the vertices in each set.
If \code{v2} is \code{NULL} (default), set \code{v2} to \code{v1}. Can
alternatively be a vector if integers corresponding to vertex indices.}
}
\value{
Adjacency matrix
}
\description{
Make adjacency matrix between two sets of vertices on the same mesh.
}
\keyword{internal}
