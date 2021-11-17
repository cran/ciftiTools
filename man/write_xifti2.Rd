% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/write_xifti2.R
\name{write_xifti2}
\alias{write_xifti2}
\title{Write a \code{"xifti"} object to GIFTI and NIFTI files}
\usage{
write_xifti2(
  xifti,
  brainstructures = NULL,
  cortexL_fname = NULL,
  cortexR_fname = NULL,
  subcortVol_fname = NULL,
  subcortLabs_fname = NULL,
  ROI_brainstructures = "all",
  ROIcortexL_fname = NULL,
  ROIcortexR_fname = NULL,
  ROIsubcortVol_fname = NULL,
  write_dir = NULL,
  verbose = FALSE
)
}
\arguments{
\item{xifti}{A \code{"xifti"} object.}

\item{brainstructures}{Character vector indicating which brain structure(s)
to obtain: \code{"left"} (left cortex), \code{"right"} (right cortex)
and/or \code{"subcortical"} (subcortex and cerebellum). Can also be
\code{"all"} (obtain all three brain structures). Default:
\code{c("left","right")} (cortex only).

If a brain structure is indicated but does not exist in the CIFTI file, a
warning will occur and that brain structure will be skipped.}

\item{cortexL_fname}{(Optional) GIFTI file names
(*.[func/label].gii) to save the [left/right] cortex data to. If not provided,
defaults to \code{"*[L/R].\[func/label\].gii"}, where * is the file name
component of \code{cifti_fname}. Will be written in \code{write_dir}.

dtseries and dscalar files should use "func", whereas dlabel files should
use "label".}

\item{cortexR_fname}{(Optional) GIFTI file names
(*.[func/label].gii) to save the [left/right] cortex data to. If not provided,
defaults to \code{"*[L/R].\[func/label\].gii"}, where * is the file name
component of \code{cifti_fname}. Will be written in \code{write_dir}.

dtseries and dscalar files should use "func", whereas dlabel files should
use "label".}

\item{subcortVol_fname}{(Optional) NIFTI file names to save
the subcortical [volume/labels] to. If not provided, defaults to
\code{"*[/.labels].nii"}, where * is the file name component of
\code{cifti_fname}. Will be written in \code{write_dir}.}

\item{subcortLabs_fname}{(Optional) NIFTI file names to save
the subcortical [volume/labels] to. If not provided, defaults to
\code{"*[/.labels].nii"}, where * is the file name component of
\code{cifti_fname}. Will be written in \code{write_dir}.}

\item{ROI_brainstructures}{Which ROIs should be obtained? \code{"all"}
(default) to obtain ROIs for each of the \code{brainstructures}. \code{NULL}
to not obtain any ROIs. This should be a subset of \code{brainstructures}.}

\item{ROIcortexL_fname}{(Optional) GIFTI file names
(*.[func/label].gii) to save the [left/right] cortex ROI to. If not provided,
defaults to \code{"*ROI_[L/R].\[func/label\].gii"}, where * is the file name component
of \code{cifti_fname}. The cortical ROIs typically represent the medial wall
mask, with values of 1 for in-ROI (non-medial wall) vertices and 0 for
out-of-ROI (medial wall) vertices. Will be written in \code{write_dir}.

dtseries and dscalar files should use "func", whereas dlabel files should
use "label".}

\item{ROIcortexR_fname}{(Optional) GIFTI file names
(*.[func/label].gii) to save the [left/right] cortex ROI to. If not provided,
defaults to \code{"*ROI_[L/R].\[func/label\].gii"}, where * is the file name component
of \code{cifti_fname}. The cortical ROIs typically represent the medial wall
mask, with values of 1 for in-ROI (non-medial wall) vertices and 0 for
out-of-ROI (medial wall) vertices. Will be written in \code{write_dir}.

dtseries and dscalar files should use "func", whereas dlabel files should
use "label".}

\item{ROIsubcortVol_fname}{(Optional) NIFTI file names to save
the subcortical ROI to. If not provided, defaults to
\code{"*ROI.nii"}, where * is the file name component of
\code{cifti_fname}. The subcortical ROI typically represents the volumetric
mask for the entire subcortical structure, with values of 1 for in-ROI
(in subcortex) voxels and 0 for out-of-ROI (not in subcortex) voxels. Will
be written in \code{write_dir}.}

\item{write_dir}{Where should the separated files be placed? \code{NULL}
(default) will write them to the current working directory.

\code{write_dir} must already exist, or an error will occur.}

\item{verbose}{Should occasional updates be printed? Default: \code{FALSE}.}
}
\value{
List of written files
}
\description{
Write metric or label GIFTIs for the cortical surface data and NIFTIs for the
subcortical labels and mask in a \code{"xifti"} object. Each present
brainstructure will be written; if a brainstructure is absent the
corresponding file is not written.
}