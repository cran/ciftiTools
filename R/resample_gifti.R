#' Resample a GIFTI file (with its ROI)
#'
#' Perform spatial resampling of GIFTI data on the cortical surface (metric
#'  and label), or of GIFTI surface geometry data itself.
#'
#' @param original_fname The GIFTI file to resample.
#' @param target_fname Where to save the resampled file.
#' @param hemisphere \code{"left"} (default) or \code{"right"}. An error will
#'  occur if the hemisphere indicated in the GIFTI metadata does not match.
#' @param file_type \code{"metric"}, \code{"label"}, \code{"surf"}, or \code{NULL}
#'  (default) to infer from \code{original_fname}.
#' @param original_res The resolution of the original file. If \code{NULL}
#'  (default), infer from the file. Alternatively, provide
#'  \code{sphere_original_fname} which will override \code{original_res}.
#'
#'  In general, \code{original_res} should be used when the original file is
#'  in registration with the spheres created by the Workbench command
#'  \code{-surface-create-sphere}, and \code{sphere_original_fname} should be
#'  used when it is not compatible.
#' @param resamp_res Target resolution for resampling (number of
#'  cortical surface vertices per hemisphere). Alternatively, provide
#'  \code{sphere_target_fname} which will override \code{resamp_res}.
#'
#'  In general, \code{resamp_res} should be used when the target file will be
#'  in registration with the spheres created by the Workbench command
#'  \code{-surface-create-sphere}, and \code{sphere_target_fname} should be
#'  used when it is not compatible.
#' @param resamp_method \code{"barycentric"} (default) or \code{"adaptive"}
#'  resampling. These options correspond to the Workbench command options
#'  \code{"BARYCENTRIC"} and \code{"ADAP_BARY_AREA"}, respectively.
#'
#'  While adaptive resampling is recommended for metric or label
#'  data, it requires that \code{area_original_fname} be provided.
#' @param area_original_fname,area_target_fname File paths to the surfaces to
#'  use for vertex area correction during adaptive resampling. (Ignored if
#'  resampling with the barycentric method.) \code{area_original_fname} should
#'  match the current resolution of the data, and \code{area_target_fname}
#'  should match \code{resamp_res}. If \code{area_target_fname} is not provided,
#'  \code{area_original_fname} will be resampled with the barycentric method,
#'  and the result will be used as \code{area_target_fname}.
#' @param ROIcortex_original_fname The name of the ROI file corresponding to
#'  \code{original_fname}. Leave as \code{NULL} (default) if this doesn't exist
#'  or shouldn't be resampled.
#' @param ROIcortex_target_fname The name of the resampled ROI file. Only
#'  applicable if \code{ROIcortex_original_fname} is provided.
#' @param sphere_original_fname,sphere_target_fname File paths to the sphere
#'  surfaces in the original and target resolutions. If possible, the simpler
#'  arguments \code{original_res} and \code{resamp_res} can be used instead. But
#'  those depend on the surface being compatible with that created by
#'  \code{-surface-create-sphere}, which isn't always true. Therefore
#'  \code{sphere_original_fname} and \code{sphere_target_fname} can be used if
#'  needed.
#' @param read_dir Directory to append to the path of every file name in
#'  \code{original_fname} and \code{ROIcortex_original_fname}. If \code{NULL}
#'  (default), do not append any directory to the path.
#' @param write_dir Directory to append to the path of every file name in
#'  \code{target_fname} and \code{ROIcortex_target_fname}. If \code{NULL}
#'  (default), do not append any directory to the path.
#'
#' @return The resampled GIFTI file name, invisibly
#'
#' @importFrom gifti readgii
#'
#' @family gifting
#' @export
#'
#' @section Connectome Workbench:
#' This function interfaces with the \code{"-metric-resample"}, \code{"-label-resample"},
#'  and/or \code{"-surface-resample"} Workbench commands, depending on the input.
#'
resample_gifti <- function(
  original_fname, target_fname, hemisphere=c("left", "right"),
  file_type=NULL, original_res=NULL,
  resamp_res=NULL, resamp_method=c("barycentric", "adaptive"),
  area_original_fname=NULL, area_target_fname=NULL,
  ROIcortex_original_fname=NULL, ROIcortex_target_fname=NULL,
  sphere_original_fname=NULL, sphere_target_fname=NULL,
  read_dir=NULL, write_dir=NULL) {

  # ----------------------------------------------------------------------------
  # Check arguments. -----------------------------------------------------------
  # ----------------------------------------------------------------------------

  # File type
  if (is.null(file_type)) {
    if (grepl("func.gii", original_fname, fixed=TRUE)) {
      file_type <- "metric"
    } else if (grepl("label.gii", original_fname, fixed=TRUE)) {
      file_type <- "label"
    } else if (grepl("surf.gii", original_fname, fixed=TRUE)) {
      file_type <- "surface"
    } else {
      stop(paste(
        "Could not infer file type of ", original_fname,
        ". Please set the file_type argument."
      ))
    }
  }
  file_type <- match.arg(file_type, c("metric", "label", "surface"))

  # Original & target file names
  original_fname <- format_path(original_fname, read_dir, mode=4)
  stopifnot(file.exists(original_fname))
  target_fname <- format_path(target_fname, write_dir, mode=2)

  # Hemisphere
  hemisphere <- match.arg(hemisphere, c("left", "right"))
  if (file_type == "surface") {
    surf <- make_surf(original_fname, hemisphere)
  }

  # Original ROI & target ROI file names
  do_ROI <- !is.null(ROIcortex_original_fname)
  if (do_ROI) {
    ROIcortex_original_fname <- format_path(
      ROIcortex_original_fname, read_dir, mode=4
    )
    stopifnot(file.exists(ROIcortex_original_fname))
    if (is.null(ROIcortex_target_fname)) {
      ROIcortex_target_fname <- cifti_component_suffix(
        paste0("ROIcortex", switch(hemisphere, left="L", right="R")),
        switch(file_type, metric="func", label="func", surface="surf")
      )
    }
    ROIcortex_target_fname <- format_path(
      ROIcortex_target_fname, write_dir, mode=2
    )
  }

  # Check resamp_method and area files
  resamp_method <- match.arg(resamp_method, c("barycentric", "adaptive"))
  if (resamp_method == "adaptive") {
    if (is.null(area_original_fname)) {
      stop("Adaptive sampling requires `area_original_fname`.")
    }
  }

  # Compatible args
  if (do_ROI & file_type=="surface") {
    stop("do_ROI AND file_type=='surface', but surface files do not use ROI.")
  }

  # ----------------------------------------------------------------------------
  # Spheres & Area surface. ----------------------------------------------------
  # ----------------------------------------------------------------------------

  # Spheres

  tdir <- tempdir()

  if (is.null(sphere_original_fname)) {

    # Check `original_res`
    if (is.null(original_res)) {
      gii <- readgii(original_fname)
      original_res <- switch(file_type,
        metric = nrow(gii$data[[1]]),
        label = nrow(gii$data[[1]]),
        surface = nrow(gii$data$pointset)
      )
    }
    stopifnot(is.numeric(original_res) && original_res > 0)

    sphereL_original_fname <- format_path(paste0("sphereL_", original_res, ".surf.gii"), tdir, mode=2)
    sphereR_original_fname <- format_path(paste0("sphereR_", original_res, ".surf.gii"), tdir, mode=2)
    write_spheres(sphereL_original_fname, sphereR_original_fname, original_res)

    sphere_original_fname <- switch(hemisphere,
      left = sphereL_original_fname,
      right = sphereR_original_fname
    )

    # Check that the sphere is of compatible resolution.
    sphere_original_res <- nrow(readgii(sphere_original_fname)$data$pointset)
    if (sphere_original_res != original_res) {
      stop(paste(
        "Unable to create sphere with matching resolution to input GIFTI.",
        "Please provide `sphere_original_fname` (and `sphere_target_fname`)",
        "instead of `original_res` (and `resamp_res`)."
      ))
    }
  }

  if (is.null(sphere_target_fname)) {

    # Check `resamp_res`.
    if (is.null(resamp_res)) {
      stop("Provide either `sphere_target_fname` or `resamp_res`.")
    }
    stopifnot(is.numeric(resamp_res) && resamp_res > 0)

    sphereL_target_fname <- format_path(paste0("sphereL_", resamp_res, ".surf.gii"), tdir, mode=2)
    sphereR_target_fname <- format_path(paste0("sphereR_", resamp_res, ".surf.gii"), tdir, mode=2)
    write_spheres(sphereL_target_fname, sphereR_target_fname, resamp_res)

    sphere_target_fname <- switch(hemisphere,
      left = sphereL_target_fname,
      right = sphereR_target_fname
    )
  }

  # Area surface for adaptive resampling
  if (resamp_method == "adaptive") {
    if (is.null(area_target_fname)) {
      area_target_fname <- paste0(
        tempfile(), "_area_target_", hemisphere, ".surf.gii"
      )
      resample_gifti(
        original_fname=area_original_fname,
        target_fname=area_target_fname,
        hemisphere=hemisphere,
        file_type="surf",
        original_res=original_res,
        resamp_res=resamp_res,
        sphere_original_fname=sphere_original_fname,
        sphere_target_fname=sphere_target_fname
      )
    }
  }

  # ----------------------------------------------------------------------------
  # Make and run command. ------------------------------------------------------
  # ----------------------------------------------------------------------------

  resamp_method_wb <- switch(resamp_method,
    adaptive="ADAP_BARY_AREA",
    barycentric="BARYCENTRIC"
  )

  if (file_type=="surface" && resamp_res>original_res) {
    ciftiTools_warn("Upsampling a surface is not recommended, if avoidable.")
  }

  cmd_name <- switch(file_type,
    metric="-metric-resample",
    label="-label-resample",
    surface="-surface-resample"
  )

  cmd <- paste(
    cmd_name,
    sys_path(original_fname), sys_path(sphere_original_fname),
    sys_path(sphere_target_fname), resamp_method_wb, sys_path(target_fname)
  )
  if (resamp_method_wb=="ADAP_BARY_AREA") {
    cmd <- paste(
      cmd,
      ifelse(endsWith(area_original_fname, "shape.gii"), "-area-metrics", "-area-surfs"),
      sys_path(area_original_fname),
      sys_path(area_target_fname)
    )
  }
  if (do_ROI) {
    cmd <- paste(
      cmd,
      "-current-roi", sys_path(ROIcortex_original_fname),
      "-valid-roi-out", sys_path(ROIcortex_target_fname)
    )
  }
  run_wb_cmd(cmd)

  if (do_ROI) {
    out <- c(target_fname, ROIcortex_target_fname)
    names(out) <- c(file_type, "ROI")
  } else {
    out <- target_fname
    names(out) <- file_type
  }
  invisible(out)
}

#' Generate GIFTI sphere surface files
#'
#' This function generates a pair of GIFTI vertex-matched left and right spheres
#'  in the target resolution. These are required for resampling CIFTI and GIFTI
#'  files.
#'
#' @param sphereL_fname File path to left-hemisphere spherical GIFTI to be
#'  created
#' @param sphereR_fname File path to right-hemisphere spherical GIFTI to be
#'  created
#' @inheritParams resamp_res_Param_required
#' @param write_dir (Optional) directory to place the sphere files in. If
#'  \code{NULL} (default), do not append any directory to the sphere file paths.
#'
#' @return The names of the written sphere files, invisibly
#'
#' @keywords internal
#'
write_spheres <- function(
  sphereL_fname, sphereR_fname, resamp_res,
  write_dir=NULL) {

  sphereL_fname <- format_path(sphereL_fname, write_dir, mode=2)
  sphereR_fname <- format_path(sphereR_fname, write_dir, mode=2)

  resamp_res <- format(resamp_res, scientific=FALSE)

  # Make left
  run_wb_cmd(
    paste("-surface-create-sphere", resamp_res[1], sys_path(sphereL_fname)),
  )

  # Make right
  if (length(resamp_res)==1 || length(unique(resamp_res))==1) {
    run_wb_cmd(
      paste("-surface-flip-lr", sys_path(sphereL_fname), sys_path(sphereR_fname)),
    )
  } else {
    sphereR_temp <- paste0(tempfile(), ".surf.gii")
    run_wb_cmd(
      paste("-surface-create-sphere", resamp_res[2], sphereR_temp),
    )
    run_wb_cmd(
      paste("-surface-flip-lr", sphereR_temp, sys_path(sphereR_fname)),
    )
  }

  # Set structure
  run_wb_cmd(
    paste("-set-structure", sys_path(sphereL_fname), "CORTEX_LEFT"),
  )
  run_wb_cmd(
    paste("-set-structure", sys_path(sphereR_fname), "CORTEX_RIGHT"),
  )

  invisible(list(sphereL_fname=sphereL_fname, sphereR_fname=sphereR_fname))
}

#' @rdname resample_gifti
#' @export
resampleGIfTI <- function(
  original_fname, target_fname, hemisphere,
  file_type=NULL, original_res=NULL, resamp_res,
  ROIcortex_original_fname=NULL, ROIcortex_target_fname=NULL,
  read_dir=NULL, write_dir=NULL){

  resample_gifti(
    original_fname, target_fname, hemisphere,
    file_type, original_res, resamp_res,
    ROIcortex_original_fname, ROIcortex_target_fname,
    read_dir, write_dir
  )
}

#' @rdname resample_gifti
#' @export
resamplegii <- function(
  original_fname, target_fname, hemisphere,
  file_type=NULL, original_res=NULL, resamp_res,
  ROIcortex_original_fname=NULL, ROIcortex_target_fname=NULL,
  read_dir=NULL, write_dir=NULL){

  resample_gifti(
    original_fname, target_fname, hemisphere,
    file_type, original_res, resamp_res,
    ROIcortex_original_fname, ROIcortex_target_fname,
    read_dir, write_dir
  )
}
