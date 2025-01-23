## ----include=FALSE------------------------------------------------------------
library(knitr)
knitr::opts_chunk$set(autodep = TRUE, cache = FALSE)

## -----------------------------------------------------------------------------
# Check if package installed. If not, install it.
if(!require('ciftiTools', quietly=TRUE)){
  install.packages('ciftiTools')
  # devtools::install_github('mandymejia/ciftiTools') # development version
}

## -----------------------------------------------------------------------------
library(ciftiTools)

## ----eval=FALSE---------------------------------------------------------------
# # Replace '~/Applications' with the actual path to the
# #   Connectome Workbench folder on your computer. If
# #   successful, the Workbench executable path will be printed.
# ciftiTools.setOption('wb_path', '~/Applications')

## -----------------------------------------------------------------------------
xii <- load_parc() # Loads the Schaefer 100 parcellation.
xii # Summary of the `"xifti"` object.

## -----------------------------------------------------------------------------
library(rgl)
rgl::setupKnitr()

# Sometimes the first OpenGL window does not render properly.
rgl::open3d(); rgl::close3d()

# These are also required.
library(manipulateWidget)
library(ggpubr)

## ----fig.cap="Schaefer 100 parcellation", rgl=TRUE, format="jpg", fig.height=2.1, fig.width=2.5----
# Normally `cex.title` doesn't need to be set, as it defaults to a good choice.
#   But when knitting static images this way, the default becomes a bit too big
#   based on how knitting works.
view_xifti_surface(xii, idx=1, title='Schaefer 100', cex.title=1.2)

## -----------------------------------------------------------------------------
xii <- add_surf(xii, "midthickness", "midthickness")
xii

## -----------------------------------------------------------------------------
xii <- remove_xifti(xii, c("cortex_right", "surf_right"))
xii

## ----fig.cap="Plotting the FPole parcel", rgl=TRUE, format="jpg", fig.height=2, fig.width=1.3----
label_to_viz <- "17networks_LH_DefaultB_FPole_1"
key_idx <- which(rownames(xii$meta$cifti$labels$parcels)==label_to_viz)
key <- xii$meta$cifti$labels$parcels$Key[key_idx]
xii <- transform_xifti(xii, function(v){ifelse(v==key, v, 0)})
view_xifti_surface(xii)

## -----------------------------------------------------------------------------
citation("ciftiTools")

