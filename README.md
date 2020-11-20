
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ciftiTools

<!-- badges: start -->

[![Travis build
status](https://travis-ci.com/mandymejia/ciftiTools.svg?branch=master)](https://travis-ci.com/mandymejia/ciftiTools)
[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/mandymejia/ciftiTools?branch=master&svg=true)](https://ci.appveyor.com/project/mandymejia/ciftiTools)
[![Coveralls test
coverage](https://coveralls.io/repos/github/mandymejia/ciftiTools/badge.svg)](https://coveralls.io/github/mandymejia/ciftiTools)
<!-- badges: end -->

Tools for reading and visualizing CIFTI brain imaging files.

CIFTI files contain brain imaging data in “gray-ordinates”, which
represent the gray matter as cortical surface vertices (left and right)
and subcortical voxels (cerebellum, basal ganglia, and other deep gray
matter).`ciftiTools` uses the Connectome Workbench to read CIFTI files
into R and apply common pre-processing steps (e.g. smoothing,
resampling). It also provides tools for visualizing the cortical surface
with GIFTI files, and for visualizing the subcortical volume.

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("mandymejia/ciftiTools")
```

## Vignette

See this link to view the tutorial:
<https://htmlpreview.github.io/?https://github.com/mandymejia/ciftiTools/blob/master/vignettes/ciftiTools_vignette.html>
