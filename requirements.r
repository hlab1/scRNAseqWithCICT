# requirements.r

install.packages(c(
  "caret",
  "Rtsne",
  "umap",
  "devtools",
  "hutils",
  "sn",
  "tidyr",
  "infotheo",
  "mpmi",
  "stringr",
  "RCy3",
  "igraph",
  "data.table",
  "plyr",
  "tidyverse",
  "dplyr",
  "zeallot",
  "fitdistrplus",
  "gendist",
  "ggplot2",
  "SimilarityMeasures",
  "PerformanceAnalytics",
  "Lmoments",
  "sna",
  "tidygraph",
  "graphlayouts",
  "statnet",
  "knitr"
))

# Alias functions
select <- dplyr::select
mutate <- dplyr::mutate
arrange <- dplyr::arrange
filter <- dplyr::filter
group_by <- dplyr::group_by
`%<-z%` <- zeallot::`%<-%`

# Define constants
PSUDEO_ZERO <- 0.000001
