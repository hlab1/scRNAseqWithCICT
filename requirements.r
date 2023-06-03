# requirements.r

install.packages("BiocManager")
BiocManager::install("minet")
BiocManager::install('RCy3')
BiocManager::install(c('impute', 'preprocessCore', 'GO.db', 'AnnotationDbi'))

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
  "knitr",
  "arules",
  "doFuture",
  "here",
  "GGally",
  "doParallel",
  "doSNOW",
  "WGCNA"
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
