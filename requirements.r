# requirements.r

# BioConductor
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
  "WGCNA",
  "precrec"
))

# h2o
if ("package:h2o" %in% search()) { detach("package:h2o", unload=TRUE) }
if ("h2o" %in% rownames(installed.packages())) { remove.packages("h2o") }
pkgs <- c("RCurl","jsonlite")
for (pkg in pkgs) {
    if (! (pkg %in% rownames(installed.packages()))) { install.packages(pkg) }
}
install.packages("h2o", type="source", repos=(c("http://h2o-release.s3.amazonaws.com/h2o/latest_stable_R")))



# Alias functions
select <- dplyr::select
mutate <- dplyr::mutate
arrange <- dplyr::arrange
filter <- dplyr::filter
group_by <- dplyr::group_by
`%<-z%` <- zeallot::`%<-%`

# Define constants
PSUDEO_ZERO <- 0.000001
