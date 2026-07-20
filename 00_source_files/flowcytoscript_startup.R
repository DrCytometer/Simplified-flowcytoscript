# start up for simplified flowcytoscript

# dummy variables for report
experimental.system <- "experimental system"
tissue.type <- "tissue type"
clustering.method <- NULL
p.value.message <- 0
trex.condition <- "None"
fcs.mfi.stats.dir <- "./marker_stats/"
fcs.channel.label <- "None"
analysis.calc.time <- 0
setup.time <- 0
fcs.channel <- NULL
fcs.cluster.label <- NULL
setup.start.time <- Sys.time()
input.file.type <- 2
run.crossentropy <- 2
clustering.method <- 2
species.options <- 1
dmrd.knn <- NULL

# Welcome and check R status
if (Be.Chatty == TRUE){
  cat( 
"Welcome to flowcytoscript!\n
This simplified version of the Liston Lab flow cytometry analysis\n
pipeline will try to take care of as much as possible.\n
\n"
  )
  Sys.sleep(message.delay.time)
  cat(  
"We're going to have you tell us what your groups are,\n
which markers you want to analyze, and how many cells\n
you want to work with.\n 
\n"
  )
  Sys.sleep(message.delay.time)
  cat( 
"After that, we'll try to cluster\n
your data, and provide you with visualizations in the forms\n
of tSNE, UMAP, PCA, heatmaps and barcharts.\n
\n"
  )
  Sys.sleep(message.delay.time)
  cat(
"For best results, make sure your R and RStudio are up-to-date.\n
\n"
  )
  
  Sys.sleep(message.delay.time)
  
  # check packages, give warning, install missing packages--------------
  
  cat(
"Now we're going to try to install any of the required packages\n
that you don't already have installed.\n
\n"
  )
  Sys.sleep(message.delay.time)
}

required.packages <- c(
  "digest", "dunn.test", "ggplot2", "ggridges", "ggrepel", "ggtext",
  "RColorBrewer", "Rtsne", "uwot", "dplyr", "tidyr", "RcppHNSW",
  "parallel", "data.table", "remotes", "pak",
  "coda", "emmeans", "EmbedSOM", "ConsensusClusterPlus", "flowCore",
  "flowWorkspace", "readxl", "scattermore", "parallelly"
)

# CRAN packages, plus Bioconductor packages prefixed "bioc::" so a single
# pak::pkg_install() call can resolve and install both in one dependency pass.
cran.packages <- c(
  "digest", "dunn.test", "ggplot2", "ggridges", "ggrepel", "ggtext",
  "RColorBrewer", "Rtsne", "uwot", "dplyr", "tidyr", "RcppHNSW",
  "parallel", "data.table", "remotes", "coda", "emmeans",
  "EmbedSOM", "readxl", "scattermore", "devtools", "parallelly"
)

bioconductor.packages <- c("bioc::ConsensusClusterPlus", "bioc::flowCore", "bioc::flowWorkspace")

# FastPG (the Phenograph engine) is handled separately below, since it's a
# GitHub package that always compiles from source and commonly fails on
# macOS due to missing OpenMP support in Apple's default compiler. It has
# its own try/fallback logic rather than being a hard requirement here.

if (!requireNamespace("pak", quietly = TRUE)) install.packages("pak")

if ( length(setdiff(required.packages, rownames(installed.packages()))) !=0 ){
  cat( "It looks like this might be your first time using flowcytoscript.\n
We need to install some packages before we get started.\n
This might take a few minutes.\n
If you're prompted to update, please do so.\n")
  
  pak::pkg_install(setdiff(cran.packages, rownames(installed.packages())))
  
  if ( length(setdiff(c("ConsensusClusterPlus", "flowCore", "flowWorkspace"), rownames(installed.packages()))) !=0 ){
    pak::pkg_install(bioconductor.packages)
  }
  
}

# Phenograph clustering engine: FastPG (preferred) with an Rphenograph fallback ----
# FastPG is fast but must compile C++ code with OpenMP, which fails on many Mac
# systems because Apple's default compiler doesn't ship OpenMP support out of the
# box. Rather than requiring every Mac user to install Homebrew's libomp and hand-
# edit their Makevars, we try FastPG first and silently fall back to Rphenograph
# (a slower, pure R/Rcpp implementation of the same Louvain-clustering approach,
# with no OpenMP dependency) if FastPG can't be installed. FlowSOM users never
# need either of these packages, so this never blocks that path.

fcs.phenograph.engine <- NA_character_

if (Be.Chatty == TRUE){
  cat("\nSetting up the Phenograph clustering engine (trying FastPG first,\n
with an automatic fallback if that isn't available on your system)...\n\n")
}

fastpg.ok <- suppressWarnings(tryCatch({
  if (!requireNamespace("FastPG", quietly = TRUE)) {
    pak::pkg_install("sararselitsky/FastPG")
  }
  requireNamespace("FastPG", quietly = TRUE)
}, error = function(e) FALSE, warning = function(w) FALSE))

if (isTRUE(fastpg.ok)) {
  fcs.phenograph.engine <- "FastPG"
} else {
  if (Be.Chatty == TRUE){
    cat("\nFastPG could not be installed. This is a known, common issue on\n
macOS (Apple's default compiler doesn't include OpenMP, which FastPG needs\n
to compile) and isn't something wrong with your R installation.\n
Falling back to Rphenograph instead: same style of Phenograph clustering,\n
just slower, and it needs no OpenMP-capable compiler.\n\n")
  }
  
  rphenograph.ok <- suppressWarnings(tryCatch({
    if (!requireNamespace("Rphenograph", quietly = TRUE)) {
      pak::pkg_install("JinmiaoChenLab/Rphenograph")
    }
    requireNamespace("Rphenograph", quietly = TRUE)
  }, error = function(e) FALSE, warning = function(w) FALSE))
  
  if (isTRUE(rphenograph.ok)) {
    fcs.phenograph.engine <- "Rphenograph"
  } else if (Be.Chatty == TRUE){
    cat("\nRphenograph could not be installed either. Phenograph clustering\n
won't be available this session, but FlowSOM will still work fine.\n\n")
  }
}

if (length(setdiff(required.packages, rownames(installed.packages()))) !=0 ){
  cat("\n")
  cat("Installation appears to have failed for one or more packages.\n
You'll need to get help before proceeding.\n
      \n")
} else {
  cat("\n")
  cat("That's all done. Now, on to the analysis!\n
      \n")
}

stopifnot("Not all packages were installed" = 
            length(setdiff(required.packages, rownames(installed.packages()))) == 0 )
