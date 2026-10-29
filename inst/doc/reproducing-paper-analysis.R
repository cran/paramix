## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5,
  dev.args = list(bg = "white")
)
if (isNamespaceLoaded("thematic")) {
  getFromNamespace("thematic_off", "thematic")()
}

## ----copy-pipeline, eval = FALSE----------------------------------------------
# # Choose a destination directory
# destination <- file.path("~", "Downloads", "paramixdemo")
# 
# if (!dir.exists(destination)) {
#   dir.create(destination, recursive = TRUE)
# }
# 
# # Copy the bundled analysis files
# pipeline_dir <- system.file("analysis", package = "paramix")
# analysis_files <- list.files(
#   pipeline_dir,
#   full.names = TRUE,
#   recursive = FALSE,
#   include.dirs = TRUE
# )
# 
# file.copy(from = analysis_files, to = destination, recursive = TRUE)

