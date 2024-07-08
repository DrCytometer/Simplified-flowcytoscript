# set up for simplified flowcytoscript

# Welcome and check R status
cat("This script will try to help you update R\n
    \n")
Sys.sleep(2)
cat(
  "If you can do this yourself, that may work better, particularly\n
  for non-Windows users.\n
  \n"
)

install.packages('installr')
suppressPackageStartupMessages( library(installr) )
updateR()


