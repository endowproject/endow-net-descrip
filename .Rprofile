# .Rprofile
#
# Defines the path variables that the net_descriptives/ scripts expect.
# R runs this automatically at startup when launched from this folder
# (e.g. by opening the RStudio project, or via setwd() + restarting R).
#
# NOTE: a project-level .Rprofile is used INSTEAD of ~/.Rprofile. If you rely
# on settings in your home .Rprofile, uncomment the next line to keep them:
# if (file.exists("~/.Rprofile")) source("~/.Rprofile")

EndowGitHub      <- normalizePath(getwd()) # repo root
EndowDatabase    <- file.path(EndowGitHub, "FauxData") # networks/*_networks.csv, people_observations/*.csv
EndowDerivedData <- file.path(EndowGitHub, "FauxData") # *.rdata / *.rds intermediate objects
EndowDropbox     <- file.path(EndowGitHub, "FauxDataOutputs") # Figures/ and Tables/ are written here