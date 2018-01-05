############################################################
####
####    SET-UP THE FILES AND FOLDERS TO CREATE R LIBRARY
####
############################################################

library(devtools)
library(roxygen2)
library(Rcpp)

# Input the name of the package here
mypackage <- "seminribm"

# path to the model code:
code.model.folder <- "./src/"

# Retrieve all relevant C++ files:
cppfiles <- system(paste0("ls ",code.model.folder,"*.cpp"),intern = TRUE)
hfiles   <- system(paste0("ls ",code.model.folder,"*.h"),  intern = TRUE)
cppfiles <- cppfiles[!grepl(pattern = "main",x = cppfiles)]   # <-- Remove "main" files
c.path <- c(cppfiles, hfiles)

# R files:
Rfiles <- paste0(code.model.folder,
                 c('run_seminribm.R'))

# This generates all the necessary files 
# when creating an R package from scrath 
# with the goal of interfacing with C++ (using Rcpp)
#?Rcpp::Rcpp.package.skeleton
Rcpp.package.skeleton(name =  mypackage,
                      example_code = FALSE,
                      author = "David Champredon",
                      cpp_files = c.path, 
                      code_files = Rfiles,
                      force = TRUE
)

# Documentation
devtools::document(pkg = mypackage, clean = TRUE)
