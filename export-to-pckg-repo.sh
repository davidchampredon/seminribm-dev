### Export to LIVE package Github repository
### Does NOT export the DESCRIPTION, NAMESPACE and LICENSE files

### R files:
cp seminribm/R/*.R ../seminribm/R

### C++ source files:
cp seminribm/src/* ../seminribm/src

### Documentation:
cp seminribm/man/* ../seminribm/man
cp vignettes/*.Rmd ../seminribm/vignettes

### Data added to the package:
cp data/*.RData ../seminribm/data
cp data/data.R ../seminribm/R
