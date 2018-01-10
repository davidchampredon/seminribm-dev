rm -rf seminribm*
rm -rf ./lib
Rscript build_library.R
cp Makevars seminribm/src 
R CMD build seminribm
#R CMD check seminribm
mkdir ./lib
R CMD INSTALL -l ./lib seminribm
echo --------------------------------------------------------
echo -- done --

