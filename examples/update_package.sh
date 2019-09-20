rm src/*.o
rm src/*.dll
Rscript.exe -e 'Rcpp::compileAttributes()'
Rscript.exe -e 'devtools::document()'
R CMD INSTALL --no-multiarch --with-keep.source ../bcf-1