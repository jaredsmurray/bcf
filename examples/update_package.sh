rm src/*.o
rm src/*.so
rm src/*.dll
Rscript.exe -e 'Rcpp::compileAttributes()'
Rscript.exe -e 'pkgbuild::compile_dll()'
Rscript.exe -e 'devtools::document()'
R CMD INSTALL --no-multiarch --with-keep.source ../bcf-1

# Rscript.exe -e 'pkgdown::build_site()'

# pkgdown::build_site()
# devtools::check()
# rhub::check()

# Rscript.exe examples/simple_example.R

# Rscript.exe examples/test_pred2.R


# Rscript.exe examples/simple_example.R 2>&1 | tee examples/bcf_run_log.txt 