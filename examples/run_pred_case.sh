# source scripts/remove_all_compiled_files.sh
#source scripts/remove_some_compiled_files.sh

source scripts/install_bcf2.sh 2>&1 | tee scripts/bcf_install_log.txt 
Rscript.exe examples/test_pred2.R 2>&1 | tee scripts/bcf_run_log.txt 