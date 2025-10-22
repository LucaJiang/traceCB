# run pipeline of gmm

id_merge=$(sbatch --parsable -w hpcnode022 shell/run_merge.sh)
id_ld=$(sbatch --parsable --dependency=afterok:${id_merge} -w hpcnode022 shell/run_ld.sh)
sbatch --dependency=afterok:${id_ld} -w hpcnode022 shell/run_gmm.sh

sbatch --dependency=afterok:1220916 -w hpcnode022 shell/run_gmm.sh
sbatch --dependency=afterok:1220164 -w hpcnode022 coloc/run_coloc.sh

id_merge=$(sbatch --parsable -w hpcnode022 shell/run_merge.sh)
id_ld=$(sbatch --parsable -w hpcnode022 shell/run_ld.sh)
sbatch -w hpcnode022 shell/run_gmm.sh
