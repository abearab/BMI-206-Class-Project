qsub -cwd -pe smp 6 -l mem_free=2G -l scratch=5G -l h_rt=24:00:00 run_T_nocvar.sh
