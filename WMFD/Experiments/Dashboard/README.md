# ARI dashboard (noise=drop)
Contient :
- Scripts : grid_experiments.py, gptree_cluster_refined.py, step1b_metric_wmfd.py, step1c_clustering_kmedoids.py, plot_dashboard.py
- Résultats agrégés : results_all.csv
- Figure : dashboard_final.png

Reproduire rapidement (PowerShell) :
python grid_experiments.py --scripts_dir . --workdir A_noise_dense `
  --k_list 1,2,3,4,5 --L_list 10,20,30,40,60,80,100 --n_list 8,16,32,64,128 `
  --noise_list 0.10,0.25,0.50,0.75 --noise_mode drop `
  --plevel_list 0.30,0.40,0.50,0.60,0.70 --repeats 3 --metric wmfd `
  --seed 123 --out results_all.csv
python plot_dashboard.py --results results_all.csv `
  --out dashboard_final.png --noise_ref 0.50 --noise_low 0.10 --noise_high 0.75
