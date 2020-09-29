# Analysis of RNA Polymerase II stak modeling

Computation of the sampling precision proceeds through the following three commands:

`python ./run_analysis_trajectories.py ../modeling/ example`

`python ./run_extract_models.py 0`

```
python -m IMP.sampcon.exhaust -n rnapol \
--rmfA ./model_analysis/A_models_clust0.rmf3 \
--rmfB ./model_analysis/B_models_clust0.rmf3 \
--scoreA ./model_analysis/A_models_clust0.txt \ 
--scoreB ./model_analysis/B_models_clust0.txt \
-d density_rnapol.txt \
-m cpu_omp -c 4 
-g 2.0 -gp
```
