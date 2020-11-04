## The IMP analysis tutorial for RNA Polymerase III

Last updated: Nov. 2020

## Purpose

This tutorial walks the user through analyzing, validating and depositing integrative models produced by the Integrative Modeling Platform. 

## Running the tutorial.

The tutorial contains code to perform a complete integrative modeling pipeline from model building through archiving the models and protocols in mmCIF format.  Examples for a partially sampled and more extensively sampled protocols are provided. 

### 0. Installation
Install the tutorial and IMP as described in the [main page README](../README.md)

### 1. Building the RNA Polymerase III model, scoring function and sampling

This step is optional. One may skip directly to analysis (Step 2) using the provided output files

* Navigate to `./modeling`

* Change line 7 in `./run_rnapolii_modeling.sh` to point to your python installation

* run `./run_rnapolii_modeling.sh output_dir N n_steps`
  * `output_dir` is the prefix of the output directory
  * `N` is the number of independent sampling runs to do (at least 2, up to the number of processors you are willing to sacrifice) 
  * `n_steps` is the number of Monte-Carlo steps to run. 

Each 1000 steps takes a few to 20 minutes depending on your machine. The script will produce a significant amount of output to screen and create directories `output_dir0`, `output_dir1`, `...` which contain the output that will be analyzed. 

To sufficiently sample the system, at  `N>=4` and `n_steps>=50000` should be performed. Analysis can still be done on smaller samples. 

### 2. Score-based clustering

First, we do a quick assessment of the output for sampling mixing by clustering based on score values. 

* Navigate to `./analysis`
* run `python ./run_analysis_trajectories.py ../modeling/ output_dir`
  * `../modeling` is the path to the output directories
  * `output_dir` is the prefix of the output directories (same as above or `example` if you wish to analyze the supplied output. 
* Look at `model_analysis/plot_clustering_scores.png` (below)

<img src="analysis/model_analysis/plot_clustering_scores.png" width="350">

Here, we see the output for the `example` output, distinct, non-overlapping clusters of solutions.  This can be indicative of poor sampling.  A more complete sampling is seen for the same plot from `./analysis_extensive`

<img src="analysis_extensive/model_analysis/plot_clustering_scores.png" width="350">
