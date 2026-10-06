Adaptive Causal Variable Selection via HSIC and Polarization Decoupling (HSICADL)

This repository contains all code and data required to reproduce the numerical simulations, real data empirical processing, and empirical analysis presented in the paper "Adaptive Causal Variable Selection via HSIC and Polarization Decoupling".

Prerequisites & Dependencies
The implementation of comparison methods relies on the `GOAL` R package. Please refer to (https://doi.org/10.1093/bib/bbab331).

Repository Structure
01_code_functions.R: Core function definitions required by subsequent simulation and empirical pipelines.
dhsic.test.R: Auxiliary independence testing functions based on Hilbert-Schmidt Independence Criterion (HSIC).
simulation.R: Script for running parallelized Monte Carlo numerical simulations.
depmap_l1000_causal_data.csv: Preprocessed empirical genomic dataset from the DepMap L1000 benchmark.
empirical.R: Code for real data application and empirical causal inference analysis.

Quick Start & Reproducibility Guide
To reproduce the results or run quick verification tests, execute the scripts in the following order:

If you want to quickly check the pipeline and verify that the algorithm runs correctly on generated data without running full 1,000 Monte Carlo replications:
1. Open `simulation.R`.
2. Modify line 8 to set `n_sim <- 5` (reduces replications to 5).
3. Modify line 159 by replacing `seq_len(nrow(all_config))` with `1` (runs only the first configuration setting).
4. Run `simulation.R`. Execution should complete within a few minutes.

```R
# Line 8 in simulation.R
n_sim <- 5

# Line 159 in simulation.R
for (run_cfg_idx in 1) { ... }
