# GHVI and MAP Algorithms for JADCE

This repository is the official implementation of [Activity Detection for Massive Connectivity in Cell-Free Networks With Unknown Large-Scale Fading, Channel Statistics, Noise Variance, and Activity Probability: A Bayesian Approach](https://ieeexplore.ieee.org/document/10418889).

## Function Description

### main.m function
Main program for activity detection in cell-free systems

### MAP_GH_cellfree.m function
MAP algorithm for activity status and system parameter estimation in cell-free systems

### VIAD_GH_cellfree.m function
GHVI algorithm for activity status and system parameter estimation in cell-free systems

### DominantAPSelection.m function
Script to select dominant AP for each device

### PFAPMD.m function
Generate PFA and PMD by selecting different thresholds for hyper parameter z

### PFAPMDNMSE_cellfree.m function
Generate PFA, PMD, and NMSE by selecting different threshold for dominant channel energy

## Copyright
If you use this code, please cite [our published TSP article](https://ieeexplore.ieee.org/document/10418889).