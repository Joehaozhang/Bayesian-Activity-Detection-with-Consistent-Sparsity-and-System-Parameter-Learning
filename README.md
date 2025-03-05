# GHVI and MAP Algorithms for JADCE

This repository is the official implementation of [Activity Detection for Massive Connectivity in Cell-Free Networks With Unknown Large-Scale Fading, Channel Statistics, Noise Variance, and Activity Probability: A Bayesian Approach](https://ieeexplore.ieee.org/document/10418889).

## Abstract
Activity detection is an important task in the next generation grant-free multiple access. While there are a number of existing algorithms designed for this purpose, they mostly require precise information about the network, such as large-scale fading coefficients, small-scale fading channel statistics, noise variance at the access points, and user activity probability. Acquiring these information would take a significant overhead and their estimated values might not be accurate. This problem is even more severe in cell-free networks as there are many of these parameters to be acquired. Therefore, this paper sets out to investigate the activity detection problem without the above-mentioned information. In order to handle so many unknown parameters, this paper employs the Bayesian approach, where the unknown variables are endowed with prior distributions which effectively act as regularizations. Together with the likelihood function, a maximum a posteriori (MAP) estimator and a variational inference algorithm are derived. Extensive simulations demonstrate that the proposed methods, even without the knowledge of these system parameters, perform better than existing state-of-the-art methods, such as covariance-based and approximate message passing methods.

## Function Description

### main
Main program for activity detection in cell-free systems

### MAP_GH_cellfree
MAP algorithm for activity status and system parameter estimation in cell-free systems

### VIAD_GH_cellfree
GHVI algorithm for activity status and system parameter estimation in cell-free systems

### DominantAPSelection
Script to select dominant AP for each device

### PFAPMD
Generate PFA and PMD by selecting different thresholds for hyper parameter z

### PFAPMDNMSE_cellfree
Generate PFA, PMD, and NMSE by selecting different threshold for dominant channel energy

## Copyright
If you use this code, please cite [our published TSP article](https://ieeexplore.ieee.org/document/10418889).