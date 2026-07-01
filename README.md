# DMADFO

[![DOI](https://doi.org/10.5281/zenodo.21104547.svg)](https://doi.org/10.5281/zenodo.21104547)
[![FWF Grant DOI](https://img.shields.io/badge/FWF%20Grant-10.55776%2FPAT2747625-blue)](https://doi.org/10.55776/PAT2747625)

**DMADFO** performs approximate unconstrained noisy derivative-free minimization of a possibly nonsmooth function of many continuous variables.

The method does **not** require gradients or subgradients, and noisy function evaluations are allowed.

DMADFO is a large-scale diagonal extension of **MADFO** designed for noisy derivative-free optimization problems. It replaces the full affine-scaling matrix used in MADFO with a diagonal conjugacy-based scaling mechanism, while preserving the main components of the original MADFO framework.

The solver is intended for high-dimensional optimization problems where derivative information is unavailable, unreliable, or too expensive to compute.

## Contents

This repository contains:

* MATLAB source code for **DMADFO Version 1.0**;
* a driver program, `driver.m`, showing how to use the solver;
* the accompanying paper describing the implemented method and numerical test results.

## Usage

To run the example driver in MATLAB, use:

```matlab
driver
```

The file `driver.m` provides a basic example of how to call the solver.

## Author

The program was written by:

**Morteza Kimiaei**
Faculty of Mathematics
University of Vienna

Please inform the author at `morteza.kimiaei@univie.ac.at` if you make serious use of this code.

## DOI and funding

The archived software release is available on Zenodo:

[https://doi.org/10.5281/zenodo.21104547](https://doi.org/10.5281/zenodo.21104547)

This project was funded by the **Austrian Science Fund (FWF)** under:

[https://doi.org/10.55776/PAT2747625](https://doi.org/10.55776/PAT2747625)

## References

### MADFO

M. Kimiaei and A. Neumaier,
*Effective matrix adaptation strategy for noisy derivative-free optimization*,
Mathematical Programming Computation 16, 459--501, 2024.
https://doi.org/10.1007/s12532-024-00261-z

### DMADFO

M. Kimiaei and S. Babaie-Kafaki,
*Diagonal Hessian Approximation Based on Conjugacy Condition for Noisy Derivative-Free Optimization Problems in High Dimensions*,
arXiv:2606.20304, 2026.
https://arxiv.org/abs/2606.20304

This paper describes the implemented method and reports numerical test results.

## Citation

Please cite the DMADFO paper when using this package in scientific work:

```bibtex
@article{kimiaei_babaiekafaki_2026_dmadfo,
  title   = {Diagonal Hessian Approximation Based on Conjugacy Condition for Noisy Derivative-Free Optimization Problems in High Dimensions},
  author  = {Kimiaei, Morteza and Babaie-Kafaki, Saman},
  journal = {arXiv preprint arXiv:2606.20304},
  year    = {2026},
  url     = {https://arxiv.org/abs/2606.20304}
}
```
