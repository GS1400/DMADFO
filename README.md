DMADFO performs approximate unconstrained noisy derivative-free minimization
of a not necessarily smooth function of many continuous variables. No gradients
or subgradients are required, and noisy function evaluations are allowed.

DMADFO is a large-scale diagonal extension of MADFO designed for noisy
derivative-free optimization problems. It replaces the full affine-scaling
matrix used in MADFO by a diagonal conjugacy-based scaling mechanism, while
preserving the main components of the original MADFO framework.

The program was written by Morteza Kimiaei (University of Vienna).
Please inform the author at morteza.kimiaei@univie.ac.at if you make serious use
of this code.

This directory contains:
* the MATLAB source code for Version 1.0 (June 18, 2026),
* a driver program, driver.m, showing its use,
* and the accompanying paper.

Reference for MADFO:

M. Kimiaei and A. Neumaier,
Effective matrix adaptation strategy for noisy derivative-free optimization.
Mathematical Programming Computation 16, 459--501 (2024).
https://doi.org/10.1007/s12532-024-00261-z

Reference for DMADFO:

M. Kimiaei and S. Babaie-Kafaki,
Diagonal Hessian Approximation Based on Conjugacy Condition for Noisy
Derivative-Free Optimization Problems in High Dimensions.

This paper describes the implemented method and reports numerical test results.
Please cite this paper when using this package in scientific work.

