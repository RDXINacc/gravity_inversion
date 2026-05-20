# gravity_inversion
# Gravity Inversion for 2D and 3D Sedimentary Basins

This repository contains codes and example workflows for 2D and 3D gravity forward modeling and inversion of synthetic basins and a real sedimentary basin. The implementations focus on basin geometry and density reconstruction using B-spline parameterization, Differential Evolution, and gradient-based refinement.

## Features

- 2D forward gravity modeling for depth-dependent density profiles
- 3D forward modeling on regular grids for synthetic and real basins
- Nonlinear gravity inversion with:
  - Differential Evolution (global search)
  - L-BFGS-B (local refinement)
- Constant-density and vertically varying density basin models
- Noise sensitivity and control-point sensitivity analysis
- Cost-function topography analysis in a 2D PCA subspace
- Real basin modeling and inversion using field gravity data
