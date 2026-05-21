<!--
***************************************************************
*** 3D SPoDEA programe that includes a set of *.ipynb files to compute basement depth of the complex 3D sedimentary basin.
*** Source Code is mainly written for research purposes. The codes are
*** having copyrights and required proper citations whenever it is used.
*** Originated by:
***       Prof. Chandra Prakash Dubey (email:p.dubey48@gmail.com)
***       Mr. Anurag Mondal
***       Geology and Geophysics Department
***       IIT Kharagpur
***       West Bengal, India
****************************************************************
-->

# 3D SPoDEA — Gravity Inversion for 2D and 3D Sedimentary Basins

This is a help file for a description of all Data, Source Code, and Subroutine used for the implementation of our present paper **'3D Gravity Inversion for Heterogeneous Sedimentary Basin with cubic b-Spline Polynomial Approximation using Differential Evolution Algorithm.'**

> *(Copy all set of files including data in one folder)*

---

## 1. Subroutines

| File | Description |
|------|-------------|
| `gravity3d.py` | 3D gravity forward modeling (constant density) |
| `gravity3d_variable_density.py` | 3D gravity forward modeling (variable density) |
| `guass_fft_general3d.py` | Gauss-FFT based 3D gravity computation |

---

## 2. Source Codes

| Notebook | Description |
|----------|-------------|
| `model_1_constant_density_single_basin_pca.ipynb` | Calculates the inversion of gravity anomaly for a synthetic single sedimentary basin having fixed density contrast with and without noise case (Model 1) |
| `model_2_constant_density_twin_basin_pca.ipynb` | Calculates the inversion of gravity anomaly for a synthetic twin sedimentary basin having fixed density contrast with and without noise case (Model 2) |
| `model_3_variable_density_single_basin_pca.ipynb` | Calculates the inversion of gravity anomaly for a synthetic single sedimentary basin having varying density contrast with and without noise case (Model 3) |
| `model_4_variable_density_twin_basin_pca.ipynb` | Calculates the inversion of gravity anomaly for a synthetic twin sedimentary basin having varying density contrast with and without noise case (Model 4) |
| `noise_and_control_point_sensitivity_gravity_inversion_2D.ipynb` | Noise sensitivity and control-point sensitivity analysis for 2D gravity inversion |

---

## Citation

If you use this code in your research, please cite the corresponding paper and give proper attribution to the authors.

**Authors:**
- Prof. Chandra Prakash Dubey — [p.dubey48@gmail.com](mailto:p.dubey48@gmail.com)
- Mr. Anurag Mondal
- Geology and Geophysics Department, IIT Kharagpur, West Bengal, India
