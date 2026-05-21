# Gravity Inversion for 2D and 3D Sedimentary Basins

This is a help file for a description of all Data, Source Code, and Subroutine used for the implementation of our present paper **'3D Gravity Inversion for Heterogeneous Sedimentary Basin with cubic b-Spline Polynomial Approximation using Differential Evolution Algorithm.'**


---

## 1. Subroutines

| File | Description |
|------|-------------|
| `poly_gravity.py` | 2D gravity forward modeling (constant density) |
| `gravity3d1.py` | 3D gravity forward modeling (constant density) |
| `gravity3d_variable_density1.py` | 3D gravity forward modeling (variable density) |

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
- Mr. Rajdeep Das — [phmgleb@gmail.com](mailto:phmgleb@gmail.com)
- Department of Applied Geophysics, IIT (ISM) Dhanbad, Jharkhand, India
