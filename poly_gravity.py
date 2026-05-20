from numba import njit
import numpy as np

@njit
def poly_gravity(x_obs, z_obs, x, z, roh, t, c):
    G = 6.67408e-11  # Gravitational constant
    n_poly = len(x)
    x = np.append(x, x[0])
    z = np.append(z, z[0])
    grav = np.zeros_like(x_obs)

    for i in range(len(x_obs)):
        value = 0.0
        for j in range(n_poly):
            for k in range(len(t)):
                x_t = x[j] * (1 - t[k]) + x[j + 1] * t[k]
                z_t = z[j] * (1 - t[k]) + z[j + 1] * t[k]
                ax1 = x_t - x_obs[i]
                ax2 = z_t - z_obs
                integrand = -2 * roh * G * np.arctan2(ax1, ax2) * (z[j + 1] - z[j])
                value += c[k] * integrand
        grav[i] = 1e5 * value  # convert to mGal
    return grav
