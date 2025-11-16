import numpy as np
import matplotlib.pyplot as plt

# -----------------------------
# Dipole Antenna Parameters
# -----------------------------
frequency = 1e9          # 1 GHz
c = 3e8                  # Speed of light
wavelength = c / frequency
k = 2 * np.pi / wavelength  # Wave number

L = wavelength / 2       # Dipole length (half-wavelength)

theta = np.linspace(0, np.pi, 500)  # Angle from z-axis

# -----------------------------
# Radiation Pattern (E-theta)
# -----------------------------
# Formula for thin half-wave dipole:
E_theta = np.sin((k * L / 2) * np.cos(theta)) / np.sin(theta)
E_theta = np.nan_to_num(E_theta)  # Replace NaNs (0/0) with 0

# Normalize
E_theta = E_theta / np.max(np.abs(E_theta))

# -----------------------------
# Plot 2D Radiation Pattern
# -----------------------------
plt.figure(figsize=(6,6))
plt.polar(theta, np.abs(E_theta))
plt.title("Dipole Antenna Radiation Pattern (2D)")
plt.show()

# -----------------------------
# Plot 3D Radiation Pattern
# -----------------------------
from mpl_toolkits.mplot3d import Axes3D

phi = np.linspace(0, 2*np.pi, 100)
theta_grid, phi_grid = np.meshgrid(theta, phi)

# Radiation magnitude in 3D
E_3D = np.sin((k * L / 2) * np.cos(theta_grid)) / np.sin(theta_grid)
E_3D = np.nan_to_num(E_3D)
E_3D = E_3D / np.max(np.abs(E_3D))

# Convert to Cartesian for 3D plotting
X = E_3D * np.sin(theta_grid) * np.cos(phi_grid)
Y = E_3D * np.sin(theta_grid) * np.sin(phi_grid)
Z = E_3D * np.cos(theta_grid)

fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, Z, facecolors=plt.cm.viridis(E_3D), rstride=2, cstride=2, alpha=0.8)
ax.set_title("Dipole Antenna Radiation Pattern (3D)")
plt.show()
