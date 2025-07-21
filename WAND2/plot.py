#TODO Requires filter for when Bragg condition is met

import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
from intensities2detint import intensities2detint

cif_path = '/epics/crystals/SiO2/EntryWithCollCode176.cif'
hkl_path = '/epics/crystals/SiO2/EntryWithCollCode176.hkl'

R = 0.7  # m
wavelength = 1.486  # A
width = 120  # degree
cyl_center = (0,0)
ray_origin = np.array([0.0, 0.0, 0.0])
min_intensity = 100

angle_deg = 7.5
angle_rad = np.deg2rad(angle_deg)
zmax = R * np.tan(angle_rad)
zmin = -zmax

gamma_axis = [0,0,1]
delta_axis = [0,-1,0]

lst = intensities2detint(
    cif_path, hkl_path, wavelength, min_intensity,
    R, cyl_center, ray_origin, zmin, zmax, gamma_axis, delta_axis
)
print(lst)

data = np.array(lst)
theta = data[:,0]
z = data[:,1]
intensity = data[:,2]

# Create heatmap grid
nx, ny = 500, 200  # resolution of the heatmap
theta_grid = np.linspace(0, 360, nx)
z_grid = np.linspace(zmin, zmax, ny)
X, Y = np.meshgrid(theta_grid, z_grid)

heatmap = np.zeros_like(X)

# Bin the data into the grid
for t, zz, inten in zip(theta, z, intensity):
    i = np.searchsorted(theta_grid, t)
    j = np.searchsorted(z_grid, zz)
    if 0 <= i < nx and 0 <= j < ny:
        heatmap[j, i] += inten

# Gaussian blur to create smooth halos
heatmap_blurred = gaussian_filter(heatmap, sigma=3)

# Plot heatmap
plt.figure(figsize=(10,4))
plt.imshow(
    heatmap_blurred,
    extent=[0,360, zmin, zmax],
    origin='lower',
    aspect='auto',
    cmap='viridis'
)
plt.colorbar(label='Intensity')
plt.xlabel("θ (degrees)")
plt.ylabel("z")
plt.title("Unfolded Cylinder Heatmap with Gaussian Halos")
plt.tight_layout()
plt.show()

