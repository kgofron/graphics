# WAND2 beamline

### Overview

The folder contains six Python scripts, each contributing to a workflow for simulating, processing, and visualizing diffraction data, likely for a cylindrical detector geometry (as used in neutron or X-ray diffraction experiments). The codebase uses scientific libraries such as NumPy, pandas, matplotlib, and domain-specific libraries (e.g., gi.repository.Hkl).

* * * * *

### File-by-file Analysis

#### 1. detectorposcalc.py

-   Purpose: Provides geometric calculations for mapping real-space angles and positions to detector coordinates, specifically for a cylindrical detector.

-   Key Functions:

-   rotation_matrix: Computes a 3D rotation matrix for a given axis and angle.

-   cylinder_intersection: Finds the intersection of a ray with a cylinder.

-   real2det: Converts sample rotation angles (omega, chi, phi, theta) into unfolded detector coordinates (theta, z) using the above utilities.

#### 2. dfhkl2dfhklaxes.py

-   Purpose: Converts a DataFrame of reflection indices (h, k, l) and intensities into a DataFrame with corresponding instrument axis values (omega, chi, phi, tth) using a geometry engine.

-   Key Functions:

-   dfhkl2dfhklaxes: For each reflection above a minimum intensity, calculates the instrument angles required to access that reflection, using a provided geometry engine and detector/sample objects.

#### 3. hkl2dfhkl.py

-   Purpose: Reads a .hkl file (reflection list) and extracts both lattice parameters and a DataFrame of reflection data.

-   Key Functions:

-   hkl2dfhkl: Parses the file for lattice constants and reflection data, returning both as structured data.

#### 4\. intensities2detint.py

-   Purpose: Orchestrates the conversion from crystallographic data (CIF and HKL files) to a list of detector coordinates and intensities, suitable for visualization.

-   Key Functions:

-   intensities2detint: Runs an external tool to generate HKL data from a CIF file, processes it through the above modules, and returns a list of (theta, z, intensity, h, k, l) tuples for valid detector hits.

#### 5\. interactive_plot.py

-   Purpose: Provides an interactive matplotlib-based visualization of the simulated detector data, allowing users to adjust parameters (e.g., omega window, z offset, intensity threshold) and explore the resulting heatmap.

-   Key Features:

-   Uses sliders and text boxes for interactive parameter adjustment.

-   Displays a heatmap of intensity on the unfolded detector surface.

-   Shows tooltips with reflection indices (hkl) on hover.

#### 6\. plot.py

-   Purpose: Generates a static heatmap visualization of the detector data.

-   Key Features:

-   Bins the simulated detector hits into a 2D grid.

-   Applies Gaussian blur to simulate detector resolution.

-   Plots the result as a heatmap.

* * * * *

### Workflow Summary

1.  Data Preparation:

-   hkl2dfhkl.py reads reflection data and lattice parameters.

-   dfhkl2dfhklaxes.py computes instrument angles for each reflection.

1.  Geometric Mapping:

-   detectorposcalc.py maps these angles to detector coordinates.

1.  Simulation Orchestration:

-   intensities2detint.py ties together the above steps, starting from a CIF file and producing detector hit data.

1.  Visualization:

-   plot.py and interactive_plot.py visualize the simulated detector data, with the latter providing interactive exploration.

* * * * *

### Scientific Context

This codebase is designed for simulating and visualizing diffraction patterns on a cylindrical detector, likely for single-crystal or powder diffraction experiments. It automates the process from crystallographic input files to detector hit maps, supporting both static and interactive analysis.