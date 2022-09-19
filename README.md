# NPRE598CP

## Project 1. Ion motion in a magnetic dipole (Earth’s magnetic field)

* Create a routine that calculates the Earth’s magnetic field (assumed to be a dipole field without solar wind interaction). Calculate in spherical coordinates (axial symmetry), then express in GSM coordinates for convenience: `Bx, By, Bz = bfield.dipoleEarth(x, y, z)`
* Sample the magnetic field in the r-theta-plane, and create a function that interpolates the vector field at a generic point P(x,y,z) for computational efficiency.
* Compare computational time with and without interpolation.
* Solve the ion trajectory (Lorentz equation) using the interpolator; initialize ions with Maxwellian distribution at the equatorial plane, and integrate the trajectory over time. Plot in 3D (plotly).
* (Optional) Replace the distribution with a more realistic model.
* Do a Monte Carlo simulation to visualize trapped ions and ions that escape. Compare pathways of escaped ions (e-, H+, N+, O+).
* Create a function that calculates the criteria of trapping a particle in a magnetic dipole.
* Compare the simulation with the calculation.
