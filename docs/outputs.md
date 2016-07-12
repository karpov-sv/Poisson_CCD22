# Output Files

This page documents the output files produced by the Poisson program and read by the various python scripts.

All output files are either in HDF5 format (with extension `.hdf5`) or plain text (with extension `.dat`).

## Charge density

File names ending with `_rho.hdf5` contain a 3D array of charge density values expressed `rho/eps` in units of Volts/(micron^2), where `eps` is the permittivity of silicon.

## Electric potential

File names ending with `_phi.hdf5` contain a 3D array of electric potential values in units of Volts.

## Electric field

File names ending with `_Ex.hdf5`, `_Ey.hdf5`, or `_Ez.hdf5` contain 3D arrays of electric field component (x,y,z) values in units of Volts/micron.

## Mobile charges

File names ending with `_Elec.hdf5` or `_Hole.hdf5` contain 3D arrays of mobile electron and hole number densities in units of 1/(micron^3).
