# Poisson_CCD22 - 'hole17' branch
Poisson solver for LSST CCDs
Description of stand-alone Poisson solver.
Craig Lage - UC Davis - 3-Feb-17

Description: This code is a simple grid-based Poisson's equation solver intended to simulate pixel distortion effects in thick fully-depleted CCD's to be used in the LSST digital camera.  The code builds a 3D rectilinear grid to represent a portion of the CCD, assigns the appropriate charge densities and applied potentials, then solves Poisson's equation using multi-grid methods.  A 360^3 grid, which is adequate for most purposes solves in less than one minute on a typical laptop.  Plotting routines are available to plot the potentials, E-Fields, Pixel shapes, and electron paths.  Installation instructions and instructions on how to run the examples are in the README.pdf file.


Changes:

This 'hole17' branch is a major revision. It is still under development and may contain bugs. However, it was
used to create all of the plots in the document docs/PACCD_Paper_2Feb17.pdf, so it is reasonably mature.

Mobile Charges: The treatment of mobile charges, both holes and electrons, is the biggest change in the
version. Mobile charges are now dealt with using the Quasi-Fermi level formalism. This self consistently
solves for the potential and mobile charge densities in regions containing mobile charges. The hole-containing
regions in these devices are a spatially continuous region which are in quasi-equilibrium, so a single hole
Quasi-Fermi level serves to set the holes densities throughout the device. This is set in the .cfg file by the
parameter qfh. Setting this parameter to a value of 0.4 Volts will set the potential in the hole containing
region to be near 0 Volts, which is appropriate for the ITL devices. For the electrons, a different value of
the electron Quasi-Fermi level is set in each pixel based on the number of electrons in this pixel. Because
the relation between the electron Quasi-Fermi level and the number of electrons in the pixel is a compelx
non-linear relationship, the code builds a look-up table for this. Building this look-up table is controlled
by the BuildQFe, QFeMin, QFeMax, and NQFe parameters. This only need to be done each time the
CCD parameters (voltages, dopings, etc.) are changed. If BuildQFeLookup = 0, the code will look for the
*QFe.dat file to read in.

Geometry of gate and field oxide regions: The code now has a more realistic gate and field oxide
geometry, controlled by the GateOxide, FieldOxide, and FieldOxideTaper parameters.
4 Examples


Hopefully you find this useful.  Comments and questions are encouraged and should be addressed to: cslage@ucdavis.edu

