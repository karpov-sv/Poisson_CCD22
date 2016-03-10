# Poisson_CCD22
Poisson solver for LSST CCDs
Description of stand-alone Poisson solver.
Craig Lage - UC Davis - 10-Mar-16

Description: This code is a simple grid-based Poisson's equation solver intended to simulate pixel distortion effects in thick fully-depleted CCD's to be used in the LSST digital camera.  The code builds a 3D rectilinear grid to represent a portion of the CCD, assigns the appropriate charge densities and applied potentials, then solves Poisson's equation using multi-grid methods.  A 360^3 grid, which is adequate for most purposes solves in less than one minute on a typical laptop.  Plotting routines are available to plot the potentials, E-Fields, Pixel shapes, and electron paths.  Installation instructions and instructions on hom to run the examples are in the README.pdf file.

Hopefully you find this useful.  Comments and questions are encouraged and should be addressed to: cslage@ucdavis.edu

