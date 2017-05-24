# Poisson_CCD22 - 'hole18' branch
Poisson solver for LSST CCDs
Description of stand-alone Poisson solver.
Craig Lage - UC Davis - 24-May-17

Description: This code is a simple grid-based Poisson's equation solver intended to simulate pixel distortion effects in thick fully-depleted CCD's to be used in the LSST digital camera.  The code builds a 3D rectilinear grid to represent a portion of the CCD, assigns the appropriate charge densities and applied potentials, then solves Poisson's equation using multi-grid methods.  A 360^3 grid, which is adequate for most purposes solves in less than one minute on a typical laptop.  Plotting routines are available to plot the potentials, E-Fields, Pixel shapes, and electron paths.  Installation instructions and instructions on how to run the examples are in the README.pdf file.


Changes:

This 'hole18' branch has a number of changes from the 'hole17' branch.  It has not been tested as thoroughly as I would like.  however, it has run with all of the .cfg files in the data directory and produced the plots found there.  It has been modified to include the results of physical measurements (cross sectional SEM's, threshold voltage measurements, and SIMS measurements of dopant profiles).  With these modifications, it has produced the most accurate simulations to date of the ITL CCD detector.  Some of the changes include:

(1) To better fit the measured dopant profiles, it now supports multiple Gaussians for the profiles, as well as a surface charge in the channel region.

(2) In now supports non-square pixel shapes.  This may cause some problems with running legacy plotting files, since parameters like PixelSize are now replaced with PixelSizeX and PixelSizeY.

(3) There have been a lot of changes to support a list of VoltageRegions.  This allos simulation of regions like the I/Os and the Serial regions, which are geometrically more complicated than the regular imaging array.

(4) Changes have been made to the Multi-Grid strategy to make the coarser grids mimic the finest grid as closely as possible.  this helps speed convergence.

Data files:

Five different .cfg files are supplied and have been tested.  The .cfg files and the resultsing plots are in the data directory, but not the large data output files.  The files are:

(0) data/run0/bf.cfg - This is a 'quick and dirty' file that generates a 9x9 grid of pixels with the center pixle having 80,000 electrons, at lower resolution.  It should run in 10-15 minutes.

(1) data/run1/bf.cfg - This is the same as run0, but at higher resolution.  It will take several hours to run.

(2) data/run2/bf.cfg - Similar to run1, but a 13x13 pixel array, and calculates the pixel shapes, which is slow.  The file data/run2/BF_0_Vertices.dat is the most accurate simulation I have of the distorted pixel shapes, and produced the plot data/run2/plots/Area_Corr_19May17_0.pdf, which is the best fit I have to the correlation measurements.

(3) data/run3/trans.cfg - this is a simulation of the output transistor, and was used to simulate the threshold volatge of this device for comparison to measurements.

(4) data/run4/io.cfg - this is a simulation of the entire IO region, including the output transistor, the reset gate, and the end of the serial chain.  It is still being developed, but seems to give sensible results.

(5) data/run5/sat.cfg - this is what was used to predict the onset of blooming in the array as a function of the parallel gate voltages.

Other comments:

I have not updated the README.pdf file for this version, so these notes are the only documentation on this version of the code.


Hopefully you find this useful.  Comments and questions are encouraged and should be addressed to: cslage@ucdavis.edu

