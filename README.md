# Poisson_CCD22
Poisson solver for LSST CCDs
Description of stand-alone Poisson solver.
Craig Lage - UC Davis - 11-Nov-16

Description: This code is a simple grid-based Poisson's equation solver intended to simulate pixel distortion effects in thick fully-depleted CCD's to be used in the LSST digital camera.  The code builds a 3D rectilinear grid to represent a portion of the CCD, assigns the appropriate charge densities and applied potentials, then solves Poisson's equation using multi-grid methods.  A 360^3 grid, which is adequate for most purposes solves in less than one minute on a typical laptop.  Plotting routines are available to plot the potentials, E-Fields, Pixel shapes, and electron paths.  Installation instructions and instructions on how to run the examples are in the README.pdf file.

Changes: This version is a major revision from the Mar, 2016 version, with several new options and control parameters added.  The major revisions are as follows:

(1) File Formatting: All output files now have extensions of either ".hdf5" or ".dat" to make it clearer
what type of file they are.

(2) Changes to Formatting of Pts files: The "*_Pts.dat", files, which have the information on electron
tracking, have been re-formatted based on input and code from David Kirkby. The files now have a unique
ID for each electron, and phase information of where the electron is in its track.

(3) Non-Linear Z-axis: This option "stretches" the z-axis to put more grid cells near the bottom of the silicon
where things are changing rapidly, and fewer grid cells at the top where things change slowly. This allows
higher resolution without increasing compute time significantly. It is controlled by the parameter "NZExp",
which is the stretch factor at z=0. A value of NZExp = 1.0 reverts to a linear z-axis. A value of NZExp =
10.0 is recommended. This option was also in the "hole7" branch, but was hard coded at a value of NZExp
= 10.0. An adjustment to the Z-axis to more correctly display the difference in dielectric constant between
the gate oxide and the silicon is also included. More details are in Appendix D of the white paper BF_WP_11Nov16.pdf.

(4) Ability to read in and continue a long simulation: The control parameters "Continuation" and
"LastContinuationStep" allow you to read in the mobile charge locations and continue a simulation from
where it left off. This is useful for simuations which have crashed or run out of time.

(5)Option of free holes in the channel stop region: One hypothesis still being explored is that the
channel stop region is not fully depleted, and that there are free holes present in this region. This option
was hard coded into the earlier "hole7" branch, but is now available as an option. It is controlled by
the "UndepletedChannelStop" parameter (0 = no free holes; 1= free holes), and the "Vchannelstop" and
"HoleConvergenceVoltage" parameters. The output file "*_Hole.hdf5" contains the locations of the free
holes. It has been found that the inclusion of free holes in the channel stop region gives better agreement
with measurements of the "brighter-fatter effect".  Rather than run a fixed number of iterations, as in Hole7,
it now iterates until the difference in voltage from Vchannelstop is less than the limit given
in the parameter HoleConvergenceVoltage. The code for generating the free hole density is slower and less robust than
I would like, but seems to run OK in most cases.

(6) Better control of electron tracking: Two parameters, "EquilibrateSteps" and "BottomSteps" have
been added to give better control of the electron tracking. The first (which was previously hard coded
at 100) counts the number of diffusion steps after the electron reaches the collecting well before we start
storing the charge location, and the second (which was previously hard coded at 1000) counts the number
of diffusion steps after the electron reaches the collecting well after we start storing the charge location.

(7) The code now outputs files called *_grid.dat (*=x,y,z) to document the grids used.  This is based
on a suggestion from David Kirkby, and makes the plotting much easier and more correct.  This revision
chnages the format from what David had suggested in order to include the cell boundaries as well as
the cell center.  This is needed to accomodate the non-linear Z-axis discussed above.

Hopefully you find this useful.  Comments and questions are encouraged and should be addressed to: cslage@ucdavis.edu

