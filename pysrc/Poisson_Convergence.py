#!/usr/bin/env python

#Author: Craig Lage, NYU;
#Date: 10-Sep17

#This program plots the Poisson equation solutions from the C++ Poisson solver
import matplotlib
matplotlib.use("PDF")
from pylab import *
import os, sys, time, h5py
sys.path.append(os.path.realpath('./pysrc'))
from pysubs import *  # These are the plotting subroutines

#****************MAIN PROGRAM*****************

# First, read the .cfg file

configfile = sys.argv[1]
run = int(sys.argv[2])
ConfigData = ReadConfigFile(configfile)

outputfilebase = ConfigData["outputfilebase"]
outputfiledir = ConfigData["outputfiledir"]
ScaleFactor = ConfigData["ScaleFactor"]
GridsPerPixelX = ConfigData["GridsPerPixelX"]
# Create the output directory if it doesn't exist
if not os.path.isdir(outputfiledir+"/plots"):
    os.mkdir(outputfiledir+"/plots")


print "Making Potential and Charge Density Convergence plots\n"
figure()

suptitle("Convergence Tests")
subplot(1,2,1)
title("Phi-Collect Gate", fontsize=12)
subplot(1,2,2)
title("Phi-Collect Gate", fontsize=12)
zmax = 64 * ScaleFactor
colors = ['black','red','green','blue','magenta','cyan']

NumMultis = int(3 + log2(ScaleFactor))

for i in range(NumMultis):

    # This holds all of the data
    dat = Array3dHDF5(outputfiledir, outputfilebase, run, multi=i)
    nxx = dat.nx - 1
    nyy = dat.ny - 1
    nzz = dat.nz - 1
    nxcenter = nxx/2
    nycenter = nyy/2
    subplot(2,2,1)
    plot(dat.z, dat.phi[nxcenter, nycenter,:], label = "Grid %d"%i, color = colors[i], marker = 'x')
    ylim(0.0, 25.0)
    xlim(0.0,5.0)
    subplot(2,2,2)
    zm = zmax / pow(2,i)
    plot(dat.z[0:zm], dat.phi[nxcenter, nycenter,0:zm], label = "Grid %d"%i, color = colors[i], marker = 'x')
    ylim(-15.0, 25.0)
    nxcenter = nxcenter + ScaleFactor * GridsPerPixelX / 2 / pow(2,i)
    subplot(2,2,3)
    plot(dat.z, dat.phi[nxcenter, nycenter,:], label = "Grid %d"%i, color = colors[i], marker = 'x')
    ylim(-60.0, 20.0)
    subplot(2,2,4)
    plot(dat.z[0:zm], dat.phi[nxcenter, nycenter,0:zm], label = "Grid %d"%i, color = colors[i], marker = 'x')
    ylim(-6.0, 9.0)
subplot(2,2,3)
legend(loc = 'upper right')
savefig(outputfiledir+"/plots/"+outputfilebase+"_Grid_Comparison.pdf")


print "Making 1D Electron Density plots\n"
figure()

suptitle("Convergence Tests")
subplot(1,1,1)
title("Electron Density - Collect Gate", fontsize=12)
colors = ['black','red','green','blue','magenta','cyan']
for i in range(NumMultis):

    # This holds all of the data
    dat = Array3dHDF5(outputfiledir, outputfilebase, run, multi=i)
    nxx = dat.nx - 1
    nyy = dat.ny - 1
    nzz = dat.nz - 1
    zmax = dat.Elec.shape[2]
    nxcenter = nxx/2
    nycenter = nyy/2
    plot(dat.z[0:zmax], dat.Elec[nxcenter, nycenter,:] / (dat.dx[0] * dat.dy[0]) / dat.dz[0:zmax], label = "Grid %d"%i, color = colors[i], marker = '*')
    #ylim(-60.0, 20.0)
legend(loc = 'upper right')
savefig(outputfiledir+"/plots/"+outputfilebase+"_Elec_Comparison.pdf")
