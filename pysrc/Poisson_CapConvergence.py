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
Nzelec = ConfigData["Nzelec"]

# Create the output directory if it doesn't exist
if not os.path.isdir(outputfiledir+"/plots"):
    os.mkdir(outputfiledir+"/plots")

if ConfigData["SimulationRegionLowerLeft"][0] < 50.0:
    print "Field Oxide"
    DeltaPhi_calc = ConfigData["ChannelStopSurfaceCharge"] * ConfigData["FieldOxide"] * 1.0E-4 * 1.6E-19 / (8.85E-14 * 4.3)
else:
    print "Gate Oxide"
    DeltaPhi_calc = ConfigData["ChannelSurfaceCharge"] * ConfigData["GateOxide"] * 1.0E-4 * 1.6E-19 / (8.85E-14 * 4.3)    

print "Making 1D Potential and Charge Density plots\n"
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
    zm = zmax / pow(2,i)
    for k in range(zm):
        if dat.rho[nxcenter,nycenter,k] > 0.1:
            # Finds where the charge is
            break
    DeltaPhi = dat.phi[nxcenter, nycenter,k] - dat.phi[nxcenter,nycenter,0]
    print "For i = %d, Nz = %d, DeltaPhi = %f"%(i,nzz,DeltaPhi)
    subplot(1,2,1)
    plot(dat.z, dat.phi[nxcenter, nycenter,:], label = "Grid %d"%i, color = colors[i], marker = 'x')
    ylim(0.0, 20.0)
    xlim(0.0,2.0)
    subplot(1,2,2)
    plot(dat.z[0:zm], dat.phi[nxcenter, nycenter,0:zm], label = "Grid %d"%i, color = colors[i], marker = 'x')
    text(2.0, -5.0-1.0*i, "i = %d, Nz = %d, DeltaPhi = %f"%(i,nzz,DeltaPhi), fontsize = 8)
    ylim(-15.0, 20.0)
    
subplot(1,2,2)
text(2.0, -5.0-1.0*NumMultis, "Calculated DeltaPhi = %f"%DeltaPhi_calc, fontsize = 8)
print "Calculated DeltaPhi = %f"%DeltaPhi_calc
subplot(1,2,1)
legend(loc = 'upper left')
savefig(outputfiledir+"/plots/"+outputfilebase+"_Grid_Comparison.pdf")

