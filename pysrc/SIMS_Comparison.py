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

# This holds all of the data
dat = Array3dHDF5(outputfiledir, outputfilebase, run)

ScaleFactor = ConfigData["ScaleFactor"]
GridsPerPixelX = ConfigData["GridsPerPixelX"]
GridsPerPixelY = ConfigData["GridsPerPixelY"]
ChargeFactor = 1.6E-19 * 1.0E6 / (11.7 * 8.85E-12)/((dat.x[1]-dat.x[0])*(dat.y[1]-dat.y[0])) #(QE*MICRON_PER_M/(EPSILON_0*EPSILON_SI))/(dx*dy)


[boron_depth, boron_conc, phos_depth, phos_conc] = ReadSIMSData()

nxx = dat.nx - 1
nyy = dat.ny - 1
nzz = dat.nz - 1

nxcenter = nxx/2
nycenter = nyy/2
nxcenter2 = nxcenter

NumPixelsPlotted = 4
nycenter2 = nycenter
nymin = nycenter - (NumPixelsPlotted * ScaleFactor * GridsPerPixelY)/2
nymax = nycenter + (NumPixelsPlotted * ScaleFactor * GridsPerPixelY)/2

nxmin = nxcenter - (NumPixelsPlotted * ScaleFactor * GridsPerPixelX)/2
nxmax = nxcenter + (NumPixelsPlotted * ScaleFactor * GridsPerPixelX)/2

nzmin = 0
nzmax = 16 * ScaleFactor

# Create the output directory if it doesn't exist
if not os.path.isdir(outputfiledir+"/plots"):
    os.mkdir(outputfiledir+"/plots")

rcParams['contour.negative_linestyle'] = 'solid'
rcParams.update({'font.size': 6})


print "Making SIMS Comparison plots\n"
figure()

suptitle("Charge Density Comparison. Grid = %d*%d*%d."%(nxx,nyy,nzz),fontsize = 18)
plotcounter = 1
subplots_adjust(hspace=0.3, wspace=0.3)
numzs = 160

subplot(1,2,1)
title("Rho-Collect Gate", fontsize=12)

# Determine location of bottom of silicon
for nz in range(dat.nz):
    if abs(dat.rho[nxcenter2, nycenter2, nz]) > 1.0E-12:
        break

Chan_z0 = dat.z[nz]

plot(phos_depth + Chan_z0, phos_conc, label = "SIMS Profiles", color = 'red', zorder = -1)
scatter(dat.z[0:numzs], (dat.rho[nxcenter2,nycenter2,0:numzs]+dat.rho[nxcenter2-1,nycenter2,0:numzs]+dat.rho[nxcenter2,nycenter2-1,0:numzs]+dat.rho[nxcenter2-1,nycenter2-1,0:numzs])/4.0, label = "Poisson Model", color='green', zorder = 1)

legend(loc = "upper right")
xlabel("Z-Dimension (microns)", fontsize=12)
ylabel('Charge Density ($V/ \mu m^2$)',fontsize=12)
ylim(-20.0, 100.0)
xlim(0.0,2.0)
nxcenter3 = nxcenter2 + 3 * GridsPerPixelX * ScaleFactor / 2
nycenter3 = nycenter2
nycenter4 = nycenter2 + GridsPerPixelY * ScaleFactor / 2

subplot(1,2,2)
title("Rho-ChanStop", fontsize=12)

# Determine location of bottom of silicon
for nz in range(dat.nz):
    if abs(dat.rho[nxcenter3, nycenter3, nz]) > 1.0E-12:
        break

ChanStop_z0 = dat.z[nz]

plot(boron_depth + ChanStop_z0, boron_conc, label = "SIMS Profiles", color = 'red', zorder = -1)
scatter(dat.z[0:numzs], ((dat.rho[nxcenter3,nycenter3,0:numzs]+dat.rho[nxcenter3-1,nycenter3,0:numzs]+dat.rho[nxcenter3,nycenter3-1,0:numzs]+dat.rho[nxcenter3-1,nycenter3-1,0:numzs])/4.0), color = 'green', label = "Poisson Model", zorder = 1)


legend(loc = "upper right")
xlabel("Z-Dimension (microns)", fontsize=12)
ylim(-100.0, 20.0)
xlim(0.0,4.0)
savefig(outputfiledir+"/plots/"+outputfilebase+"_%d_SIMS_Comparison_1.pdf"%run)


figure()

suptitle("Charge Density Comparison. Grid = %d*%d*%d."%(nxx,nyy,nzz),fontsize = 18)
plotcounter = 1
subplots_adjust(hspace=0.3, wspace=0.3)
numzs = 160

subplot(1,2,1)
title("Rho-Collect Gate", fontsize=12)

# Determine location of bottom of silicon
for nz in range(dat.nz):
    if abs(dat.rho[nxcenter2, nycenter2, nz]) > 1.0E-12:
        break

Chan_z0 = dat.z[nz]

plot(phos_depth + Chan_z0, phos_conc, label = "SIMS Profiles", color = 'red', zorder = -1)
scatter(dat.z[0:numzs], (dat.rho[nxcenter2,nycenter2,0:numzs]+dat.rho[nxcenter2-1,nycenter2,0:numzs]+dat.rho[nxcenter2,nycenter2-1,0:numzs]+dat.rho[nxcenter2-1,nycenter2-1,0:numzs])/4.0, label = "Poisson Model", color='green', zorder = 1)

legend(loc = "upper right")
xlabel("Z-Dimension (microns)", fontsize=12)
ylabel('Charge Density ($V/ \mu m^2$)',fontsize=12)
ylim(-20.0, 1500.0)
xlim(0.0,2.0)
nxcenter3 = nxcenter2 + 3 * GridsPerPixelX * ScaleFactor / 2
nycenter3 = nycenter2
nycenter4 = nycenter2 + GridsPerPixelY * ScaleFactor / 2

subplot(1,2,2)
title("Rho-ChanStop", fontsize=12)

# Determine location of bottom of silicon
for nz in range(dat.nz):
    if abs(dat.rho[nxcenter3, nycenter3, nz]) > 1.0E-12:
        break

ChanStop_z0 = dat.z[nz]

plot(boron_depth + ChanStop_z0, boron_conc, label = "SIMS Profiles", color = 'red', zorder = -1)
scatter(dat.z[0:numzs], ((dat.rho[nxcenter3,nycenter3,0:numzs]+dat.rho[nxcenter3-1,nycenter3,0:numzs]+dat.rho[nxcenter3,nycenter3-1,0:numzs]+dat.rho[nxcenter3-1,nycenter3-1,0:numzs])/4.0), color = 'green', label = "Poisson Model", zorder = 1)


legend(loc = "upper right")
xlabel("Z-Dimension (microns)", fontsize=12)
ylim(-200.0, 20.0)
xlim(0.0,4.0)
savefig(outputfiledir+"/plots/"+outputfilebase+"_%d_SIMS_Comparison_2.pdf"%run)



