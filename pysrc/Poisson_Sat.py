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


nxx = dat.nx - 1
nyy = dat.ny - 1
nzz = dat.nz - 1

nxcenter = nxx/2
nycenter = nyy/2
nxcenter2 = nxcenter

NumPixelsPlottedX = 3
NumPixelsPlottedY = 4    
nycenter2 = nycenter
nymin = nycenter - (NumPixelsPlottedY * ScaleFactor * GridsPerPixelY)/2
nymax = nycenter + (NumPixelsPlottedY * ScaleFactor * GridsPerPixelY)/2

nxmin = nxcenter - (NumPixelsPlottedX * ScaleFactor * GridsPerPixelX)/2
nxmax = nxcenter + (NumPixelsPlottedX * ScaleFactor * GridsPerPixelX)/2

nzmin = 0
nzmax = 16 * ScaleFactor
"""
print "Channel CG"
for i in range(dat.Elec.shape[2]):
    print "k = %d, z=%.3f, phi = %.3f, holes = %.3f, elec = %.3f, rho = %.3f, Ez = %.3f"%(i,dat.z[i],dat.phi[nxcenter, nycenter,i],dat.Hole[nxcenter, nycenter,i],dat.Elec[nxcenter, nycenter,i],dat.rho[nxcenter, nycenter,i],dat.Ez[nxcenter, nycenter,i])

print "Channel BG"
for i in range(dat.Elec.shape[2]):
    print "k = %d, z=%.3f, phi = %.3f, holes = %.3f, elec = %.3f, rho = %.3f, Ez = %.3f"%(i,dat.z[i],dat.phi[nxcenter, nycenter+ScaleFactor*GridsPerPixelY/2,i],dat.Hole[nxcenter, nycenter+ScaleFactor*GridsPerPixelY/2,i],dat.Elec[nxcenter, nycenter+ScaleFactor*GridsPerPixelY/2,i],dat.rho[nxcenter, nycenter+ScaleFactor*GridsPerPixelY/2,i],dat.Ez[nxcenter, nycenter+ScaleFactor*GridsPerPixelY/2,i])
"""

# Create the output directory if it doesn't exist
if not os.path.isdir(outputfiledir+"/plots"):
    os.mkdir(outputfiledir+"/plots")

rcParams['contour.negative_linestyle'] = 'solid'
rcParams.update({'font.size': 6})

print "Making 1D Potential and Charge Density plots\n"
figure()

suptitle("1D Potential and Charge Density Slices. Grid = %d*%d*%d."%(nxx,nyy,nzz),fontsize = 18)
plotcounter = 1
subplots_adjust(hspace=0.3, wspace=0.3)
phinumzs = 160
numzs = 160
elecnumzs = ConfigData["Nzelec"] * ConfigData["ScaleFactor"]

centers = [(0,0), (-3,3), (0,3), (3,3), (-3,-3), (0,-3), (3,-3)]
fluxes = []
fluxes.append(0)
for kk in range(6):
    kkk = 3 * kk
    fluxes.append(ConfigData["CollectedCharge_0_%d"%kkk])

subplot(2,2,1)
title("Phi-Collect Gate", fontsize=12)

for m,(dnx, dny) in enumerate(centers):
    nxcenter2 = nxcenter + dnx * ScaleFactor * GridsPerPixelX
    nycenter2 = nycenter + dny * ScaleFactor * GridsPerPixelY
    plot(dat.z[0:phinumzs],(dat.phi[nxcenter2,nycenter2,0:phinumzs]+dat.phi[nxcenter2-1,nycenter2,0:phinumzs]+dat.phi[nxcenter2,nycenter2-1,0:phinumzs]+dat.phi[nxcenter2-1,nycenter2-1,0:phinumzs])/4.0, label = "%d"%fluxes[m])
legend(loc = "upper right")
#xlabel("Z-Dimension (microns)")
ylabel('$\phi(x,y,z)$ [V]',fontsize=12)
ylim(-10.0, 30.0)
xlim(0.0,4.0)

subplot(2,2,2)
title("Phi-Barrier Gate", fontsize=12)

for m,(dnx, dny) in enumerate(centers):
    nxcenter2 = nxcenter + dnx * ScaleFactor * GridsPerPixelX
    nycenter2 = nycenter + dny * ScaleFactor * GridsPerPixelY + ScaleFactor * GridsPerPixelY / 2
    plot(dat.z[0:phinumzs],(dat.phi[nxcenter2,nycenter2,0:phinumzs]+dat.phi[nxcenter2-1,nycenter2,0:phinumzs]+dat.phi[nxcenter2,nycenter2-1,0:phinumzs]+dat.phi[nxcenter2-1,nycenter2-1,0:phinumzs])/4.0, label = "%d"%fluxes[m])
legend(loc = "upper right")
#xlabel("Z-Dimension (microns)")
ylabel('$\phi(x,y,z)$ [V]',fontsize=12)
ylim(-10.0, 30.0)
xlim(0.0,4.0)

subplot(2,2,3)
title("Rho-Collect Gate", fontsize=12)
nxcenter2 = nxcenter
nycenter2 = nycenter
plot(dat.z[0:numzs], (dat.rho[nxcenter2,nycenter2,0:numzs]+dat.rho[nxcenter2-1,nycenter2,0:numzs]+dat.rho[nxcenter2,nycenter2-1,0:numzs]+dat.rho[nxcenter2-1,nycenter2-1,0:numzs])/4.0, label = "Fixed charge", color='green')
for m,(dnx, dny) in enumerate(centers):
    nxcenter2 = nxcenter + dnx * ScaleFactor * GridsPerPixelX
    nycenter2 = nycenter + dny * ScaleFactor * GridsPerPixelY
    plot(dat.z[0:elecnumzs], -ChargeFactor * (dat.Elec[nxcenter2,nycenter2,0:elecnumzs]+dat.Elec[nxcenter2-1,nycenter2,0:elecnumzs]+dat.Elec[nxcenter2,nycenter2-1,0:elecnumzs]+dat.Elec[nxcenter2-1,nycenter2-1,0:elecnumzs])/4.0 / dat.dz[0:elecnumzs], label = "%d"%fluxes[m])

legend(loc = "upper right")
xlabel("Z-Dimension (microns)", fontsize=12)
ylabel('$\\rho(x,y,z)/\epsilon_{Si}$ [V/um$^2$]',fontsize=12)
ylim(-80.0, 80.0)
xlim(0.0,4.0)

slicez = 8 * ScaleFactor
ddny = 3 * ScaleFactor * GridsPerPixelY / 2
subplot(2,2,4)
title("Potential Barrier, z = %.2f"%dat.z[slicez])
for m,(dnx, dny) in enumerate(centers):
    nxcenter2 = nxcenter + dnx * ScaleFactor * GridsPerPixelX
    nycenter2 = nycenter + dny * ScaleFactor * GridsPerPixelY
    nymin = nycenter2 - ddny
    nymax = nycenter2 + ddny
    deltay = dat.y[nycenter2]
    plot(dat.y[nymin:nymax]-deltay,dat.phi[nxcenter2,nymin:nymax, slicez], label = "%d"%fluxes[m])
ylim(-5.0, 15.0)
xlim(dat.y[nycenter-ddny]-dat.y[nycenter], dat.y[nycenter+ddny]-dat.y[nycenter])
xlabel("Y-Dimension (microns)")
ylabel("Potential(Volts)")
legend(loc = "lower right")


savefig(outputfiledir+"/plots/"+outputfilebase+"_1D_Multi_%d.pdf"%run)

