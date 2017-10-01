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
nxcenter2 = int(((ConfigData["PixelRegionLowerLeft_0"][0] + ConfigData["PixelRegionUpperRight_0"][0]) / 2.0 - ConfigData["SimulationRegionLowerLeft"][0]) / (GridsPerPixelX * ScaleFactor))
nycenter2 = int(((ConfigData["PixelRegionLowerLeft_0"][1] + ConfigData["PixelRegionUpperRight_0"][1]) / 2.0 - ConfigData["SimulationRegionLowerLeft"][1]) / (GridsPerPixelY * ScaleFactor))

nxmin = 0
nxmax = nxx
nymin = 0
nymax = nyy

nzmin = 0
nzmax = 16 * ScaleFactor

# Create the output directory if it doesn't exist
if not os.path.isdir(outputfiledir+"/plots"):
    os.mkdir(outputfiledir+"/plots")

rcParams['contour.negative_linestyle'] = 'solid'
rcParams.update({'font.size': 6})

print "Making array edge potential plots\n"
figure()
suptitle("Array Edge Potentials. Grid = %d*%d*%d."%(nxx,nyy,nzz),fontsize = 18)
subplot(2,2,1)
title("Front Edge")
ylim(-20.0, 20.0)
for slicez in [0,1,2,3,10]:
    plot(dat.x[:],dat.phi[:,0,slicez], label = '$z_0=%.1f$'%dat.z[slicez])
plt.xlabel('$x$ [um]')
plt.ylabel('$\phi(x,y_F,z_0)$ [V]')
legend()

subplot(2,2,2)
title("Back Edge")
ylim(-80.0, 20.0)
for slicez in [0,1,2,3,10]:
    plot(dat.x[:],dat.phi[:,dat.ny-1,slicez], label = '$z_0=%.1f$'%dat.z[slicez])
plt.xlabel('$x$ [um]')
plt.ylabel('$\phi(x,y_B,z_0)$ [V]')
legend()

subplot(2,2,3)
title("Left Edge")
ylim(-75.0, 25.0)
for slicez in [0,1,2,3,10]:
    plot(dat.y[:],dat.phi[0,:,slicez], label = '$z_0=%.1f$'%dat.z[slicez])
plt.xlabel('$y$ [um]')
plt.ylabel('$\phi(x_L,y,z_0)$ [V]')
legend()

subplot(2,2,4)
title("Right Edge")
ylim(-75.0, 25.0)
for slicez in [0,1,2,3,10]:
    plot(dat.y[:],dat.phi[dat.nx-1,:,slicez], label = '$z_0=%.1f$'%dat.z[slicez])
plt.xlabel('$y$ [um]')
plt.ylabel('$\phi(x_R,y,z_0)$ [V]')
legend()

savefig(outputfiledir+"/plots/"+outputfilebase+"_Edge_Potentials_%d.pdf"%run)

print "Making 1D potential and Charge Density plots\n"
figure()

suptitle("1D Potential and Charge Density Slices. Grid = %d*%d*%d."%(nxx,nyy,nzz),fontsize = 18)
plotcounter = 1
subplots_adjust(hspace=0.3, wspace=0.3)
phinumzs = 160
numzs = 160
elecnumzs = ConfigData["Nzelec"] * ConfigData["ScaleFactor"]

subplot(2,3,1)
title("Phi-Collect Gate", fontsize=12)

nxcenter3 = nxcenter2 + 2 * GridsPerPixelX*ScaleFactor

plot(dat.z[0:phinumzs],(dat.phi[nxcenter3,nycenter2,0:phinumzs]+dat.phi[nxcenter3-1,nycenter2,0:phinumzs]+dat.phi[nxcenter3,nycenter2-1,0:phinumzs]+dat.phi[nxcenter3-1,nycenter2-1,0:phinumzs])/4.0, label = "Empty Well")
legend(loc = "lower left")
#xlabel("Z-Dimension (microns)")
ylabel('$\phi(x,y,z)$ [V]',fontsize=12)
ylim(-5.0, 25.0)
xlim(0.0,4.0)

subplot(2,3,4)
title("Rho-Collect Gate", fontsize=12)
plot(dat.z[0:numzs], dat.rho[nxcenter2,nycenter2,0:numzs], label = "Fixed charge", color='green')
plot(dat.z[0:elecnumzs], -ChargeFactor * dat.Elec[nxcenter3,nycenter2,0:elecnumzs] / dat.dz[0:elecnumzs], label = "Empty well", color='magenta')

legend(loc = "upper right")
xlabel("Z-Dimension (microns)", fontsize=12)
ylabel('$\\rho(x,y,z)/\epsilon_{Si}$ [V/um$^2$]',fontsize=12)
ylim(-80.0, 250.0)
xlim(0.0,4.0)
nxcenter3 = nxcenter2 + GridsPerPixelX * ScaleFactor / 2
nycenter3 = nycenter2
nycenter4 = nycenter2 + GridsPerPixelY * ScaleFactor / 2
subplot(2,3,2)
title("Phi-ChanStop", fontsize=12)
plot(dat.z[0:phinumzs],dat.phi[nxcenter3,nycenter3,0:phinumzs], label = "CollectGate")
plot(dat.z[0:phinumzs],dat.phi[nxcenter3,nycenter4,0:phinumzs], label = "Barrier Gate")
legend(loc = "lower left")
#xlabel("Z-Dimension (microns)")
#ylabel('$\phi(x,y,z)$ [V]',fontsize=9)
ylim(-10.0, 12.0)
xlim(0.0,10.0)
subplot(2,3,5)
title("Rho-ChanStop", fontsize=12)

plot(dat.z[0:numzs], dat.rho[nxcenter3,nycenter3,0:numzs], color = 'green', label = "Fixed charge")

plot(dat.z[0:elecnumzs], ChargeFactor / dat.dz[0:elecnumzs] * dat.Hole[nxcenter3,nycenter3,0:elecnumzs], label = "Holes Collect Gate", color = 'red')

plot(dat.z[0:elecnumzs], ChargeFactor / dat.dz[0:elecnumzs] * dat.Hole[nxcenter3,nycenter4,0:elecnumzs], label = "Holes Barrier Gate", color = 'orange', )

legend(loc = "lower right")
xlabel("Z-Dimension (microns)", fontsize=12)
#ylabel('$\\rho(x,y,z)/\epsilon_{Si}$ [V/um$^2$]',fontsize=9)
ylim(-100.0, 100.0)
xlim(0.0,10.0)

nxcenter3 = nxcenter2
nycenter3 = nycenter2 + GridsPerPixelY * ScaleFactor / 2
subplot(2,3,3)
title("Phi-Barrier Gate", fontsize=12)
plot(dat.z[0:phinumzs],dat.phi[nxcenter3,nycenter3,0:phinumzs])
#xlabel("Z-Dimension (microns)")
#ylabel('$\phi(x,y,z)$ [V]',fontsize=9)
ylim(-10.0, 12.0)
xlim(0.0,4.0)
subplot(2,3,6)
title("Rho-Barrier Gate", fontsize=12)

plot(dat.z[0:numzs], dat.rho[nxcenter3,nycenter3,0:numzs], color = 'green', label = "Fixed charge")

plot(dat.z[0:elecnumzs], ChargeFactor / dat.dz[0:elecnumzs] * dat.Hole[nxcenter3,nycenter3,0:elecnumzs], color = 'red', label = "Holes")

legend(loc = "lower left")
xlabel("Z-Dimension (microns)", fontsize=12)
#ylabel('$\\rho(x,y,z)/\epsilon_{Si}$ [V/um$^2$]',fontsize=9)
ylim(-80.0, 250.0)
xlim(0.0,4.0)
savefig(outputfiledir+"/plots/"+outputfilebase+"_1D_Potentials_%d.pdf"%run)

print "Making summary plots\n"
figure()
suptitle("Potentials. Grid = %d*%d*%d."%(nxx,nyy,nzz),fontsize = 18)
[yy,xx] = meshgrid(dat.y[nymin:nymax],dat.x[nxmin:nxmax])
subplots_adjust(hspace=0)

subplot(2,1,1)
plot(dat.y[:],dat.phi[nxcenter,:,0], label="z=%.2f"%dat.z[0])
plot(dat.y[:],dat.phi[nxcenter,:,10], label="z=%.2f"%dat.z[10])
plot(dat.y[:],dat.phi[nxcenter,:,int(dat.nz/8)], label="z=%.2f"%dat.z[int(dat.nz/8)])
plot(dat.y[:],dat.phi[nxcenter,:,int(dat.nz/4)], label="z=%.2f"%dat.z[int(dat.nz/4)])
plot(dat.y[:],dat.phi[nxcenter,:,int(dat.nz/2)], label="z=%.2f"%dat.z[int(dat.nz/2)])
plot(dat.y[:],dat.phi[nxcenter,:,int(3*dat.nz/4)], label="z=%.2f"%dat.z[int(3*dat.nz/4)])
zmaxx = min(100.0, dat.z[dat.nz-1])
plot(dat.y[:],dat.phi[nxcenter,:,dat.nz-1], label="z=%.2f"%zmaxx)
ylabel('Potential(V)')
xlabel('Y (micron)')
legend()

slicez = 0
subplot(2,1,2, aspect=1)
title("Potentials at z=0", fontsize = 12)
levels = linspace(-60.0, 20.0, 161)
contour(yy,xx,dat.phi[nxmin:nxmax,nymin:nymax,slicez],levels,lw=0.1)
contourf(yy,xx,dat.phi[nxmin:nxmax,nymin:nymax,slicez],levels)
ylabel("X-Dimension (microns)")
yticks([40.0, 70.0])
xlabel("Y-Dimension (microns)")
cb=colorbar(orientation='horizontal')
cb.set_ticks([-60.0,-45.0,-30.0,-15.0,0.0,15.0])
cb.set_label('$\phi(x,y)$ [V]')
savefig(outputfiledir+"/plots/"+outputfilebase+"_Bottom_Potential_Vbb_%.0f.pdf"%ConfigData["Vbb"])


# Next, plots of the pixel boundaries
print "Making pixel plots\n"
figure()
rcParams['contour.negative_linestyle'] = 'solid'
#rcParams.update({'font.size': 18})

suptitle("CCD Pixel Plots. Grid = %d*%d*%d."%(nxx,nyy,nzz),fontsize = 24)
plotcounter = 1
subplots_adjust(hspace=0.3, wspace=0.1)

filename = outputfiledir+"/"+outputfilebase+'_'+str(run)+"_Pts.dat"
file = open(filename,"r")
lines = file.readlines()
file.close()
if len(lines) < 2:
    print "No data in Pts file.  Quitting"
    sys.exit()
redsx=[]
redsy=[]
blacksx=[]
blacksy=[]
plottedxin = -1000.0
plottedyin = -1000.0
lines.remove(lines[0])
for line in lines:
    values = line.split()
    phase = int(values[2])
    #zout = float(values[5])        
    if phase == 0:
        xin = float(values[3])
        yin = float(values[4])
    elif phase == 2:
    #elif zout < 1.20:
        xout = float(values[3])
        yout = float(values[4])
        if isnan(xout) or isnan(yout):
            print "xin = %.3f, yin = %.3f is a nan"
            continue
        pixxout = int(xout/10.0)
        pixyout = int(yout/10.0)
        if (pixxout + pixyout) % 2 == 0:
            redsx.append(xin)
            redsy.append(yin)
        else:
            blacksx.append(xin)
            blacksy.append(yin)
        continue
    else:
        continue

subplot(1,1,1,aspect=1)
title("Pixel Boundaries",fontsize = 12)
if ConfigData["PixelBoundaryTestType"] == 0:
    spotsize = 10.0 * ConfigData["PixelBoundaryStepSize"][0] * ConfigData["PixelBoundaryStepSize"][1]
else:
    spotsize = 0.1
scatter(redsx,redsy,s=spotsize,color="red")
scatter(blacksx,blacksy,s=spotsize,color="black")

for linex in linspace(120.0,200.0,9):
    plot((linex,linex),(20.0,70.0),linewidth=1.0, color='blue')

xlabel("X(microns)",fontsize = 18)
xticks([40.0, 70.0])
ylabel("Y(microns)",fontsize = 18)
xlim(ConfigData["PixelBoundaryLowerLeft"][0], ConfigData["PixelBoundaryUpperRight"][0])
ylim(ConfigData["PixelBoundaryLowerLeft"][1], ConfigData["PixelBoundaryUpperRight"][1])


savefig(outputfiledir+"/plots/"+outputfilebase+"_Pixels_%d.pdf"%run)
