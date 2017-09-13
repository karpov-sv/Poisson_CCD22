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

NumPixelsPlotted = 1
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
ylim(-20.0, 20.0)
for slicez in [0,1,2,3,10]:
    plot(dat.x[:],dat.phi[:,dat.ny-1,slicez], label = '$z_0=%.1f$'%dat.z[slicez])
plt.xlabel('$x$ [um]')
plt.ylabel('$\phi(x,y_B,z_0)$ [V]')
legend()

subplot(2,2,3)
title("Left Edge")
ylim(-20.0, 20.0)
for slicez in [0,1,2,3,10]:
    plot(dat.y[:],dat.phi[0,:,slicez], label = '$z_0=%.1f$'%dat.z[slicez])
plt.xlabel('$y$ [um]')
plt.ylabel('$\phi(x_L,y,z_0)$ [V]')
legend()

subplot(2,2,4)
title("Right Edge")
ylim(-20.0, 20.0)
for slicez in [0,1,2,3,10]:
    plot(dat.y[:],dat.phi[dat.nx-1,:,slicez], label = '$z_0=%.1f$'%dat.z[slicez])
plt.xlabel('$y$ [um]')
plt.ylabel('$\phi(x_R,y,z_0)$ [V]')
legend()

savefig(outputfiledir+"/plots/"+outputfilebase+"_Edge_Potentials_%d.pdf"%run)

print "Making 1D Potential and Charge Density plots\n"
figure()

suptitle("1D Potential and Charge Density Slices. Grid = %d*%d*%d."%(nxx,nyy,nzz),fontsize = 18)
plotcounter = 1
subplots_adjust(hspace=0.3, wspace=0.3)
phinumzs = 160 * ScaleFactor
numzs = 160 * ScaleFactor
elecnumzs = ConfigData["Nzelec"] * ConfigData["ScaleFactor"]

subplot(2,3,1)
title("Phi-Collect Gate", fontsize=12)
nxcenter3 = nxcenter2

plot(dat.z[0:phinumzs],(dat.phi[nxcenter3,nycenter2,0:phinumzs]+dat.phi[nxcenter3-1,nycenter2,0:phinumzs]+dat.phi[nxcenter3,nycenter2-1,0:phinumzs]+dat.phi[nxcenter3-1,nycenter2-1,0:phinumzs])/4.0, label = "Filled Well")
legend(loc = "lower left")
#xlabel("Z-Dimension (microns)")
ylabel('$\phi(x,y,z)$ [V]',fontsize=12)
ylim(-5.0, 25.0)
xlim(0.0,4.0)

subplot(2,3,4)
title("Rho-Collect Gate", fontsize=12)
plot(dat.z[0:numzs], (dat.rho[nxcenter2,nycenter2,0:numzs]+dat.rho[nxcenter2-1,nycenter2,0:numzs]+dat.rho[nxcenter2,nycenter2-1,0:numzs]+dat.rho[nxcenter2-1,nycenter2-1,0:numzs])/4.0, label = "Fixed charge", color='green')
plot(dat.z[0:elecnumzs], -ChargeFactor * (dat.Elec[nxcenter3,nycenter2,0:elecnumzs]+dat.Elec[nxcenter3-1,nycenter2,0:elecnumzs]+dat.Elec[nxcenter3,nycenter2-1,0:elecnumzs]+dat.Elec[nxcenter3-1,nycenter2-1,0:elecnumzs])/4.0 / dat.dz[0:elecnumzs], label = "Empty well", color='magenta')

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
plot(dat.z[0:phinumzs],(dat.phi[nxcenter3,nycenter3,0:phinumzs]+dat.phi[nxcenter3-1,nycenter3,0:phinumzs]+dat.phi[nxcenter3,nycenter3-1,0:phinumzs]+dat.phi[nxcenter3-1,nycenter3-1,0:phinumzs])/4.0, label = "CollectGate")
plot(dat.z[0:phinumzs],(dat.phi[nxcenter3,nycenter4,0:phinumzs]+dat.phi[nxcenter3-1,nycenter4,0:phinumzs]+dat.phi[nxcenter3,nycenter4-1,0:phinumzs]+dat.phi[nxcenter3-1,nycenter4-1,0:phinumzs])/4.0, label = "Barrier Gate")
legend(loc = "lower left")
#xlabel("Z-Dimension (microns)")
#ylabel('$\phi(x,y,z)$ [V]',fontsize=9)
ylim(-10.0, 12.0)
xlim(0.0,10.0)
subplot(2,3,5)
title("Rho-ChanStop", fontsize=12)

plot(dat.z[0:numzs], ((dat.rho[nxcenter3,nycenter3,0:numzs]+dat.rho[nxcenter3-1,nycenter3,0:numzs]+dat.rho[nxcenter3,nycenter3-1,0:numzs]+dat.rho[nxcenter3-1,nycenter3-1,0:numzs])/4.0), color = 'green', label = "Fixed charge")

plot(dat.z[0:elecnumzs], ChargeFactor / dat.dz[0:elecnumzs] * ((dat.Hole[nxcenter3,nycenter3,0:elecnumzs]+dat.Hole[nxcenter3-1,nycenter3,0:elecnumzs]+dat.Hole[nxcenter3,nycenter3-1,0:elecnumzs]+dat.Hole[nxcenter3-1,nycenter3-1,0:elecnumzs])/4.0), label = "Holes Collect Gate", color = 'red')

plot(dat.z[0:elecnumzs], ChargeFactor / dat.dz[0:elecnumzs] * ((dat.Hole[nxcenter3,nycenter4,0:elecnumzs]+dat.Hole[nxcenter3-1,nycenter4,0:elecnumzs]+dat.Hole[nxcenter3,nycenter4-1,0:elecnumzs]+dat.Hole[nxcenter3-1,nycenter4-1,0:elecnumzs])/4.0), label = "Holes Barrier Gate", color = 'orange', )

legend(loc = "lower right")
xlabel("Z-Dimension (microns)", fontsize=12)
#ylabel('$\\rho(x,y,z)/\epsilon_{Si}$ [V/um$^2$]',fontsize=9)
ylim(-100.0, 100.0)
xlim(0.0,10.0)

nxcenter3 = nxcenter2
nycenter3 = nycenter2 + GridsPerPixelY * ScaleFactor / 2
subplot(2,3,3)
title("Phi-Barrier Gate", fontsize=12)
plot(dat.z[0:phinumzs],(dat.phi[nxcenter3,nycenter3,0:phinumzs]+dat.phi[nxcenter3-1,nycenter3,0:phinumzs]+dat.phi[nxcenter3,nycenter3-1,0:phinumzs]+dat.phi[nxcenter3-1,nycenter3-1,0:phinumzs])/4.0)
#xlabel("Z-Dimension (microns)")
#ylabel('$\phi(x,y,z)$ [V]',fontsize=9)
ylim(-10.0, 12.0)
xlim(0.0,4.0)
subplot(2,3,6)
title("Rho-Barrier Gate", fontsize=12)

plot(dat.z[0:numzs], ((dat.rho[nxcenter3,nycenter3,0:numzs]+dat.rho[nxcenter3-1,nycenter3,0:numzs]+dat.rho[nxcenter3,nycenter3-1,0:numzs]+dat.rho[nxcenter3-1,nycenter3-1,0:numzs])/4.0), color = 'green', label = "Fixed charge")

plot(dat.z[0:elecnumzs], ChargeFactor / dat.dz[0:elecnumzs] * ((dat.Hole[nxcenter3,nycenter3,0:elecnumzs]+dat.Hole[nxcenter3-1,nycenter3,0:elecnumzs]+dat.Hole[nxcenter3,nycenter3-1,0:elecnumzs]+dat.Hole[nxcenter3-1,nycenter3-1,0:elecnumzs])/4.0), color = 'red', label = "Holes")

legend(loc = "lower left")
xlabel("Z-Dimension (microns)", fontsize=12)
#ylabel('$\\rho(x,y,z)/\epsilon_{Si}$ [V/um$^2$]',fontsize=9)
ylim(-80.0, 250.0)
xlim(0.0,4.0)
savefig(outputfiledir+"/plots/"+outputfilebase+"_1D_Potentials_%d.pdf"%run)

print "Making 1D potential Plots #2 \n"
figure()
suptitle("1D Potentials in Storage Region. Grid = %d*%d*%d."%(nxx,nyy,nzz),fontsize = 18)
subplots_adjust(hspace=0.3, wspace=0.3)
slicez = 8 * ScaleFactor
subplot(1,2,1, aspect = 1)
title("Phi, z = %.2f"%dat.z[slicez])
levels = linspace(-20.0, 20.0, 21)
[yy,xx] = meshgrid(dat.y[nymin:nymax],dat.x[nxmin:nxmax])
contour(xx,yy,dat.phi[nxmin:nxmax,nymin:nymax,slicez],levels,lw=0.1)
contourf(xx,yy,dat.phi[nxmin:nxmax,nymin:nymax,slicez],levels)
xlabel("X-Dimension (microns)")
ylabel("Y-Dimension (microns)")
plot([dat.x[nxmin+1],dat.x[nxmax-1]],[dat.y[nycenter2],dat.y[nycenter2]],ls = "-", color="k")
plot([dat.x[nxcenter2],dat.x[nxcenter2]],[dat.y[nymin+1],dat.y[nymax-1]],ls = "-", color="k")
colorbar(orientation='horizontal').set_label('$\phi(x,y,z)$ [V]')

subplot(1,2,2)
title("Phi-Collect Gate, z = %.2f"%dat.z[slicez])
plot(dat.x[nxmin:nxmax],dat.phi[nxmin:nxmax, nycenter2, slicez], label = "XSlice, y = %.2f"%dat.y[nycenter2])
plot(dat.y[nymin:nymax],dat.phi[nxcenter2,nymin:nymax, slicez], label = "YSlice, x = %.2f"%dat.x[nxcenter2])
ylim(-10.0, 20.0)
xlim(dat.x[nxmin],dat.x[nxmax])
xlabel("X,Y-Dimension (microns)")
ylabel("Potential(Volts)")
legend()
savefig(outputfiledir+"/plots/"+outputfilebase+"_1D_Potentials_2_%d.pdf"%run)

print "Making 1D potential Plots #3 \n"
slicezs = [4,6,8,10,12,14,16,18]
figure()
suptitle("1D Potentials in Storage Region. Grid = %d*%d*%d."%(nxx,nyy,nzz),fontsize = 18)
subplots_adjust(hspace=0.3, wspace=0.3)

for i,slicez in enumerate(slicezs):
    subplot(2,4,i+1)
    title("Phi, z = %.2f"%dat.z[slicez])
    plot(dat.x[nxmin:nxmax],dat.phi[nxmin:nxmax, nycenter2, slicez], label = "XSlice, y = %.2f"%dat.y[nycenter2])
    plot(dat.y[nymin:nymax],dat.phi[nxcenter2,nymin:nymax, slicez], label = "YSlice, x = %.2f"%dat.x[nxcenter2])
    ylim(-10.0, 20.0)
    xlim(dat.x[nxmin],dat.x[nxmax])
    xlabel("X,Y-Dimension (microns)")
    ylabel("Potential(Volts)")
    legend()
savefig(outputfiledir+"/plots/"+outputfilebase+"_1D_Potentials_3_%d.pdf"%run)

print "Making 1D potential Plots #5 \n"
slicezs = [4,8,12,16]
figure()
suptitle("1D Potentials in Storage Region. Grid = %d*%d*%d."%(nxx,nyy,nzz),fontsize = 18)
subplots_adjust(hspace=0.3, wspace=0.3)

for i,slicez in enumerate(slicezs):
    subplot(2,4,i+1)
    title("Phi, z = %.2f"%dat.z[slicez])
    plot(dat.x[nxmin:nxmax],dat.phi[nxmin:nxmax, nycenter2, slicez], label = "XSlice, y = %.2f"%dat.y[nycenter2])
    plot(dat.y[nymin:nymax],dat.phi[nxcenter2,nymin:nymax, slicez], label = "YSlice, x = %.2f"%dat.x[nxcenter2])
    ylim(-10.0, 20.0)
    xlim(dat.x[nxmin],dat.x[nxmax])
    xlabel("X,Y-Dimension (microns)")
    ylabel("Potential(Volts)")
    legend()
for i,slicez in enumerate(slicezs):
    subplot(2,4,i+5)
    title("Holes, z = %.2f"%dat.z[slicez])
    plot(dat.x[nxmin:nxmax],dat.Hole[nxmin:nxmax, nycenter2, slicez], label = "XSlice, y = %.2f"%dat.y[nycenter2])
    plot(dat.y[nymin:nymax],dat.Hole[nxcenter2,nymin:nymax, slicez], label = "YSlice, x = %.2f"%dat.x[nxcenter2])
    #ylim(-5.0, 0.0)
    xlim(dat.x[nxmin],dat.x[nxmax])
    xlabel("X,Y-Dimension (microns)")
    ylabel("Holes")
    legend()
savefig(outputfiledir+"/plots/"+outputfilebase+"_1D_Potentials_5_%d.pdf"%run)

print "Making 1D potential Plots #4 \n"
figure()

suptitle("1D Potentials in Isolation Regions. Grid = %d*%d*%d."%(nxx,nyy,nzz),fontsize = 18)
subplots_adjust(hspace=0.3, wspace=0.3)
slicez = 16 * ScaleFactor
nxcenter3 = nxcenter2 + GridsPerPixelX * ScaleFactor / 2
nycenter3 = nycenter2 + GridsPerPixelY * ScaleFactor / 2
subplot(1,2,1, aspect = 1)
title("Phi, z = %.2f"%dat.z[slicez])
levels = linspace(-20.0, 20.0, 21)
[yy,xx] = meshgrid(dat.y[nymin:nymax],dat.x[nxmin:nxmax])
contour(xx,yy,dat.phi[nxmin:nxmax,nymin:nymax,slicez],levels,lw=0.1)
contourf(xx,yy,dat.phi[nxmin:nxmax,nymin:nymax,slicez],levels)
xlabel("X-Dimension (microns)")
ylabel("Y-Dimension (microns)")
plot([dat.x[nxmin+1],dat.x[nxmax-1]],[dat.y[nycenter3],dat.y[nycenter3]],ls = "-", color="k")
plot([dat.x[nxcenter3],dat.x[nxcenter3]],[dat.y[nymin+1],dat.y[nymax-1]],ls = "-", color="k")
#colorbar()

subplot(1,2,2)
title("Phi-Collect Gate, z = %.2f"%dat.z[slicez])
plot(dat.x[nxmin:nxmax],dat.phi[nxmin:nxmax, nycenter3, slicez], label = "XSlice, y = %.2f"%dat.y[nycenter3])
plot(dat.y[nymin:nymax],dat.phi[nxcenter3,nymin:nymax, slicez], label = "YSlice, x = %.2f"%dat.x[nxcenter3])
ylim(-10.0, 10.0)
xlim(dat.x[nxmin],dat.x[nxmax])
xlabel("X,Y-Dimension (microns)")
ylabel("Potential(Volts)")
legend()
savefig(outputfiledir+"/plots/"+outputfilebase+"_1D_Potentials_4_%d.pdf"%run)


print "Making summary plots\n"
figure()
suptitle("CCD Charge Collection. Grid = %d*%d*%d."%(nxx,nyy,nzz),fontsize = 18)
plotcounter = 1
subplots_adjust(hspace=0.3, wspace=0.3)
[yy,xx] = meshgrid(dat.y[nymin:nymax],dat.x[nxmin:nxmax])

slicez = 0
subplot(2,2,1, aspect = 1)
title("Phi, z = 0.0")
levels = linspace(-10.0, 10.0, 31)
contour(xx,yy,dat.phi[nxmin:nxmax,nymin:nymax,slicez],levels,lw=0.1)
contourf(xx,yy,dat.phi[nxmin:nxmax,nymin:nymax,slicez],levels)
xlabel("X-Dimension (microns)")
ylabel("Y-Dimension (microns)")
colorbar().set_label('$\phi(x,y,z)$ [V]')

subplot(2,2,2, aspect = 1)
nxmin2 = nxmin -  GridsPerPixelX * ScaleFactor / 2
nymin2 = nymin -  GridsPerPixelY * ScaleFactor / 2
rho0 = dat.rho[nxmin2,nymin2,slicez+1]

title("Rho, z = %.2f - %.2f"%(dat.z[nzmin],dat.z[nzmax]))
levels = linspace(-2.0,15.0,34)
plotarray = array(dat.rho[nxmin:nxmax,nymin:nymax,nzmin:nzmax].sum(axis=2)/(nzmax-nzmin))
contour(xx,yy,plotarray, levels, lw=0.1)
contourf(xx,yy,plotarray, levels)
xlabel("X-Dimension (microns)")
ylabel("Y-Dimension (microns)")
colorbar().set_label('$\\rho(x,y,z) / \epsilon_{Si}$ [V/um$^2$]')

slicez = 8 * ScaleFactor
subplot(2,2,3, aspect = 1)
title("Phi, z = %.2f"%dat.z[slicez])
levels = linspace(-20.0, 20.0, 21)
contour(xx,yy,dat.phi[nxmin:nxmax,nymin:nymax,slicez],levels,lw=0.1)
contourf(xx,yy,dat.phi[nxmin:nxmax,nymin:nymax,slicez],levels)
xlabel("X-Dimension (microns)")
ylabel("Y-Dimension (microns)")
plot([dat.x[nxmin+1],dat.x[nxmax-1]],[dat.y[nycenter2],dat.y[nycenter2]],ls = "-", color="k")
plot([dat.x[nxcenter2],dat.x[nxcenter2]],[dat.y[nymin+1],dat.y[nymax-1]],ls = "-", color="k")
colorbar().set_label('$\phi(x,y,z)$ [V]')

slicez = 16 * ScaleFactor
subplot(2,2,4, aspect = 1)
title("Phi, z = %.2f"%dat.z[slicez])
contour(xx,yy,dat.phi[nxmin:nxmax,nymin:nymax,slicez],levels,lw=0.1)
contourf(xx,yy,dat.phi[nxmin:nxmax,nymin:nymax,slicez],levels)
xlabel("X-Dimension (microns)")
ylabel("Y-Dimension (microns)")
colorbar().set_label('$\phi(x,y,z)$ [V]')

savefig(outputfiledir+"/plots/"+outputfilebase+"_Summary_1_%d.pdf"%run)

figure()
suptitle("CCD Charge Collection. Grid = %d*%d*%d."%(nxx,nyy,nzz),fontsize = 18)
subplots_adjust(hspace=0.3, wspace=0.3)

levels = linspace(-20.0, 20.0, 41)

subplot(1,2,1)
title("Phi and (-)E in Gate Region. y = %.2f"%dat.y[nycenter2])
xlabel("X-Dimension (microns)")
ylabel("Z-Dimension (microns)")
nzmin = 0
nzmax = 16 * ScaleFactor
[zz,xx] = meshgrid(dat.z[nzmin:nzmax],dat.x[nxmin:nxmax])
contour(xx,zz,dat.phi[nxmin:nxmax,nycenter2,nzmin:nzmax],levels,lw=0.1)
contourf(xx,zz,dat.phi[nxmin:nxmax,nycenter2,nzmin:nzmax],levels)
colorbar().set_label('$\phi(x,y,z)$ [V]')
if ConfigData["LogEField"] == 1:
    nzmin = 1
    [zz,xx] = meshgrid(dat.z[nzmin:nzmax],dat.x[nxmin:nxmax])
    quiver(xx, zz, dat.Ex[nxmin:nxmax,nycenter2,nzmin:nzmax], dat.Ez[nxmin:nxmax,nycenter2,nzmin:nzmax], color='b', scale = 150.0)#, scale_units="width")
[zz,xx] = meshgrid(dat.z[nzmin:nzmax],dat.x[nxmin:nxmax])
ylim(zz[0,0], zz[-1,-1])
xlim(xx[0,0], xx[-1,-1])


subplot(1,2,2)
title("Phi and (-)E in Gate Region. x = %.2f"%dat.x[nxcenter2])
xlabel("Y-Dimension (microns)")
ylabel("Z-Dimension (microns)")
nzmin = 0
[zz,yy] = meshgrid(dat.z[nzmin:nzmax],dat.y[nymin:nymax])
contour(yy,zz,dat.phi[nxcenter2,nymin:nymax,nzmin:nzmax],levels,lw=0.1)
contourf(yy,zz,dat.phi[nxcenter2,nymin:nymax,nzmin:nzmax],levels)
colorbar().set_label('$\phi(x,y,z)$ [V]')
if ConfigData["LogEField"] == 1:
    nzmin = 1
    [zz,yy] = meshgrid(dat.z[nzmin:nzmax],dat.y[nymin:nymax])
    quiver(yy, zz, dat.Ey[nxcenter2,nymin:nymax,nzmin:nzmax], dat.Ez[nxcenter,nymin:nymax,nzmin:nzmax], color='b', scale = 150.0)#, scale_units="width")
nzmin = 0
[zz,yy] = meshgrid(dat.z[nzmin:nzmax],dat.y[nymin:nymax])
ylim(zz[0,0], zz[-1,-1])
xlim(yy[0,0], yy[-1,-1])

savefig(outputfiledir+"/plots/"+outputfilebase+"_Summary_2_%d.pdf"%run)

