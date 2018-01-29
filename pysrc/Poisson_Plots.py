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

NumPixelsPlotted = 4
nycenter2 = nycenter
nxmin = nxcenter - (NumPixelsPlotted * ScaleFactor * GridsPerPixelX)/2
nxmax = nxcenter + (NumPixelsPlotted * ScaleFactor * GridsPerPixelX)/2
nymin = nycenter - (NumPixelsPlotted * ScaleFactor * GridsPerPixelY)/2
nymax = nycenter + (NumPixelsPlotted * ScaleFactor * GridsPerPixelY)/2

nzmin = 0
nzmax = 16 * ScaleFactor

print "Channel"
for i in range(dat.Elec.shape[2]):
    print "k = %d, z=%.3f, phi = %.3f, holes = %.3f, elec = %.3f, rho = %.3f, Ez = %.3f"%(i,dat.z[i],dat.phi[nxcenter, nycenter,i],dat.Hole[nxcenter, nycenter,i],dat.Elec[nxcenter, nycenter,i],dat.rho[nxcenter, nycenter,i],dat.Ez[nxcenter, nycenter,i])

nxcenter3 = nxcenter2 + GridsPerPixelX * ScaleFactor / 2
print "Channel Stop"
for i in range(dat.Elec.shape[2]):
    print "k = %d, z=%.3f, phi = %.3f, holes = %.3f, elec = %.3f, rho = %.3f, Ez = %.3f"%(i,dat.z[i],dat.phi[nxcenter3, nycenter,i],dat.Hole[nxcenter3, nycenter,i],dat.Elec[nxcenter3, nycenter,i],dat.rho[nxcenter3, nycenter,i],dat.Ez[nxcenter3, nycenter,i])
    
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
plot(dat.z[0:phinumzs],(dat.phi[nxcenter2,nycenter2,0:phinumzs]+dat.phi[nxcenter2-1,nycenter2,0:phinumzs]+dat.phi[nxcenter2,nycenter2-1,0:phinumzs]+dat.phi[nxcenter2-1,nycenter2-1,0:phinumzs])/4.0, label = "80K electrons")

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
plot(dat.z[0:elecnumzs], -ChargeFactor * dat.Elec[nxcenter2,nycenter2,0:elecnumzs] / dat.dz[0:elecnumzs], label = "80K electrons", color='blue')
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

print "Making 1D potential Plots #2 \n"
figure()
suptitle("1D Potentials in Storage Region. Grid = %d*%d*%d."%(nxx,nyy,nzz),fontsize = 18)
subplots_adjust(hspace=0.3, wspace=0.3)
slicez = 16 * ScaleFactor
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
slicezs = [4,8,12,14,16,18,20,24]
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
slicezs = array([8,9,10,12]) * ScaleFactor
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
nxmin2 = nxmin - 9 * GridsPerPixelX * ScaleFactor / 2
nymin2 = nymin - 9 * GridsPerPixelY * ScaleFactor / 2
rho0 = dat.rho[nxmin2,nymin2,slicez+1]

title("Rho, z = %.2f - %.2f"%(dat.z[nzmin],dat.z[nzmax]))
levels = linspace(-40.0,20.0,61)
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
    if phase == 0:
        xin = float(values[3])
        yin = float(values[4])
    elif phase == 3 or phase == 4:
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

xlabel("X(microns)",fontsize = 18)
ylabel("Y(microns)",fontsize = 18)
xlim(ConfigData["PixelBoundaryLowerLeft"][0], ConfigData["PixelBoundaryUpperRight"][0])
ylim(ConfigData["PixelBoundaryLowerLeft"][1], ConfigData["PixelBoundaryUpperRight"][1])


savefig(outputfiledir+"/plots/"+outputfilebase+"_Pixels_%d.pdf"%run)


if ConfigData["LogPixelPaths"] != 0 and run % ConfigData["LogPixelPaths"] == 0:
    # Last, plots of the electron paths
    print "Making array electron path plots\n"
    # Plotting the paths along a line through the center
    yline = (ConfigData["PixelBoundaryLowerLeft"][1] + ConfigData["PixelBoundaryUpperRight"][1]) / 2.0
    xline = (ConfigData["PixelBoundaryLowerLeft"][0] + ConfigData["PixelBoundaryUpperRight"][0]) / 2.0

    vertical_zoom = 1
    figure()
    suptitle("Electron Path Plot - Vertical Zoom = %d"%vertical_zoom, fontsize = 24)
    subplots_adjust(wspace=0.2)

    for line in lines:
        values = line.split()
        phase = int(values[2])
        if phase == 0:
            xin = float(values[3])
            yin = float(values[4])
            if (yin > yline - .10) and (yin < yline + .10):
                YPlotThisID = True
                xpaths=[]
                zxpaths=[]
            else:
                YPlotThisID = False
            if (xin > xline - .10) and (xin < xline + .10):
                XPlotThisID = True
                ypaths=[]
                zypaths=[]
            else:
                XPlotThisID = False
            continue
        if XPlotThisID or YPlotThisID:
            xout = float(values[3])
            yout = float(values[4])
            zout = float(values[5])
            if isnan(xout) or isnan(yout) or isnan(zout):
                continue
            if YPlotThisID:
                xpaths.append(xout)
                zxpaths.append(zout)
                if phase == 3 or phase == 4:
                    pixxin = int(xin/10.0)
                    if pixxin % 2 == 0:
                        color = "red"
                    else:
                        color = "black"
                    subplot(1,2,1,aspect=vertical_zoom)
                    plot(xpaths, zxpaths, color = color, linewidth = 0.1)

            if XPlotThisID:
                ypaths.append(yout)
                zypaths.append(zout)
                if phase == 3 or phase == 4:
                    pixyin = int(yin/10.0)
                    if pixyin % 2 == 0:
                        color = "red"
                    else:
                        color = "black"
                    subplot(1,2,2,aspect=vertical_zoom)                        
                    plot(ypaths, zypaths, color = color, linewidth = 0.1)

    subplot(1,2,1,aspect=vertical_zoom)
    ylabel("Z(microns)")
    xlabel("X (microns)")
    ylim(0.0,110.0)
    xlim(ConfigData["PixelBoundaryLowerLeft"][0], ConfigData["PixelBoundaryUpperRight"][0])
    subplot(1,2,2,aspect=vertical_zoom)
    ylabel("Z(microns)")
    xlabel("Y (microns)")
    ylim(0.0,110.0)
    xlim(ConfigData["PixelBoundaryLowerLeft"][1], ConfigData["PixelBoundaryUpperRight"][1])
    savefig(outputfiledir+"/plots/"+outputfilebase+"_Paths_%d.pdf"%run)

