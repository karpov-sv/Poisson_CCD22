#!/usr/bin/env python

#Author: Craig Lage, NYU;
#Date: 10-Nov-16

#This program plots the Poisson equation solutions from the C++ Poisson solver
import matplotlib
matplotlib.use("PDF")
from pylab import *
import os, sys, time, h5py

#****************SUBROUTINES*****************
class Array3dHDF5(object):
    def __init__(self, dir, filebase, LogEField, run):
        phifile = dir+'/'+filebase+'_'+str(run)+'_phi' + '.hdf5'
        rhofile = dir+'/'+filebase+'_'+str(run)+'_rho' + '.hdf5'
        xfile = dir+'/'+'grid_x.dat'
        yfile = dir+'/'+'grid_y.dat'
        zfile = dir+'/'+'grid_z.dat'        

        xgrid = loadtxt(xfile, skiprows=1)
        ygrid = loadtxt(yfile, skiprows=1)
        zgrid = loadtxt(zfile, skiprows=1)        

        self.nx=xgrid.shape[0]
        self.ny=ygrid.shape[0]
        self.nz=zgrid.shape[0]

        self.xmin=xgrid[0,1]
        self.ymin=ygrid[0,1]
        self.zmin=zgrid[0,1]

        self.xmax=xgrid[self.nx-1,3]
        self.ymax=ygrid[self.ny-1,3]
        self.zmax=zgrid[self.nz-1,3]

        self.x=xgrid[:,2]
        self.y=ygrid[:,2]
        self.z=zgrid[:,2]
        self.dz=(zgrid[:,3] - zgrid[:,1])        

        hdfphi = h5py.File(phifile,'r')
        self.phi=array(hdfphi[hdfphi.items()[0][0]])
        hdfrho = h5py.File(rhofile,'r')
        self.rho=array(hdfrho[hdfrho.items()[0][0]])
        if LogEField == 1:
            Exfile = dir+'/'+filebase+'_'+str(run)+'_Ex' + '.hdf5'
            Eyfile = dir+'/'+filebase+'_'+str(run)+'_Ey' + '.hdf5'
            Ezfile = dir+'/'+filebase+'_'+str(run)+'_Ez' + '.hdf5'
            hdfEx = h5py.File(Exfile,'r')
            self.Ex=array(hdfEx[hdfEx.items()[0][0]])
            hdfEy = h5py.File(Eyfile,'r')
            self.Ey=array(hdfEy[hdfEy.items()[0][0]])
            hdfEz = h5py.File(Ezfile,'r')
            self.Ez=array(hdfEz[hdfEz.items()[0][0]])

        elecfile = dir+'/'+filebase+'_'+str(run)+'_Elec' + '.hdf5'
        hdfelec = h5py.File(elecfile,'r')
        holefile = dir+'/'+filebase+'_'+str(run)+'_Hole' + '.hdf5'
        hdfhole = h5py.File(holefile,'r')
        self.elec=array(hdfelec[hdfelec.items()[0][0]])
        self.hole=array(hdfhole[hdfhole.items()[0][0]])

def ReadConfigFile(filename):
    # This reads the Poisson simulator config file for
    # the settings that were run
    # and returns a dictionary with the values
    ConfigData = {}
    try:
        file = open(filename,'r')
        lines=file.readlines()
        file.close()
    except IOError:
        print "Configuration file %s not found"%filename
        return False, ConfigData 

    try:
        for line in lines:
            ThisLine=line.strip().split()
            ThisLineLength=len(ThisLine)
            if ThisLineLength < 3:
                continue
            if list(ThisLine[0])[0]=='#' or ThisLine[0]=='\n':
                continue
            try:
                ParamName = ThisLine[0]
                ThisLine.remove(ThisLine[0])
                for counter,item in enumerate(ThisLine):
                    if list(item)[0] == '#':
                        del ThisLine[counter:] # Strip the rest of the line as a comment
                        continue
                    if item == '=':
                        ThisLine.remove(item)
                        continue
                if len(ThisLine) == 0:
                    continue
                elif len(ThisLine) == 1:
                    ThisParam = ThisLine[0]
                    try: ConfigData[ParamName] = int(ThisParam)
                    except ValueError:
                        try:
                            ConfigData[ParamName] = float(ThisParam)
                        except ValueError:
                            try:
                                ConfigData[ParamName] = ThisParam
                            except ValueError:
                                return False, ConfigData 
                else:
                    ThisParam = []
                    for item in ThisLine:
                        try: ThisParam.append(int(item))
                        except ValueError:
                            try: ThisParam.append(float(item))
                            except ValueError:
                                ThisParam.append(item)
                    ConfigData[ParamName] = ThisParam
            except (IOError, ValueError):
                continue
    except Exception as e:
        print "Error reading configuration file %s. Exception of type %s and args = \n"%(filename,type(e).__name__), e.args 
        return False, ConfigData 

    return True, ConfigData

#****************MAIN PROGRAM*****************

# First, read the .cfg file

configfile = sys.argv[1]
run = int(sys.argv[2])
cfg_success, ConfigData = ReadConfigFile(configfile)
if not cfg_success:
    print "Configuration file issue. Quitting"
    sys.exit()
outputfilebase = ConfigData["outputfilebase"]
outputfiledir = ConfigData["outputfiledir"]

# This holds all of the data
dat = Array3dHDF5(outputfiledir, outputfilebase, ConfigData["LogEField"], run)

ScaleFactor = ConfigData["ScaleFactor"]
GridsPerPixel = ConfigData["GridsPerPixelX"]
ChargeFactor = 1.6E-19 * 1.0E6 / (11.7 * 8.85E-12)/((dat.x[1]-dat.x[0])*(dat.y[1]-dat.y[0])) #(QE*MICRON_PER_M/(EPSILON_0*EPSILON_SI))/(dx*dy)


nxx = dat.nx - 1
nyy = dat.ny - 1
nzz = dat.nz - 1

# A couple of things to customize the plots
PlotEField = bool(ConfigData["PlotEField"])
EdgePlot = bool(ConfigData["EdgePlot"])

nxcenter = nxx/2
nycenter = nyy/2
nxcenter2 = nxcenter
nycenter2 = nycenter

nxmin = 0
nxmax = nxx
nymin = 0
nymax = nyy

nzmin = 0
nzmax = nzz

"""
print "Vdd contact"
xplot = 41 * ScaleFactor
yplot = 90 * ScaleFactor
for i in range(dat.elec.shape[2]):
    print "k = %d, z=%.3f, phi = %.3f, holes = %.3f, elec = %.3f, rho = %.3f, Ez = %.3f"%(i,dat.z[i],dat.phi[xplot,yplot,i],dat.hole[xplot,yplot,i]*ChargeFactor/dat.dz[i],dat.elec[xplot,yplot,i]*ChargeFactor/dat.dz[i],dat.rho[xplot,yplot,i],dat.Ez[xplot,yplot,i])

print "Vdd Region"
xplot = 41 * ScaleFactor
yplot = 72 * ScaleFactor
for i in range(dat.elec.shape[2]):
    print "k = %d, z=%.3f, phi = %.3f, holes = %.3f, elec = %.3f, rho = %.3f, Ez = %.3f"%(i,dat.z[i],dat.phi[xplot,yplot,i],dat.hole[xplot,yplot,i]*ChargeFactor/dat.dz[i],dat.elec[xplot,yplot,i]*ChargeFactor/dat.dz[i],dat.rho[xplot,yplot,i],dat.Ez[xplot,yplot,i])
"""
print "OG"
xplot = 41 * ScaleFactor
yplot = 48 * ScaleFactor
for i in range(dat.elec.shape[2]):
    print "k = %d, z=%.3f, phi = %.3f, holes = %.3f, elec = %.3f, rho = %.3f, Ez = %.3f"%(i,dat.z[i],dat.phi[xplot,yplot,i],dat.hole[xplot,yplot,i]*ChargeFactor/dat.dz[i],dat.elec[xplot,yplot,i]*ChargeFactor/dat.dz[i],dat.rho[xplot,yplot,i],dat.Ez[xplot,yplot,i])

print "OD"
xplot = 30 * ScaleFactor
yplot = 48 * ScaleFactor
for i in range(dat.elec.shape[2]):
    print "k = %d, z=%.3f, phi = %.3f, holes = %.3f, elec = %.3f, rho = %.3f, Ez = %.3f"%(i,dat.z[i],dat.phi[xplot,yplot,i],dat.hole[xplot,yplot,i]*ChargeFactor/dat.dz[i],dat.elec[xplot,yplot,i]*ChargeFactor/dat.dz[i],dat.rho[xplot,yplot,i],dat.Ez[xplot,yplot,i])
"""
print "Field Region"
xplot = 47 * ScaleFactor
yplot = 17 * ScaleFactor
for i in range(dat.elec.shape[2]):
    print "k = %d, z=%.3f, phi = %.3f, holes = %.3f, elec = %.3f, rho = %.3f, Ez = %.3f"%(i,dat.z[i],dat.phi[xplot,yplot,i],dat.hole[xplot,yplot,i]*ChargeFactor/dat.dz[i],dat.elec[xplot,yplot,i]*ChargeFactor/dat.dz[i],dat.rho[xplot,yplot,i],dat.Ez[xplot,yplot,i])

print "Hole Region"
xplot = 20
for i in range(dat.elec.shape[2]):
    print "k = %d, z=%.3f, phi = %.3f, holes = %.3f, rho = %.3f, Ez = %.3f"%(i,dat.z[i],dat.phi[xplot,nycenter,i],dat.hole[xplot,nycenter,i],dat.rho[xplot,nycenter,i],dat.Ez[xplot,nycenter,i])

print "Center"
xplot = nxcenter
for i in range(dat.elec.shape[2]):
    print "k = %d, z=%.3f, phi = %.3f, holes = %.3f, rho = %.3f, Ez = %.3f"%(i,dat.z[i],dat.phi[xplot,nycenter,i],dat.hole[xplot,nycenter,i],dat.rho[xplot,nycenter,i],dat.Ez[xplot,nycenter,i])

print "Right Edge"
#xplot = 560 # Top and Bottom
xplot = 150 # Left and Right
for i in range(dat.elec.shape[2]):
    print "k = %d, z=%.3f, phi = %.3f, elec = %.3f, rho = %.3f, Ez = %.3f"%(i,dat.z[i],dat.phi[xplot,nycenter,i],dat.elec[xplot,nycenter,i],dat.rho[xplot,nycenter,i],dat.Ez[xplot,nycenter,i])
"""
# Create the output directory if it doesn't exist
if not os.path.isdir(outputfiledir+"/plots"):
    os.mkdir(outputfiledir+"/plots")

rcParams['contour.negative_linestyle'] = 'solid'
rcParams.update({'font.size': 6})

"""
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
if EdgePlot:
    ylim(-75.0, 25.0)
else:
    ylim(-20.0, 20.0)
for slicez in [0,1,2,3,10]:
    plot(dat.y[:],dat.phi[0,:,slicez], label = '$z_0=%.1f$'%dat.z[slicez])
plt.xlabel('$y$ [um]')
plt.ylabel('$\phi(x_L,y,z_0)$ [V]')
legend()

subplot(2,2,4)
title("Right Edge")
if EdgePlot:
    ylim(-75.0, 25.0)
else:
    ylim(-20.0, 20.0)
for slicez in [0,1,2,3,10]:
    plot(dat.y[:],dat.phi[dat.nx-1,:,slicez], label = '$z_0=%.1f$'%dat.z[slicez])
plt.xlabel('$y$ [um]')
plt.ylabel('$\phi(x_R,y,z_0)$ [V]')
legend()

savefig(outputfiledir+"/plots/"+outputfilebase+"_Edge_Potentials_%d.pdf"%run)

"""

print "Making 1D Potential and Charge Density plots\n"
figure()

suptitle("1D Potential and Charge Density Slices. Grid = %d*%d*%d."%(nxx,nyy,nzz),fontsize = 18)
plotcounter = 1
subplots_adjust(hspace=0.3, wspace=0.3)
phinumzs = 160 * ScaleFactor
numzs = 160 * ScaleFactor
elecnumzs = ConfigData["Nzelec"] * ConfigData["ScaleFactor"]

subplot(2,2,1)
title("Phi - S/D Regions", fontsize=12)
nxcenter2 = 30 * ScaleFactor
nycenter2 = 48 * ScaleFactor
plot(dat.z[0:phinumzs],dat.phi[nxcenter2,nycenter2,0:phinumzs], label = "OD")
nxcenter2 = 54 * ScaleFactor
nycenter2 = 48 * ScaleFactor
plot(dat.z[0:phinumzs],dat.phi[nxcenter2,nycenter2,0:phinumzs], label = "OS")
xlabel("Z-Dimension (microns)")
ylabel('$\phi(z)$ [V]',fontsize=12)
ylim(-1.0, 1.0)
xlim(0.0,10.0)
legend(loc = 'upper right')

subplot(2,2,2)
title("Phi - Gate Region", fontsize=12)
nxcenter2 = 41 * ScaleFactor
nycenter2 = 48 * ScaleFactor
plot(dat.z[0:phinumzs],dat.phi[nxcenter2,nycenter2,0:phinumzs], label = "Gate")
xlabel("Z-Dimension (microns)")
ylabel('$\phi(z)$ [V]',fontsize=12)
ylim(-2.0, 2.0)
xlim(0.0,10.0)
legend(loc = 'upper right')

subplot(2,2,3)
title("Rho-S/D Regions", fontsize=12)
nxcenter2 = 30 * ScaleFactor
nycenter2 = 48 * ScaleFactor
plot(dat.z[0:numzs], dat.rho[nxcenter2,nycenter2,0:numzs], label = "Fixed charge", color='green')
plot(dat.z[0:elecnumzs], -ChargeFactor * dat.elec[nxcenter2,nycenter2,0:elecnumzs] / dat.dz[0:elecnumzs], label = "Electrons", color='blue')
plot(dat.z[0:elecnumzs], ChargeFactor * dat.hole[nxcenter2,nycenter2,0:elecnumzs] / dat.dz[0:elecnumzs], label = "Holes", color='red')
legend(loc = "lower right")
xlabel("Z-Dimension (microns)", fontsize=12)
ylabel('$\\rho(x,y,z)/\epsilon_{Si}$ [V/um$^2$]',fontsize=12)
ylim(-80.0, 80.0)
xlim(0.0,4.0)

subplot(2,2,4)
title("Rho-Gate Region", fontsize=12)
nxcenter2 = 41 * ScaleFactor
nycenter2 = 48 * ScaleFactor
plot(dat.z[0:numzs], dat.rho[nxcenter2,nycenter2,0:numzs], label = "Fixed charge", color='green')
plot(dat.z[0:elecnumzs], -ChargeFactor * dat.elec[nxcenter2,nycenter2,0:elecnumzs] / dat.dz[0:elecnumzs], label = "Electrons", color='blue')
plot(dat.z[0:elecnumzs], ChargeFactor * dat.hole[nxcenter2,nycenter2,0:elecnumzs] / dat.dz[0:elecnumzs], label = "Holes", color='red')
legend(loc = "lower right")
xlabel("Z-Dimension (microns)", fontsize=12)
ylabel('$\\rho(x,y,z)/\epsilon_{Si}$ [V/um$^2$]',fontsize=12)
ylim(-80.0, 80.0)
xlim(0.0,4.0)
savefig(outputfiledir+"/plots/"+outputfilebase+"_1D_Potentials_%d.pdf"%run)


print "Making summary plots\n"
figure()
suptitle("Potentials on Bottom. Grid = %d*%d*%d."%(nxx,nyy,nzz),fontsize = 18)
[yy,xx] = meshgrid(dat.y[nymin:nymax],dat.x[nxmin:nxmax])

slicez = 0
subplot(1,2,1, aspect = 1)
levels = linspace(-10.0, 26.0, 37)
title ("Z = %.2f"%dat.z[slicez])
contour(xx,yy,dat.phi[nxmin:nxmax,nymin:nymax,slicez],levels,lw=0.1)
contourf(xx,yy,dat.phi[nxmin:nxmax,nymin:nymax,slicez],levels)
xlabel("X-Dimension (microns)")
ylabel("Y-Dimension (microns)")
cb=colorbar(orientation='horizontal')
cb.set_ticks([-10.0,0.0,10.0,20.0,30.0])
cb.set_label('$\phi(x,y)$ [V]')

slicez = 32
subplot(1,2,2, aspect = 1)
levels = linspace(-5.0, 15.0, 41)
title ("Z = %.2f"%dat.z[slicez])
contour(xx,yy,dat.phi[nxmin:nxmax,nymin:nymax,slicez],levels,lw=0.1)
contourf(xx,yy,dat.phi[nxmin:nxmax,nymin:nymax,slicez],levels)
xlabel("X-Dimension (microns)")
ylabel("Y-Dimension (microns)")
cb=colorbar(orientation='horizontal')
cb.set_ticks([-5.0,0.0,5.0,10.0,15.0])
cb.set_label('$\phi(x,y)$ [V]')

savefig(outputfiledir+"/plots/"+outputfilebase+"_Summary_1_%d.pdf"%run)

print "Making summary plots\n"
figure()
suptitle("Potentials on Bottom. Grid = %d*%d*%d."%(nxx,nyy,nzz),fontsize = 18)
[yy,xx] = meshgrid(dat.y[nymin:nymax],dat.x[nxmin:nxmax])

slicez = 0
subplot(1,2,1, aspect = 1)
levels = linspace(-10.0, 26.0, 37)
title ("Z = %.2f"%dat.z[slicez])
contour(xx,yy,dat.phi[nxmin:nxmax,nymin:nymax,slicez],levels,lw=0.1)
contourf(xx,yy,dat.phi[nxmin:nxmax,nymin:nymax,slicez],levels)
xlabel("X-Dimension (microns)")
ylabel("Y-Dimension (microns)")
cb=colorbar(orientation='horizontal')
cb.set_ticks([-10.0,0.0,10.0,20.0,30.0])
cb.set_label('$\phi(x,y)$ [V]')

slicez = 9
subplot(1,2,2, aspect = 1)
levels = linspace(-10.0, 20.0, 61)
title ("Z = %.2f"%dat.z[slicez])
contour(xx,yy,dat.phi[nxmin:nxmax,nymin:nymax,slicez],levels,lw=0.1)
contourf(xx,yy,dat.phi[nxmin:nxmax,nymin:nymax,slicez],levels)
xlabel("X-Dimension (microns)")
ylabel("Y-Dimension (microns)")
cb=colorbar(orientation='horizontal')
cb.set_ticks([-10.0,0.0,10.0,20.0])
cb.set_label('$\phi(x,y)$ [V]')

savefig(outputfiledir+"/plots/"+outputfilebase+"_Summary_2_%d.pdf"%run)

figure()
suptitle("Vertical Cross-Section. Grid = %d*%d*%d."%(nxx,nyy,nzz),fontsize = 18)
subplots_adjust(hspace=0.3, wspace=0.3)
nxcenter = 38 * ScaleFactor
nymin = 5 * ScaleFactor
nymax = 120 * ScaleFactor
nycenter = 42 * ScaleFactor
nxmin = 16 * ScaleFactor
nxmax = 64 * ScaleFactor
subplot(1,2,1)
title("Phi(x,z) y = %.2f"%dat.y[nycenter])
xlabel("X-Dimension (microns)")
ylabel("Z-Dimension (microns)")
nzmin = 0
nzmax = 160 * ScaleFactor
[zz,xx] = meshgrid(dat.z[nzmin:nzmax],dat.x[nxmin:nxmax])
contour(xx,zz,dat.phi[nxmin:nxmax,nycenter,nzmin:nzmax],levels,lw=0.1)
contourf(xx,zz,dat.phi[nxmin:nxmax,nycenter,nzmin:nzmax],levels)
cb=colorbar(orientation='horizontal')
cb.set_ticks([-60.0,-45.0,-30.0,-15.0,0.0,15.0])
cb.set_label('$\phi(x,z)$ [V]')
[zz,xx] = meshgrid(dat.z[nzmin:nzmax],dat.x[nxmin:nxmax])
ylim(zz[0,0], zz[-1,-1])
xlim(xx[0,0], xx[-1,-1])

subplot(1,2,2)
title("Phi(y,z) x = %.2f"%dat.x[nxcenter])
xlabel("Y-Dimension (microns)")
ylabel("Z-Dimension (microns)")
nzmin = 0
nzmax = 160 * ScaleFactor

[zz,yy] = meshgrid(dat.z[nzmin:nzmax],dat.y[nymin:nymax])
contour(yy,zz,dat.phi[nxcenter,nymin:nymax,nzmin:nzmax],levels,lw=0.1)
contourf(yy,zz,dat.phi[nxcenter,nymin:nymax,nzmin:nzmax],levels)
cb=colorbar(orientation='horizontal')
cb.set_ticks([-60.0,-45.0,-30.0,-15.0,0.0,15.0])
cb.set_label('$\phi(y,z)$ [V]')
[zz,yy] = meshgrid(dat.z[nzmin:nzmax],dat.y[nymin:nymax])
ylim(zz[0,0], zz[-1,-1])
xlim(yy[0,0], yy[-1,-1])

savefig(outputfiledir+"/plots/"+outputfilebase+"_Summary_3_%d.pdf"%run)


nzmax = 40 * ScaleFactor

figure()
suptitle("Vertical Cross-Section. Grid = %d*%d*%d."%(nxx,nyy,nzz),fontsize = 18)
subplots_adjust(hspace=0.3, wspace=0.3)
nxcenter = 38 * ScaleFactor
nymin = 5 * ScaleFactor
nymax = 120 * ScaleFactor
nycenter = 42 * ScaleFactor
nxmin = 16 * ScaleFactor
nxmax = 64 * ScaleFactor
subplot(1,2,1)
title("Elec(x,z) y = %.2f"%dat.y[nycenter])
xlabel("X-Dimension (microns)")
ylabel("Z-Dimension (microns)")
nzmin = 0
nzmax = 160 * ScaleFactor
[zz,xx] = meshgrid(dat.z[nzmin:nzmax],dat.x[nxmin:nxmax])
contour(xx,zz,dat.elec[nxmin:nxmax,nycenter,nzmin:nzmax],lw=0.1)
contourf(xx,zz,dat.elec[nxmin:nxmax,nycenter,nzmin:nzmax])
cb=colorbar(orientation='horizontal')
cb.set_ticks([-60.0,-45.0,-30.0,-15.0,0.0,15.0])
cb.set_label('elec(x,z) [V]')
[zz,xx] = meshgrid(dat.z[nzmin:nzmax],dat.x[nxmin:nxmax])
ylim(zz[0,0], zz[-1,-1])
xlim(xx[0,0], xx[-1,-1])

subplot(1,2,2)
title("Elec(y,z) x = %.2f"%dat.x[nxcenter])
xlabel("Y-Dimension (microns)")
ylabel("Z-Dimension (microns)")
nzmin = 0
nzmax = 160 * ScaleFactor

[zz,yy] = meshgrid(dat.z[nzmin:nzmax],dat.y[nymin:nymax])
contour(yy,zz,dat.elec[nxcenter,nymin:nymax,nzmin:nzmax],lw=0.1)
contourf(yy,zz,dat.elec[nxcenter,nymin:nymax,nzmin:nzmax])
cb=colorbar(orientation='horizontal')
cb.set_ticks([-60.0,-45.0,-30.0,-15.0,0.0,15.0])
cb.set_label('elec(y,z) [V]')
[zz,yy] = meshgrid(dat.z[nzmin:nzmax],dat.y[nymin:nymax])
ylim(zz[0,0], zz[-1,-1])
xlim(yy[0,0], yy[-1,-1])

savefig(outputfiledir+"/plots/"+outputfilebase+"_Summary_4_%d.pdf"%run)

figure()
suptitle("Vertical Cross-Section. Grid = %d*%d*%d."%(nxx,nyy,nzz),fontsize = 18)
subplots_adjust(hspace=0.3, wspace=0.3)
nxcenter = 38 * ScaleFactor
nymin = 5 * ScaleFactor
nymax = 120 * ScaleFactor
nycenter = 42 * ScaleFactor
nxmin = 16 * ScaleFactor
nxmax = 64 * ScaleFactor
subplot(1,2,1)
title("Hole(x,z) y = %.2f"%dat.y[nycenter])
xlabel("X-Dimension (microns)")
ylabel("Z-Dimension (microns)")
nzmin = 0
nzmax = 160 * ScaleFactor
[zz,xx] = meshgrid(dat.z[nzmin:nzmax],dat.x[nxmin:nxmax])
contour(xx,zz,dat.hole[nxmin:nxmax,nycenter,nzmin:nzmax],lw=0.1)
contourf(xx,zz,dat.hole[nxmin:nxmax,nycenter,nzmin:nzmax])
cb=colorbar(orientation='horizontal')
cb.set_ticks([-60.0,-45.0,-30.0,-15.0,0.0,15.0])
cb.set_label('hole(x,z) [V]')
[zz,xx] = meshgrid(dat.z[nzmin:nzmax],dat.x[nxmin:nxmax])
ylim(zz[0,0], zz[-1,-1])
xlim(xx[0,0], xx[-1,-1])

subplot(1,2,2)
title("Hole(y,z) x = %.2f"%dat.x[nxcenter])
xlabel("Y-Dimension (microns)")
ylabel("Z-Dimension (microns)")
nzmin = 0
nzmax = 160 * ScaleFactor

[zz,yy] = meshgrid(dat.z[nzmin:nzmax],dat.y[nymin:nymax])
contour(yy,zz,dat.hole[nxcenter,nymin:nymax,nzmin:nzmax],lw=0.1)
contourf(yy,zz,dat.hole[nxcenter,nymin:nymax,nzmin:nzmax])
cb=colorbar(orientation='horizontal')
cb.set_ticks([-60.0,-45.0,-30.0,-15.0,0.0,15.0])
cb.set_label('hole(y,z) [V]')
[zz,yy] = meshgrid(dat.z[nzmin:nzmax],dat.y[nymin:nymax])
ylim(zz[0,0], zz[-1,-1])
xlim(yy[0,0], yy[-1,-1])

savefig(outputfiledir+"/plots/"+outputfilebase+"_Summary_5_%d.pdf"%run)
"""
figure()
suptitle("Hole Concentration (log scale). Grid = %d*%d*%d."%(nxx,nyy,nzz),fontsize = 18)
subplots_adjust(hspace=0.3, wspace=0.3)
levels = linspace(-2.0, 4.0, 161)
subplot(1,2,1)
title("Holes (Vertical dimension magnified). y = %.2f"%dat.y[nycenter2])
xlabel("X-Dimension (microns)")
ylabel("Z-Dimension (microns)")
nzmin = 0
nzmax = 160 * ScaleFactor
[zz,xx] = meshgrid(dat.z[nzmin:nzmax],dat.x[nxmin:nxmax])
contour(xx,zz,clip(log10(dat.hole[nxmin:nxmax,nycenter2,nzmin:nzmax]+1.0E-9),-2.0,10.0),levels,lw=0.1)
contourf(xx,zz,clip(log10(dat.hole[nxmin:nxmax,nycenter2,nzmin:nzmax]+1.0E-9),-2.0,10.0),levels)
cb=colorbar(orientation='horizontal')
cb.set_ticks([-1.0,0.0,1.0,2.0,3.0])
cb.set_label('log10(holes(x,z))')
ylim(zz[0,0], zz[-1,-1])
xlim(xx[0,0], xx[-1,-1])

subplot(1,2,2,aspect=1)
title("Holes (No vertical magnification). y = %.2f"%dat.y[nycenter2])
xlabel("X-Dimension (microns)")
ylabel("Z-Dimension (microns)")
nzmin = 0
nzmax = 160 * ScaleFactor
[zz,xx] = meshgrid(dat.z[nzmin:nzmax],dat.x[nxmin:nxmax])
contour(xx,zz,clip(log10(dat.hole[nxmin:nxmax,nycenter2,nzmin:nzmax]+1.0E-9),-2.0,10.0),levels,lw=0.1)
contourf(xx,zz,clip(log10(dat.hole[nxmin:nxmax,nycenter2,nzmin:nzmax]+1.0E-9),-2.0,10.0),levels)
cb=colorbar(orientation='horizontal')
cb.set_ticks([-1.0,0.0,1.0,2.0,3.0])
cb.set_label('log10(holes(x,z))')
ylim(zz[0,0], zz[-1,-1])
xlim(xx[0,0], xx[-1,-1])
savefig(outputfiledir+"/plots/"+outputfilebase+"_Summary_3_%d.pdf"%run)

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

if EdgePlot:
    for linex in linspace(120.0,200.0,9):
        plot((linex,linex),(20.0,70.0),linewidth=1.0, color='blue')

xlabel("X(microns)",fontsize = 18)
ylabel("Y(microns)",fontsize = 18)
xlim(ConfigData["PixelBoundaryLowerLeft"][0], ConfigData["PixelBoundaryUpperRight"][0])
ylim(ConfigData["PixelBoundaryLowerLeft"][1], ConfigData["PixelBoundaryUpperRight"][1])


savefig(outputfiledir+"/plots/"+outputfilebase+"_Pixels_%d.pdf"%run)


if ConfigData["LogPixelPaths"] != 0 and run % ConfigData["LogPixelPaths"] == 0:
    # Last, plots of the electron paths
    print "Making array electron path plots\n"
    # Plotting the paths along a line through the center
    yline = (ConfigData["PixelBoundaryLowerLeft"][1] + ConfigData["PixelBoundaryUpperRight"][1] + ConfigData["PixelBoundaryStepSize"][1]) / 2.0
    xline = (ConfigData["PixelBoundaryLowerLeft"][0] + ConfigData["PixelBoundaryUpperRight"][0] + ConfigData["PixelBoundaryStepSize"][0]) / 2.0
    #print xline, yline
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
                if phase == 2:
                #if zout < 1.20:
                    #print "Finished this path", xin, yin 
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
                if phase == 2:
                #if zout < 1.20:                    
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

"""
