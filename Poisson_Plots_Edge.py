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
GridsPerPixelX = ConfigData["GridsPerPixelX"]
GridsPerPixelY = ConfigData["GridsPerPixelY"]
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

nxmin = 0
nxmax = nxx
nymin = 0
nymax = nyy

nzmin = 0
nzmax = 16 * ScaleFactor

# We need to customize the plot locations for this particular .cfg file

lowerleft = array(ConfigData["SimulationRegionLowerLeft"])

pixelcenter = (array(ConfigData["PixelRegionUpperRight_0"]) + array(ConfigData["PixelRegionLowerLeft_0"])) / 2.0
pixelnxcenter = int((pixelcenter[0] - lowerleft[0]) / ConfigData["PixelSizeX"] * GridsPerPixelX * ScaleFactor)
pixelnycenter = int((pixelcenter[1] - lowerleft[1]) / ConfigData["PixelSizeY"] * GridsPerPixelY * ScaleFactor)
#print pixelcenter, pixelnxcenter, pixelnycenter

serialcenter = (array(ConfigData["FixedRegionUpperRight_9"]) + array(ConfigData["FixedRegionLowerLeft_9"])) / 2.0
serialnxcenter = int((serialcenter[0] - lowerleft[0]) / ConfigData["PixelSizeX"] * GridsPerPixelX * ScaleFactor)
serialnycenter = int((serialcenter[1] - lowerleft[1]) / ConfigData["PixelSizeY"] * GridsPerPixelY * ScaleFactor)
serialcscenter = (array(ConfigData["FixedRegionUpperRight_8"]) + array(ConfigData["FixedRegionLowerLeft_8"])) / 2.0
serialcsnycenter = int((serialcscenter[1] - lowerleft[1]) / ConfigData["PixelSizeY"] * GridsPerPixelY * ScaleFactor)
#print serialcenter, serialnxcenter, serialnycenter
#print serialcscenter, serialcsnycenter

scuppercenter = (array(ConfigData["FixedRegionUpperRight_6"]) + array(ConfigData["FixedRegionLowerLeft_6"])) / 2.0
scuppernxcenter = int((scuppercenter[0] - lowerleft[0]) / ConfigData["PixelSizeX"] * GridsPerPixelX * ScaleFactor)
scuppernycenter = int((scuppercenter[1] - lowerleft[1]) / ConfigData["PixelSizeY"] * GridsPerPixelY * ScaleFactor)
#print scuppercenter, scuppernxcenter, scuppernycenter

vbbringcenter = (array(ConfigData["FixedRegionUpperRight_0"]) + array(ConfigData["FixedRegionLowerLeft_0"])) / 2.0
vbbringnxcenter = int((vbbringcenter[0] - lowerleft[0]) / ConfigData["PixelSizeX"] * GridsPerPixelX * ScaleFactor)
vbbringnycenter = int((vbbringcenter[1] - lowerleft[1]) / ConfigData["PixelSizeY"] * GridsPerPixelY * ScaleFactor)
#print vbbringcenter, vbbringnxcenter, vbbringnycenter

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
ylim(-75.0, 20.0)
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

print "Making 1D Potential and Charge Density plots\n"
figure()

nxcenter2 = pixelnxcenter
nycenter2 = pixelnycenter

suptitle("1D Potential and Charge Density Slices. Pixel Region",fontsize = 18)
plotcounter = 1
subplots_adjust(hspace=0.3, wspace=0.3)
phinumzs = 160
numzs = 160
elecnumzs = ConfigData["Nzelec"] * ConfigData["ScaleFactor"]

subplot(2,3,1)
title("Phi-Collect Gate", fontsize=12)
plot(dat.z[0:phinumzs],(dat.phi[nxcenter2,nycenter2,0:phinumzs]+dat.phi[nxcenter2-1,nycenter2,0:phinumzs]+dat.phi[nxcenter2,nycenter2-1,0:phinumzs]+dat.phi[nxcenter2-1,nycenter2-1,0:phinumzs])/4.0, label = "Empty well")
legend(loc = "lower left")
#xlabel("Z-Dimension (microns)")
ylabel('$\phi(x,y,z)$ [V]',fontsize=12)
ylim(-5.0, 25.0)
xlim(0.0,4.0)

subplot(2,3,4)
title("Rho-Collect Gate", fontsize=12)
plot(dat.z[0:numzs], (dat.rho[nxcenter2,nycenter2,0:numzs]+dat.rho[nxcenter2-1,nycenter2,0:numzs]+dat.rho[nxcenter2,nycenter2-1,0:numzs]+dat.rho[nxcenter2-1,nycenter2-1,0:numzs])/4.0, label = "Fixed charge", color='green')
plot(dat.z[0:elecnumzs], -ChargeFactor * (dat.elec[nxcenter2,nycenter2,0:elecnumzs]+dat.elec[nxcenter2-1,nycenter2,0:elecnumzs]+dat.elec[nxcenter2,nycenter2-1,0:elecnumzs]+dat.elec[nxcenter2-1,nycenter2-1,0:elecnumzs])/4.0 / dat.dz[0:elecnumzs], label = "Empty well", color='magenta')

legend(loc = "lower left")
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
ylim(-10.0, 10.0)
xlim(0.0,10.0)
subplot(2,3,5)
title("Rho-ChanStop", fontsize=12)

plot(dat.z[0:numzs], ((dat.rho[nxcenter3,nycenter3,0:numzs]+dat.rho[nxcenter3-1,nycenter3,0:numzs]+dat.rho[nxcenter3,nycenter3-1,0:numzs]+dat.rho[nxcenter3-1,nycenter3-1,0:numzs])/4.0), color = 'green', label = "Fixed charge")

plot(dat.z[0:elecnumzs], ChargeFactor / dat.dz[0:elecnumzs] * ((dat.hole[nxcenter3,nycenter3,0:elecnumzs]+dat.hole[nxcenter3-1,nycenter3,0:elecnumzs]+dat.hole[nxcenter3,nycenter3-1,0:elecnumzs]+dat.hole[nxcenter3-1,nycenter3-1,0:elecnumzs])/4.0), label = "Holes Collect Gate", color = 'red')

plot(dat.z[0:elecnumzs], ChargeFactor / dat.dz[0:elecnumzs] * ((dat.hole[nxcenter3,nycenter4,0:elecnumzs]+dat.hole[nxcenter3-1,nycenter4,0:elecnumzs]+dat.hole[nxcenter3,nycenter4-1,0:elecnumzs]+dat.hole[nxcenter3-1,nycenter4-1,0:elecnumzs])/4.0), label = "Holes Barrier Gate", color = 'orange', )

legend(loc = "lower left")
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
ylim(-10.0, 10.0)
xlim(0.0,4.0)
subplot(2,3,6)
title("Rho-Barrier Gate", fontsize=12)

plot(dat.z[0:numzs], ((dat.rho[nxcenter3,nycenter3,0:numzs]+dat.rho[nxcenter3-1,nycenter3,0:numzs]+dat.rho[nxcenter3,nycenter3-1,0:numzs]+dat.rho[nxcenter3-1,nycenter3-1,0:numzs])/4.0), color = 'green', label = "Fixed charge")

plot(dat.z[0:elecnumzs], ChargeFactor / dat.dz[0:elecnumzs] * ((dat.hole[nxcenter3,nycenter3,0:elecnumzs]+dat.hole[nxcenter3-1,nycenter3,0:elecnumzs]+dat.hole[nxcenter3,nycenter3-1,0:elecnumzs]+dat.hole[nxcenter3-1,nycenter3-1,0:elecnumzs])/4.0), color = 'red', label = "Holes")

legend(loc = "lower left")
xlabel("Z-Dimension (microns)", fontsize=12)
#ylabel('$\\rho(x,y,z)/\epsilon_{Si}$ [V/um$^2$]',fontsize=9)
ylim(-80.0, 250.0)
xlim(0.0,4.0)
savefig(outputfiledir+"/plots/"+outputfilebase+"_1D_Potentials_%d.pdf"%run)


figure()

nxcenter2 = serialnxcenter
nycenter2 = serialnycenter
nycenter3 = serialcsnycenter

suptitle("1D Potential and Charge Density Slices. Serial Region",fontsize = 18)
plotcounter = 1
subplots_adjust(hspace=0.3, wspace=0.3)
phinumzs = 160
numzs = 160
elecnumzs = ConfigData["Nzelec"] * ConfigData["ScaleFactor"]

subplot(2,2,1)
title("Phi-Collect Gate", fontsize=12)
plot(dat.z[0:phinumzs],(dat.phi[nxcenter2,nycenter2,0:phinumzs]+dat.phi[nxcenter2-1,nycenter2,0:phinumzs]+dat.phi[nxcenter2,nycenter2-1,0:phinumzs]+dat.phi[nxcenter2-1,nycenter2-1,0:phinumzs])/4.0, label = "Empty well")
legend(loc = "lower left")
#xlabel("Z-Dimension (microns)")
ylabel('$\phi(x,y,z)$ [V]',fontsize=12)
ylim(-5.0, 25.0)
xlim(0.0,4.0)

subplot(2,2,3)
title("Rho-Collect Gate", fontsize=12)
plot(dat.z[0:numzs], (dat.rho[nxcenter2,nycenter2,0:numzs]+dat.rho[nxcenter2-1,nycenter2,0:numzs]+dat.rho[nxcenter2,nycenter2-1,0:numzs]+dat.rho[nxcenter2-1,nycenter2-1,0:numzs])/4.0, label = "Fixed charge", color='green')
plot(dat.z[0:elecnumzs], -ChargeFactor * (dat.elec[nxcenter2,nycenter2,0:elecnumzs]+dat.elec[nxcenter2-1,nycenter2,0:elecnumzs]+dat.elec[nxcenter2,nycenter2-1,0:elecnumzs]+dat.elec[nxcenter2-1,nycenter2-1,0:elecnumzs])/4.0 / dat.dz[0:elecnumzs], label = "Empty well", color='magenta')

legend(loc = "lower left")
xlabel("Z-Dimension (microns)", fontsize=12)
ylabel('$\\rho(x,y,z)/\epsilon_{Si}$ [V/um$^2$]',fontsize=12)
ylim(-80.0, 250.0)
xlim(0.0,4.0)
nxcenter3 = nxcenter2 
nycenter4 = nycenter3

subplot(2,2,2)
title("Phi-ChanStop", fontsize=12)
plot(dat.z[0:phinumzs],(dat.phi[nxcenter3,nycenter3,0:phinumzs]+dat.phi[nxcenter3-1,nycenter3,0:phinumzs]+dat.phi[nxcenter3,nycenter3-1,0:phinumzs]+dat.phi[nxcenter3-1,nycenter3-1,0:phinumzs])/4.0, label = "Average Gate")
legend(loc = "lower left")
#xlabel("Z-Dimension (microns)")
#ylabel('$\phi(x,y,z)$ [V]',fontsize=9)
ylim(-10.0, 10.0)
xlim(0.0,10.0)
subplot(2,2,4)
title("Rho-ChanStop", fontsize=12)

plot(dat.z[0:numzs], ((dat.rho[nxcenter3,nycenter3,0:numzs]+dat.rho[nxcenter3-1,nycenter3,0:numzs]+dat.rho[nxcenter3,nycenter3-1,0:numzs]+dat.rho[nxcenter3-1,nycenter3-1,0:numzs])/4.0), color = 'green', label = "Fixed charge")

plot(dat.z[0:elecnumzs], ChargeFactor / dat.dz[0:elecnumzs] * ((dat.hole[nxcenter3,nycenter3,0:elecnumzs]+dat.hole[nxcenter3-1,nycenter3,0:elecnumzs]+dat.hole[nxcenter3,nycenter3-1,0:elecnumzs]+dat.hole[nxcenter3-1,nycenter3-1,0:elecnumzs])/4.0), label = "Holes Average Gate", color = 'red')

legend(loc = "lower left")
xlabel("Z-Dimension (microns)", fontsize=12)
#ylabel('$\\rho(x,y,z)/\epsilon_{Si}$ [V/um$^2$]',fontsize=9)
ylim(-100.0, 100.0)
xlim(0.0,10.0)
savefig(outputfiledir+"/plots/"+outputfilebase+"_1D_Serial_%d.pdf"%run)


figure()

nxcenter2 = scuppernxcenter
nycenter2 = scuppernycenter
nxcenter3 = vbbringnxcenter
nycenter3 = vbbringnycenter

suptitle("1D Potential and Charge Density Slices. Fixed Regions",fontsize = 18)
plotcounter = 1
subplots_adjust(hspace=0.3, wspace=0.3)
phinumzs = 160
numzs = 160
elecnumzs = ConfigData["Nzelec"] * ConfigData["ScaleFactor"]

subplot(2,2,1)
title("Phi-Scupper", fontsize=12)
plot(dat.z[0:phinumzs],(dat.phi[nxcenter2,nycenter2,0:phinumzs]+dat.phi[nxcenter2-1,nycenter2,0:phinumzs]+dat.phi[nxcenter2,nycenter2-1,0:phinumzs]+dat.phi[nxcenter2-1,nycenter2-1,0:phinumzs])/4.0, label = "Scupper")
legend(loc = "lower left")
#xlabel("Z-Dimension (microns)")
ylabel('$\phi(x,y,z)$ [V]',fontsize=12)
ylim(-80.0, 30.0)
xlim(0.0,100.0)

subplot(2,2,3)
title("Rho-Scupper", fontsize=12)
plot(dat.z[0:numzs], (dat.rho[nxcenter2,nycenter2,0:numzs]+dat.rho[nxcenter2-1,nycenter2,0:numzs]+dat.rho[nxcenter2,nycenter2-1,0:numzs]+dat.rho[nxcenter2-1,nycenter2-1,0:numzs])/4.0, label = "Fixed charge", color='green')
plot(dat.z[0:elecnumzs], -ChargeFactor * (dat.elec[nxcenter2,nycenter2,0:elecnumzs]+dat.elec[nxcenter2-1,nycenter2,0:elecnumzs]+dat.elec[nxcenter2,nycenter2-1,0:elecnumzs]+dat.elec[nxcenter2-1,nycenter2-1,0:elecnumzs])/4.0 / dat.dz[0:elecnumzs], label = "Electrons", color='magenta')

legend(loc = "lower left")
xlabel("Z-Dimension (microns)", fontsize=12)
ylabel('$\\rho(x,y,z)/\epsilon_{Si}$ [V/um$^2$]',fontsize=12)
ylim(-80.0, 80.0)
xlim(0.0,10.0)

subplot(2,2,2)
title("Phi-Vbb Ring", fontsize=12)
plot(dat.z[0:phinumzs],(dat.phi[nxcenter3,nycenter3,0:phinumzs]+dat.phi[nxcenter3-1,nycenter3,0:phinumzs]+dat.phi[nxcenter3,nycenter3-1,0:phinumzs]+dat.phi[nxcenter3-1,nycenter3-1,0:phinumzs])/4.0, label = "Vbb Ring")
legend(loc = "lower left")
#xlabel("Z-Dimension (microns)")
#ylabel('$\phi(x,y,z)$ [V]',fontsize=9)
ylim(-80.0, 30.0)
xlim(0.0,100.0)
subplot(2,2,4)
title("Rho-Vbb Ring", fontsize=12)

plot(dat.z[0:numzs], ((dat.rho[nxcenter3,nycenter3,0:numzs]+dat.rho[nxcenter3-1,nycenter3,0:numzs]+dat.rho[nxcenter3,nycenter3-1,0:numzs]+dat.rho[nxcenter3-1,nycenter3-1,0:numzs])/4.0), color = 'green', label = "Fixed charge")

plot(dat.z[0:elecnumzs], ChargeFactor / dat.dz[0:elecnumzs] * ((dat.hole[nxcenter3,nycenter3,0:elecnumzs]+dat.hole[nxcenter3-1,nycenter3,0:elecnumzs]+dat.hole[nxcenter3,nycenter3-1,0:elecnumzs]+dat.hole[nxcenter3-1,nycenter3-1,0:elecnumzs])/4.0), label = "Holes Average Gate", color = 'red')

legend(loc = "lower left")
xlabel("Z-Dimension (microns)", fontsize=12)
#ylabel('$\\rho(x,y,z)/\epsilon_{Si}$ [V/um$^2$]',fontsize=9)
ylim(-50.0, 50.0)
xlim(0.0,10.0)
savefig(outputfiledir+"/plots/"+outputfilebase+"_1D_Fixed_Regions_%d.pdf"%run)


print "Making summary plots\n"
nxmin = 0
nxmax = nxx
nymin = 0
nymax = nyy

figure()
suptitle("Potentials on Bottom. Grid = %d*%d*%d."%(nxx,nyy,nzz),fontsize = 18)
[yy,xx] = meshgrid(dat.y[nymin:nymax],dat.x[nxmin:nxmax])

slicez = 0
subplot(1,1,1, aspect = 1)
levels = linspace(-60.0, 20.0, 161)
contour(xx,yy,dat.phi[nxmin:nxmax,nymin:nymax,slicez],levels,lw=0.1)
contourf(xx,yy,dat.phi[nxmin:nxmax,nymin:nymax,slicez],levels)
xlabel("X-Dimension (microns)")
ylabel("Y-Dimension (microns)")
cb=colorbar(orientation='horizontal')
cb.set_ticks([-60.0,-45.0,-30.0,-15.0,0.0,15.0,30.0])
cb.set_label('$\phi(x,y)$ [V]')
savefig(outputfiledir+"/plots/"+outputfilebase+"_Summary_1_%d.pdf"%run)


