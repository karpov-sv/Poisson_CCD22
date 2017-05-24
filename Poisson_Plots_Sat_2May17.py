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

if EdgePlot:
    NumPixelsPlotted = 8 # Special for edge plots
    nycenter2 = nycenter - 16 * GridsPerPixelX * ScaleFactor
    nymax = nycenter + 21 * ScaleFactor * GridsPerPixelY
    nymin = nycenter - 21 * ScaleFactor * GridsPerPixelY

else:
    NumPixelsPlottedX = 3
    NumPixelsPlottedY = 4    
    nycenter2 = nycenter
    nymin = nycenter - (NumPixelsPlottedY * ScaleFactor * GridsPerPixelY)/2
    nymax = nycenter + (NumPixelsPlottedY * ScaleFactor * GridsPerPixelY)/2

nxmin = nxcenter - (NumPixelsPlottedX * ScaleFactor * GridsPerPixelX)/2
nxmax = nxcenter + (NumPixelsPlottedX * ScaleFactor * GridsPerPixelX)/2

nzmin = 0
nzmax = 16 * ScaleFactor

print "Channel CG"
for i in range(dat.elec.shape[2]):
    print "k = %d, z=%.3f, phi = %.3f, holes = %.3f, elec = %.3f, rho = %.3f, Ez = %.3f"%(i,dat.z[i],dat.phi[nxcenter, nycenter,i],dat.hole[nxcenter, nycenter,i],dat.elec[nxcenter, nycenter,i],dat.rho[nxcenter, nycenter,i],dat.Ez[nxcenter, nycenter,i])

print "Channel BG"
for i in range(dat.elec.shape[2]):
    print "k = %d, z=%.3f, phi = %.3f, holes = %.3f, elec = %.3f, rho = %.3f, Ez = %.3f"%(i,dat.z[i],dat.phi[nxcenter, nycenter+ScaleFactor*GridsPerPixelY/2,i],dat.hole[nxcenter, nycenter+ScaleFactor*GridsPerPixelY/2,i],dat.elec[nxcenter, nycenter+ScaleFactor*GridsPerPixelY/2,i],dat.rho[nxcenter, nycenter+ScaleFactor*GridsPerPixelY/2,i],dat.Ez[nxcenter, nycenter+ScaleFactor*GridsPerPixelY/2,i])


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
    plot(dat.z[0:elecnumzs], -ChargeFactor * (dat.elec[nxcenter2,nycenter2,0:elecnumzs]+dat.elec[nxcenter2-1,nycenter2,0:elecnumzs]+dat.elec[nxcenter2,nycenter2-1,0:elecnumzs]+dat.elec[nxcenter2-1,nycenter2-1,0:elecnumzs])/4.0 / dat.dz[0:elecnumzs], label = "%d"%fluxes[m])

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
ylim(-10.0, 10.0)
xlim(dat.y[nycenter-ddny]-dat.y[nycenter], dat.y[nycenter+ddny]-dat.y[nycenter])
xlabel("Y-Dimension (microns)")
ylabel("Potential(Volts)")
legend(loc = "lower right")


savefig(outputfiledir+"/plots/"+outputfilebase+"_1D_Multi_%d.pdf"%run)

