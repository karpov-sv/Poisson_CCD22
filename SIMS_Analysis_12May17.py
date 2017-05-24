#!/usr/bin/env python

#Author: Craig Lage, NYU;
#Date: 10-Nov-16

#This program plots the Poisson equation solutions from the C++ Poisson solver
import matplotlib
matplotlib.use("PDF")
from pylab import *
import os, sys, time, h5py, xlrd

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

def ReadSIMSData():

    epsilon_si = 11.7 * 8.85E-14
    qe = 1.6E-19
    boron_wb = xlrd.open_workbook('C0HVL528x02_B site.xls')
    boron_data = boron_wb.sheet_by_name('Processed data')
    boron_conc = []
    boron_depth = []
    for i in range(boron_data.nrows):
        try:
            if type(boron_data.row(i)[0].value) is float:
                if boron_data.row(i)[1].value < 1.0E17: # One data point is garbage
                    boron_depth.append(boron_data.row(i)[0].value)
                    boron_conc.append(boron_data.row(i)[1].value)
        except:
            continue

    phos_wb = xlrd.open_workbook('C0HVL528L05_P Site.xls')
    phos_data = phos_wb.sheet_by_name('Processed data')
    phos_conc = []
    phos_depth = []
    for i in range(phos_data.nrows):
        try:
            if type(phos_data.row(i)[0].value) is float:
                phos_depth.append(phos_data.row(i)[0].value)
                phos_conc.append(phos_data.row(i)[1].value)
        except:
            continue

    boron_conc = -array(boron_conc) * qe / epsilon_si * 1.0E-8 # Convert to V/um/2
    phos_conc = array(phos_conc) * qe / epsilon_si * 1.0E-8 # Convert to V/um/2
    boron_depth = array(boron_depth)
    phos_depth = array(phos_depth)    
    return [boron_depth, boron_conc, phos_depth, phos_conc]


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


print "Making 1D potential and Charge Density plots\n"
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
savefig(outputfiledir+"/plots/"+outputfilebase+"_%d_SIMS_Comparison.pdf"%run)



