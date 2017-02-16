#!/usr/bin/env python

#Author: Craig Lage, UC Davis; 
#Date: 16-Feb-17

#This program builds a fits file from the Poisson_CCD output files
from pylab import *
import sys, time, subprocess, h5py
import pyfits as pf

#****************SUBROUTINES*****************

class Array2dSet:
    def __init__(self,xmin,xmax,nx,ymin,ymax,ny,nstamps):
        # This packages up a set of nstamps postage stamp images,
        # each image of which is nx * ny pixels
        self.nx=nx
        self.ny=ny
        self.nstamps=nstamps

        self.xmin=xmin
        self.ymin=ymin
        
        self.xmax=xmax
        self.ymax=ymax
        
        self.dx=(xmax-xmin)/nx
        self.dy=(ymax-ymin)/ny
        
        self.x=linspace(xmin+self.dx/2,xmax-self.dx/2,nx)
        self.y=linspace(ymin+self.dy/2,ymax-self.dy/2,ny)

        self.data=zeros([nx,ny,nstamps], dtype = int32)
        self.xoffset=zeros([nstamps])
        self.yoffset=zeros([nstamps])
        self.imax=zeros([nstamps])

class Array3dHDF5Elec(object):
    def __init__(self, dir, filebase, n):
        elecfile = dir+'/'+filebase+'_'+str(n)+'_Elec.hdf5'
        hdfelec = h5py.File(elecfile,'r')
        Dimension = hdfelec[hdfelec.items()[0][0]].attrs[u'Dimension']
        self.nx=Dimension[0]
        self.ny=Dimension[1]
        self.nz=Dimension[2]
        
        Lower_Left = hdfelec[hdfelec.items()[0][0]].attrs[u'Lower_Left']
        self.xmin=Lower_Left[0]
        self.ymin=Lower_Left[1]
        self.zmin=Lower_Left[2]

        Upper_Right = hdfelec[hdfelec.items()[0][0]].attrs[u'Upper_Right']
        self.xmax=Upper_Right[0]
        self.ymax=Upper_Right[1]
        self.zmax=Upper_Right[2]
        
        self.dx=(self.xmax-self.xmin)/self.nx
        self.dy=(self.ymax-self.ymin)/self.ny
        self.dz=(self.zmax-self.zmin)/self.nz
        self.volume = self.dx * self.dy * self.dz
        
        self.x=linspace(self.xmin+self.dx/2,self.xmax-self.dx/2,self.nx)
        self.y=linspace(self.ymin+self.dy/2,self.ymax-self.dy/2,self.ny)
        self.z=linspace(self.zmin+self.dz/2,self.zmax-self.dz/2,self.nz)

        self.elec=array(hdfelec[hdfelec.items()[0][0]])

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

def FillSpotlist(ConfigData, run):
    # Note sigmas and offset are in pixels, not microns.
    # Postage Stamp size
    nx = ConfigData['PixelBoundaryNx']
    ny = ConfigData['PixelBoundaryNy']

    outputfiledir = ConfigData['outputfiledir']
    outputfilebase = ConfigData['outputfilebase']
    GridsPerPixel = ConfigData['GridsPerPixel'] * ConfigData['ScaleFactor']
    PixelSize = ConfigData['PixelSize']
    stampxmin = -(int(nx/2)+0.5)
    stampxmax = -stampxmin
    stampymin = -(int(ny/2)+0.5)
    stampymax = -stampymin

    spotlist = Array2dSet(stampxmin,stampxmax,nx,stampymin,stampymax,ny,1)

    dat = Array3dHDF5Elec(outputfiledir, outputfilebase, run)
    
    for i in range(nx):
        nxmin = int((ConfigData['PixelBoundaryLowerLeft'][0] - dat.xmin) / dat.dx) + GridsPerPixel * i
        nxmax = nxmin + GridsPerPixel
        for j in range(ny):
            nymin = int((ConfigData['PixelBoundaryLowerLeft'][1] - dat.ymin) / dat.dy) + GridsPerPixel * j
            nymax = nymin + GridsPerPixel
            electrons_in_pixel = dat.elec[nxmin:nxmax,nymin:nymax,:].sum()
            spotlist.data[i,j,0] = int(electrons_in_pixel)
    return spotlist

#****************MAIN PROGRAM*****************

# First, read the .cfg file

configfile = sys.argv[1]
run = int(sys.argv[2])
fitsfile = sys.argv[3]
cfg_success, ConfigData = ReadConfigFile(configfile)
if not cfg_success:
    print "Configuration file issue. Quitting"
    sys.exit()

spotlist = FillSpotlist(ConfigData, run)

primary_hdu = pf.PrimaryHDU(data=spotlist.data[:,:,0])
hdulist = pf.HDUList([primary_hdu])
hdulist.writeto(fitsfile, clobber=True)
