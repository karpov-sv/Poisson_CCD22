#!/usr/bin/env python

#Author: Craig Lage, NYU;
#Date: 26-Jan-15

#This program plots the Poisson equation solutions from the C++ Poisson solver
import matplotlib
matplotlib.use("PDF")
from pylab import *
from matplotlib.colors import from_levels_and_colors
import os, sys, time, h5py

#****************SUBROUTINES*****************
class Array3dHDF5Elec(object):
    def __init__(self, dir, filebase, n):
        elecfile = dir+'/'+filebase+'_'+str(n)+'_Elec' + '.hdf5'
        hdfelec = h5py.File(elecfile,'r')
        holefile = dir+'/'+filebase+'_'+str(n)+'_Hole' + '.hdf5'
        hdfhole = h5py.File(holefile,'r')
        self.elec=array(hdfelec[hdfelec.items()[0][0]])
        self.hole=array(hdfhole[hdfhole.items()[0][0]])
        #self.comb = log10(self.hole.clip(1.0,1.0E30)) - log10(self.elec.clip(1.0,1.0E30))
        phifile = dir+'/'+filebase+'_'+str(run)+'_phi' + '.hdf5'
        rhofile = dir+'/'+filebase+'_'+str(run)+'_rho' + '.hdf5'
        hdfrho = h5py.File(rhofile,'r')
        self.rho=array(hdfrho[hdfrho.items()[0][0]])
        hdfphi = h5py.File(phifile,'r')
        self.phi=array(hdfphi[hdfphi.items()[0][0]])        
        xfile = dir+'/'+'grid_x.dat'
        yfile = dir+'/'+'grid_y.dat'
        zfile = dir+'/'+'grid_z.dat'        

        xgrid = loadtxt(xfile, skiprows=1)
        ygrid = loadtxt(yfile, skiprows=1)
        zgrid = loadtxt(zfile, skiprows=1)        

        self.nx=xgrid.shape[0]
        self.ny=ygrid.shape[0]

        Dimension = hdfelec[hdfelec.items()[0][0]].attrs[u'Dimension']
        self.nz=Dimension[2] # Note the Elec and Hole grids are smaller than the phi,etc. grids.

        self.xmin=xgrid[0,1]
        self.ymin=ygrid[0,1]
        self.zmin=zgrid[0,1]

        self.xmax=xgrid[self.nx-1,3]
        self.ymax=ygrid[self.ny-1,3]
        self.zmax=zgrid[self.nz-1,3]

        self.x=xgrid[:,2]
        self.y=ygrid[:,2]
        self.z=zgrid[0:self.nz,2]

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

def ChargeDepth(filename, nxcenter, nycenter, ScaleFactor, GridsPerPixelX, GridsPerPixelY):
    file = open(filename,"w")
    nzxmin = nxcenter - (ScaleFactor * GridsPerPixelX)/2
    nzxmax = nzxmin + (ScaleFactor * GridsPerPixelX)
    nzymin = nycenter - (ScaleFactor * GridsPerPixelY)/2
    nzymax = nzymin + (ScaleFactor * GridsPerPixelY)
    ncenter = 0.0
    nzcenter = 0.0

    Total_elec = 0.0
    Below_kmin = 0.0
    for nx in range(nzxmin, nzxmax):
        for ny in range(nzymin, nzymax):
            for nz in range(nzz):
                if dat.elec[nx,ny,nz] > 0.0:
                    ncenter += dat.elec[nx,ny,nz]
                    nzcenter += dat.z[nz] * dat.elec[nx,ny,nz]
    if ncenter > 0:
        meanz = nzcenter / ncenter
    else:
        meanz = 0.0
    print "Electrons in center Pixel = %.1f, Mean z = %.3f microns\n"%(ncenter, meanz)
    print "Total Electrons = %.1f\n"%dat.elec.sum()
    file.write("Electrons in center Pixel = %.1f, Mean z = %.3f microns\n"%(ncenter, meanz))
    file.close()

def BuildPlotArray(plotdata, axis, nxmin, nxmax, nymin, nymax, nzmin, nzmax, ForceZero, cmap):
    plotarray = ma.masked_array(zeros([nxmax-nxmin,nymax-nymin]))
    dxx = zeros([nxmax-nxmin])
    dyy = zeros([nymax-nymin])
    for i in range(nxmax-nxmin):
        if axis == 0:
            dxx[i] = dat.z[nxmin+i]
        elif axis == 1:
            dxx[i] = dat.x[nxmin+i]            
        elif axis == 2:             
            dxx[i] = dat.x[nxmin+i]            
    for j in range(nymax-nymin):
        if axis == 0:
            dyy[j] = dat.y[nymin+j]
        elif axis == 1:
            dyy[j] = dat.z[nymin+j]
        elif axis == 2:             
            dyy[j] = dat.y[nymin+j]            

    for i in range(nxmax-nxmin):
        for j in range(nymax-nymin):
            for k in range(nzmax-nzmin):
                if axis == 0: 
                    # Y-Z slice
                    plotarray[i,j] += plotdata[nzmin+k,nymin+j,nxmin+i]
                elif axis == 1: 
                    # X-Z slice
                    plotarray[i,j] += plotdata[nxmin+i,nzmin+k,nymin+j]                
                elif axis == 2: 
                    # X-Y slice
                    plotarray[i,j] += plotdata[nxmin+i,nymin+j,nzmin+k]                

    num_levels = 100
    if ForceZero:
        pdata = plotarray.clip(0.0,1.0E20)
        ndata = plotarray.clip(-1.0E30,0.0)
        plotarray = (pdata/pdata.max()-ndata/ndata.min())
    vmax = plotarray.max()*1.01 + 0.01
    vmin = plotarray.min()*1.01 - 0.01
    levels = np.linspace(vmin, vmax, num_levels)    

    for i in range(nxmax-nxmin):
        for j in range(nymax-nymin):
            if axis == 0:
                if abs(dat.rho[(nzmin+nzmax)/2,nymin+j,nxmin+i]) < 1.0E-12:
                    plotarray[i,j] = vmax+1.0E6
                if abs(Vph - dat.phi[(nzmin+nzmax)/2,nymin+j,nxmin+i]) < 1.0E-12 or abs(Vpl - dat.phi[(nzmin+nzmax)/2,nymin+j,nxmin+i]) < 1.0E-12:
                    plotarray[i,j] = vmin-1.0E6
            elif axis == 1:
                if abs(dat.rho[nxmin+i, (nzmin+nzmax)/2,nymin+j]) < 1.0E-12:
                    plotarray[i,j] = vmax+100000.0
                if abs(Vph - dat.phi[nxmin+i,(nzmin+nzmax)/2,nymin+j]) < 1.0E-12 or abs(Vpl - dat.phi[nxmin+i,(nzmin+nzmax)/2,nymin+j]) < 1.0E-12:
                    plotarray[i,j] = vmin-1.0E6
    if axis == 0:
        for j in range(nymax-nymin):
            plotarray[1,j] =  vmin-1.0E6
    if axis == 1:
        for i in range(nxmax-nxmin):
            plotarray[i,1] =  vmin-1.0E6

    cmap.set_under("green")
    cmap.set_over("yellow")
                
    [pyy,pxx] = meshgrid(dyy, dxx)# Data grid for plots
    return [plotarray, pxx, pyy, levels, cmap]


def BuildPlotSlice(plotdata, axis, nxmin, nxmax, nymin, nymax, nzmin, nzmax):
    plotslice = zeros([nxmax-nxmin])
    xpoints = zeros([nxmax-nxmin])
    for i in range(nxmax-nxmin):
        if axis == 0: 
            xpoints[i] = dat.x[nxmin+i]            
        elif axis == 1:
            xpoints[i] = dat.y[nxmin+i]            
        elif axis == 2:             
            xpoints[i] = dat.z[nxmin+i]

    for i in range(nxmax-nxmin):
        for j in range(nymax-nymin):
            for k in range(nzmax-nzmin):
                if axis == 0: 
                    # X Slice
                    plotslice[i] += plotdata[nxmin+i,nymin+j,nzmin+k]                
                elif axis == 1: 
                    # Y Slice
                    plotslice[i] += plotdata[nymin+j,nxmin+i,nzmin+k]                
                elif axis == 2: 
                    # Z Slice
                    plotslice[i] += plotdata[nymin+j,nzmin+k,nxmin+i]                
    return [plotslice, xpoints]

#****************MAIN PROGRAM*****************

# First, read the .cfg file

configfile = sys.argv[1]
run = int(sys.argv[2])
NumPixelsPlotted = int(sys.argv[3])
    
cfg_success, ConfigData = ReadConfigFile(configfile)
if not cfg_success:
    print "Configuration file issue. Quitting"
    sys.exit()
outputfilebase = ConfigData["outputfilebase"]
outputfiledir = ConfigData["outputfiledir"]
Vpl = ConfigData["Vparallel_lo"]
Vph = ConfigData["Vparallel_hi"]
# This holds all of the data
dat = Array3dHDF5Elec(outputfiledir, outputfilebase, run)

ScaleFactor = ConfigData["ScaleFactor"]
GridsPerPixelX = ConfigData["GridsPerPixelX"]
GridsPerPixelY = ConfigData["GridsPerPixelY"]

ZMult = 2.0
kmax = int(dat.nz*0.9)

nxx = dat.nx - 1
nyy = dat.ny - 1
nzz = dat.nz - 1
nxcenter = nxx/2
nycenter = nyy/2
dnx = ScaleFactor*GridsPerPixelX/2
dny = ScaleFactor*GridsPerPixelY/2
csxcenter = nxcenter + dnx
bgycenter = nycenter + dny
nxmin = nxcenter - NumPixelsPlotted * dnx
nxmax = nxmin + NumPixelsPlotted * 2 * dnx
nymin = nycenter - NumPixelsPlotted * dny
nymax = nymin + NumPixelsPlotted * 2 * dny

ChargeDepth(outputfiledir+'/charge.txt', nxcenter, nycenter, ScaleFactor, GridsPerPixelX, GridsPerPixelY)

carriers = ['Electron', 'Hole', 'Mobile', 'Fixed']
plotdatas = [dat.elec, dat.hole, dat.hole-dat.elec, dat.rho]
cmap0 = cm.get_cmap("jet")
cmap1 = cm.get_cmap("seismic")
cmaps = [cmap0, cmap0, cmap1, cmap1]
ForceZeros = [False, False, True, True]
xslicemins = [nxcenter-dnx, csxcenter-1, nxcenter-dnx, nxcenter-dnx]
xslicemaxs = [nxcenter+dnx, csxcenter+2, nxcenter+dnx, nxcenter+dnx]
# Create the output directory if it doesn't exist
if not os.path.isdir(outputfiledir+"/plots"):
    os.mkdir(outputfiledir+"/plots")

for i, plotdata in enumerate(plotdatas):

    fig = figure(figsize = (12,12))
    suptitle("%s Charge Distribution"%carriers[i], fontsize = 36)

    ax1=axes([0.10,0.40,0.50,0.50],aspect=1)
    ax1.set_title("X-Y Slice")
    ax1.set_xticks([])
    ax1.set_yticks([])
    [plotarray, dxx, dyy, levels, my_cmap] = BuildPlotArray(plotdata, 2, nxmin, nxmax, nymin, nymax, 0, kmax, ForceZeros[i], cmaps[i])
    ax1.contourf(dxx, dyy, plotarray, levels = levels, cmap = my_cmap, extend='both')

    ax2=axes([0.10,0.20,0.50,0.20], aspect=ZMult)
    ax2.set_title("X-Z Slice")
    ax2.set_xticks([])
    ax2.set_yticks([0.0,1.0,2.0])
    ax2.set_ylabel("Z (Microns)")
    [plotarray, dxx, dyy, levels, my_cmap] = BuildPlotArray(plotdata, 1, nxmin, nxmax, 0, kmax, nymin, nymax, ForceZeros[i], cmaps[i])
    ax2.contourf(dxx, dyy, plotarray, levels = levels, cmap = my_cmap, extend='both')

    ax3=axes([0.10,0.10,0.50,0.10])
    ax3.set_title("X-Cut, Y = %.2f"%dat.y[nycenter])
    ax3.set_xlabel("X (Microns)")
    ax3.set_ylabel("Charge Density \n(arb. units)")
    [plotslice, xpoints] = BuildPlotSlice(plotdata, 0, nxmin, nxmax, nycenter-1, nycenter+2, 0, kmax)
    ax3.plot(xpoints, plotslice)
    
    ax4=axes([0.60,0.40,0.20,0.50], aspect=1.0/ZMult)
    ax4.set_title("Y-Z Slice")
    ax4.set_xticks([0.0,1.0,2.0])
    ax4.set_yticks([])
    ax4.set_xlabel("Z (Microns)")
    [plotarray, dxx, dyy, levels, my_cmap] = BuildPlotArray(plotdata, 0, 0, kmax, nymin, nymax, xslicemins[i], xslicemaxs[i], ForceZeros[i], cmaps[i])
    ax4.contourf(dxx, dyy, plotarray, levels = levels, cmap = my_cmap, extend='both')        
    ax4.invert_xaxis()

    ax5=axes([0.80,0.40,0.10,0.50])
    ax5.set_title("Y-Cut, X = %.2f"%dat.x[nxcenter])
    ax5.yaxis.tick_right()
    ax5.yaxis.set_label_position("right")
    for tick in ax5.get_xticklabels():
        tick.set_rotation(90)
    ax5.set_ylabel("Y (Microns)")
    ax5.set_xlabel("Charge Density \n(arb. units)")
    [plotslice, xpoints] = BuildPlotSlice(plotdata, 1, nymin, nymax, nxcenter-1, nxcenter+2, 0, kmax)
    ax5.plot(plotslice, xpoints)
    ax5.set_xlim(ax5.get_xlim()[::-1])

    ax6=axes([0.65,0.20,0.10,0.10])
    ax6.set_title("Z-Cut")
    ax6.set_ylabel("Log Charge Density \n(arb. units)")
    ax6.set_xticks([0.0,1.0,2.0])
    ax6.set_xlabel("Z ( Microns)")
    ax6.yaxis.set_label_position("right")
    ax6.set_ylim(-4.0, 1.0)
    ax6.text(3.0,-8.0,"Z-Axis has a %.1fX Scale Multiplier"%ZMult, fontsize = 16)
    if i == 2:
        [plotslicen, xpoints] = BuildPlotSlice(plotdata, 2, 0, kmax, nxcenter-1, nxcenter+2, nymin, nymax)
        [plotslicep, xpoints] = BuildPlotSlice(plotdata, 2, 0, kmax, csxcenter-1, csxcenter+2, nymin, nymax)        
        ax6.plot(xpoints, log10(plotslicen/plotslicen.min()+1.0E-5), color='blue')
        ax6.plot(xpoints, log10(plotslicep/plotslicep.max()+1.0E-5), color='red')        
        ax6.text(3.0,-9.0,"+ Charge in Red, - Charge in Blue", fontsize = 16)        
    elif i == 3:
        [plotslicep, xpoints] = BuildPlotSlice(plotdata, 2, 0, kmax, nxcenter-1, nxcenter+2, nymin, nymax)
        [plotslicen, xpoints] = BuildPlotSlice(plotdata, 2, 0, kmax, csxcenter-1, csxcenter+2, nymin, nymax)        
        ax6.plot(xpoints, log10(plotslicen/plotslicen.min()+1.0E-5), color='blue')
        ax6.plot(xpoints, log10(plotslicep/plotslicep.max()+1.0E-5), color='red')        
        ax6.text(3.0,-9.0,"+ Charge in Red, - Charge in Blue", fontsize = 16)        
    else:
        [plotslice, xpoints] = BuildPlotSlice(plotdata, 2, 0, kmax, nxmin, nxmax, nymin, nymax)
        ax6.plot(xpoints, log10(plotslice/plotslice.max()+1.0E-5))
    ax6.set_xlim(ax6.get_xlim()[::-1])
    ax6.text(3.0,-10.0,"Oxide in yellow, Gates in green", fontsize = 16)        
    savefig(outputfiledir+"/plots/%sDistribution_27Jan17_%d.pdf"%(carriers[i],run))
    close(fig)
