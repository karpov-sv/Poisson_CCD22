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
from scipy import stats

#****************MAIN PROGRAM*****************

nfiles = int(sys.argv[1])
# First, read the .cfg file

configfile = sys.argv[2]
run = int(sys.argv[3])
ConfigData = ReadConfigFile(configfile)

mainfiledir = ConfigData["outputfiledir"]

figure()

for n in range(nfiles):

    # First, read the .cfg file
    configfile = 'data/satrun_%d/sat.cfg'%n
    run = 0
    ConfigData = ReadConfigFile(configfile)
    outputfilebase = ConfigData["outputfilebase"]
    outputfiledir = ConfigData["outputfiledir"]
    file = open(outputfiledir+'/barrier.txt', 'w')
    file.write("   Flux      Barrier      \n")
    # Create the output directory if it doesn't exist
    if not os.path.isdir(outputfiledir+"/plots"):
        os.mkdir(outputfiledir+"/plots")

    # This holds all of the data
    dat = Array3dHDF5(outputfiledir, outputfilebase, run)

    ScaleFactor = ConfigData["ScaleFactor"]
    GridsPerPixel = ConfigData["GridsPerPixelX"]
    Vph = float(ConfigData["Vparallel_hi"])    
    Vpl = float(ConfigData["Vparallel_lo"])    
    qfh = float(ConfigData["qfh"])    
    dn = GridsPerPixel*ScaleFactor
    nxx = dat.nx - 1
    nyy = dat.ny - 1
    nzz = dat.nz - 1

    nxc = nxx/2
    nyc = nyy/2

    kmin = 4

    centers = [(nxc-3*dn,nyc+3*dn),(nxc,nyc+3*dn),(nxc+3*dn,nyc+3*dn),(nxc-3*dn,nyc-3*dn),(nxc,nyc-3*dn),(nxc+3*dn,nyc-3*dn)]
    fluxes = []
    for kk in range(6):
        fluxes.append(ChargeDepth(dat,centers[kk][0], centers[kk][1], dn/2, dn/2, dat.elec.shape[2]))
    bheights = []
    linefluxes = []
    linebheights = []

    for m, (nxcenter,nycenter) in enumerate(centers):
        elecmax = 0
        for i in range(dat.elec.shape[2]):
            if dat.elec[nxcenter,nycenter,i] > 0.0 and i > elecmax: elecmax = i

        bheight = 20.0
        #print nxcenter, nycenter, elecmax    
        for k in range(4, elecmax):
            if dat.elec[nxcenter, nycenter,k] < 1.0:
                continue

            beginphi = 0.0
            minphi = 20.0
            endphi = 0.0
            ReachedEnd = False
            for i in range(nycenter,nycenter+dn):
                if dat.elec[nxcenter,i,k] < 1.0 and ReachedEnd:
                    if dat.phi[nxcenter,i,k] < minphi:
                        minphi = dat.phi[nxcenter,i,k]
                        imin = i
                    if dat.elec[nxcenter,i+1,k] > 1.0:
                        endphi = dat.phi[nxcenter,i+1,k]
                        iend = i
                        break
                if dat.elec[nxcenter,i,k] > 1.0 and not ReachedEnd:
                    if dat.elec[nxcenter,i+1,k] < 1.0:
                        beginphi = dat.phi[nxcenter,i,k]
                        ibegin = i
                        ReachedEnd = True
            if ReachedEnd == False: bheight = 0.0
            bheight = min(bheight,min(abs(beginphi - minphi),abs(endphi-minphi)))
            #print nxcenter, nycenter, beginphi, minphi,endphi, ibegin,imin,iend
            beginphi = 0.0
            minphi = 20.0
            endphi = 0.0
            ReachedEnd = False
            for i in range(nycenter-dn,nycenter):
                if dat.elec[nxcenter,i,k] < 1.0 and ReachedEnd:
                    if dat.phi[nxcenter,i,k] < minphi:
                        minphi = dat.phi[nxcenter,i,k]
                    if dat.elec[nxcenter,i+1,k] > 1.0:
                        endphi = dat.phi[nxcenter,i+1,k]
                        break
                if dat.elec[nxcenter,i,k] > 1.0 and not ReachedEnd:
                    if dat.elec[nxcenter,i+1,k] < 1.0:
                        beginphi = dat.phi[nxcenter,i,k]
                        ReachedEnd = True
            if ReachedEnd == False: bheight = 0.0
            bheight = min(bheight,min(abs(beginphi - minphi),abs(endphi-minphi)))
            #print nxcenter, nycenter, beginphi, minphi,endphi, ibegin,imin,iend
        bheights.append(bheight)
        if bheight > 0:
            linebheights.append(bheight * 1000.0)
            linefluxes.append(fluxes[m])
    slope, intercept, r_value, p_value, std_err = stats.linregress(linefluxes, linebheights)
    value = (75.0 - intercept) / slope
    plot(fluxes, array(bheights)*1000.0, label = "Vph = %.1f"%Vph)

    for kk in range(6):
        file.write(" %.1f     %.4f\n"%(fluxes[kk], bheights[kk]))
    file.write("Flux at 5kT/q = %f\n"%value)
    file.close()


title("Interwell Barrier Height vs Pixel Charge")
plot([0,300000],[75.0, 75.0], color = 'blue', ls = '--')
text(200000, 85.0,"5*kT/q")
xlabel("Pixel Charge (e-)")
ylabel("Barrier Height (mV)")
legend()
savefig(outputfiledir+"/plots/Barrier_Height_Vpl_%.1f_Vph_%.1f.pdf"%(Vpl,Vph))

