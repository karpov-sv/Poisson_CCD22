#!/usr/bin/env python

#Author: Craig Lage, UC Davis; 
#Date: 3-Mar-17

# This program converts the GalSim PhotonArray fits file to a text file
# that can be read into Poisson_CCD
import sys
import pyfits as pf

#****************SUBROUTINES*****************

def ConvertPhotonFile(fitsfilename):
    pyfits_version = float(pf.version.__version__)
    datfilename = fitsfilename.strip('.fits') + '.dat'

    # Read the fits file
    with pf.open(fitsfilename) as fits:
        data = fits[1].data
    N = len(data)
    if pyfits_version > '3.0':
        names = data.columns.names
    else: # pragma: no cover
        names = data.dtype.names
    x = data['x']
    y = data['y']
    flux = data['flux']
    if 'dxdz' in names:
        dxdz = data['dxdz']
        dydz = data['dydz']
    if 'wavelength' in names:
        wavelength = data['wavelength']

    # Write the dat file
    outfile = open(datfilename, 'w')
    outfile.write('ID \t X(pixels) \t Y(pixels) \t dxdz      \t dydz      \t lambda(nm)\n')
    for i in range(N):
        xi = x[i] # in pixels
        yi = y[i] #in pixels

        if 'dxdz' in names:
            dxdzi = dxdz[i]
            dydzi = dydz[i]
        else:
            dxdzi = 0.0
            dydzi = 0.0

        if 'wavelength' in names:
            wavelengthi = wavelength[i]
        else:
            wavelengthi = 0.0

        outfile.write('%d \t %.6f \t %.6f \t %.6f \t %.6f \t %.6f\n'%(i,xi,yi,dxdzi,dydzi,wavelengthi))
    outfile.close()
    return 

#****************MAIN PROGRAM*****************

# First, read the .fits file

fitsfilename = sys.argv[1]
ConvertPhotonFile(fitsfilename)
