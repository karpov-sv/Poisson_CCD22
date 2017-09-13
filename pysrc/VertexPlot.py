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
Nx = ConfigData["PixelBoundaryNx"]
Ny = ConfigData["PixelBoundaryNy"]
XCenter = ConfigData["FilledPixelCoords_0_0"][0]
YCenter = ConfigData["FilledPixelCoords_0_0"][1]
PixelSizeX = ConfigData["PixelSizeX"]
PixelSizeY = ConfigData["PixelSizeY"]
NxCenter = int((XCenter - ConfigData["PixelBoundaryLowerLeft"][0]) / PixelSizeX)
NyCenter = int((YCenter - ConfigData["PixelBoundaryLowerLeft"][1]) / PixelSizeY)
Area_0 = 100.0
NumAngles = 4 * ConfigData["NumVertices"] + 4
NumElec = ConfigData["CollectedCharge_0_0"]

PlotDelta = int(sys.argv[3])

filename = outputfiledir + '/' + outputfilebase +'_%d_Area'%run + '.dat'

[area,sim] = ReadAreaFile(filename, Nx, Ny, NxCenter, NyCenter, Area_0)

figure()
subplot(1,1,1,aspect = 1)
title("Pixel Area: %d e-"%NumElec)

filename = outputfiledir + '/' + outputfilebase +'_%d_Vertices'%run + '.dat'

(vx, vy) = ReadVertexFile(filename, Nx, Ny, NumAngles)

LineXMin = XCenter - (PlotDelta + 0.5) * PixelSizeX
LineXMax = XCenter + (PlotDelta + 0.5) * PixelSizeX
LineYMin = YCenter - (PlotDelta + 0.5) * PixelSizeY
LineYMax = YCenter + (PlotDelta + 0.5) * PixelSizeY

title("Pixel Vertices: %d e-"%NumElec)

for i in range(NxCenter-PlotDelta, NxCenter+PlotDelta+1):
    plot([XCenter+PixelSizeX*(i-NxCenter-0.5), XCenter+PixelSizeX*(i-NxCenter-0.5)], [LineYMin, LineYMax], color = 'black', ls = '--')
for j in range(NyCenter-PlotDelta, NyCenter+PlotDelta+1):
    plot([LineXMin, LineXMax], [YCenter+PixelSizeY*(j-NyCenter-0.5), YCenter+PixelSizeY*(j-NyCenter-0.5)], color = 'black', ls = '--')


for i in range(NxCenter-PlotDelta, NxCenter+PlotDelta+1):
    for j in range(NyCenter-PlotDelta, NyCenter+PlotDelta+1):
        if i == NxCenter and j == NyCenter:
            textcolor = 'red'
        else:
            textcolor = 'black'
        text(XCenter+PixelSizeX*(i-NxCenter-0.2), YCenter+PixelSizeY*(j-NyCenter-0.1), "%.4f"%area[i,j], color = textcolor, fontsize = 12/PlotDelta, fontweight = 'bold')
        x = []
        y = []
        for k in range(NumAngles):
            x.append(vx[i,j,k])
            y.append(vy[i,j,k])
        x.append(vx[i,j,0])
        y.append(vy[i,j,0])
        plot(x, y, lw = 0.5)

savefig(outputfiledir+"/plots/PixelVertices_%d.pdf"%run)
