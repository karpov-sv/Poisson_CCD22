#!/usr/bin/env python

#Author: Craig Lage, NYU;
#Date: 16-Nov-15

#This program plots the pixel area plots from the Poisson CCD solver
import matplotlib
matplotlib.use("PDF")
from pylab import *
import sys, time, h5py

#****************SUBROUTINES*****************

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

def ReadAreaFile(filename, nx, ny):
    area = zeros([nx, ny])
    file = open(filename, 'r')
    lines = file.readlines()
    file.close()
    for line in lines:
        items = line.split()
        if items[0] == "Nx":
            continue
        i = int(items[0])
        j = int(items[1])
        area[i,j] = float(items[2])
    return area

def ReadVertexFile(filename, nx, ny, NumAngles):
    vx = zeros([nx, ny, NumAngles])
    vy = zeros([nx, ny, NumAngles])
    file = open(filename, 'r')
    lines = file.readlines()
    file.close()
    FirstPass = True
    for line in lines:
        items = line.split()
        if items[0] == "X0":
            continue
        x = int(float(items[0]))
        y = int(float(items[1]))
        if FirstPass:
            lastx = x; lasty = y; i = 0; j = 0; k = 0; FirstPass = False
        if y != lasty:
            k = 0
            j += 1
            lasty = y
        if x != lastx:
            k = 0
            j = 0
            i += 1
            lastx = x
            lasty = y

        #print i,j,k
        vx[i,j,k] = float(items[3])
        vy[i,j,k] = float(items[4])
        k += 1
    return (vx, vy)


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

area = ReadAreaFile(filename, Nx, Ny)
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
