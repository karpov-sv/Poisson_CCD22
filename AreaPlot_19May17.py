#!/usr/bin/env python

#Author: Craig Lage, UC Davis
#Date: 19-May-17

#This program plots the covariance plots simulations against measured data
import matplotlib
matplotlib.use("Agg")
from pylab import *
import sys, time, h5py

#****************SUBROUTINES*****************

def ReadConfigFile(filename):
    # This reads the config file for the necessary settings
    # and returns a dictionary with the values
    file = open(filename,'r')
    lines=file.readlines()
    file.close()
    ConfigData = {}
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
                        try: ConfigData[ParamName] = float(ThisParam)
                        except ValueError:
                            try:
                                ConfigData[ParamName] = ThisParam
                            except ValueError:
                                print "Error reading .cfg file"
                else:
                    ThisParam = []
                    for item in ThisLine:
                        try: ThisParam.append(int(item))
                        except ValueError:
                            try: ThisParam.append(float(item))
                            except ValueError:
                                ThisParam.append(item)
                    ConfigData[ParamName] = ThisParam
            except:
                continue
    except:
        print "Error reading .cfg file"

    return ConfigData

def ReadCorrData(filename):
    # This reads the correlation data file
    # and returns arrays with the values and sigmas
    systematic_error = 0.0001
    data = zeros([6,6])
    sigma = zeros([6,6])
    file = open(filename,'r')
    lines=file.readlines()
    file.close()
    lines.remove(lines[0]) # Strip the title line
    try:
        for line in lines:
            items = line.split()
            if items[0] == 'Slope':
                break
            i = int(items[0])
            j = int(items[1])        
            data[i,j] = float(items[2])
            sigma[i,j] = float(items[3]) + systematic_error        
    except:
        print "Error reading data file"
        sys.exit()
    return [data, sigma]

def ReadAreaFile(filename, nx, ny, nxcenter, nycenter, Area_0):
    # This reads the correlation data file
    # and returns an array with the expected correlations
    area = zeros([nx, ny])
    file = open(filename, 'r')
    lines = file.readlines()
    file.close()
    lines.remove(lines[0]) # Strip the title line    
    for line in lines:
        items = line.split()
        i = int(items[0])
        j = int(items[1])
        area[i,j] = float(items[2])

    sim = zeros([6,6])
    num = zeros([6,6], dtype = int)    
    for i in range(Nx):
        for j in range(Ny):
            ii = abs(i - NxCenter)
            jj = abs(j - NxCenter)
            sim[jj,ii] += (area[i,j] - Area_0) / Area_0
            num[jj,ii] += 1
    for i in range(6):
        for j in range(6):
            if num[i,j] >0:
                sim[i,j] /= float(num[i,j])
    return [area,sim]


#****************MAIN PROGRAM*****************

# First, read the .cfg file

configfile = sys.argv[1]
step = int(sys.argv[2])
ConfigData = ReadConfigFile(configfile)
outputfilebase = ConfigData["outputfilebase"]
outputfiledir = ConfigData["outputfiledir"]
Nx = ConfigData["PixelBoundaryNx"]
Ny = ConfigData["PixelBoundaryNy"]
NxCenter = 4
NyCenter = 4
Area_0 = 99.9972

areafilename = outputfiledir + '/' + outputfilebase +'_%d_Area'%step+'.dat'
datafilename = outputfiledir + '/corr_meas.txt'
[area,sim] = ReadAreaFile(areafilename, Nx, Ny, NxCenter, NyCenter, Area_0)
[data,sigma] = ReadCorrData(datafilename)

NumElec = 80000#step * 1000

x = []
y = []
xfit = []
yfit = []
xneg = []
yneg = []
xdata = []
ydata = []
yerr = []
ysigma = []
xdataneg = []
ydataneg = []
yerrneg = []
ysigmaneg = []

FOM = 0.0
Num_FOM = 0.0

fullfile = open(outputfiledir+"/full_corr.txt","w")
for i in range(5):
    for j in range(5):
        rsquared = float(i**2 + j**2)
        if i+j == 0 or rsquared > 18:
            continue
        FOM_contrib = ((data[i,j] - sim[i,j]) / sigma[i,j])**2
        print "i = %d, j = %d, data = %.4f, sim = %.4f, FOM_contrib = %.4f"%(i,j,data[i,j],sim[i,j],FOM_contrib)
        FOM += FOM_contrib
        Num_FOM += 1.0
        yvalue = sim[i,j]
        if yvalue < 0:
            xneg.append(rsquared)
            yneg.append(-yvalue)
        else:
            x.append(rsquared)
            y.append(yvalue)
        yvalue = data[i,j]
        if yvalue < 0:
            xdataneg.append(rsquared)
            ydataneg.append(-yvalue)
            yerrneg.append(sigma[i,j])
        else:
            xdata.append(rsquared)
            ydata.append(yvalue)
            yerr.append(sigma[i,j])
        if rsquared > 1.1 and i < 3 and j < 3 and yvalue > 0:
            fullfile.write("i = %d, j = %d, C = %.6f\n"%(i,j,yvalue))
            xfit.append(rsquared)
            yfit.append(yvalue)


chisquared = FOM / Num_FOM
print "Chi-squared = %.2f"%chisquared 
figure()
title("Covariance vs Distance: %d e-"%NumElec)
xscale('log')
yscale('log')
xlim(0.8, 100.0)
ylim(1E-5, 1E-1)
yerr = array(yerr)
ylower = np.maximum(1.1E-5, ydata - yerr)
yerr_lower = ydata - ylower
yerrneg = array(yerrneg)
ylowerneg = np.maximum(1.1E-5, ydataneg - yerrneg)
yerr_lowerneg = ydataneg - ylowerneg
errorbar(xdata,ydata, yerr = [yerr_lower, 2.0*yerr] , ls = 'None',marker = '.', ms = 10, color = 'green', label = 'Data')
scatter(array(x), array(y), marker = 'x', s = 50, color = 'blue', label = 'Sim')
if len(xdataneg) > 0:
    errorbar(xdataneg,ydataneg, yerr = [yerr_lowerneg, 2.0*yerrneg] , ls = 'None',marker = '.', ms = 10, color = 'cyan', label = 'Data-Neg')
if len(xneg) > 0:
    scatter(array(xneg), array(yneg), marker = 'x', s = 50, color = 'magenta', label = 'Sim-Neg')

from scipy import stats
slope, intercept, r_value, p_value, std_err = stats.linregress(log10(xfit),log10(yfit))
xplot=linspace(0.0, 2.0, 100)
yplot = slope * xplot + intercept
plot(10**xplot, 10**yplot, color='red', lw = 2, ls = '--')
text(6.0, 0.0316, "Slope = %.3f"%slope)
text(6.0, 0.0125, "C10 = %.5g"%sim[0,1])
text(6.0, 0.0079, "C01 = %.5g"%sim[1,0])
text(6.0, 0.005, "$\chi^2$ = %.2f"%chisquared)
xlabel("$i^2 + j^2$")
ylabel("$\delta$ Area / Area")
xlim(0.8,100)
legend()
savefig(outputfiledir+"/plots/Area_Corr_19May17_%d.pdf"%step)


file = open(outputfiledir+"/corr.txt","w")
file.write("Slope = %.3f\n"%slope)
file.write("C10 = %.5g\n"%sim[1,0])
file.write("C01 = %.5g\n"%sim[0,1])
file.write("Chi-squared = %.2f\n"%chisquared)
fullfile.write("Slope = %.3f\n"%slope)
fullfile.write("C10 = %.5g\n"%sim[1,0])
fullfile.write("C01 = %.5g\n"%sim[0,1])
fullfile.write("Chi-squared = %.2f\n"%chisquared)
file.close()
fullfile.close()

