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
NxCenter = 4
NyCenter = 4
Area_0 = 99.9972

areafilename = outputfiledir + '/' + outputfilebase +'_%d_Area.dat'%run
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
savefig(outputfiledir+"/plots/Area_Covariance_%d.pdf"%run)


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

