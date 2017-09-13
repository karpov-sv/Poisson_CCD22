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


q = 1.6E-19   # MKS
me = 9.11E-31   # MKS
k = 1.38E-23   # MKS
T = 273.0
E=1000.0
Tox = 1.0E-5
Eps_ox = 8.85E-14 * 4.2
mu = mu_si(E, T)
print "Mobility = %f"%mu
Vds = 0.5
L = 5.0E-4
Vgates = []
Is = []
for run in range(13):
    ConfigData = ReadConfigFile('data/run%d/trans.cfg'%run)

    Vgate = ConfigData["FixedRegionVoltage_11"]
    file = open('data/run%d/charge.txt'%run,'r')
    lines = file.readlines()
    file.close
    Ne = float(lines[0].split()[4].strip(','))
    print run, Vgate, Ne
    I = Vds * q * mu * Ne / L**2
    Vgates.append(Vgate)
    Is.append(I)

Qss = ConfigData["ChannelSurfaceCharge"]
[Vgs, Ids] = Read_STA3800_IV_Data("STA3800_meas.xls")

#Delta_V = -0.7
#Delta_Q = Eps_ox / Tox * Delta_V / q

figure()
ax1=axes([0.2,0.1,0.6,0.6])
ax1.set_title("STA3800C Id-Vg")
ax1.plot(Vgs, array(Ids)*1000.0, color = 'green', lw = 2, label='Measured')
ax1.scatter(array(Vgates), array(Is)*1000.0, marker = 'x', color='red', label='Sim - Qss = %g'%Qss)
#ax1.scatter(array(Vgates)+Delta_V, array(Is)*1000.0, marker = 'x', color='red', label='Sim - QS=1.4E12')
ax1.set_xlabel("Vgs (volts)")
ax1.set_ylabel("Ids(mA)")
ax1.set_ylim(0,1.0)
ax1.legend(loc='upper left')
ax1.text(-23, 0.3, 'Vds = 0.5V')
#ax1.text(-23, 0.5, 'DeltaV = %.1f, DeltaQ = %g'%(Delta_V, Delta_Q))
savefig("IdVg_12Sep17.pdf")
