#!/usr/bin/env python

#Author: Craig Lage, NYU; 
#Date: 14-Sep-17

#This program manages the running of multiple Poisson tests
from pylab import *
import os, sys, time, subprocess
from mpi4py import MPI
sys.path.append(os.path.realpath('./pysrc'))
from pysubs import *
#****************MAIN PROGRAM*****************
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
run = rank

NumTrans = 13 # Number of transistor runs
NumSats = 24 # Number of saturation runs
NumSpots = size - 3 # Number of BF spots
NumRuns = 80 # Number of spot fluxes

if run == 0:
    # This processor runs the pixel runs
    cmd = '~/Software/Poisson_CCD_Hole20/src/Poisson data/smallpixel/smallpixel.cfg'
    print "Launching %s"%cmd
    cpp_job = subprocess.Popen(cmd, shell=True)
    subprocess.Popen.wait(cpp_job)
    time.sleep(1.0)
    poiplot = subprocess.Popen('python pysrc/Poisson_Small.py data/smallpixel/smallpixel.cfg 0', shell=True) 
    subprocess.Popen.wait(poiplot)
    poiplot2 = subprocess.Popen('python pysrc/Poisson_Convergence.py data/smallpixel/smallpixel.cfg 0', shell=True) 
    subprocess.Popen.wait(poiplot2)

    cmd = '~/Software/Poisson_CCD_Hole20/src/Poisson data/smallcapg/smallcap.cfg'
    print "Launching %s"%cmd
    cpp_job = subprocess.Popen(cmd, shell=True)
    subprocess.Popen.wait(cpp_job)
    time.sleep(1.0)
    poiplot = subprocess.Popen('python pysrc/Poisson_CapConvergence.py data/smallcapg/smallcap.cfg 0', shell=True) 
    subprocess.Popen.wait(poiplot)

    cmd = '~/Software/Poisson_CCD_Hole20/src/Poisson data/smallcapf/smallcap.cfg'
    print "Launching %s"%cmd
    cpp_job = subprocess.Popen(cmd, shell=True)
    subprocess.Popen.wait(cpp_job)
    time.sleep(1.0)
    poiplot = subprocess.Popen('python pysrc/Poisson_CapConvergence.py data/smallcapf/smallcap.cfg 0', shell=True) 
    subprocess.Popen.wait(poiplot)

    for pixrun in range(4):
        cmd = '~/Software/Poisson_CCD_Hole20/src/Poisson data/pixel%d/pixel.cfg'%pixrun
        print "Launching %s"%cmd
        cpp_job = subprocess.Popen(cmd, shell=True)
        subprocess.Popen.wait(cpp_job)
        time.sleep(1.0)
        if pixrun < 3:
            poiplot = subprocess.Popen('python pysrc/Poisson_Plots.py data/pixel%d/pixel.cfg 0'%pixrun, shell=True) 
            subprocess.Popen.wait(poiplot)
            chargeplot = subprocess.Popen('python pysrc/ChargePlots.py data/pixel%d/pixel.cfg 0 2'%pixrun, shell=True) 
            subprocess.Popen.wait(chargeplot)
        if pixrun > 0 and pixrun < 3:
            pixplot = subprocess.Popen('python pysrc/VertexPlot.py data/pixel%d/pixel.cfg 0 2'%pixrun, shell=True) 
            subprocess.Popen.wait(pixplot)
            areaplot = subprocess.Popen('python pysrc/Area_Covariance_Plot.py data/pixel%d/pixel.cfg 0'%pixrun, shell=True) 
            subprocess.Popen.wait(areaplot)
        if pixrun == 3:
            poiplot = subprocess.Popen('python pysrc/Poisson_Plots.py data/pixel%d/pixel.cfg 80'%pixrun, shell=True) 
            subprocess.Popen.wait(poiplot)
            chargeplot = subprocess.Popen('python pysrc/ChargePlots.py data/pixel%d/pixel.cfg 80 2'%pixrun, shell=True) 
            subprocess.Popen.wait(chargeplot)
            pixplot = subprocess.Popen('python pysrc/VertexPlot.py data/pixel%d/pixel.cfg 80 2'%pixrun, shell=True) 
            subprocess.Popen.wait(pixplot)
            areaplot = subprocess.Popen('python pysrc/Area_Covariance_Plot.py data/pixel%d/pixel.cfg 80'%pixrun, shell=True) 
            subprocess.Popen.wait(areaplot)
        if pixrun == 1:
            simsplot = subprocess.Popen('python pysrc/SIMS_Comparison.py data/pixel%d/pixel.cfg 0'%pixrun, shell=True) 
            subprocess.Popen.wait(simsplot)

elif run == 1:
    # This processor runs edge, IO, and sat runs

    cmd = '~/Software/Poisson_CCD_Hole20/src/Poisson data/edgerun/edge.cfg'
    print "Launching %s"%cmd
    cpp_job = subprocess.Popen(cmd, shell=True)
    subprocess.Popen.wait(cpp_job)
    time.sleep(1.0)
    poiplot = subprocess.Popen('python pysrc/Poisson_Edge.py data/edgerun/edge.cfg 0', shell=True) 
    subprocess.Popen.wait(poiplot)
    chargeplot = subprocess.Popen('python pysrc/Pixel_Shift.py data/edgerun/edge.cfg 0', shell=True) 
    subprocess.Popen.wait(chargeplot)

    cmd = '~/Software/Poisson_CCD_Hole20/src/Poisson data/iorun/io.cfg'
    print "Launching %s"%cmd
    cpp_job = subprocess.Popen(cmd, shell=True)
    subprocess.Popen.wait(cpp_job)
    time.sleep(1.0)
    poiplot = subprocess.Popen('python pysrc/Poisson_IO.py data/iorun/io.cfg 0', shell=True) 
    subprocess.Popen.wait(poiplot)
    chargeplot = subprocess.Popen('python pysrc/ChargePlots_IO.py data/iorun/io.cfg 0', shell=True) 
    subprocess.Popen.wait(chargeplot)

    # This processor runs the saturation runs, sequentially
    Vph = linspace(0.0,5.0,6)
    Vpl = linspace(-2.0,-8.0,4)

    for satrun in range(NumSats):
        newdir = subprocess.Popen('mkdir -p data/satrun_%d'%satrun, shell=True) 
        subprocess.Popen.wait(newdir)
        incfgfile = 'data/satrun/sat.cfg'
        outcfgfile = 'data/satrun_%d/sat.cfg'%satrun
        ConfigData = ReadConfigFile(incfgfile)
        Vphrun = Vph[satrun%6]
        Vplrun = Vpl[(satrun - satrun%6) / 6]
        New_Sat_Cfg_File(incfgfile, outcfgfile, satrun, Vphrun, Vplrun)
        time.sleep(1.0)
        
        cmd = '~/Software/Poisson_CCD_Hole20/src/Poisson data/satrun_%d/sat.cfg'%satrun
        print "Launching %s"%cmd
        cpp_job = subprocess.Popen(cmd, shell=True)
        subprocess.Popen.wait(cpp_job)
        time.sleep(1.0)
        poiplot = subprocess.Popen('python pysrc/Poisson_Sat.py data/satrun_%d/sat.cfg 0'%satrun, shell=True) 
        subprocess.Popen.wait(poiplot)
        chargeplot = subprocess.Popen('python pysrc/ChargePlots.py data/satrun_%d/sat.cfg 0 9'%satrun, shell=True) 
        subprocess.Popen.wait(chargeplot)
        barplot = subprocess.Popen('python pysrc/Barrier.py data/satrun_%d/sat.cfg 0'%satrun, shell=True) 
        subprocess.Popen.wait(barplot)
    
    time.sleep(1.0)
    satplot = subprocess.Popen('python pysrc/Plot_SatLevel.py data/satrun/sat.cfg 0 %d'%NumSats, shell=True) 
    subprocess.Popen.wait(satplot)

elif run == 2:
    # This processor runs the transistor runs, sequentially, and the treering run
    Vgates = linspace(-20.0,10.0,NumTrans)

    for transrun in range(NumTrans):
        newdir = subprocess.Popen('mkdir -p data/transrun_%d'%transrun, shell=True) 
        subprocess.Popen.wait(newdir)
        incfgfile = 'data/transrun/trans.cfg'
        outcfgfile = 'data/transrun_%d/trans.cfg'%transrun
        ConfigData = ReadConfigFile(incfgfile)
        New_Trans_Cfg_File(incfgfile, outcfgfile, transrun, Vgates[transrun])
        time.sleep(1.0)
        
        cmd = '~/Software/Poisson_CCD_Hole20/src/Poisson data/transrun_%d/trans.cfg'%transrun
        print "Launching %s"%cmd
        cpp_job = subprocess.Popen(cmd, shell=True)
        subprocess.Popen.wait(cpp_job)
        time.sleep(1.0)
        poiplot = subprocess.Popen('python pysrc/Poisson_Trans.py data/transrun_%d/trans.cfg 0'%transrun, shell=True) 
        subprocess.Popen.wait(poiplot)
        chargeplot = subprocess.Popen('python pysrc/ChargePlots_Trans.py data/transrun_%d/trans.cfg 0'%transrun, shell=True) 
        subprocess.Popen.wait(chargeplot)
    
    time.sleep(1.0)
    mosplot = subprocess.Popen('python pysrc/MOSFET_Calculated_IV.py data/transrun/trans.cfg 0', shell=True) 
    subprocess.Popen.wait(mosplot)
    time.sleep(1.0)
    cmd = '~/Software/Poisson_CCD_Hole20/src_25sep17/Poisson data/treering/treering.cfg'
    print "Launching %s"%cmd
    cpp_job = subprocess.Popen(cmd, shell=True)
    subprocess.Popen.wait(cpp_job)
    time.sleep(1.0)
    pixplot = subprocess.Popen('python pysrc/VertexPlot.py data/treering/treering.cfg 0 4', shell=True) 
    subprocess.Popen.wait(pixplot)

else:
    # All remaining processors process BF spots
    bfrun = run - 3
    newdir = subprocess.Popen('mkdir -p data/bfrun_%d'%bfrun, shell=True) 
    subprocess.Popen.wait(newdir)
    incfgfile = 'data/bfrun/bf.cfg'
    outcfgfile = 'data/bfrun_%d/bf.cfg'%bfrun
    ConfigData = ReadConfigFile(incfgfile)
    New_BF_Cfg_File(incfgfile, outcfgfile, bfrun)
    time.sleep(1.0)
    
    cmd = '~/Software/Poisson_CCD_Hole20/src/Poisson data/bfrun_%d/bf.cfg'%bfrun
    print "Launching %s"%cmd
    cpp_job = subprocess.Popen(cmd, shell=True)
    subprocess.Popen.wait(cpp_job)
    time.sleep(1.0)
    poiplot = subprocess.Popen('python pysrc/Poisson_Plots.py data/bfrun_%d/bf.cfg 80'%bfrun, shell=True) 
    subprocess.Popen.wait(poiplot)
    chargeplot = subprocess.Popen('python pysrc/ChargePlots.py data/bfrun_%d/bf.cfg 80 9'%bfrun, shell=True) 
    subprocess.Popen.wait(chargeplot)


    if run == size - 1:
        # Use this one to plot the spot results
        # First wait until everyone is done
        FinishedSpots = 0
        while FinishedSpots < NumSpots:
            FinishedSpots = int(subprocess.Popen('ls data/bfrun_*/*%d_CC.dat | wc -l'%NumRuns,shell=True, stdout = subprocess.PIPE).communicate()[0].strip('\n'))
            time.sleep(30.0)

        bfplot = subprocess.Popen('python pysrc/Plot_BF_Spots.py data/bfrun/bf.cfg 0 %d %d'%(NumSpots,NumRuns), shell=True) 
        subprocess.Popen.wait(bfplot)

