from pylab import *
import os, sys, time, h5py
sys.path.append(os.path.realpath('./pysrc'))
from pysubs import *  # These are the plotting subroutines
from subprocess import *

# Craig Lage - 14Sep17
# This program copies over the output to a smaller directory for transferring without the
# large data files.

maindir = "data"
transdir = "transdata"
baredir = "baredata"

# Create the new directories if they don't exist
if not os.path.isdir(transdir):
    os.mkdir(transdir)
if not os.path.isdir(baredir):
    os.mkdir(baredir)

datadirs = os.listdir(maindir)
for datadir in datadirs:
    subdir_0 = False
    #print datadir
    if '_' in list(datadir): 
        if list(datadir)[-2:] != ['_','0']:
            continue
        else:
            subdir_0 = True
    os.mkdir(transdir+'/'+datadir)
    files = os.listdir(maindir+'/'+datadir)
    for file in files:
        if (list(file)[-3:] == ['c','f','g']) or (list(file)[-3:] == ['t','x','t']):
            cmd = 'cp %s/%s/%s  %s/%s/'%(maindir,datadir,file,transdir,datadir)
            cp3 = Popen(cmd, shell=True)
            Popen.wait(cp3)
        else:
            continue
        try:
            cmd = 'cp -r %s/%s/plots  %s/%s/'%(maindir,datadir,transdir,datadir)
            cp4 = Popen(cmd, shell=True)
            Popen.wait(cp4)
        except:
            continue
    if subdir_0:
        continue
    if not os.path.isdir(baredir+'/'+datadir):
        os.mkdir(baredir+'/'+datadir)
    files = os.listdir(maindir+'/'+datadir)
    for file in files:
        if list(file)[-3:] == ['c','f','g']:
            cmd = 'cp %s/%s/%s  %s/%s/'%(maindir,datadir,file,baredir,datadir)
            cp3 = Popen(cmd, shell=True)
            Popen.wait(cp3)
        else:
            continue
