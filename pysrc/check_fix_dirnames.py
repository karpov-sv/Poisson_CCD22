from pylab import *
import os, sys, time, h5py
sys.path.append(os.path.realpath('./pysrc'))
from pysubs import *  # These are the plotting subroutines
from subprocess import *

# Craig Lage - 14Sep17
# This program checks that the outputfiledir names are correct and fixes them if not.

maindir = "baredata"
datadirs = os.listdir(maindir)

linelist = range(48,82)

# First check if the outputfiledir matches the directory it's in
for datadir in datadirs:
    files = os.listdir(maindir+'/'+datadir)

    for file in files:
        if list(file)[-3:] == ['c','f','g']:
            break
    incfgfile = maindir+'/'+datadir+'/'+file
    ConfigData = ReadConfigFile(incfgfile)
    dirname = ConfigData['outputfiledir'].split('/')[-1]
    if datadir == dirname:
        print "File %s names match"%dirname
    else:
        try:
            print "File %s names don't match. Fixing."%incfgfile
            lines = OpenFile(incfgfile)
            for i,line in enumerate(lines):
                if len(line.split()) == 0:
                    continue
                if line.split()[0] == 'outputfiledir':
                    #print "Found it in line %d"%i
                    lines[i] = 'outputfiledir = '+maindir+'/'+datadir+'\n'
            newcfgfile = open(incfgfile, 'w')
            for line in lines:
                newcfgfile.write(line)
            newcfgfile.close()
        except:
            print "File %s names don't match, but fix failed."%dirname
            continue

# Now check that the silicon parameters match between directories
files = os.listdir(maindir+'/pixel1')
for file in files:
    if list(file)[-3:] == ['c','f','g']:
        firstcfg = file
        break
print "Comparison file = %s"%file
firstcfgfile = maindir+'/pixel1/'+file
firstlines = OpenFile(firstcfgfile)

for datadir in datadirs:
    if datadir == 'smallcapg' or datadir == 'smallcapf':
        continue
    files = os.listdir(maindir+'/'+datadir)
    for file in files:
        if list(file)[-3:] == ['c','f','g']:
            break
    print "Checking file %s"%file
    incfgfile = maindir+'/'+datadir+'/'+file
    lines = OpenFile(incfgfile)
    for i in linelist:
        if len(lines[i]) < 2:
            continue
        if lines[i].split()[1] != '=':
            continue

        param = lines[i].split('=')[1].split('#')[0].strip()
        firstparam = firstlines[i].split('=')[1].split('#')[0].strip()
        if param != firstparam:
            print "Param %s in file %s doesn't match param %s in %s"%(param, file, firstparam,firstcfg)
