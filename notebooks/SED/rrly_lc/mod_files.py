import sys
list = open(sys.argv[1])
fnum = 1
for file in list:
  fh = open(file.rstrip())
  ofh = open("%i_normed.txt"%(fnum), "w")
  for l in fh:
    if l.startswith("#"):
      els = l.rstrip().split()
      period = float(els[2])
      norms = []
      for n in els[3:]:
        norms.append(float(n))
      norms.append(float(n))
      ofh.write(l.rstrip()+"\n")
    else:
      els = l.rstrip().split()
      line = []
      line.append(float(els[0])*period)
      for i in range(5):
        line.append(float(els[1+i])*norms[i])
      line.append(float(els[5])*norms[5])
      ofh.write(" ".join(els)+"\n")
  ofh.close()
  fh.close()
  fnum += 1
