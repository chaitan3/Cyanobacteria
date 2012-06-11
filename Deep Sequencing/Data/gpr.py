import numpy as np
import matplotlib.pylab as plt
import csv
f = open("/home/chaitan3/GSE-Raw/GSM724228_F60.gpr")
while f.readline()[0]!='"':
  pass
while f.readline()[0]=='"':
  pass
imeanA = list()
imean = list()
isd = list()
rows = csv.reader(f, delimiter='\t')
genes = 0
for r in rows:
  imeanA.append([int(r[9])])
  imean.append(int(r[9]))
  isd.append(int(r[10]))
imeanA = np.array(imeanA)
imean = np.array(imean)
isd = np.array(isd)
m = np.linalg.lstsq(imeanA, isd)[0][0]
print m
isdd = isd - m*imean
print np.mean(isdd), np.std(isdd)
plt.plot(imean, isd,'.')
x = 100000
plt.plot([0,x],[0,m*x])
plt.show()
