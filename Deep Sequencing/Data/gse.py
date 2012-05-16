#Generate pearson coefficients from microarray data with respect to psbA1
#write to total.mat
import csv
from scipy.stats import pearsonr
from scipy.io import savemat
from numpy import array, zeros

f = open("GSE29264-GPL13535_series_matrix.txt")
while f.readline()[0]=='!':
  pass
val = list()
rows = csv.reader(f, delimiter='\t')
genes = 0
for r in rows:
  if r[0][0] == '!':
    break
  val.append(zeros(len(r)))
  genes += 1
  for i in range(1,len(r)):
    val[-1][i-1] = r[i]
p = []
test = 34490
for i in range(0, genes):
  p.append((pearsonr(val[i], val[test-1])[0],i+1))
p.sort()
savemat('total',{'p':p})
#s=p[0:100]
#s.sort()
#print s

