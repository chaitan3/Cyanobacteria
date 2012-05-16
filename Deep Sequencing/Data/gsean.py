#Create a subset of opposite transcripters with psbA1, use data by gse.py(total.mat)
#Make a map from id to start stop positions, write to subset/full.perl
from scipy.io import loadmat
import csv
p=loadmat('total.mat')['p']
l=p[0:,1].copy()
l.sort()
#print l

p[:,[1,0]]=p[:,[0,1]]
ps = dict(p)

f = open("GPL13535-11452.txt")
while f.readline()[0]=='#':
  pass
val = list()
rows = csv.reader(f, delimiter='\t')
genes = 0
for r in rows:
  val.append(dict())
  val[genes]['chr']=r[1]
  val[genes]['db']=r[2]
  val[genes]['strand']=r[3]
  val[genes]['start']=r[4]
  val[genes]['stop']=r[5]
  genes += 1
print genes
f.close()

#stop start chr db strand
f = open('full.perl','w')
w = csv.writer(f)
for i in l:
  v = val[int(i-1)]
  v['pearson'] = ps[i]
  #print val[int(i-1)]
  w.writerow([v['stop'],v['start'], v['chr'], v['db'], v['strand'], v['pearson']])
f.close()
