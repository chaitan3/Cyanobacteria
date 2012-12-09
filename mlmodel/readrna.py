from pylab import *

f = open('../Deep Sequencing/Data/GSE29264_RNA_Sequencing_Chr_Plus.wig')
f.readline()
f.readline()
lines = f.readlines()
f.close()
rna = []
for line in lines:
  if len(line) > 0:
    rna.append(int(line.split('\t')[1]))
# log    
rna = log2(array(rna)+1)
