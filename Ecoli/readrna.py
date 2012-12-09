from pylab import *

f = open('GSM613993_Early_Exponential.basecount.txt')
lines = f.readlines()
f.close()
rna = []
for line in lines:
  if len(line) > 0:
    rna.append(int(line.split('\t')[1]))
# log    
rna = log2(array(rna)+1)

#talk about
# negative data generation
# training and test sets
# sets implication, comparative and actual results
# RNA sequencing multiple time points data situation
