import csv
#import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from Bio import SeqIO

neg = 229169
gene = 34489
time = [60, 68, 72, 80]

f = open("GSE29264-GPL13535_series_matrix.txt")
while f.readline()[0]=='!':
  pass
val = list()
rows = csv.reader(f, delimiter='\t')
val.append([])
for r in rows:
  if r[0][0] == '!':
    break
  val.append([0]*(len(r)-1))
  for i in range(1,len(r)):
    val[-1][i-1] = float(r[i])
print 'Microarray imported'

c = 0
p = []
db = open("../../microRNA/db.fasta")
for seq in SeqIO.parse(db, "fasta"):
	c += 1
	print seq.id
	s = seq.id.split(',')
	s[0] = s[0][1:]
	strand = int(s[2])
	if strand == 1:
		start = int(s[0]) + int(s[3])
		end = int(s[0]) + int(s[4])
		probe = start/12
	else:
		start = int(s[1]) - int(s[4])
		end = int(s[1]) - int(s[3])
		probe = start/12 + neg
	print probe
	p.append(pearsonr(val[probe], val[gene])[0])

	print c, pearsonr(val[probe], val[gene])
#	plt.plot(time, val[probe], label=str(c))

print sorted(enumerate(p),key=lambda x:x[1],reverse=True)

#plt.plot(time, val[gene],'r--',label='psbA')
#plt.legend()
#plt.show()
