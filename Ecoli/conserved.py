from Bio import SeqIO

win = 20
spread = range(0, win)
pseq = SeqIO.parse(open("genome.fasta"), "fasta").next().seq

f = open('genome.ptt')
for i in range(0, 3):
  f.readline()
lines = f.readlines()
f.close()

p = []
for i in spread:
  p.append({'A':0,'T':0,'G':0,'C':0})
c = 0
genes = []
for line in lines:
  vals = line.split('\t')
  val = vals[0].split('..')
  start = int(val[0])-1
  end = int(val[1])-1
  if start < win:
    continue
  if vals[1] == '+':
    #1 is for strand, check if it is an mRNA
    genes.append([start, end, 1])
    window = pseq[start - win:start]
    c += 1
    for i in spread:
      p[i][window[i]] += 1

for i in spread:
  for j in p[i].keys():
    p[i][j] = log(1 + p[i][j]*1.0/c)
    
#~ maxs = []
#~ for i in spread:
  #~ print sorted(p[i].items(), key=lambda codon: codon[1], reverse=True), c, sum(p[i].values())
  #~ maxs.append(p[i]['a'] + p[i]['t'])
#~ plot(spread, maxs)
#~ show()
hmm = p[8:14]

