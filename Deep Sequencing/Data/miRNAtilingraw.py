from math import log
import commands
from Bio import SeqIO

neg = 229169
gene = 34489
time = [60, 68, 72, 80]
base = '/home/chaitan3/GSE-Raw/GSM7242'
num = 28
seqio= open("../../microRNA/miRNA7.fasta")
acc = []

c = 0
for seq in SeqIO.parse(seqio, "fasta"):
	c += 1
	s = seq.id.split(',')
	s[0] = s[0][1:]
	strand = int(s[2])
	if strand == 1:
		start = int(s[0]) + int(s[3])
		end = int(s[0]) + int(s[4])
		probe = start/12
		pref = 'F'
	else:
		start = int(s[1]) - int(s[4])
		end = int(s[1]) - int(s[3])
		probe = start/12# + neg
		pref = 'R'
	print probe, c
	fg = []
	bg = []
#	probe = 34489
	for i in range(0, len(time)):
		filename = base + str(num + i) + '_' + pref + str(time[i]) + '.gpr'
		row = pref + str(probe)
		out = commands.getoutput('grep ' + row + ' ' + filename)
		fields = out.split('\t')
		#4-Name, 9-F635 Mean, 14-B635 Mean, 44-Log Ratio, 47-F-B635,48-F-B532
		fg.append(float(fields[48]))
		bg.append(float(fields[47]))
	if min(fg) > 100:
		print "qualifies"
		acc.append(c)
		for i in range(0, len(time)):
			print fg[i], bg[i], log(fg[i]/bg[i],2)
print acc
		
	
	
	
