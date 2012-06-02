
import subprocess
from Bio import SeqIO

neg = 229169
gene = 34489
time = [60, 68, 72, 80]
base = '/home/chaitan3/GSE-Raw/GSM7242'
num = 28

c = 0
for seq in SeqIO.parse("../../microRNA/miRNA7.fasta", "fasta"):
	c += 1
	print seq.id
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
		probe = start/12 + neg
		pref = 'R'
	print probe, c
	
	for i in range(0, len(time)):
		filename = base + str(num + i) + '_' + pref + str(time[i])
		row = pref + str(probe)
		print filename
		
		print subprocess.check_output('grep ' + row + ' ' + filename)
	
	
