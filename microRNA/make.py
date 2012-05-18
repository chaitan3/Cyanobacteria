from Bio import SeqIO

utr_size = 180
window_size = 30

gene = 413873
utr = gene-utr_size
search = gene-600

r = SeqIO.parse("../Genome/circular-elongatus.fasta",'fasta')
s = r.next()
mirnas = []
for i in range(0, 151, 10):
	mirnas.append(s[utr+i:utr+i+window_size])
SeqIO.write(mirnas,'utr.fasta','fasta')
SeqIO.write([s[search:gene]],'search.fasta','fasta')


