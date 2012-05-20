from Bio import SeqIO

utr5_size = 47
utr3_size = 63


start = 413873
end = 414955




r = SeqIO.parse("../Genome/circular-elongatus.fasta",'fasta')
s = r.next()

reg = [
#3'UTR
s[415031:415062],
s[415078:415122],
s[415151:415177],
s[415217:415253],
s[415291:415310],
s[415352:415377]

]
SeqIO.write(reg,'db.fasta','fasta')

SeqIO.write([s[end+utr3_size:end+500]],'search.fasta','fasta')


