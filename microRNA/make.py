from Bio import SeqIO





r = SeqIO.parse("../Genome/circular-elongatus.fasta",'fasta')
s = r.next()

#psbA

#~ utr5_size = 47
#~ utr3_size = 63
#~ 
#~ 
#~ start = 413873
#~ end = 414955

#~ reg = [
#~ s[413937:413843:-1]
#~ s[415031:415062],
#~ s[415078:415122],
#~ s[415151:415177],
#~ s[415217:415253],
#~ s[415291:415310],
#~ s[415352:415377]
#~ ]

#kaiA
utr5_size = 66
utr3_size = 11


start = 1241411
end = 1240557

reg = [
#Anti sense
s[1240836:1240876],
s[1240925:1240957],
s[1240995:1241023],
#3'UTR
s[1240527:1240501:-1]
]
SeqIO.write(reg,'kaiA-db.fasta','fasta')

SeqIO.write([s[start+utr5_size:start:-1]],'kaiA-utr5.fasta','fasta')
SeqIO.write([s[end:end-utr3_size:-1]],'kaiA-utr3.fasta','fasta')
