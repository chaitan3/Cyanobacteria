from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
#~ gene = 'kaiA'
#~ gid='637799647'
#~ seq = SeqIO.parse(open("../Genome/circular-elongatus.fasta"), "fasta").next()
#~ seq = seq.seq[1241411:1241517].reverse_complement()

#~ gene = 'kaiB'
#~ gid='637799646'
#~ seq = SeqIO.parse(open("../Genome/circular-elongatus.fasta"), "fasta").next()
#~ seq = seq.seq[1240467:1240535].reverse_complement()

#~ gene = 'psbA1'
#~ gid='637798830'
#~ seq = SeqIO.parse(open("../Genome/circular-elongatus.fasta"), "fasta").next()
#~ seq = seq.seq[413831:413873]#.reverse_complement()

#~ gene = 'nir'
#~ gid='637799669'
#~ seq = SeqIO.parse(open("../Genome/circular-elongatus.fasta"), "fasta").next()
#~ seq = seq.seq[1264801:1264857].reverse_complement()

gene = 'pex'
gid='637799089'
seq = SeqIO.parse(open("../Genome/circular-elongatus.fasta"), "fasta").next()
seq = seq.seq[672462:672657].reverse_complement()



recs = SeqRecord(seq, id='Elongatus 7942, '+gid,description=gene+'UTR5 region for binding')
SeqIO.write(recs, gene+"-utr5.fasta", "fasta")
