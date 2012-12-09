from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein
from Bio import SeqIO

mcut = log2(2)

recs = []
x = arange(0, len(rna))
y = zeros(len(x))
for i in arange(nbp, len(rna), nbp+nba):
  if c.predict(transform([arange(i-nbp,i+nba)])[0]):
    print transform([arange(i-nbp,i+nba)])
    r = arange(i, i+nbp)
    m = mean(rna[r])
    window = r + nbp
    while window[-1] < len(rna):
      if abs(mean(rna[window]) - m) < mcut:
        window += nbp
      else:
        break
    if window[0] - r[0] > nbp:
      print r[0]
      recs.append(SeqRecord(Seq(str(pseq[r[0]:window[0]])), id=str(r[0])+','+ \
      str(window[0]-1)+',1',description=''))
SeqIO.write(recs, "all-transcripts.fasta", "fasta")
      
      
