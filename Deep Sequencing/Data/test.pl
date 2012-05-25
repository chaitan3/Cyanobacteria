use Bio::SeqIO;
use Bio::Seq;
$seq_in = Bio::SeqIO->new(-file=>'../../Genome/circular-elongatus.fasta');
$seq = $seq_in->next_seq();
print substr($seq->seq,2,5);
$db_seq = Bio::Seq->new(-seq=>substr($seq->seq,2,5),-id=>$seq->id);
$db = Bio::SeqIO->new(-file=>'>db.fasta');
$db->write_seq($db_seq);
