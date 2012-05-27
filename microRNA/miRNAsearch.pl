use Bio::SeqIO;

$gene_in = Bio::SeqIO->new(-file=>'utr5.fasta');
$gene = $gene_in->next_seq();

$seq_in = Bio::SeqIO->new(-file=>'../Deep Sequencing/Data/Plus7.fasta');
$seq = $seq_in->next_seq();
while(my $seq = $seq_in->next_seq()) {
	$db = Bio::SeqIO->new(-file=>'>tmp.fasta', -flush=>1);
	$db->write_seq($seq);
	$db->write_seq($gene);
	
	open('RNAduplex <tmp.fasta') || die "RNA failed";
}



