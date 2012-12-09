use Bio::SeqIO;

$gene_in = Bio::SeqIO->new(-file=>'gene.fasta');
$gene = $gene_in->next_seq();
$seq_in = Bio::SeqIO->new(-file=>'db.fasta');

$c=0;
mkdir 'rna';
chdir 'rna';
while(my $seq = $seq_in->next_seq()) {
	
	mkdir $c;
	chdir $c;
	open (MYFILE, '>tmp.fasta');
	$seq1=$seq->seq;
	$seq2=$gene->seq;
	print MYFILE "$seq1\&$seq2\n";
	close (MYFILE); 
	my $r=`RNAcofold -p <tmp.fasta`;
	open (MY, '>rna.out');
	print MY $r;
	chdir '..';
	$c++;
}
chdir '..';



