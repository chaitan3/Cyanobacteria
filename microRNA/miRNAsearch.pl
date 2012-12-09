use Bio::SeqIO;

$strand=shift;
if ($strand eq '+') {
	$name = 'Plus';
	$comp = 1;
}
else {
	$name = 'Minus';
	$comp = -1;
}
$infile="$name.fasta";

$gene_in = Bio::SeqIO->new(-file=>'gene.fasta');
$gene = $gene_in->next_seq();
$seq_in = Bio::SeqIO->new(-file=>"../Deep Sequencing/Data/$infile");
#$seq_in = Bio::SeqIO->new(-file=>'miRNA7+5.fasta');

$seq = $seq_in->next_seq();
my @recs; 
while(my $seq = $seq_in->next_seq()) {
	$db = Bio::SeqIO->new(-file=>'>tmp.fasta', -flush=>1,-format=>'raw');
	$db->write_seq($seq);
	$db->write_seq($gene);
	my $r=`RNAduplex <tmp.fasta`;
	$r=~ s/>[^[:space:]]*//g;
	$r=~ s/\n//g;
	$r=~ s/\)/ /g;
	$r=~ s/\(/ /g;
	my @a = split(' ', $r);
	$score = $a[-1];
	my @b = split(',', $a[-2]);
	my @c = split(',', $a[-4]);
	my $match=$b[1]-$b[0];
	$id = $seq->display_id;
	$seq->display_id("$id,$c[0],$c[1],$b[0],$b[1]");
	my @last = ($seq->display_id, $seq->length, $match, $score, $seq, $c[0], $c[1]);
  push @recs, \@last;
}

print "Done\n";
$mi = Bio::SeqIO->new(-file=>'>db.fasta', -flush=>1);
@sorted = sort {$a->[3] <=> $b->[3]} @recs;
$size = @recs;
for ($i=0; $i<$size; $i++)  {
	$newseq = Bio::Seq->new( -display_id => $sorted[$i][0],-seq => $sorted[$i][-3]->subseq($sorted[$i][-2],$sorted[$i][-1]));
	$mi->write_seq($newseq);
	foreach $j (@{$sorted[$i]}) {
		print "$j ";
	}
	print "\n";
}





