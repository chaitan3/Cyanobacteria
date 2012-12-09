#Find overlapping regions of subsets of opposite transcripters with genes. read subset.perl. write final/curated.out
#Find overlapping asRNA regions with psbA1, read full.perl, write stdout

$strand=shift;
if ($strand eq '+') {
	$name = 'Plus';
	$comp = 1;
}
else {
	$name = 'Minus';
	$comp = -1;
}
$outfile=">$name.fasta";
print "$outfile\n";
use Bio::SeqFeature::Generic;
use Bio::SeqFeature::Collection;
use Bio::FeatureIO;
use Bio::SeqIO;
use Bio::Seq;
$seq_in = Bio::SeqIO->new(-file=>'../../mlmodel/all-transcripts.fasta');
$db = Bio::SeqIO->new(-file=>$outfile);

select((select(STDOUT), $|=1)[0]);


my $fin = Bio::FeatureIO->new(-file=>'../../Genome/circular-elongatus.ptt', -format=>'ptt');
my $psb;
my @genes;
while (my $f = $fin->next_feature) {
  push @genes, $f;
}
my $col = Bio::SeqFeature::Collection->new();
$col->add_features(\@genes);
my $count = 0;
while(my $seq = $seq_in->next_seq()) {
  $flag = 1;
  my @s = split(',', $seq->display_id);
  $feat = new Bio::SeqFeature::Generic ( -start => $s[0], -end => $s[1],
                                -strand => $s[2]);
  my @feats = $col->features_in_range(-start => $s[0], -end => $s[1], -strand => $s[2]
  -contain => false, -strandmatch => 'strong');
  my $size = $#feats+1;
  print "$s[0]\n";
  if ($size > 0) {
    $count = $count + 1;
    $db->write_seq($seq);
  }
} 
print $count;
close rna;





