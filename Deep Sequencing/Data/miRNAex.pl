#Find overlapping regions of subsets of opposite transcripters with genes. read subset.perl. write final/curated.out
#Find overlapping asRNA regions with psbA1, read full.perl, write stdout

use Bio::SeqFeature::Generic;
use Bio::FeatureIO;
use Bio::SeqIO;
use Bio::Seq;
$seq_in = Bio::SeqIO->new(-file=>'../../Genome/circular-elongatus.fasta');
$seq = $seq_in->next_seq();
$db = Bio::SeqIO->new(-file=>'>db.fasta');

select((select(STDOUT), $|=1)[0]);


my $fin = Bio::FeatureIO->new(-file=>'elongatus.ptt', -format=>'ptt');
my $psb;
my @genes;
while (my $f = $fin->next_feature) {
  push @genes, $f;
}
print "ptt imported\n";

open (rna, 'GSE29264_RNA_Sequencing_Chr_Plus.wig') or die $!;
$s = 0;$e = 0;
foreach $line (<rna>) {
  ($pos, $sig) = split('\t', $line);
  if ($s == $e) {
    if ($sig > 7) {
      $s = $pos;
    }
  }
  elsif ($s > $e) {
    if ($sig < 7) {
	  $e = $pos;
	}
  }
  else {
    $feat = new Bio::SeqFeature::Generic ( -start => $s, -end => $e,
                                -strand => 1);
    $flag = 1;
    foreach (@genes) {
      if ($feat->overlaps($_) and ($feat->strand == $_->strand)) {
        $flag = 0;
      }
    }
    if ($flag) {
   	  print $s;
   	  print ",";
   	  print $e;
   	  print "\n";
   	  $db_seq = Bio::Seq->new(-seq=>substr($seq->seq,$s-1,$e-$s),-id=>$seq->id);
      $db->write_seq($db_seq);
    }
    $s = $e;
  }
} 
close rna;





