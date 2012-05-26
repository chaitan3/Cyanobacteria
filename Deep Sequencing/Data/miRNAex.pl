#Find overlapping regions of subsets of opposite transcripters with genes. read subset.perl. write final/curated.out
#Find overlapping asRNA regions with psbA1, read full.perl, write stdout

$cut=int(shift);
$strand=shift;
if ($strand=='+') {
	$name = 'Plus';
	$comp = 1;
}
else {
	$name = 'Minus';
	$comp = 0;
}
$file="GSE29264_RNA_Sequencing_Chr_$name.wig";
$outfile=">$name$cut.fasta";
print "$outfile\n";
use Bio::SeqFeature::Generic;
use Bio::FeatureIO;
use Bio::SeqIO;
use Bio::Seq;
use Bio::Location::Simple;
$seq_in = Bio::SeqIO->new(-file=>'../../Genome/circular-elongatus.fasta');
$seq = $seq_in->next_seq();
$db = Bio::SeqIO->new(-file=>$outfile);

select((select(STDOUT), $|=1)[0]);


my $fin = Bio::FeatureIO->new(-file=>'elongatus.ptt', -format=>'ptt');
my $psb;
my @genes;
while (my $f = $fin->next_feature) {
  push @genes, $f;
}

open (rna, $file) or die $!;
$s = 0;$e = 0;
foreach $line (<rna>) {
  ($pos, $sig) = split('\t', $line);
  if ($s == $e) {
    if ($sig > $cut) {
      $s = $pos;
    }
  }
  elsif ($s > $e) {
    if ($sig < $cut) {
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
	  if ($comp) {
		$location = Bio::Location::Simple->new(-start  => $s,
                                          -end   => $e-1,
                                          -strand => "-1");
	    $pseq = $seq->subseq($location);
      }
	  else {
	    $pseq = $seq->subseq($s,$e-1);
   	  }
	  $db_seq = Bio::Seq->new(-seq=>$pseq,-id=>$seq->id);
      $db->write_seq($db_seq);
    }
    $s = $e;
  }
} 
close rna;





