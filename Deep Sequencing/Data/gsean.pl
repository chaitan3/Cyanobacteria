#Find overlapping regions of subsets of opposite transcripters with genes. read subset.perl. write final/curated.out
#Find overlapping asRNA regions with psbA1, read full.perl, write stdout

use Text::CSV;
use Bio::SeqFeature::Generic;
use Bio::FeatureIO;


my $fin = Bio::FeatureIO->new(-file=>'elongatus.ptt', -format=>'ptt');
my $psb;
my @genes;
while (my $f = $fin->next_feature) {
  push @genes, $f;
  if (($f->get_tag_values('Synonym'))[0] eq 'Synpcc7942_0424') {
    $psb = $f;
  }
}
print $psb->location->start;
print "\n";

my @transcripts;
my $csv = Text::CSV->new();
open (CSV, "<", 'full.perl') or die $!;
while (<CSV>) {
  $csv->parse($_);
  my @r = $csv->fields();
  if ($r[4] eq '+') {
    $strand=1;
  }
  else {
    $strand=-1;
  }
  $feat = new Bio::SeqFeature::Generic ( -start => $r[1], -end => $r[0],
                                -strand => $strand,
                                -score  => 1000,
                                -tag => {
                                  chrome=>$r[3],
                                  pearson=>$r[5]
                                }
                                     );
  push @transcripts, $feat;
  #~ #antisense RNA
  if ($feat->overlaps($psb) and ($feat->strand == -1) and (($feat->get_tag_values('pearson'))[0] lt 0)) {
    print $feat->location->start;
    print " ";
    print $feat->location->end;
    print " ";
    print $feat->strand;
    print " ";
    @a=$feat->get_tag_values('pearson');
    print $a[0];
    print "\n";
  }
  #~ foreach (@genes) {
    #~ if ($feat->overlaps($_) and ($feat->strand == $_->strand) and (($feat->get_tag_values('chrome'))[0] eq 'CP000100')) {
      #~ print $feat->location->start;
      #~ print " ";
      #~ print $feat->location->end;
      #~ print " ";
      #~ print $feat->strand;
      #~ print " ";
      #~ @a=$feat->get_tag_values('pearson');
      #~ print $a[0];
      #~ print "\n";
      #~ #COG, Length, PID
      #~ $pro = ($_->get_tag_values('Product'))[0];
      #~ $id = ($_->get_tag_values('Synonym'))[0];
      #~ print "$id\n$pro\n\n";
    #~ }
  #~ }
} 
close CSV;




