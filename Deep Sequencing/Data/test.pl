use Bio::SeqIO;
use Bio::Seq;
$seq_in = Bio::SeqIO->new(-file=>'../../Genome/circular-elongatus.fasta');
$seq = $seq_in->next_seq();
use Bio::Location::Simple;
my $location = Bio::Location::Simple->new(-start  => 46,
                                          -end   => 63,
                                          -strand => "-1");


print $seq->subseq(46,63);
print "\n";
print $seq->subseq($location);
print "\n";
