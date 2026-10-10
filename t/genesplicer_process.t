use strict;
use warnings;
use Test::More;
use File::Temp qw(tempdir);
use FindBin;
BEGIN {
  package Bio::EnsEMBL::Utils::Sequence;
  sub import {}
  $INC{'Bio/EnsEMBL/Utils/Sequence.pm'} = __FILE__;
  package Bio::EnsEMBL::Variation::Utils::VariationEffect;
  sub import {}
  $INC{'Bio/EnsEMBL/Variation/Utils/VariationEffect.pm'} = __FILE__;
  package Bio::EnsEMBL::Variation::Utils::BaseVepPlugin;
  sub new { my ($class,$config,@params)=@_; bless {config=>$config,params=>\@params},$class }
  sub params { $_[0]->{params} }
  $INC{'Bio/EnsEMBL/Variation/Utils/BaseVepPlugin.pm'} = __FILE__;
}
require "$FindBin::Bin/../GeneSplicer.pm";
my $dir = tempdir(CLEANUP=>1);
my $bin = "$dir/predictor";
mkdir "$dir/training data";
mkdir "$dir/temporary files";
open my $fh, '>', $bin or die $!;
print {$fh} <<'PREDICTOR';
#!/usr/bin/perl
use strict;
use warnings;
if (!@ARGV) { print "USAGE predictor input training -f output\n"; exit 1 }
my ($input,$training,$flag,$output) = @ARGV;
die 'argument paths were split' unless -f $input && -d $training && $flag eq '-f';
open my $in,'<',$input or die $!;
<$in>; my $seq=<$in>; close $in; chomp $seq;
if ($seq eq 'FAIL') { print STDERR "controlled prediction failure\n"; exit 7 }
if ($seq eq 'MISSING') { exit 0 }
open my $out,'>',$output or die $!;
print {$out} "10 11 2.5 High donor\n" if $seq ne 'EMPTY';
close $out;
if ($seq eq 'PARTIAL') { print STDERR "partial prediction failure\n"; exit 9 }
PREDICTOR
close $fh;
chmod 0755, $bin;
my $plugin = GeneSplicer->new({}, "binary=$bin", "training=$dir/training data", "tmpdir=$dir/temporary files");
for my $seq ('FAIL','PARTIAL') {
  my $result=eval { $plugin->results_from_seq($seq) };
  like($@, qr/genesplicer failed.*prediction failure/s, "$seq raises with predictor diagnostic");
  ok(!defined($result), "$seq does not return predictions");
  my @left=glob("$dir/temporary\\ files/*");
  is(scalar(@left),0,"$seq temporary files removed");
}
is_deeply($plugin->results_from_seq('OK'), [{end5=>10,end3=>11,score=>2.5,confidence=>'High',type=>'donor'}], 'valid prediction with spaces in paths');
is_deeply($plugin->results_from_seq('EMPTY'), [], 'successful empty file');
is_deeply($plugin->results_from_seq('MISSING'), [], 'successful no-file prediction preserves behavior');
$plugin->{_bin} = "$dir/nonexistent";
eval { $plugin->results_from_seq('OK') };
like($@,qr/(exec|open3)/, 'execution failure raises');
done_testing();
