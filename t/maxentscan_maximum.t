use strict;
use warnings;
use Test::More;
use FindBin;
# The scorer is controlled so this tests selection independently of model data.
BEGIN {
  package Bio::EnsEMBL::Utils::Sequence;
  sub import {}
  $INC{'Bio/EnsEMBL/Utils/Sequence.pm'} = __FILE__;
  package Bio::EnsEMBL::Variation::Utils::VariationEffect;
  sub import {}
  $INC{'Bio/EnsEMBL/Variation/Utils/VariationEffect.pm'} = __FILE__;
  package Bio::EnsEMBL::Variation::Utils::BaseVepPlugin;
  sub new { bless {}, shift }
  $INC{'Bio/EnsEMBL/Variation/Utils/BaseVepPlugin.pm'} = __FILE__;
}
require "$FindBin::Bin/../MaxEntScan.pm";
{
  package TestMaxEntScan;
  our @ISA = ('MaxEntScan');
  sub score5 { shift @{$_[0]->{scores}} }
  sub score3 { shift @{$_[0]->{scores}} }
}
for my $method ('get_max_donor', 'get_max_acceptor') {
  my $width = $method eq 'get_max_donor' ? 9 : 23;
  for my $case ([0,-1,0,1], [-1,0,0,2], [0,0,0,1], [-2,-1,-1,2], [1,0,1,1]) {
    my ($first, $second, $want, $frame) = @$case;
    my $plugin = bless {scores => [$first, $second]}, 'TestMaxEntScan';
    my $got = $plugin->$method(('A' x $width).'C');
    is($got->[2], $want, "$method [$first,$second] score");
    is($got->[1], $frame, "$method [$first,$second] frame");
    is($got->[0], $frame == 1 ? 'A' x $width : ('A' x ($width-1)).'C', 'selected sequence');
  }
  my $plugin = bless {scores => []}, 'TestMaxEntScan';
  is_deeply($plugin->$method('A'), [undef,undef,undef], "$method empty window");
}
done_testing();
