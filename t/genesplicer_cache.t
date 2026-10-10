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
my $plugin = bless {cache => []}, 'GeneSplicer';
my $calls = 0;
{
  no warnings qw(redefine once);
  local *GeneSplicer::md5_hex = sub { ++$calls; Digest::MD5::md5_hex($_[0]) };
  for my $size (0, 1, 50) {
    $plugin->{cache} = [map { {hex=>Digest::MD5::md5_hex("SEQ$_"),results=>[$_]} } 1..$size];
    $calls = 0;
    is($plugin->results_from_cache('missing'), undef, "cache miss at size $size");
    is($calls, 1, "one digest for a miss at size $size");
    if ($size) {
      $calls = 0;
      is_deeply($plugin->results_from_cache('SEQ1'), [1], 'cached result returned');
      is($calls, 1, 'one digest for a hit');
    }
  }
  $plugin->{cache} = [{hex=>Digest::MD5::md5_hex('empty'),results=>[]}];
  is_deeply($plugin->results_from_cache('empty'), [], 'cached empty results remain a hit');
}
done_testing();
