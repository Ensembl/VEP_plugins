=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

 Will McLaren <wm2@ebi.ac.uk>
    
=cut

=head1 NAME

 MaxEntScan

=head1 SYNOPSIS

 mv MaxEntScan.pm ~/.vep/Plugins
 perl variant_effect_predictor.pl -i variants.vcf --plugin MaxEntScan,[path_to_maxentscan_dir]

=head1 DESCRIPTION

 This is a plugin for the Ensembl Variant Effect Predictor (VEP) that
 runs MaxEntScan (http://genes.mit.edu/burgelab/maxent/Xmaxentscan_scoreseq.html)
 to get splice site predictions.

 The plugin copies most of the code verbatim from the score5.pl and score3.pl
 scripts provided in the MaxEntScan download. To run the plugin you must get and
 unpack the archive from http://genes.mit.edu/burgelab/maxent/download/; the path
 to this unpacked directory is then the param you pass to the --plugin flag.

 The plugin executes the logic from one of the scripts depending on which
 splice region the variant overlaps:

 score5.pl : last 3 bases of exon    --> first 6 bases of intron
 score3.pl : last 20 bases of intron --> first 3 bases of exon

 The plugin reports the reference, alternate and difference (REF - ALT) maximum
 entropy scores.

=cut

package MaxEntScan;

use strict;
use warnings;

use Digest::MD5 qw(md5_hex);

use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use Bio::EnsEMBL::Variation::Utils::VariationEffect qw(overlap);

use Bio::EnsEMBL::Variation::Utils::BaseVepPlugin;
use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

# how many seq/score pairs to cache in memory
our $CACHE_SIZE = 50;

sub new {
  my $class = shift;

  my $self = $class->SUPER::new(@_);
  
  # we need sequence, so no offline mode unless we have FASTA
  die("ERROR: cannot function in offline mode without a FASTA file\n") if $self->{config}->{offline} && !$self->{config}->{fasta};

  my $params = $self->params;

  my $dir = shift @$params;
  die("ERROR: MaxEntScan directory not specified\n") unless $dir;
  die("ERROR: MaxEntScan directory not found\n") unless -d $dir;
  $self->{_dir} = $dir;

  # tell VEP this plugin uses a cache
  $self->{has_cache} = 1;

  ## setup from score5.pl
  $self->{'score5_me2x5'} = $self->score5_makescorematrix($dir.'/me2x5');
  $self->{'score5_seq'}   = $self->score5_makesequencematrix($dir.'/splicemodels/splice5sequences');

  ## setup from score3.pl
  $self->{'score3_metables'} = $self->score3_makemaxentscores;

  return $self;
}

sub feature_types {
  return ['Transcript'];
}

sub get_header_info {
  return {
    MaxEntScan_ref => "MaxEntScan reference sequence score",
    MaxEntScan_alt => "MaxEntScan alternate sequence score",
    MaxEntScan_diff => "MaxEntScan score difference",
  };
}

sub run {
  my ($self, $tva) = @_;

  my $vf = $tva->variation_feature;
  my $tr = $tva->transcript;
  my $tr_strand = $tr->strand;

  return {} unless $vf->{start} == $vf->{end} && $tva->feature_seq =~ /^[ACGT]$/;

  foreach my $intron(@{$tr->get_all_Introns}) {
    # get coords depending on strand
    # MaxEntScan does different predictions for 5 and 3 prime
    # and we need to feed it different bits of sequence for each
    #
    # 5prime, 3 bases of exon, 6 bases of intron:
    # ===------
    #
    # 3prime, 20 bases of intron, 3 bases of exon
    # --------------------===

    my ($five_start, $five_end, $three_start, $three_end);

    if($tr_strand > 0) {
      ($five_start, $five_end)   = ($intron->start - 3, $intron->start + 5);
      ($three_start, $three_end) = ($intron->end - 19, $intron->end + 3);
    }

    else {
      ($five_start, $five_end)   = ($intron->end - 5, $intron->end + 3);
      ($three_start, $three_end) = ($intron->start - 3, $intron->start + 19);
    }

    if(overlap($vf->start, $vf->end, $five_start, $five_end)) {
      my ($ref_seq, $alt_seq) = @{$self->get_seqs($tva, $five_start, $five_end)};

      return {} unless $ref_seq =~ /^[ACGT]+$/ && $alt_seq =~ /^[ACGT]+$/;

      my $ref_score = $self->score5($ref_seq);
      my $alt_score = $self->score5($alt_seq);

      return {
        MaxEntScan_ref => $ref_score,
        MaxEntScan_alt => $alt_score,
        MaxEntScan_diff => $ref_score - $alt_score,
      }
    }

    if(overlap($vf->start, $vf->end, $three_start, $three_end)) {
      my ($ref_seq, $alt_seq) = @{$self->get_seqs($tva, $three_start, $three_end)};

      return {} unless $ref_seq =~ /^[ACGT]+$/ && $alt_seq =~ /^[ACGT]+$/;

      my $ref_score = $self->score3($ref_seq);
      my $alt_score = $self->score3($alt_seq);

      return {
        MaxEntScan_ref => $ref_score,
        MaxEntScan_alt => $alt_score,
        MaxEntScan_diff => $ref_score - $alt_score,
      }
    }
  }

  return {};
}

sub get_seqs {
  my ($self, $tva, $start, $end) = @_;
  my $vf = $tva->variation_feature;

  my $ref_seq = $vf->{slice}->sub_Slice(
    $start,
    $end,
    $tva->transcript->strand
  )->seq;

  my $alt_seq = $ref_seq;
  substr($alt_seq, $vf->{start} - $start, ($vf->{end} - $vf->{start}) + 1) = $tva->feature_seq;

  return [$ref_seq, $alt_seq];
}

sub score5 {
  my $self = shift;
  my $seq = shift;
  my $hex = md5_hex($seq);

  # check cache
  if($self->{cache}) {
    my ($res) = grep {$_->{hex} eq $hex} @{$self->{cache}->{score5}};

    return $res->{score} if $res; 
  }

  my $a = $self->score5_scoreconsensus($seq);
  die("ERROR: No score5_scoreconsensus\n") unless defined($a);

  my $b = $self->score5_getrest($seq);
  die("ERROR: No score5_getrest\n") unless defined($b);

  my $c = $self->{'score5_seq'}->{$b};
  die("ERROR: No score5_seq for $b\n") unless defined($c);

  my $d = $self->{'score5_me2x5'}->{$c};
  die("ERROR: No score5_me2x5 for $c\n") unless defined($d);

  my $score = $self->log2($a * $d);

  # cache it
  push @{$self->{cache}->{score5}}, { hex => $hex, score => $score };
  shift @{$self->{cache}->{score5}} while scalar @{$self->{cache}->{score5}} > $CACHE_SIZE;

  return $score;
}

sub score3 {
  my $self = shift;
  my $seq = shift;
  my $hex = md5_hex($seq);

  # check cache
  if($self->{cache}) {
    my ($res) = grep {$_->{hex} eq $hex} @{$self->{cache}->{score3}};

    return $res->{score} if $res; 
  }

  my $a = $self->score3_scoreconsensus($seq);
  die("ERROR: No score3_scoreconsensus\n") unless defined($a);

  my $b = $self->score3_getrest($seq);
  die("ERROR: No score3_getrest\n") unless defined($b);

  my $c = $self->score3_maxentscore($b, $self->{'score3_metables'});
  die("ERROR: No score3_maxentscore for $b\n") unless defined($c);

  my $score = $self->log2($a * $c);

  # cache it
  push @{$self->{cache}->{score3}}, { hex => $hex, score => $score };
  shift @{$self->{cache}->{score3}} while scalar @{$self->{cache}->{score3}} > $CACHE_SIZE;

  return $score;
}


## methods copied from score5.pl
################################

sub score5_makesequencematrix {
  my $self = shift;
  my $file = shift;
  my %matrix;
  my $n=0;
  open(SCOREF, $file) || die "Can't open $file!\n";
  while(<SCOREF>) { 
    chomp;
    $_=~ s/\s//;
    $matrix{$_} = $n;
    $n++;
  }
  close(SCOREF);
  return \%matrix;
}

sub score5_makescorematrix {
  my $self = shift;
  my $file = shift;
  my %matrix;
  my $n=0;
  open(SCOREF, $file) || die "Can't open $file!\n";
  while(<SCOREF>) { 
    chomp;
    $_=~ s/\s//;
    $matrix{$n} = $_;
    $n++;
  }
  close(SCOREF);
  return \%matrix;
}

sub score5_getrest {
  my $self = shift;
  my $seq = shift;
  my @seqa = split(//,uc($seq));
  return $seqa[0].$seqa[1].$seqa[2].$seqa[5].$seqa[6].$seqa[7].$seqa[8];
}

sub score5_scoreconsensus {
  my $self = shift;
  my $seq = shift;
  my @seqa = split(//,uc($seq));
  my %bgd; 
  $bgd{'A'} = 0.27; 
  $bgd{'C'} = 0.23; 
  $bgd{'G'} = 0.23; 
  $bgd{'T'} = 0.27;  
  my %cons1;
  $cons1{'A'} = 0.004;
  $cons1{'C'} = 0.0032;
  $cons1{'G'} = 0.9896;
  $cons1{'T'} = 0.0032;
  my %cons2;
  $cons2{'A'} = 0.0034; 
  $cons2{'C'} = 0.0039; 
  $cons2{'G'} = 0.0042; 
  $cons2{'T'} = 0.9884;
  my $addscore = $cons1{$seqa[3]}*$cons2{$seqa[4]}/($bgd{$seqa[3]}*$bgd{$seqa[4]}); 
  return $addscore;
}

sub log2 {
  my ($self, $val) = @_;
  return log($val)/log(2);
}


## methods copied from score3.pl
################################

sub score3_hashseq {
  #returns hash of sequence in base 4
  # $self->score3_hashseq('CAGAAGT') returns 4619
  my $self = shift;
  my $seq = shift;
  $seq = uc($seq);
  $seq =~ tr/ACGT/0123/;
  my @seqa = split(//,$seq);
  my $sum = 0;
  my $len = length($seq);
  my @four = (1,4,16,64,256,1024,4096,16384);
  my $i=0;
  while ($i<$len) {
    $sum+= $seqa[$i] * $four[$len - $i -1] ;
    $i++;
  }
  return $sum;
}

sub score3_makemaxentscores {
  my $self = shift;
  my $dir = $self->{'_dir'}."/splicemodels/";
  my @list = ('me2x3acc1','me2x3acc2','me2x3acc3','me2x3acc4',
    'me2x3acc5','me2x3acc6','me2x3acc7','me2x3acc8','me2x3acc9');
  my @metables;
  my $num = 0 ;
  foreach my $file (@list) {
    my $n = 0;
    open (SCOREF,"<".$dir.$file) || die "Can't open $file!\n";
    while(<SCOREF>) {
      chomp;
      $_=~ s/\s//;
      $metables[$num]{$n} = $_;
      $n++;
    }
    close(SCOREF);
    #print STDERR $file."\t".$num."\t".$n."\n";
    $num++;
  }
  return \@metables;
}

sub score3_maxentscore {
  my $self = shift;
  my $seq = shift;
  my $table_ref = shift;
  my @metables = @$table_ref;
  my @sc;
  $sc[0] = $metables[0]{$self->score3_hashseq(substr($seq,0,7))};
  $sc[1] = $metables[1]{$self->score3_hashseq(substr($seq,7,7))};
  $sc[2] = $metables[2]{$self->score3_hashseq(substr($seq,14,7))};
  $sc[3] = $metables[3]{$self->score3_hashseq(substr($seq,4,7))};
  $sc[4] = $metables[4]{$self->score3_hashseq(substr($seq,11,7))};
  $sc[5] = $metables[5]{$self->score3_hashseq(substr($seq,4,3))};
  $sc[6] = $metables[6]{$self->score3_hashseq(substr($seq,7,4))};
  $sc[7] = $metables[7]{$self->score3_hashseq(substr($seq,11,3))};
  $sc[8] = $metables[8]{$self->score3_hashseq(substr($seq,14,4))};
  my $finalscore = $sc[0] * $sc[1] * $sc[2] * $sc[3] * $sc[4] / ($sc[5] * $sc[6] * $sc[7] * $sc[8]);
  return $finalscore;
}

sub score3_getrest {
  my $self = shift;
  my $seq = shift;
  my $seq_noconsensus = substr($seq,0,18).substr($seq,20,3);
  return $seq_noconsensus;
}

sub score3_scoreconsensus {
  my $self = shift;
  my $seq = shift;
  my @seqa = split(//,uc($seq));
  my %bgd; 
  $bgd{'A'} = 0.27; 
  $bgd{'C'} = 0.23; 
  $bgd{'G'} = 0.23; 
  $bgd{'T'} = 0.27;  
  my %cons1;
  $cons1{'A'} = 0.9903;
  $cons1{'C'} = 0.0032;
  $cons1{'G'} = 0.0034;
  $cons1{'T'} = 0.0030;
  my %cons2;
  $cons2{'A'} = 0.0027; 
  $cons2{'C'} = 0.0037; 
  $cons2{'G'} = 0.9905; 
  $cons2{'T'} = 0.0030;
  my $addscore = $cons1{$seqa[18]} * $cons2{$seqa[19]}/ ($bgd{$seqa[18]} * $bgd{$seqa[19]}); 
  return $addscore;
}

1;

