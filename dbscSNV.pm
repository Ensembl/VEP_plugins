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

 dbscSNV

=head1 SYNOPSIS

  mv dbscSNV.pm ~/.vep/Plugins
  perl variant_effect_predictor.pl -i variations.vcf --plugin dbscSNV,/path/to/dbscSNV.txt.gz

=head1 DESCRIPTION

  A VEP plugin that retrieves data for splicing variants from a tabix-indexed
  dbscSNV file.

  Please cite the dbscSNV publication alongside the VEP if you use this resource:
  http://nar.oxfordjournals.org/content/42/22/13534

  The tabix utility must be installed in your path to use this plugin. The dbscSNV
  data file can be downloaded from
  https://sites.google.com/site/jpopgen/dbNSFP.

  The file must be processed and indexed by tabix before use by this plugin:

  > wget ftp://dbscsnv:dbscsnv@dbscsnv.softgenetics.com/dbscSNV.zip
  > unzip dbscSNV.zip
  > head -n1 dbscSNV.chr1 > h
  > cat dbscSNV.chr* | grep -v ^chr | cat h - | bgzip -c > dbscSNV.txt.gz
  > tabix -s 1 -b 2 -e 2 -c c dbscSNV.txt.gz

  Note that in the last command we tell tabix that the header line starts with "c";
  this may change to the default of "#" in future versions of dbscSNV.

  Tabix also allows the data file to be hosted on a remote server. This plugin is
  fully compatible with such a setup - simply use the URL of the remote file:

  --plugin dbscSNV,http://my.files.com/dbscSNV.txt.gz

  Note that transcript sequences referred to in dbscSNV may be out of sync with
  those in the latest release of Ensembl; this may lead to discrepancies with
  scores retrieved from other sources.
 
=cut

package dbscSNV;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);

use Bio::EnsEMBL::Variation::Utils::BaseVepPlugin;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

sub new {
  my $class = shift;
  
  my $self = $class->SUPER::new(@_);
  
  # test tabix
  die "ERROR: tabix does not seem to be in your path\n" unless `which tabix 2>&1` =~ /tabix$/;
  
  # get dbscSNV file
  my $file = $self->params->[0];
  
  # remote files?
  if($file =~ /tp\:\/\//) {
    my $remote_test = `tabix -f $file 1:1-1 2>&1`;
    if($remote_test && $remote_test !~ /get_local_version/) {
      die "$remote_test\nERROR: Could not find file or index file for remote annotation file $file\n";
    }
  }

  # check files exist
  else {
    die "ERROR: dbscSNV file $file not found\n" unless -e $file;
    die "ERROR: Tabix index file $file\.tbi not found - perhaps you need to create it first?\n" unless -e $file.'.tbi';
  }
  
  $self->{file} = $file;
  
  # get headers
  open HEAD, "tabix -fh $file 1:1-1 2>&1 | ";
  while(<HEAD>) {
    next unless /^c/;
    chomp;
    $self->{headers} = [split];
  }
  close HEAD;
  
  die "ERROR: Could not read headers from $file\n" unless defined($self->{headers}) && scalar @{$self->{headers}};
  
  # check alt and Ensembl_transcriptid headers
  foreach my $h(qw(alt Ensembl_gene)) {
    die "ERROR: Could not find required column $h in $file\n" unless grep {$_ eq $h} @{$self->{headers}};
  }
  
  $self->{cols} = {
    'ada_score' => 1,
    'rf_score'  => 1
  };
  
  return $self;
}

sub feature_types {
  return ['Transcript'];
}

sub variation_feature_types {
  return ['VariationFeature'];
}

sub get_header_info {
  return {
    ada_score => 'dbscSNV ADA score',
    rf_score  => 'dbscSNV RF score'
  }
}

sub run {
  my ($self, $tva) = @_;
  
  my $vf = $tva->variation_feature;
  
  return {} unless $vf->{start} eq $vf->{end};
  return {} unless grep {$_->SO_term =~ /splic/} @{$tva->get_all_OverlapConsequences};
  
  # get allele, reverse comp if needed
  my $allele = $tva->variation_feature_seq;
  reverse_comp(\$allele) if $vf->{strand} < 0;
  
  return {} unless $allele =~ /^[ACGT]$/;
  
  # get gene stable ID
  my $g_id = $tva->transcript->{_gene_stable_id} || $tva->transcript->gene->stable_id;
    
  my $pos_string = sprintf("%s:%i-%i", $vf->{chr}, $vf->{start}, $vf->{end});
  
  my @dbscsnv_data;
  
  # cached?
  if(defined($self->{cache}) && defined($self->{cache}->{$pos_string})) {
    @dbscsnv_data = @{$self->{cache}->{$pos_string}};
  }
  
  # read from file
  else {
    open TABIX, sprintf("tabix -f %s %s |", $self->{file}, $pos_string);
    
    while(<TABIX>) {
      chomp;
      s/\r$//g;
      my @split = split /\t/;
      
      # parse data into hash of col names and values
      my %data = map {$self->{headers}->[$_] => $split[$_]} (0..(scalar @{$self->{headers}} - 1));
      
      push @dbscsnv_data, \%data;
    }
    
    close TABIX;
  }
  
  # overwrite cache
  $self->{cache} = {$pos_string => \@dbscsnv_data};
  
  my $data;
  
  foreach my $tmp_data(@dbscsnv_data) {
    # compare allele and transcript
    next unless
      defined($tmp_data->{alt}) &&
      $tmp_data->{alt} eq $allele; # &&
#       defined($tmp_data->{Ensembl_gene}) &&
#       $tmp_data->{Ensembl_gene} =~ /$g_id($|;)/;
    
    $data = $tmp_data;
    last;
  }
  
  return {} unless scalar keys %$data;
  
  $DB::single = 1;
  
  # get required data
  my %return =
    map {$_ => $data->{$_}}
    grep {$data->{$_} ne '.'}              # ignore missing data
    grep {defined($self->{cols}->{$_})}  # only include selected cols
    keys %$data;
  
  return \%return;
}

1;

