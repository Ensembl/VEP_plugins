=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2017] EMBL-European Bioinformatics Institute

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

 Ensembl <http://www.ensembl.org/info/about/contact/index.html>

=cut

=head1 NAME

 REVEL

=head1 SYNOPSIS

 mv REVEL.pm ~/.vep/Plugins
 ./vep -i variations.vcf --plugin REVEL,/path/to/revel/data.tsv.gz

=head1 DESCRIPTION

 This is a plugin for the Ensembl Variant Effect Predictor (VEP) that
 adds the REVEL score for missense variants to VEP output.

 Please cite the CADD publication alongside the VEP if you use this resource:
 https://www.ncbi.nlm.nih.gov/pubmed/27666373

 The tabix utility must be installed in your path to use this plugin.

=cut

package REVEL;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

sub version {
    return '1.0';
}

sub feature_types {
    return ['Feature','Intergenic'];
}

sub get_header_info {
  return { REVEL => 'Rare Exome Variant Ensemble Learner '};
}
sub new {
  my $class = shift;
  my $self = $class->SUPER::new(@_);

  # test tabix
  die "ERROR: tabix does not seem to be in your path\n" unless `which tabix 2>&1` =~ /tabix$/;

  # get REVEL file
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

  $self->{headers} = `tabix -fh $file 1:1-1`;
  return $self;
}

sub run {
  my ($self, $tva) = @_;

  # only for missense variants
  return {} unless grep {$_->SO_term eq 'missense_variant'} @{$tva->get_all_OverlapConsequences};

  my $vf = $tva->variation_feature;
  my $allele = $tva->variation_feature_seq;
  my $reg_start = $vf->{start};
  my $reg_end = $vf->{end};
  my $pos_string = sprintf("%s:%i-%i", $vf->{chr}, $reg_start, $reg_end);

  # clear cache if it looks like the coords are the same
  # but allele type is different
  delete $self->{cache} if
    defined($self->{cache}->{$pos_string}) &&
    scalar keys %{$self->{cache}->{$pos_string}} &&
    !defined($self->{cache}->{$pos_string}->{$allele});

  my %revel_data;

  # cached?
  if(defined($self->{cache}) && defined($self->{cache}->{$pos_string})) {
    %revel_data = %{$self->{cache}->{$pos_string}};
  }

  # read from file(s)
  else {
    my $file = $self->{file};
    open TABIX, sprintf("tabix -f %s %s |", $file, $pos_string);

    while(<TABIX>) {
      chomp;
      s/\r$//g;
      my ($c, $s, $ref, $alt, $refaa, $altaa, $revel_value) = split /\t/;

      next unless $reg_start == $vf->{start} && $reg_end  == $vf->{end};
      $revel_data{$alt} = {
        REVEL   => $revel_value,
      };
    }

    close TABIX;
  }

  # overwrite cache
  $self->{cache} = {$pos_string => \%revel_data};

  return defined($revel_data{$allele}) ? $revel_data{$allele} : {};
}
1;
