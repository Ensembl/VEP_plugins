=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2023] EMBL-European Bioinformatics Institute

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

 CADD

=head1 SYNOPSIS

 mv CADD.pm ~/.vep/Plugins
 ./vep -i variations.vcf --plugin CADD,/FULL_PATH_TO_CADD_FILE/whole_genome_SNVs.tsv.gz,/FULL_PATH_TO_CADD_FILE/InDels.tsv.gz

=head1 DESCRIPTION

 A VEP plugin that retrieves CADD scores for variants from one or more
 tabix-indexed CADD data files.
 
 Please cite the CADD publication alongside the VEP if you use this resource:
 https://www.ncbi.nlm.nih.gov/pubmed/24487276
 
 The tabix utility must be installed in your path to use this plugin. The CADD
 data files can be downloaded from
 http://cadd.gs.washington.edu/download

 The plugin works with all versions of available CADD files. The plugin only
 reports scores and does not consider any additional annotations from a CADD
 file. It is therefore sufficient to use CADD files without the additional
 annotations. 
 
=cut

package CADD;

use strict;
use warnings;

use Bio::EnsEMBL::Variation::Utils::Sequence qw(get_matched_variant_alleles);

use Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin);

# List here all columns in headers that should be included
my %INCLUDE_COLUMNS = (
    "PHRED" => {
      "name" => "CADD_PHRED",
      "description" => 'PHRED-like scaled CADD score'
    },
    "RawScore" => {
      "name" => "CADD_RAW",
      "description" => 'Raw CADD score'
    }
);

sub new {
  my $class = shift;
  
  my $self = $class->SUPER::new(@_);

  # Test if tabix exists
  die "\nERROR: tabix does not seem to be in your path\n" unless `which tabix 2>&1` =~ /tabix$/;

  $self->expand_left(0);
  $self->expand_right(0);

  $self->get_user_params();

  # Check files in arguments
  my @params = @{$self->params};

  die "\nERROR: No CADD files specified\nTip: Add a file after command, example:\nvep ... --plugin CADD,/FULL_PATH_TO_CADD_FILE/whole_genome_SNVs.tsv.gz\n" unless @params > 0;
  $self->add_file($_) for @params;

  my $assembly = $self->{config}->{assembly};

  $self->{header} = ();

  foreach my $file (@params) {
    open IN, "tabix -f -h ".$file." 1:1-1 |";

    my @lines = <IN>;

    my @assembly_header_matches =  grep { /$assembly/ } @lines if (defined($assembly));

    if (!@assembly_header_matches && $assembly) {
      die "\nERROR: Assembly is " . $assembly .
          " but CADD file does not contain " .
          $assembly . " in header.\n";
    }

    while (my $line = shift @lines) {

      next if (rindex $line, "#Chrom", 0);
      chomp $line;
      $self->{$file} = $line;
    }

    # Make sure it has a known prefix in header
    die "'#Chrom' was not found on header" unless $self->{$file};

    my $file_check = 0;
    # Conditional header
    for (split /\t/, $self->{$file}){
      next unless (exists($INCLUDE_COLUMNS{$_}));
      $file_check = 1;
      $self->{header}{$INCLUDE_COLUMNS{$_}{"name"}} = $INCLUDE_COLUMNS{$_}{"description"};
    }

    die "\nERROR: $file does not have a known column to be included" unless $file_check;

  }

  close IN;

  return $self;
}

sub feature_types {
  return ['Feature','Intergenic'];
}

sub get_header_info {
  my $self = shift;

  return $self->{header}
}

sub run {
  my ($self, $tva) = @_;
  
  my $vf = $tva->variation_feature;
  
  # get allele
  my $allele = $tva->variation_feature_seq;
  
  return {} unless $allele =~ /^[ACGT-]+$/;

  my @data =  @{$self->get_data($vf->{chr}, $vf->{start} - 2, $vf->{end})};

  foreach (@data) {
    my $matches = get_matched_variant_alleles(
      {
        ref    => $vf->ref_allele_string,
        alts   => [$allele],
        pos    => $vf->{start},
        strand => $vf->strand
      },
      {
       ref  => $_->{ref},
       alts => [$_->{alt}],
       pos  => $_->{start},
      }
    );
    return $_->{result} if (@$matches);
  }
  return {};
}

sub parse_data {
  my ($self, $line, $file) = @_;

  my @headers = split /\t/, $self->{$file};
  my @values = split /\t/, $line;

  my %data = map {$headers[$_] => $values[$_]} (0..(@headers - 1));

  my $c = $data{"#Chrom"};
  my $s = $data{"Pos"};
  my $ref = $data{"Ref"};
  my $alt = $data{"Alt"};

  # Conditional result
  my %result = ();
  foreach (keys %INCLUDE_COLUMNS){
    next unless (exists($data{$_}));
    $result{$INCLUDE_COLUMNS{$_}{"name"}} = $data{$_};
  }

  # do VCF-like coord adjustment for mismatched subs
  my $end = ($s + length($ref)) - 1;
  if(length($alt) != length($ref)) {
    my $first_ref = substr($ref, 0, 1);
    my $first_alt = substr($alt, 0, 1);
    if ($first_ref eq $first_alt) {
      $s++;
      $ref = substr($ref, 1);
      $alt = substr($alt, 1);
      $ref ||= '-';
      $alt ||= '-';
    }
  }

  return {
    ref => $ref,
    alt => $alt,
    start => $s,
    end => $end,
    result => \%result
  };
}

sub get_start {
  return $_[1]->{start};
}

sub get_end {
  return $_[1]->{end};
}

1;
