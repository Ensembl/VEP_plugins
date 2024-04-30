=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2024] EMBL-European Bioinformatics Institute

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

 Geno2MP

=head1 SYNOPSIS

 cp Geno2MP.pm ${HOME}/.vep/Plugins
 ./vep -i variations.vcf --plugin Geno2MP,file=/path/to/Geno2MP/data.vcf.gz

 # Return more columns from Geno2MP VCF file
 ./vep -i variations.vcf --plugin Geno2MP,file=/path/to/Geno2MP/data.vcf.gz,cols=HPO_CT:FXN:nhomalt_male_aff:nhomalt_male_unaff

 # Build and return Geno2MP URL based on GRCh37 variant location
 ./vep -i variations.vcf --plugin Geno2MP,file=/path/to/Geno2MP/data.vcf.gz,url=1

=head1 DESCRIPTION

 A VEP plugin that adds information from Geno2MP, a web-accessible database of
 rare variant genotypes linked to phenotypic information.

 Parameters can be set using a key=value system:
   file : VCF file containing Geno2MP data
   cols : colon-delimited list of Geno2MP columns to return from INFO fields
          (by default it only returns the column HPO_CT)
   url  : build and return URL to Geno2MP variant page (boolean; 0 by default);
          the variant location in Geno2MP website is based on GRCh37 coordinates

 Please cite Geno2MP alongside the VEP if you use this resource:
 Geno2MP, NHGRI/NHLBI University of Washington-Center for Mendelian Genomics (UW-CMG), Seattle, WA
 (URL: http://geno2mp.gs.washington.edu [date (month, yr) accessed]).

=cut
package Geno2MP;

use strict;
use warnings;

use Bio::EnsEMBL::Variation::Utils::Sequence qw(get_matched_variant_alleles);
use Bio::EnsEMBL::Variation::Utils::VEP qw(parse_line);

use Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin;
use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin);

sub prepare_header_info {
  my $self = shift;
  
  # Header information based on Geno2MP file (for user-selected columns only)
  $self->{header_info} = {};
  open IN, "tabix -f -h " . $self->{_files}[0] . " 1:1-1 |";
  my @header = <IN>;
  while (my $line = shift @header) {
    next unless $line =~ /##INFO/;
    
    # Get column name and description for user-selected columns
    my ($col, $description) = $line =~ /.*ID=(.*?),.*,Description="(.*)".*/g;
    next unless $self->{cols}->{$col};
    
    if ($col =~ "HPO_CT") {
      $col = $self->{hpo_ct};
      $description = "Number of phenotypic profiles in Geno2MP";
    }
    my $key = $self->{label} . "_" . $col;
    $self->{header_info}->{$key} = $description;
  }
  close IN;
  
  # Geno2MP URL header information
  if ($self->{url}) {
    my $key = $self->{label} . "_URL";
    $self->{header_info}->{$key} = "Link to Geno2MP variant page"
  }
}

sub new {
  my $class = shift;

  my $self = $class->SUPER::new(@_);
  my $config = $self->{config};

  $self->expand_left(0);
  $self->expand_right(0);

  my $param_hash = $self->params_to_hash();

  # Add Geno2MP file
  my $file = $param_hash->{file};
  die "\nERROR: No Geno2MP file specified\nTry using 'file=path/to/Geno2MP/file.vcf.gz'\n"
    unless defined($file);
  $self->add_file($file);
  
  # Build URL linking to Geno2MP variant page
  $self->{url} = defined($param_hash->{url}) ? $param_hash->{url} : 0;
  
  # Retrieve specific columns from VCF
  my @cols = $param_hash->{cols} ? split(/\:/, $param_hash->{cols}) : "HPO_CT";
  $self->{cols}->{$_} = 1 for @cols;

  # Prepare header information from Geno2MP file header
  $self->{label} = "Geno2MP";
  $self->{hpo_ct} = "HPO_count";
  $self->prepare_header_info();
  
  return $self;
}

sub feature_types {
  return ['Feature', 'Intergenic'];
}

sub get_header_info {
  my $self = shift;
  return $self->{header_info};
}

sub run {
  my ($self, $tva) = @_;
  my $vf = $tva->variation_feature;

  # get allele
  my $alt_alleles = $tva->base_variation_feature->alt_alleles;
  my $ref_allele  = $vf->ref_allele_string;

  my @data = @{ $self->get_data($vf->{chr}, $vf->{start} - 2, $vf->{end}) };
  return {} unless(@data);

  foreach my $var (@data) {
    my $matches = get_matched_variant_alleles(
      {
        ref    => $ref_allele,
        alts   => $alt_alleles,
        pos    => $vf->{start},
        strand => $vf->strand
      },
      {
       ref  => $var->{ref},
       alts => [$var->{alt}],
       pos  => $var->{start},
      }
    );
    return $var->{result} if (@$matches);
  }

  return {};
}

sub parse_data {
  my ($self, $line) = @_;
  
  # Parse VCF fields
  my ($chrom, $pos, $id, $ref, $alt, $qual, $filter, $info) = split /\t/, $line;
  
  # Parse INFO fields
  my (%res, $grch37);
  my @attribs = split/\;/, $info;
  foreach my $att (@attribs){
    my($col, $val) = split/\=/, $att;
    # Save GRCh37 variant location if available
    $grch37 = $val if $col eq "GRCh37";

    if ($self->{cols}->{$col}) {
      $col = $self->{hpo_ct} if $col =~ "HPO_CT";
      my $key = $self->{label} . "_" . $col;
      $res{$key} = $val;
    }
  }
  
  # Add URL based on GRCh37 variant location (Geno2P webpage only supports GRCh37)
  if ($self->{url}) {
    # Guess between SNP or indels (other types not supported in Geno2MP)
    my $type = (length($alt) == 1 && length($ref) == 1) ? "snp" : "indel";
    
    # Get GRCh37 variant location from INFO if available
    my ($c, $p) = $grch37 ? split(/_/, $grch37) : ($chrom, $pos);
    my $key  = $self->{label} . "_URL";
    my $base = "https://geno2mp.gs.washington.edu/Geno2MP";
    $res{$key} = "$base/#/variant/$c/$p/$ref%3E$alt/$type";
  }

  return {
    ref => $ref,
    alt => $alt,
    start => $pos,
    result => \%res
  }
}

sub get_start {
  return $_[1]->{start};
}

sub get_end {
  return $_[1]->{end};
}

1;
