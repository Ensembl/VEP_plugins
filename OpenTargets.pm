=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2025] EMBL-European Bioinformatics Institute

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

 OpenTargets

=head1 SYNOPSIS

 mv OpenTargets.pm ~/.vep/Plugins

 # print Open Targets Genetics scores and respective gene identifiers (default)
 ./vep -i variations.vcf --plugin OpenTargets,file=path/to/data.tsv.bz
 
 # print all information from Open Targets Genetics
 ./vep -i variations.vcf --plugin OpenTargets,file=path/to/data.tsv.bz,cols=all

=head1 DESCRIPTION

 An Ensembl VEP plugin that integrates data from Open Targets Platform
 (https://platform.opentargets.org), a tool that highlights variant-centric
 statistical evidence to allow both prioritisation of candidate causal variants
 at trait-associated loci and identification of potential drug targets.

 Data from Open Targets Genetics includes two types of associations: GWAS based
 and QTL base. For GWAS studies, the data contains a locus-to-gene (L2G) scores
 to predict causal genes at GWAS loci.

 The tabix utility must be installed in your path to use this plugin. The Open
 Targets Genetics file and respective index (TBI) file can be downloaded from:
 https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/latest/etc/export/vep_data/

 Before using the plugin, split the Open Targets file in two:
   1) a file containing all GWAS associations (and corresponding tabix index file)
    (zcat open_targets_vep.tsv.bgz | head -n 1 | sed 's/.*/#&/'; zcat open_targets_vep.tsv.bgz | awk -F $'\t' '$7 == "gwas"' ) | cut -d $'\t' -f 1-16,17-19 | bgzip -c > open_targets_vep_gwas.tsv.bgz
    tabix -s 1 -b 2 -e 2 open_targets_vep_gwas.tsv.bgz
   2) a file containing all QTL associations (and corresponding tabix index file)
    (zcat open_targets_vep.tsv.bgz | head -n 1 | sed 's/.*/#&/'; zcat open_targets_vep.tsv.bgz | awk -F $'\t' '$7 ~ /qtl$/' ) | cut -d $'\t' -f 1-16,20,21 | bgzip -c > open_targets_vep_qtl.tsv.bgz
    tabix -s 1 -b 2 -e 2 open_targets_vep_qtl.tsv.bgz

  To use the original file (with all associations), comment out the header line and regenerate the tabix-index to include the header.
    zcat open_targets_vep.tsv.bgz | sed '1 s/.*/#&/' | bgzip -c > open_targets_vep.tsv.bgz.tmp; mv open_targets_vep.tsv.bgz.tmp open_targets_vep.tsv.bgz
    tabix -s 1 -b 2 -e 2 open_targets_vep_qtl.tsv.bgz

 Options are passed to the plugin as key=value pairs:
   file : (mandatory) Tabix-indexed file from Open Targets platform. File should contain GWAS- or QTL-annotations, or both.
   lead_variants_only : (optional) boolean to filter out non-lead variant associations, 1 (filtered) or 0 (unfiltered). (default: 1)
   cols : (optional) Colon-separated list of columns to return from the plugin file.
          (default: "gwasLocusToGeneScore:gwasGeneId:gwasDiseases" on GWAS annotations and
          "qtlGeneId:qtlBiosampleName" for QTL annotations); use 'all' to print all data

 Please cite the Open Targets Genetics publication alongside Ensembl VEP if 
 you use this resource: https://doi.org/10.1093/nar/gkaa84

=cut

package OpenTargets;

use strict;
use warnings;

use File::Basename;
use Bio::SeqUtils;
use Bio::EnsEMBL::Variation::Utils::Sequence qw(get_matched_variant_alleles);

use Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin;
use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin);

sub _plugin_name {
  return 'OpenTargets';
}

sub _prefix_cols {
  my $cols = shift;
  my $prefix = _plugin_name() . '_';
  my @prefixed_cols = map { $_ =~ /^$prefix/ ? $_ : $prefix . $_ } @$cols;
  return \@prefixed_cols;
}

sub _get_colnames {
  my $self = shift;

  # Open file header
  open IN, "tabix -H " . $self->{_files}[0] . " |"
    or die "ERROR: cannot open tabix file for " . $self->{_files}[0];

  # Get last line from header
  my $last;
  $last = $_ while <IN>;
  close IN;

  $last =~ s/(^#|\n$)//g
    or die "ERROR: no header found in tabix file for " . $self->{_files}[0];

  # Parse column names from header
  my @cols = split /\t/, $last;
  @cols = splice @cols, 4; # first columns only identify the variant
  return \@cols;
}

sub _get_studytypes {
  my $self = shift;
  my $file = shift;
  my $gwas = my $qtl = 0;

  if( grep(/^gwas/, @{ $self->{colnames} }) ){
    $gwas = 1;
  }

  if( grep(/^qtl/, @{ $self->{colnames} }) ){
    $qtl = 1;
  }

  return ($gwas, $qtl);
}

sub _parse_colnames {
  my $self = shift;
  my $cols = shift;

  # Parse output columns filter
  if ($cols eq "all") {
    $self->{cols} = $self->{colnames};
  } else {
    my @cols = split(/:/, $cols);
    $self->{cols} = \@cols;

    # Check validity of all columns
    my @invalid_cols;
    for my $col (@{ $self->{cols} }) {
      push(@invalid_cols, $col) unless grep(/^$col$/, @{ $self->{colnames} });
    }

    die "\n  ERROR: The following columns were not found in file header: ",
      join(", ", @invalid_cols), "\n" if @invalid_cols;
  }
}

sub new {
  my $class = shift;  
  my $self = $class->SUPER::new(@_);

  $self->expand_left(0);
  $self->expand_right(0);
  $self->get_user_params();

  my $param_hash = $self->params_to_hash();

  # Check file
  my $file = $param_hash->{file};
  die "\n  ERROR: No file specified\nTry using 'file=path/to/data.tsv.bz'\n"
     unless defined($file);
  $self->add_file($file);

  # Get file columns
  $self->{colnames} = $self->_get_colnames();

  # Get study types
  _get_studytypes($self);

  (my $gwas, my $qtl) = _get_studytypes($self);

  $self->{studytypes} = {};
  $self->{studytypes}{GWAS} = $gwas;
  $self->{studytypes}{QTL} = $qtl;

  die "\n  ERROR: Input file does not contain GWAS nor QTL data (missing columns).\n"
     unless $self->{studytypes}->{GWAS} || $self->{studytypes}->{QTL};

  # parse lead-variant param
  if (exists($param_hash->{lead_variants_only})) {
    $self->{lead_variants_only} = $param_hash->{lead_variants_only};
  }
  else{
    $self->{lead_variants_only} = 1;
  }

  # Parse column names filter
  my $cols = $param_hash->{cols};
  if(!$cols){
    $cols = "";
    if($self->{studytypes}->{GWAS}){
      $cols .= ":gwasLocusToGeneScore:gwasGeneId:gwasDiseases";
    }
    if($self->{studytypes}->{QTL}){
      $cols .= ":qtlGeneId:qtlBiosampleName";
    }

    $cols =~ s/^://;
  }
  $self->_parse_colnames($cols);

  # Prefix all column names
  $self->{colnames} = _prefix_cols $self->{colnames};
  $self->{cols}     = _prefix_cols $self->{cols};

  return $self;
}

sub feature_types {
  return ['Feature', 'Intergenic'];
}

sub get_header_info {
  my $self = shift;
  my %header;
  my @keys = @{ $self->{colnames} };

  my $description = "column from " . basename $self->{_files}[0];
  my @vals = map { $description } @keys;
  @header{ @keys } = @vals;

  # Custom headers
  $header{_plugin_name() . '_gwasLocusToGeneScore'}  = "Locus-to-gene (L2G) scores to predict causal genes at GWAS loci; " . $description;

  # Filter by user-selected columns
  %header = map { $_ => $header{$_} } @{ $self->{cols} };

  return \%header;
}

sub _filter_selected_cols {
  my $res = shift;
  my $cols = shift;

  my %tmp = %$res;
  %tmp = map { $_ => $res->{$_} } @$cols;
  return \%tmp;
}

sub _join_results {
  my $self = shift;
  my $all_results = shift;
  my $res = shift;

  # Filter user-selected columns
  $res = _filter_selected_cols($res, $self->{cols});

  if ($self->{config}->{output_format} eq 'json' || $self->{config}->{rest}) {
    # Group results for JSON and REST
    my $name = _plugin_name();
    $all_results->{$name} = [] unless defined $all_results->{$name};
    push(@{ $all_results->{$name} }, $res);
  } else {
    # Create array of results per key
    for (keys %$res) {
      $all_results->{$_} = [] if !$all_results->{$_};
      push(@{ $all_results->{$_} }, $res->{$_} || "NA");
    }
  }
  return $all_results;
}

sub run {
  my ($self, $tva) = @_;
  
  my $vf = $tva->variation_feature;
  my $allele = $tva->base_variation_feature->alt_alleles;
  my @data = @{$self->get_data($vf->{chr}, $vf->{start} - 2, $vf->{end})};

  my @filtered_data;
  if($self->{lead_variants_only}) {
    # Filter @data to only include lead variants
    @filtered_data = grep { $_->{result}->{'OpenTargets_isLead'} eq 'true' } @data;
  }
  else{
    @filtered_data = @data;
  }

  my $all_results = {};
  foreach (@filtered_data) {
    my $matches = get_matched_variant_alleles(
      {
        ref    => $vf->ref_allele_string,
        alts   => $allele,
        pos    => $vf->{start},
        strand => $vf->strand
      },
      {
       ref  => $_->{ref},
       alts => [$_->{alt}],
       pos  => $_->{start},
      }
    );
    $all_results = $self->_join_results($all_results, $_->{result}) if @$matches;
  }
  return $all_results;
}

sub parse_data {
  my ($self, $line) = @_;
  my ($chrom, $pos, $ref, $alt, @vals) = split /\t/, $line;
  my $start = $pos;
  my $end   = $pos;

  # VCF-like adjustment of mismatched substitutions for comparison with VEP
  if(length($alt) != length($ref)) {
    my $first_ref = substr($ref, 0, 1);
    my $first_alt = substr($alt, 0, 1);
    if ($first_ref eq $first_alt) {
      $start++;
      $ref = substr($ref, 1);
      $alt = substr($alt, 1);
      $ref ||= '-';
      $alt ||= '-';
    }
  }
  
  my %res;
  @res{ @{ $self->{colnames} } } = @vals;

  return {
    ref => $ref,
    alt => $alt,
    start => $start,
    result => \%res
  };
}

sub get_start {
  return $_[1]->{start};
}

sub get_end {
  return $_[1]->{end};
}

1;
