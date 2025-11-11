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

 EVE

=head1 SYNOPSIS

 cp EVE.pm ${HOME}/.vep/Plugins
 ./vep -i variations.vcf --plugin EVE,file=/path/to/eve/data.vcf.gz # By default, Class75 is used.
 ./vep -i variations.vcf --plugin EVE,file=/path/to/eve/data.vcf.gz,class_number=60

=head1 DESCRIPTION

 This is a plugin for the Ensembl Variant Effect Predictor (VEP) that
 adds information from EVE (evolutionary model of variant effect).

 This plugin only report EVE scores for input variants
 and does not merge input lines to report on adjacent variants.
 It is only available for GRCh38.

 Please cite EVE publication alongside Ensembl VEP if you use this resource:
 https://www.nature.com/articles/s41586-021-04043-8

###################################################
# Bash script to merge all VCFs from EVE dataset. #
###################################################

### BEGIN

# EVE input file can be downloaded from https://evemodel.org/api/proteins/bulk/download/ 
# Input: VCF files by protein (vcf_files_missense_mutations inside zip folder)
# Output: Compressed Merged VCF file (vcf.gz) + index file (.tbi)

DATA_FOLDER='/<PATH-TO>/vcf_files_missense_mutations' # Fill this line
OUTPUT_FOLDER='/<PATH-TO>/eve_plugin' # Fill this line
OUTPUT_NAME='eve_merged.vcf' # Default output name

# Get header from first VCF
cat `ls ${DATA_FOLDER}/*vcf | head -n1` > header

# Get variants from all VCFs and add to a single-file
ls ${DATA_FOLDER}/*vcf | while read VCF; do grep -v '^#' ${VCF} >> variants; done

# Merge Header + Variants in a single file
cat header variants | \
awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1V -k2,2n"}' > ${OUTPUT_FOLDER}/${OUTPUT_NAME};

# Remove temporary files
rm header variants

# Compress and index
bgzip ${OUTPUT_FOLDER}/${OUTPUT_NAME};

# If not installed, use: sudo apt install tabix
tabix ${OUTPUT_FOLDER}/${OUTPUT_NAME}.gz;

### END

=cut
package EVE;

use strict;
use warnings;

use Bio::EnsEMBL::Variation::Utils::Sequence qw(get_matched_variant_alleles);

use Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin;
use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin);
use Bio::DB::HTS::Tabix;

sub new {
  my $class = shift;
  my $self = $class->SUPER::new(@_);
  my $config = $self->{config};

  $self->expand_left(0);
  $self->expand_right(0);

  my $assembly = $config->{assembly} || $config->{human_assembly};
  die "\nAssembly is not GRCh38, EVE only works with GRCh38.\n" if ($assembly ne "GRCh38");

  my $param = $self->params_to_hash();

  my $eve_file    = $param->{file} || $param->{eve_file};
  my $popeve_file = $param->{popeve_file};

  die "\nERROR: No input specified\nUse 'file=<eve.vcf.gz>' and/or 'popeve_file=<popeve.vcf.gz>'\n"
    unless ($eve_file || $popeve_file);

  if ($eve_file) {
    $self->add_file($eve_file);
    $self->{has_eve} = 1;
    $self->{eve_file} = $eve_file;
  }
  if ($popeve_file) {
    $self->add_file($popeve_file);
    $self->{has_pop} = 1;
    $self->{pop_file}  = $popeve_file;
  }

  my @valid_class_numbers = (10, 20, 25, 30, 40, 50, 60, 70, 75, 80, 90);
  if (defined($param->{class_number})) {
    my $class_number = $param->{class_number};
    die "\nERROR: This class_number: '$class_number' does not exist.\nTry any of: " . join(', ', @valid_class_numbers)
      unless grep(/^$class_number$/, @valid_class_numbers);
    $self->{class_number} = $class_number;
  } else {
    $self->{class_number} = 75;
  }

  return $self;
}

sub feature_types {
  return ['Transcript'];
}

sub _vcf_info_descriptions {
  my ($file) = @_;
  my %desc;

  return \%desc unless $file;

  eval {
    my $tbx    = Bio::DB::HTS::Tabix->new( filename => $file );
    my $header = $tbx->header;

    foreach my $line (split /\n/, $header) {
      next unless $line =~ /^##INFO<(.+)>/ || $line =~ /^##INFO=<(.+)>/;
      my $body = $1;
      my ($id)   = $body =~ /\bID=([^,>]+)/;
      my ($d)    = $body =~ /Description="([^"]*)"/;

      $desc{$id} = $d if defined $id && defined $d;
    }
  };
  # if tabix/header parsing fails, return an empty hash and let get_header_info use defaults
  return \%desc;
}

sub get_header_info {
  my $self = shift;

  # try to read INFO descriptions from the EVE and popEVE VCF headers
  my %eve_desc;
  my %pop_desc;

  if ($self->{has_eve} && $self->{eve_file}) {
    my $h = _vcf_info_descriptions($self->{eve_file});
    %eve_desc = %{$h} if $h;
  }

  if ($self->{has_pop} && $self->{pop_file}) {
    my $h = _vcf_info_descriptions($self->{pop_file});
    %pop_desc = %{$h} if $h;
  }

  my $class_key = "Class" . $self->{class_number};

  # prefer INFO descriptions pulled from VCF but fall back to hard-coded descriptions
  my $h = {
    EVE_SCORE => ($eve_desc{EVE} || "Score from EVE model"),
    EVE_CLASS => ($eve_desc{$class_key} || "Classification (Benign, Uncertain, or Pathogenic) when setting $self->{class_number}% as uncertain"),
  };

  if ($self->{has_pop}) {
    $h->{popEVE_SCORE}              = ($pop_desc{popEVE}              || "Score from popEVE");
    $h->{popEVE_EVE}                = ($pop_desc{EVE}                 || "Raw EVE score (popEVE file)");
    $h->{popEVE_ESM1v}              = ($pop_desc{ESM1v}               || "Raw ESM1v score (popEVE file)");
    $h->{popEVE_pop_adjusted_EVE}   = ($pop_desc{'pop-adjusted_EVE'}  || $pop_desc{'pop_adjusted_EVE'} || "Population-adjusted EVE");
    $h->{popEVE_pop_adjusted_ESM1v} = ($pop_desc{'pop-adjusted_ESM1v'}|| $pop_desc{'pop_adjusted_ESM1v'} || "Population-adjusted ESM1v");
    $h->{popEVE_gap_frequency}      = ($pop_desc{gap_frequency}       || "Gap frequency");
    $h->{popEVE_gene}               = ($pop_desc{gene}                || "Gene symbol");
    $h->{popEVE_protein}            = ($pop_desc{protein}             || "Protein accession");
    $h->{popEVE_mutant}             = ($pop_desc{mutant}              || "Protein-level change");
  }

  return $h;
}

sub run {
  my ($self, $tva) = @_;
  my $vf = $tva->variation_feature;

  return {} unless grep {$_->SO_term eq 'missense_variant'} @{$tva->get_all_OverlapConsequences};

  my $alt_alleles = $tva->base_variation_feature->alt_alleles;
  my $ref_allele  = $vf->ref_allele_string;

  my @data = @{
    $self->get_data(
      $vf->{chr},
      $vf->{start} - 2,
      $vf->{end}
    )
  };

  return {} unless @data;

  my %out;

  foreach my $variant (@data) {
    my $matches = get_matched_variant_alleles(
      { ref => $ref_allele, alts => $alt_alleles, pos => $vf->{start}, strand => $vf->strand },
      { ref => $variant->{ref},  alts => [ $variant->{alt} ],          pos => $variant->{start} }
    );
    next unless @$matches;

    # MERGE instead of returning immediately
    # merge results from every matching record within the window, instead of returning on the first match
    # this allows EVE (allele in codon format: XXX) and popEVE (allele in SNV format: X) to both annotate the same input variant
    @out{ keys %{ $variant->{result} } } = values %{ $variant->{result} };
  }

  return %out ? \%out : {};
}

sub parse_data {
  my ($self, $line) = @_;
  chomp $line;  # ensure INFO regexes see clean line endings

  my ($chrom, $pos, $id, $ref, $alt, $qual, $filter, $info) = split /\t/, $line, 8;

  # source detection
  my $is_pop = ($info =~ /(;\s*)?(popEVE|protein|gene|mutant|gap_frequency|ESM1v|pop-adjusted_EVE|pop-adjusted_ESM1v)=/);

  if ($is_pop) {
    # -------- popEVE branch: extract popEVE_* fields --------
    my ($score)     = $info =~ /(?:^|;)popEVE=([^;]+)/;
    my ($raw_eve)   = $info =~ /(?:^|;)EVE=([^;]+)/;
    my ($esm1v)     = $info =~ /(?:^|;)ESM1v=([^;]+)/;
    my ($padj_eve)  = $info =~ /(?:^|;)pop[-_]adjusted_EVE=([^;]+)/;
    my ($padj_esm)  = $info =~ /(?:^|;)pop[-_]adjusted_ESM1v=([^;]+)/;
    my ($gap)       = $info =~ /(?:^|;)gap_frequency=([^;]+)/;
    my ($gene)      = $info =~ /(?:^|;)gene=([^;]+)/;
    my ($protein)   = $info =~ /(?:^|;)protein=([^;]+)/;
    my ($mutant)    = $info =~ /(?:^|;)mutant=([^;]+)/;

    my %res;
    $res{popEVE_SCORE}              = $score     if defined $score;
    $res{popEVE_EVE}                = $raw_eve   if defined $raw_eve; # avoids collision with the EVE col from the EVE file
    $res{popEVE_ESM1v}              = $esm1v     if defined $esm1v;
    $res{popEVE_pop_adjusted_EVE}   = $padj_eve  if defined $padj_eve;
    $res{popEVE_pop_adjusted_ESM1v} = $padj_esm  if defined $padj_esm;
    $res{popEVE_gap_frequency}      = $gap       if defined $gap;
    $res{popEVE_gene}               = $gene      if defined $gene;
    $res{popEVE_protein}            = $protein   if defined $protein;
    $res{popEVE_mutant}             = $mutant    if defined $mutant;

    my $end = $pos + (length($ref||'') ? length($ref)-1 : 0);
    return { ref=>$ref, alt=>$alt, start=>$pos, end=>$end, result=>\%res };
  }

  # -------- EVE branch --------
  my ($EVE_SCORE) = $info =~ /EVE=([^;]+)/;
  my $class_number = $self->{class_number};
  my ($EVE_CLASS) = $info =~ /Class$class_number=([^;]+)/;

  my $end = $pos + (length($ref||'') ? length($ref)-1 : 0);
  return {
    ref => $ref,
    alt => $alt,
    start => $pos,
    end => $end,
    result => {
      EVE_SCORE => $EVE_SCORE,
      EVE_CLASS => $EVE_CLASS
    }
  };
}

sub get_start {
  return $_[1]->{start};
}

sub get_end {
  # safe fallback if 'end' wasn't set in parse_data
  my $v = $_[1];
  return defined $v->{end} ? $v->{end} : ($v->{start} + (length($v->{ref} // '') ? length($v->{ref}) - 1 : 0));
}

1;