=head1 LICENSE
Copyright [2025-2026] EMBL-European Bioinformatics Institute

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
 gnomADMt

=head1 SYNOPSIS
 mv gnomADMt.pm ~/.vep/Plugins

 # Default usage
 ./vep -i variations.vcf --plugin gnomADMt,file=/path/to/gnomad-chrM.vcf.gz
 # Only return (haplotype-)grouped data for the specified haplogroups
 ./vep -i variations.vcf --plugin gnomADMt,file=/path/to/gnomad-chrM.vcf.gz,hap_filter=HV:A:B
 # Only return population-grouped data for the specified populations
 ./vep -i variations.vcf --plugin gnomADMt,file=/path/to/gnomad-chrM.vcf.gz,fields=pop_AC_hom:pop_AC_het:pop_AF_hom:pop_AF_het:pop_AN,pop_filter=eas:fin:afr

=head1 DESCRIPTION
 An Ensembl VEP plugin that retrieves allele frequencies for mitochondrial variants
 from the gnomAD genomes mitochondrial annotations files, available here:
   https://gnomad.broadinstitute.org/downloads#v3-mitochondrial-dna

 Options are passed to the plugin as key=value pairs:
   file : (mandatory) Tabix-indexed VCF file from gnomAD
   fields : (optional) Colon-separated list of information from mitochondrial variants to
            output (default: 'hap_AC_hom:hap_AC_het:hap_AF_hom:hap_AF_het:hap_AN');
            keyword 'all' can be used to print all fields;
            Available fields are the names of INFO fields.
   hap_filter : (optional) Colon-separated list of haplogroups to limit returned haplogroup field data to.
                The order of the list is preserved in the output. Example: 'HV:A:B'.
                By default, data for all haplogroups is returned for all haplogroup fields (/hap_.+/).
   pop_filter : (optional) Colon-separated list of populations to limit returned population field data to.
                The order of the list is preserved in the output. Example: 'eas:fin:afr'.
                By default, data for all populations is returned for all population fields (/pop_.+/).

 To download the gnomad Mitochondrial allele annotations file in VCF format (GRCh38 based) and the corresponding index file:
    wget https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1/vcf/genomes/gnomad.genomes.v3.1.sites.chrM.vcf.bgz --no-check-certificate
    wget https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1/vcf/genomes/gnomad.genomes.v3.1.sites.chrM.vcf.bgz.tbi --no-check-certificate

 If you use this plugin, please see the terms and data information:
   https://gnomad.broadinstitute.org/terms

 You must have the Bio::DB::HTS module or the tabix utility must be installed
 in your path to use this plugin.
=cut

package gnomADMt;

use strict;
use warnings;

use Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin;
use Bio::EnsEMBL::Variation::Utils::Sequence qw(get_matched_variant_alleles);

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin);

sub _validate_list_selection {
  # Validate string $selected can be split into a list of valid items,
  # each of which must be found in arrayref $available
  my $selected  = shift;
  my $available = shift;
  my $separator = shift || ':';
  my $selection_type = shift || 'field';

  # return all available selection options when using 'all'
  return $available if $selected eq 'all';
  my @selected_items = split(/$separator/, $selected);

  # check if the selected items exist
  my @valid;
  my @invalid;
  for my $item (@selected_items) {
    if ( grep { $_ eq $item } @$available ) {
      push(@valid, $item);
    } else {
      push(@invalid, $item);
    }
  }

  die "ERROR: all ".$selection_type."s given are invalid. Available ".$selection_type."s are:\n" .
    join(", ", @$available)."\n" unless @valid;
  warn "gnomADMt plugin: WARNING: the following ".$selection_type." are not valid and were ignored: ",
    join(", ", @invalid), "\n" if @invalid;

  return \@valid;
}

sub _parse_vcf_description_groups {
  my $description = shift;
  my $groupname = shift;

  my @groups = ();
  if( $description =~ /$groupname order: \[(.+?)\]/ ){
    my $group_string = $1;
    @groups = $group_string =~ /'(\w+)'/g;
  }
  else{
    warn "gnomADMt plugin: WARNING: $groupname order not found in VCF description.\n";
  }

  return \@groups;
}

sub _rewrite_vcf_description_groups {
  my $description = shift;
  my $groupname = shift;
  my $new_groups = shift;

  my $match_pattern = qr/($groupname order: )\[.+?\]/;
  my $subsitute_pattern = '['.join(',', map { "'$_'" } @$new_groups).']';
  $description =~ s/$match_pattern/$1$subsitute_pattern/;

  return $description;
}

sub new {
  my $class = shift;

  my $self = $class->SUPER::new(@_);

  $self->expand_left(0);
  $self->expand_right(0);

  $self->get_user_params();
  my $params = $self->params_to_hash();

  my $file = $params->{file};

  # get INFO fields names from VCF
  my $vcf_file = Bio::EnsEMBL::IO::Parser::VCF4Tabix->open($file);
  my $info = $vcf_file->get_metadata_by_pragma('INFO');
  my $info_ids = [ map { $_->{ID} } @$info ];
  my $info_descriptions = { map { $_->{ID} => $_->{Description} } @$info };

  # check if AF_hom and AF_het fields exists
  if (!defined $info_descriptions->{AF_hom} || !defined $info_descriptions->{AF_het}) {
    die "ERROR: Provided file does not contain AF_hom or AF_het fields. Please provide the file as downloaded from GnomAD.\n"
  }

  $self->add_file($file);

  # Process the requested fields
  my $requested_fields_str = defined $params->{fields} ? $params->{fields} : 'hap_AC_hom:hap_AC_het:hap_AF_hom:hap_AF_het:hap_AN';
  my $fields = _validate_list_selection($requested_fields_str, $info_ids);

  $self->{fields} = $fields;
  $self->{fields_descriptions} = { map { $_ => $info_descriptions->{$_} } @$fields };

  if( defined $params->{hap_filter} ){
    # Check any haplogroup fields are requested.
    # Otherwise, show a warning that haplogroup filtering can only be applied to haplogroup fields.
    if(!grep /hap_.+/, @{$self->{fields}}){
      warn "gnomADMt plugin: WARNING: Haplogroup filtering can only be applied to haplogroup fields.\n";
    }
    else{
      # Parse haplogroup output order
      my $available_haplogroups = _parse_vcf_description_groups($info_descriptions->{'hap_AN'}, 'haplogroup');
      my $haplogroups = _validate_list_selection($params->{hap_filter}, $available_haplogroups, ':', 'haplogroup');

      my @available_haplogroup_idxs = (0..scalar(@$available_haplogroups)-1);
      my @haplogroup_idxs = ();
      for my $haplogroup (@$haplogroups){
        push(@haplogroup_idxs, grep { $available_haplogroups->[$_] eq $haplogroup } @available_haplogroup_idxs);
      }
      $self->{report_haplogroups} = $haplogroups;
      $self->{report_haplogroup_idxs} = \@haplogroup_idxs;
    }
  }

  if( defined $params->{pop_filter} ){
    # Check any population fields are requested.
    # Otherwise, show a warning that population filtering can only be applied to population fields.
    if(!grep /pop_.+/, @{$self->{fields}}){
      warn "gnomADMt plugin: WARNING: Population filtering can only be applied to population fields.\n";
    }
    else{
      # Parse population output order
      my $available_populations = _parse_vcf_description_groups($info_descriptions->{'pop_AN'}, 'population');
      my $populations = _validate_list_selection($params->{pop_filter}, $available_populations, ':', 'population');

      my @available_population_idxs = (0..scalar(@$available_populations)-1);
      my @population_idxs = ();
      for my $population (@$populations){
        push(@population_idxs, grep { $available_populations->[$_] eq $population } @available_population_idxs);
      }
      $self->{report_populations} = $populations;
      $self->{report_population_idxs} = \@population_idxs;
    }
  }

  my $prefix = 'gnomAD';
  $self->{prefix} = $prefix;

  return $self;
}

sub feature_types {
  return ['Feature', 'Intergenic'];
}

sub get_header_info {
  my $self = shift;

  my $prefix = $self->{prefix} . '_';

  my %field_descriptions = map {
    $_ => $self->{fields_descriptions}->{$_}
  } @{$self->{fields}};

  # Report new haplogroup output order based on haplogroup filtering option
  if( defined($self->{report_haplogroups}) ){
    foreach my $field (keys %field_descriptions) {
      if($field =~ /hap_.+/){
        $field_descriptions{$field} = _rewrite_vcf_description_groups($field_descriptions{$field}, 'haplogroup', $self->{report_haplogroups});
      }
    }
  }

  # Report new population output order based on population filtering option
  if( defined($self->{report_populations}) ){
    foreach my $field (keys %field_descriptions) {
      if($field =~ /pop_.+/){
        $field_descriptions{$field} = _rewrite_vcf_description_groups($field_descriptions{$field}, 'population', $self->{report_populations});
      }
    }
  }

  my %header_info = map {
    ''.$prefix.$_ => $field_descriptions{$_};
  } @{$self->{fields}};

  return \%header_info;
}

sub run {
  my ($self, $vfoa) = @_;

  my $vf = $vfoa->variation_feature;

  (my $vf_chr = $vf->{chr}) =~ s/^chr//;

  #Normalize alternative chromosome names
  if( $vf_chr eq 'MT' ){
    $vf_chr = 'M';
  }

  my ($vf_start, $vf_end) = ($vf->{start}, $vf->{end});

  $vf_end = $vf_start if $vf_start > $vf_end;

  my @data = @{ $self->get_data($vf_chr, $vf_start -2, $vf_end) };

  return {} unless @data;

  # Match the relevant alternative allele at the given position
  my $gnomad_freqs = {};
  my $vfoa_variation_feature_seq = $vfoa->variation_feature_seq;

  for my $data_candidate (@data) {
    my $candidate_alt_allele = $data_candidate->{'alt'};

    my $matches = get_matched_variant_alleles(
      {
        ref    => $vf->ref_allele_string,
        alts   => [ $vfoa_variation_feature_seq ],
        pos    => $vf->{start},
        strand => $vf->strand
      },
      {
        ref  => $data_candidate->{ref},
        alts => [ $candidate_alt_allele ],
        pos  => $data_candidate->{start},
      }
    );

    if (@$matches){
      $gnomad_freqs = $data_candidate->{'result'};
      last;
    }
  }

  return $gnomad_freqs;
}

sub _parse_vcf {
  my ($self, $line) = @_;

  my ($chrom, $start, $id, $ref, $alt, $qual, $filter, $info) = split /\t/, $line;
  # VCF-like adjustment of mismatched substitutions for comparison with VEP
  if(length($alt) != length($ref)) {
    my $first_ref = substr($ref, 0, 1);
    my $first_alt = substr($alt, 0, 1);
    if ($first_ref eq $first_alt) {
      $start++;
      $ref = substr($ref, 1) || "-";
      $alt = substr($alt, 1) || "-";
    }
  }

  my %data = (
    'chromosome' => $chrom,
    'start'      => $start,
    'identifier' => $id,
    'alleles'    => $ref.'/'.$alt,
    'ref'        => $ref,
    'alt'        => $alt,
    'quality'    => $qual,
    'filter'     => $filter,
  );
  # fetch data from all INFO fields
  for my $field ( split /;/, $info ) {
    my ($key, $value) = split /=/, $field;
    $data{$key} = $value;
  }

  return \%data;
}

sub parse_data {
  my ($self, $line) = @_;

  my $vcf_data = $self->_parse_vcf($line);

  # Process haplogroup output values (and order) to match haplogroup filtering
  if(defined($self->{report_haplogroups})){
    foreach my $info_field (keys %{$vcf_data}){
      if ($info_field =~ m/^hap_/ && defined($vcf_data->{$info_field})){
        my $separator = ($info_field eq 'hap_hl_hist')?',':'|';
        my $split_pattern = quotemeta($separator);
        my @values = split(/$split_pattern/, $vcf_data->{$info_field});
        my @filtered_values = @values[@{$self->{report_haplogroup_idxs}}];
        $vcf_data->{$info_field} = join($separator, @filtered_values);
      }
    }
  }

  # Process population output values (and order) to match population filtering
  if(defined($self->{report_populations})){
    foreach my $info_field (keys %{$vcf_data}){
      if ($info_field =~ m/^pop_/ && defined($vcf_data->{$info_field})){
        my $separator = ($info_field eq 'pop_hl_hist')?',':'|';
        my $split_pattern = quotemeta($separator);
        my @values = split(/$split_pattern/, $vcf_data->{$info_field});
        my @filtered_values = @values[@{$self->{report_population_idxs}}];
        $vcf_data->{$info_field} = join($separator, @filtered_values);
      }
    }
  }

  # Filter VCF data down to selected fields
  # and prefix the field names
  my %keys = map {
    $_ => join('_', $self->{prefix}, $_)
  }
  @{$self->{fields}};

  my $result = {map { $keys{$_} => $vcf_data->{$_} } @{$self->{fields}}};

  return {
    ref => $vcf_data->{ref},
    alt => $vcf_data->{alt},
    start => $vcf_data->{start},
    result => $result
  }
}

sub get_start {
  return $_[1]->{start};
}

sub get_end {
  return $_[1]->{end};
}

1;

