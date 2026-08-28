=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2026] EMBL-European Bioinformatics Institute

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

 NearestGene

=head1 SYNOPSIS

 mv NearestGene.pm ~/.vep/Plugins
 ./vep -i variations.vcf --cache --plugin NearestGene

=head1 DESCRIPTION

 This is a plugin for the Ensembl Variant Effect Predictor (VEP) that
 finds the nearest gene(s) to a non-genic variant. More than one gene
 may be reported if the genes overlap the variant or if genes are
 equidistant or if option 'both_directions' is used.

 Various key=value parameters can be altered by passing them to the plugin command:

   limit           : limit the number of genes returned (default: 1)
   range           : initial search range in bp (default: 1000)
   max_range       : maximum search range in bp (default: 50000)
   both_directions : return the nearest genes upstream and downstream of the variant
                     this option overwrites the limit to 1
                     note that the max_range affects the search range in both directions
                     results include distance and direction, with direction relative to
                     the forward strand
   gff3            : path to a tabix-indexed GFF3 file to use for gene locations instead
                     of a database connection. Gene identifiers are read from gene_id, ID,
                     Name or gene_name attributes, in that order.
   regulatory      : set to 0 to restrict the plugin to Intergenic consequences and prevent
                     VEP from automatically enabling regulatory and motif annotation
                     (default: 1)
   parent_only      : set to 1 with gff3 to consider only gene parent feature types
                     supported by VEP's GFF parser, excluding artifact genes and genes
                     represented only by readthrough transcripts (default: 0)

 Parameters are passed e.g.:

 --plugin NearestGene,limit=3,max_range=50000
 --plugin NearestGene,max_range=50000,both_directions=1
 --plugin NearestGene,gff3=/path/to/genes.gff3.gz
 --plugin NearestGene,gff3=/path/to/genes.gff3.gz,regulatory=0
 --plugin NearestGene,gff3=/path/to/genes.gff3.gz,parent_only=1

 By default, this plugin requires a database connection. It cannot be run in offline mode
 i.e. using the --offline flag. When the gff3 option is provided, gene locations are read
 from the tabix-indexed GFF3 file and the plugin can be run without a database connection.

=cut

package NearestGene;

use strict;
use warnings;

use Bio::EnsEMBL::Variation::Utils::BaseVepPlugin;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

my $char_sep = "|";

my %CONFIG = (
  limit => 1,
  range => 1000,
  max_range => 50000,
);

my $TABIX_CACHE_EXPANSION = 1_000_000;
my $TABIX_CACHE_SIZE = 2;

sub new {
  my $class = shift;
  
  my $self = $class->SUPER::new(@_);

  my $params = $self->params;

  # get output format
  $char_sep = ":" if ($self->{config}->{output_format} eq 'vcf');

  foreach my $param(@$params) {
    my ($key, $val) = split('=', $param);

    die("ERROR: Failed to parse parameter $param\n") unless defined($key) && defined($val);
    if($key eq "both_directions") {
      $self->{both_directions} = 1;
    }
    elsif($key eq "gff3") {
      $self->{gff3_file} = $val;
    }
    elsif($key eq "regulatory") {
      die("ERROR: regulatory must be either 0 or 1\n") unless $val =~ /^[01]$/;
      $self->{include_regulatory} = $val;
    }
    elsif($key eq "parent_only") {
      die("ERROR: parent_only must be either 0 or 1\n") unless $val =~ /^[01]$/;
      $self->{parent_only} = $val;
    }
    else {
      $CONFIG{$key} = $val;
    }
  }

  if($self->{gff3_file}) {
    $self->_enable_tabix_gff3;
    # BaseVepTabixPlugin does not clamp expanded starts to 1; left expansion
    # near a chromosome boundary would therefore produce an invalid region.
    $self->expand_left(0);
    $self->expand_right($TABIX_CACHE_EXPANSION);
    $self->cache_size($TABIX_CACHE_SIZE);
    $self->add_file($self->{gff3_file});
  }

  return $self;
}

sub feature_types {
  my $self = shift;

  return ['Intergenic'] if defined($self->{include_regulatory}) && !$self->{include_regulatory};
  return ['Intergenic','MotifFeature','RegulatoryFeature'];
}

sub variant_feature_types {
  return ['BaseVariationFeature'];
}

sub get_header_info {
  return {
    NearestGene => "Ensembl identifier of nearest gene; with both_directions, also includes distance and forward-strand direction"
  };
}

sub _enable_tabix_gff3 {
  require Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin;

  no strict 'refs';
  unshift @NearestGene::ISA, 'Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin'
    unless grep { $_ eq 'Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin' } @NearestGene::ISA;
}

sub run {
  my ($self, $vfoa) = @_;
  
  if ($self->{config}->{offline} && !$self->{gff3_file}) {
      die("ERROR: the plugin NearestGene does not work in --offline mode\n");
  }

  my $vf = $vfoa->base_variation_feature;
  my $loc_string = sprintf("%s:%i-%i", $vf->{chr} || $vf->seq_region_name, $vf->{start}, $vf->{end});
  
  if(!exists($self->{_cache}) || !exists($self->{_cache}->{$loc_string})) {
    if($self->{gff3_file}) {
      my @result;
      my $seq_region_name = $vf->{chr} || $vf->seq_region_name;

      @result = $self->_nearest_gff3_genes($seq_region_name, $vf->{start}, $vf->{end});

      $self->{_cache}->{$loc_string} = scalar @result ? join(",", @result) : undef;
      return $self->{_cache}->{$loc_string} ? { NearestGene => $self->{_cache}->{$loc_string} } : {};
    }

    $self->{config}->{ga} = $self->{config}->{reg}->get_adaptor($self->{config}->{species}, $self->{config}->{core_type}, 'gene');
    $self->{ga} ||= $self->{config}->{ga};
    die("ERROR: Could not get gene adaptor\n") unless $self->{ga};

    my %opts = map {'-'.$_ => $CONFIG{$_}} keys %CONFIG;
    $opts{-feature} = $vf;

    my @result;

    if($self->{both_directions}) {
      # Overwrite the limit - we want to return only one gene on each direction
      $opts{-limit} = 1;

      # Get upstream genes
      $opts{-upstream} = "upstream";

      my $list_of_genes = $self->{ga}->fetch_all_by_outward_search(%opts);

      for my $gene_result (@{$list_of_genes}){
        my $gene_id = @{$gene_result}[0]->stable_id;
        my $distance = @{$gene_result}[1];
        push(@result, join($char_sep, $gene_id, $distance, 'upstream'))
      }

      # Get downstream genes
      delete $opts{-upstream};
      $opts{-downstream} = "downstream";

      my $list_of_genes_2 = $self->{ga}->fetch_all_by_outward_search(%opts);

      for my $gene_result (@{$list_of_genes_2}){
        my $gene_id = @{$gene_result}[0]->stable_id;
        my $distance = @{$gene_result}[1];
        push(@result, join($char_sep, $gene_id, $distance, 'downstream'))
      }
    }
    else {
      # Default behaviour
      @result = map {$_->[0]->stable_id} @{
        $self->{ga}->fetch_all_by_outward_search(%opts)
      };
    }

    $self->{_cache}->{$loc_string} = scalar @result ? join(",", @result) : undef;
  }

  return $self->{_cache}->{$loc_string} ? { NearestGene => $self->{_cache}->{$loc_string} } : {};
}

sub _gff3_gene_id {
  my ($attributes) = @_;

  my $parsed_attributes = _gff3_attributes($attributes);

  foreach my $key(qw(gene_id ID Name gene_name)) {
    next unless exists $parsed_attributes->{$key};

    my $stable_id = $parsed_attributes->{$key};
    $stable_id =~ s/^gene://;

    return $stable_id;
  }

  return undef;
}

sub _gff3_attributes {
  my ($attributes) = @_;

  my %parsed_attributes;
  foreach my $attribute(split(/;/, $attributes)) {
    my ($key, $value) = split(/=/, $attribute, 2);
    next unless defined $key && defined $value;
    $value = _gff3_unescape($value);
    $parsed_attributes{$key} = exists $parsed_attributes{$key}
      ? $parsed_attributes{$key}.','.$value
      : $value;
  }

  return \%parsed_attributes;
}

sub _is_gff3_gene {
  my ($feature_type, $attributes) = @_;

  return 1 if $feature_type eq 'gene' || $feature_type =~ /gene$/;
  return 1 if $attributes =~ /(?:^|;)ID=gene:/;

  return 0;
}

sub _gff3_unescape {
  my ($value) = @_;

  $value =~ s/%([0-9A-Fa-f]{2})/chr(hex($1))/eg;
  return $value;
}

sub _nearest_gff3_genes {
  my ($self, $seq_region_name, $vf_start, $vf_end) = @_;

  my $slice_start = $self->{both_directions} ? $vf_start - $CONFIG{max_range} : $vf_start;
  $slice_start = 1 if $slice_start < 1;

  my $slice_end = $vf_end + $CONFIG{max_range};
  my $records = $self->get_data($seq_region_name, $slice_start, $slice_end);
  my $genes = $self->_gff3_parent_only_filtered_genes($seq_region_name, $records);

  return $self->_nearest_gff3_genes_both_directions($genes, $vf_start, $vf_end) if $self->{both_directions};

  my @matches = sort {
    $a->[1] <=> $b->[1] || $a->[0]->{start} <=> $b->[0]->{start} || $a->[0]->{stable_id} cmp $b->[0]->{stable_id}
  } grep {
    $_->[0]->{start} >= $vf_end && $_->[1] <= $CONFIG{max_range}
  } map {
    [$_, _distance_to_feature($_, $vf_start, $vf_end)]
  } @$genes;

  if(@matches) {
    my $last_index = _last_index($CONFIG{limit}, scalar @matches);
    my $max_distance = $matches[$last_index]->[1];
    @matches = grep { $_->[1] <= $max_distance } @matches;
  }

  return map { $_->[0]->{stable_id} } @matches;
}

sub _gff3_parent_only_filtered_genes {
  my ($self, $seq_region_name, $records) = @_;

  my @genes = grep { $_->{record_type} eq 'gene' } @$records;
  return \@genes unless $self->{parent_only};

  my %has_non_readthrough_transcript;
  foreach my $record(grep { $_->{record_type} eq 'transcript' } @$records) {
    next if $record->{readthrough};
    $has_non_readthrough_transcript{$_} = 1 for @{$record->{gene_ids}};
  }

  @genes = grep {
    $self->_gff3_gene_passes_parent_only_filter(
      $seq_region_name,
      $_,
      $has_non_readthrough_transcript{$_->{stable_id}},
    )
  } @genes;

  return \@genes;
}

sub _gff3_gene_passes_parent_only_filter {
  my ($self, $seq_region_name, $gene, $has_non_readthrough_transcript) = @_;

  my $stable_id = $gene->{stable_id};
  return $self->{_gff3_parent_only_filter_cache}->{$stable_id}
    if exists $self->{_gff3_parent_only_filter_cache}->{$stable_id};

  my $gff_source = $self->_vep_gff_source;
  my $feature_type = $gene->{feature_type} || '';

  return $self->{_gff3_parent_only_filter_cache}->{$stable_id} = 0
    unless $gff_source->include_feature_types->{$feature_type}
      && $gff_source->_record_is_gene({ type => $feature_type });

  return $self->{_gff3_parent_only_filter_cache}->{$stable_id} = 0
    if ($gene->{biotype} || '') eq 'artifact';

  return $self->{_gff3_parent_only_filter_cache}->{$stable_id} = 1
    if $has_non_readthrough_transcript;

  my $records = $self->get_data(
    $seq_region_name,
    $gene->{start},
    $gene->{end},
  );

  my @transcripts = grep {
    $_->{record_type} eq 'transcript' &&
    grep { $_ eq $stable_id } @{$_->{gene_ids}}
  } @$records;

  my $passes = !@transcripts || grep { !$_->{readthrough} } @transcripts;
  return $self->{_gff3_parent_only_filter_cache}->{$stable_id} = $passes ? 1 : 0;
}

sub _nearest_gff3_genes_both_directions {
  my ($self, $genes, $vf_start, $vf_end) = @_;

  my (@upstream, @downstream);

  foreach my $gene(@$genes) {
    my $distance = _distance_to_feature($gene, $vf_start, $vf_end);
    next if $distance > $CONFIG{max_range};

    if($gene->{end} < $vf_start) {
      push @upstream, [$gene, $distance];
    }
    elsif($gene->{start} > $vf_end) {
      push @downstream, [$gene, $distance];
    }
    else {
      push @upstream, [$gene, $distance];
      push @downstream, [$gene, $distance];
    }
  }

  @upstream = sort { $a->[1] <=> $b->[1] || $b->[0]->{end} <=> $a->[0]->{end} || $a->[0]->{stable_id} cmp $b->[0]->{stable_id} } @upstream;
  @downstream = sort { $a->[1] <=> $b->[1] || $a->[0]->{start} <=> $b->[0]->{start} || $a->[0]->{stable_id} cmp $b->[0]->{stable_id} } @downstream;

  return (
    (scalar @upstream ? join($char_sep, $upstream[0]->[0]->{stable_id}, $upstream[0]->[1], 'upstream') : ()),
    (scalar @downstream ? join($char_sep, $downstream[0]->[0]->{stable_id}, $downstream[0]->[1], 'downstream') : ()),
  );
}

sub _distance_to_feature {
  my ($gene, $vf_start, $vf_end) = @_;

  return $vf_start - $gene->{end} if $gene->{end} < $vf_start;
  return $gene->{start} - $vf_end if $gene->{start} > $vf_end;
  return 0;
}

sub _last_index {
  my ($limit, $count) = @_;

  return -1 unless $count;
  return $count - 1 if $count < $limit;
  return $limit - 1;
}

sub parse_data {
  my ($self, $line) = @_;

  return undef if $line =~ /^#/;
  chomp $line;
  return undef unless $line;

  my @cols = split(/\t/, $line);
  return undef unless scalar @cols >= 9;

  my ($seq_region_name, $feature_type, $start, $end, $attributes) = @cols[0, 2, 3, 4, 8];
  my $parsed_attributes = _gff3_attributes($attributes);

  if(_is_gff3_gene($feature_type, $attributes)) {
    my $stable_id = _gff3_gene_id($attributes);
    return undef unless defined $stable_id && length $stable_id;

    return {
      record_type     => 'gene',
      seq_region_name => $seq_region_name,
      stable_id       => $stable_id,
      feature_type    => $feature_type,
      biotype         => $parsed_attributes->{biotype},
      start           => $start,
      end             => $end,
    };
  }

  return undef unless $self->{parent_only};

  my @gene_ids = map {
    my $gene_id = $_;
    $gene_id =~ s/^gene://;
    $gene_id;
  } grep {
    /^gene:/
  } split(/,/, $parsed_attributes->{Parent} || '');

  return undef unless @gene_ids;
  return undef unless $parsed_attributes->{transcript_id} ||
    ($parsed_attributes->{ID} || '') =~ /^transcript:/;

  my %tags = map { $_ => 1 } split(/,/, $parsed_attributes->{tag} || '');

  return {
    record_type     => 'transcript',
    seq_region_name => $seq_region_name,
    gene_ids        => \@gene_ids,
    readthrough     => $tags{readthrough_transcript} || $parsed_attributes->{readthrough_tra} ? 1 : 0,
    start           => $start,
    end             => $end,
  };
}

sub _vep_gff_source {
  my $self = shift;

  require Bio::EnsEMBL::VEP::AnnotationSource::File::GFF;

  return $self->{_vep_gff_source} ||= bless(
    { gff_type => 'transcript' },
    'Bio::EnsEMBL::VEP::AnnotationSource::File::GFF',
  );
}

sub get_start {
  return $_[1]->{start};
}

sub get_end {
  return $_[1]->{end};
}

1;
