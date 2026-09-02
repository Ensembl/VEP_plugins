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

 NearestExonJB

=head1 SYNOPSIS

 mv NearestExonJB.pm ~/.vep/Plugins
 ./vep -i variations.vcf --cache --plugin NearestExonJB

=head1 DESCRIPTION

 This is a plugin for the Ensembl Variant Effect Predictor (VEP) that
 finds the nearest exon junction boundary to a coding sequence variant. More than 
 one boundary may be reported if the boundaries are equidistant or if using option
 --intronic.

 The plugin will report the Ensembl identifier of the exon, the distance to the
 exon boundary, the boundary type (start or end of exon) and the total
 length in nucleotides of the exon.

 Various key=value parameters can be altered by passing them to the plugin command:

   max_range : maximum search range in bp (default: 10000)

   intronic  : set to 1 to check nearest exons for intronic variants (default: 0)
               returns the nearest exon upstream and downstream without considering
               the max_range.

   gff3      : path to a tabix-indexed GFF3 file to use for exon locations instead
               of transcript exon objects. Exon identifiers are read from exon_id,
               ID or Name attributes, in that order.

   parent_only: set to 1 with gff3 to exclude exons from readthrough transcripts,
                artifact genes and gene parent feature types not supported by VEP's
                GFF parser (default: 0)


 Parameters are passed e.g.:

   --plugin NearestExonJB,max_range=50000
   --plugin NearestExonJB,max_range=50000,intronic=1
   --plugin NearestExonJB,intronic=1
   --plugin NearestExonJB,gff3=/path/to/genes.gff3.gz
   --plugin NearestExonJB,gff3=/path/to/genes.gff3.gz,parent_only=1

=cut

package NearestExonJB;

use strict;
use warnings;

use Bio::EnsEMBL::Variation::Utils::BaseVepPlugin;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

my $char_sep = "|";

my %CONFIG = (
  max_range => 10000,
  intronic => 0
);

my $TABIX_CACHE_EXPANSION = 1_000_000;
my $TABIX_CACHE_SIZE = 2;

sub new {
  my $class = shift;

  my $self = $class->SUPER::new(@_);

  my $params = $self->params;

  # get output format
  $char_sep = "+" if ($self->{config}->{output_format} eq 'vcf');

  foreach my $param(@$params) {
    my ($key, $val) = split('=', $param);
    die("ERROR: Failed to parse parameter $param\n") unless defined($key) && defined($val);
    if($key eq 'gff3') {
      $self->{gff3_file} = $val;
    }
    elsif($key eq 'parent_only') {
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
  return ['Transcript'];
}

sub get_header_info {
  my $header;

  if($CONFIG{intronic} == 1) {
    $header = 'Nearest Exon Junction Boundary. Format:';
  }
  else {
    $header = 'Nearest Exon Junction Boundary (coding sequence variants only). Format:';
  }

  $header .= join($char_sep, qw(ExonID distance start/end length) );

  return {
    NearestExonJB => $header,
  }
}

sub _enable_tabix_gff3 {
  require Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin;

  no strict 'refs';
  unshift @NearestExonJB::ISA, 'Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin'
    unless grep { $_ eq 'Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin' } @NearestExonJB::ISA;
}

# BaseVepTabixPlugin only recognises TBI indexes in some VEP releases.
# Tabix itself also supports CSI indexes, so accept either format here.
sub check_file {
  my ($self, $file) = @_;

  die("ERROR: No file specified\n") unless $file;

  if($file !~ /tp\:\/\//) {
    die "ERROR: Data file $file not found\n" unless -e $file;
    die "ERROR: Tabix index file $file.tbi or $file.csi not found - perhaps you need to create it first?\n"
      unless -e $file.'.tbi' || -e $file.'.csi';
  }

  return 1;
}

sub run {
  my ($self, $tva) = @_;

  my $vf = $tva->base_variation_feature;
  my $trv = $tva->base_transcript_variation;

  my $loc_string = sprintf("%s:%s-%i-%i", $trv->transcript_stable_id, $vf->{chr} || $vf->seq_region_name, $vf->{start}, $vf->{end});

  if(!exists($self->{_cache}) || !exists($self->{_cache}->{$loc_string})) {
    my $exons = $self->{gff3_file}
      ? $self->_gff3_overlapped_exons($tva, $trv, $vf)
      : $trv->_overlapped_exons; # intronic variants do not overlap any exon
    my %dists;
    my $min = $CONFIG{max_range};

    # For option --intronic, fetch the list of exons with different method
    # Do not take into account the max_range
    if(scalar @{$exons} == 0 && $CONFIG{intronic} == 1) {
      my $intron_numbers = $trv->intron_number();
      my $consequences = join(",", map { $_->SO_term } @{$tva->get_all_OverlapConsequences});

      if(defined $intron_numbers && $consequences =~ /intron/) {
        $exons = $self->{gff3_file}
          ? $self->_gff3_transcript_exons($tva, $trv, $vf)
          : $trv->_sorted_exons;
        my ($intron_number, $total_number) = split(/\//, $intron_numbers);
 
        # Get the number of exons before and after the intron
        my $exon_before = $intron_number;
        my $exon_after = $intron_number + 1;

        my @exons_tmp;
        # In the reverse strand we get the last two exons from the list
        if($tva->transcript->strand < 0) {
          push(@exons_tmp, $exons->[-$exon_before]);
          push(@exons_tmp, $exons->[-$exon_after]);
        }
        else {
          push(@exons_tmp, $exons->[$exon_before -1]);
          push(@exons_tmp, $exons->[$exon_after -1]);
        }

        $exons = \@exons_tmp;
      }
    }

    foreach my $exon (@$exons) {
      my $exon_id = $self->{gff3_file} ? $exon->{stable_id} : $exon->stable_id;
      my $exon_start = $self->{gff3_file} ? $exon->{start} : $exon->seq_region_start;
      my $exon_end = $self->{gff3_file} ? $exon->{end} : $exon->seq_region_end;
      my $exon_length = $self->{gff3_file} ? $exon->{length} : $exon->length;

      my $startD = abs ($vf->start - $exon_start);
      my $endD = abs ($vf->end - $exon_end);
      if ($startD < $endD){
        $dists{$exon_id}{$startD} = 'start';
        $dists{$exon_id}{len} = $exon_length;
        $min = $startD if $min > $startD;
      } elsif ($startD > $endD){
        $dists{$exon_id}{$endD} = 'end';
        $dists{$exon_id}{len} = $exon_length;
        $min = $endD if $min > $endD;
      } else {
        $dists{$exon_id}{$startD} = "start_end";
        $dists{$exon_id}{len} = $exon_length;
        $min = $startD if $min > $startD;
      }
    }

    my @finalRes;
    # For option --intronic, return the closest exons (upstream/dowsntream) from the intron
    if(scalar @{$exons} == 2 && $CONFIG{intronic} == 1) {
      foreach my $exon (keys %dists) {
        my $inner_hash = $dists{$exon};
        my $length_value;
        my $type;
        my $distance_value;

        for my $internal_key (keys %{$inner_hash}) {
          if($internal_key eq "len") {
            $length_value = $inner_hash->{$internal_key};
          }
          else {
            $type = $inner_hash->{$internal_key};
            $distance_value = $internal_key;
          }
        }

        my $string = $exon . $char_sep . $distance_value . $char_sep . $type . $char_sep . $length_value;
        push(@finalRes, $string);
      }
    }
    else {
        # This is the default behaviour of the plugin
        foreach my $exon (keys %dists){
        if (exists $dists{$exon}{$min}) {
          push(@finalRes, $exon.$char_sep.$min.$char_sep.$dists{$exon}{$min}.$char_sep.$dists{$exon}{len})
        }
      }
    }

    $self->{_cache}->{$loc_string} = scalar @finalRes ? join(",", @finalRes) : undef;
  }
  return $self->{_cache}->{$loc_string} ? { NearestExonJB => $self->{_cache}->{$loc_string} } : {};
}

sub _gff3_overlapped_exons {
  my ($self, $tva, $trv, $vf) = @_;

  return $self->_gff3_exons_for_region(
    $vf->{chr} || $vf->seq_region_name,
    $vf->{start},
    $vf->{end},
    $trv->transcript_stable_id
  );
}

sub _gff3_transcript_exons {
  my ($self, $tva, $trv, $vf) = @_;

  my $transcript = $tva->transcript;

  return $self->_gff3_exons_for_region(
    $vf->{chr} || $vf->seq_region_name,
    $transcript->seq_region_start,
    $transcript->seq_region_end,
    $trv->transcript_stable_id
  );
}

sub _gff3_exons_for_region {
  my ($self, $seq_region_name, $start, $end, $transcript_id) = @_;

  my $records = $self->get_data($seq_region_name, $start, $end);

  return [] if $self->{parent_only} &&
    !$self->_gff3_transcript_passes_parent_only_filter($records, $transcript_id);

  return [
    sort {
      $a->{start} <=> $b->{start} || $a->{end} <=> $b->{end} || $a->{stable_id} cmp $b->{stable_id}
    }
    grep {
      $_->{record_type} eq 'exon' &&
      exists $_->{parents}->{$transcript_id}
    } @$records
  ];
}

sub _gff3_transcript_passes_parent_only_filter {
  my ($self, $records, $transcript_id) = @_;

  return $self->{_gff3_parent_only_filter_cache}->{$transcript_id}
    if exists $self->{_gff3_parent_only_filter_cache}->{$transcript_id};

  my ($transcript) = grep {
    $_->{record_type} eq 'transcript' &&
    $_->{stable_id} eq $transcript_id
  } @$records;

  # Exon-only GFF3 files have no metadata to filter against.
  return $self->{_gff3_parent_only_filter_cache}->{$transcript_id} = 1
    unless $transcript;

  return $self->{_gff3_parent_only_filter_cache}->{$transcript_id} = 0
    if $transcript->{readthrough};

  return $self->{_gff3_parent_only_filter_cache}->{$transcript_id} = 1
    unless @{$transcript->{gene_ids}};

  my %genes = map {
    $_->{stable_id} => $_
  } grep {
    $_->{record_type} eq 'gene'
  } @$records;

  my $gff_source = $self->_vep_gff_source;
  my $passes = grep {
    my $gene = $genes{$_};
    my $feature_type = $gene ? ($gene->{feature_type} || '') : '';
    $gene &&
    $gff_source->include_feature_types->{$feature_type} &&
    $gff_source->_record_is_gene({ type => $feature_type }) &&
    ($gene->{biotype} || '') ne 'artifact'
  } @{$transcript->{gene_ids}};

  return $self->{_gff3_parent_only_filter_cache}->{$transcript_id} = $passes ? 1 : 0;
}

sub _gff3_attributes {
  my ($attributes) = @_;

  my %attributes;
  foreach my $attribute(split(/;/, $attributes)) {
    my ($key, $value) = split(/=/, $attribute, 2);
    next unless defined $key && defined $value;
    $attributes{$key} = _gff3_unescape($value);
  }

  return \%attributes;
}

sub _gff3_exon_id {
  my ($attributes) = @_;

  foreach my $key(qw(exon_id ID Name)) {
    next unless exists $attributes->{$key};

    my $stable_id = $attributes->{$key};
    $stable_id =~ s/^exon://;

    return $stable_id;
  }

  return undef;
}

sub _gff3_stable_id {
  my ($attributes, $prefix, @keys) = @_;

  foreach my $key(@keys) {
    next unless exists $attributes->{$key};

    my $stable_id = $attributes->{$key};
    $stable_id =~ s/^\Q$prefix\E://;

    return $stable_id;
  }

  return undef;
}

sub _gff3_parent_transcripts {
  my ($attributes) = @_;

  my %parents;
  foreach my $parent(split(/,/, $attributes->{Parent} || '')) {
    $parent =~ s/^transcript://;
    $parents{$parent} = 1 if length $parent;
  }

  return \%parents;
}

sub _gff3_parent_genes {
  my ($attributes) = @_;

  my @parents;
  foreach my $parent(split(/,/, $attributes->{Parent} || '')) {
    next unless $parent =~ s/^gene://;
    push @parents, $parent if length $parent;
  }

  return \@parents;
}

sub _is_gff3_gene {
  my ($feature_type, $attributes) = @_;

  return 1 if $feature_type eq 'gene' || $feature_type =~ /gene$/;
  return 1 if ($attributes->{ID} || '') =~ /^gene:/;

  return 0;
}

sub _gff3_unescape {
  my ($value) = @_;

  $value =~ s/%([0-9A-Fa-f]{2})/chr(hex($1))/eg;
  return $value;
}

sub parse_data {
  my ($self, $line) = @_;

  return undef if $line =~ /^#/;
  chomp $line;
  return undef unless $line;

  my @cols = split(/\t/, $line);
  return undef unless scalar @cols >= 9;

  my $attributes = _gff3_attributes($cols[8]);

  if($cols[2] eq 'exon') {
    my $stable_id = _gff3_exon_id($attributes);
    my $parents = _gff3_parent_transcripts($attributes);

    return undef unless defined $stable_id && length $stable_id;
    return undef unless scalar keys %$parents;

    return {
      record_type     => 'exon',
      seq_region_name => $cols[0],
      stable_id       => $stable_id,
      start           => $cols[3],
      end             => $cols[4],
      length          => $cols[4] - $cols[3] + 1,
      parents         => $parents,
    };
  }

  return undef unless $self->{parent_only};

  if(_is_gff3_gene($cols[2], $attributes)) {
    my $stable_id = _gff3_stable_id($attributes, 'gene', qw(gene_id ID Name gene_name));
    return undef unless defined $stable_id && length $stable_id;

    return {
      record_type     => 'gene',
      seq_region_name => $cols[0],
      stable_id       => $stable_id,
      feature_type    => $cols[2],
      biotype         => $attributes->{biotype},
      start           => $cols[3],
      end             => $cols[4],
    };
  }

  my $stable_id = _gff3_stable_id($attributes, 'transcript', qw(transcript_id ID));
  my $gene_ids = _gff3_parent_genes($attributes);
  return undef unless defined $stable_id && length $stable_id;
  return undef unless @$gene_ids || ($attributes->{ID} || '') =~ /^transcript:/;

  my %tags = map { $_ => 1 } split(/,/, $attributes->{tag} || '');

  return {
    record_type     => 'transcript',
    seq_region_name => $cols[0],
    stable_id       => $stable_id,
    gene_ids        => $gene_ids,
    readthrough     => $tags{readthrough_transcript} || $attributes->{readthrough_tra} ? 1 : 0,
    start           => $cols[3],
    end             => $cols[4],
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
