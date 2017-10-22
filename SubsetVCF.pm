=head1 NAME

 SubsetVCF

=head1 SYNOPSIS

 ./vep -i variations.vcf --plugin SubsetVCF,file=filepath.vcf.gz,name=myvfc,filter=true,fields=AC*%AN*

=head1 DESCRIPTION

 A VEP plugin to retrieve all requested info fields from a tabix-indexed vcf file.
 When a matching variant is found, info values for all alternatives are returned.
 Additionally, matched alleles are re-formatted to be concordant with the alleles found by VEP.

 Returns for each overlapping allele:
     <name>_POS: Original variant position
     <name>_REF: Original reference
     <name>_ALT: All alternatives with matching alleles minimized
  <name>_<info>: List of requested info values
 
 Parameters:
     name: short name added used as a prefix
     file: path to tabix-index vcf file
   filter: only consider variants marked as 'PASS', 1 or 0 
   fields: info fields to be returned
            '%' can delimit multiple fields
            '*' can be used as a wildcard 

=cut

package SubsetVCF;

use strict;
use warnings;

use Storable qw(dclone);
use Data::Dumper;

use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use Bio::EnsEMBL::Variation::Utils::Sequence qw(get_matched_variant_alleles);
use Bio::EnsEMBL::Variation::Utils::VEP qw(parse_line get_all_consequences);
use Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin;
use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin);

sub simple_vf {
  my ($vf, $filter) = @_;
  my @alleles = split /\//, $vf->{allele_string};
  my $ref = shift @alleles;
  my @line = split /\t/, $vf->{_line};
  my $ret = {
      chr    => $vf->{chr},
      pos    => $vf->{start},
      start  => $vf->{start},
      end    => $vf->{end},
      strand => $vf->{strand},
      alts   => [@alleles],
      line   => [@line],
      ref    => $ref};
  return $filter && $line[6] ne "PASS" ? {} : $ret;
}

sub parse_info {
  my ($line, $valid_fields) = @_;
  my %ret;
  for my $dat (split /;/, $line->[7]) {
    my ($field, $val) = split /=/, $dat;
    if (grep { $field eq $_ } @$valid_fields) {
      $ret{$field} = [split /,/, $val];
    }
  }
  return \%ret;
}

sub new {
  my $class = shift;
  my $self = $class->SUPER::new(@_);

  $self->expand_left(0);
  $self->expand_right(0);
  $self->get_user_params();

  # Add file via parameter hash
  my $params = $self->params_to_hash();
  $self->add_file($params->{file});
  $self->{filter} = $params->{filter};
  $self->{name} = $params->{name};

  # Mung filter to make: # AC*%AN* => AC[^,]+|AN[^,]+
  $params->{fields} =~ s/%/|/g;
  $params->{fields} =~ s/\*/[^,]+/g;

  # Get input file headers
  my %fields;
  my $info_regex = "^##INFO=<ID=($params->{fields}),.*Description=\"([^\"]+).*";
  open HEAD, "tabix -fh $params->{file} 1:1-1 2>&1 | ";
  while(my $line = <HEAD>) {
    next unless $line =~ $info_regex;
    $fields{$1} = $2;
  }
  die "Could not find any valid info fields" if not %fields;

  $self->{fields} = \%fields;
  $self->{valid_fields} = [keys %fields];
  return $self;
}

sub feature_types {
  return ['Feature', 'Intergenic'];
}

sub get_header_info {
  my $self = shift;
  my %ret;
  while (my ($field, $desc) = each %{$self->{fields}}) {
    $ret{"$self->{name}_$field"} = $desc;
  }
  $ret{"$self->{name}_POS"} = "Original POS";
  $ret{"$self->{name}_REF"} = "Original refrance allele";
  $ret{"$self->{name}_ALT"} = "All alternatives with matching alleles minimized";
  return \%ret;
}

sub run {
  my ($self, $tva) = @_;
  my $vf = simple_vf($tva->variation_feature);
  
  # get allele, reverse comp if needed
  my $allele = $tva->variation_feature_seq;
  reverse_comp(\$allele) if $vf->{strand} < 0;

  # Zero-indexing start for tabix and adding 1 to end for VEP indels
  my @data = @{$self->get_data($vf->{chr}, ($vf->{start} - 1), ($vf->{end} + 1))};

  my (%ret, $found_vf, @matches);
  for my $dat (@data) {
    next unless %$dat;
    @matches = @{get_matched_variant_alleles($vf, $dat)};
    if (@matches) {
      $found_vf = dclone $dat;
      last;
    }
  }

  if (@matches) {
    # Convert found alleles to match the VEP output
    my @found_alts = @{$found_vf->{alts}};
    for my $match (@matches) {
      $found_alts[$match->{b_index}] = $match->{a_allele};
    }

    # Parse info fields
    my %found_fields = %{parse_info($found_vf->{line}, $self->{valid_fields})};

    # Organize results
    while (my ($field, $val) = each %found_fields) {
      $ret{"$self->{name}_$field"} = [@$val];
    }

    $ret{"$self->{name}_POS"} = $found_vf->{pos};
    $ret{"$self->{name}_REF"} = $found_vf->{ref};
    $ret{"$self->{name}_ALT"} = [@found_alts];
  }
  return \%ret;
}

sub parse_data {
  my ($self, $line) = @_;
  my ($vf) = @{parse_line({format => 'vcf', minimal => 1}, $line)};
  return simple_vf($vf, $self->{filter});
}

sub get_start {
  return $_[1]->{start};
}

sub get_end {
  return $_[1]->{end};
}

1;
