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

 Phenotypes

=head1 SYNOPSIS

 mv Phenotypes.pm ~/.vep/Plugins
 perl variant_effect_predictor.pl -i variations.vcf --plugin Phenotypes

=head1 DESCRIPTION

 A VEP plugin that retrieves overlapping phenotype information
 
=cut

package Phenotypes;

use strict;
use warnings;

use Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin);

# default config
my %CONFIG = (
  exclude_sources => 'HGMD_PUBLIC&COSMIC',
  exclude_types => 'StructuralVariation&SupportingStructuralVariation',
  expand_right => 100000,
);

my @FIELDS = qw(seqname source type start end score strand frame attributes comments);

sub new {
  my $class = shift;
  
  my $self = $class->SUPER::new(@_);
  
  my $params_hash = $self->params_to_hash();
  $CONFIG{$_} = $params_hash->{$_} for keys %$params_hash;

  unless($CONFIG{file}) {
    my $pkg = __PACKAGE__;
    $pkg .= '.pm';

    my $config = $self->{config};

    my $species = $config->{species};
    my $version = $config->{db_version} || $config->{reg}->software_version;
    my $assembly = $config->{assembly};

    $CONFIG{file} = sprintf("%s_%s_%i_%s.bed.gz", $INC{$pkg}, $species, $version, $assembly);
  }

  $self->generate_phenotype_gff($CONFIG{file}) if !(-e $CONFIG{file}) || (-e $CONFIG{file}.'.lock');

  $self->add_file($CONFIG{file});

  $self->get_user_params();

  return $self;
}

sub feature_types {
  return ['Feature','Intergenic'];
}

sub get_header_info {
  my $self = shift;
  return {
    phenotypes => 'Phenotypes associated with overlapping genomic features'
  }
}

sub generate_phenotype_gff {
  my ($self, $file) = @_;

  my $config = $self->{config};
  die("ERROR: Unable to generate GFF file in offline mode\n") if $config->{offline};
  
  # test bgzip
  die "ERROR: bgzip does not seem to be in your path\n" unless `which bgzip 2>&1` =~ /bgzip$/;

  unless($config->{quiet}) {
    print STDERR "### Phenotypes plugin: Generating GFF file $file from database\n";
    print STDERR "### Phenotypes plugin: This will take some time but it will only run once per species, assembly and release\n";
  }

  my $pfa = $self->{config}->{reg}->get_adaptor($config->{species}, 'variation', 'phenotypefeature');

  print STDERR "### Phenotypes plugin: Querying database\n" unless $config->{quiet};

  my $sth = $pfa->dbc->prepare(qq{
    SELECT
      sr.name AS seqname,
      REPLACE(s.name, " ", "_") AS source,
      pf.type AS feature,
      pf.seq_region_start AS start,
      pf.seq_region_end AS end,
      NULL AS score,
      IF(pf.seq_region_strand = 1, '+', '-') AS strand,
      NULL AS frame,
      CONCAT_WS('; ',
        CONCAT('id=', pf.object_id),
        CONCAT('phenotype="', REPLACE(p.description, '"', ''), '"'),
        GROUP_CONCAT(at.code, "=", concat('"', pfa.value, '"') SEPARATOR '; ')
      ) AS attribute

      FROM
        seq_region sr,
        source s,
        phenotype p,
        phenotype_feature pf

      LEFT JOIN (phenotype_feature_attrib pfa, attrib_type `at`)
      ON (
        pf.phenotype_feature_id = pfa.phenotype_feature_id
        AND pfa.attrib_type_id = at.attrib_type_id
      )

      WHERE sr.seq_region_id = pf.seq_region_id
      AND s.source_id = pf.source_id
      AND pf.phenotype_id = p.phenotype_id

      GROUP BY pf.phenotype_feature_id
      ORDER BY pf.seq_region_id, pf.seq_region_start, pf.seq_region_end
  }, { mysql_use_result => 1});

  $sth->execute();

  print STDERR "### Phenotypes plugin: Writing to file\n" unless $config->{quiet};

  my $lock = "$file\.lock";
  open LOCK, ">$lock" or die "ERROR: Unable to write to lock file $lock\n";
  print LOCK "1\n";
  close LOCK;

  open OUT, " | bgzip -c > $file" or die "ERROR: Unable to write to file $file\n";

  while(my $row = $sth->fetchrow_arrayref()) {
    print OUT join("\t", map {defined($_) ? $_ : '.'} @$row)."\n";
  }

  close OUT;

  unlink($lock);

  $sth->finish();

  print STDERR "### Phenotypes plugin: Indexing file with tabix\n" unless $config->{quiet};

  system("tabix -p gff $file") and die("ERROR: tabix failed\n");

  print STDERR "### Phenotypes plugin: All done!\n" unless $config->{quiet};
}

sub run {
  my ($self, $tva) = @_;
  
  my $vf = $tva->variation_feature;
  
  # adjust coords to account for VCF-like storage of indels
  my ($s, $e) = ($vf->{start} - 1, $vf->{end} + 1);
  
  my $data = $self->get_data($vf->{chr}, $s, $e);

  return {} unless $data && scalar @$data;

  return {
    phenotypes => $self->{config}->{json} ? $data : join(",", map {$_->{phenotype} =~ s/\s+/\_/gr} @$data)
  };
}

sub parse_data {
  my ($self, $line) = @_;

  my @split = split /\t/, $line;
  
  my $data;
  
  # parse split data into hash
  for my $i(0..$#split) {
    $data->{$FIELDS[$i]} = $split[$i];
  }

  my $inc_sources = $self->include_sources;
  if(scalar keys %$inc_sources) {
    return undef if $data->{source} && !$inc_sources->{$data->{source}};
  }
  else {
    return undef if $data->{source} && $self->exclude_sources->{$data->{source}};
  }

  my $inc_types = $self->include_types;
  if(scalar keys %$inc_types) {
    return undef if $data->{type} && !$inc_types->{$data->{type}};
  }
  else {
    return undef if $data->{type} && $self->exclude_types->{$data->{type}};
  }

  # parse attributes
  if(defined($data->{attributes})) {
    $data->{attributes} =~ s/^\s+//g;
    
    my %attribs;
    
    foreach my $pair(split /;\s*/, $data->{attributes}) {
      my ($key, $value) = split /\=/, $pair;

      next unless defined($key) and defined($value);
      
      # remove quote marks
      $value =~ s/\"//g;

      # avoid overwriting if an attrib key duplicates a main key
      $key = 'attrib_'.$key if exists($data->{lc($key)});
      
      # lowercase key to reduce chances of mess up!
      $data->{lc($key)} = $value;
    }
    
    delete $data->{attributes};
  }

  # delete empty
  map {delete $data->{$_}} grep {$data->{$_} eq '.' || $data->{$_} eq ''} keys %$data;

  return $data;
}

sub get_start {
  return $_[1]->{start};
}

sub get_end {
  return $_[1]->{end};
}

sub exclude_sources {
  return $_[0]->_generic_inc_exc('exclude_sources');
}

sub exclude_types {
  return $_[0]->_generic_inc_exc('exclude_types');
}

sub include_sources {
  return $_[0]->_generic_inc_exc('include_sources');
}

sub include_types {
  return $_[0]->_generic_inc_exc('include_types');
}

sub _generic_inc_exc {
  my ($self, $key) = @_;

  if(!exists($self->{'_'.$key})) {
    my %exc = map {$_ => 1} split('&', $CONFIG{$key} || '');
    $self->{'_'.$key} = \%exc;
  }

  return $self->{'_'.$key};
}

1;

