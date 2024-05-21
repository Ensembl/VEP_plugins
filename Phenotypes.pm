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

Phenotypes

=head1 SYNOPSIS

 mv Phenotypes.pm ~/.vep/Plugins

 # Automatically download phenotype annotation files if needed and annotate
 # variants with phenotypes
 ./vep -i variations.vcf --plugin Phenotypes

 # Fetch only gene-associated phenotypes
 ./vep -i variations.vcf --plugin Phenotypes,include_types=Gene

 # Set directory with phenotypes annotations (phenotype annotation file is
 # automatically downloaded if not available in this directory)
 ./vep -i variations.vcf --plugin Phenotypes,dir=${HOME},include_types=Gene

 # Specify a file with phenotypes annotation (file is automatically
 # downloaded and saved with this name if it does not exist)
 ./vep -i variations.vcf --plugin Phenotypes,file=${HOME}/phenotypes.gff.gz,include_types=Gene

=head1 DESCRIPTION

 A VEP plugin that retrieves overlapping phenotype information.

 On the first run for each new version/species/assembly will
 download a GFF-format dump to ~/.vep/Plugins/

 Ensembl provides phenotype annotations mapped to a number of genomic
 feature types, including genes, variants and QTLs.

 This plugin is best used with JSON output format; the output will be
 more verbose and include all available phenotype annotation data and
 metadata.

 For other output formats, only a concatenated list of phenotype
 description strings is returned.

 Several paramters can be set using a key=value system:

 dir            : Path to directory where to look for phenotypes annotation. If
                  the required file does not exist, the file is downloaded and
                  saved in the provided directory (download requires using
                  database or cache mode).

 file           : File path to phenotypes annotation. If the file does
                  not exist, the file is downloaded and saved with this name
                  (download requires using database or cache mode).

 exclude_sources: &-separated list of phenotype sources to exclude. By default,
                  HGMD-PUBLIC and COSMIC annotations are excluded. See
                  http://www.ensembl.org/info/genome/variation/phenotype/sources_phenotype_documentation.html

 include_sources: &-separated list of phenotype sources to include. If defined,
                  exclude_sources is ignored.

 exclude_types  : &-separated list of feature types to exclude: Gene, Variation,
                  QTL, StructuralVariation, SupportingStructuralVariation,
                  RegulatoryFeature. By default, StructuralVariation and 
                  SupportingStructuralVariation annotations are always excluded
                  (due to size issues) and Variation is excluded when annotating
                  structural variants; to get these annotations in all cases, use
                  include_types=StructuralVariation&SupportingStructuralVariation&Variation

 include_types  : &-separated list of feature types to include. If defined,
                  exclude_types is ignored.

 expand_right   : Cache size in bp. By default, annotations 100000bp (100kb)
                  downstream of the initial lookup are cached.

 phenotype_feature : Boolean to report the gene/variation associated with the
                     phenotype (such as overlapping gene or structural
                     variation) and annotation source (default: 0)

 cols           : &-separated list of column and/or attribute names to output
                  from the gff file. The output fields will be ordered in the same
                  way given in cols argument. (default: 'phenotype' or 'source,phenotype,id'
                  if you set phenotype_feature=1)

 id_match       : Return results only if the identifiers matches with the 
                  variant or the gene depending on the type (default: 0)

 Example:

 --plugin Phenotypes,file=${HOME}/phenotypes.gff.gz,include_types=Gene
 --plugin Phenotypes,dir=${HOME},include_types=Gene

=cut

package Phenotypes;

use strict;
use warnings;

use Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin);

# default config
my %CONFIG = (
  exclude_sources => 'HGMD-PUBLIC&COSMIC',
  exclude_types => 'StructuralVariation&SupportingStructuralVariation',
  expand_right => 100000,
  phenotype_feature => 0,
);

my %output_format;
my $char_sep = "|";

my (%cols, @fields_order);

my @FIELDS = qw(seq_region_name source type start end score strand frame attributes comments);

sub new {
  my $class = shift;
  
  my $self = $class->SUPER::new(@_);
  
  my $params_hash = $self->params_to_hash();
  $CONFIG{$_} = $params_hash->{$_} for keys %$params_hash;

  #for REST calls report all data (use json output flag)
  $self->{config}->{output_format} ||= $CONFIG{output_format};

  # get output format
  if ($self->{config}->{output_format}) {
    $output_format{$self->{config}->{output_format}} = 1;
  }
  $char_sep = "+" if ($output_format{'vcf'});

  #DEFAULTS are not refreshed automatically by multiple REST calls unless forced
  my $refresh = 0;
  $refresh = 1 if (exists $CONFIG{species} && $CONFIG{species} ne $self->{config}{species});

  unless($CONFIG{file} && !$refresh) {
    my $pkg = __PACKAGE__;
    $pkg .= '.pm';

    my $config = $self->{config};

    my $species = $config->{species};
    my $version = $config->{db_version} || 'Bio::EnsEMBL::Registry'->software_version;
    my $assembly = $config->{assembly};

    my $dir = $CONFIG{dir};
    if(defined $dir && -d $dir){
      $dir =~ s/\/?$/\//; #ensure dir path string ends in slash
      if( $species eq 'homo_sapiens' || $species eq 'human'){
        $assembly ||= $config->{human_assembly};
        $CONFIG{file} = sprintf("%s_%s_%i_%s.gvf.gz", $dir.$pkg, $species, $version, $assembly);
      } else {
        $CONFIG{file} = sprintf("%s_%s_%i.gvf.gz", $dir.$pkg, $species, $version);
      }
    } else { #assembly value will be automatically populated by VEP script but not by REST server
      $CONFIG{file} = sprintf("%s_%s_%i_%s.gvf.gz", $INC{$pkg}, $species, $version, $assembly);
    }
    $CONFIG{species} = $species;
  }

  # set which columns to output
  if ($CONFIG{cols}){
    %cols = %{ $self->cols };
    @fields_order = split('&', $CONFIG{cols} || '');
  }
  elsif ($CONFIG{phenotype_feature}){
    %cols =  (phenotype => 1, source => 1, id => 1);
    @fields_order = ("phenotype", "source", "id");
  }
  else {
    %cols =  (phenotype => 1);
    @fields_order = ("phenotype");
  }

  $self->generate_phenotype_gff($CONFIG{file}) if !(-e $CONFIG{file}) || (-e $CONFIG{file}.'.lock');

  $self->add_file($CONFIG{file});

  $self->get_user_params();

  return $self;
}

sub feature_types {
  return ['Feature','Intergenic'];
}

sub variant_feature_types {
  return ['BaseVariationFeature'];
}

sub get_header_info {
  my $self = shift;
  return {
    PHENOTYPES => 'Phenotypes associated with overlapping genomic features'
  }
}

sub generate_phenotype_gff {
  my ($self, $file) = @_;

  my $config = $self->{config};
  die("ERROR: File not found (".$file.") and unable to generate GFF file in offline mode\n") if $config->{offline};
  die("ERROR: Not allowed to generate GFF file in rest mode\n") if $config->{rest};
  
  # test bgzip
  die "ERROR: bgzip does not seem to be in your path\n" unless `which bgzip 2>&1` =~ /bgzip$/;

  unless($config->{quiet}) {
    print STDERR "### Phenotypes plugin: Generating GFF file $file from database\n";
    print STDERR "### Phenotypes plugin: This will take some time but it will only run once per species, assembly and release\n";
  }

  my $pfa = $self->{config}->{reg}->get_adaptor($config->{species}, 'variation', 'phenotypefeature');
  die ("ERROR: no variation db found, please check that it exists for this release \n") unless defined $pfa;

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
        GROUP_CONCAT(at.code, "=", concat('"', pfa.value, '"') SEPARATOR '; '),
        GROUP_CONCAT('submitter_name="', sub.description, '"')
      ) AS attribute

      FROM
        seq_region sr,
        source s,
        phenotype p,
        phenotype_feature pf

      LEFT JOIN phenotype_feature_attrib pfa
        ON pf.phenotype_feature_id = pfa.phenotype_feature_id
      LEFT JOIN attrib_type `at`
        ON pfa.attrib_type_id = at.attrib_type_id
      LEFT JOIN submitter sub
        ON (
          pfa.value = sub.submitter_id
          AND pfa.attrib_type_id = at.attrib_type_id
          AND at.code = 'submitter_id'
        )

      WHERE sr.seq_region_id = pf.seq_region_id
      AND s.source_id = pf.source_id
      AND pf.phenotype_id = p.phenotype_id

      GROUP BY pf.phenotype_feature_id
      ORDER BY pf.seq_region_id, pf.seq_region_start, pf.seq_region_end
  }, { mysql_use_result => 1});

  $sth->execute();

  print STDERR "### Phenotypes plugin: Writing to file\n" unless $config->{quiet};
  my $file_sorted = $file;
  $file .= ".tmp";

  my $lock = "$file\.lock";
  open LOCK, ">$lock" or die "ERROR: Unable to write to lock file $lock\n";
  print LOCK "1\n";
  close LOCK;

  open OUT, " | bgzip -c > $file" or die "ERROR: Unable to write to file $file\n";
  print OUT "##gvf-version 1.10\n"; #HEADER

  while(my $row = $sth->fetchrow_arrayref()) {
    # swap start end for insertions
    @$row[3,4] = @$row[4,3] if (@$row[3] > @$row[4]);
    print OUT join("\t", map {defined($_) ? $_ : '.'} @$row)."\n";
  }

  close OUT;

  unlink($lock);

  $sth->finish();

  print STDERR "### Phenotypes plugin: Sorting file with sort\n" unless $config->{quiet};

  system("(zgrep '^#' $file; LC_ALL=C zgrep -v '^#' $file | sort -k1,1 -k4,4n ) | bgzip -c > $file_sorted") and die("ERROR: sort failed\n");

  print STDERR "### Phenotypes plugin: Indexing file with tabix\n" unless $config->{quiet};

  system("tabix -p gff $file_sorted") and die("ERROR: tabix failed\n");

  print STDERR "### Phenotypes plugin: All done!\n" unless $config->{quiet};
}

sub run {
  my ($self, $bvfo) = @_;
  
  my $vf = $bvfo->base_variation_feature;
  $self->{is_sv} = $vf->isa('Bio::EnsEMBL::Variation::StructuralVariationFeature');

  my $tr = $bvfo->transcript;
  my $gene_stable_id = defined $tr ? $tr->{_gene_stable_id} : "";

  # adjust coords for tabix
  my ($s, $e) = ($vf->{start}, $vf->{end});
  ($s, $e) = ($vf->{end}, $vf->{start}) if ($vf->{start} > $vf->{end}); # swap for insertions

  my $data = $self->get_data($vf->{chr}, $s, $e);

  return {} unless $data && scalar @$data;

  my @f_data = @$data;
  if ($CONFIG{id_match}){
    @f_data = grep { $_->{id} eq ($_->{type} eq "Gene" ? $gene_stable_id : $vf->variation_name) } @$data;
  }

  return { PHENOTYPES =>  \@f_data } if ($output_format{'json'} && !$CONFIG{phenotype_feature});

  my %tmp_res_uniq;
  my @result_str = ();
  my @result_data = ();

  foreach my $tmp_data(@f_data) {
    # subset phenotype data columns
    my %tmp = map { $_ => ($tmp_data->{$_} || '') } keys %cols;
    $tmp_data = \%tmp;

    if (!$output_format{'json'}) {
      # replace link characters with _
      foreach (keys %$tmp_data){
        $tmp_data->{$_} =~ tr/ ;,)(/\_\_\_\_\_/;
        $tmp_data->{$_} =~ s/\+/%2B/ if $output_format{'vcf'};
      }

      # report only unique set of fields
      my $record_line = join(",", values %$tmp_data);
      next if defined $tmp_res_uniq{$record_line};
      $tmp_res_uniq{$record_line} = 1;

      push(@result_str, join($char_sep, @$tmp_data{@fields_order}));
    }

    push @result_data, $tmp_data;
  }

  # output options: phenotype_feature + json OR phenotype_feature + vep|vcf|tab
  return {
    PHENOTYPES => defined($output_format{'json'}) ? \@result_data : \@result_str
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
    return undef if $data->{type} && (
                      $self->exclude_types->{$data->{type}} ||
                      # avoid annotating SVs with short variant-associated phenotypes
                      ($self->{is_sv} && $data->{type} eq 'Variation')
                    );
  }

  # parse attributes
  if(defined($data->{attributes})) {
    $data->{attributes} =~ s/^\s+//g;
    
    my %attribs;
    
    foreach my $pair(split /;\s*/, $data->{attributes}) {
      my ($key, $value) = split /\=/, $pair;

      next unless defined($key) and defined($value);
      next if $key eq 'submitter_id'; # if submitter_id exists, submitter_name should also exist and displayed
      
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

sub cols {
  return $_[0]->_generic_inc_exc('cols');
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

