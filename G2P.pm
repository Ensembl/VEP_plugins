=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016] EMBL-European Bioinformatics Institute

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

 dbNSFP

=head1 SYNOPSIS

 mv G2P.pm ~/.vep/Plugins
 perl variant_effect_predictor.pl -i variations.vcf --plugin G2P,/path/to/G2P.csv.gz

=head1 DESCRIPTION

 A VEP plugin that uses G2P allelic requirements to assess variants in genes
 for potential phenotype involvement.

 The plugin has multiple configuration options, though minimally requires only
 the CSV file of G2P data.

 Options are passed to the plugin as key=value pairs, (defaults in parentheses):

 file        : path to G2P data file, as found at
               http://www.ebi.ac.uk/gene2phenotype/gene2phenotype-webcode/cgi-bin/handler.cgi?show_downloads=all

 maf         : maximum allele frequency for inclusion (0.01)

 maf_key     : allele key to use; by default this is (minor_allele_freq), which
               is the 1000 Genomes global frequency. Choose from:
               1000 genomes: minor_allele_freq,AFR,AMR,EAS,EUR,SAS
               ESP: AA,EA
               ExAC: ExAC,ExAC_AFR,ExAC_AMR,ExAC_Adj,ExAC_EAS,ExAC_FIN,ExAC_NFE,ExAC_OTH,ExAC_SAS

 default_maf : default frequency of the input variant if no frequency data is
               found (0). This determines whether such variants are included;
               the value of 0 forces variants with no frequency data to be
               included as this is considered equivalent to having a frequency
               of 0. Set to 1 (or any value higher than maf) to exclude them.

 types       : SO consequence types to include. Separate multiple values with '&'
               (splice_donor_variant,splice_acceptor_variant,stop_gained,
               frameshift_variant,stop_lost,initiator_codon_variant,
               inframe_insertion,inframe_deletion,missense_variant,
               coding_sequence_variant,start_lost,transcript_ablation,
               transcript_amplification,protein_altering_variant)


 Example:

 --plugin G2P,file=G2P.csv.gz,maf=0.05,maf_key=ExAC_Adj,types=stop_gained&frameshift_variant
 
=cut

package G2P;

use strict;
use warnings;
use Scalar::Util qw(looks_like_number);

use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);

use Bio::EnsEMBL::Variation::Utils::BaseVepPlugin;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

my %DEFAULTS = (

  # vars must have a frequency <= to this to pass
  maf => 0.01,

  # by default we look at the global MAF
  # configure this to use e.g. a continental MAF or an ExAC one
  maf_key => 'minor_allele_freq',

  # if no MAF data is found, default to 0
  # this means absence of MAF data is considered equivalent to MAF=0
  # set to 1 to do the "opposite", i.e. exclude variants with no MAF data
  default_maf => 0,

  # only include variants with these consequence types
  # currently not ontology-resolved, exact term matches only
  types => {map {$_ => 1} qw(splice_donor_variant splice_acceptor_variant stop_gained frameshift_variant stop_lost initiator_codon_variant inframe_insertion inframe_deletion missense_variant
 coding_sequence_variant start_lost transcript_ablation transcript_amplification protein_altering_variant)}
);

sub new {
  my $class = shift;
  
  my $self = $class->SUPER::new(@_);

  my $params = $self->params_to_hash();
  my $file = '';

  # user only supplied file as first param?
  if(!keys %$params) {
    $file = $self->params->[0];
  }
  else {
    $file = $params->{file};

    # process types
    if($params->{types}) {
      $params->{types} = {map {$_ => 1} split(/[\;\&\|]/, $params->{types})};
    }

    # check maf
    if($params->{maf}) {
      die("ERROR: Invalid value for maf: ".$params->{maf}."\n") unless
        looks_like_number($params->{maf}) &&
        ($params->{maf} >= 0 && $params->{maf} <= 1)
    }
  }

  # copy in default params
  $params->{$_} //= $DEFAULTS{$_} for keys %DEFAULTS;
  $self->{user_params} = $params;

  # read data from file
  $self->{gene_data} = $self->read_gene_data_from_file($file);

  # force some config params
  $self->{config}->{individual} //= ['all'];
  $self->{config}->{symbol} = 1;

  if($params->{maf} > 0) {
    $self->{config}->{check_existing} = 1;
    $self->{config}->{check_alleles} = 1;
    $self->{config}->{gmaf} = 1;
  }

  # tell VEP we have a cache so stuff gets shared/merged between forks
  $self->{has_cache} = 1;

  return $self;
}

sub feature_types {
  return ['Transcript'];
}

sub get_header_info {
  my $self = shift;

  return {
    G2P_flag => 'Flags zygosity of valid variants for a G2P gene',
    G2P_complete => 'Indicates this variant completes the allelic requirements for a G2P gene',
  };
}

sub run {
  my ($self, $tva, $line) = @_;

  my $params = $self->{user_params};

  # only interested if we know the zygosity
  my $zyg = $line->{Extra}->{ZYG};
  return {} unless $zyg;

  # only interested in given gene set
  my $tr = $tva->transcript;
  my $gene_symbol = $tr->{_gene_symbol} || $tr->{_gene_hgnc};
  my $gene_data = $self->gene_data($gene_symbol);
  return {} unless $gene_data;

  my $ar = $gene_data->{'allelic requirement'};
  return {} unless $ar && $ar =~ /^(mono|bi)allelic$/;
  
  # limit by type
  return {} unless grep {$self->{user_params}->{types}->{$_->SO_term}} @{$tva->get_all_OverlapConsequences};
  
  # limit by MAF
  if($params->{maf} > 0) {
    return {} unless $self->get_freq($tva) < $params->{maf};
  }

  my %return = (
    G2P_flag => $zyg
  );

  my $cache = $self->{cache}->{$tr->stable_id} ||= {};

  # create a name for the VF to cache by
  my $vf = $tva->base_variation_feature;
  my $vf_name = $vf->variation_name || $vf->{start}.'_'.$vf->{allele_string};
  delete $cache->{$vf_name} if exists($cache->{$vf_name});

  # biallelic genes require >=1 hom or >=2 hets
  if($ar eq 'biallelic') {

    # homozygous, report complete
    if(uc($zyg) eq 'HOM') {
      $return{G2P_complete} = 1;
    }

    # heterozygous
    # we need to cache that we've observed one
    elsif(uc($zyg) eq 'HET') {
      if(scalar keys %$cache) {
        $return{G2P_complete} = 1;
      }
      $cache->{$vf_name} = 1;
    }
  }
  # monoallelic genes require only one allele
  elsif($ar eq 'monoallelic') {
    $return{G2P_complete} = 1;
  }
  else {
    return {};
  }

  return \%return;
}

# read G2P CSV dump
# as from http://www.ebi.ac.uk/gene2phenotype/gene2phenotype-webcode/cgi-bin/handler.cgi?show_downloads=all
sub read_gene_data_from_file {
  my $self = shift;
  my $file = shift;
  my $delimiter = shift;

  my (@headers, %gene_data);
  
  die("ERROR: No file specified or could not read from file ".($file || '')."\n") unless $file && -e $file;

  # allow file to be (b)gzipped
  if($file =~ /\.gz/) {
    open FILE, "gzip -dc $file |";
  }
  else {
    open FILE, $file;
  }

  # this regexp allows for nested ",", e.g.
  # item,description
  # cheese,"salty,delicious"
  my $re = qr/(?: "\( ( [^()""]* ) \)" |  \( ( [^()]* ) \) |  " ( [^"]* ) " |  ( [^,]* ) ) , \s* /x;

  while(<FILE>) {
    chomp;

    my @split = grep defined, "$_," =~ /$re/g;

    unless(@headers) {
      @headers = map {s/\"//g; $_} @split;
    }
    else {
      my %tmp = map {$headers[$_] => $split[$_]} (0..$#split);

      die("ERROR: Gene symbol column not found\n$_\n") unless $tmp{"gene symbol"};

      # store data as hash keyed on gene symbol
      $gene_data{$tmp{"gene symbol"}} = \%tmp;
    }
  }

  close FILE;

  return \%gene_data;
}

# return either whole gene data hash or one gene's data
# this should allow updates to this plugin to e.g. query a REST server, for example
sub gene_data {
  my ($self, $gene_symbol) = @_;
  return $gene_symbol ? $self->{gene_data}->{$gene_symbol} : $self->{gene_data};
}

# get the relevant allele frequency
sub get_freq {
  my $self = shift;
  my $tva = shift;

  my $vf     = $tva->base_variation_feature;
  my $allele = $tva->variation_feature_seq;
  reverse_comp(\$allele) if $vf->{strand} < 0;

  my $cache = $vf->{_g2p_freqs} ||= {};

  # cache it on VF...
  if(!exists($cache->{allele})) {

    my $freq    = $self->{user_params}->{default_maf};
    my $maf_key = $self->{user_params}->{maf_key};

    foreach my $ex(@{$vf->{existing} || []}) {
      
      if($maf_key eq 'minor_allele_freq') {
        if(($ex->{minor_allele} || '') eq $allele && defined($ex->{minor_allele_freq})) {
          $freq = $ex->{minor_allele_freq};
          last;
        }
      }
      else {
        foreach my $pair(split(',', $ex->{$maf_key} || '')) {
          my ($a, $f) = split(':', $pair);

          if(($a || '') eq $allele && defined($f)) {
            $freq = $f;
            last;
          }
        }
      }
    }

    $cache->{$allele} = $freq;
  }

  return $cache->{$allele};
}

1;

