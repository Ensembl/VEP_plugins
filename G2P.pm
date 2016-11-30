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

 file            : path to G2P data file, as found at
                   http://www.ebi.ac.uk/gene2phenotype/downloads

 maf_monoallelic : maximum allele frequency for inclusion for monoallelic genes (0.0001)

 maf_biallelic   : maximum allele frequency for inclusion for biallelic genes (0.005)

 maf_key         : allele key to use; by default this is (minor_allele_freq), which
                   is the 1000 Genomes global frequency. Choose from:
                   1000 genomes: minor_allele_freq,AFR,AMR,EAS,EUR,SAS
                   ESP: AA,EA
                   ExAC: ExAC,ExAC_AFR,ExAC_AMR,ExAC_Adj,ExAC_EAS,ExAC_FIN,ExAC_NFE,ExAC_OTH,ExAC_SAS

 default_maf     : default frequency of the input variant if no frequency data is
                   found (0). This determines whether such variants are included;
                   the value of 0 forces variants with no frequency data to be
                   included as this is considered equivalent to having a frequency
                   of 0. Set to 1 (or any value higher than maf) to exclude them.

 types           : SO consequence types to include. Separate multiple values with '&'
                   (splice_donor_variant,splice_acceptor_variant,stop_gained,
                   frameshift_variant,stop_lost,initiator_codon_variant,
                   inframe_insertion,inframe_deletion,missense_variant,
                   coding_sequence_variant,start_lost,transcript_ablation,
                   transcript_amplification,protein_altering_variant)

  exac_file      : ExAC data file, Visit ftp://ftp.broadinstitute.org/pub/ExAC_release/current
                   to download the latest ExAC VCF
  
  log_dir        : write stats to log files in log_dir 
  

 Example:

 --plugin G2P,file=G2P.csv.gz,maf_monoallelic=0.05,maf_key=ExAC_Adj,types=stop_gained&frameshift_variant
 
=cut

package G2P;

use strict;
use warnings;


use ExAC;

use Scalar::Util qw(looks_like_number);

use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);

use Bio::EnsEMBL::Variation::Utils::BaseVepPlugin;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

my %DEFAULTS = (

  # vars must have a frequency <= to this to pass
  maf => 0.001,
  maf_monoallelic => 0.0001,
  maf_biallelic => 0.005, 

  # by default we look at the global MAF
  # configure this to use e.g. a continental MAF or an ExAC one
  maf_key => 'minor_allele_freq',

  # if no MAF data is found, default to 0
  # this means absence of MAF data is considered equivalent to MAF=0
  # set to 1 to do the "opposite", i.e. exclude variants with no MAF data
  default_maf => 0,

  # only include variants with these consequence types
  # currently not ontology-resolved, exact term matches only
  types => {map {$_ => 1} qw(splice_donor_variant splice_acceptor_variant stop_gained frameshift_variant stop_lost initiator_codon_variant inframe_insertion inframe_deletion missense_variant coding_sequence_variant start_lost transcript_ablation transcript_amplification protein_altering_variant)},

);

my $maf_key_2_population_name = {
  minor_allele_freq => '1000GENOMES:phase_3:ALL',
  AFR => '1000GENOMES:phase_3:AFR',
  AMR => '1000GENOMES:phase_3:AMR',
  EAS => '1000GENOMES:phase_3:EAS',
  EUR => '1000GENOMES:phase_3:EUR',
  SAS => '1000GENOMES:phase_3:SAS',
  AA => 'ESP6500:African_American',
  EA => 'ESP6500:European_American',
};

my @population_wide = qw(AA EA AFR AMR EAS EUR SAS ExAC_AFR ExAC_AMR ExAC_Adj ExAC_EAS ExAC_FIN ExAC_NFE ExAC_OTH ExAC_SAS);

sub new {
  my $class = shift;
  
  my $self = $class->SUPER::new(@_);
  my $supported_maf_keys = { map {$_ => 1} qw(minor_allele_freq AFR AMR EAS EUR SAS AA EA ExAC ExAC_AFR ExAC_AMR ExAC_Adj ExAC_EAS ExAC_FIN ExAC_NFE ExAC_OTH ExAC_SAS) };

  my $params = $self->params_to_hash();
  my $file = '';

  # user only supplied file as first param?
  if (!keys %$params) {
    $file = $self->params->[0];
  }
  else {
    $file = $params->{file};

    # process types
    if ($params->{types}) {
      $params->{types} = {map {$_ => 1} split(/[\;\&\|]/, $params->{types})};
    }

    # check maf
    foreach my $maf (qw/maf_monoallelic maf_biallelic/) {
      if($params->{$maf}) {
        die("ERROR: Invalid value for maf: ".$params->{$maf} . "\n") unless
          looks_like_number($params->{$maf}) && ($params->{$maf} >= 0 && $params->{$maf} <= 1)
      }
    }

    my $use_exac = 0;
    my $population_wide = 0;
    my @maf_keys = ();

    if ($params->{maf_key}) {
      foreach my $maf_key (split(',', $params->{maf_key})) {
        die("ERROR: maf_key: " . $maf_key . " not supported. Check plugin documentation for supported maf_keys.\n") unless $supported_maf_keys->{$maf_key};
        $use_exac = 1 if ($maf_key =~ /^ExAC/);
        push @maf_keys, $maf_key;
      }
    }

    if (scalar @maf_keys == 0) {
      $params->{maf_keys} = \@population_wide;  
      $population_wide = 1;
    } else {
      $params->{maf_keys} = \@maf_keys;
    }
    
    if ($use_exac || $population_wide) {
      my $file = $params->{exac_file};
      die("ERROR: ExAC data file is required if you want to filter by ExAC frequencies") unless $file;
      my $exac_plugin = ExAC->new($self->{config}, ($file));
      $self->{config}->{exac_plugin} = $exac_plugin;
    }

  }

  if (!$params->{log_file}) {
    $params->{log_file} = 'g2p_plugin_logs';
  }
  open(my $fh, '>', $params->{log_file}) or die "Could not open file '$params->{log_file}' $!";
  close $fh;

  if (!$params->{log_dir}) {
    $params->{log_dir} = 'g2p_plugin_log_dir';
  }
  my $log_dir = $params->{log_dir};
  if (-d $log_dir) {
    my @files = <$log_dir/*>;
    if (scalar @files > 0) {
      unlink glob "'$log_dir/*.*'";
    }
    @files = <$log_dir/*>;
    die("G2P plugin log dir ($log_dir) is not empty") if (scalar @files > 0);
  } else {
    mkdir $log_dir, 0755;
  }
  # copy in default params
  $params->{$_} //= $DEFAULTS{$_} for keys %DEFAULTS;
  $self->{user_params} = $params;

  my $va = $self->{config}->{reg}->get_adaptor($self->{config}->{species}, 'variation', 'variation');
  $va->db->use_vcf(1);
  $va->db->include_failed_variations(1);
  $self->{config}->{va} = $va;
  my $pa = $self->{config}->{reg}->get_adaptor($self->{config}->{species}, 'variation', 'population');
  $self->{config}->{pa} = $pa;

  # read data from file
  $self->{gene_data} = $self->read_gene_data_from_file($file);
  $self->synonym_mappings();

  # force some config params
  $self->{config}->{individual} //= ['all'];
  $self->{config}->{symbol} = 1;

  if($params->{maf} > 0) {
    $self->{config}->{check_existing} = 1;
    $self->{config}->{failed} = 1;
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
    G2P_gene_req => 'MONO or BI depending on the context in which this gene has been explored',
  };
}

sub run {
  my ($self, $tva, $line) = @_;

  # only interested if we know the zygosity
  my $zyg = $line->{Extra}->{ZYG};
  return {} unless $zyg;

  # only interested in given gene set
  my $tr = $tva->transcript;
  my $gene_symbol = $tr->{_gene_symbol} || $tr->{_gene_hgnc};
  my $gene_data = $self->gene_data($gene_symbol);
  $self->write_report('G2P_in_vcf', $gene_symbol);
  return {} unless $gene_data;

  my @ars = ($gene_data->{'allelic requirement'}) ? @{$gene_data->{'allelic requirement'}} : ();
  return {} unless (@ars && ( grep {$_ =~ /^(mono|bi)allelic$/} @ars ));
 
  # limit by type
  my @consequence_types = map { $_->SO_term } @{$tva->get_all_OverlapConsequences};
  return {} unless grep {$self->{user_params}->{types}->{$_->SO_term}} @{$tva->get_all_OverlapConsequences};

  # limit by MAF
  my $threshold = 0; 
  my $ar_passed = {};
  my ($freqs, $existing_variant);
  foreach my $ar (@ars) {
    if ($ar eq 'monoallelic') {
      $threshold = $self->{user_params}->{maf_monoallelic};
    } else {
      $threshold = $self->{user_params}->{maf_biallelic};
    }
    ($freqs, $existing_variant) = @{$self->get_freq($tva)};
    foreach my $maf_key (keys %$freqs) {
      return {} if ($freqs->{$maf_key} > $threshold);
    }
    $ar_passed->{$ar} = 1;
  }

  my $vf = $tva->base_variation_feature;
  my $allele = $tva->variation_feature_seq;
  my $start = $vf->{start};
  my $end = $vf->{end};

  my $individual = $vf->{individual};
  my $vf_name = $vf->variation_name || $vf->{start}.'_'.$vf->{allele_string};

  my $allele_string = $vf->{allele_string};
  my @alleles = split('/', $allele_string);
  my $ref = $alleles[0]; 
  my $seq_region_name = $vf->{chr};

  my $params = $self->{user_params};
  my $refseq = $tr->{_refseq};
  my $tr_stable_id = $tr->stable_id;
  my $hgvs_t = $tva->hgvs_transcript || 'NA';
  my $hgvs_p = $tva->hgvs_protein || 'NA';
  
  my ($clin_sig, $novel, $failed, $frequencies, $existing_name) = ('NA', 'yes', 'NA', 'NA', 'NA');
  if ($existing_variant) {
    $clin_sig = $existing_variant->{clin_sig} || 'NA';
    $failed = ($existing_variant->{failed}) ? 'yes' : 'no';
    $existing_name = $existing_variant->{variation_name} || 'NA';
    $novel = 'no';
  }
 
  if (scalar keys %$freqs > 0) {
    $frequencies = join(',', map {"$_=$freqs->{$_}"} keys %$freqs);
  }   
 
  my $ar = join(',', sort keys %$ar_passed);
  my $g2p_data = {
    'zyg' => $zyg,
    'allele_requirement' => $ar,
    'frequencies' => $frequencies,
    'consequence_types' => join(',', @consequence_types),
    'refseq' => $refseq,
    'failed' => $failed,
    'clin_sig' => $clin_sig, 
    'novel' => $novel,
    'existing_name' => $existing_name,
    'hgvs_t' => $hgvs_t,
    'hgvs_p' => $hgvs_p,
    'vf_location' => "$seq_region_name:$start-$end $ref/$allele",
  };

  my %return = (
    G2P_flag => $zyg
  );


  $self->write_report('G2P_flag', $gene_symbol, $tr_stable_id, $individual, $vf_name, $g2p_data);

  $self->write_report('G2P_complete', $gene_symbol, $tr_stable_id, $individual, $vf_name, $ar, $zyg);

  my $cache = $self->{cache}->{$individual}->{$tr->stable_id} ||= {};

  delete $cache->{$vf_name} if exists($cache->{$vf_name});

  # biallelic genes require >=1 hom or >=2 hets
  my $gene_reqs = {};
  foreach my $ar (keys %$ar_passed) { 
    if($ar eq 'biallelic') {
      # homozygous, report complete
      if(uc($zyg) eq 'HOM') {
        $return{G2P_complete} = 1;
        $gene_reqs->{BI} = 1;
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
      $gene_reqs->{MONO} = 1;
    }
    else {
      return {};
    }
  }
  if ($return{G2P_complete}) {
    $return{G2P_gene_req} = join(',', sort keys %$gene_reqs);
  }

  return \%return;
}

# read G2P CSV dump
# as from http://www.ebi.ac.uk/gene2phenotype/downloads
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
    $_ =~ s/\R//g;

    my @split = grep defined, "$_," =~ /$re/g;

    unless(@headers) {
      @headers = map {s/\"//g; $_} @split;
    }
    else {
      my %tmp = map {$headers[$_] => $split[$_]} (0..$#split);
      die("ERROR: Gene symbol column not found\n$_\n") unless $tmp{"gene symbol"};
      my $gene_symbol = $tmp{"gene symbol"};
      $gene_data{$gene_symbol}->{"prev symbols"} = $tmp{"prev symbols"};
      push @{$gene_data{$gene_symbol}->{"allelic requirement"}}, $tmp{"allelic requirement"} if ($tmp{"allelic requirement"});
      $self->write_report('G2P_list', $tmp{"gene symbol"}, $tmp{"DDD category"});
    }
  }

  close FILE;

  return \%gene_data;
}

# return either whole gene data hash or one gene's data
# this should allow updates to this plugin to e.g. query a REST server, for example
sub gene_data {
  my ($self, $gene_symbol) = @_;
  my $gene_data = $self->{gene_data}->{$gene_symbol};
  if (!$gene_data) {
    my $prev_gene_symbol = $self->{prev_symbol_mappings}->{$gene_symbol};
    return $prev_gene_symbol ? $self->{gene_data}->{$prev_gene_symbol} : $self->{gene_data};
  } 
  return $gene_data;
}

sub synonym_mappings {
  my $self = shift;
  my $gene_data = $self->{gene_data};
  my $synonym_mappings = {};
  foreach my $gene_symbol (keys %$gene_data) {
    my $prev_symbols = $gene_data->{$gene_symbol}->{'prev symbols'};
    if ($prev_symbols) {
      foreach my $prev_symbol (split(';', $prev_symbols)) {
        $synonym_mappings->{$prev_symbol} = $gene_symbol;
      }
    }
  }
  $self->{prev_symbol_mappings} = $synonym_mappings;
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
  if (!exists($cache->{$allele}->{freq})) {
    foreach my $ex (@{$vf->{existing} || []}) {
      my $existing_allele_string = $ex->{allele_string};
      my $variation_name = $ex->{variation_name};
      next if ($variation_name !~ /^rs/);
      my $freqs = {};  
      my $has_exac = 0;

      foreach my $maf_key (@{$self->{user_params}->{maf_keys}}) {
        my $freq = $self->{user_params}->{default_maf};
        if ($maf_key eq 'minor_allele_freq') {
          if (defined $ex->{minor_allele_freq}) {
            if (($ex->{minor_allele} || '') eq $allele ) {
              $freq = $ex->{minor_allele_freq};
            } else {
              $freq = $self->correct_frequency($tva, $existing_allele_string, $ex->{minor_allele}, $ex->{minor_allele_freq}, $allele, $variation_name, $maf_key) || $freq;
            }
          }
        }
        else {
          foreach my $pair(split(',', $ex->{$maf_key} || '')) {
            my ($a, $f) = split(':', $pair);
            if(($a || '') eq $allele && defined($f)) {
              $freq = $f;
            } else {
              $freq = $self->correct_frequency($tva, $existing_allele_string, $a, $f, $allele, $variation_name, $maf_key) || $freq;
            }
            if ($maf_key =~ /^ExAC/ && $freq) {
              $has_exac = 1;
            }
          }
        }
        $freqs->{$maf_key} = $freq if ($freq);
      }

      if (!$has_exac) {
        my $exac_data = $self->get_ExAC_frequencies($tva); 
        foreach my $maf_key (@{$self->{user_params}->{maf_keys}}) {
          if ($maf_key =~ /^ExAC/) {
            my $exac_key = $maf_key;
            $exac_key =~ s/ExAC/ExAC_AF/;
            my $freq = $exac_data->{$exac_key};
            $freqs->{$maf_key} = $freq if ($freq);
          }
        }      
      }

      $cache->{$allele}->{freq} = $freqs;
      $cache->{$allele}->{ex_variant} = $ex;
    }
  }
  return [$cache->{$allele}->{freq}, $cache->{$allele}->{ex_variant}];
}

sub get_ExAC_frequencies {
  my $self = shift;
  my $tva = shift;
  my $exac_plugin = $self->{config}->{exac_plugin};
  my $exac_data = {};
  eval {
    $exac_data = $exac_plugin->run($tva);
  };
  warn "Problem in ExAC plugin: $@" if $@; 
  return $exac_data;
}

sub correct_frequency {
  my ($self, $tva, $allele_string, $minor_allele, $maf, $allele, $variation_name, $maf_key) = @_;

  if ($maf_key =~ /^ExAC/) {
    $maf_key =~ s/ExAC/ExAC_AF/;
    my $exac_plugin = $self->{config}->{exac_plugin};
    my $exac_data = {};
    eval {
      $exac_data = $exac_plugin->run($tva);
    };
    warn "Problem in ExAC plugin: $variation_name $allele_string $allele $@" if $@; 
    my $freq = $exac_data->{$maf_key};
    return $freq;
  }

  my @existing_alleles = split('/', $allele_string);
  if ($maf_key eq 'minor_allele_freq' && (scalar @existing_alleles == 2)) {
    my $existing_ref_allele = $existing_alleles[0];
    my $existing_alt_allele = $existing_alleles[1];
    if ( ($minor_allele eq $existing_ref_allele && ($allele eq $existing_alt_allele)) || 
         ($minor_allele eq $existing_alt_allele && ($allele eq $existing_ref_allele)) ) {
      return (1.0 - $maf);
    } 
  } else {
    my $va = $self->{config}->{va};
    my $pa = $self->{config}->{pa};
    my $variation = $va->fetch_by_name($variation_name);
    my $maf_key = $self->{user_params}->{maf_key};
    my $population_name = $maf_key_2_population_name->{$maf_key};
    if ($population_name) {
      my $population = $self->{config}->{$population_name};
      if (!$population) {
        $population = $pa->fetch_by_name($population_name);
        $self->{config}->{$population_name} = $population;
      }
      foreach (@{$variation->get_all_Alleles($population)}) {
        if ($_->allele eq $allele) {
          return $_->frequency;
        }
      }
    }
  }     
  return 0.0;
}

sub write_report {
  my $self = shift;
  my $flag = shift;
  my $log_dir = $self->{user_params}->{log_dir};
  my $log_file = "$log_dir/$$.txt";
  open(my $fh, '>>', $log_file) or die "Could not open file '$log_file' $!";
  if ($flag eq 'G2P_list') {
    my ($gene_symbol, $DDD_category) = @_;
    $DDD_category ||= 'Not assigned';
    print $fh "$flag\t$gene_symbol\t$DDD_category\n";
  } elsif ($flag eq 'G2P_in_vcf') {
    my $gene_symbol = shift;
    print $fh "$flag\t$gene_symbol\n";
  } elsif ($flag eq 'G2P_complete') {
    print $fh join("\t", $flag, @_), "\n";
  } else {
    my ($gene_symbol, $tr_stable_id, $individual, $vf_name, $data) = @_;
    $data = join(';', map {"$_=$data->{$_}"} sort keys %$data);
    print $fh join("\t", $flag, $gene_symbol, $tr_stable_id, $individual, $vf_name, $data), "\n";
  }
  close $fh;
}

sub finish {
  my $self = shift;
}

1;
