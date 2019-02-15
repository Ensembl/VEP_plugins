=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2019] EMBL-European Bioinformatics Institute

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

 G2P

=head1 SYNOPSIS

 mv G2P.pm ~/.vep/Plugins
 ./vep -i variations.vcf --plugin G2P,file=/path/to/G2P.csv.gz

=head1 DESCRIPTION

 A VEP plugin that uses G2P allelic requirements to assess variants in genes
 for potential phenotype involvement.

 The plugin has multiple configuration options, though minimally requires only
 the CSV file of G2P data.

 Options are passed to the plugin as key=value pairs, (defaults in parentheses):

 file                  : path to G2P data file, as found at http://www.ebi.ac.uk/gene2phenotype/downloads

 af_monoallelic        : maximum allele frequency for inclusion for monoallelic genes (0.0001)

 af_biallelic          : maximum allele frequency for inclusion for biallelic genes (0.005)
 all_confidence_levels : set value to 1 to include all confidence levels: confirmed, probable and possible. 
                         Default levels are confirmed and probable. 
 af_keys               : reference populations used for annotating variant alleles with observed
                         allele frequencies. Allele frequencies are stored in VEP cache files. 
                         Default populations are:
                         ESP: AA, EA
                         1000 Genomes: AFR, AMR, EAS, EUR, SAS 
                         gnomAD exomes: gnomAD, gnomAD_AFR, gnomAD_AMR, gnomAD_ASJ, gnomAD_EAS, gnomAD_FIN, gnomAD_NFE, gnomAD_OTH, gnomAD_SAS 
                         Separate multiple values with '&'
 af_from_vcf           : set value to 1 to include allele frequencies from VCF file. 
                         Specifiy the list of reference populations to include with --af_from_vcf_keys    
 af_from_vcf_keys      : reference populations used for annotating variant alleles with observed
                         allele frequencies. Allele frequencies are retrieved from VCF files. If
                         af_from_vcf is set to 1 but no populations specified with --af_from_vcf_keys
                         all available reference populations are included. 
                         TOPmed: TOPMed
                         UK10K: ALSPAC, TWINSUK
                         gnomAD exomes: gnomADe:AFR, gnomADe:ALL, gnomADe:AMR, gnomADe:ASJ, gnomADe:EAS, gnomADe:FIN, gnomADe:NFE, gnomADe:OTH, gnomADe:SAS
                         gnomAD genomes: gnomADg:AFR, gnomADg:ALL, gnomADg:AMR, gnomADg:ASJ, gnomADg:EAS, gnomADg:FIN, gnomADg:NFE, gnomADg:OTH
                         Separate multiple values with '&'
 default_af            : default frequency of the input variant if no frequency data is
                         found (0). This determines whether such variants are included;
                         the value of 0 forces variants with no frequency data to be
                         included as this is considered equivalent to having a frequency
                         of 0. Set to 1 (or any value higher than af) to exclude them.
 types                 : SO consequence types to include. Separate multiple values with '&'
                         (splice_donor_variant,splice_acceptor_variant,stop_gained,
                         frameshift_variant,stop_lost,initiator_codon_variant,
                         inframe_insertion,inframe_deletion,missense_variant,
                         coding_sequence_variant,start_lost,transcript_ablation,
                         transcript_amplification,protein_altering_variant)
  
  log_dir              : write stats to log files in log_dir 

  txt_report           : write all G2P complete genes and attributes to txt file

  html_report          : write all G2P complete genes and attributes to html file

 Example:

 --plugin G2P,file=G2P.csv,af_monoallelic=0.05,af_keys=AA&gnomAD_ASJ,types=stop_gained&frameshift_variant
 --plugin G2P,file=G2P.csv,af_monoallelic=0.05,types=stop_gained&frameshift_variant
 --plugin G2P,file=G2P.csv,af_monoallelic=0.05,af_from_vcf=1
 --plugin G2P,file=G2P.csv
 
=cut

package G2P;

use strict;
use warnings;

use Cwd;
use Scalar::Util qw(looks_like_number);
use FileHandle;
use Text::CSV;
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);

use Bio::EnsEMBL::Variation::Utils::BaseVepPlugin;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

our $CAN_USE_HTS_PM;

BEGIN {
  if (eval { require Bio::DB::HTS::Tabix; 1 }) {
    $CAN_USE_HTS_PM = 1;
  }
}


my %DEFAULTS = (

  # vars must have a frequency <= to this to pass
  af => 0.001,
  af_monoallelic => 0.0001,
  af_biallelic => 0.005, 

  af_keys => [qw(AA AFR AMR EA EAS EUR SAS gnomAD gnomAD_AFR gnomAD_AMR gnomAD_ASJ gnomAD_EAS gnomAD_FIN gnomAD_NFE gnomAD_OTH gnomAD_SAS)],

  af_from_vcf_keys => [qw(ALSPAC TOPMed TWINSUK gnomADe:AFR gnomADe:ALL gnomADe:AMR gnomADe:ASJ gnomADe:EAS gnomADe:FIN gnomADe:NFE gnomADe:OTH gnomADe:SAS gnomADg:AFR gnomADg:ALL gnomADg:AMR gnomADg:ASJ gnomADg:EAS gnomADg:FIN gnomADg:NFE gnomADg:OTH)],

  # if no MAF data is found, default to 0
  # this means absence of MAF data is considered equivalent to MAF=0
  # set to 1 to do the "opposite", i.e. exclude variants with no MAF data
  default_af => 0,

  confidence_levels => [qw(confirmed probable)],

  # only include variants with these consequence types
  # currently not ontology-resolved, exact term matches only
  types => {map {$_ => 1} qw(splice_donor_variant splice_acceptor_variant stop_gained frameshift_variant stop_lost initiator_codon_variant inframe_insertion inframe_deletion missense_variant coding_sequence_variant start_lost transcript_ablation transcript_amplification protein_altering_variant)},

);

my $af_key_2_population_name = {
  minor_allele_freq => 'global allele frequency (AF) from 1000 Genomes Phase 3 data',
  AFR => '1000GENOMES:phase_3:AFR',
  AMR => '1000GENOMES:phase_3:AMR',
  EAS => '1000GENOMES:phase_3:EAS',
  EUR => '1000GENOMES:phase_3:EUR',
  SAS => '1000GENOMES:phase_3:SAS',
  AA => 'Exome Sequencing Project 6500:African_American',
  EA => 'Exome Sequencing Project 6500:European_American',
  gnomAD => 'Genome Aggregation Database:Total',
  gnomAD_AFR => 'Genome Aggregation Database exomes:African/African American',
  gnomAD_AMR => 'Genome Aggregation Database exomes:Latino',
  gnomAD_ASJ => 'Genome Aggregation Database exomes:Ashkenazi Jewish',
  gnomAD_EAS => 'Genome Aggregation Database exomes:East Asian',
  gnomAD_FIN => 'Genome Aggregation Database exomes:Finnish',
  gnomAD_NFE => 'Genome Aggregation Database exomes:Non-Finnish European',
  gnomAD_OTH => 'Genome Aggregation Database exomes:Other (population not assigned)',
  gnomAD_SAS => 'Genome Aggregation Database exomes:South Asian',
  ALSPAC => 'UK10K:ALSPAC cohort',
  TOPMed => 'Trans-Omics for Precision Medicine (TOPMed) Program',
  TWINSUK => 'UK10K:TWINSUK cohort',
  'gnomADe:AFR' => 'Genome Aggregation Database exomes v170228',
  'gnomADe:ALL' => 'Genome Aggregation Database exomes v170228',
  'gnomADe:AMR' => 'Genome Aggregation Database exomes v170228',
  'gnomADe:ASJ' => 'Genome Aggregation Database exomes v170228',
  'gnomADe:EAS' => 'Genome Aggregation Database exomes v170228',
  'gnomADe:FIN' => 'Genome Aggregation Database exomes v170228',
  'gnomADe:NFE' => 'Genome Aggregation Database exomes v170228',
  'gnomADe:OTH' => 'Genome Aggregation Database exomes v170228',
  'gnomADe:SAS' => 'Genome Aggregation Database exomes v170228',
  'gnomADg:AFR' => 'Genome Aggregation Database genomes v170228:African/African American',
  'gnomADg:ALL' => 'Genome Aggregation Database genomes v170228:All gnomAD genomes individuals',
  'gnomADg:AMR' => 'Genome Aggregation Database genomes v170228:Latino',
  'gnomADg:ASJ' => 'Genome Aggregation Database genomes v170228:Ashkenazi Jewish',
  'gnomADg:EAS' => 'Genome Aggregation Database genomes v170228:East Asian',
  'gnomADg:FIN' => 'Genome Aggregation Database genomes v170228:Finnish',
  'gnomADg:NFE' => 'Genome Aggregation Database genomes v170228:Non-Finnish European',
  'gnomADg:OTH' => 'Genome Aggregation Database genomes v170228:Other (population not assigned)',
};

my $allelic_requirements = {
  'biallelic' => { af => 0.005, rules => {HET => 2, HOM => 1} },
  'monoallelic' => { af => 0.0001, rules => {HET => 1, HOM => 1} },
  'hemizygous' => { af => 0.0001, rules => {HET => 1, HOM => 1} },
  'x-linked dominant' => { af => 0.0001, rules => {HET => 1, HOM => 1} },
  'x-linked over-dominance' => { af => 0.0001, rules => {HET => 1, HOM => 1} },
};

my @allelic_requirement_terms = keys %$allelic_requirements;

my @population_wide = qw(minor_allele_freq AA AFR ALSPAC AMR EA EAS EUR SAS TOPMed TWINSUK gnomAD gnomAD_AFR gnomAD_AMR gnomAD_ASJ gnomAD_EAS gnomAD_FIN gnomAD_NFE gnomAD_OTH gnomAD_SAS gnomADe:AFR gnomADe:ALL gnomADe:AMR gnomADe:ASJ gnomADe:EAS gnomADe:FIN gnomADe:NFE gnomADe:OTH gnomADe:SAS gnomADg:AFR gnomADg:ALL gnomADg:AMR gnomADg:ASJ gnomADg:EAS gnomADg:FIN gnomADg:NFE gnomADg:OTH);

sub new {
  my $class = shift;

  my $self = $class->SUPER::new(@_);
# suppress warnings that the FeatureAdpators spit if using no_slice_cache
  Bio::EnsEMBL::Utils::Exception::verbose(1999);

  my $supported_af_keys = { map {$_ => 1} @population_wide }; 

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

    # check af
    foreach my $af (qw/af_monoallelic af_biallelic/) {
      if($params->{$af}) {
        die("ERROR: Invalid value for af: ".$params->{$af} . "\n") unless
          looks_like_number($params->{$af}) && ($params->{$af} >= 0 && $params->{$af} <= 1)
      }
    }

    my $assembly =  $self->{config}->{assembly};
    my $af_from_vcf_key_2_collection_id = {
      ALSPAC => {GRCh37 => 'uk10k_GRCh37', GRCh38 => 'uk10k_GRCh38'},
      TOPMed => {GRCh37 => 'topmed_GRCh37', GRCh38 => 'topmed_GRCh38'},
      TWINSUK =>  {GRCh37 => 'uk10k_GRCh37', GRCh38 => 'uk10k_GRCh38'},
      'gnomADe:AFR' => {GRCh37 => 'gnomADe_GRCh37', GRCh38 => 'gnomADe_GRCh38'},
      'gnomADe:ALL' => {GRCh37 => 'gnomADe_GRCh37', GRCh38 => 'gnomADe_GRCh38'},
      'gnomADe:AMR' => {GRCh37 => 'gnomADe_GRCh37', GRCh38 => 'gnomADe_GRCh38'},
      'gnomADe:ASJ' => {GRCh37 => 'gnomADe_GRCh37', GRCh38 => 'gnomADe_GRCh38'},
      'gnomADe:EAS' => {GRCh37 => 'gnomADe_GRCh37', GRCh38 => 'gnomADe_GRCh38'},
      'gnomADe:FIN' => {GRCh37 => 'gnomADe_GRCh37', GRCh38 => 'gnomADe_GRCh38'},
      'gnomADe:NFE' => {GRCh37 => 'gnomADe_GRCh37', GRCh38 => 'gnomADe_GRCh38'},
      'gnomADe:OTH' => {GRCh37 => 'gnomADe_GRCh37', GRCh38 => 'gnomADe_GRCh38'},
      'gnomADe:SAS' => {GRCh37 => 'gnomADe_GRCh37', GRCh38 => 'gnomADe_GRCh38'},
      'gnomADg:AFR' => {GRCh37 => 'gnomADg_GRCh37', GRCh38 => 'gnomADg_GRCh38'},
      'gnomADg:ALL' => {GRCh37 => 'gnomADg_GRCh37', GRCh38 => 'gnomADg_GRCh38'},
      'gnomADg:AMR' => {GRCh37 => 'gnomADg_GRCh37', GRCh38 => 'gnomADg_GRCh38'},
      'gnomADg:ASJ' => {GRCh37 => 'gnomADg_GRCh37', GRCh38 => 'gnomADg_GRCh38'},
      'gnomADg:EAS' => {GRCh37 => 'gnomADg_GRCh37', GRCh38 => 'gnomADg_GRCh38'},
      'gnomADg:FIN' => {GRCh37 => 'gnomADg_GRCh37', GRCh38 => 'gnomADg_GRCh38'},
      'gnomADg:NFE' => {GRCh37 => 'gnomADg_GRCh37', GRCh38 => 'gnomADg_GRCh38'},
      'gnomADg:OTH' => {GRCh37 => 'gnomADg_GRCh37', GRCh38 => 'gnomADg_GRCh38'},
    };

    my @keys = ();
    my $vcf_collection_ids = {};
    if ($params->{af_keys}) {
      push @keys, $params->{af_keys};
    } else {
      push @keys, @{$DEFAULTS{af_keys}};
    }
    if ($params->{af_from_vcf}) {
      if ($params->{af_from_vcf_keys}) {
        push @keys, $params->{af_from_vcf_keys};
      } else {
        push @keys, @{$DEFAULTS{af_from_vcf_keys}};
      }
    }
    
    my @af_keys = ();
    foreach my $af_key_set (@keys) {
      foreach my $af_key (split(/[\;\&\|]/, $af_key_set)) {
        die("ERROR: af_key: " . $af_key . " not supported. Check plugin documentation for supported af_keys.\n") unless $supported_af_keys->{$af_key};
        push @af_keys, $af_key;
        if ($af_from_vcf_key_2_collection_id->{$af_key}) {
          
          $vcf_collection_ids->{$af_from_vcf_key_2_collection_id->{$af_key}->{$assembly}} = 1;
        }
      }
    }
    $params->{af_keys} = \@af_keys;
    $params->{vcf_collection_ids} = $vcf_collection_ids;
  }

  my ($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst) = localtime(time);
  $year += 1900;
  $mon++;
  my $stamp = join('_', ($year, $mon, $mday, $hour, $min));
  my $cwd_dir = getcwd;
  my $new_log_dir = "$cwd_dir/g2p_log_dir\_$stamp";
  my $log_dir = $params->{log_dir} || $new_log_dir;
  if (!-d $log_dir) {
    mkdir $log_dir, 0755;
    $params->{log_dir} = $log_dir;
  } 

  foreach my $report_type (qw/txt_report html_report/) {
    if (!$params->{$report_type}) {
      my $file_type = ($report_type eq 'txt_report') ? 'txt' : 'html';
      $params->{$report_type} = $cwd_dir . "/$report_type\_$stamp.$file_type";
    } 
  }

  if ($params->{all_confidence_levels}) {
    push @{$params->{confidence_levels}}, 'possible', @{$DEFAULTS{confidence_levels}};
  }

  # copy in default params
  $params->{$_} //= $DEFAULTS{$_} for keys %DEFAULTS;
  $self->{user_params} = $params;

  if (!defined($self->{config}->{reg})) {
    my $reg = 'Bio::EnsEMBL::Registry';
    $reg->load_registry_from_db(
      -host       => $self->config->{host},
      -user       => $self->config->{user},
      -port       => $self->config->{port},
      -db_version => $self->config->{db_version},
      -species    => $self->config->{species},
      -no_cache   => $self->config->{no_slice_cache},
    );
    $self->{config}->{reg} = $reg;
  }

  my $va = $self->{config}->{reg}->get_adaptor($self->{config}->{species}, 'variation', 'variation');
  $va->db->use_vcf(1) if ($CAN_USE_HTS_PM);
  $va->db->include_failed_variations(1);
  $self->{config}->{va} = $va;
  my $pa = $self->{config}->{reg}->get_adaptor($self->{config}->{species}, 'variation', 'population');
  $self->{config}->{pa} = $pa;
  my $vca = $self->{config}->{reg}->get_adaptor($self->{config}->{species}, 'variation', 'VCFCollection');
  $self->{config}->{vca} = $vca;
  my $ta = $self->{config}->{reg}->get_adaptor($self->{config}->{species}, 'core', 'transcript');
  $self->{config}->{ta} = $ta;

  # read data from file
  $self->{gene_data} = $self->read_gene_data_from_file($file);
  $self->synonym_mappings();

  # force some config params
  $self->{config}->{individual} //= ['all'];
  $self->{config}->{symbol} = 1;

  $self->{config}->{check_existing} = 1;
  $self->{config}->{failed} = 1;
  $self->{config}->{af} = 1;
  $self->{config}->{af_1kg} = 1;
  $self->{config}->{af_esp} = 1;
  $self->{config}->{af_gnomad} = 1;
#  $self->{config}->{sift} = 'b';
#  $self->{config}->{polyphen} = 'b';
#  $self->{config}->{hgvsc} = 1;
#  $self->{config}->{hgvsp} = 1;

  # tell VEP we have a cache so stuff gets shared/merged between forks
  $self->{has_cache} = 1;
  $self->{cache}->{g2p_in_vcf} = {};

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
  my $zyg = $line->{Extra}->{ZYG} || $line->{ZYG};
  return {} unless $zyg;

  # only interested in given gene set
  my $tr = $tva->transcript;

  my $ensembl_gene_id = $tr->{_gene}->stable_id;
  my $gene_symbol = $tr->{_gene_symbol} || $tr->{_gene_hgnc};
  return {} unless $gene_symbol;
  my $gene_data = $self->gene_data($gene_symbol);
  if (! defined $gene_data) {
    my $ensembl_gene_id = $tr->{_gene}->stable_id;
    $gene_data = $self->gene_data($ensembl_gene_id);
  }

  if (!$self->{cache}->{g2p_in_vcf}->{$gene_symbol}) {
    $self->write_report('G2P_in_vcf', $gene_symbol);
    $self->{cache}->{g2p_in_vcf}->{$gene_symbol} = 1;
  }
  return {} unless defined $gene_data;

  my @ars = ($gene_data->{'allelic requirement'}) ? @{$gene_data->{'allelic requirement'}} : ();
  my %seen;
  @ars = grep { !$seen{$_}++ } @ars;

  return {} unless (@ars && ( grep { exists($allelic_requirements->{$_}) } @ars));
  # limit by type
  my @consequence_types = map { $_->SO_term } @{$tva->get_all_OverlapConsequences};

  return {} unless grep {$self->{user_params}->{types}->{$_->SO_term}} @{$tva->get_all_OverlapConsequences};

  # limit by MAF
  my $threshold = 0; 
  my ($freqs, $existing_variant, $ar_passed) = @{$self->get_freq($tva, \@ars)};

  return {} if (!keys %$ar_passed);

  my $vf = $tva->base_variation_feature;
  my $allele = $tva->variation_feature_seq;
  my $start = $vf->{start};
  my $end = $vf->{end};

  my $individual = $vf->{individual};
  my $vf_name = $vf->variation_name;
  if ($vf_name || $vf_name eq '.') {
    $vf_name = ($vf->{original_chr} || $vf->{chr}) . '_' . $vf->{start} . '_' . ($vf->{allele_string} || $vf->{class_SO_term});
  }
  my $allele_string = $vf->{allele_string};
  my @alleles = split('/', $allele_string);
  my $ref = $alleles[0]; 
  my $seq_region_name = $vf->{chr};

  my $params = $self->{user_params};
  my $refseq = $tr->{_refseq} || 'NA';
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

  my $pph_score   = (defined $tva->polyphen_score) ? $tva->polyphen_score : 'NA';
  my $pph_pred    = (defined $tva->polyphen_prediction) ? $tva->polyphen_prediction : 'NA';
  my $sift_score  = (defined $tva->sift_score) ? $tva->sift_score : 'NA';
  my $sift_pred   = (defined $tva->sift_prediction) ? $tva->sift_prediction : 'NA';
 
  if (scalar keys %$freqs > 0) {
    $frequencies = join(',', map {"$_=$freqs->{$_}"} keys %$freqs);
  }   
 
  my $ar = join(',', sort keys %$ar_passed);
  my $ar_in_g2pdb = join(',', sort @ars);
  my $g2p_data = {
    'zyg' => $zyg,
    'allele_requirement' => $ar,
    'ar_in_g2pdb' => $ar_in_g2pdb,
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
    'sift_score' => "$sift_score",
    'sift_prediction' => $sift_pred,
    'polyphen_score' => "$pph_score",
    'polyphen_prediction' => $pph_pred,
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
    elsif($ar eq 'monoallelic' || $ar eq 'x-linked dominant' || $ar eq 'hemizygous' || $ar eq 'x-linked over-dominance') {
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

  my $assembly =  $self->{config}->{assembly};
  die("ERROR: No file specified or could not read from file ".($file || '')."\n") unless $file && -e $file;

  my @confidence_levels = @{$self->{user_params}->{confidence_levels}};

  # determine file type
  my $file_type;
  my $fh = FileHandle->new($file, 'r');
  while (<$fh>) {
    chomp;
      if (/Model_Of_Inheritance/) {
        $file_type = 'panelapp';
      } elsif (/"allelic requirement"/) {
        $file_type = 'g2p';
      } else {
        $file_type = 'unknown';
      }
      last;
  }
  $fh->close();
  if ($file_type eq 'unknown') {
    if ($file =~ /gz$/) { 
      die("ERROR: G2P plugin can only read uncompressed data");
    } else {
      die("ERROR: Could not recognize input file format. Format must be one of panelapp, g2p or custom. Check website for details: https://www.ebi.ac.uk/gene2phenotype/g2p_vep_plugin");
    }
  }

  if ($file_type eq 'panelapp') {
    my @headers = ();
    my $csv = Text::CSV->new ({ sep_char => "\t" });
    open my $fh, "<:encoding(utf8)", "$file" or die "$file: $!";
    while ( my $row = $csv->getline( $fh ) ) {
      unless (@headers) {
        @headers = @$row;
      } else {
        my %tmp = map {$headers[$_] => $row->[$_]} (0..$#headers);
        my $gene_symbol = $tmp{"Gene Entity Symbol"};
        my $ensembl_gene_id = "";
        if ($assembly eq 'GRCh37') { 
          $ensembl_gene_id = $tmp{"EnsemblId(GRch37)"};
        } else { # GRCh38
          $ensembl_gene_id = $tmp{"EnsemblId(GRch38)"};
        }
        if ($ensembl_gene_id) {
          my @ars = ();
          my $allelic_requirement_panel_app = $tmp{"Model_Of_Inheritance"};
          if ($allelic_requirement_panel_app =~ m/MONOALLELIC|BOTH/) {
            push @ars, 'monoallelic';
          } elsif ($allelic_requirement_panel_app =~ m/BIALLELIC|BOTH/) {
            push @ars, 'biallelic';
          } elsif ($allelic_requirement_panel_app eq 'X-LINKED: hemizygous mutation in males, biallelic mutations in females') {
            push @ars, 'hemizygous';
          } elsif ($allelic_requirement_panel_app eq 'X-LINKED: hemizygous mutation in males, monoallelic mutations in females may cause disease (may be less severe, later onset than males)') {
            push @ars, 'x-linked dominant';
          } else {
            $self->write_report('log', "no allelelic_requirement for $ensembl_gene_id");
          }
          foreach my $ar (@ars) {
            push @{$gene_data{$ensembl_gene_id}->{"allelic requirement"}}, $ar;
          }
        } else {
          $self->write_report('log', "no ensembl gene id for $gene_symbol");
        }
      }
    }
    $csv->eof or $csv->error_diag();
    close $fh;
  }

  if ($file_type eq 'g2p') {
    # this regexp allows for nested ",", e.g.
    # item,description
    # cheese,"salty,delicious"
    my $re = qr/(?: "\( ( [^()""]* ) \)" |  \( ( [^()]* ) \) |  " ( [^"]* ) " |  ( [^,]* ) ) , \s* /x;

    my $fh = FileHandle->new($file, 'r');

    while(<$fh>) {
      chomp;
      $_ =~ s/\R//g;
      my @split = grep defined, "$_," =~ /$re/g;
      unless(@headers) {
        if ($file_type eq 'g2p') {
          @headers = map {s/\"//g; $_} @split;
        } else {
          @headers = @split;
        }
      }
      else {
        my %tmp = map {$headers[$_] => $split[$_]} (0..$#split);
        die("ERROR: Gene symbol column not found\n$_\n") unless $tmp{"gene symbol"};
        my $confidence_value = $tmp{"DDD category"};
        next if (!grep{$_ eq $confidence_value} @confidence_levels);
        my $gene_symbol = $tmp{"gene symbol"};
        $gene_data{$gene_symbol}->{"prev symbols"} = $tmp{"prev symbols"};
        push @{$gene_data{$gene_symbol}->{"allelic requirement"}}, $tmp{"allelic requirement"} if ($tmp{"allelic requirement"});
        $self->write_report('G2P_list', $tmp{"gene symbol"}, $tmp{"DDD category"});
      }
    }
    $fh->close;
  }
  return \%gene_data;
}

# return either whole gene data hash or one gene's data
# this should allow updates to this plugin to e.g. query a REST server, for example
sub gene_data {
  my ($self, $gene_symbol) = @_;
  my $gene_data = $self->{gene_data}->{$gene_symbol};
  if (!$gene_data) {
    my $prev_gene_symbol = $self->{prev_symbol_mappings}->{$gene_symbol};
    return $prev_gene_symbol ? $self->{gene_data}->{$prev_gene_symbol} : undef;
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

sub get_freq {
  my $self = shift;
  my $tva = shift;
  my $ars = shift;
  my $allele = $tva->variation_feature_seq;
  my $vf     = $tva->base_variation_feature;
  reverse_comp(\$allele) if $vf->{strand} < 0;
  my $vf_name = $vf->variation_name;
  if ($vf_name || $vf_name eq '.') {
    $vf_name = ($vf->{original_chr} || $vf->{chr}) . '_' . $vf->{start} . '_' . ($vf->{allele_string} || $vf->{class_SO_term});
  }
  my $cache = $self->{cache}->{$vf_name}->{_g2p_freqs} ||= {};

  if (exists $cache->{$allele}->{failed}) {
    return [$cache->{$allele}->{freq}, $cache->{$allele}->{ex_variant}, {}];
  }

  if (exists $cache->{$allele}->{freq}) {
    return [$cache->{$allele}->{freq}, $cache->{$allele}->{ex_variant}, $cache->{$allele}->{passed_ar}];
  }

  if (!$vf->{existing} || ! scalar @{$vf->{existing}}) {
    my $failed_ars = {};
    my $freqs = {};
    my $passed = $self->frequencies_from_VCF($freqs, $vf, $allele, $ars, $failed_ars);
    if (!$passed) {
      $cache->{$allele}->{failed} = 1;
      return [{}, {}, {}];
    } else {
     $cache->{$allele}->{freq} = $freqs;
     $cache->{$allele}->{ex_variant} = undef;
     # if we get to here return all allelic requirements that passed threshold filtering
     my $passed_ar = {};
     foreach my $ar (@$ars) {
      if (!$failed_ars->{$ar}) {
        $passed_ar->{$ar} = 1;
      }
     }
     $cache->{$allele}->{passed_ar} = $passed_ar;
     return [$cache->{$allele}->{freq}, $cache->{$allele}->{ex_variant}, $cache->{$allele}->{passed_ar}];
    }
  }

  my @existing_variants = @{$vf->{existing}};
  # favour dbSNP variants
  my @dbSNP_variants = grep {$_->{variation_name} =~ /^rs/} @existing_variants;
  if (@dbSNP_variants) {
    @existing_variants = @dbSNP_variants;
  }
  foreach my $ex (@existing_variants) {
    my $existing_allele_string = $ex->{allele_string};
    my $variation_name = $ex->{variation_name};
    my $freqs = {};
    my $failed_ars = {};
    foreach my $af_key (@{$self->{user_params}->{af_keys}}) {
      my $freq = $self->{user_params}->{default_af};
      if ($af_key eq 'minor_allele_freq') {
        if (defined $ex->{minor_allele_freq}) {
          if (($ex->{minor_allele} || '') eq $allele ) {
            $freq = $ex->{minor_allele_freq};
          } else {
            $freq = $self->correct_frequency($tva, $existing_allele_string, $ex->{minor_allele}, $ex->{minor_allele_freq}, $allele, $variation_name, $af_key, $vf_name) || $freq;
          }
        }
      }
      else {
        my @pairs = split(',', $ex->{$af_key} || '');
        my $found = 0;
        if (scalar @pairs == 0) {
          $found = 1; # no allele frequency for this population/af_key available
        }
        foreach my $pair (@pairs) {
          my ($a, $f) = split(':', $pair);
          if(($a || '') eq $allele && defined($f)) {
            $freq = $f;
            $found = 1;
          }
        }
        if (!$found) {
          $freq = $self->correct_frequency($tva, $existing_allele_string, undef, undef, $allele, $variation_name, $af_key, $vf_name) || $freq;
        }
      }
      if (!$self->continue_af_annotation($ars, $failed_ars, $freq)) {
        # cache failed results
        $cache->{$allele}->{failed} = 1;
        return [$cache->{$allele}->{freq}, $cache->{$allele}->{ex_variant}, {}];
      }
      $freqs->{$af_key} = $freq if ($freq);
    }
    if ($self->{user_params}->{af_from_vcf}) {
      my $passed = $self->frequencies_from_VCF($freqs, $vf, $allele, $ars, $failed_ars);
      if (!$passed) {
        $cache->{$allele}->{failed} = 1;
        return [$cache->{$allele}->{freq}, $cache->{$allele}->{ex_variant}, {}];
      }
    }
    $cache->{$allele}->{freq} = $freqs;
    $cache->{$allele}->{ex_variant} = $ex;

    # if we get to here return all allelic requirements that passed threshold filtering
    my $passed_ar = {};
    foreach my $ar (@$ars) {
      if (!$failed_ars->{$ar}) {
        $passed_ar->{$ar} = 1;
      }
    }
    $cache->{$allele}->{passed_ar} = $passed_ar;

  }

  return [$cache->{$allele}->{freq}, $cache->{$allele}->{ex_variant}, $cache->{$allele}->{passed_ar}];
}

sub correct_frequency {
  my ($self, $tva, $allele_string, $minor_allele, $af, $allele, $variation_name, $af_key, $vf_name) = @_;

  my @existing_alleles = split('/', $allele_string);
  if (!grep( /^$allele$/, @existing_alleles)) {
    return 0.0;  
  } 

  if ($af_key eq 'minor_allele_freq' && (scalar @existing_alleles == 2)) {
    my $existing_ref_allele = $existing_alleles[0];
    my $existing_alt_allele = $existing_alleles[1];
    if ( ($minor_allele eq $existing_ref_allele && ($allele eq $existing_alt_allele)) || 
         ($minor_allele eq $existing_alt_allele && ($allele eq $existing_ref_allele)) ) {
      return (1.0 - $af);
    } 
  } else {
    my $va = $self->{config}->{va};
    my $pa = $self->{config}->{pa};
    my $variation = $va->fetch_by_name($variation_name);
    my $af_key = $self->{user_params}->{af_keys};
    my $population_name = $af_key_2_population_name->{$af_key};
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

sub frequencies_from_VCF {
  my $self = shift;
  my $freqs = shift;
  my $vf = shift;
  my $vf_allele = shift;
  my $ars = shift;
  my $failed_ars = shift;
  return 1 if (!defined $self->{user_params}->{af_from_vcf});
  return 1 if (!$CAN_USE_HTS_PM);
  my $vca = $self->{config}->{vca};
  my $collections = $vca->fetch_all;
  foreach my $vc (@$collections) {
    next if (! $self->{user_params}->{vcf_collection_ids}->{$vc->id});
    my $alleles = $vc->get_all_Alleles_by_VariationFeature($vf);
    foreach my $allele (@$alleles) {
      if ($allele->allele eq $vf_allele) {
        my $af_key = $allele->population->name;
        my $freq = $allele->frequency;
        return 0 if (!$self->continue_af_annotation($ars, $failed_ars, $freq));
        $freqs->{$af_key} = $freq;
      }
    }
  }
  return 1;
}

sub continue_af_annotation {
  my $self = shift;
  my $ars = shift;
  my $failed_ars = shift;
  my $freq = shift;
  foreach my $ar (@$ars) {
    if (!$failed_ars->{$ar})  {
      if (defined $allelic_requirements->{$ar}) {
        my $threshold = $allelic_requirements->{$ar}->{af};
        if ($freq > $threshold) {
          $failed_ars->{$ar} = 1;
        }
      }
    } 
  }
  return (scalar @$ars != scalar keys %$failed_ars);
} 



sub write_report {
  my $self = shift;
  my $flag = shift;
  my $log_dir = $self->{user_params}->{log_dir};
  my $log_file = "$log_dir/$$.txt";
  open(my $fh, '>>', $log_file) or die "Could not open file '$flag $log_file' $!";
  if ($flag eq 'G2P_list') {
    my ($gene_symbol, $DDD_category) = @_;
    $DDD_category ||= 'Not assigned';
    print $fh "$flag\t$gene_symbol\t$DDD_category\n";
  } elsif ($flag eq 'G2P_in_vcf') {
    my $gene_symbol = shift;
    print $fh "$flag\t$gene_symbol\n";
  } elsif ($flag eq 'G2P_complete') {
    print $fh join("\t", $flag, @_), "\n";
  } elsif ($flag eq 'log') {
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
  $self->generate_report;
}

sub generate_report {
  my $self = shift;
  my $result_summary = $self->parse_log_files;
  my $chart_txt_data = $self->chart_and_txt_data($result_summary);
  my $chart_data = $chart_txt_data->{chart_data};
  my $txt_data = $chart_txt_data->{txt_data};
  my $canonical_transcripts = $chart_txt_data->{canonical_transcripts};
  $self->write_txt_output($txt_data);
  $self->write_charts($result_summary, $chart_data, $canonical_transcripts);
}

sub write_txt_output {
  my $self = shift;
  my $txt_output_data = shift; 
  my $txt_output_file = $self->{user_params}->{txt_report};
  my $fh_txt = FileHandle->new($txt_output_file, 'w');
  foreach my $individual (sort keys %$txt_output_data) {
    foreach my $gene_symbol (keys %{$txt_output_data->{$individual}}) {
      foreach my $ar (keys %{$txt_output_data->{$individual}->{$gene_symbol}}) {
        foreach my $tr_stable_id (keys %{$txt_output_data->{$individual}->{$gene_symbol}->{$ar}}) {
          my $is_canonical = $txt_output_data->{$individual}->{$gene_symbol}->{$ar}->{$tr_stable_id}->{is_canonical};
          my $canonical_tag = ($is_canonical) ? 'is_canonical' : 'not_canonical';
          my $req =  $txt_output_data->{$individual}->{$gene_symbol}->{$ar}->{$tr_stable_id}->{REQ};
          my $variants = join(';', @{$txt_output_data->{$individual}->{$gene_symbol}->{$ar}->{$tr_stable_id}->{variants}});
          print $fh_txt join("\t", $individual, $gene_symbol, $tr_stable_id, $canonical_tag, "OBS=$ar", "REQ=$req", $variants), "\n";
        }
      }
    }
  }
  $fh_txt->close();
}

sub write_charts {
  my $self = shift;
  my $result_summary = shift;
  my $chart_data = shift;
  my $canonical_transcripts = shift;

  my $count_g2p_genes = keys %{$result_summary->{g2p_list}};
  my $count_in_vcf_file = keys %{$result_summary->{in_vcf_file}};
  my $count_complete_genes = scalar keys %{$result_summary->{complete_genes}};

  my @charts = ();
  my @frequencies_header = (); 

  foreach my $short_name (sort @{$self->{user_params}->{af_keys}}) {
    my $text = $af_key_2_population_name->{$short_name};
    push @frequencies_header, "<a style=\"cursor: pointer\" data-placement=\"top\" data-toggle=\"tooltip\" data-container=\"body\" title=\"$text\">$short_name</a>";
  }

  my $count = 1;
  my @new_header = (
    'Variant location and alleles (REF/ALT)',
    'Variant name', 
    'Existing name', 
    'Zygosity', 
    'All allelic requirements from G2P DB',
    'Consequence types', 
    'ClinVar annotation', 
    'SIFT', 
    'PolyPhen', 
    'Novel variant', 
    'Has been failed by Ensembl', 
    @frequencies_header,
    'HGVS transcript', 
    'HGVS protein', 
    'RefSeq IDs', 
  );

  my $html_output_file = $self->{user_params}->{html_report};
  my $fh_out = FileHandle->new($html_output_file, 'w');
  print $fh_out stats_html_head(\@charts);
  print $fh_out "<div class='main_content container'>";

 
  print $fh_out "<h1>G2P report</h1>";
  print $fh_out "<p>Input and output files:</p>";

  print $fh_out "<dl class='dl-horizontal'>";
  print $fh_out "<dt>G2P list</dt>";
  print $fh_out "<dd>" . $self->{user_params}->{file} .  "</dd>";
  print $fh_out "<dt>Log directory</dt>";
  print $fh_out "<dd>" . $self->{user_params}->{log_dir} .  "</dd>";
  print $fh_out "<dt>HTML report</dt>";
  print $fh_out "<dd>" . $self->{user_params}->{html_report} .  "</dd>";
  print $fh_out "<dt>TXT report</dt>";
  print $fh_out "<dd>" . $self->{user_params}->{txt_report} .  "</dd>";
  print $fh_out "</dl>";

  print $fh_out "<p>Counts:</p>";
  print $fh_out "<dl class='dl-horizontal text-overflow'>";
  print $fh_out "<dt>$count_g2p_genes</dt>";
  print $fh_out "<dd>G2P genes</dd>";
  print $fh_out "<dt>$count_in_vcf_file</dt>";
  print $fh_out "<dd>G2P genes in input VCF file</dd>";
  print $fh_out "<dt>$count_complete_genes</dt>";
  print $fh_out "<dd>G2P complete genes in input VCF file</dd>";
  print $fh_out "</dl>";


  print $fh_out "<h1>Summary of G2P complete genes per individual</h1>";
  print $fh_out "<p>G2P complete gene: A sufficient number of variant hits for the observed allelic requirement in at least one of the gene's transcripts. Variants are filtered by frequency.</p>";
  print $fh_out "<p>Frequency thresholds and number of required variant hits for each allelic requirement:</p>";

  print $fh_out "<table class='table table-bordered'>";
  print $fh_out "<thead>";
  print $fh_out "<tr><th>Allelic requirement</th><th>Frequency threshold for filtering</th><th>Variant counts by zygosity</th></tr>";
  print $fh_out "</thead>";
  print $fh_out "<tbody>";
  foreach my $ar (sort keys %$allelic_requirements) {
    my $af = $allelic_requirements->{$ar}->{af};
    my $rules =  $allelic_requirements->{$ar}->{rules};
    my $rule = join(' OR ', map {"$_ >= $rules->{$_}"} keys %$rules);
    print $fh_out "<tr><td>$ar</td><td>$af</td><td>$rule</td></tr>";
  }
  print $fh_out "</tbody>";
  print $fh_out "</table>";

my $switch =<<SHTML;
<form>
<div class="checkbox">
  <label>
    <input class="target" type="checkbox"> Show only canonical transcript
  </label>
</div>
</form>
SHTML

  print $fh_out $switch;

  foreach my $individual (sort keys %$chart_data) {
    foreach my $gene_symbol (keys %{$chart_data->{$individual}}) {
      foreach my $ar (keys %{$chart_data->{$individual}->{$gene_symbol}}) {
        print $fh_out "<ul>\n";
        foreach my $transcript_stable_id (keys %{$chart_data->{$individual}->{$gene_symbol}->{$ar}}) {
          my $class = ($canonical_transcripts->{$transcript_stable_id}) ? 'is_canonical' : 'not_canonical';
          print $fh_out "<li><a class=\"$class\" href=\"#$individual\_$gene_symbol\_$ar\_$transcript_stable_id\">" . "$individual &gt; $gene_symbol &gt; $ar &gt; $transcript_stable_id" . "</a> </li>\n";
        }
        print $fh_out "</ul>\n";
      }
    }
  }

  foreach my $individual (sort keys %$chart_data) {
    foreach my $gene_symbol (keys %{$chart_data->{$individual}}) {
      foreach my $ar (keys %{$chart_data->{$individual}->{$gene_symbol}}) {
        foreach my $transcript_stable_id (keys %{$chart_data->{$individual}->{$gene_symbol}->{$ar}}) {
          my $class = ($canonical_transcripts->{$transcript_stable_id}) ? 'is_canonical' : 'not_canonical';
          print $fh_out "<div class=\"$class\">";
          my $name = "$individual\_$gene_symbol\_$ar\_$transcript_stable_id";
          my $title = "$individual &gt; $gene_symbol &gt; $ar &gt; $transcript_stable_id";
          print $fh_out "<h3><a name=\"$name\"></a>$title <a title=\"Back to Top\" data-toggle=\"tooltip\" href='#top'><span class=\"glyphicon glyphicon-arrow-up\" aria-hidden=\"true\"></span></a></h3>\n";
          print $fh_out "<div class=\"table-responsive\" style=\"width:100%\">\n";
          print $fh_out "<TABLE  class=\"table table-bordered table-condensed\" style=\"margin-left: 2em\">";
          print $fh_out "<thead>\n";
          print $fh_out "<tr>" . join('', map {"<th>$_</th>"} @new_header) . "</tr>\n";
          print $fh_out "</thead>\n";
          print $fh_out "<tbody>\n";
          foreach my $vf_data (@{$chart_data->{$individual}->{$gene_symbol}->{$ar}->{$transcript_stable_id}}) {
            my $data_row = $vf_data->[0];
            my @tds = ();
            foreach my $cell (@$data_row) {
              my $value = $cell->[0];
              my $class = $cell->[1];
              if ($class) {
                push @tds, "<td class=\"$class\">$value</td>";
              } else {
                push @tds, "<td>$value</td>";
              }
            }
            print $fh_out "<tr>", join('', @tds), "</tr>\n";
          }
          print $fh_out "</tbody>\n";
          print $fh_out "</TABLE>\n";
          print $fh_out "</div>\n";
          print $fh_out "</div>\n";
        }
      }
    }
  }
  print $fh_out stats_html_tail();
}

sub chart_and_txt_data {
  my $self = shift;
  my $result_summary = shift;
  my $individuals = $result_summary->{individuals};
  my $complete_genes = $result_summary->{complete_genes};
  my $acting_ars = $result_summary->{acting_ars};
  my $new_order = $result_summary->{new_order};

#  my @frequencies_header = sort keys $af_key_2_population_name;
  my @frequencies_header = sort @{$self->{user_params}->{af_keys}};

  my $assembly = $self->{config}->{assembly};
  my $transcripts = {};
  my $canonical_transcripts = {};
  my $transcript_adaptor = $self->{config}->{ta};
  my $chart_data = {};
  my $txt_output_data = {};

  my $prediction2bgcolor = {
    'probably damaging' => 'danger',
    'deleterious' => 'danger',
    'possibly damaging' => 'warning',
    'unknown'  => 'warning',
    'benign' => 'success',
    'tolerated' => 'success',
  };

  foreach my $individual (sort keys %$new_order) {
    foreach my $gene_symbol (keys %{$new_order->{$individual}}) {
      foreach my $ar (keys %{$new_order->{$individual}->{$gene_symbol}}) {
        foreach my $transcript_stable_id (keys %{$new_order->{$individual}->{$gene_symbol}->{$ar}}) {
          foreach my $vf_name (keys %{$new_order->{$individual}->{$gene_symbol}->{$ar}->{$transcript_stable_id}}) {
            my $data = $individuals->{$individual}->{$gene_symbol}->{$vf_name}->{$transcript_stable_id};

            my $hash = {};
            foreach my $pair (split/;/, $data) {
              my ($key, $value) = split('=', $pair, 2);
              $value ||= '';
              $hash->{$key} = $value;
            }
            my $vf_location = $hash->{vf_location};
            my $existing_name = $hash->{existing_name};
            if ($existing_name ne 'NA') {
              $existing_name = "<a href=\"http://$assembly.ensembl.org/Homo_sapiens/Variation/Explore?v=$existing_name\">$existing_name</a>";
            }
            my $refseq = $hash->{refseq};
            my $failed = $hash->{failed};
            my $clin_sign = $hash->{clin_sig};
            my $novel = $hash->{novel};
            my $hgvs_t = $hash->{hgvs_t};
            my $hgvs_p = $hash->{hgvs_p};
            my $allelic_requirement = $hash->{allele_requirement};
            my $observed_allelic_requirement = $hash->{ar_in_g2pdb};
            my $consequence_types = $hash->{consequence_types};
            my $zygosity = $hash->{zyg};
            my $sift_score = $hash->{sift_score} || '0.0';
            my $sift_prediction = $hash->{sift_prediction};
            my $sift = 'NA';
            my $sift_class = '';
            if ($sift_prediction ne 'NA') {
              $sift = "$sift_prediction(" . "$sift_score)";
              $sift_class = $prediction2bgcolor->{$sift_prediction};
            }
            my $polyphen_score = $hash->{polyphen_score} || '0.0';
            my $polyphen_prediction = $hash->{polyphen_prediction};
            my $polyphen = 'NA';
            my $polyphen_class = '';
            if ($polyphen_prediction ne 'NA') {
              $polyphen = "$polyphen_prediction($polyphen_score)";
              $polyphen_class =  $prediction2bgcolor->{$polyphen_prediction};
            }
            
            my %frequencies_hash = ();
            if ($hash->{frequencies} ne 'NA') {
              %frequencies_hash = split /[,=]/, $hash->{frequencies};
            }
            my @frequencies = ();
            my @txt_output_frequencies = ();
            foreach my $population (@frequencies_header) {
              my $frequency = $frequencies_hash{$population} || '';
              push @frequencies, ["$frequency"];
              if ($frequency) {
                push @txt_output_frequencies, "$population=$frequency";
              }
            }
            my $is_canonical = 0;
            if ($hash->{is_canonical}) {
              $is_canonical = ($hash->{is_canonical} eq 'yes') ? 1 : 0;
            } else {
              if ($transcripts->{$transcript_stable_id}) {
                $is_canonical = 1 if ($canonical_transcripts->{$transcript_stable_id});
              } else {
                my $transcript = $transcript_adaptor->fetch_by_stable_id($transcript_stable_id);
                if ($transcript) {
                  $is_canonical = $transcript->is_canonical();
                  $transcripts->{$transcript_stable_id} = 1;
                  $canonical_transcripts->{$transcript_stable_id} = 1 if ($is_canonical);
                }
              }
            }
            my ($location, $alleles) = split(' ', $vf_location);
            $location =~ s/\-/:/;
            $alleles =~ s/\//:/;

            push @{$chart_data->{$individual}->{$gene_symbol}->{$ar}->{$transcript_stable_id}}, [[
              [$vf_location], 
              [$vf_name], 
              [$existing_name], 
              [$zygosity], 
              [$observed_allelic_requirement],
              [$consequence_types], 
              [$clin_sign], 
              [$sift, $sift_class], 
              [$polyphen, $polyphen_class], 
              [$novel], 
              [$failed], 
              @frequencies,
              [$hgvs_t], 
              [$hgvs_p], 
              [$refseq] 
            ], $is_canonical];

            my $txt_output_variant = "$location:$alleles:$zygosity:$consequence_types:SIFT=$sift:PolyPhen=$polyphen";
            if (@txt_output_frequencies) {
              $txt_output_variant .= ':' . join(',', @txt_output_frequencies);
            }
            $txt_output_data->{$individual}->{$gene_symbol}->{$ar}->{$transcript_stable_id}->{is_canonical} = $is_canonical;
            $txt_output_data->{$individual}->{$gene_symbol}->{$ar}->{$transcript_stable_id}->{REQ} = $observed_allelic_requirement;
            push @{$txt_output_data->{$individual}->{$gene_symbol}->{$ar}->{$transcript_stable_id}->{variants}}, $txt_output_variant;
          }
        }
      }
    }
  }
  return {txt_data => $txt_output_data, chart_data => $chart_data, canonical_transcripts => $canonical_transcripts};
}

sub parse_log_files {
  my $self = shift;

  my $log_dir = $self->{user_params}->{log_dir}; 
  my @files = <$log_dir/*>;

  my $genes = {};
  my $individuals = {};
  my $complete_genes = {};
  my $g2p_list = {};
  my $in_vcf_file = {};
  my $cache = {};
  my $acting_ars = {};

  my $new_order = {};

  foreach my $file (@files) {
    my $fh = FileHandle->new($file, 'r');
    while (<$fh>) {
      chomp;
      if (/^G2P_list/) {
        my ($flag, $gene_symbol, $DDD_category) = split/\t/;
        $g2p_list->{$gene_symbol} = 1;
      } elsif (/^G2P_in_vcf/) {
        my ($flag, $gene_symbol) = split/\t/;
        $in_vcf_file->{$gene_symbol} = 1;
      } elsif (/^G2P_complete/) {
        my ($flag, $gene_symbol, $tr_stable_id, $individual, $vf_name, $ars, $zyg) = split/\t/;
        foreach my $ar (split(',', $ars)) {
          if ($ar eq 'biallelic') {
            # homozygous, report complete
            if (uc($zyg) eq 'HOM') {
              $complete_genes->{$gene_symbol}->{$individual}->{$tr_stable_id} = 1;
              $acting_ars->{$gene_symbol}->{$individual}->{$ar} = 1;
              $new_order->{$individual}->{$gene_symbol}->{$ar}->{$tr_stable_id}->{$vf_name} = 1;
            }
            # heterozygous
            # we need to cache that we've observed one
            elsif (uc($zyg) eq 'HET') {
              if (scalar keys %{$cache->{$individual}->{$tr_stable_id}} >= 1) {
                $complete_genes->{$gene_symbol}->{$individual}->{$tr_stable_id} = 1;
                $acting_ars->{$gene_symbol}->{$individual}->{$ar} = 1;
                $new_order->{$individual}->{$gene_symbol}->{$ar}->{$tr_stable_id}->{$vf_name} = 1;
                # add first observed het variant to the list
                foreach my $vf (keys %{$cache->{$individual}->{$tr_stable_id}}) {
                  $new_order->{$individual}->{$gene_symbol}->{$ar}->{$tr_stable_id}->{$vf} = 1;
                }
              }
              $cache->{$individual}->{$tr_stable_id}->{$vf_name}++;
            }
          }
          # monoallelic genes require only one allele
          elsif ($ar eq 'monoallelic' || $ar eq 'x-linked dominant' || $ar eq 'hemizygous' || $ar eq 'x-linked over-dominance') {
            $complete_genes->{$gene_symbol}->{$individual}->{$tr_stable_id} = 1;
            $acting_ars->{$gene_symbol}->{$individual}->{$ar} = 1;
            $new_order->{$individual}->{$gene_symbol}->{$ar}->{$tr_stable_id}->{$vf_name} = 1;
          }
        }
      } elsif (/^G2P_flag/) {
        my ($flag, $gene_symbol, $tr_stable_id, $individual, $vf_name, $g2p_data) = split/\t/;
        $genes->{$gene_symbol}->{"$individual\t$vf_name"}->{$tr_stable_id} = $g2p_data;
        $individuals->{$individual}->{$gene_symbol}->{$vf_name}->{$tr_stable_id} = $g2p_data;
      } else {

      }
    }
    $fh->close();
  }
  return {
    genes => $genes,
    individuals => $individuals,
    complete_genes => $complete_genes,
    g2p_list => $g2p_list,
    in_vcf_file => $in_vcf_file,
    acting_ars => $acting_ars,
    new_order => $new_order,
  };
}


sub stats_html_head {
    my $charts = shift;

    my $html =<<SHTML;
<html>
<head>
  <title>VEP summary</title>
  <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css" integrity="sha384-BVYiiSIFeK1dGmJRAkycuHAHRg32OmUcww7on3RYdg4Va+PmSTsz/K68vbdEjh4u" crossorigin="anonymous">
  <style>
    a.inactive {
      color: grey;
      pointer-events:none;
    }
  </style>
</head>
<body>
SHTML
  return $html;
}

sub stats_html_tail {
  my $script =<<SHTML;
  <script src="https://ajax.googleapis.com/ajax/libs/jquery/1.12.4/jquery.min.js"></script>
  <script src="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/js/bootstrap.min.js"></script>
  <script type="text/javascript" src="http://www.google.com/jsapi"></script>
  <script>
    \$( "input[type=checkbox]" ).on( "click", function(){
      if (\$('.target').is(':checked')) {
        \$( "div.not_canonical" ).hide();
        \$("a.not_canonical").addClass("inactive");
      } else {
        \$( "div.not_canonical" ).show();
        \$("a.not_canonical").removeClass("inactive");
      }
    } );
  \$(document).ready(function(){
    \$('[data-toggle="tooltip"]').tooltip(); 
  });
  </script>
SHTML
  return "\n</div>\n$script\n</body>\n</html>\n";
}

sub sort_keys {
  my $data = shift;
  my $sort = shift;
  print $data, "\n";
  my @keys;

  # sort data
  if(defined($sort)) {
    if($sort eq 'chr') {
      @keys = sort {($a !~ /^\d+$/ || $b !~ /^\d+/) ? $a cmp $b : $a <=> $b} keys %{$data};
    }
    elsif($sort eq 'value') {
      @keys = sort {$data->{$a} <=> $data->{$b}} keys %{$data};
    }
    elsif(ref($sort) eq 'HASH') {
      @keys = sort {$sort->{$a} <=> $sort->{$b}} keys %{$data};
    }
  }
  else {
    @keys = keys %{$data};
  }

  return \@keys;
}

1;
