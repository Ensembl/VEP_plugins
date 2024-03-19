=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2023] EMBL-European Bioinformatics Institute

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

 GWAS

=head1 SYNOPSIS

 mv GWAS.pm ~/.vep/Plugins
 ./vep -i variations.vcf --plugin GWAS,file=/FULL_PATH_TO/gwas_catalog_v1.0.2-associations_e107_r2022-09-14.tsv
 ./vep -i variations.vcf --plugin GWAS,type=sstate,file=/FULL_PATH_TO/17463246-GCST000028-EFO_0001360.h.tsv.gz

=head1 DESCRIPTION

 A VEP plugin that retrieves relevant NHGRI-EBI GWAS Catalog data given the file. 
 
 This plugin supports both the curated data that is found in the download section of the NHGRI-EBI GWAS Catalog website and the 
 summary statistics file. By default the plugin assumes the file provided is the curated file but you can pass "type=sstate" 
 to say you want to annotate with a summary statistics file.

 Please cite the following publication alongside the VEP if you use this resource:
 https://pubmed.ncbi.nlm.nih.gov/30445434/

 Pre-requisites:

 For curated NHGRI-EBI GWAS Catalog file -
 GWAS files can be downloaded from - https://www.ebi.ac.uk/gwas/api/search/downloads/alternative
 
 For summary statistics file -
 The plugin can process the harmonised version of the summary statistics file. Which can be downloaded from the FTP site - 
 http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics
 
 They are under directory with related to their specific GCST id. For example -
 http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST000001-GCST001000/GCST000028/harmonised/17463246-GCST000028-EFO_0001360.h.tsv.gz
 
 Please keep the filename format as it is because filename is parsed to get information.
 
 When run for the first time for either type of file, the plugin will create a processed file that have genomic locations and indexed and 
 put it under the --dir location determined by Ensembl VEP. If db=1 option is used, depending on the file size it might take hour(s) to create 
 the processed file. Subsequent runs will be faster as the plugin will be using the already generated processed file. This option is not used by 
 default and the variant information is generally taken directly from the file provided.

 Options are passed to the plugin as key=value pairs:

 file     : (mandatory) Path to GWAS curated or summary statistics file
 type     : type of the file. Valid values are "curated" and "sstate" (summary statistics). Default is "curated".
 verbose  : display info level messages. Valid values are 0 or 1. Default is 0.
 db       : get variant information from Ensembl database during creation of processed file. Valid values are 0 or 1. If Default is 0 (variant
            information is retrieved from curated file)

=cut

package GWAS;

use strict;
use warnings;
use File::Basename;
use Storable qw(dclone);

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Variation::Utils::Sequence qw(get_matched_variant_alleles);

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin);

sub new {
  my $class = shift;

  my $self = $class->SUPER::new(@_);
  
  $self->expand_left(0);
  $self->expand_right(0);
  
  my $param_hash = $self->params_to_hash();

  die "ERROR: please supply a file of NHGRI-EBI GWAS Catalog data using the 'file' parameter to use the GWAS plugin\n" unless defined $param_hash->{file};
  $self->{file} = $param_hash->{file};
  
  $self->{type} = $param_hash->{type} || "curated";
  die "ERROR: provided type ($self->{type}) is not recognized\n" unless(
    $self->{type} eq "curated" || $self->{type} eq "sstate");

  $self->{verbose} = $param_hash->{verbose} || 0;
  $self->{db} = $param_hash->{db} || 0;
  
  # processed file is assumed to be present under --dir 
  my $config = $self->{config};
  my $dir = $config->{dir};
  my $input_filename = basename($self->{file});
  $self->{processed_file} = $dir . "/" . $input_filename;
  $self->{processed_file} .= ".gz" unless $self->{processed_file} =~ /gz$/;
  
  # process for cureated file 
  if ($self->{type} eq "curated") {
    # create processed file with genomic location and index - only run if already not created
    unless (-e $self->{processed_file}){
      # create the processed file from the input file given
      $self->create_curated_processed_file($self->{file});
    }
  }
  # process for summary statistics file
  else {
    # parse the header to know the desired column location
    $self->{"sstate_colmap"} = $self->parse_sstate_header();
    
    # get the chr and pos column numbers in sstate file
    my ($chr, $pos);
    foreach (keys %{ $self->{"sstate_colmap"} }){ 
      $chr = $_ + 1 if $self->{"sstate_colmap"}->{$_} eq "hm_chrom";
      $pos = $_ + 1 if $self->{"sstate_colmap"}->{$_} eq "hm_pos";
    };
    
    # create processed file for sstate - only run if already not created
    unless (-e $self->{processed_file}){
      $self->create_state_processed_file($chr, $pos);
    }
    
    # get some information from the filename
    basename($self->{"file"}) =~ m/(.+?)-(.+?)-(.+?)\.(.+)/;
    $self->{"pmid"} = $1;
    $self->{"study"} = $2;
    $self->{"accession"} = $3;
  }
    
  $self->add_file($self->{"processed_file"});
  
  return $self;
}

sub feature_types {
  return ['Feature', 'Intergenic'];
}

sub get_header_info {
  my $self = shift;

  my %header;
  $header{"GWAS_associated_gene"} = "Gene(s) reported by author in the NHGRI-EBI GWAS Catalog";
  $header{"GWAS_risk_allele"} = "Allele associated with the variant that is strongly associated with the trait";
  $header{"GWAS_p_value"} = "P-value reported for the variant";
  $header{"GWAS_study"} = "GWAS study id";
  $header{"GWAS_pmid"} = "Pubmed identifier of the paper published";
  $header{"GWAS_accessions"} = "URI with trait mapped to a ontology term";
  $header{"GWAS_beta_coef"} = "Beta co-efficient reported for the variant";
  $header{"GWAS_odds_ratio"} = "Odds ratio reported for the variant";

  return \%header;
}

sub run {
  my ($self, $tva) = @_;

  my $vf = $tva->variation_feature();
  return {} unless $vf;
  
  my @data =  @{$self->get_data($vf->{chr}, $vf->{start} - 2, $vf->{end})};
  
  foreach (@data) {
    # if db=0 is used we do not check ref allele 
    $_->{ref} = $vf->ref_allele_string if ($_->{ref} eq "" || $_->{ref} eq "N");

    my $matches = get_matched_variant_alleles(
      {
        ref    => $vf->ref_allele_string,
        alts   => $tva->base_variation_feature->alt_alleles,
        pos    => $vf->{"start"},
        strand => $vf->strand
      },
      {
        ref  => $_->{"ref"},
        alts => [$_->{"result"}->{"GWAS_risk_allele"}],
        pos  => $_->{"start"},
      }
    );
    return $_->{"result"} if (@$matches);
  }
  
  return {};
}

sub get_vfs_from_db {
  my ($self, $id) = @_;

  my $reg = 'Bio::EnsEMBL::Registry';
  my $config = $self->{config};

  if($config->{host}) {
      $reg->load_registry_from_db(
          -host       => $config->{host},
          -user       => $config->{user},
          -pass       => $config->{password},
          -port       => $config->{port},
          -species    => $config->{species},
          -db_version => $config->{db_version},
          -no_cache   => $config->{no_slice_cache},
      );
  }

  my $va = $reg->get_adaptor($config->{species}, 'variation', 'variation');
  return [] unless defined $va;
  
  my $v = $va->fetch_by_name($id);
  return [] unless defined $v;
  
  my $locations = [];
  foreach my $vf (@{ $v->get_all_VariationFeatures() }) {
    my $location = {
      "seq"   => $vf->seq_region_name(),
      "start" => $vf->seq_region_start(),
      "end"   => $vf->seq_region_end(),
      "ref"   => $vf->ref_allele_string()
    };
    
    push @{ $locations }, $location;
  }
  
  return $locations;
}

sub get_vfs_from_file {
  my ($self, $chr, $start, $end) = @_;

  my @chrs = map { local $_ = $_; s/\s+//g; $_ } split(/[;,x]/, $chr);
  my @starts = map { local $_ = $_; s/\s+//g; $_ } split(/[;,x]/, $start);
  my @ends = map { local $_ = $_; s/\s+//g; $_ } split(/[;,x]/, $end);

  return [] if (scalar @chrs != scalar @starts && scalar @chrs != scalar @ends);

  my $locations = [];
  for (0..$#chrs) {
    my $location = {
      "seq"   => $chrs[$_],
      "start" => $starts[$_],
      "end"   => $ends[$_],
      "ref"   => "N"
    };

    push @{ $locations }, $location;
  }
  
  return $locations;
}

sub parse_curated_file {
  my ($self, $content) = @_;

  my $pubmed_id      = $content->{'PUBMEDID'};
  my $study          = $content->{'STUDY ACCESSION'};
  my $phenotype      = $content->{'DISEASE/TRAIT'};
  my $gene           = ($content->{'REPORTED GENE(S)'} =~ /\?/) ? '' : $content->{'REPORTED GENE(S)'};
  my $rs_risk_allele = ($content->{'STRONGEST SNP-RISK ALLELE'} =~ /\?/) ? '' : $content->{'STRONGEST SNP-RISK ALLELE'};
  my $rs_id          = $content->{'SNPS'};
  my $pvalue         = ($content->{'P-VALUE'} ne '') ? $content->{'P-VALUE'} : '';
  my $ratio          = $content->{'OR OR BETA'};
  my $ratio_info     = $content->{'95% CI (TEXT)'};
  my @accessions     = split/\,/, $content->{'MAPPED_TRAIT_URI'};
  my $chr            = $content->{'CHR_ID'};
  my $start          = $content->{'CHR_POS'};;
  my $end            = $content->{'CHR_POS'};; 

  warn "INFO: 'DISEASE/TRAIT' entry is empty for '$rs_id'\n" if (($phenotype eq '') && $self->{verbose});
  return {} if ($phenotype eq '');

  $gene =~ s/\s+//g;
  $gene =~ s/–/-/g;
  $gene =~ s/[^\x00-\x7F]//g; # Remove non ASCII characters
  $gene = '' if $gene eq '-' or $gene eq 'NR'; #Skip uninformative entries, missing data in original curation see the NHGRI-EBI GWAS Catalog curation

  my %data = (
    'GWAS_associated_gene' => $gene,
    'GWAS_p_value' => $pvalue,
    'GWAS_accessions'   => \@accessions,
    'GWAS_pmid' => $pubmed_id,
    'GWAS_study' => $study 
  );

  # post process the ratio data
  if (defined($ratio)) {
    if ($ratio =~ /(\d+)?(\.\d+)$/) {
      my $pre  = $1;
      my $post = $2;
      $ratio = (defined($pre)) ? "$pre$post" : "0$post";
      $ratio = 0 if ($ratio eq '0.00');
    } else {
      $ratio = undef;
    }
  }

  # add ratio/coef
  if (defined($ratio)) {
    # parse the ratio info column to extract the unit information (we are not interested in the confidence interval)
    if ($ratio_info =~ /^\s*(\[.+\])?\s*(.+)$/) {
      my $unit = $2;
          $unit =~ s/\(//g;
          $unit =~ s/\)//g;
          $unit =~ s/µ/micro/g;
      if ($unit =~ /decrease|increase/) {
        $data{'GWAS_beta_coef'} = "$ratio $unit";
      }
      else {
        $data{'GWAS_odds_ratio'} = $ratio;
      }
    }
    else {
      $data{'GWAS_odds_ratio'} = $ratio;
    }
  }
  $data{'GWAS_odds_ratio'} = "" unless defined $data{'GWAS_odds_ratio'};
  $data{'GWAS_beta_coef'} = "" unless defined $data{'GWAS_beta_coef'};

  # parse the ids
  my @ids;
  $rs_id ||= "";
  while ($rs_id =~ m/(rs[0-9]+)/g) {
    push(@ids, $1);
  }

  # if we did not get any rsIds, skip this row (this will also get rid of the header)
  warn "INFO: Could not parse any rsIds from string '$rs_id'\n" if (!scalar(@ids) && $self->{verbose});
  return {} if (!scalar(@ids));

  my @phenotypes;
  map {
    my $id = $_;

    my $t_data = dclone \%data;
    
    my $vfs = $self->{db} ? $self->get_vfs_from_db($id) : $self->get_vfs_from_file($chr, $start, $end);
    $t_data->{"id"} = $id;
    $t_data->{"vfs"} = $vfs;
    
    my $risk_allele;
    map {
      if ($_ =~ /$id/) {
        my $risk_allele_with_id = $_;
        $risk_allele = ( split("-", $risk_allele_with_id) )[1];
      }
    } split(/[;,x]/, $rs_risk_allele);
    $t_data->{"GWAS_risk_allele"} = $risk_allele;

    push(@phenotypes, $t_data);
  } @ids;

  my %result = ('phenotypes' => \@phenotypes);
  return \%result;
}

sub create_curated_processed_file {
  my ($self, $input_file) = @_;

  my $temp_processed_file = $self->{"processed_file"} . "_temp";
  open(my $temp_processed_FH, '>', $temp_processed_file) || die ("Could not open " . $self->{processed_file} . " for writing: $!\n");

  # open the input file for reading
  my $input_FH;
  if($input_file =~ /gz$/) {
    open($input_FH, "zcat " . $input_file . " |") || die ("Could not open $input_file for reading: $!\n");
  }
  else {
    open($input_FH, '<', $input_file) || die ("Could not open $input_file for reading: $!\n");
  }

  my %headers;
  # read through the file and parse out the desired fields
  while (<$input_FH>) {
    chomp;

    my @row_data = split(/\t/,$_);

    # header
    if(/^DATE\s+ADDED\s+TO\s+CATALOG/) {
      $headers{uc($row_data[$_])} = $_ for 0..$#row_data;
    }
    else {
      die ("ERROR: Could not find header data\n") unless %headers;

      my %content;
      $content{$_} = $row_data[$headers{$_}] for keys %headers;
      my $data = $self->parse_curated_file(\%content);
  
      foreach my $phenotype (@{ $data->{"phenotypes"} }){
        foreach my $vf (@{ $phenotype->{"vfs"} }){
          next unless $vf;
          
          next unless $phenotype->{"GWAS_risk_allele"};
          my $GWAS_risk_allele = $phenotype->{"GWAS_risk_allele"} || "";
          
          my $GWAS_associated_gene = $phenotype->{"GWAS_associated_gene"} || "";
          my $GWAS_p_value = $phenotype->{"GWAS_p_value"} || "";
          my $GWAS_study = $phenotype->{"GWAS_study"} || "";
          my $GWAS_pmid = $phenotype->{"GWAS_pmid"} || "";
          my $accessions = join(",", @{ $phenotype->{"GWAS_accessions"} }) || "";
          my $GWAS_beta_coef = $phenotype->{"GWAS_beta_coef"} || "";
          my $GWAS_odds_ratio = $phenotype->{"GWAS_odds_ratio"} || "";

          my $line = join("\t", (
            $vf->{"seq"}, $vf->{"start"}, $vf->{"end"}, $vf->{"ref"}, 
            $GWAS_associated_gene,
            $GWAS_risk_allele,
            $GWAS_p_value,
            $GWAS_study,
            $GWAS_pmid,
            $accessions,
            $GWAS_beta_coef,
            $GWAS_odds_ratio,
          ));
          
          print $temp_processed_FH $line . "\n";
        }
      }
    }
  }
  
  close($temp_processed_FH);
  
  system("sort -k1,1 -k2,2n -k3,3n $temp_processed_file | bgzip -c > $self->{processed_file}") == 0
    or die "Failed to sort and compress $temp_processed_file: $?\n";
  system("rm $temp_processed_file") == 0
    or die "Failed to delete $temp_processed_file: $?\n";
  system("tabix -s 1 -b 2 -e 3 -f " . $self->{processed_file}) == 0
    or die "Failed to create index " . $self->{processed_file} . ": $?\n";;
}

sub parse_sstate_header {
  my ($self) = @_;
  
  my $header = `head -n 1 $self->{"file"}` or die "Cannot get header from " . $self->{"file"} . ": $?\n";
  $header =~ s/^#//;
  
  my @cols = (split("\t", $header));
  my @required_cols = qw/hm_chrom hm_pos p_value hm_beta hm_odds_ratio hm_effect_allele/;
  
  my $colmap = {};
  map {
    my $index = $_;
    my @matched = (grep {/^$cols[$index]$/} @required_cols);
    $colmap->{$index} = $matched[0] if @matched;
  } 0 .. $#cols;

  return $colmap;
}

sub create_state_processed_file {
  my ($self, $c, $p) = @_;
  
  my $file = $self->{"file"};
  my $processed_file = $self->{"processed_file"};
  
  my $filtered_file = $file . ".filtered";
  system("gunzip -cf $file | tail -n +2 $file | awk '\$$c != \"NA\" && \$$p != \"NA\" {print}' > $filtered_file") == 0
    or die "Failed to filter $file: $?\n";
  
  my $sorted_file = $file . ".sorted";
  system("sort -k${c},${c} -k${p},${p}n $filtered_file > $sorted_file") == 0
    or die "Failed to sort $filtered_file: $?\n";

  my $zipped_file = $sorted_file . ".gz";
  system("bgzip -c $sorted_file > $zipped_file") == 0
    or die "Failed to zip $sorted_file: $?\n";

  system("tabix -s $c -b $p -e $p -f $zipped_file") == 0
    or die "Failed to create index for $zipped_file: $?\n";
    
    
  system("mv $zipped_file $processed_file") == 0
    or die "Failed to rename files $zipped_file to $processed_file: $?\n";
    
  system("mv $zipped_file.tbi $processed_file.tbi") == 0
    or die "Failed to rename files $zipped_file.tbi to $processed_file: $?\n";
}

sub parse_data {
  my ($self, $line) = @_;
  my ($c, $s, $e, $ref, $a_gene, $r_allele, $p_val, $study, $pmid, $acc, $beta, $odds);
  
  if ($self->{"type"} eq "curated"){
    ($c, $s, $e, $ref, $a_gene, $r_allele, $p_val, $study, $pmid, $acc, $beta, $odds) = split /\t/, $line;
  }
  else {
    $study = $self->{"study"};
    $pmid = $self->{"pmid"};
    $acc = $self->{"accession"};
    
    my @cols = split /\t/, $line;
    my $parsed_data = {};
    map {
      if ($self->{"sstate_colmap"}->{$_}){
        my $col_name = $self->{"sstate_colmap"}->{$_};
        $parsed_data->{$col_name} = $cols[$_] ;
      }
    } 0 .. $#cols;
    
    $c = $parsed_data->{"hm_chrom"};
    $s = $parsed_data->{"hm_pos"};
    $e = $s;
    $r_allele = $parsed_data->{"hm_effect_allele"} || "";
    $p_val = $parsed_data->{"p_value"} || "";
    $beta = $parsed_data->{"hm_beta"} || "";
    $odds = $parsed_data->{"hm_odds_ratio"} || "";
    
    $a_gene = "";
    $ref = "";
  }
  
  return {
    chr => $c,
    start => $s,
    end => $e,
    ref => $ref,
    result => {
      GWAS_associated_gene => $a_gene,
      GWAS_risk_allele => $r_allele,
      GWAS_p_value => $p_val,
      GWAS_study => $study,
      GWAS_pmid => $pmid,
      GWAS_accessions => $acc,
      GWAS_beta_coef => $beta,
      GWAS_odds_ratio => $odds
    }
  };
}

1;
