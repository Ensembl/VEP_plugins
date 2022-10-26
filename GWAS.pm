=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2022] EMBL-European Bioinformatics Institute

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
 ./vep -i variations.vcf --plugin GWAS,type=sstate,file=/FULL_PATH_TO/17463246-GCST000028-EFO_0001360.h.tsv.sorted.gz

=head1 DESCRIPTION

 A VEP plugin that retrieves GWAS catalog data given the file. 
 
 This plugin supports both the curated data that is found in the download section of the GWAS catalog website and the 
 summary statistics file. By default the plugin assumes the file provided is the curated file but you can pass "type=sstate" 
 to say you want to annotate with a summary statistics file

 Please cite the following publication alongside the VEP if you use this resource:
 https://pubmed.ncbi.nlm.nih.gov/30445434/

 Pre-requisites:

 For curated GWAS catalog file -
 GWAS files can be downloaded from - https://www.ebi.ac.uk/gwas/api/search/downloads/alternative
 When run for the first time the plugin will create a processed file that have genomic locations and indexed and put it under 
 the --dir determined by Ensembl VEP. It may take >1 hour to create the processed file depending on the file size. But subsequent 
 runs will be faster as the plugin will use the already generated processed file.
 
 For summary statistics file -
 The plugin can understand the harmonised version of the summary statistics file. Which can be downloaded from the FTP site - 
 http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics
 
 They are under the directory with specific GCST id. For example -
 http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST000001-GCST001000/GCST000028/harmonised/17463246-GCST000028-EFO_0001360.h.tsv.gz
 
 Please keep the filename format same as it is because filename is parsed to get information.
 
 Follow the following steps -
 1. filter out lines that don't have genomic locations -
 tail -n +2 17463246-GCST000028-EFO_0001360.h.tsv | awk '$3 != "NA" && $4 != "NA" {print}' > 17463246-GCST000028-EFO_0001360.h.tsv.filtered

 replace $3 and $4 with column location of the hm_chrom and hm_pos
 
 2. sort the file. Also keep the header and put # in front - 
 (head -n 1 17463246-GCST000028-EFO_0001360.h.tsv && sort -k3,3 -k4,4n 17463246-GCST000028-EFO_0001360.h.tsv.filtered) > 17463246-GCST000028-EFO_0001360.h.tsv.sorted
 sed -i '1 s/^/#/' 17463246-GCST000028-EFO_0001360.h.tsv.sorted

 3. zip the file -
 bgzip -c 17463246-GCST000028-EFO_0001360.h.tsv.sorted > 17463246-GCST000028-EFO_0001360.h.tsv.sorted.gz

 4. index the file -
 tabix -s 3 -b 4 -e 4 -c "#" -f 17463246-GCST000028-EFO_0001360.h.tsv.sorted.gz

 Options are passed to the plugin as key=value pairs:

 file			: (mandatory) Path to GWAS curated or summary statistics file
 type     : type of the file. Valid values are "curated" and "sstate".

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

  die "ERROR: file is not specified which is a mandatory parameter\n" unless defined $param_hash->{file};
  $self->{file} = $param_hash->{file};
  
  $self->{type} = $param_hash->{type} || "curated";
  die "ERROR: provided type ($self->{type}) is not recognized\n" unless(
    $self->{type} eq "curated" || $self->{type} eq "sstate");
  
  # process for cureated file 
  if ($self->{type} eq "curated") {
    my $config = $self->{config};
    
    my $dir = $config->{dir};
    my $input_filename = basename($self->{file});
    $self->{processed_file} = $dir . "/" . $input_filename;
    $self->{processed_file} .= ".gz" unless $self->{processed_file} =~ /gz$/;
    
    # create processed file with genomic location and index - only run if already not created
    unless (-e $self->{processed_file}){
      my $reg = 'Bio::EnsEMBL::Registry';

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
      
      $self->{"va"} = $reg->get_adaptor('human', 'variation', 'variation');
      
      # parse the raw file
      my $data = $self->parse_curated_file($self->{file});
      # create the processed file
      $self->create_processed_file($data);
    }
    
    $self->add_file($self->{"processed_file"});
  }
  # process for summary statistics file
  else {
    $self->add_file($self->{"file"});
    
    # get some information from the filename
    basename($self->{"file"}) =~ m/(.+?)-(.+?)-(.+?)\.(.+)/;
    $self->{"pmid"} = $1;
    $self->{"study"} = $2;
    $self->{"accession"} = $3;
    
    # parse the header to know the desired column location
    $self->{"sstate_colmap"} = $self->parse_sstate_header();
  }
  
  return $self;
}

sub feature_types {
  return ['Feature', 'Intergenic'];
}

sub get_header_info {
  my $self = shift;

  my %header;
  $header{"GWAS_associated_gene"} = "Gene(s) reported by author in GWAS catalog";
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
    $_->{ref} = $vf->ref_allele_string if $_->{ref} eq "";
    my $matches = get_matched_variant_alleles(
      {
        ref    => $vf->ref_allele_string,
        alts   => [$tva->variation_feature_seq],
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

sub get_vfs_from_id {
  my ($self, $id) = @_;
  
  return [] unless defined $self->{va};
  
  my $v = $self->{va}->fetch_by_name($id);
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

sub parse_curated_file {
  my ($self, $input_file) = @_;

  # open the input file for reading
  my $input_FH;
  if($input_file =~ /gz$/) {
    open($input_FH, "zcat " . $input_file . " |") || die ("Could not open $input_file for reading: $!\n");
  }
  else {
    open($input_FH, '<', $input_file) || die ("Could not open $input_file for reading: $!\n");
  }

  my (%headers, @phenotypes);
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

      my $pubmed_id      = $content{'PUBMEDID'};
      my $study          = $content{'STUDY ACCESSION'};
      my $phenotype      = $content{'DISEASE/TRAIT'};
      my $gene           = ($content{'REPORTED GENE(S)'} =~ /\?/) ? '' : $content{'REPORTED GENE(S)'};
      my $rs_risk_allele = ($content{'STRONGEST SNP-RISK ALLELE'} =~ /\?/) ? '' : $content{'STRONGEST SNP-RISK ALLELE'};
      my $rs_id          = $content{'SNPS'};
      my $pvalue         = ($content{'P-VALUE'} ne '') ? $content{'P-VALUE'} : '';
      my $ratio          = $content{'OR OR BETA'};
      my $ratio_info     = $content{'95% CI (TEXT)'};
      my @accessions     = split/\,/, $content{'MAPPED_TRAIT_URI'};

      warn "WARNING: 'DISEASE/TRAIT' entry is empty for '$rs_id'\n" if ($phenotype eq '');
      next if ($phenotype eq '');

      $gene =~ s/\s+//g;
      $gene =~ s/–/-/g;
      $gene =~ s/[^\x00-\x7F]//g; # Remove non ASCII characters
      $gene = '' if $gene eq '-' or $gene eq 'NR'; #Skip uninformative entries, missing data in original curation see GWAS catalog curation

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
      warn "WARNING: Could not parse any rsIds from string '$rs_id'\n" if (!scalar(@ids));
      next if (!scalar(@ids));

      map {
        my $id = $_;

        my $t_data = dclone \%data;
        
        my $vfs = $self->get_vfs_from_id($id);
        $t_data->{"id"} = $id;
        $t_data->{"vfs"} = $vfs;
        
        my $risk_allele;
        map {
          if ($_ =~ /$id/) {
            my $risk_allele_with_id = $_;
            $risk_allele = ( split("-", $risk_allele_with_id) )[1];
          }
        } split(";", $rs_risk_allele);
        $t_data->{"GWAS_risk_allele"} = $risk_allele;

        push(@phenotypes, $t_data);
      } @ids;
    }
  }
  close($input_FH);

  my %result = ('phenotypes' => \@phenotypes);
  return \%result;
}

sub create_processed_file {
  my ($self, $data) = @_;
  
  my $temp_processed_file = $self->{"processed_file"} . "_temp";
  open(my $temp_processed_FH, '>', $temp_processed_file) || die ("Could not open " . $self->{processed_file} . " for writing: $!\n");
  
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
  
  my $header = `tabix $self->{"file"} -H` or die "Cannot get header from " . $self->{"file"} . ": $?\n";
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
