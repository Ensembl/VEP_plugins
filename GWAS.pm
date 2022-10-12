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
 ./vep -i variations.vcf --plugin GWAS,/FULL_PATH_TO/gwas_catalog_v1.0.2-associations_e107_r2022-09-14.tsv

=head1 DESCRIPTION

 A VEP plugin that retrieves GWAS catalog data given the file.

 Please cite the following publication alongside the VEP if you use this resource:
 https://pubmed.ncbi.nlm.nih.gov/30445434/

 Pre-requisites:

 GWAS files can be downloaded from -
 https://www.ebi.ac.uk/gwas/api/search/downloads/alternative

 Options are passed to the plugin as key=value pairs:

 file			: (mandatory) Path to GWAS tsv file

=cut

package GWAS;

use strict;
use warnings;
use Path::Tiny qw(path);
use Storable qw(dclone);

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

sub new {
  my $class = shift;

  my $self = $class->SUPER::new(@_);

  die "ERROR: file is not specified which is a mandatory parameter\n" unless defined $self->params->[0];
  $self->{file} = $self->params->[0];

  $self->{output_file} = $self->{config}->{output_file};

  $self->{data} = $self->parse_input_file($self->{file}, $self->{output_file});

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
  $header{"GWAS_study"} = "Pubmed identifier of the paper published";
  $header{"GWAS_accessions"} = "URI with trait mapped to a ontology term";
  $header{"GWAS_beta_coef"} = "Beta co-efficient reported for the variant";
  $header{"GWAS_odds_ratio"} = "Odds ratio reported for the variant";

  return \%header;
}

sub run {
  my ($self, $tva) = @_;

  my $vf = $tva->variation_feature();
  my $variant_name = $vf->name();

  my $phenotypes = $self->{data}->{phenotypes};

  my $result;
  if ( exists $phenotypes->{$variant_name} ) {
    foreach my $phenotype ( @{ $phenotypes->{$variant_name} } ){
      map {
        my $field = $_;

        if ($_ eq "GWAS_risk_allele") {
          map {
            if ($_ =~ /$variant_name/) {
              my $risk_allele_with_id = $_;
              $phenotype->{$field} = ( split("-", $risk_allele_with_id) )[1];
            }
          } split(";", $phenotype->{$_});
        }

        $result->{$_} = defined $result->{$_} ? $phenotype->{$_} : "," . $phenotype->{$_};
      } keys %{ $phenotype };
    }

    return $result;
  }

  return {};
}


sub parse_input_file {
  my ($self, $input_file, $output_file) = @_;

  my $input_filename = path($input_file)->basename();
  my $workdir = path($output_file)->parent->stringify();

  my ($err_file, $err_FH);
  $err_file = $workdir . '/log_' . $input_filename . ".err";
  open($err_FH, ">", $err_file) || die ("Could not open $err_file for writing: $!\n");

  # Open the input file for reading
  my $input_FH;
  if($input_file =~ /gz$/) {
    open($input_FH, "zcat " . $input_file . " |") || die ("Could not open $input_file for reading: $!\n");
  }
  else {
    open($input_FH, '<', $input_file) || die ("Could not open $input_file for reading: $!\n");
  }

  my (%headers, $phenotypes);
  # Read through the file and parse out the desired fields
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
      my $study          = $content{'STUDY'};
      my $phenotype      = $content{'DISEASE/TRAIT'};
      my $gene           = ($content{'REPORTED GENE(S)'} =~ /\?/) ? '' : $content{'REPORTED GENE(S)'};
      my $rs_risk_allele = ($content{'STRONGEST SNP-RISK ALLELE'} =~ /\?/) ? '' : $content{'STRONGEST SNP-RISK ALLELE'};
      my $rs_id          = $content{'SNPS'};
      my $pvalue         = ($content{'P-VALUE'} ne '') ? $content{'P-VALUE'} : '';
      my $ratio          = $content{'OR OR BETA'};
      my $ratio_info     = $content{'95% CI (TEXT)'};
      my @accessions     = split/\,/, $content{'MAPPED_TRAIT_URI'};

      print $err_FH "WARNING: 'DISEASE/TRAIT' entry is empty for '$rs_id'\n" if ($phenotype eq '');
      next if ($phenotype eq '');

      my $risk_frequency = '';
      if ($rs_risk_allele =~ /^\s*$rs_id-+\s*(\w+)\s*$/i) {
        # $rs_risk_allele = $1;
        $risk_frequency = $content{'RISK ALLELE FREQUENCY'};
      }

      $gene =~ s/\s+//g;
      $gene =~ s/–/-/g;
      $gene =~ s/[^\x00-\x7F]//g; # Remove non ASCII characters
      $gene = '' if $gene eq '-' or $gene eq 'NR'; #Skip uninformative entries, missing data in original curation see GWAS catalog curation

      my %data = (
        'GWAS_associated_gene' => $gene,
        'GWAS_risk_allele' => $rs_risk_allele,
        'GWAS_p_value' => $pvalue,
        'GWAS_accessions'   => \@accessions
      );

      # Post process the ratio data
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

      # Add ratio/coef
      if (defined($ratio)) {
        # Parse the ratio info column to extract the unit information (we are not interested in the confidence interval)
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

      # Parse the ids
      my @ids;
      $rs_id ||= "";
      while ($rs_id =~ m/(rs[0-9]+)/g) {
        push(@ids, $1);
      }

      # Setting pubmed prefix fixed here - would be better to be get from the API
      $data{'GWAS_study'} = "PMID:" . $pubmed_id if (defined($pubmed_id));

      # If we did not get any rsIds, skip this row (this will also get rid of the header)
      print $err_FH "WARNING: Could not parse any rsIds from string '$rs_id'\n" if (!scalar(@ids));
      next if (!scalar(@ids));

      map {
        $phenotypes->{$_} = [] unless defined $phenotypes->{$_};

        my $t_data = dclone \%data;
        $t_data->{"id"} = $_;

        push($phenotypes->{$_}, $t_data);
      } @ids;
    }
  }
  close($input_FH);
  close($err_FH);

  my %result = ('phenotypes' => $phenotypes);
  return \%result;
}

1;
