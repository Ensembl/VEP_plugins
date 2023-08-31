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

 TrioAnalysis

=head1 SYNOPSIS

 mv TrioAnalysis.pm ~/.vep/Plugins
 ./vep -i variations.vcf --plugin TrioAnalysis,ped=samples.ped

=head1 DESCRIPTION

 A VEP plugin that identifies de novo variants in a VCF file.
 
 Options are passed to the plugin as key=value pairs:

 ped                : Path to PED file (mandatory)
                      The file is tab or white-space delimited with six mandatory columns:
                        - family ID
                        - individual ID
                        - paternal ID
                        - maternal ID
                        - sex
                        - phenotype

 report_dir         : write files in report_dir (optional)


 The plugin can then be run:
 ./vep -i variations.vcf --plugin TrioAnalysis,ped=samples.ped
 ./vep -i variations.vcf --plugin TrioAnalysis,ped=samples.ped,report_dir=path/to/dir


=cut

package TrioAnalysis;

use strict;
use warnings;
use Cwd;

use Bio::EnsEMBL::Variation::VariationFeature;
use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepFilterPlugin);


sub _parse_ped_file {
  my $self = shift;
  my $file = shift;

  open FILE, $file;
  while(<FILE>) {
    chomp;
    my ($family_id, $ind_id, $paternal_id, $maternal_id, $sex, $pheno) = split/\s+|\t/;
    
    if(!defined ($family_id && $ind_id && $paternal_id && $maternal_id && $sex && $pheno)) {
      die "ERROR: PED file requires dix columns: family ID, individual ID, paternal ID, maternal ID, sex, phenotype\n";
    }
    
    $self->{linkage}->{$ind_id}->{paternal_id} = $paternal_id;
    $self->{linkage}->{$ind_id}->{maternal_id} = $maternal_id;
    $self->{linkage}->{$ind_id}->{sex} = $sex;
    $self->{linkage}->{$ind_id}->{pheno} = $pheno;
    $self->{linkage}->{$ind_id}->{child} = 1 if ($pheno == 2);
  }
  close FILE;
}

sub new {
  my $class = shift;

  my $self = $class->SUPER::new(@_);

  my $param_hash = $self->params_to_hash();

  if (defined $param_hash->{ped} && -e $param_hash->{ped}) {
    $self->_parse_ped_file($param_hash->{ped});
  }
  else {
    die "ERROR: Please provide a valid PED file\n";
  }

  if (defined $param_hash->{report_dir}) {
    if (!-d $param_hash->{report_dir}) {
      my $return = mkdir $param_hash->{report_dir}, 0755;
      die("ERROR: Couldn't create report_dir ", $param_hash->{report_dir}, " $!\n") if (!$return);
    }
    $self->{report_dir} = $param_hash->{report_dir};
  }
  else {
    my $cwd_dir = getcwd;
    $self->{report_dir} = "$cwd_dir";
  }

  # force some config params
  $self->{config}->{individual_zyg} = ['all'];

  # Report files
  $self->{report_de_novo} = "variants_de_novo.txt";
  $self->{report_all} = "variants_child_and_both_parents.txt";
  # Write header
  write_header($self->{report_dir}.'/'.$self->{report_de_novo}, 1);
  write_header($self->{report_dir}.'/'.$self->{report_all});

  return $self;
}

sub feature_types {
  return ['Feature', 'Intergenic'];
}

sub get_header_info {
  my $self = shift;

  my %header;

  $header{'TrioAnalysis'} = '';

  return \%header;
}

sub run {
  my ($self, $tva, $line) = @_;
  
  my $vf = $tva->variation_feature;
  my $zyg = defined($line->{Extra}) ? $line->{Extra}->{ZYG} : $line->{ZYG};
  
  my $chr = $vf->{chr};
  my $start = $vf->{start};

  my $result;
  my $list_of_ind;

  foreach my $geno_ind (@{$zyg}) {
    my ($ind, $geno) = split(':', $geno_ind);

    # Check if VCF and PED file have the same individual IDs
    if(!$self->{linkage}->{$ind}) {
      die "ERROR: VCF and PED individuals do not match. Please check the individual IDs in PED file: ", join(',', keys %{$self->{linkage}}), "\n";
    }

    # HOM : homozygous
    # HET : heterozygous
    # HOMREF : homozygous reference (not a variant)
    if($geno eq 'HOM' || $geno eq 'HET') {
      if($self->{linkage}->{$ind} && $self->{linkage}->{$ind}->{child}) {
        $list_of_ind->{'child'} = $geno_ind;
      }
      elsif($self->{linkage}->{$ind}) {
        push @{$list_of_ind->{'parent'}}, $ind;
      }
    }
  }

  if(scalar(keys %{$list_of_ind}) == 1) {
    if(defined $list_of_ind->{'child'}) {
      $result = 'only_in_child';
      write_report($self->{report_dir}.'/'.$self->{report_de_novo}, $line, $tva, $list_of_ind->{'child'});
    }
    elsif(scalar(@{$list_of_ind->{'parent'}}) == 1) {
      $result = 'only_in_one_parent';
    }
    else {
      $result = 'only_in_both_parents';
    }
  }
  elsif(scalar(keys %{$list_of_ind}) == 2) {
    if(scalar(@{$list_of_ind->{'parent'}}) == 2) {
      $result = 'in_child_and_both_parents';
      write_report($self->{report_dir}.'/'.$self->{report_all} , $line, $tva);
    }
    else {
      $result = 'in_child_and_one_parent';
    }
  }
  else {
    $result = 'not_found';
  }

  return $result ? { TrioAnalysis => $result } : {};
}

sub write_header {
  my $file = shift;
  my $flag = shift;
  
  open(my $fh, '>', $file) or die "Could not open file $file $!\n";
  
  my $line = "Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tProtein_position
            Amino_acids\tCodons";

  if($flag) {
    $line .= "\tde_novo_sample\tde_novo_zygosity";
  }

  print $fh $line . "\n";

  close($fh);
}

sub write_report {
  my $file = shift;
  my $line = shift;
  my $tva = shift;
  my $geno_ind = shift;

  my $cdna = $tva->transcript_variation->cdna_start() ? $tva->transcript_variation->cdna_start() : '-';
  my $peptide_start = defined($tva->transcript_variation->translation_start) ? $tva->transcript_variation->translation_start : '-';
  my $aa_string = $tva->pep_allele_string ? $tva->pep_allele_string : '-';
  my $codon = $tva->transcript_variation->codons ? $tva->transcript_variation->codons : '-';
  my $existing_variation;
  
  my $ind;
  my $geno;
  
  open(my $fh, '>>', $file) or die "Could not open file $file $!\n";

  my $out_line = $line->{Uploaded_variation}."\t".$line->{Location}."\t".$line->{Allele}."\t".$line->{Gene}."\t"
  .$line->{Feature}."\t".$line->{Feature_type}."\t".join(',', @{$line->{Consequence}})."\t".$cdna."\t".$peptide_start.
  "\t".$aa_string."\t".$codon;

  if($geno_ind) {
    ($ind, $geno) = split(':', $geno_ind);
    $out_line .= "\tde_novo_sample=".$ind."\tde_novo_zyg=".$geno;
  }

  print $fh $out_line."\n";

  close($fh);
}

1;
