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

 DeNovo

=head1 SYNOPSIS

 mv DeNovo.pm ~/.vep/Plugins
 ./vep -i variations.vcf --plugin DeNovo,ped=samples.ped
  ./vep -i variations.vcf --plugin DeNovo,ped=samples.ped,full_report=1

=head1 DESCRIPTION

 A VEP plugin that identifies de novo variants in a VCF file.
 The plugin is not compatible with JSON output format.
 
 Options are passed to the plugin as key=value pairs:

 ped                : Path to PED file (mandatory)
                      The file is tab or white-space delimited with five mandatory columns:
                        - family ID
                        - individual ID
                        - paternal ID
                        - maternal ID
                        - sex
                        - phenotype (optional)

 report_dir         : Write files in report_dir (optional)

 full_report        : Set to 1 to report all types of variants (optional)
                      By default, the plugin only reports de novo variants.


 The plugin can then be run:
 ./vep -i variations.vcf --plugin DeNovo,ped=samples.ped
 ./vep -i variations.vcf --plugin DeNovo,ped=samples.ped,report_dir=path/to/dir
 ./vep -i variations.vcf --plugin DeNovo,ped=samples.ped,report_dir=path/to/dir,full_report=1


=cut

package DeNovo;

use strict;
use warnings;
use Cwd;
use List::MoreUtils qw(uniq);

use Bio::EnsEMBL::Variation::VariationFeature;
use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepFilterPlugin);


sub _parse_ped_file {
  my $self = shift;
  my $file = shift;

  open(my $fh, '<', $file) or die "Could not open file $file $!\n";

  while(<$fh>) {
    chomp;
    my ($family_id, $ind_id, $paternal_id, $maternal_id, $sex, $pheno) = split/\s+|\t/;

    if(!defined ($family_id && $ind_id && $paternal_id && $maternal_id && $sex)) {
      die "ERROR: PED file requires five columns: family ID, individual ID, paternal ID, maternal ID, sex\n";
    }

    $self->{linkage}->{$ind_id}->{family_id} = $family_id;
    $self->{linkage}->{$ind_id}->{sex} = $sex;
    $self->{linkage}->{$ind_id}->{pheno} = $pheno if(defined $pheno);
    $self->{linkage}->{$ind_id}->{parents}->{mother} = $maternal_id if ($maternal_id);
    $self->{linkage}->{$ind_id}->{parents}->{father} = $paternal_id if ($paternal_id);

    $self->{family_list}->{$family_id} = 1;
  }
  close $fh;
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

    # Check if directory is empty
    opendir(DIR, $self->{report_dir}) or die "Cannot open directory ", $self->{report_dir}, "$!";
    my @files = readdir(DIR);
    foreach my $file (@files) {
      if($file =~ /.txt$/) {
        die("ERROR: Report directory is not empty please empty directory before running the plugin: ", $self->{report_dir}, "\n");
      }
    }
  }
  else {
    my $cwd_dir = getcwd;
    $self->{report_dir} = "$cwd_dir";
  }

  if (defined $param_hash->{full_report}) {
    $self->{full_report} = 1;
  }

  # force some config params
  $self->{config}->{individual_zyg} = ['all'];

  # Report files for each family
  my $family_list = $self->{family_list};
  foreach my $family (keys %{$family_list}) {
    $self->{report_de_novo}->{$family} = $family . "_variants_de_novo.txt";
    # Write header
    write_header($self->{report_dir}.'/'.$self->{report_de_novo}->{$family}, 1);

    if($self->{full_report}) {
      $self->{report_both_parents}->{$family} = $family . "_variants_both_parents.txt";
      $self->{report_all}->{$family} = $family . "_variants_child_and_both_parents.txt";
      # Write header
      write_header($self->{report_dir}.'/'.$self->{report_both_parents}->{$family});
      write_header($self->{report_dir}.'/'.$self->{report_all}->{$family});
    }
  }

  return $self;
}

sub feature_types {
  return ['Feature', 'Intergenic'];
}

sub get_header_info {
  my $self = shift;

  my %header;
  my $header_text;

  if ($self->{full_report}) {
    $header_text = 'De novo variants identified. The output includes the following flags: ' .
                    'de_novo_alt - variant only found in the proband, ' .
                    'only_in_one_parent - variant found in one parent, ' .
                    'only_in_both_parents - variant found in both parents, ' .
                    'in_child_and_one_parent - variant found in proband and one parent, ' .
                    'in_child_hom_and_one_parent_het - variant found in proband (homozygous) and one parent (heterozygous), ' .
                    'in_child_and_both_parents - variant found in proband and both parents, ' .
                    'in_child_het_and_both_parents_hom - variant found in proband (heterozygous) and both parents (homozygous),' .
                    'not_found - variant not found in any sample';
  }
  else {
    $header_text = 'De novo variants identified. The output includes the following flags: ' .
                   'de_novo_alt - variant only found in the proband, ' .
                   'not_found - variant not found in any sample';
  }

  $header{'DeNovo'} = $header_text;

  return \%header;
}

sub run {
  my ($self, $tva, $line) = @_;

  my $vf = $tva->variation_feature;
  my $zyg = defined($line->{Extra}) ? $line->{Extra}->{ZYG} : $line->{ZYG};
  # only interested if we know the zygosity
  return {} if(!$zyg);

  my $chr = $vf->{chr};
  my $start = $vf->{start};
  my @alleles = split /\//, $vf->allele_string;
  my $ref_allele = shift @alleles;
  my $alt_allele = $tva->variation_feature_seq;
  my $n_alt_alleles = scalar(@alleles);

  my $vf_genotype = $vf->{genotype_ind} ? $vf->{genotype_ind} : undef;

  my @result = ();
  my $list_of_ind;

  foreach my $geno_ind (@{$zyg}) {
    my ($ind, $geno) = split(':', $geno_ind);

    # Check if VCF and PED file have the same individual IDs
    if(!$self->{linkage}->{$ind}) {
      die "ERROR: VCF and PED individuals do not match. Please check the individual IDs in PED file: ", join(',', keys %{$self->{linkage}}), "\n";
    }

    my $family_id = $self->{linkage}->{$ind}->{family_id};

    # HOM : homozygous
    # HET : heterozygous
    # HOMREF : homozygous reference (not a variant)
    if($geno eq 'HOM' || $geno eq 'HET') {
      if($self->{linkage}->{$ind} && $self->{linkage}->{$ind}->{parents}) {
        $list_of_ind->{$family_id}->{'child'} = $geno_ind;
      }
      elsif($self->{linkage}->{$ind}) {
        push @{$list_of_ind->{$family_id}->{'parent'}}, $geno_ind;
      }
    }
  }

  foreach my $family (keys %{$self->{family_list}}) {
    if(scalar(keys %{$list_of_ind->{$family}}) == 1) {
      # First check if variant is multi-allelic
      # The checks are different for these variants
      if($n_alt_alleles == 2) {
        my ($child_ind, $child_geno) = split(':', $list_of_ind->{$family}->{'child'}) if($list_of_ind->{$family}->{'child'});
        my $parent_geno = $list_of_ind->{$family}->{'parent'} if($list_of_ind->{$family}->{'parent'});
        my $multi_allelic_r = $self->_check_multi_allelic($ref_allele, $alt_allele, $vf_genotype, $child_ind, $parent_geno);
        push @result, "$family:$multi_allelic_r" if defined $multi_allelic_r;
        if(defined $multi_allelic_r && $multi_allelic_r eq 'de_novo_alt') {
          write_report($self->{report_dir}.'/'.$self->{report_de_novo}->{$family}, $line, $tva, $list_of_ind->{$family}->{'child'});
        }
        elsif(defined $multi_allelic_r && $multi_allelic_r eq 'in_child_and_both_parents') {
          write_report($self->{report_dir}.'/'.$self->{report_all}->{$family}, $line, $tva);
        }
        elsif(defined $multi_allelic_r && $multi_allelic_r eq 'only_in_both_parents') {
          write_report($self->{report_dir}.'/'.$self->{report_both_parents}->{$family}, $line, $tva);
        }
      }
      # If variant is not multi-allelic continue with normal checks
      elsif(defined $list_of_ind->{$family}->{'child'}) {
        push @result, "$family:de_novo_alt";
        write_report($self->{report_dir}.'/'.$self->{report_de_novo}->{$family}, $line, $tva, $list_of_ind->{$family}->{'child'});
      }
      elsif(scalar(@{$list_of_ind->{$family}->{'parent'}}) == 1) {
        push @result, "$family:only_in_one_parent" if $self->{full_report};
      }
      elsif($self->{full_report}) {
        push @result, "$family:only_in_both_parents";
        write_report($self->{report_dir}.'/'.$self->{report_both_parents}->{$family}, $line, $tva);
      }
    }
    elsif(scalar(keys %{$list_of_ind->{$family}}) == 2) {
      # check the parents zygosity
      my ($child_ind, $child_geno) = split(':', $list_of_ind->{$family}->{'child'});
      my $parent_geno = $list_of_ind->{$family}->{'parent'};
      my $parents_het = 0; # number of parents that are heterozygous
      my $parents_hom = 0; # number of parents that are homozygous
      foreach my $p (@{$parent_geno}) {
        my ($p_ind, $p_geno) = split(':', $p);
        if($p_geno eq 'HET') {
          $parents_het += 1;
        }
        if($p_geno eq 'HOM') {
          $parents_hom += 1;
        }
      }

      if(scalar(@{$parent_geno}) == 2) {
        # found in the child and both parents
        if($child_geno eq 'HET' && $parents_hom == 2 && $n_alt_alleles == 1) {
          push @result, "$family:de_novo_alt";
          write_report($self->{report_dir}.'/'.$self->{report_de_novo}->{$family}, $line, $tva, $list_of_ind->{$family}->{'child'});
        }
        elsif($n_alt_alleles == 1) {
          if($self->{full_report}) {
            push @result, "$family:in_child_and_both_parents";
            write_report($self->{report_dir}.'/'.$self->{report_all}->{$family}, $line, $tva);
          }
        }
        # Multi-allelic variants are different
        else {
          my $multi_allelic_r = $self->_check_multi_allelic($ref_allele, $alt_allele, $vf_genotype, $child_ind, $parent_geno);
          push @result, "$family:$multi_allelic_r" if defined $multi_allelic_r;
          if(defined $multi_allelic_r && $multi_allelic_r eq 'de_novo_alt') {
            write_report($self->{report_dir}.'/'.$self->{report_de_novo}->{$family}, $line, $tva, $list_of_ind->{$family}->{'child'});
          }
          elsif(defined $multi_allelic_r && $multi_allelic_r eq 'in_child_and_both_parents') {
            write_report($self->{report_dir}.'/'.$self->{report_all}->{$family}, $line, $tva);
          }
          elsif(defined $multi_allelic_r && $multi_allelic_r eq 'only_in_both_parents') {
            write_report($self->{report_dir}.'/'.$self->{report_both_parents}->{$family}, $line, $tva);
          }
        }
      }
      else {
        # found in the child and one parent
        if($child_geno eq 'HOM' && $parents_het == 1 && $n_alt_alleles == 1) {
          push @result, "$family:in_child_hom_and_one_parent_het" if $self->{full_report};
        }
        elsif($n_alt_alleles == 1) {
          push @result, "$family:in_child_and_one_parent" if $self->{full_report};
        }
        # Multi-allelic variants are different
        else {
          my $multi_allelic_r = $self->_check_multi_allelic($ref_allele, $alt_allele, $vf_genotype, $child_ind, $parent_geno);
          push @result, "$family:$multi_allelic_r" if defined $multi_allelic_r;
          if(defined $multi_allelic_r && $multi_allelic_r eq 'de_novo_alt') {
            write_report($self->{report_dir}.'/'.$self->{report_de_novo}->{$family}, $line, $tva, $list_of_ind->{$family}->{'child'});
          }
          elsif(defined $multi_allelic_r && $multi_allelic_r eq 'in_child_and_both_parents') {
            write_report($self->{report_dir}.'/'.$self->{report_all}->{$family}, $line, $tva);
          }
          elsif(defined $multi_allelic_r && $multi_allelic_r eq 'only_in_both_parents') {
            write_report($self->{report_dir}.'/'.$self->{report_both_parents}->{$family}, $line, $tva);
          }
        }
      }
    }
    else {
      push @result, "$family:not_found";
    }
  }

  my $final = join("&", @result) if(scalar(@result) > 0);

  return $final ? { DeNovo => $final } : {};
}

# This method checks multi-allelic variants
# Checking the genotype returned by --individual_zyg is not enough to determine if variant is de novo
# Example:
#           1 46352728 A/C/G
#           0|2 1|1 0|1
#           HET HOM HET
sub _check_multi_allelic {
  my ($self, $ref_allele, $alt_allele, $vf_genotype, $child_ind, $parent_geno) = @_;

  my $found_child = 0;
  my $n_parents = 0;
  my $result;

  # $vf_genotype has the genotype of all individuals from all families
  # we want to select the individuals for this specific family
  my $parents_alleles;
  my $child_alleles;

  if(defined $child_ind) {
    foreach my $x (@{$vf_genotype->{$child_ind}}) {
      if($x ne $ref_allele && $x eq $alt_allele) {
        $child_alleles->{$x} = 1;
      }
    }
  }

  foreach my $p (@{$parent_geno}) {
    my ($parent_ind, $parent_geno) = split(':', $p);

    foreach my $x (uniq (@{$vf_genotype->{$parent_ind}})) {
      if($x ne $ref_allele && $x eq $alt_allele) {
        $parents_alleles->{$x} += 1;
      }
    }
  }

  if(defined $child_alleles->{$alt_allele}) {
    $found_child = $child_alleles->{$alt_allele};
  }
  if(defined $parents_alleles->{$alt_allele}) {
    $n_parents = $parents_alleles->{$alt_allele};
  }

  if($found_child && !$n_parents) {
    $result = 'de_novo_alt';
  }
  elsif($found_child && $n_parents == 2) {
    $result = 'in_child_and_both_parents';
  }
  elsif($found_child && $n_parents == 1) {
    $result = 'in_child_and_one_parent';
  }
  elsif(!$found_child && $n_parents == 1) {
    $result = 'only_in_one_parent';
  }
  elsif(!$found_child && $n_parents == 2) {
    $result = 'only_in_both_parents';
  }
  else {
    $result = 'not_found';
  }

  return !$self->{full_report} && ($result ne 'de_novo_alt' || $result ne 'not_found') ? undef : $result;
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

  my $cdna = '-';
  my $peptide_start = '-';
  my $codon = '-';
  my $aa_string = '-';
  my $gene = '-';
  my $feature = '-';
  my $feature_type = '-';

  if(!$tva->isa('Bio::EnsEMBL::Variation::IntergenicVariationAllele')) {
    $cdna = defined($tva->transcript_variation) && defined($tva->transcript_variation->cdna_start()) ? $tva->transcript_variation->cdna_start() : '-';
    $peptide_start = defined($tva->transcript_variation) && defined($tva->transcript_variation->translation_start) ? $tva->transcript_variation->translation_start : '-';
    $codon = defined($tva->transcript_variation) && defined($tva->transcript_variation->codons) ? $tva->transcript_variation->codons : '-';
    $aa_string = defined($tva->pep_allele_string) ? $tva->pep_allele_string : '-';
    $gene = $line->{Gene};
    $feature = $line->{Feature};
    $feature_type = $line->{Feature_type};
  }

  my $existing_variation;

  my $ind;
  my $geno;

  open(my $fh, '>>', $file) or die "Could not open file $file $!\n";

  my $out_line = $line->{Uploaded_variation}."\t".$line->{Location}."\t".$line->{Allele}."\t".$gene."\t"
  .$feature."\t".$feature_type."\t".join(',', @{$line->{Consequence}})."\t".$cdna."\t".$peptide_start.
  "\t".$aa_string."\t".$codon;

  if($geno_ind) {
    ($ind, $geno) = split(':', $geno_ind);
    $out_line .= "\tde_novo_sample=".$ind."\tde_novo_zyg=".$geno;
  }

  print $fh $out_line."\n";

  close($fh);
}

1;
