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

 Mastermind

=head1 SYNOPSIS

 mv Mastermind.pm ~/.vep/Plugins
 ./vep -i variations.vcf --plugin Mastermind,/path/to/data.vcf.gz

=head1 DESCRIPTION

 This is a plugin for the Ensembl Variant Effect Predictor (VEP) that
 reports variants that have clinical evidence cited in the medical literature. 
 It is available for both GRCh37 and GRCh38.
 
 The output includes three unique numbers for each variant (MMCNT1, MMCNT2, MMCNT3)
 and one value (MMID3) to be used to build an URL which shows all articles
 from MMCNT3. To build the URL, substitute the 'gene:key' in the following link with the
 value from MMID3:  
 https://mastermind.genomenon.com/detail?disease=all%20diseases&gene=gene&mutation=gene:key 
  
 More information can be found at: https://www.genomenon.com/cvr/


 The following steps are necessary before running this plugin:
 
 Download and Registry (free): 
 https://www.genomenon.com/cvr/ 
 
 GRCh37 VCF:
 unzip mastermind_cited_variants_reference-XXXX.XX.XX-grch37-vcf.zip
 bgzip mastermind_cited_variants_reference-XXXX.XX.XX-GRCh37-vcf
 tabix -p vcf mastermind_cited_variants_reference-XXXX.XX.XX.GRCh37-vcf.gz

 GRCh38 VCF:
 unzip mastermind_cited_variants_reference-XXXX.XX.XX-grch38-vcf.zip
 bgzip mastermind_cited_variants_reference-XXXX.XX.XX-GRCh38-vcf
 tabix -p vcf mastermind_cited_variants_reference-XXXX.XX.XX.GRCh38-vcf.gz
  
 
 The plugin can then be run with:
 ./vep -i variations.vcf --plugin Mastermind,/path/to/mastermind_cited_variants_reference-XXXX.XX.XX.GRChXX-vcf.gz


=cut

package Mastermind;

use strict;
use warnings;
 
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp); 

use Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin); 

sub new {
  my $class = shift;

  my $self = $class->SUPER::new(@_);

  $self->expand_left(0);
  $self->expand_right(0);

  $self->get_user_params();
  
  my $file = $self->params->[0]; 
  $self->{file} = $file; 
    
  return $self;
}

sub feature_types {
  return ['Transcript'];
}

sub get_header_info {

  return{
    'MMCNT1' => 'Count of Mastermind articles with cDNA matches for this specific variant',
    'MMCNT2' => 'Count of Mastermind articles with variants either explicitly matching at the cDNA level or given only at protein level',
    'MMCNT3' => 'Count of Mastermind articles including other DNA-level variants resulting in the same amino acid change',
    'MMID3'  => 'Mastermind variant identifiers, as gene:key, for MMCNT3',
  };

}

sub run {
  my ($self, $tva) = @_;
    
  my $vf = $tva->variation_feature; 
  my $tv = $tva->transcript_variation; 
  my $chr = $vf->{chr};

  my $chr_syn;
  my @new_chr_array;
  my $new_chr;
  
  $self->parse_chromosome_synonyms($self->config->{'synonyms'}) if $self->config->{cache} && (not defined($self->{config}->{_chromosome_synonyms}));
  
  if(defined($self->{syn_cache}->{$chr})) {
    $new_chr = $self->{syn_cache}->{$chr}; 
  }
  else {
    if($self->config->{database}) {
      my $srs_adaptor = $vf->slice->adaptor->db->get_SeqRegionSynonymAdaptor();
      $chr_syn = $srs_adaptor->get_synonyms( $vf->slice->get_seq_region_id($vf->slice) );
      @new_chr_array = map {$_->{name}} (grep {$_->{name} =~ 'NC_'} @{$chr_syn});
    }
    elsif($self->config->{cache}) {
      $chr_syn = $self->config->{_chromosome_synonyms}->{($vf->{chr})};
      @new_chr_array = grep(/NC_/, keys($chr_syn))
    }

    $new_chr = shift(@new_chr_array);
    $self->{syn_cache}->{$chr} = $new_chr;
  }
  
  my @alleles = split /\//, $vf->allele_string;  
  my $ref_allele = shift @alleles; 
  my $alt_allele = shift @alleles; 
  
  my $end = $vf->{end};
  my $start = $vf->{start};
  ($start, $end) = ($end, $start) if $start > $end;
  my $consequence = $vf->consequence_type;  

  my @data = @{$self->get_data($new_chr, $start, $end)}; 
  
  my $result_data;  
  
  foreach my $data_value (@data) {

    if($data_value->{data}) {
      my $ref_allele_comp = $ref_allele; 
      my $alt_allele_comp = $alt_allele; 
      reverse_comp(\$ref_allele_comp); 
      reverse_comp(\$alt_allele_comp);  
     
      # Ref and alt alleles from mastermind file  
      my $mm_ref = $data_value->{ref};
      my $mm_alt = $data_value->{alt};

      if( ($ref_allele eq $mm_ref && $alt_allele eq $mm_alt) || ($ref_allele_comp eq $mm_ref && $alt_allele_comp eq $mm_alt) ) { 
        
        my $peptide_start = defined($tv->translation_start) ? $tv->translation_start : undef;  
        my $aa_alterations = $data_value->{aa};

        foreach my $aa_alteration (@$aa_alterations) {
          $aa_alteration =~ s/.*\:[A-Za-z]+//;
          $aa_alteration =~ s/[A-Za-z]+|\*//;

          # checks if citation refers to an UTR variant 
          my $mm_utr = $aa_alteration =~/UTR/;
          if($mm_utr) { 
            $result_data = $data_value->{result};
          }
          # If there's a protein alteration then it only adds citations for the exact alteration cited  
          elsif(defined($aa_alteration) && defined($peptide_start) && $peptide_start == $aa_alteration) {
            $result_data = $data_value->{result};
          } 
        } 
      } 

    } 

  }

  my $result = defined($result_data) ? $result_data : {};
 
  return $result; 
} 

sub parse_data {
  my ($self, $line) = @_;

  my ($chr, $start, $id, $ref, $alt, $x, $xx, $data) = split /\t/, $line;

  my ($hgvs, $gene, $mmcnt1, $mmcnt2, $mmcnt3, $mmid3, $mmuri3) = split /;/, $data; 
  
  my $mm_data = $mmcnt1 . ';' . $mmcnt2 . ';' . $mmcnt3 . ';' . $mmid3; 

  $mmcnt1 =~ s/MMCNT1=//;
  $mmcnt2  =~ s/MMCNT2=//;
  $mmcnt3  =~ s/MMCNT3=//;
  $mmid3 =~ s/MMID3=//;

  my %mm_hash;
  $mm_hash{'MMCNT1'} = $mmcnt1;
  $mm_hash{'MMCNT2'} = $mmcnt2;
  $mm_hash{'MMCNT3'} = $mmcnt3;
  $mm_hash{'MMID3'}  = $mmid3;

  my @aa_alteration = split /,/, $mmid3;   

  return {
    chr    => $chr, 
    start  => $start, 
    ref    => $ref,
    alt    => $alt,
    aa     => \@aa_alteration,  
    data   => $mm_data,
    result => \%mm_hash, 
  };
}

sub parse_chromosome_synonyms {
  my $self = shift;
  my $file = shift;

  if($file) {
    open INPUT, $file or throw("ERROR: Could not read synonyms file $file: $!");

    my $synonyms = $self->config->{_chromosome_synonyms} ||= {};

    while(<INPUT>) {
      chomp;
      my @split = split(/\s+/, $_);

      my $ref = shift @split;

      foreach my $syn(@split) {
        $synonyms->{$ref}->{$syn} = 1;
        $synonyms->{$syn}->{$ref} = 1;
      }
    }

    close INPUT; 
  }

  return $self->config->{_chromosome_synonyms} ||= {};
}

1;
