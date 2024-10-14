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

    Xiaolei Zhang
    Ensembl <http://www.ensembl.org/info/about/contact/index.html>

=cut

=head1 NAME

    UTRAnnotator

=head1 SYNOPSIS
 
    mv UTRAnnotator.pm ~/.vep/Plugins
    vep -i variations.vcf --plugin UTRAnnotator,file=/path/to/uORF_starts_ends_GRCh38_PUBLIC.txt

    # skip annotation for variants with a 80% or higher overlap of the UTR
    vep -i variations.vcf --plugin UTRAnnotator,file=/path/to/uORF_starts_ends_GRCh38_PUBLIC.txt,max_overlap=80

=head1 DESCRIPTION
 
    A VEP plugin that annotates the effect of 5' UTR variant especially for variant creating/disrupting upstream ORFs.
    Available for both GRCh37 and GRCh38.

    Options are passed to the plugin as key=value pairs:

    file        : (Required) Path to UTRAnnotator data file:
                  - Download 'uORF_5UTR_GRCh37_PUBLIC.txt' or 'uORF_5UTR_GRCh38_PUBLIC.txt'
                    from https://github.com/Ensembl/UTRannotator
                  - Download from http://sorfs.org
    max_overlap : (Optional) Maximum percentage of overlap between variant and
                  UTR for UTR annotation (default: 100)

    Citation

    About the role of 5'UTR variants in human genetic disease:

    Whiffin, N., Karczewski, K.J., Zhang, X. et al. Characterising the loss-of-function impact of 5’ untranslated region variants in 15,708 individuals. Nat Commun 11, 2523 (2020). https://doi.org/10.1038/s41467-019-10717-9

    About UTRAnnotator:

    The original UTRAnnotator plugin is written by Xiaolei Zhang et al. Later adopted by Ensembl VEP plugins with some changes.
    You can find the original plugin here -
    https://github.com/ImperialCardioGenetics/UTRannotator 

    Please cite the UTRannotator publication alongside the Ensembl VEP if you use this resource -
    Annotating high-impact 5'untranslated region variants with the UTRannotator Zhang, X., Wakeling, M.N., Ware, J.S, Whiffin, N. Bioinformatics; doi: https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btaa783/5905476

=cut


package UTRAnnotator;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);

use List::Util qw(min max);
use Scalar::Util qw(looks_like_number);

use strict;

sub feature_types {
  return ['Transcript'];
}

sub variant_feature_types {
  return ['VariationFeature', 'StructuralVariationFeature'];
}

sub new {
  my $class = shift;
  my $self = $class->SUPER::new(@_);

  printf "\nWarning: no FASTA file specified (--fasta). VEP running will take a long time." if ($self->config->{fasta} eq "");

  my $param_hash = $self->params_to_hash();

  if (-e $param_hash->{file}){
    open my $fh, "<", $param_hash->{file} or die $!;
    my %uORF_evidence;

    while (<$fh>) {
      chomp;
      my ($chr, $pos, $gene, $strand, $type, $stop_pos) = split /\t/;
      next if $chr eq "chr" && $pos eq "pos";

      die "$chr does not exist; Check your reference file." unless ($chr =~ m/^(chr*)/);
      die "Position $pos is not a number; Check your reference file." unless looks_like_number($pos);

      my $key = $chr . ":" . $pos; # chr has 'chr' prefix
      $uORF_evidence{$key} = 1;
    }
    close $fh;
    $self->{uORF_evidence} = \%uORF_evidence;
  } else {
      printf "Warning: small ORF file not found. For human, you could use the curated list of uORFs found in the repository (https://github.com/Ensembl/UTRannotator):\n" .
      "'uORF_starts_ends_GRCh37_PUBLIC.txt' for GRCh37 or 'uORF_starts_ends_GRCh38_PUBLIC.txt' for GRCh38.\n";
  }

  # maximum percentage UTR overlap
  $self->{max_overlap} = defined $param_hash->{max_overlap} ? $param_hash->{max_overlap} : 100;

  # Kozak-Strength Class
  my %kozak_strength;
  $kozak_strength{1}='Weak';
  $kozak_strength{2}='Moderate';
  $kozak_strength{3}='Strong';

  $self->{kozak_strength} = \%kozak_strength;

  return $self;
}

sub get_header_info {

  my $self->{_header_info} = {
        '5UTR_consequence' => 'Variant consequence from UTRAnnotator',
        '5UTR_annotation' => 'Variant annotation from UTRAnnotator',
        Existing_uORFs => 'The number of existing uORFs with a stop codon within the 5 prime UTR',
        Existing_OutOfFrame_oORFs => 'The number of existing out-of-frame overlapping ORFs (OutOfFrame oORF) at the 5 prime UTR',
        Existing_InFrame_oORFs => 'The number of existing inFrame overlapping ORFs (inFrame oORF) at the 5 prime UTR',
  };
  return $self->{_header_info};
}

sub run {
  my ($self, $tva) = @_;

  #only annotate the effect if the variant is 5_prime_UTR_variant
  return {} unless grep {$_->SO_term eq '5_prime_UTR_variant'}  @{$tva->get_all_OverlapConsequences};

  #retrieve the variant info
  my $vf = $tva->can('variation_feature') ? $tva->variation_feature : $tva->structural_variation_feature;
  my $chr     = ($vf->{chr}   || $vf->seq_region_name);
  my $pos     = ($vf->{start} || $vf->seq_region_start);
  my $pos_end = ($vf->{end}   || $vf->seq_region_end);
  my @alleles = split /\//, $vf->allele_string;
  my $ref = shift @alleles;
  my $alt = $tva->can('variation_feature_seq') ? $tva->variation_feature_seq : '';
  my %variant = (
        "chr" => $chr,
        "pos" => $pos,
        "ref" => $ref,
        "alt" => $alt,
  );

  #retrieve the UTR info: transcript id, strand, five prime UTR sequence, start and end genomic coordinates.
  my $t = $tva->transcript;
  my $transcript_id = (defined $t? $t->stable_id: undef);

  #retrieve the gene symbol of the transcript
  my $symbol = $t->{_gene_symbol} || $t->{_gene_hgnc};

  #retrieve the strand of the transcript
  my $tr_strand = $t->strand + 0;
  my $cds = $t->translateable_seq();

  #retrieve the five prime utr sequence
  my $five_prime_seq = (defined $t->five_prime_utr? $t->five_prime_utr->seq(): undef);

  #Type: five_prime_feature - Bio::EnsEMBL::Feature
  my $UTRs = $t->get_all_five_prime_UTRs();
  my @five_utr_starts;
  my @five_utr_ends;
  foreach  my $utr (@$UTRs){
    my $UTR_start = $utr->start(); # this will return the absolute starting positions in chromosome of the UTR exons
    my $UTR_end = $utr->end(); # this will return the absolute ending positions in chromosome of the UTR exons

    #avoid storing UTRs if variant goes over the maximum allowed UTR overlap percentage
    my $UTR_length = $UTR_end - $UTR_start + 1;
    my $UTR_overlap_start = max($UTR_start, $pos);
    my $UTR_overlap_end   = min($UTR_end, $pos_end);
    my $UTR_overlap_perc  = 100 * ($UTR_overlap_end - $UTR_overlap_start) / $UTR_length;
    next if $UTR_overlap_perc >= $self->{max_overlap};

    push(@five_utr_starts, $UTR_start);
    push(@five_utr_ends,   $UTR_end);
  }
  return {} unless (@five_utr_starts || @five_utr_ends);

  my @sorted_starts = sort {$a <=> $b} @five_utr_starts;
  my @sorted_ends = sort {$a <=> $b} @five_utr_ends;

  my %UTR_info = (
        "gene" => $symbol,
        "start" => \@sorted_starts,
        "end" => \@sorted_ends,
        "seq" => $five_prime_seq,   #the exon sequences in 5'UTR
        "strand" => $tr_strand,
        "cds_seq" => $cds,
  );

  $self->{ref_coding} = $self->get_ref_coding($ref);
  $self->{alt_coding} = $self->get_alt_coding($alt,$tr_strand);
  
  my %uAUG_gained = {};
  my %uSTOP_lost = {};
  my %uAUG_lost = {};
  my %uSTOP_gained = {};
  my %uFrameshift = {};
  my ($mut_pos, $end_pos) = $self->get_allele_exon_pos($tr_strand, $pos, $self->{ref_coding}, \%UTR_info);

  if(defined($mut_pos) && defined($end_pos)) {
    my @sequence = split //, $five_prime_seq;
    my $mut_utr_seq = $self->mut_utr_sequence(
        \@sequence,$mut_pos,
        $self->{ref_coding},
        $self->{alt_coding},
        $tr_strand
    );

    my @utr_sequence = split //, $UTR_info{seq};
    my %existing_uORF = %{$self->existing_uORF(\@utr_sequence)};

    my @mut_utr_seq = split //,$mut_utr_seq;
    my %existing_utr_uORF = %{$self->existing_uORF(\@mut_utr_seq)};

    my @start = @{$self->get_ATG_pos(\@utr_sequence)};

    %uAUG_gained = %{$self->uAUG_gained(\%variant,\%UTR_info, $mut_pos, $mut_utr_seq, \%existing_utr_uORF)};
    %uSTOP_lost = %{$self->uSTOP_lost(\%variant,\%UTR_info, $mut_pos, $mut_utr_seq, \%existing_uORF, \%existing_utr_uORF)};
    %uAUG_lost = %{$self->uAUG_lost(\%variant,\%UTR_info, $mut_pos, $mut_utr_seq, \%existing_uORF, \@start)};
    %uSTOP_gained = %{$self->uSTOP_gained(\%variant,\%UTR_info, $mut_pos, $end_pos, $mut_utr_seq, \%existing_uORF, \%existing_utr_uORF, \@start)};
    %uFrameshift = %{$self->uFrameshift(\%variant,\%UTR_info, $mut_pos, $end_pos, $mut_utr_seq, \%existing_uORF, \%existing_utr_uORF, \@start)};

  }

  my %five_prime_flag = (
        "uAUG_gained" => $uAUG_gained{'uAUG_gained_flag'},
        "uSTOP_lost" => $uSTOP_lost{'uSTOP_lost_flag'},
        "uAUG_lost" => $uAUG_lost{'uAUG_lost_flag'},
        "uSTOP_gained" => $uSTOP_gained{'uSTOP_gained_flag'},
        "uFrameshift" =>$uFrameshift{'uFrameShift_flag'},
  );

  my %five_prime_annotation = (
        "uAUG_gained" => $uAUG_gained{"uAUG_gained_effect"},
        "uSTOP_lost" => $uSTOP_lost{"uSTOP_lost_effect"},
        "uAUG_lost" => $uAUG_lost{"uAUG_lost_effect"},
        "uSTOP_gained" => $uSTOP_gained{'uSTOP_gained_effect'},
        "uFrameshift" =>$uFrameshift{'uFrameShift_effect'},
  );


  my $output_five_prime_flag;
  my $output_five_prime_annotation;

  foreach my $flag (keys %five_prime_flag){
    if($five_prime_flag{$flag}){
      $output_five_prime_flag = $five_prime_flag{$flag};
      $output_five_prime_annotation = $five_prime_annotation{$flag};
    }
  };

  my %utr_effect;

  if ($self->{config}->{output_format} eq "json" || $self->{config}->{rest}){

    %utr_effect = (
        "5UTR_consequence" => $output_five_prime_flag,
        "5UTR_annotation" => $output_five_prime_annotation,
    );

  } else {

    my @consequences;

        foreach my $value (keys %{$output_five_prime_annotation}){
          push @consequences, join(":", map { "$_=$output_five_prime_annotation->{$value}{$_}" } keys %{$output_five_prime_annotation->{$value}});
        };

    my $consequence;

    $consequence = join(",", @consequences);

    %utr_effect = (
        "5UTR_consequence" => defined($output_five_prime_flag) ? $output_five_prime_flag : "-",
        "5UTR_annotation" => defined($consequence) && $consequence ne '' ? $consequence : "-",
    );

  }

  my $existing_uORF_num = $self->count_number_ATG($five_prime_seq);
  my $output ={%utr_effect, %$existing_uORF_num};
  return $output? $output: {};
}

##################
# Main functions #
##################


sub uAUG_gained {
  # Description: annotate if a five_prime_UTR_variant creates ATG

  my ($self, $variant_info,$UTR_info, $mut_pos, $mut_utr_seq, $existing_utr_uorf) = @_;

  my $pos = $variant_info->{pos};
  my $ref = $variant_info->{ref};
  my $alt = $variant_info->{alt};

  my @sequence = split //, $UTR_info->{seq};
  my $strand = $UTR_info->{strand};

  #return annotators
  my $uAUG_gained_DistanceToCDS = "";  # the distance between the gained uAUG to CDS
  my $uAUG_gained_KozakContext = "";  # the Kozak context sequence of the gained uAUG
  my $uAUG_gained_KozakStrength = ""; # the Kozak strength of the gained uAUG
  my $uAUG_gained_type = ""; # the type of uORF created - any of the following: uORF, inframe_oORF,OutOfFrame_oORF
  my $uAUG_gained_DistanceFromCap = ""; # the distance between the gained uAUG to the start of the five prime UTR
  my $uAUG_gained_DistanceToStop = ""; #the distance between the gained uAUG to stop codon (could be in CDS)

  #indicate whether the variant creates a ATG
  my $flag = 0;
  my $current_kozak = "";
  my $current_kozak_strength ="";

  #the relative position of input variant in the UTR sequence

  my %result = (); #result is a hash table with two elements: $flag and $output_effects
  my $output_flag = "";
  my %output_effects;

  my $ref_coding = $self->{ref_coding};
  my $alt_coding = $self->{ref_coding};

  my @mut_utr_seq = split //,$mut_utr_seq;
  my $mut_utr_length = @mut_utr_seq;

  #get the nt sequence that might have the pos_A for a new ATG codon
  my @mut_seq = @mut_utr_seq[$mut_pos-2..$mut_pos+length($alt_coding)+1];
  my @mut_atg_pos = @{$self->get_ATG_pos(\@mut_seq)};

  # check whether there is a ATG codon is to check the length of the array
  #if there is a ATG codon, flag it and output the relative positive of A
  my $pos_A;
  if (@mut_atg_pos){
    #get the pos of A (in ATG) in the UTR sequence
    $pos_A=$mut_pos-2+$mut_atg_pos[0];
    $flag=1;
  }

  if ($flag){

    ################################################################################
    #annotator 1: get the distance to the start codon of the main ORF
    ################################################################################
    $uAUG_gained_DistanceToCDS = $mut_utr_length-$pos_A;
    $uAUG_gained_DistanceFromCap = $pos_A;

    ################################################################################
    #annotator 2: determine kozak context;
    ################################################################################
    if ((($pos_A-3)>=0)&&($mut_utr_seq[($pos_A+3)])){
      $current_kozak = $mut_utr_seq[($pos_A-3)].$mut_utr_seq[($pos_A-2)].$mut_utr_seq[$pos_A-1]."ATG".$mut_utr_seq[$pos_A+3];
    }
    else{
      $current_kozak = '-';
    }
    #get the strength of kozak context
    if ($current_kozak !~ /-/){
      my @split_kozak = split //, $current_kozak;
      $current_kozak_strength = 1;
      if ((($split_kozak[0] eq 'A')||($split_kozak[0] eq 'G'))&&($split_kozak[6] eq 'G')){
        $current_kozak_strength = 3;
      }
      elsif ((($split_kozak[0] eq 'A')||($split_kozak[0] eq 'G'))||($split_kozak[6] eq 'G')){
        $current_kozak_strength = 2;
      }
    }

    $uAUG_gained_KozakContext=$current_kozak;
    $uAUG_gained_KozakStrength=$self->{kozak_strength}{$current_kozak_strength}? $self->{kozak_strength}{$current_kozak_strength}:$current_kozak_strength;


    ################################################################################
    #annotator 3: Type of new ORF with respect the main ORF: uORF/Overlapping_inFrame/Overlapping_outofframe
    ################################################################################

    if(exists($existing_utr_uorf->{$pos_A})){ #if there is stop codon within 5'UTR
      $uAUG_gained_type = "uORF";
    }
    else{
      if(($mut_utr_length - $pos_A) % 3){
        $uAUG_gained_type = "OutOfFrame_oORF";
      }
      else{
        $uAUG_gained_type = "inFrame_oORF";
      }
    }

    my @overlapping_seq = split //, $mut_utr_seq.$UTR_info->{cds_seq};
    my %existing_uORF = %{$self->existing_uORF(\@overlapping_seq)};
    if(exists($existing_uORF{$pos_A})){
      my @stop_pos_array = sort{$a<=>$b}@{$existing_uORF{$pos_A}};
      my $stop_pos = $stop_pos_array[0];
      $uAUG_gained_DistanceToStop = $stop_pos-$pos_A;
    } else {
      $uAUG_gained_DistanceToStop = "NA";
    }


    my %uORF_effect = (
            "KozakContext" => $uAUG_gained_KozakContext,
            "KozakStrength" => $uAUG_gained_KozakStrength,
            "DistanceToCDS" => $uAUG_gained_DistanceToCDS,
            "type" => $uAUG_gained_type,
            "DistanceToStop" => $uAUG_gained_DistanceToStop,
            "CapDistanceToStart" => $uAUG_gained_DistanceFromCap,
    );

    $output_flag = "5_prime_UTR_premature_start_codon_gain_variant";
    my $size = (keys %output_effects) + 1;
    $output_effects{$size} = \%uORF_effect;
  }

  $result{'uAUG_gained_flag'} = $output_flag;
  $result{'uAUG_gained_effect'} = \%output_effects;

  return \%result;
}

sub uSTOP_gained {
  # Description: annotate whether a five_prime_UTR_variant creates new stop codon. It only evaluate SNVs.

  my ($self, $variant_info,$UTR_info, $mut_pos, $end_pos, $mut_utr_seq, $existing_ref_uORF, $mut_uORF, $start) = @_;

  my $chr = $variant_info->{chr};
  my $pos = $variant_info->{pos};
  my $ref = $variant_info->{ref};
  my $alt = $variant_info->{alt};

  my @sequence = split //, $UTR_info->{seq};
  my $strand = $UTR_info->{strand};
  my $utr_length = @sequence;

  #return annotators
  my $uSTOP_gained_ref_StartDistanceToCDS = "";  # the distance between the uAUG of the disrupting uORF to CDS
  my $uSTOP_gained_KozakContext = "";  # the Kozak context sequence of the disrupting uORF
  my $uSTOP_gained_KozakStrength = ""; # the Kozak strength of the the disrupting uORF
  my $uSTOP_gained_ref_type = ""; # the type of uORF being disrupted - any of the following: uORF, inframe_oORF,OutOfFrame_oORF
  my $uSTOP_gained_newSTOPDistanceToCDS = ""; # the distance between the gained uSTOP to the start of the CDS
  my $uSTOP_gained_evidence = ""; # the translation evidence of the disrupted uORF

  #indicate whether the variant creates a stop codon
  my $flag = 0;
  my $current_kozak = "";
  my $current_kozak_strength ="";

  #the relative position of input variant in the UTR sequence

  my %result = (); #result is a hash table with two elements: $flag and $output_effects
  my $output_flag = 0;
  my %output_effects;

  my $ref_coding = $self->{ref_coding};
  my $alt_coding = $self->{alt_coding};

  #only evaluate SNVs and MNVs, for deletion and insertion it would be evaluated in the uframeshift
  return{} unless(length($ref_coding) eq length($alt_coding));
  # the length of the 5'UTR won't change so as the coordinates of the nts

  my @mut_utr_seq = split //,$mut_utr_seq;
  my $mut_utr_length = @mut_utr_seq;

  my @start = @{$start};

  if(@start){
    for (my $i=0;$i<@start;$i++){
      $flag=0;
      my $start_pos = $start[$i];

      # to set the range for searching new STOP codon: not including start_codon, (1) for uORFs - start_codon...stop_codon (2) for oORFs - start_codon...end of 5'UTR

      my $check_point = $start_pos;
      # the checking end point of eligible area
      #For uORF: start_pos .. check_point
      #For overlapping ORF: start_pos .. 3' end of 5'UTR sequence
      #if the variant is entirely in the uORF (within 5'UTR
      if(exists($existing_ref_uORF->{$start_pos})) {
        my @stops = sort {$a <=> $b} @{$existing_ref_uORF->{$start_pos}};
        $check_point = $stops[0]-1;
      } else{
        #if the existing uORF is an oORF
        $check_point = $utr_length-1;
      }

      # only check the ones between start and stop codons
      # ignore the cases at the boundary of start and stop codon since the effect on start/stop would be evaluated anyway
      next if(($mut_pos<$start_pos+3)|($end_pos>$check_point));
      #check whether there are new stop codon induced by this mutation
      my @mut_stops;
      next unless(exists($mut_uORF->{$start_pos}));
      @mut_stops = sort {$a <=> $b} @{$mut_uORF->{$start_pos}};

      my $mut_stop = $mut_stops[0];
      if($mut_stop<$check_point){
        $flag=1;
      }

      if($flag){
        #getting the Kozak context and Kozak strength of the start codon
        if ((($start_pos-3)>=0)&&($sequence[($start_pos+3)])){
          $current_kozak = $sequence[($start_pos-3)].$sequence[($start_pos-2)].$sequence[$start_pos-1]."ATG".$sequence[$start_pos+3];
        }
        else{
          $current_kozak = '-';
        }

        if ($current_kozak !~ /-/){
          my @split_kozak = split //, $current_kozak;
          $current_kozak_strength = 1;
          if ((($split_kozak[0] eq 'A')||($split_kozak[0] eq 'G'))&&($split_kozak[6] eq 'G')){
            $current_kozak_strength = 3;
          }
          elsif ((($split_kozak[0] eq 'A')||($split_kozak[0] eq 'G'))||($split_kozak[6] eq 'G')){
            $current_kozak_strength = 2;
          }
        }

        $uSTOP_gained_KozakContext=$current_kozak;
        $uSTOP_gained_KozakStrength=$self->{kozak_strength}{$current_kozak_strength}? $self->{kozak_strength}{$current_kozak_strength}:$current_kozak_strength;

        #the annotation of the original uORF
        if (exists($existing_ref_uORF->{$start_pos})){ #if there is stop codon within 5'UTR
          $uSTOP_gained_ref_type = "uORF"
        } elsif (($utr_length-$start_pos) % 3){
          $uSTOP_gained_ref_type = "OutOfFrame_oORF";
        } else {
          $uSTOP_gained_ref_type = "InFrame_oORF";
        }

        $uSTOP_gained_ref_StartDistanceToCDS = $utr_length - $start_pos;
        $uSTOP_gained_newSTOPDistanceToCDS = $mut_utr_length - $mut_stop;

        #find evidence in added reference file
        $uSTOP_gained_evidence = (exists $self->{uORF_evidence})? $self->find_uorf_evidence($UTR_info,$chr,$start_pos): "NA";

        my %uORF_effect = (
                "ref_type" => $uSTOP_gained_ref_type,
                "ref_StartDistanceToCDS" => $uSTOP_gained_ref_StartDistanceToCDS,
                "newSTOPDistanceToCDS" => $uSTOP_gained_newSTOPDistanceToCDS,
                "KozakContext" => $uSTOP_gained_KozakContext,
                "KozakStrength" => $uSTOP_gained_KozakStrength,
                "Evidence" => $uSTOP_gained_evidence,
              );

        $output_flag = "5_prime_UTR_uORF_stop_codon_gain_variant";
        my $size = (keys %output_effects) + 1;
        $output_effects{$size} = \%uORF_effect;
      }

    }
  }

  $result{'uSTOP_gained_flag'} = $output_flag;
  $result{'uSTOP_gained_effect'} = \%output_effects;
  return \%result;

}

sub uSTOP_lost {

  # Description: annotate if a five_prime_UTR_varint removes a stop codon of an existing uORF (given that uORF does not not change)

  my ($self, $variant_info, $UTR_info, $mut_pos, $mut_utr_seq, $existing_uORF, $mut_uORF) = @_;

  my $chr = $variant_info->{chr};
  my $pos = $variant_info->{pos};
  my $ref = $variant_info->{ref};
  my $alt = $variant_info->{alt};

  my @sequence = split //, $UTR_info->{seq};
  my $strand = $UTR_info->{strand};

  #return annotators
  my $uSTOP_lost_AltStop = "";  #whether there is an alternative stop codon downstream
  my $uSTOP_lost_AltStopDistanceToCDS = ""; # the distance between alternative stop codon and the CDS
  my $uSTOP_lost_FrameWithCDS = ""; # the frame of the uORF with respect to CDS
  my $uSTOP_lost_KozakContext = ""; # the Kozak context sequence of the uORF with the lost stop codon
  my $uSTOP_lost_KozakStrength = ""; # the Kozak strength of the uORF with the lost stop codon
  my $uSTOP_lost_evidence = ""; # whether this uORF has any evidence of being translated
  ###

  my $current_kozak = "";
  my $current_kozak_strength = "";

  #indicate whether the variant ever disrupts a stop codon
  my $output_flag = "";
  #indicate whether the variant ever disrupts the uORF that being analyzed
  my $flag_uORF=0;

  #the relative position of input variant in the UTR sequence
  my %result = (); #result is a hash table with two elements: $flag and $output_effects
  my %output_effects;

  my @stop_codons = ("TAA","TGA","TAG");

  my $ref_coding = $self->{ref_coding};
  my $alt_coding = $self->{alt_coding};

  my @mut_utr_seq = split //,$mut_utr_seq;
  my $length = @mut_utr_seq;

  my @start = sort {$a <=> $b} (keys %{$existing_uORF});

  if(%{$existing_uORF}){
    for (my $i=0;$i<@start;$i++){
      $flag_uORF=0;
      my $start_pos = $start[$i];
      my @stops = sort {$a <=> $b} @{$existing_uORF->{$start_pos}};
      my $stop_pos=$stops[0];

      next if ($mut_pos-$stop_pos>2);

      #for snps and deletion
      next if (defined($ref_coding) & $mut_pos+length($ref_coding)-1<$stop_pos);

      next if ($mut_pos<$stop_pos & !defined($ref_coding));   #for insertion

      #for deletion, it definitely disrupting the stop codon.
      if (length($alt_coding) eq 0) {
        $flag_uORF=1;
      } else {
        my $mut_codon = $mut_utr_seq[$stop_pos].$mut_utr_seq[$stop_pos+1].$mut_utr_seq[$stop_pos+2];

        next if (grep( /^$mut_codon$/, @stop_codons));
        $flag_uORF=1;
      }

      if($flag_uORF){
        #getting the Kozak context and Kozak strength of the start codon

        if ((($start_pos-3)>=0)&&($sequence[($start_pos+3)])){
          $current_kozak = $sequence[($start_pos-3)].$sequence[($start_pos-2)].$sequence[$start_pos-1]."ATG".$sequence[$start_pos+3];
        }
        else{
          $current_kozak = '-';
        }

        if ($current_kozak !~ /-/){
          my @split_kozak = split //, $current_kozak;
          $current_kozak_strength = 1;
          if ((($split_kozak[0] eq 'A')||($split_kozak[0] eq 'G'))&&($split_kozak[6] eq 'G')){
            $current_kozak_strength = 3;
          }
          elsif ((($split_kozak[0] eq 'A')||($split_kozak[0] eq 'G'))||($split_kozak[6] eq 'G')){
            $current_kozak_strength = 2;
          }
        }

        $uSTOP_lost_KozakContext=$current_kozak;
        $uSTOP_lost_KozakStrength=$self->{kozak_strength}{$current_kozak_strength}? $self->{kozak_strength}{$current_kozak_strength}:$current_kozak_strength;

        # the sequence before mut_pos should not be changed. Thus start_pos shall still correspond to a start codon;
        # if the sequence is indeed very short as such ATGTGA

        my @mut_stops;
        if(exists($mut_uORF->{$start_pos})){
          @mut_stops = sort {$a <=> $b} @{$mut_uORF->{$start_pos}}
        }

        if (@mut_stops>0){
          $uSTOP_lost_AltStop = "True";
          $uSTOP_lost_AltStopDistanceToCDS = $length-$mut_stops[0];
        } #if there is no alternative stop codon
        else{
          $uSTOP_lost_AltStop = "False";
          $uSTOP_lost_AltStopDistanceToCDS = "NA";
        }
        if (($length-$start_pos) % 3){
          $uSTOP_lost_FrameWithCDS = "outOfFrame";
        }
        else {
          $uSTOP_lost_FrameWithCDS = "inFrame";
        }
        #find evidence from sorf

        $uSTOP_lost_evidence= (exists $self->{uORF_evidence})? $self->find_uorf_evidence($UTR_info,$chr,$start_pos): "NA";

        my %uORF_effect = (
          "AltStop" => $uSTOP_lost_AltStop,
          "AltStopDistanceToCDS" => $uSTOP_lost_AltStopDistanceToCDS,
          "FrameWithCDS" => $uSTOP_lost_FrameWithCDS,
          "KozakContext" => $uSTOP_lost_KozakContext,
          "KozakStrength" => $uSTOP_lost_KozakStrength,
          "Evidence" => $uSTOP_lost_evidence,
        );

        $output_flag = "5_prime_UTR_uORF_stop_codon_loss_variant";
        my $size = (keys %output_effects) + 1;
        $output_effects{$size} = \%uORF_effect;
      }
    }
  }

  $result{'uSTOP_lost_flag'} = $output_flag;
  $result{'uSTOP_lost_effect'} = \%output_effects;

  return \%result;
}

sub uAUG_lost {
  # Description: annotate if a five_prime_UTR_varint removes a start codon of an existing uORF

  my ($self, $variant_info, $UTR_info, $mut_pos, $mut_utr_seq, $existing_ref_uORF, $start) = @_;

  my $chr = $variant_info->{chr};
  my $pos = $variant_info->{pos};
  my $ref = $variant_info->{ref};
  my $alt = $variant_info->{alt};

  my @sorted_starts = @{$UTR_info->{start}};
  my @sequence = split //, $UTR_info->{seq};
  my $utr_length = @sequence;
  my $strand = $UTR_info->{strand};

  #return annotators
  my $uAUG_lost_type = "";   # the uORF type with the reference allele - uORF, inframe_oORF, outOfFrame_oORF
  my $uAUG_lost_DistanceToCDS = ""; # the distance of the lost uAUG to CDS
  my $uAUG_lost_DistanceToSTOP = ""; # the distance of the lost uAUG to stop codon (could be in CDS)
  my $uAUG_lost_KozakContext = ""; # the Kozak context sequence of the lost uAUG
  my $uAUG_lost_KozakStrength = ""; # the strength of KozakContext of the lost uAUG
  my $uAUG_lost_evidence = ""; # whether this uAUG has any evidence of being translated

  my $current_kozak = "";
  my $current_kozak_strength = "";

  #indicate whether the variant ever disrupts a stop codon
  my $output_flag = "";
  #indicate whether the variant ever disrupts the uORF that being analyzed
  my $flag_uORF=0;

  #the relative position of input variant in the UTR sequence

  my %result = (); #result is a hash table with two elements: $flag and $output_effects
  my %output_effects;

  my $ref_coding = $self->{ref_coding};
  my $alt_coding = $self->{alt_coding};

  my @mut_utr_seq = split //,$mut_utr_seq;

  my @start = @{$start};

  for (my $i=0;$i<@start;$i++){
    $flag_uORF=0;
    my $start_pos = $start[$i];

    next if ($mut_pos-$start_pos>2);

    if (length($ref_coding) ne 0){
      next if ($mut_pos+length($ref_coding)-1<$start_pos);
    }
    next if ($mut_pos<$start_pos);   #for insertion

    #for deletion, it definitely disrupting the uORF.

    if (length($alt_coding) eq 0) {
      $flag_uORF=1;
    } else {
      my $mut_codon = $mut_utr_seq[$start_pos].$mut_utr_seq[$start_pos+1].$mut_utr_seq[$start_pos+2];
      next if ($mut_codon eq "ATG");
      $flag_uORF=1;
    }

    if($flag_uORF){

      #getting the Kozak context and Kozak strength of the lost start codon
      if ((($start_pos-3)>=0)&&($sequence[($start_pos+3)])){
        $current_kozak = $sequence[($start_pos-3)].$sequence[($start_pos-2)].$sequence[$start_pos-1]."ATG".$sequence[$start_pos+3];
      } else {
        $current_kozak = '-';
      }

      if ($current_kozak !~ /-/){
        my @split_kozak = split //, $current_kozak;
        $current_kozak_strength = 1;
        if ((($split_kozak[0] eq 'A')||($split_kozak[0] eq 'G'))&&($split_kozak[6] eq 'G')){
          $current_kozak_strength = 3;
        } elsif ((($split_kozak[0] eq 'A')||($split_kozak[0] eq 'G'))||($split_kozak[6] eq 'G')){
          $current_kozak_strength = 2;
        }
      }
      $uAUG_lost_KozakContext=$current_kozak;
      $uAUG_lost_KozakStrength=$self->{kozak_strength}{$current_kozak_strength}? $self->{kozak_strength}{$current_kozak_strength}:$current_kozak_strength;

      #check what kind of uORF does that correspond to?
      #first check whether it's overlapping with CDS

      if (exists($existing_ref_uORF->{$start_pos})){ #if there is stop codon within 5'UTR
        $uAUG_lost_type = "uORF"
      } elsif(($utr_length-$start_pos) % 3){
        $uAUG_lost_type = "OutOfFrame_oORF";
      } else {
        $uAUG_lost_type = "InFrame_oORF";
      }

      $uAUG_lost_DistanceToCDS = $utr_length - $start_pos;

      my @overlapping_seq = split //, $UTR_info->{seq}.$UTR_info->{cds_seq};
      my %existing_overlapping_uORF = %{$self->existing_uORF(\@overlapping_seq)};

      if(exists($existing_overlapping_uORF{$start_pos})){
        my @stop_pos_array = sort{$a<=>$b}@{$existing_overlapping_uORF{$start_pos}};
        my $stop_pos = $stop_pos_array[0];
        $uAUG_lost_DistanceToSTOP = $stop_pos-$start_pos;
      } else {
        $uAUG_lost_DistanceToSTOP = "NA"
      }

      $uAUG_lost_evidence= (exists $self->{uORF_evidence})? $self->find_uorf_evidence($UTR_info,$chr,$start_pos): "NA";

      my %uORF_effect = (
          "type" => $uAUG_lost_type,
          "CapDistanceToStart" =>$start_pos,
          "DistanceToCDS" => $uAUG_lost_DistanceToCDS,
          "DistanceToStop" => $uAUG_lost_DistanceToSTOP,
          "KozakContext" => $uAUG_lost_KozakContext,
          "KozakStrength" => $uAUG_lost_KozakStrength,
          "Evidence" => $uAUG_lost_evidence,
        );

      $output_flag = "5_prime_UTR_premature_start_codon_loss_variant";
      my $size = (keys %output_effects) + 1;
      $output_effects{$size} = \%uORF_effect;
    }
  }

  $result{'uAUG_lost_flag'} = $output_flag;
  $result{'uAUG_lost_effect'} = \%output_effects;

  return \%result;
}

sub uFrameshift {

  # Description: annotate if a five_prime_UTR_varint create a frameshift in existing uORFs

  my ($self, $variant_info, $UTR_info, $mut_pos, $end_pos, $mut_utr_seq, $existing_uORF, $mut_uORF, $start) = @_;

  my $chr = $variant_info->{chr};
  my $pos = $variant_info->{pos};
  my $ref = $variant_info->{ref};
  my $alt = $variant_info->{alt};

  my @sequence = split //, $UTR_info->{seq};
  my $strand = $UTR_info->{strand};
  my $utr_length = @sequence;

  #return annotators
  my $uFrameshift_ref_type = ""; # the type of uORF with the reference allele
  my $uFrameshift_ref_type_length = ""; # the length of the uORF with of the reference allele.
  my $uFrameshift_StartDistanceToCDS = ""; # the distance between the start codon of the disrupted uORF and CDS
  my $uFrameshift_alt_type = ""; # the type of uORF with the alternative allele
  my $uFrameshift_alt_type_length = ""; # the length of the uORF with of the alternative allele
  my $uFrameshift_KozakContext = ""; # the Kozak context sequence of the disrupted uORF
  my $uFrameshift_KozakStrength = ""; # the Kozak strength of the disrupted uORF
  my $uFrameshift_evidence = ""; # whehter there is translation evidence for the disrupted uORF

  #indicate whether the variant ever introduce a frameshift variant
  my $flag_uORF;
  my $output_flag = "";

  my $current_kozak="";
  my $current_kozak_strength="";

  my %result = (); #result is a hash table with two elements: $flag and $output_effects
  my %output_effects;

  my $ref_coding = $self->{ref_coding};
  my $alt_coding = $self->{alt_coding};

  #skip alleles with same length
  return {} unless(length($ref_coding) ne length($alt_coding));

  #if it's a deletion at the boundary of exon and intron, we would skip the annotation

  my @mut_utr_seq = split //,$mut_utr_seq;
  my $length = @mut_utr_seq;

  my @start = @{$start};

  #check for each uORF
  if(@start){
    for (my $i=0;$i<@start;$i++){
      $flag_uORF=0;
      my $start_pos = $start[$i];

      my $check_point = $start_pos;
      # the checking end point of eligible area
      #For uORF: start_pos .. check_point
      #For overlapping ORF: start_pos .. 3' end of 5'UTR sequence

      #if the variant is entirely in the uORF (within 5'UTR)

      if(exists($existing_uORF->{$start_pos})) {
        my @stops = sort {$a <=> $b} @{$existing_uORF->{$start_pos}};
        $check_point = $stops[0]-1;
      } else {
        #if the existing uORF is an oORF
        $check_point = $utr_length-1;
      }

      # only check the ones between start and stop codons and
      # ignore the cases at the boundary of start and stop codon since the effect on start/stop would be evaluated anyway
      if(($mut_pos>=$start_pos+3)&($end_pos<=$check_point)){
        #snp and insertion could only be annotated here
        $flag_uORF = abs(length($ref_coding)-length($alt_coding))%3;
      }

      if($flag_uORF){
        #getting the Kozak context and Kozak strength of the start codon
        if ((($start_pos-3)>=0)&&($sequence[($start_pos+3)])){
          $current_kozak = $sequence[($start_pos-3)].$sequence[($start_pos-2)].$sequence[$start_pos-1]."ATG".$sequence[$start_pos+3];
        }
        else{
          $current_kozak = '-';
        }

        if ($current_kozak !~ /-/){
          my @split_kozak = split //, $current_kozak;
          $current_kozak_strength = 1;
          if ((($split_kozak[0] eq 'A')||($split_kozak[0] eq 'G'))&&($split_kozak[6] eq 'G')){
            $current_kozak_strength = 3;
          }
          elsif ((($split_kozak[0] eq 'A')||($split_kozak[0] eq 'G'))||($split_kozak[6] eq 'G')){
            $current_kozak_strength = 2;
          }
        }

        $uFrameshift_KozakContext=$current_kozak;
        $uFrameshift_KozakStrength=$self->{kozak_strength}{$current_kozak_strength}? $self->{kozak_strength}{$current_kozak_strength}:$current_kozak_strength;

        #the annotation of the original uORF
        my @ref_overlapping_seq = split //, $UTR_info->{seq}.$UTR_info->{cds_seq};
        my %ref_existing_oORF = %{$self->existing_uORF(\@ref_overlapping_seq)};

        if (exists($existing_uORF->{$start_pos})){ #if there is stop codon within 5'UTR
          $uFrameshift_ref_type = "uORF";
        } elsif (($utr_length-$start_pos) % 3) {
          $uFrameshift_ref_type = "OutOfFrame_oORF";
        } else{
          $uFrameshift_ref_type = "InFrame_oORF";
        }

        if (exists($ref_existing_oORF{$start_pos})){
          my @stops = sort {$a <=> $b} @{$ref_existing_oORF{$start_pos}};
          $uFrameshift_ref_type_length = $stops[0]-$start_pos+3;
        } else {
          $uFrameshift_ref_type_length = "NA";
        }

        $uFrameshift_StartDistanceToCDS = $utr_length - $start_pos;

        my @alt_overlapping_seq = split //, $mut_utr_seq.$UTR_info->{cds_seq};
        my %alt_existing_oORF = %{$self->existing_uORF(\@alt_overlapping_seq)};

        if (exists($mut_uORF->{$start_pos})){
          $uFrameshift_alt_type = "uORF";
        } elsif(($length-$start_pos)%3){
          $uFrameshift_alt_type = "OutOfFrame_oORF";
        } else {
          $uFrameshift_alt_type = "InFrame_oORF";
        }

        #get the length
        if (exists($alt_existing_oORF{$start_pos})){
          my @stops = sort {$a <=> $b} @{$alt_existing_oORF{$start_pos}};
          $uFrameshift_alt_type_length = $stops[0]-$start_pos+3;
        } else {
          $uFrameshift_alt_type_length = "NA";
        }

        #find evidence in added reference file
        $uFrameshift_evidence= (exists $self->{uORF_evidence})? $self->find_uorf_evidence($UTR_info,$chr,$start_pos): "NA";

        my %uORF_effect = (
            "ref_type" => $uFrameshift_ref_type,
            "ref_type_length" => $uFrameshift_ref_type_length,
            "ref_StartDistanceToCDS" => $uFrameshift_StartDistanceToCDS,
            "alt_type" => $uFrameshift_alt_type,
            "alt_type_length" => $uFrameshift_alt_type_length,
            "KozakContext" => $uFrameshift_KozakContext,
            "KozakStrength" => $uFrameshift_KozakStrength,
            "Evidence" => $uFrameshift_evidence,
        );

        $output_flag = "5_prime_UTR_uORF_frameshift_variant";
        my $size = (keys %output_effects) + 1;
        $output_effects{$size} = \%uORF_effect;
      }
    }
  }

  $result{'uFrameShift_flag'} = $output_flag;
  $result{'uFrameShift_effect'} = \%output_effects;

  return \%result;
}

#########
# Utils #
#########

=head2 count_number_ATG

  Arg [1]    : $seq
  Description: Count the number of existing ATGs in the five prime UTR sequence.
  Return: returns the number of uORFs and the number of oORFs
  Returntype : hashref

=cut

sub count_number_ATG {
  my ($self,$seq) = @_;

  my @sequence = split //, $seq;
  my $length = @sequence;

  my @atg_pos = @{$self->get_ATG_pos(\@sequence)};
  my @mes_pos = @{$self->get_stopcodon_pos(\@sequence)};

  my $inframe_stop_num=0;
  my $outofframe_atg_num=0;
  my $inframeORF_num=0;

  foreach my $atg (@atg_pos){

    my $flag=0;

    #indicate whether there is a stop codon with respect to this ATG
    foreach my $mes (@mes_pos){
      unless (($mes-$atg) % 3 || ($mes-$atg) < 0 ){
        $flag=1;
        last;
      };
    }
    #if there is no stop codon, then look at whether it's Out_of_frame or Inframe
    if($flag==1){
      $inframe_stop_num++;
    } else {
      (($length-$atg) % 3)? $outofframe_atg_num++ : $inframeORF_num++;
    }

  }

  my %existing_uORF_num = (
    "Existing_uORFs" => $inframe_stop_num,
    "Existing_OutOfFrame_oORFs" => $outofframe_atg_num,
    "Existing_InFrame_oORFs" => $inframeORF_num,
  );

  return \%existing_uORF_num;

}


=head2 existing_uORF

  Arg [1]    : $seq
  Description: obtaining the relative coordinates of start and end pos of existing uORF in the five prime UTR sequence.
  Return: returns as the key the position of the first nucleotide of start codon and the value is all the positions of the first nucleotide of the stop codon.
  Returntype : hashref

=cut


sub existing_uORF {
  my ($self,$seq) = @_;

  my @atg_pos = @{$self->get_ATG_pos($seq)};
  my @mes_pos = @{$self->get_stopcodon_pos($seq)};
  my %uORF;

  foreach my $atg (@atg_pos){
    my @stop_pos;
    foreach my $mes (@mes_pos){
      push @stop_pos, $mes unless ((($mes-$atg) % 3) || ($mes-$atg)<0);
    }

    $uORF{$atg} = [@stop_pos] if (@stop_pos);
  }

  return \%uORF;
}


=head2 get_ATG_pos

  Arg [1]    : $seq
  Description:  get all the relative position of A in ATG of the five prime UTR sequence
  Return: returns ATG position in sequence.
  Returntype : listref

=cut


sub get_ATG_pos {
  my ($self,$seq) = @_;

  my @sequence = @{$seq};
  my $length = @sequence;
  my $seq_str=join '', @sequence;
  my @atg_pos = grep {(substr ($seq_str, $_,3) eq 'ATG')} 0..($length);

  return \@atg_pos;
}


=head2 get_stopcodon_pos

  Arg [1]    : $seq
  Description:  get all the relative position of stop codons in the five prime UTR sequence
  Return: returns stop-codon position in sequence.
  Returntype : listref

=cut


sub get_stopcodon_pos {
  my ($self,$seq) = @_;

  my @sequence = @{$seq};
  my $length = @sequence;

  my @met_pos;
  if($length){
    for (my $seq_n=0; $seq_n<$length; $seq_n++){
      if ((($sequence[$seq_n] eq 'T')&&($sequence[$seq_n+1] eq 'A')&&($sequence[$seq_n+2] eq 'A'))
          ||(($sequence[$seq_n] eq 'T')&&($sequence[$seq_n+1] eq 'A')&&($sequence[$seq_n+2] eq 'G'))
          ||(($sequence[$seq_n] eq 'T')&&($sequence[$seq_n+1] eq 'G')&&($sequence[$seq_n+2] eq 'A'))){
           push @met_pos,$seq_n;
      }
    }
  }

  return \@met_pos;
}

sub get_ref_coding {
  my ($self,$ref) = @_;

  my $ref_coding = $ref;

  if($ref_coding eq "-"){   #for insertion
    $ref_coding = "";
  }

  return $ref_coding;
}

sub get_alt_coding {

  my ($self,$alt,$strand) = @_;
  my $alt_coding = $alt;

  $alt_coding = "" if($alt_coding eq "-");

  reverse_comp(\$alt_coding) if($strand < 0);

  return $alt_coding;
}


=head2 mut_utr_sequence

  Description:  get the mutated five prime UTR sequence
  Return: Mutated UTR sequence
  Returntype : String

=cut


sub mut_utr_sequence {
  my ($self,$seq,$mut_pos,$ref_coding,$alt_coding,$strand) = @_;

  my @sequence = @{$seq};
  my $length = @sequence;

  #order from 5' to 3'
  my $mut_sequence;
  $mut_sequence = (join "", @sequence[0..$mut_pos-1]).$alt_coding.(join "", @sequence[$mut_pos+length($ref_coding)..$length-1]);

  return $mut_sequence;
}

sub transform_hash_to_string {
  my ($self,$hash) = @_;

  my %hash = %{$hash};

  my $output_str = "";
  foreach my $key (sort keys %hash){
    $output_str = $output_str.$key.":".$hash{$key}.",";
  };

  chop($output_str) if($output_str);

  return $output_str;
}

sub utr_exon_position {
  my ($self,$UTR_info) = @_;

  my @sorted_starts = @{$UTR_info->{start}};
  my @sorted_ends = @{$UTR_info->{end}};
  my $num_exons = @sorted_starts;
  my $strand = $UTR_info->{strand};
  my %chr_position;
  my $utr_position = 0;

  if ($strand == 1){
    #create a map from chromosome position to UTR position
    for (my $m=0; $m<$num_exons; $m++){
      for (my $p=$sorted_starts[$m]; $p<=$sorted_ends[$m]; $p++){
        $chr_position{$p}=$utr_position;
        $utr_position++;
      }
    }
  }

  if ($strand == -1){
    #create a map from chromosome position to UTR position
    #the exons were arranged in increasing order above
    for (my $m=$num_exons-1; $m>=0; $m--){
      for (my $p=$sorted_ends[$m]; $p>=$sorted_starts[$m]; $p--){
        $chr_position{$p}=$utr_position;
        $utr_position++;
      }
    }
  }

  return \%chr_position;
}

sub chr_position {
  my ($self,$UTR_info) = @_;

  my @sorted_starts = @{$UTR_info->{start}};
  my @sorted_ends = @{$UTR_info->{end}};
  my $num_exons = @sorted_starts;
  my $strand = $UTR_info->{strand};
  my %utr_position;
  my $utr_position = 0;

  if ($strand == 1){
  #create a map from chromosome position to UTR position
    for (my $m=0; $m<$num_exons; $m++){
      for (my $p=$sorted_starts[$m]; $p<=$sorted_ends[$m]; $p++){
        $utr_position{$utr_position}=$p;
        $utr_position++;
      }
    } 
  }

  if ($strand == -1){
    #create a map from chromosome position to UTR position
    #the exons were arranged in increasing order above
    for (my $m=$num_exons-1; $m>=0; $m--){
      for (my $p=$sorted_ends[$m]; $p>=$sorted_starts[$m]; $p--){
        $utr_position{$utr_position}=$p;
        $utr_position++;
      }
    }
  }

  return \%utr_position;
}

sub get_allele_exon_pos {
  # return the 5' start and end position at the exon of the ref allele
  my ($self, $strand, $pos, $ref_coding, $UTR_info) = @_;
  my %utr_exon_pos = %{$self->utr_exon_position($UTR_info)};

  my $ref_start;
  my $ref_end;

  if ($strand == 1) {
    $ref_start = $utr_exon_pos{$pos};

    if(length($ref_coding)){
      $ref_end = $utr_exon_pos{$pos + length($ref_coding) - 1};
    } else {
      $ref_end = $ref_start;
    }
  }
  else {
    $ref_start = $utr_exon_pos{$pos+length($ref_coding)-1};

    if(length($ref_coding)){
      $ref_end = $utr_exon_pos{$pos};
    } else {
      $ref_end = $ref_start;
    }
  }

  return ($ref_start, $ref_end);
}

sub find_uorf_evidence {
  my ($self,$UTR_info,$chr,$start_pos)=@_;

  my %utr_pos = %{$self->chr_position($UTR_info)};
  my $start_chr_pos = $utr_pos{$start_pos};

  my $query = ($chr=~/chr/i)? $chr . ":" . $start_chr_pos : "chr" . $chr . ":" . $start_chr_pos;
  my $evidence = $self->{uORF_evidence}->{$query}? "True" : "False";

  return $evidence;
}

1;
