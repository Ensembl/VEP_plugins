=head1 CONTACT

 Abhishek Niroula <abhishek.niroula@med.lu.se>
 Mauno Vihinen <mauno.vihinen@med.lu.se>

=cut

=head1 NAME

 PON_P2

=head1 SYNOPSIS

 mv PON_P2.pm ~/.vep/Plugins
 ./vep -i variations.vcf --plugin PON_P2,pyscript=/path/to/python/script/ponp2.py,hg=hg37

=head1 DESCRIPTION
 This plugin for Ensembl Variant Effect Predictor (VEP) computes the predictions of PON-P2
 for amino acid substitutions in human proteins.

 PON-P2 is developed and maintained by Protein Structure and Bioinformatics Group
 at Lund University and is available at http://structure.bmc.lu.se/PON-P2/.

 If you use this data, please cite the following publication
 Niroula, A., Vihinen, M. Harmful somatic amino acid substitutions affect key pathways in cancers.
 BMC Med Genomics 8, 53 (2015). https://doi.org/10.1186/s12920-015-0125-x

 There are two ways to run the plugin:

 1. To compute the predictions from the PON-P2 API, use python script 'ponp2.py' (*)
    and select the reference genome (acceptable values are: hg37 and hg38):
    --plugin PON_P2,pyscript=/path/to/python/script/ponp2.py,hg=hg37
    (*) To run this mode, you will require a python script and its dependencies (Python,
    python suds). The python file can be downloaded from http://structure.bmc.lu.se/PON-P2/vep.html/
    and the complete path to this file must be supplied while using this plugin.

 2. To fetch the predictions from a file containing pre-calculated predictions for somatic variations 
 please use the following key=value option (only available for GRCh37):

   file :         COSMIC text file with pre-calculated predictions downloaded from
                  http://structure.bmc.lu.se/PON-P2/cancer30.html/

   The following steps are necessary before using the file:
   > (head -n 1 COSMIC.txt && tail -n +2 COSMIC.txt | sort -t $'\t' -k1,1 -k2,2n) > cosmic_sorted.txt
   > sed -i 's/Chromosome/#Chromosome/' cosmic_sorted.txt
   > bgzip cosmic_sorted.txt
   > tabix -s 1 -b 2 -e 2 cosmic_sorted.txt.gz

   --plugin PON_P2,file=path/to/cosmic_sorted.txt.gz 

=cut

package PON_P2;


use strict;
use warnings;

use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin;
use Bio::EnsEMBL::Variation::Utils::Sequence qw(get_matched_variant_alleles);
use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin);
use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

sub feature_types {
  return ['Transcript'];
}

sub get_header_info {
  return {
    PON_P2 => "PON-P2 prediction and score for amino acid substitutions"
  };
}

sub new {
  my $class = shift;
  my $self = $class->SUPER::new(@_);

  # get parameters
  #my $command = $self->params->[0];
  #my $Hg = $self->params->[1];
  my $params = $self->params_to_hash();
  my $pyscript;
  my $hg;
  my $from_file;

  if (!%{$params}) {
    $pyscript = $self->params->[0];
    $hg = $self->params->[1];
  } else {
    if (exists $params->{"file"}) {
      $from_file = $params->{file};
      die "\nERROR: No PON_P2 file specified\nTry using 'file=path/to/file.txt.gz'\n" unless -e $from_file;
      $self->add_file($from_file);
      $self->{from_file} = $from_file;
    } else {
      $pyscript = $params->{pyscript} if (defined ($params->{pyscript})); 
      $hg = $params->{hg} if (defined($params->{hg}))
    }
  }

  if(!$from_file) {
    die 'ERROR: Path to python script not specified! Specify path to python script e.g. --plugin PON_P2,/path/to/python/client/for/ponp2.py,[hg37/hg38]\n' unless defined($pyscript);
    die 'ERROR: Reference genome not specified! Specify the reference genome after the path to python file e.g. --plugin PON_P2,/path/to/python/client/for/ponp2.py,[hg37/hg38]\n' unless defined($hg);
    die "ERROR: Wrong reference genome specified! It should be either 'hg37' or 'hg38'\n" unless grep(/^hg3[78]$/, $hg);
    die 'ERROR: Incorrect path to ponp2.py\n' unless -e $pyscript;
  }

  $self->{command} = $pyscript;
  $self->{Hg} = $hg;

  return $self;
}


sub run {
  my ($self, $tva) = @_;

  # only for missense variants
  return {} unless grep {$_->SO_term eq 'missense_variant'} @{$tva->get_all_OverlapConsequences};

  ## Now get the variation features
  my $vf = $tva->variation_feature;

  ## If not snp return
  return {} unless $vf->{start} == $vf->{end};

  my $chr = $vf->{chr};
  my $start = $vf->{start};
  my $end = $vf->{end};
  my $ref_allele = $vf->ref_allele_string;
  my $alt_allele = $tva->base_variation_feature->alt_alleles;

  ## Check for single nucleotide substitution
  return {} unless $ref_allele =~ /^[ACGT]$/;

  if($self->{from_file}) {
    my @data = @{ $self->get_data($chr, $start, $end) };
    return {} unless(@data);

    my $result = 0;

    foreach (@data) {

      my $matches = get_matched_variant_alleles(
        {
          ref    => $ref_allele,
          alts   => $alt_allele,
          pos    => $start,
          strand => $vf->strand
        },
        {
           ref  => $_->{ref},
           alts => [$_->{alt}],
           pos  => $_->{start},
        }
      );

      if (@{$matches}) {
        my $tv = $tva->transcript_variation;
        my $peptide_start = defined($tv->translation_start) ? $tv->translation_start : undef;
        my @vf_aa = split '/', $tva->pep_allele_string;

        if ($_->{aa} eq $vf_aa[0].$peptide_start.$vf_aa[1] || $_->{aa} eq $vf_aa[1].$peptide_start.$vf_aa[0]) {
          my $pred = $_->{pon_p2_pred};
          my $class = $_->{pon_p2_class};
          $result = 1;
          return {
            PON_P2 => "$pred($class)"
          };
        }
        else {
          return {};
        }
      }
    }
    return {} if(!$result) ;
  }
  else {
    foreach my $alt (@{$alt_allele}) {
      my $command = $self->{command};
      my $Hg = $self->{Hg};
      my $V = $chr."_".$start."_".$ref_allele."_".$alt;

      ## Call pon-p2 python script here
      my $ponp2Res = `python2 $command $V $Hg` or return {};
      $ponp2Res =~ s/\R//g;

      my ($pred, $prob) =split /\t/, $ponp2Res;

      ## Can PON-P2 predict?
      return {} if $pred eq "cannot";

      ## Return predictions
      return $pred && $prob ? {
        PON_P2 => "$pred($prob)",
      } : {};
    }
  }
}

sub parse_data {
  my ($self, $line) = @_;

  my @all_data = split /\t/, $line;

  return {
      start => $all_data[1],
      end => $all_data[2],
      ref => $all_data[3],
      alt => $all_data[4],
      aa => $all_data[6],
      pon_p2_pred => $all_data[7],
      pon_p2_class => $all_data[9]
  };
}

sub get_start {
  return $_[1]->{start};
}

sub get_end {
  return $_[1]->{end};
}

1;
