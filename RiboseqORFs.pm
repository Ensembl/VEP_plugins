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

 RiboseqORFs

=head1 SYNOPSIS

 ./vep -i variations.vcf --plugin RiboseqORFs,file=/path/to/Ribo-seq_ORFs.bed.gz

=head1 DESCRIPTION

 This is a VEP plugin that uses a standardized catalog of human Ribo-seq ORFs to
 re-calculate consequences for variants located in these translated regions.

 This plugin reports new consequences based on the evidence from the Ribo-seq
 ORF annotation and supporting publications. The human Ribo-seq ORF data can be
 downloaded from: https://ftp.ebi.ac.uk/pub/databases/gencode/riboseq_orfs/data

 After downloading the annotation, please bgzip and tabix it:
   bgzip Ribo-seq_ORFs.bed
   tabix Ribo-seq_ORFs.bed.gz

 For optimal performance when running this plugin in VEP, please use a FASTA
 file (`--fasta`). A FASTA file is always required in offline mode.

 Please cite the publication for the Ribo-seq ORF annotation alongside the VEP
 if you use this resource: https://doi.org/10.1038/s41587-022-01369-0

 The tabix utility must be installed in your path to use this plugin.

=cut

package RiboseqORFs;

use strict;
use warnings;

use Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin;
use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin);

use Bio::EnsEMBL::VEP::Utils qw(format_coords);

sub new {
  my $class = shift;
  my $self  = $class->SUPER::new(@_);

  $self->expand_left(0);
  $self->expand_right(0); 
  $self->get_user_params();

  # Check if using FASTA file
  if(!defined($self->{config}->{fasta})) {
    # Error when using offline mode
    die("ERROR: Cannot generate sequences in offline mode without FASTA file (--fasta)\n")
      if defined $self->{config}->{offline};

    # Warn when using database instead of FASTA file
    warn("WARNING: RiboseqORFs plugin will fetch sequences from Ensembl database; use a FASTA file (--fasta) for optimal performance\n")
      if defined $self->{config}->{cache} and !defined $self->{config}->{quiet};
  }

  # Add plugin data file
  my $param_hash = $self->params_to_hash();
  my $file = $param_hash->{file};
  die "\n  ERROR: No file specified\nTry using 'file=path/to/Ribo-seq_ORFs.bed.gz'\n"
     unless defined($file);
  $self->add_file($file);

  return $self;
}

sub feature_types {
  return [ 'Feature', 'Intergenic' ];
}

sub get_header_info {
  return {
    RiboseqORFs_id                => 'Ribo-seq ORF identifier',
    RiboseqORFs_consequences      => 'Recalculated consequences based on Ribo-seq ORF evidence',
    RiboseqORFs_cDNA_position     => 'Recalculated position of base pair in cDNA sequence',
    RiboseqORFs_CDS_position      => 'Recalculated position of base pair in coding sequence',
    RiboseqORFs_protein_position  => 'Recalculated position of amino acid in protein',
    RiboseqORFs_amino_acids       => 'Reference and recalculated variant codon sequence',
    RiboseqORFs_codons            => 'Reference and recalculated variant amino acids',
    RiboseqORFs_impact            => 'Subjective impact classification of recalculated consequence type',
    RiboseqORFs_publications      => 'Associated publications',
  };
}

sub _transcripts_match {
  my ($tva, $transcripts) = @_;

  # Get transcript ID for Ensembl
  my $tr = $tva->transcript->{stable_id};

  my @all_trs = split(/,/, $transcripts);
  return grep { $tr eq $_ } @all_trs;
}

sub _recreate_TranscriptVariationAllele_with_ORF {
  my $tva = shift;
  my $orf = shift;

  my $slice  = $tva->transcript->slice;
  my $strand = $_->{strand} eq '+' ? 1 : -1;

  # Create one exon for each ORF block
  my $exons;
  my @block_start = split(',', $orf->{chromStarts});
  my @block_size  = split(',', $orf->{blockSize});
  my @exon_frames = split(',', $orf->{exonFrames});

  for my $k ( 0 .. $orf->{blockCount} - 1 ) {
    my $exon_start = $orf->{chromStart} + $block_start[$k];
    my $exon_end   = $exon_start + $block_size[$k];

    push @$exons, Bio::EnsEMBL::Exon->new(
                    -SLICE     => $slice,
                    -START     => $exon_start + 1,
                    -END       => $exon_end,
                    -STRAND    => $strand,
                    -PHASE     => $exon_frames[$k],
                  );
  }
  $exons = [ reverse @$exons ] if $strand == -1;

  # Create new transcript (the ORF) with previous exons (ORF blocks)
  my $new_tr = Bio::EnsEMBL::Transcript->new(
    -STABLE_ID => $orf->{name},
    -BIOTYPE   => $orf->{transcript_biotype},
    -SLICE     => $slice,
    -START     => $orf->{chromStart},
    -END       => $orf->{chromEnd},
    -STRAND    => $strand,
    -VERSION   => 1,
    -EXONS     => $exons,
  );

  # Create respective translation
  $new_tr->translation(Bio::EnsEMBL::Translation->new(
    -START_EXON => $exons->[0],
    -END_EXON   => $exons->[-1],
    -SEQ_START  => 1,
    -SEQ_END    => length($new_tr->spliced_seq),
    -VERSION    => 1
  ));

  # Check if translation matches the one from annotation
  my $seq_pred  = $new_tr->translation->seq;
  my $seq_annot = $orf->{sequence};
  my $var       = $tva->variation_feature->name;
  die "Unexpected translation sequence for ${var}:\n" .
      "  - Predicted: ${seq_pred}\n" .
      "  - Annotated: ${seq_annot}\n" if $seq_pred ne $seq_annot;

  # Create TranscriptVariation object
  my $tv = Bio::EnsEMBL::Variation::TranscriptVariation->new(
    -transcript        => $new_tr,
    -variation_feature => $tva->variation_feature,
  );

  # Fix start_lost consequence not being returned when the variant starts upstream of the translation
  $tv->translation_start(1) unless defined $tv->translation_start;
  $tv->cdna_start(1) unless defined $tv->cdna_start;
  $tv->cds_start(1) unless defined $tv->cds_start;
  $tv->cdna_start_unshifted(1) unless defined $tv->cdna_start_unshifted;

  # Create new TranscriptVariationAllele object with the current allele info
  my $ref_allele = Bio::EnsEMBL::Variation::TranscriptVariationAllele->new(
    -base_variation_feature_overlap => $tv,
    -variation_feature_seq          => $tva->base_variation_feature->ref_allele_string,
    -is_reference                   => 1,
  );
  $tv->add_TranscriptVariationAllele($ref_allele);

  my $tva_orf = Bio::EnsEMBL::Variation::TranscriptVariationAllele->new(
    -base_variation_feature_overlap => $tv,
    -variation_feature_seq          => $tva->variation_feature_seq,
    -is_reference                   => 0,
  );
  $tv->add_TranscriptVariationAllele($tva_orf);

  my @ocs    = sort {$a->rank <=> $b->rank} @{$tv->get_all_alternate_TranscriptVariationAlleles->[0]->get_all_OverlapConsequences};
  my $impact = $ocs[0]->impact if @ocs;

  my $hash = {
    RiboseqORFs_id               => $new_tr->stable_id,
    RiboseqORFs_consequences     => [ map { $_->SO_term } @ocs ],
    RiboseqORFs_impact           => $impact,
    RiboseqORFs_codons           => $tva_orf->display_codon_allele_string,
    RiboseqORFs_amino_acids      => $tva_orf->pep_allele_string,
    RiboseqORFs_protein_position => format_coords($tv->translation_start, $tv->translation_end),
    RiboseqORFs_cDNA_position    => format_coords($tv->cdna_start, $tv->cdna_end),
    RiboseqORFs_CDS_position     => format_coords($tv->cds_start, $tv->cds_end),
    RiboseqORFs_publications     => [ split(",", $orf->{ref_studies}) ],
  };
  delete $hash->{RiboseqORFs_codons}      unless $hash->{RiboseqORFs_codons};
  delete $hash->{RiboseqORFs_amino_acids} unless $hash->{RiboseqORFs_amino_acids};
  return $hash;
}

sub run {
  my ($self, $tva) = @_;
  my $vf   = $tva->variation_feature;

  my ($vf_start, $vf_end) = ($vf->{start}, $vf->{end});
     ($vf_start, $vf_end) = ($vf_end, $vf_start) if $vf_start > $vf_end;

  my @data = @{$self->get_data($vf->{chr}, $vf_start, $vf_end)};

  for (@data) {
    next unless _transcripts_match($tva, $_->{'all_transcript_ids'});

    my $res = _recreate_TranscriptVariationAllele_with_ORF($tva, $_);

    # Group results for JSON and REST
    if ($self->{config}->{output_format} eq 'json' || $self->{config}->{rest}) {
      $res = { "RiboseqORFs" => $res };
    }
    return $res;
  }
  return {};
}

sub parse_data {
  my ($self, $line) = @_;

  my %output_hash;
  my @header = qw{
    chrom chromStart chromEnd name score strand thickStart thickEnd reserved
    blockCount blockSize chromStarts name2 cdsStartStat cdsEndStat exonFrames
    type geneName geneName2 geneType transcript_biotype sequence
    all_transcript_ids all_gene_ids replicated ref_studies
  };
  @output_hash{ @header } = split /\t/, $line;
  return \%output_hash;
}


sub get_start {
  return $_[1]->{start};
}

sub get_end {
  return $_[1]->{end};
}

1;
