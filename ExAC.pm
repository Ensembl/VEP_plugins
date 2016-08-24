=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016] EMBL-European Bioinformatics Institute

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

 Will McLaren <wm2@ebi.ac.uk>

=cut

=head1 NAME

 ExAC

=head1 SYNOPSIS

 mv ExAC.pm ~/.vep/Plugins
 perl variant_effect_predictor.pl -i variations.vcf --plugin ExAC,/path/to/ExAC/ExAC.r0.3.sites.vep.vcf.gz
 perl variant_effect_predictor.pl -i variations.vcf --plugin ExAC,/path/to/ExAC/ExAC.r0.3.sites.vep.vcf.gz,AC
 perl variant_effect_predictor.pl -i variations.vcf --plugin ExAC,/path/to/ExAC/ExAC.r0.3.sites.vep.vcf.gz,,AN
 perl variant_effect_predictor.pl -i variations.vcf --plugin ExAC,/path/to/ExAC/ExAC.r0.3.sites.vep.vcf.gz,AC,AN



=head1 DESCRIPTION

 A VEP plugin that retrieves ExAC allele frequencies.

 Visit ftp://ftp.broadinstitute.org/pub/ExAC_release/current to download the latest ExAC VCF.

 Note that the currently available version of the ExAC data file (0.3) is only available
 on the GRCh37 assembly; therefore it can only be used with this plugin when using the
 VEP on GRCh37. See http://www.ensembl.org/info/docs/tools/vep/script/vep_other.html#assembly

 The tabix utility must be installed in your path to use this plugin.

 The plugin takes 3 command line arguments. Second and third arguments are not mandatory. If AC specified as second
 argument Allele counts per population will be included in output. If AN specified as third argument Allele specific
 chromosome counts will be included in output.


=cut

package ExAC;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use Bio::EnsEMBL::Variation::Utils::VEP qw(parse_line);

use Bio::EnsEMBL::Variation::Utils::BaseVepPlugin;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

sub new {
  my $class = shift;

  my $self = $class->SUPER::new(@_);

  # test tabix
  die "ERROR: tabix does not seem to be in your path\n" unless `which tabix 2>&1` =~ /tabix$/;

  # get ExAC file
  my $file = $self->params->[0];

  # get AC,AN options
  if (exists($self->params->[1]) && $self->params->[1] eq 'AC'){
    $self->{display_ac} = 1;
  }
  else {
    $self->{display_ac} = 0;
  }

  if (exists($self->params->[2]) && $self->params->[2] eq 'AN'){
    $self->{display_an} = 1;
  }
  else {
    $self->{display_an} = 0;
  }

  # remote files?
  if($file =~ /tp\:\/\//) {
    my $remote_test = `tabix -f $file 1:1-1 2>&1`;
    if($remote_test && $remote_test !~ /get_local_version/) {
      die "$remote_test\nERROR: Could not find file or index file for remote annotation file $file\n";
    }
  }

  # check files exist
  else {
    die "ERROR: ExAC file $file not found; you can download it from ftp://ftp.broadinstitute.org/pub/ExAC_release/current\n" unless -e $file;
    die "ERROR: Tabix index file $file\.tbi not found - perhaps you need to create it first?\n" unless -e $file.'.tbi';
  }

  # initialize current iteration of ExAC file
  $self->{exac_state} = {chr => "", pos => 0};

  $self->{file} = $file;

  return $self;
}

sub feature_types {
  return ['Feature','Intergenic'];
}

sub get_header_info {
  my $self = shift;

  if(!exists($self->{header_info})) {
    open IN, "tabix -f -h ".$self->{file}." 1:1-1 |";

    my %headers = ();
    my @lines = <IN>;

    while(my $line = shift @lines) {
      if($line =~ /ID\=AC(\_[A-Zdj]+)?\,.*\"(.+)\"/) {
        my ($pop, $desc) = ($1, $2);

        $desc =~ s/Counts?/frequency/i;
        $pop ||= '';

        my $field_name = 'ExAC_AF'.$pop;
        $headers{$field_name} = 'ExAC '.$desc;

        if ($self->{display_ac}){
          $field_name = 'ExAC_AC'.$pop;
          $headers{$field_name} = 'ExAC'.$pop.' Allele count';
        }
        if ($self->{display_an}){
          $field_name = 'ExAC_AN'.$pop;
          $headers{$field_name} = 'ExAC'.$pop.' Allele number';
        }

        # store this header on self
        push @{$self->{headers}}, 'AC'.$pop;
      }
    }

    close IN;

    die "ERROR: No valid headers found in ExAC VCF file\n" unless scalar keys %headers;

    $self->{header_info} = \%headers;
  }

  return $self->{header_info};
}


# given a pair of alleles, return the shortest Ensembl-like representation and
# adjust the starting coordindate if necessary.
# operates on references to the alleles and position!
sub fix_alleles {
  my ($ref, $alt, $start) = @_;

  # if the first base is the same, trim it and adjust starting coord
  if (substr($$ref, 0, 1) eq substr($$alt, 0, 1)) {
    $$ref = substr($$ref, 1) || "-";
    $$alt = substr($$alt, 1) || "-";
    $$start += 1;
  }

  # remove any identical sequence of bases from the end of the alleles
  while (substr($$ref, -1) eq substr($$alt, -1)) {
    $$ref = substr($$ref, 0, -1) || "-";
    $$alt = substr($$alt, 0, -1) || "-";
  }
}


# read ExAC file up to the position of the current variant feature and store
# the parsed ExAC data in the cache.
sub process_exac {
  my ($self, $vf) = @_;

  # if we're on a new chromosome or at a position we've already passed, open up a new file descriptor
  if ($vf->{chr} ne $self->{exac_state}->{chr} || $vf->{start} < $self->{exac_state}->{last_vf_start}) {

    close($self->{exac_state}->{fp}) if defined $self->{exac_state}->{fp};

    # we can skip directly to the first requested position (helpful for parallelization)
    open $self->{exac_state}->{fp}, sprintf("tabix -f %s %s:%s |", $self->{file}, $vf->{chr}, $vf->{start} - 1);

    # update state and empty cache
    $self->{exac_state}->{chr} = $vf->{chr};
    $self->{exac_state}->{pos} = 0;
    $self->{cache} = [];
  }

  $self->{exac_state}->{last_vf_start} = $vf->{start};

  # keep only cache entries that are at or ahead of the current variant feature
  # subtract one from starting position to account for indels
  # (a A/- Ensembl variant at position 10 might be a CA/C variant at position 9
  # in ExAC)
  $self->{cache} = [ grep { $_->{pos} >= $vf->{start} - 1 } @{$self->{cache}} ];

  # main ExAC parsing loop
  # iterate over ExAC file until we reach the position of the current variation
  # feature or we run out of ExAC variants
  while ($self->{exac_state}->{pos} < $vf->{start} && !eof($self->{exac_state}->{fp})) {

    $_ = readline $self->{exac_state}->{fp};
    chomp;
    s/\r$//;

    my @fields = split /\s+/;
    my $chr = $fields[0];
    my $start = $fields[1];
    my $ref = $fields[3];
    my @alts = split /,/, $fields[4];
    my $info = $fields[7];

    # if the variant is not at the position of the variation feature (accounting
    # for indels), don't parse it further
    next unless $start >= $vf->{start} - 1;

    $self->{exac_state}->{chr} = $chr;
    $self->{exac_state}->{pos} = $start;

    # map of unfixed ExAC alternate allele to hash of frequency data
    my %data = map { $_ => {} } @alts;

    # iterate over required headers
    foreach my $h(@{$self->{headers} || []}) {
      my $total_ac = 0;

      if ($info =~ /$h=([0-9,]+)/) {

        # grab AC
        my @ac = split /,/, $1;
        next unless scalar @ac == scalar @alts;

        # now sed header to get AN
        my $anh = $h;
        $anh =~ s/AC/AN/;

        my $afh = $h;
        $afh =~ s/AC/AF/;

        # get AC from header
        my $ach = $h;

        if ($info =~ /$anh=([0-9,]+)/) {

          # grab AN
          my @an = split /,/, $1;
          next unless @an;
          my $an;

          foreach my $a(@alts) {

            my $ac = shift @ac;
            $an = shift @an if @an;

            # no dividing by 0
            next unless $an;

            $total_ac += $ac;
            if ($self->{display_ac}) {
              $data{$a}->{'ExAC_'.$ach} = $ac;
            }
            if ($self->{display_an}) {
              $data{$a}->{'ExAC_'.$anh} = $an;
            }

            $data{$a}->{'ExAC_'.$afh} = sprintf("%.3g", $ac / $an);
          }
        }
      }
    }

    # store data in cache using fixed alleles (fixed for each ref/alt pair)
    foreach my $a(@alts) {
      my $a_ref = $ref;
      my $a_alt = $a;
      my $a_pos = $start;
      fix_alleles(\$a_ref, \$a_alt, \$a_pos);

      push @{$self->{cache}}, {chr => $chr, pos => $a_pos, ref => $a_ref, alt => $a_alt, data => $data{$a}};
    }
  }
}

sub run {
  my ($self, $tva) = @_;


  # make sure headers have been loaded
  $self->get_header_info();

  my $vf = $tva->variation_feature;

  # get allele, reverse comp if needed
  my $allele;

  $allele = $tva->variation_feature_seq;
  reverse_comp(\$allele) if $vf->{strand} < 0;

  # continue iteration of ExAC until we get to our variation feature position
  $self->process_exac($vf);

  my @alleles = split /\//, $vf->allele_string;
  my $ref = shift @alleles;

  my $data = {};

  # search cache for data corresponding to fixed alleles (for each ref/alt pair)
  for my $a(@alleles) {
    my $a_ref = $ref;
    my $a_alt = $a;
    my $a_pos = $vf->{start};
    fix_alleles(\$a_ref, \$a_alt, \$a_pos);

    # find cache hits
    my @hits = grep {
      $_->{chr} eq $vf->{chr} &&
      $_->{pos} == $a_pos &&
      $_->{ref} eq $a_ref &&
      $_->{alt} eq $a_alt
    } @{$self->{cache}};

    # there should never be more than one
    if (scalar @hits) {
      my $hit = shift @hits;

      # store it using the unfixed Ensembl allele, used for lookup at the end
      $data->{$a} = $hit->{data};
    }
  }

  return defined($data->{$allele}) ? $data->{$allele} : {};
}

1;

