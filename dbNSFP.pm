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

 dbNSFP

=head1 SYNOPSIS

 mv dbNSFP.pm ~/.vep/Plugins
 perl variant_effect_predictor.pl -i variations.vcf --plugin dbNSFP,/path/to/dbNSFP.gz,col1,col2

=head1 DESCRIPTION

 A VEP plugin that retrieves data for missense variants from a tabix-indexed
 dbNSFP file.
 
 Please cite the dbNSFP publication alongside the VEP if you use this resource:
 http://www.ncbi.nlm.nih.gov/pubmed/21520341
 
 The tabix utility must be installed in your path to use this plugin. The dbNSFP
 data file can be downloaded from
 https://sites.google.com/site/jpopgen/dbNSFP.
 
 The file must be processed and indexed by tabix before use by this plugin:
 
 > wget ftp://dbnsfp:dbnsfp@dbnsfp.softgenetics.com/dbNSFPv3.2a.zip
 > unzip dbNSFPv3.2a.zip
 > head -n1 dbNSFP3.2a_variant.chr1 > h
 > cat dbNSFP3.2a_variant.chr* | grep -v ^#chr | sort -k1,1 -k2,2n - | cat h - | bgzip -c > dbNSFP.gz
 > tabix -s 1 -b 2 -e 2 dbNSFP.gz
 
 When running the plugin you must list at least one column to retrieve from the
 dbNSFP file, specified as parameters to the plugin e.g.
 
 --plugin dbNSFP,/path/to/dbNSFP.gz,LRT_score,GERP++_RS
 
 Tabix also allows the data file to be hosted on a remote server. This plugin is
 fully compatible with such a setup - simply use the URL of the remote file:
 
 --plugin dbNSFP,http://my.files.com/dbNSFP.gz,col1,col2
 
 Note that transcript sequences referred to in dbNSFP may be out of sync with
 those in the latest release of Ensembl; this may lead to discrepancies with
 scores retrieved from other sources.

 If the dbNSFP README file is found in the same directory as the data file,
 column descriptions will be read from this and incorporated into the VEP output
 file header.
 
=cut

package dbNSFP;

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
  
  # get dbNSFP file
  my $file = $self->params->[0];
  $self->add_file($file);
  
  # get headers
  open HEAD, "tabix -fh $file 1:1-1 2>&1 | ";
  while(<HEAD>) {
    next unless /^\#/;
    chomp;
    $self->{headers} = [split];
  }
  close HEAD;
  
  die "ERROR: Could not read headers from $file\n" unless defined($self->{headers}) && scalar @{$self->{headers}};
  
  # check alt and Ensembl_transcriptid headers
  foreach my $h(qw(alt Ensembl_transcriptid)) {
    die "ERROR: Could not find required column $h in $file\n" unless grep {$_ eq $h} @{$self->{headers}};
  }
  
  # get required columns
  my $i = 1;
  while(defined($self->params->[$i])) {
    my $col = $self->params->[$i];
    die "ERROR: Column $col not found in header for file $file. Available columns are:\n".join(",", @{$self->{headers}})."\n" unless grep {$_ eq $col} @{$self->{headers}};
    
    $self->{cols}->{$self->params->[$i]} = 1;
    $i++;
  }
  
  die "ERROR: No columns selected to fetch. Available columns are:\n".join(",", @{$self->{headers}})."\n" unless defined($self->{cols}) && scalar keys %{$self->{cols}};
  
  return $self;
}

sub feature_types {
  return ['Transcript'];
}

sub get_header_info {
  my $self = shift;

  if(!exists($self->{_header_info})) {

    # look for readme
    my $file_dir = $self->files->[0];

    my %rm_descs;

    # won't work for remote
    if($file_dir !~ /tp\:\/\//) {

      # get just dir
      $file_dir =~ s/\/[^\/]+$/\//;

      if(opendir DIR, $file_dir) {
        my ($readme_file) = grep {/dbnsfp.*readme\.txt/i} readdir DIR;
        closedir DIR;

        if(open RM, $file_dir.$readme_file) {
          my ($col, $reading);

          # parse dbNSFP readme
          # relevant lines look like:
          #
          # 1   column1_name: description blah blah
          #     blah blah blah
          # 2   column2_name: description blah blah
          #     blah blah blah

          while(<RM>) {
            chomp;
            s/\r$//g;

            if(/^\d+\s/) {
              $reading = 1;

              m/^\d+\s+(.+?)\:\s+(.+)/;
              $col = $1;

              $rm_descs{$col} = '(from dbNSFP) '.$2;
            }
            elsif($reading && /\w/) {
              s/^\s+//;
              $rm_descs{$col} .= ' '.$_;
            }
            else {
              $reading = 0;
            }
          }

          close RM;

          # remove multiple spaces
          $rm_descs{$_} =~ s/\s+/ /g for keys %rm_descs;
        }
      }
    }

    $self->{_header_info} = {map {$_ => $rm_descs{$_} || ($_.' from dbNSFP file '.$self->{file})} keys %{$self->{cols}}};
  }
  
  return $self->{_header_info};
}

sub run {
  my ($self, $tva) = @_;
  
  # only for missense variants
  return {} unless grep {$_->SO_term eq 'missense_variant'} @{$tva->get_all_OverlapConsequences};
  
  my $vf = $tva->variation_feature;
  
  return {} unless $vf->{start} eq $vf->{end};
  
  # get allele, reverse comp if needed
  my $allele = $tva->variation_feature_seq;
  reverse_comp(\$allele) if $vf->{strand} < 0;
  
  return {} unless $allele =~ /^[ACGT]$/;
  
  # get transcript stable ID
  my $tr_id = $tva->transcript->stable_id;

  my $data;
  my $pos;
  
  foreach my $tmp_data(@{$self->get_data($vf->{chr}, $vf->{start} - 1, $vf->{end})}) {

    # compare allele and transcript
    if (exists $tmp_data->{'pos(1-based)'}) {
      $pos = $tmp_data->{'pos(1-based)'}
    } else {
      # for dbNSFP version 2.9.1
      $pos = $tmp_data->{'pos(1-coor)'}
    }
    next unless
      ${pos} == $vf->{start} &&
      defined($tmp_data->{alt}) &&
      $tmp_data->{alt} eq $allele;
    
    # make a clean copy as we're going to edit it
    %$data = %$tmp_data;

    # convert data with multiple transcript values
    # if($data->{Ensembl_transcriptid} =~ m/\;/) {

    #   # find the "index" of this transcript
    #   my @tr_ids = split(';', $data->{Ensembl_transcriptid});
    #   my $tr_index;

    #   for my $i(0..$#tr_ids) {
    #     $tr_index = $i;
    #     last if $tr_ids[$tr_index] =~ /^$tr_id(\.\d+)?$/;
    #   }

    #   next unless defined($tr_index);

    #   # now alter other fields
    #   foreach my $key(keys %$data) {
    #     if($data->{$key} =~ m/\;/) {
    #       my @split = split(';', $data->{$key});
    #       die("ERROR: Transcript index out of range") if $tr_index > $#split;
    #       $data->{$key} = $split[$tr_index];
    #     } 
    #   }
    # }
    last;
  }
  
  return {} unless scalar keys %$data;
  
  # get required data
  my %return =
    map {$_ => ($data->{$_} =~ s/[|]/&/gr) } # replace | with & to prevent conflict with existing field sep
    map {$data->{$_} =~ s/\;/\,/g; $_ }
    grep {$data->{$_} ne '.'}              # ignore missing data
    grep {defined($self->{cols}->{$_})}  # only include selected cols
    keys %$data;
  
  return \%return;
}

sub parse_data {
  my ($self, $line) = @_;

  $line =~ s/\r$//g;

  my @split = split /\t/, $line;
  
  # parse data into hash of col names and values
  my %data = map {$self->{headers}->[$_] => $split[$_]} (0..(scalar @{$self->{headers}} - 1));

  return \%data;
}

sub get_start {  
  return $_[1]->{'pos(1-based)'};
}

sub get_end {
  return $_[1]->{'pos(1-based)'};
}

1;

