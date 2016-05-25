=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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
 
 > wget ftp://dbnsfp:dbnsfp@dbnsfp.softgenetics.com/dbNSFPv3.0b2a.zip
 > unzip dbNSFPv3.0b2a.zip
 > cat dbNSFP*chr* | bgzip -c > dbNSFP.gz
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

use Bio::EnsEMBL::Variation::Utils::BaseVepPlugin;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

sub new {
  my $class = shift;
  
  my $self = $class->SUPER::new(@_);
  
  # test tabix
  die "ERROR: tabix does not seem to be in your path\n" unless `which tabix 2>&1` =~ /tabix$/;
  
  # get dbNSFP file
  my $file = $self->params->[0];
  
  # remote files?
  if($file =~ /tp\:\/\//) {
    my $remote_test = `tabix -f $file 1:1-1 2>&1`;
    if($remote_test && $remote_test !~ /get_local_version/) {
      die "$remote_test\nERROR: Could not find file or index file for remote annotation file $file\n";
    }
  }

  # check files exist
  else {
    die "ERROR: dbNSFP file $file not found\n" unless -e $file;
    die "ERROR: Tabix index file $file\.tbi not found - perhaps you need to create it first?\n" unless -e $file.'.tbi';
  }
  
  $self->{file} = $file;
  
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

sub version {
  return 71;
}

sub feature_types {
  return ['Transcript'];
}

sub get_header_info {
  my $self = shift;

  if(!exists($self->{_header_info})) {

    # look for readme
    my $file_dir = $self->{file};

    my %rm_descs;

    # won't work for remote
    if($file_dir !~ /tp\:\/\//) {

      # get just dir
      $file_dir =~ s/\/[^\/]+$/\//;

      if(opendir DIR, $file_dir) {
        my ($readme_file) = grep {/dbnsfp.*readme/i} readdir DIR;
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
    
  my $pos_string = sprintf("%s:%i-%i", $vf->{chr}, $vf->{start}, $vf->{end});
  
  my @dbnsfp_data;
  
  # cached?
  if(defined($self->{cache}) && defined($self->{cache}->{$pos_string})) {
    @dbnsfp_data = @{$self->{cache}->{$pos_string}};
  }
  
  # read from file
  else {
    open TABIX, sprintf("tabix -f %s %s |", $self->{file}, $pos_string);

    while(<TABIX>) {
      chomp;
      s/\r$//g;
      my @split = split /\t/;
      
      # parse data into hash of col names and values
      my %data = map {$self->{headers}->[$_] => $split[$_]} (0..(scalar @{$self->{headers}} - 1));
      
      push @dbnsfp_data, \%data;
    }
    
    close TABIX;
  }
  
  # overwrite cache
  $self->{cache} = {$pos_string => \@dbnsfp_data};
  
  my $data;
  
  foreach my $tmp_data(@dbnsfp_data) {
    # compare allele and transcript
    next unless
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

1;

