=head1 LICENSE

 Copyright (c) 1999-2013 The European Bioinformatics Institute and                                                   
 Genome Research Limited.  All rights reserved.                                                                      

 This software is distributed under a modified Apache license.                                                       
 For license details, please see

   http://www.ensembl.org/info/about/code_licence.html                                                               

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
 
 dbNSFP publication: http://www.ncbi.nlm.nih.gov/pubmed/21520341
 
 The tabix utility must be installed in your path to use this plugin. The dbNSFP
 data file can be downloaded from
 http://dbnsfp.houstonbioinformatics.org/dbNSFPzip/.
 
 The file must be processed and indexed by tabix before use by this plugin:
 
 > wget http://dbnsfp.houstonbioinformatics.org/dbNSFPzip/dbNSFP2.0.zip
 > unzip dbNSFP2.0.zip
 > cat dbNSFP2.0_variant.chr* | bgzip -c > dbNSFP.gz
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
    my $remote_test = `tabix $file 1:1-1 2>&1`;
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
  open HEAD, "tabix -h $file 1:1-1 2>&1 | ";
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
  
  my %tmp = map {$_ => $_.' from dbNSFP file '.$self->{file}} keys %{$self->{cols}};
  
  return \%tmp;
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
    open TABIX, sprintf("tabix %s %s |", $self->{file}, $pos_string);
    
    while(<TABIX>) {
      chomp;
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
      $tmp_data->{alt} eq $allele &&
      defined($tmp_data->{Ensembl_transcriptid}) &&
      $tmp_data->{Ensembl_transcriptid} =~ /$tr_id($|;)/;
    
    $data = $tmp_data;
    last;
  }
  
  return {} unless scalar keys %$data;
  
  # get required data
  my %return =
    map {$_ => $data->{$_}}
    grep {$data->{$_} ne '.'}              # ignore missing data
    grep {defined($self->{cols}->{$_})}  # only include selected cols
    keys %$data;
  
  return \%return;
}

1;

