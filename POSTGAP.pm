=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     https://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

Please email comments or questions to the public Ensembl
developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

Questions may also be sent to the Ensembl help desk at
<https://www.ensembl.org/Help/Contact>.
    
=cut

=head1 NAME

 POSTGAP - Add POSTGAP data fields to the VEP output

=head1 SYNOPSIS

 mv POSTGAP.pm ~/.vep/Plugins
 ./vep -i variations.vcf --plugin POSTGAP,/path/to/PostGap.gz,col1,col2

=head1 DESCRIPTION

A VEP plugin that retrieves data for variants from a tabix-indexed POSTGAP file.

Please refer to the POSTGAP github for more information:
https://github.com/Ensembl/postgap

The Bio::DB::HTS perl library or tabix utility must be installed in your path
to use this plugin. The POSTGAP data file can be downloaded from
https://storage.googleapis.com/postgap-data.

The file must be processed and indexed by tabix before use by this plugin.
POSTGAP has coordinates for both GRCh38 and GRCh37; the file must be
processed differently according to the assembly you use.

> wget https://storage.googleapis.com/postgap-data/postgap.txt.gz
> gunzip postgap.txt.gz

# GRCh38
> (grep ^"ld_snp_rsID" postgap.txt; grep -v ^"ld_snp_rsID" postgap.txt | sort -k4,4 -k5,5n ) | bgzip > postgap_GRCh38.txt.gz
> tabix -s 4 -b 5 -e 5 -c l postgap_GRCh38.txt.gz

# GRCh37
> (grep ^"ld_snp_rsID" postgap.txt; grep -v ^"ld_snp_rsID" postgap.txt | sort -k2,2 -k3,3n ) | bgzip > postgap_GRCh37.txt.gz
> tabix -s 2 -b 3 -e 3 -c l postgap_GRCh37.txt.gz

Note that in the last command we tell tabix that the header line starts with "l";
this may change to the default of "#" in future versions of POSTGAP.

Tabix also allows the data file to be hosted on a remote server. This plugin is
fully compatible with such a setup - simply use the URL of the remote file:

--plugin POSTGAP,http://my.files.com/postgap.txt.gz

Note that gene sequences referred to in POSTGAP may be out of sync with
those in the latest release of Ensembl; this may lead to discrepancies with
scores retrieved from other sources.

=cut

package POSTGAP;

use strict;
use warnings;

use DBI;
use POSIX;

use Bio::EnsEMBL::VEP::Utils qw(convert_arrayref);
use Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin);

my @fields_order;

my $out_txt = 1;
my $out_vcf = 0;
my $out_json = 0;
my $char_sep = "|";

sub new {
  my $class = shift;

  my $self = $class->SUPER::new(@_);

  $self->expand_left(0);
  $self->expand_right(0);

  my $params_hash = $self->params_to_hash();

  # get POSTGAP file
  my $file = $self->params->[0];
  $self->add_file($file);

  # get output format
  $out_vcf  = 1 if ($self->{config}->{output_format} eq "vcf");
  $out_json = 1 if ($self->{config}->{output_format} eq "json");
  $out_txt = 0 if ($out_vcf || $out_json);

  $char_sep = "+" if $out_vcf;

  # get headers
  open HEAD, "tabix -fh $file 1:1-1 2>&1 | ";
  while(<HEAD>) {
    chomp;
    $self->{headers} = [split];
  }
  close HEAD;
  die "ERROR: Could not read headers from $file\n" unless defined($self->{headers}) && scalar @{$self->{headers}};

  # get required columns
  my $all =0;
  my $i = 4;
  while(defined($self->params->[$i])) {
    my $col = $self->params->[$i];
    if($col eq 'ALL') {
      $self->{cols} = {map {$_ => $i++}
          grep {!defined($self->{cols}->{$_})} #only the extra columns
          @{$self->{headers}}};
      last;
    }
    die "ERROR: Column $col not found in header for file $file. Available columns are:\n".join(",", @{$self->{headers}})."\n" unless grep {$_ eq $col} @{$self->{headers}};

    $self->{cols}->{$self->params->[$i]} = $i;
    $i++;
  }

  #default columns always reported
  $self->{cols}->{'disease_efo_id'} = 1;
  $self->{cols}->{'disease_name'} = 2;
  $self->{cols}->{'score'} = 3;

  # get the order of the output fields into an array, $i is the total number of columns +1
  @fields_order = map { $_ }
  sort {
    (defined($self->{cols}->{$a}) ? $self->{cols}->{$a} : $i)
    <=>
    (defined($self->{cols}->{$b}) ? $self->{cols}->{$b} : $i)
    ||
    $a cmp $b
  }
  keys %{$self->{cols}};

  return $self;
}

sub feature_types {
  return ['Feature','Intergenic'];
}

sub get_header_info {
  my $self = shift;
  
  my $header = '"POSTGAP data for variation - phenotype association. Format: Allele';
  $header .= $char_sep.join($char_sep, @fields_order ).'">';

  return { 
    POSTGAP => $header,
  }
}

sub run {
  my $self = shift;
  my $tva = shift;

  my $vf = $tva->variation_feature;

  # adjust coords for tabix
  my ($start, $end) = ($vf->{start} - 1, $vf->{end});

  my $data = $self->get_data($vf->{chr}, $start, $end);

  return {} unless $data && scalar @$data;

  my @result =();
  my %result_uniq;
  my @result_str = ();

  foreach my $tmp_data(@{$data}) {
    # get required data
    my %tmp_return =
      map {$_ => $tmp_data->{$_}}
      grep {defined($self->{cols}->{$_})}  # only include selected cols
      keys %$tmp_data;

      # get only 3 significant digits
      $tmp_return{score} = sprintf("%.3f", $tmp_return{score}) if defined $tmp_return{score};

      # report only unique set of fields
      my $record_line = join(",", values %tmp_return);
      next if defined $result_uniq{$record_line};
      $result_uniq{$record_line} = 1;

      push(@result_str, $vf->allele_string.$char_sep.join($char_sep, @tmp_return{@fields_order}));
      push(@result, \%tmp_return);
  }
  
  return {
    POSTGAP => $self->{config}->{output_format} eq "json" ? \@result : \@result_str
  }
}

sub parse_data {
  my ($self, $line) = @_;

  my @split = split /\t/, $line;

  # parse data into hash of col names and values
  my %data = map {$self->{headers}->[$_] => $split[$_]} (0..(scalar @{$self->{headers}} - 1));

  return \%data;
}

1;
