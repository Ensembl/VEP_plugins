=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2019] EMBL-European Bioinformatics Institute

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

 MPRA - massively parallel reporter assays (MPRA) saturation mutagenesis on 21 regulatory elements (11 enhancers, 10 promoters) - Add MPRA data fields to the VEP output

=head1 SYNOPSIS

 mv MPRA.pm ~/.vep/Plugins
 ./vep -i variations.vcf --plugin MPRA,/path/to/MPRA_data.gz,col1,col2

=head1 DESCRIPTION

A VEP plugin that retrieves data for variants from a tabix-indexed MPRA file (1-based file).

Parameters can be set using a key=value system:
file           : required - a tabix indexed file of the MPRA data corresponding to desired assembly. 
pvalue         : p-value threshold  (default: 0.00001)
cols           : colon delimited list of additional data types to be returned from the MPRA data
                (only 'Value', 'P-Value', and 'Element' reported by default)

MPRA data was obtained for 20 disease-associated regulatory elements and one ultraconserved enhancer in different cell lines:
- ten promoters (of ​TERT, LDLR, HBB, HBG, HNF4A, MSMB, PKLR, F9, FOXE1 and GP1BB​) and
- ten enhancers (of ​SORT1, ZRS, BCL11A, IRF4, IRF6, MYC (2x), RET, TCF7L2 ​and ZFAND3​) and
- one ultraconserved enhancer (​UC88​).

Please refer to the MPRA web server and biorxiv manuscript for more information:
https://mpra.gs.washington.edu/satMutMPRA/
https://www.biorxiv.org/content/10.1101/505362v1.full

The Bio::DB::HTS perl library or tabix utility must be installed in your path
to use this plugin. The MPRA data file can be downloaded from
https://mpra.gs.washington.edu/satMutMPRA/.

MPRA data can be downloaded for both GRCh38 and GRCh37 from the web server (https://mpra.gs.washington.edu/satMutMPRA/):
'Download' section, select 'GRCh37' or 'GRCh38' for 'Genome release' and 'Download All Elements'.

The file must be processed and indexed by tabix before use by this plugin.

# GRCh38
> (grep ^"Chr" GRCh38_ALL.tsv; grep -v ^"Chr" GRCh38_ALL.tsv | sort -k1,1 -k2,2n ) | bgzip > MPRA_GRCh38_ALL.gz
> tabix -s 1 -b 2 -e 2 -c C MPRA_GRCh38_ALL.gz

# GRCh37
> (grep ^"Chr" GRCh37_ALL.tsv; grep -v ^"Chr" GRCh37_ALL.tsv | sort -k1,1 -k2,2n ) | bgzip > MPRA_GRCh37_ALL.gz
> tabix -s 1 -b 2 -e 2 -c C MPRA_GRCh37_ALL.gz

Note that in the last command we tell tabix that the header line starts with "Chr";
this may change to the default of "#" in future versions of MPRA.

When running the plugin by default 'Value', 'P-Value', and 'Element'
information is returned e.g.

--plugin MPRA,file=/path/to/MPRA_GRCh38_ALL.gz

You may include all columns with ALL; this fetches a all data per variant
(e.g. Tags, DNA, RNA, Value, P-Value, Element):

--plugin MPRA,file=/path/to/MPRA_GRCh38_ALL.gz,cols=ALL

You may want to select only a specific subset of additional information to be
reported, you can do this by specifying the extra columns as parameters to the plugin e.g.

--plugin MPRA,file=/path/to/MPRA_GRCh38_ALL.gz,cols=Tags:DNA

If a requested column is not found, the error message will report the
complete list of available columns in the MPRA file. For a detailed description
of the available information please refer to the manuscript or online web server.

Tabix also allows the data file to be hosted on a remote server. This plugin is
fully compatible with such a setup - simply use the URL of the remote file:

--plugin MPRA,file=http://my.files.com/mpra.gz

Note that gene sequences referred to in MPRA may be out of sync with
those in the latest release of Ensembl; this may lead to discrepancies with
information retrieved from other sources.

=cut

package MPRA;

use strict;
use warnings;

use Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin);

my @fields_order;

my $out_txt = 1;
my $out_vcf = 0;
my $out_json = 0;
my $char_sep = "|";

# default config
my %DEFAULTS = (
  pvalue => 0.00001,
);


sub new {
  my $class = shift;

  my $self = $class->SUPER::new(@_);

  $self->expand_left(0);
  $self->expand_right(0);

  my $params_hash = $self->params_to_hash();
  $DEFAULTS{$_} = $params_hash->{$_} for keys %$params_hash;

  # get MPRA file
  my $file = $DEFAULTS{file};
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

  my $nDefaultColumns = 3;
  # get required columns
  my $i = 1;
  ## allow specific columns to be retrieved from the VCF
  if($params_hash->{cols}){
    my @cols = split/\:/, $params_hash->{cols};
    foreach my $col (@cols){
      if($col eq 'ALL') {
        $i = $nDefaultColumns + 1;
        $self->{cols} = {map {$_ => $i++}
            grep {!defined($self->{cols}->{$_})} #only the extra columns
            @{$self->{headers}}};
        last; #if ALL is used, then the loop will exit after all existing header elements have been selected
      }
      die "ERROR: Column $col not found in header for file $file. Available columns are:\n".join(",", @{$self->{headers}})."\n" unless grep {$_ eq $col} @{$self->{headers}};

      $self->{cols}->{$col} = $nDefaultColumns + $i;
      $i++;
    }
  }

  #default columns always reported
  $self->{cols}->{'Value'} = 1;
  $self->{cols}->{'P-Value'} = 2;
  $self->{cols}->{'Element'} = 3;

  $i += $nDefaultColumns; #ensure that $i is higher than the number of selected columns

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

  my $header = 'MPRA data for variation in 21 regulatory features. Format: ';
  $header .= join($char_sep, @fields_order );

  return {
    MPRA => $header,
  }
}

sub run {
  my $self = shift;
  my $tva = shift;

  my $vf = $tva->variation_feature;
  my $ref_allele = $vf->ref_allele_string;
  my %alt_alleles = map {$_ => 1} @{$vf->alt_alleles};

  my ($start, $end) = ($vf->{start}, $vf->{end});
  # adjust coords for insertions
  ($start, $end) = ($vf->{end}, $vf->{start}) if ($vf->{start} > $vf->{end});

  my $data = $self->get_data($vf->{chr}, $start, $end);

  return {} unless $data && scalar @$data;

  my @result =();
  my %result_uniq;
  my @result_str = ();

  foreach my $tmp_data(@{$data}) {
    next unless $tmp_data->{'P-Value'} < $DEFAULTS{pvalue};
    next unless $ref_allele eq $tmp_data->{'Ref'} && exists($alt_alleles{$tmp_data->{'Alt'}});

    # get required data
    my %tmp_return =
      map {$_ => $tmp_data->{$_}}
      grep {defined($self->{cols}->{$_})}  # only include selected cols
      keys %$tmp_data;

      # report only unique set of fields
      my $record_line = join(",", values %tmp_return);
      next if defined $result_uniq{$record_line};
      $result_uniq{$record_line} = 1;

      push(@result_str, join($char_sep, @tmp_return{@fields_order}));
      push(@result, \%tmp_return);
  }

  if (scalar @result > 0){
    return {
      MPRA => $self->{config}->{output_format} eq "json" ? \@result : \@result_str
    }
  } else {
    return {};
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
