=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2017] EMBL-European Bioinformatics Institute

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

 ExAC

=head1 SYNOPSIS

 mv ExAC.pm ~/.vep/Plugins
 ./vep -i variations.vcf --plugin ExAC,/path/to/ExAC/ExAC.r0.3.sites.vep.vcf.gz
 ./vep -i variations.vcf --plugin ExAC,/path/to/ExAC/ExAC.r0.3.sites.vep.vcf.gz,AC
 ./vep -i variations.vcf --plugin ExAC,/path/to/ExAC/ExAC.r0.3.sites.vep.vcf.gz,,AN
 ./vep -i variations.vcf --plugin ExAC,/path/to/ExAC/ExAC.r0.3.sites.vep.vcf.gz,AC,AN



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
use Bio::EnsEMBL::Variation::Utils::Sequence qw(get_3prime_seq_offset);

use Bio::EnsEMBL::Variation::Utils::VEP qw(parse_line get_slice);

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
    print STDERR "$remote_test\n";
    # if($remote_test && $remote_test !~ /get_local_version/) {
    #   die "$remote_test\nERROR: Could not find file or index file for remote annotation file $file\n";
    # }
  }

  # check files exist
  else {
    die "ERROR: ExAC file $file not found; you can download it from ftp://ftp.broadinstitute.org/pub/ExAC_release/current\n" unless -e $file;
    die "ERROR: Tabix index file $file\.tbi not found - perhaps you need to create it first?\n" unless -e $file.'.tbi';
  }
  
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

sub run {
  my ($self, $tva) = @_;
  # make sure headers have been loaded
  $self->get_header_info();

  my $vf = $tva->variation_feature;
  
  # get allele, reverse comp if needed
  my $allele;
  
  $allele = $tva->variation_feature_seq;
  reverse_comp(\$allele) if $vf->{strand} < 0;
  
  # adjust coords to account for VCF-like storage of indels
  my ($s, $e) = ($vf->{start} - 1, $vf->{end} + 1);
  
  my $pos_string = sprintf("%s:%i-%i", $vf->{chr}, $s, $e);

  
  # clear cache if it looks like the coords are the same
  # but allele type is different
  delete $self->{cache} if
    defined($self->{cache}->{$pos_string}) &&
    scalar keys %{$self->{cache}->{$pos_string}} &&
    !defined($self->{cache}->{$pos_string}->{$allele});
  
  my $data = {};
  
  # cached?
  if(defined($self->{cache}) && defined($self->{cache}->{$pos_string})) {
    $data = $self->{cache}->{$pos_string};
  }
  
  # read from file
  else {
    open TABIX, sprintf("tabix -f %s %s |", $self->{file}, $pos_string);
    
    while(<TABIX>) {
      chomp;
      s/\r$//g;
      
      # parse VCF line into a VariationFeature object
      my ($vcf_vf) = @{parse_line({format => 'vcf', minimal => 1}, $_)};
      
      # check parsed OK
      next unless $vcf_vf && $vcf_vf->isa('Bio::EnsEMBL::Variation::VariationFeature');
      my @vcf_alleles = split /\//, $vcf_vf->allele_string;
      my $ref_allele  = shift @vcf_alleles;
      my $vcf_vf_start = $vcf_vf->{start};
      my $vcf_vf_end = $vcf_vf->{end};
      my ($input_ref_allele, $input_alt_allele) = split('/', $vf->allele_string);

#      next unless $vcf_vf->{start} == $vf->{start} && $vcf_vf->{end} == $vf->{end};
      my $found_variant_after_shifting = 0;
      if (($vcf_vf_start != $vf->{start} || $vcf_vf_end != $vf->{end})) {
        # look at alleles individually
        # create copy for each ref/alt
        # remove same bases for ref and alt adjust vf coordinates
        # find out if all the alts start with the same base
        # 3prime_align input and ExAC variant to check if they are equivalent
        # deal with ExAV variants that look like this: 2  241696840 . ATCCTCCTCCTCC ATCCTCCTCC,ATCCTCC,ATCC,ATCCTCCTCCTCCTCC,A

        my $seq_to_check;
        my $offset = 0;
        ## sequence to compare is the reference allele for deletion
        $seq_to_check = $input_ref_allele if ($input_alt_allele eq '-') ;
        ## sequence to compare is the alt allele
        $seq_to_check = $input_alt_allele if ($input_ref_allele eq '-');
        if ($seq_to_check) {
          ## 3' flanking sequence to check
          my ($ref_seq, $ref_start, $ref_end) = _get_flank_seq($vf); 
          my $downstream_seq = substr($ref_seq ,$ref_end);
          my $three_prime_allele;
          ($three_prime_allele, $offset ) = get_3prime_seq_offset($seq_to_check, $downstream_seq);
        }
        my $input_vf_start = $vf->{start} + $offset;
        my $input_vf_end = $vf->{end} + $offset;

        foreach my $alt_allele (@vcf_alleles) {
          my ($short_string, $long_string ) = ($ref_allele, $alt_allele);
          if (length($ref_allele) >= length($alt_allele)) {
            ($short_string, $long_string ) = ($alt_allele, $ref_allele);
          }
          my $var_class = ''; 

          my $vf_copy = { %$vcf_vf };
          bless $vf_copy, ref($vcf_vf);

          my $substring = substr($long_string, 0, length($short_string)); 
          my ($new_ref_allele, $new_alt_allele);
          if ($substring eq $short_string) {
            if ($short_string eq $ref_allele) {
              $var_class = 'insertion';
              $new_alt_allele =  substr($long_string, length($short_string));
              $vf_copy->{allele_string} = "-/$new_alt_allele";
            } else {
              $var_class = 'deletion';
              $new_ref_allele = substr($long_string, length($short_string));
              $vf_copy->{allele_string} = "$new_ref_allele/-";
            }
            my $start = $vcf_vf_start + length($short_string);   
            $vf_copy->{start} = $start;
          }

          if ($var_class eq 'insertion' || $var_class eq 'deletion') {
            $vf_copy->{slice} ||= get_slice($self->{config}, $vf_copy->{chr}, undef, 1);
            my ($ref_seq, $ref_start, $ref_end) = _get_flank_seq($vf_copy); 

            my $seq_to_check;
            ## sequence to compare is the reference allele for deletion
            $seq_to_check = $new_ref_allele  if $var_class eq 'deletion' ;

            ## sequence to compare is the alt allele
            $seq_to_check = $new_alt_allele if $var_class eq 'insertion';

            ## 3' flanking sequence to check
            my $downstream_seq = substr($ref_seq ,$ref_end);
            my $three_prime_allele;
            my $offset = 0;

            ($three_prime_allele, $offset ) = get_3prime_seq_offset($seq_to_check, $downstream_seq);

            $vf_copy->{start} = $vf_copy->{start} + $offset;           
            $vf_copy->{end} = $vf_copy->{end} + $offset;           
          }

          if ($vf_copy->{start} == $input_vf_start && $vf_copy->{end} == $input_vf_end) {
            $allele = $alt_allele;
            $found_variant_after_shifting = 1;
            last;
          }
        }       
      } 
      if (!($vcf_vf->{start} == $vf->{start} && $vcf_vf->{end} == $vf->{end})) {
        next if (!$found_variant_after_shifting);
      } 
      
      # iterate over required headers
      HEADER:
      foreach my $h(@{$self->{headers} || []}) {
        my $total_ac = 0;
        
        if(/$h\=([0-9\,]+)/) {
          
          # grab AC
          my @ac = split /\,/, $1;
          next unless scalar @ac == scalar @vcf_alleles;
          
          # now sed header to get AN
          my $anh = $h;
          $anh =~ s/AC/AN/;
          
          my $afh = $h;
          $afh =~ s/AC/AF/;

          # get AC from header
          my $ach = $h;

          if(/$anh\=([0-9\,]+)/) {
            
            # grab AN
            my @an = split /\,/, $1;
            next unless @an;
            my $an;

            foreach my $a(@vcf_alleles) {
              my $ac = shift @ac;
              $an = shift @an if @an;

              $total_ac += $ac;
              if ($self->{display_ac}){
                $data->{$a}->{'ExAC_'.$ach} = $ac;
              }
              if ($self->{display_an}){
                $data->{$a}->{'ExAC_'.$anh} = $an;
              }

              $data->{$a}->{'ExAC_'.$afh} = sprintf("%.3g", $ac / $an) if $an;
            }
            
            # use total to get ref allele freq
            if ($self->{display_ac}){
             $data->{$ref_allele}->{'ExAC_'.$ach} = $total_ac;
            }
            if ($self->{display_an}){
              $data->{$ref_allele}->{'ExAC_'.$anh} = $an;
            }
            $data->{$ref_allele}->{'ExAC_'.$afh} = sprintf("%.3g", 1 - ($total_ac / $an)) if $an;
          }
        }
      }
    }
    
    close TABIX;
  }
  
  # overwrite cache
  $self->{cache} = {$pos_string => $data};
  return defined($data->{$allele}) ? $data->{$allele} : {};
}

sub _get_flank_seq{
  my $vf = shift;
  # Get the underlying slice and sequence
  my $ref_slice = $vf->{slice};
  my $add_length = 100;  ## allow at least 100 for 3'shifting
  my @allele = split(/\//, $vf->allele_string());
  foreach my $al (@allele) { ## alleles be longer
    if(length($al) > $add_length){
      $add_length = length $al ;
    }
  }
  my $seq_start =  $vf->start() - $add_length;
  my $seq_end   =  $vf->end() + $add_length;

  ## variant position relative to flank
  my $ref_start = $add_length;
  my $ref_end   = $add_length + $vf->end() - $vf->start();

  # Should we be at the beginning of the sequence, adjust the coordinates to not cause an exception
  if ($seq_start < 0) {
    $ref_start += $seq_start;
    $ref_end   += $seq_start;
    $seq_start  = 0;
  }
  my $flank_seq = $ref_slice->subseq($seq_start + 1, $seq_end, 1);
  return ($flank_seq, $ref_start, $ref_end );
}

1;
