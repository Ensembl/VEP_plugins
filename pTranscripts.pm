=head1 LICENSE

GNU General Public License, version 3 (GPL-3.0)

=head1 CONTACT

 Questions may also be sent to <duarte.molha@ogt.co.uk>.

=cut

=head1 NAME

pTranscripts

=head1 SYNOPSIS

    mv pTranscripts.pm ~/.vep/Plugins
    perl variant_effect_predictor.pl -i variations.vcf --plugin pTranscripts
 
 
=head1 DESCRIPTION

 A VEP plugin that checks to see if an intergenic variant 
 has overlaps any prediction Transcripts. Adds that feature type 
 to the file for intergenic variants if true.

=cut

package pTranscripts;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

use strict;
use warnings;


sub version {
    return '2.6';
}


sub get_header_info {
    return {
        pTranscript => "Checks for overlaping prediction transcripts in Intergenic Variants"
    };
}


sub feature_types {
   return ['Intergenic'];
}


sub run {
    my ($self, $bvfoa, $line_hash) = @_;
    
    my $bvf = $bvfoa->{base_variation_feature_overlap}->{base_variation_feature};
    #plugin to run on intergenic variations
	
	################################### pTrans ######################################################
	my $slice = $self->{config}->{sa}->fetch_by_region( 'chromosome', $bvf->{chr}, $bvf->{start}, $bvf->{end}, $bvf->{strand} );
    
	if ($slice){
		my $ptrans = $slice->get_all_PredictionTranscripts();
		my @PT = @{ $ptrans };
		my $pt = ""; 
		my @ptrans_data = ();
		
		if (@PT){
			foreach my $pTrans (@PT){
				my $pt_ID = $pTrans->stable_id();
				my $type  = $pTrans->analysis()->logic_name();
				my @exons = $pTrans->get_all_Exons();
                push @ptrans_data, $pt_ID."(".scalar (@exons).")"; 
			}
            $line_hash->{Feature_type} = 'pTranscript';
            $line_hash->{Feature} = join "|", @ptrans_data;
		}
	}
    
    return {};
}

1;
