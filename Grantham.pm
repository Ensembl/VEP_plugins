=head1 LICENSE
                                                                                                    GNU General Public License, version 3 (GPL-3.0)
  
Based on plugin of Blosom62.pm done my Graham Ritchie <grsr@ebi.ac.uk> at EBI
                                                                                                                                                                   
=head1 CONTACT                                                                                                       
 Duarte Molha <duartemolha@gmail.com>
    
=cut

=head1 NAME

 Grantham

=head1 SYNOPSIS

 mv OGT_Grantham.pm ~/.vep/Plugins
 perl variant_effect_predictor.pl -i variations.vcf --plugin OGT_Grantham

=head1 DESCRIPTION

 This is a plugin for the Ensembl Variant Effect Predictor (VEP) that
 looks up the GRANTHAM substitution matrix score for the reference
 and alternative amino acids predicted for a missense mutation. It adds
 one new entry to the VEP's Extra column, Grantham which is the 
 associated score. 
 Grantham Matrix score derived from this paper: 
 Grantham, R. Amino Acid Difference Formula to Help Explain Protein Evolution, Science 1974 Sep 6;185(4154):862-4.

=cut

package Grantham;
use strict;
use warnings;
use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

my @grantham_matrix = qw(
  0 112 111 126 195  91 107  60  86  94  96 106  84 113  27  99  58 148 112  64
112   0  86  96 180  43  54 125  29  97 102  26  91  97 103 110  71 101  77  96
111  86   0  23 139  46  42  80  68 149 153  94 142 158  91  46  65 174 143 133
126  96  23   0 154  61  45  94  81 168 172 101 160 177 108  65  85 181 160 152
195 180 139 154   0 154 170 159 174 198 198 202 196 205 169 112 149 215 194 192
 91  43  46  61 154   0  29  87  24 109 113  53 101 116  76  68  42 130  99  96
107  54  42  45 170  29   0  98  40 134 138  56 126 140  93  80  65 152 122 121
 60 125  80  94 159  87  98   0  98 135 138 127 127 153  42  56  59 184 147 109
 86  29  68  81 174  24  40  98   0  94  99  32  87 100  77  89  47 115  83  84
 94  97 149 168 198 109 134 135  94   0   5 102  10  21  95 142  89  61  33  29
 96 102 153 172 198 113 138 138  99   5   0 107  15  22  98 145  92  61  36  32
106  26  94 101 202  53  56 127  32 102 107   0  95 102 103 121  78 110  85  97
 84  91 142 160 196 101 126 127  87  10  15  95   0  28  87 135  81  67  36  21
113  97 158 177 205 116 140 153 100  21  22 102  28   0 114 155 103  40  22  50
 27 103  91 108 169  76  93  42  77  95  98 103  87 114   0  74  38 147 110  68
 99 110  46  65 112  68  80  56  89 142 145 121 135 155  74   0  58 177 144 124
 58  71  65  85 149  42  65  59  47  89  92  78  81 103  38  58   0 128  92  69
148 101 174 181 215 130 152 184 115  61  61 110  67  40 147 177 128   0  37  88
112  77 143 160 194  99 122 147  83  33  36  85  36  22 110 144  92  37   0  55
 64  96 133 152 192  96 121 109  84  29  32  97  21  50  68 124  69  88  55   0
);

my @AAs = qw(A R N D C Q E G H I L K M F P S T W Y V);

sub new {
    my $class = shift;

    my $self = $class->SUPER::new(@_);

    # construct a hash representing the matrix for quick lookups

    my $num = @AAs;

    for (my $i = 0; $i < $num; $i++) {
        for (my $j = 0; $j < $num; $j++) {
            $self->{matrix}->{$AAs[$i]}->{$AAs[$j]} = $grantham_matrix[($i * $num) + $j];
        }
    }

    return $self;
}

sub version {
    return '2.3';
}
sub feature_types {
    return ['Transcript'];
}

sub get_header_info {
    return {
        Grantham => "Grantham Matrix score - Grantham, R. Amino Acid Difference Formula to Help Explain Protein Evolution, Science 1974 Sep 6;185(4154):862-4.",
    };
}

sub run {
    my ($self, $tva) = @_;

    if ($tva->pep_allele_string && $tva->pep_allele_string =~ /^([A-Z])\/([A-Z])$/) {
        
        my $score = $self->{matrix}->{$1}->{$2};
        
        if (defined $score) {
            return {
                Grantham => $score
            };
        }
    }

    return {};
}

1;