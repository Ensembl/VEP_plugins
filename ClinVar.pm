=head
################################################################################
# ClinVar Plugin                                                              #
#                                                                             #
# This plugin is designed for the Ensembl Variant Effect Predictor (VEP) to   #
# annotate variants using data from the ClinVar database.                     #
#                                                                             #
# LICENSE:                                                                    #
# This plugin is open-source and distributed under the MIT license. Use and   #
# redistribution are permitted provided that proper credit is given to the    #
# authors.                                                                    #
#                                                                             #
# CITATION:                                                                   #
# If you use this plugin, please cite the ClinVar database:                   #
# Landrum MJ, Lee JM, Benson M, et al. ClinVar: public archive of             #
# interpretations of clinically relevant variants. Nucleic Acids Research.    #
# 2016 Jan 4;44(D1):D862-8. doi:10.1093/nar/gkv1222.                          #
#                                                                             #
# USAGE:                                                                      #
# To use this plugin, provide the path to a tabix-indexed ClinVar VCF file.   #
# Example command:                                                            #
#                                                                             #
# perl vep --input_file input.vcf --output_file output.txt \                  #
#     --plugin ClinVar,ClinVar.vcf.gz --cache --dir_cache ~/.vep              #
#                                                                             #
# REQUIREMENTS:                                                               #
# - Ensembl VEP                                                               #
# - Tabix for indexing the ClinVar VCF file                                   #
#                                                                             #
# Developed to facilitate clinical variant interpretation by leveraging the   #
# ClinVar database within the VEP framework.                                  #
################################################################################

# ClinVar Plugin

This is a plugin for the Ensembl Variant Effect Predictor (VEP) to annotate variants with data from a ClinVar VCF file.

## Requirements

- **VEP**: Ensure you have VEP installed. Refer to the [VEP Installation Guide](https://www.ensembl.org/info/docs/tools/vep/index.html) if necessary.
- **ClinVar VCF File**: Obtain the latest ClinVar VCF file from [ClinVar](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/).
- **tabix**: Used to index the ClinVar VCF file for faster lookups.

## Usage

**Prepare the ClinVar VCF File**:
   Compress and index the ClinVar VCF file using `bgzip` and `tabix`:

   ```bash
   bgzip ClinVar.vcf
   tabix -p vcf ClinVar.vcf.gz
   ```

Run VEP with the plugin using the following command:

```bash
perl vep \
    --input_file input.vcf \
    --output_file output.txt \
    --plugin ClinVar,ClinVar.vcf.gz \
    --cache \
    --dir_cache ~/.vep
```
=cut

package ClinVar;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use Bio::EnsEMBL::Variation::Utils::Sequence qw(get_matched_variant_alleles);
use Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin;
use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin);

sub new {
    my $class = shift;
    my $self = $class->SUPER::new(@_);
    return $self;
}

sub feature_types {
    return ['Variation'];
}

sub get_header_info {
    return {
        ClinVar => 'ClinVar annotations from a VCF file',
    };
}

sub run {
    my ($self, $tva) = @_;
    # Plugin logic to annotate the variant
    return {};
}

1;
