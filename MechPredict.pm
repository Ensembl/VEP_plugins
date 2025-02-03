=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute  
Copyright [2016-2024] EMBL-European Bioinformatics Institute  

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

MechPredict  

=head1 SYNOPSIS  

 mv MechPredict.pm ~/.vep/Plugins  
 ./vep -i input.vcf --plugin MechPredict,file=mechpredict_data.tsv  

=head1 DESCRIPTION  

This is a plugin for the Ensembl Variant Effect Predictor (VEP) that annotates missense variants with predicted 
dominant-negative (DN), gain-of-function (GOF), or loss-of-function (LOF) mechanisms derived from a
Support Vector Classification (SVC) model (Badonyi et al., 2024).

Note:
- The plugin requires MechPredict_input.tsv, a pre-processed prediction dataset in TSV format. 
- The wrangled file should contain gene-level probabilities for the three mechanism categories.  
- The plugin adds the following fields to the VEP output:  
  - `MechPredict_pDN`: Probability of a **dominant-negative (DN) mechanism**  
  - `MechPredict_pGOF`: Probability of a **gain-of-function (GOF) mechanism**  
  - `MechPredict_pLOF`: Probability of a **loss-of-function (LOF) mechanism**  
  - `MechPredict_interpretation`: Statement of the most likely mechanism based on empirically-derived cutoffs from Badonyi et al., 2024. 

Usage:
1. The raw data from the Badonyi et al., 2024 manuscript can be pre-processed using the folllwoing command: 
  - GOF: https://osf.io/h45ns
  - DN: https://osf.io/xfy38
  - LOF https://osf.io/dj4qg

2. The plugin input data can then be prepared from the raw data using:  
   ```bash
    cut --complement -f4 pdn_svm_poly_2023-07-25.tsv | awk '{print $1 " " $2 "\t" $0}' | sort > pdn_mod.tsv && \
    cut --complement -f4 pgof_svm_poly_2023-07-25.tsv | awk '{print $1 " " $2 "\t" $0}' | sort > pgof_mod.tsv && \
    cut --complement -f4 plof_svm_poly_2023-07-28.tsv | awk '{print $1 " " $2 "\t" $0}' | sort > plof_mod.tsv && \
    join -t $'\t' -1 1 -2 1 pdn_mod.tsv pgof_mod.tsv | join -t $'\t' -1 1 -2 1 - plof_mod.tsv | cut --complement -f1,5,6,8,9 | sed '1i gene uniprot_id pDN pGOF pLOF' > MechPredict_input.tsv && \
    rm pdn_mod.tsv pgof_mod.tsv plof_mod.tsv && \
    awk 'BEGIN {print "gene\tuniprot_id\tmechanism\tprobability"} NR>1 {print $1, $2, "DN", $3; print $1, $2, "GOF", $4; print $1, $2, "LOF", $5;}' OFS='\t' MechPredict_input.tsv > MechPredict_input_pivot.tsv && \
    mv MechPredict_input_pivot.tsv MechPredict_input.tsv
   ```

3. VEP can be run with the MechPredict plugin as follows:  
   ```bash
   ./vep -i variations.vcf --plugin MechPredict,file=/path/to/mechpredict_data.tsv
    ```

Citation:
Badonyi M, Marsh JA (2024) Proteome-scale prediction of molecular mechanisms underlying dominant genetic diseases. 
PLoS ONE 19(8): e0307312. https://doi.org/10.1371/journal.pone.0307312

=cut

package MechPredict; 

# -- Setup ---------------------------------------------------------------------

# Enable variable checking
use strict;
use warnings;

# Gives access to core VEP plugin methods - provides a set of methods that can be reused by the plugin
use Bio::EnsEMBL::Variation::Utils::BaseVepPlugin; 

# Allows the plugin to reuse (inherit) existing methods - sets up object-oriented inheritance
use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin); 

# Permits debugging - built-in module for printing complex data structures
use Data::Dumper; 

# -- Initialise plugin ---------------------------------------------------------

# Define constructor called new
# This is called when VEP loads the plugin
# It serves to initialise core plugin config: sets up the plugin, reads user params, returns a reusable object
sub new { 

    # Retrieves the class name from the first argument passed to the subroutine - VEP passes the module name (MechPredict)
    my $class = shift; 

    # Calls the parent class (BaseVepPlugin), assigning its properties to new
    my $self = $class->SUPER::new(@_); 

    # Parse the parameters passed to the plugin
    # params is a method from BaseVepPlugin that returns the parameters passed to the plugin
    # Extract plugin parameters and store the values in an array
    # @(...) is a way to dereference an array reference like unlisting in R
    my @params = @{ $self->params };

    # Declare a hash to store extracted key-value pairs from the args
    my %params;
    foreach my $param (@params) {
        my ($key, $value) = split('=', $param, 2);  # Split "file=/path/to/file"
        $params{$key} = $value if defined $key and defined $value;
    }

    # Debugging
    # print "DEBUG: Extracted params = ", Dumper(\%params), "\n";

    # Ensure the file parameter was provided
    my $file = $params{file} || die "Error: No data file supplied for the plugin.\n";

    # Store file path in the object for later use
    $self->{file} = $file;

    # Supply path to read_tsv sub and store the result in the data slot of self
    $self->{data} = $self->read_tsv($file);

    # Returns the initialised plugin object - new is now populated and can be utilised by the plugin
    return $self;
}

# -- Define key subroutines ----------------------------------------------------

# Define subroutine for reading in the .tsv file
sub read_tsv {

    # Retreive the plugin object 
    # $self is the plugin object, $file is the file path
    # @_ is a special array in perl - it contains the arguments passed to the subroutine
    # In this case, two are passed - the plugin object and the file path
    # So, these need to be retrieved from the array and assigned to $self and $file
    my ( $self, $file ) = @_; 

    # Declare empty hash to store data from the .tsv
    my %data;    

    # Open file for reading, assing fh as the file handle for file
    # If the file cannot be opened, then exit
    open my $fh, "<", $file or die "Could not open file '$file': $!";

    # Loop over each line of the file - more mem efficient than reading the entire file into memory with a foreach
    while (<$fh>) {

        # Remove trailing \n chars
        chomp;

        # There are 4 cols in the .tsv file, so assign result of split to 4 usefully named variables
        # $_ means the current line of the file
        my ($gene, $uniprot_id, $mechanism, $probability) = split("\t", $_);

        # Store line in the data hash, grouping the data by gene - gene will be the array reference (key, in other words)
        push @{ $data{$gene} },
          {
            uniprot_id  => $uniprot_id,
            mechanism   => $mechanism,
            probability => $probability
          };
    }
    close $fh;

    # Debugging
    # print "DEBUG: Read ", scalar(keys %data), " genes from the file.\n";

    # Return the data hash
    return \%data;
}

# Defines the feature types that the plugin will run on
# In this case, Transcript, because the plugin will need to check if the variant is missense
# Transcript: MechPredict plugin is called once per transcript-variant pair
sub feature_types {
    return ['Transcript'];
}

# Define the VEP header annotation fields
sub get_header_info {
    return {
        MechPredict_pDN =>
'Probability that the gene is associated with a dominant-negative (DN) mechanism, as predicted by an SVC binary classifier model (Badonyi et al., 2024).',
        MechPredict_pGOF =>
'Probability that the gene is associated with a gain-of-function (GOF) mechanism, as predicted by an SVC binary classifier model (Badonyi et al., 2024).',
        MechPredict_pLOF =>
'Probability that the gene is associated with a loss-of-function (LOF) mechanism, as predicted by an SVC binary classifier model (Badonyi et al., 2024).', 
        MechPredict_interpretation => 
'Interpretation of the probabilities based on thresholds reccomended by Badonyi et al., 2024.'
    };
}

# -- Main logic ----------------------------------------------------------------

# This subroutine will be executed once per transcript that a variant overlaps
sub run {

    # VEP creates a TranscriptVariationAllele (tva) object for each transcript-variant pair
    # This object contains all the information needed to annotate the variant
    # Pull the tva object from the plugin object
    my ( $self, $tva ) = @_;

    # Debugging
    # print "DEBUG: Processing variant...\n";
    # print "DEBUG: TranscriptVariationAllele = ", Dumper($tva), "\n";

    # Debugging
    # print("--------------------\n");

    # Get transcript ID
    my $transcript = $tva->transcript;
    return {} unless $transcript;

    # Get gene name
    # Usually, ->get_Gene->external_name would be used, but this doesn't work in offline mode, so pull from cached value
    # eval is used to catch any errors that may occur if the gene object doesn't exist
    # Maybe add an if statement here based on whether this is offline/cached or not? Then would work both ways
    # _gene_symbol is a cached value in offline mode
    my $gene_name = eval { $transcript->{_gene_symbol} };
        return {} unless $gene_name;

    # Debugging
    # print "DEBUG: Processing transcript ", $transcript->stable_id, "\n";
    # print "DEBUG: Processing gene ", $gene_name, "\n";

    # Check if the variant has a missense consequence
    # Get all consequences
    my @consequences = @{ $tva->get_all_OverlapConsequences() };

    # Check if any consequence is missense_variant
    my $is_missense = grep { $_->SO_term eq 'missense_variant' } @consequences;

    if (!$is_missense) {
        # Debugging
        # print "DEBUG: No missense variant for transcript ", $transcript->stable_id, "\n";
        # print "DEBUG: Consequences for ", $transcript->stable_id, " = ", join(", ", map { $_->SO_term } @consequences), "\n";
        return {};
    }

    # Debugging
    # print "DEBUG: Missense variant found for transcript ", $transcript->stable_id, "\n";

    # Check whether the gene_name from the user's vcf can be found in the MechPredict prediction data
    # The first col of the input data was stored as the key in the read_tsv sub: @{$data{$gene}}
    # As such, can look up the gene_name in the keys of the data slot of self
    # Skip annotation if the user's gene isn't in the MechPredict prediction data

    # Debugging
    print "DEBUG: Searching for gene '$gene_name' in dataset...\n";

    if ( !exists $self->{data}{$gene_name} ) {
        # Debugging
        # print "DEBUG: Gene $gene_name not found in MechPredict dataset.\n";
        return {};
    }

    # Pull out MechPredict prediction data for gene_name
    my $data = $self->{data}{$gene_name};

    # Debugging
    # print "DEBUG: Prediction data = ", Dumper($data), "\n";

    # Initialise empty variables to hold the values
    my ($pdn, $pgof, $plof) = (undef, undef, undef);

    # Iterate through the array of hashes
    foreach my $entry (@$data) {
        if ($entry->{mechanism} eq 'DN') {
            $pdn = $entry->{probability};
        }
        elsif ($entry->{mechanism} eq 'GOF') {
            $pgof = $entry->{probability};
        }
        elsif ($entry->{mechanism} eq 'LOF') {
            $plof = $entry->{probability};
        }
    }

    # Debugging
    # print "DEBUG: Extracted probabilities - pDN: ", ($pdn // 'NA'), 
    #     ", pGOF: ", ($pgof // 'NA'), 
    #     ", pLOF: ", ($plof // 'NA'), "\n";

    # Define hash containing probability thresholds for each mechanism
    my %thresholds = (
        pdn  => 0.61,    # Probability of dominant-negative mechanism
        pgof => 0.63,    # Probability of gain-of-function mechanism
        plof => 0.64     # Probability of loss-of-function mechanism
    );

    # Create interpretation field
    my $interpretation = "";

    # Compare values to thresholds and populate interpretation
    # I do not capture the cases where: 
    #   Two probabilities are high at the same time
    #   All probabilities are below their respective thresholds
    #   All probabilities are above their respective thresholds
    # Instead, these end up categoried as "No conclusive dominant mechanism detected"
    if (   $pdn >= $thresholds{pdn}
        && $pgof < $thresholds{pgof}
        && $plof < $thresholds{plof} )
    {
        $interpretation = "gene_predicted_as_associated_with_dominant_negative_mechanism";
    }
    elsif ($pgof >= $thresholds{pgof}
        && $pdn < $thresholds{pdn}
        && $plof < $thresholds{plof} )
    {
        $interpretation = "gene_predicted_as_associated_with_gain_of_function_mechanism";
    }
    elsif ($plof >= $thresholds{plof}
        && $pgof < $thresholds{pgof}
        && $pdn < $thresholds{pdn} )
    {
        $interpretation = "gene_predicted_as_associated_with_loss_of_function mechanism";
    }
    else {
        $interpretation = "no_conclusive_mechanism_detected";
    }

    # Debugging
    # print "DEBUG: Interpretation = $interpretation\n";

    # Add 4 fields to the VEP output
    return {
        MechPredict_pDN => $pdn, # Probability of dominant-negative mechanism
        MechPredict_pGOF => $pgof, # Probability of gain-of-function mechanism
        MechPredict_pLOF => $plof, # Probability of loss-of-function mechanism
        MechPredict_interpretation => $interpretation
    };
}

1;