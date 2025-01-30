package MechPredict;

# Enable variable checking
use strict;
use warnings;
use Bio::EnsEMBL::Variation::Utils::BaseVepPlugin;
use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);
use Data::Dumper; # For debugging

# Define constructor called new
# This is called when VEP loads the plugin
# It serves to initialise core plugin config: sets up the plugin, reads user params, returns a reusable object
sub new {
    my $class = shift;
    my $self = $class->SUPER::new(@_);

    # Extract plugin parameters
    my @params = @{ $self->params };  # VEP stores parameters as an array

    # Create a hash to store extracted key-value pairs
    my %params;
    foreach my $param (@params) {
        my ($key, $value) = split('=', $param, 2);  # Split "file=/path/to/file"
        $params{$key} = $value if defined $key and defined $value;
    }

    print "DEBUG: Extracted params = ", Dumper(\%params), "\n";  # Debugging

    # Ensure the file parameter was provided
    my $file = $params{file} || die "Error: No data file supplied for the plugin.\n";

    # Store file path for later use
    $self->{file} = $file;

    return $self;
}

# Define subroutine for reading in the .tsv file
sub read_tsv {
    my ( $self, $file ) = @_
      ; # $self is the plugin object, $file is the file path and @_ contains all arguments passed to the function
    my %data;    # Define empty hash to store data from the .tsv

    # Open file for reading, assing fh as the file handle for file
    # If the file cannot be opened, then exit
    open my $fh, "<", $file or die "Could not open file '$file': $!";

    # Loop over each line of the file
    while (<$fh>) {

        # Remove trailing \n chars
        chomp;

        # Split line into 4 variables
        my ( $gene, $uniprot_id, $mechanism, $probability ) = split("\t");

# Store line in the data hash, grouping the data by gene - gene will be the array reference (key, in other words)
        push @{ $data{$gene} },
          {
            uniprot_id  => $uniprot_id,
            mechanism   => $mechanism,
            probability => $probability
          };
    }
    close $fh;

    return \%data;
}

# Define the feature type that this data will map to
# This ensures that MechPredict plugin is called once per tx-gene pair
# This is because MechPredict should check if a variant is missense before applying gene-level annotation
sub feature_types {
    return ['Transcript'];
}

# Define the VEP output fields - these are the new fields that the plugin will add to the VEP output
sub get_header_info {
    return {
        MechPredict_pDN =>
'Probability that the gene is associated with a dominant-negative (DN) mechanism, as predicted by an SVC binary classifier model (Badonyi et al., 2024).',
        MechPredict_pGOF =>
'Probability that the gene is associated with a gain-of-function (GOF) mechanism, as predicted by an SVC binary classifier model (Badonyi et al., 2024).',
        MechPredict_pLOF =>
'Probability that the gene to be associated with a loss-of-function (LOF) mechanism, as predicted by an SVC binary classifier model (Badonyi et al., 2024).'
    };
}

# Run function - main logic
# This subroutine processes each variant and returns annotation results
sub run {
    my ( $self, $tva ) = @_;    # pull TranscriptVariationAllele obj from vep

  # Check if the variant is missense - scores only relevant to missense variants
  # Get all consequences associated with this tx-variant pair
    my @consequences = @{ $tva->get_all_OverlapConsequences() };

    # Loop through consequences and check whether any are missense_variant
    my $is_missense = 0;
    for my $consequence (@consequences) {

# As soon as there is a missense_variant, set is_missense to 1 and break loop using last
        if ( $consequence->SO_term eq 'missense_variant' ) {
            $is_missense = 1;
            last;
        }
    }

    # Skip annotation if no missense_variant found
    return {} unless $is_missense;

# Get the gene associated with the tx-variant pair - external_name is the HGNC symbol
    my $gene_name = $tva->transcript->gene->external_name;

    # Skip annotation if no gene_name associated
    return {} unless $gene_name;

# Check whether the gene_name from the user's vcf can be found in the MechPredict prediction data
# The first col of my data was stored as the key in the read_tsv sub: @{$data{$gene}}
# As such, can look up the gene_name in the keys of the data slot of self
# Skip annotation if the user's gene isn't in the MechPredict prediction data
    if ( !exists $self->{data}{$gene_name} ) {
        return {};
    }

    # Pull out MechPredict prediction data for gene_name
    my $data = $self->{data}{$gene_name};

    # Pull out probability values for each mechanism
    my $pdn  = $data->{pDN};
    my $pgof = $data->{pGOF};
    my $plof = $data->{pLOF};

    # Define hash containing probability thresholds for each mechanism
    my %thresholds = (
        pdn  => 0.61,    # Probability of dominant-negative mechanism
        pgof => 0.63,    # Probability of gain-of-function mechanism
        plof => 0.64     # Probability of loss-of-function mechanism
    );

    # Create interpretation field
    my $interpretation = "";

    # Compare values to thresholds and populate interpretation
    if (   $pdn >= $thresholds{pdn}
        && $pgof < $thresholds{pgof}
        && $plof < $thresholds{plof} )
    {
        $interpretation =
          "Gene likely associated with a dominant-negative mechanism";
    }
    elsif ($pgof >= $thresholds{pgof}
        && $pdn < $thresholds{pdn}
        && $plof < $thresholds{plof} )
    {
        $interpretation =
          "Gene likely associated with a gain-of-function mechanism";
    }
    elsif ($plof >= $thresholds{plof}
        && $pgof < $thresholds{pgof}
        && $pdn < $thresholds{pdn} )
    {
        $interpretation =
          "Gene likely associated with a loss-of-function mechanism";
    }
    else {
        $interpretation = "No conclusive dominant mechanism detected";
    }

    # Add 4 fields to the VEP output
    return {
        MechPredict_pDN => $data->{pDN}
        ,    # Probability of dominant-negative mechanism
        MechPredict_pGOF => $data->{pGOF}
        ,    # Probability of gain-of-function mechanism
        MechPredict_pLOF => $data->{pLOF}
        ,    # Probability of loss-of-function mechanism
        MechPredict_interpretation => $interpretation
    };

}

1;
