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

# -- Define subroutines --------------------------------------------------------

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
'Probability that the gene to be associated with a loss-of-function (LOF) mechanism, as predicted by an SVC binary classifier model (Badonyi et al., 2024).'
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
    # print "\n🔍 DEBUG: Processing variant...\n";
    # print "DEBUG: TranscriptVariationAllele = ", Dumper($tva), "\n";

    # Check if the variant is missense - scores only relevant to missense variants
    # Get all consequences associated with this tx-variant pair
    my @consequences = @{ $tva->get_all_OverlapConsequences() };

    # Debugging
    # print "DEBUG: Consequences = ", join(", ", map { $_->SO_term } @consequences), "\n";

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

    # Get transcript
    my $transcript = $tva->transcript;
    return {} unless $transcript;

    # Debugging
    # print "DEBUG: Transcript ID = ", $transcript->stable_id, "\n";

    # Extract gene object from transcript
    # Usually, ->get_Gene->external_name would be used, but this doesn't work in offline mode, so pull from cached value
    # eval is used to catch any errors that may occur if the gene object doesn't exist
    # Maybe add an if statement here based on whether this is offline/cached or not? Then would work both ways
    # _gene_symbol is a cached value in offline mode
    my $gene_name = eval { $transcript->{_gene_symbol} }; 

    # Debugging
    # if (!$gene_name) {
    #     print "DEBUG: No gene name found for transcript ", $transcript->stable_id, "\n";
    #     return {};
    # }

    # Check whether the gene_name from the user's vcf can be found in the MechPredict prediction data
    # The first col of the input data was stored as the key in the read_tsv sub: @{$data{$gene}}
    # As such, can look up the gene_name in the keys of the data slot of self
    # Skip annotation if the user's gene isn't in the MechPredict prediction data

    # Debugging
    # print "DEBUG: Searching for gene '$gene_name' in dataset...\n";

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
        $interpretation = "Gene likely associated with a dominant-negative mechanism";
    }
    elsif ($pgof >= $thresholds{pgof}
        && $pdn < $thresholds{pdn}
        && $plof < $thresholds{plof} )
    {
        $interpretation = "Gene likely associated with a gain-of-function mechanism";
    }
    elsif ($plof >= $thresholds{plof}
        && $pgof < $thresholds{pgof}
        && $pdn < $thresholds{pdn} )
    {
        $interpretation = "Gene likely associated with a loss-of-function mechanism";
    }
    else {
        $interpretation = "No conclusive dominant mechanism detected";
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