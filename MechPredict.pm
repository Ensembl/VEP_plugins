package MechPredict;

# Enable variable checking
use strict;
use warnings;
# use Data::Dumper; # For de-bugging

# Import the Base VEP plugin class and tell Perl that the plugin inherits from the Base VEP plugin class
use Bio::EnsEMBL::Variation::Utils::BaseVepPlugin; # Gives access to core VEP plugin methods
use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin); # Allows the plugin to reuse (inherit) existing methods

# Define constructor called new
# This is called when VEP loads the plugin
# It serves to initialise core plugin config: sets up the plugin, reads user params, returns a reusable object
sub new {

    # Get class name - in this case, it's MechPredict, the name of the plugin
    my $class = shift; 

    # print "DEBUG: Raw arguments passed to new():\n";
    # print join(", ", @_), "\n";  # Print raw arguments

    # Calls the constructor of BaseVepPlugin, inheriting the structure
    my $self = $class->SUPER::new(@_); 

    # # Check what is actually stored in $self
    # print "DEBUG: Object contents after SUPER::new(): ", Dumper($self), "\n";

    # Retrieve user-specified parameters
    # In this case, $params->{file} will be the .tsv file that the use must supply
    # Manually extract parameters as params_to_hash() failed
    my %params = @_;
    # print "DEBUG: Manually extracted params = ", join(", ", map { "$_ => $params{$_}" } keys %params), "\n";

    # Exit plugin if the user doesn't supply the data file
    my $file = $params{file} || die "Error: No data file supplied for the plugin.\n";

    # Read and store the TSV data
    $self->{data} = $self->read_tsv($file);

    return $self;
}

# Define subroutine for reading in the .tsv file
sub read_tsv {
    my ($self, $file) = @_; # $self is the plugin object, $file is the file path and @_ contains all arguments passed to the function 
    my %data; # Define empty hash to store data from the .tsv
    
    # Open file for reading, assing fh as the file handle for file
    # If the file cannot be opened, then exit
    open my $fh, "<", $file or die "Could not open file '$file': $!";

    # Loop over each line of the file
    while (<$fh>) { 

        # Remove trailing \n chars
        chomp; 

        # Split line into 4 variables
        my ($gene, $uniprot_id, $mechanism, $probability) = split("\t"); 

        # Store line in the data hash, grouping the data by gene - gene will be the array reference (key, in other words)
        push @{$data{$gene}}, { uniprot_id => $uniprot_id, mechanism => $mechanism, probability => $probability };
    }
    close $fh;
    
    return \%data;
}

# Define the feature type that this data will map to - MechPredict provides gene-level annotation
# This ensures that MechPredict plugin is called once per gene
sub feature_types {
  return ['Gene'];
}

# Define the VEP output fields - these are the new fields that the plugin will add to the VEP output
sub get_header_info {
  return { 
    MechPredict_pDN => 'Probability that the gene is likely to be associated with a dominant-negative (DN) mechanism, as predicted by an SVC binary classifier model (Badonyi et al., 2024).',
    MechPredict_pGOF => 'Probability that the gene is likely to be associated with a gain-of-function (GOF) mechanism, as predicted by an SVC binary classifier model (Badonyi et al., 2024).',
    MechPredict_pLOF => 'Probability that the gene is likely to be associated with a loss-of-function (LOF) mechanism, as predicted by an SVC binary classifier model (Badonyi et al., 2024).'
  };
}

# Main logic





1;