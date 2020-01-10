# Contribution Guide

We welcome contributions from outside the Ensembl team and our full contribution guidelines are available [here](https://github.com/Ensembl/ensembl/blob/master/CONTRIBUTING.md)

Please submit new plugins and updates as pull requests against the MASTER branch unless otherwise discussed.

## Plugin development

Guidance on writing plugins for Ensembl VEP is given [here](https://www.ensembl.org/info/docs/tools/vep/script/vep_plugins.html)

To make development of plugins easier, we suggest you use the [Bio::EnsEMBL::Variation::Utils::BaseVepPlugin](https://github.com/Ensembl/ensembl-variation/blob/master/modules/Bio/EnsEMBL/Variation/Utils/BaseVepPlugin.pm)  module as your base class, as this provides default implementations of all the necessary methods which can be overridden as required.
The documentation in this module provides details of all required methods and a simple example of a plugin implementation. 
Also see the [Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin](https://github.com/Ensembl/ensembl-variation/blob/master/modules/Bio/EnsEMBL/Variation/Utils/BaseVepTabixPlugin.pm) for reading tabix-indexed files.


## Documentation

* Please add NAME and DESCRIPTION sections to the header - these are used to automatically create documentation pages. 
* The documentation should include enough information for others to use the plugin, including- 
  * full instructions on how to download and format any data
  * full instructions on how to download and install code
* Also acknowledge any data/tool provider with any text they specify on their website, including the project URL and latest publication where possible. 

## Licence
Please add a licence your contact details. Ensembl code is licensed under our Apache 2.0 license. Our expectation is that contributing code is made available under the same license. Any copyright assertion to other organistions should be declared in the modifying file and in the root LICENSE section. 


