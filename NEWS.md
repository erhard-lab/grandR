
# grandR 0.2.1

This is the first release of grandR.


## Major updates
* Toxicity plots
* Highlighting genes in the shiny web interface

## Smaller updates
* Compute BIC values for a list of models per gene and save it as an analysis table
* Additional parameters to PlotScatter
* Updated all calls to aes_string by tidy evaluation (for ggplot2 3.0 compatibility)
* Semantics for concentration fields
* Fixed issues when merging two data sets that were independently processed and have gene name duplettes
* FitKineticsGeneNtr now also computes residuals

# grandR 0.2.0

This is the first release of grandR.

## Available functionality

*   Functions for loading GRAND-SLAM data
*   Function for preprocessing
*   Functions for kinetic modeling of progressive labeling time courses
*   Functions to infer synthesis and half-lives from snapshot data
*   Many helper functions (plotting, read simulation, etc.)
*   Web-interface via shiny