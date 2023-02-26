
# grandR 0.2.1

This is the release of grandR including improvements for the revised version of the manuscript.

## Major updates
* Enabled interface with Seurat
* Implemented estimation for pulse-chase experiments
* New vignettes for single cell data and pulse-chase data
* 4sU dropout plots
* Highlighting genes in the shiny web interface
* Cheatsheet

## Smaller updates
* Compute BIC values for a list of models per gene and save it as an analysis table
* Additional parameters to PlotScatter
* Updated all calls to aes_string by tidy evaluation (for ggplot2 3.0 compatibility)
* Semantics for concentration fields
* Fixed issues when merging two data sets that were independently processed and have gene name dublettes
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
