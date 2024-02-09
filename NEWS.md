# grandR 0.2.4

## Shiny improvement
* If there are global plots, there is now an additional page showing the highlighted genes
* Defer now allows to specify desired figure width and height
* Fixed bug that caused all global plots to be executed before rendering the main table
* Html files that are in the same ordner are linked (under menu Reports)
* Now multiple tables can be provided (as named list)
* Additiona, static plots can be shown in floating windows

## Smaller updates
* fixed subsetting of grandR to only a single column
* improved handling of factors in coldata during merging
* Updated ReadGRAND3 to most recent Grand3 version (mode output)
* Improved ComputeColumnStatistics (percentage calculations)
* Differential gene expression analysis for a subset of the genes (including normalization to only those genes)
* FormatCorrelation can add RMSDs
* Statistics for Plot4sUDropout()

# grandR 0.2.3

## Smaller updates
* compute.M for LFC
* size factors to normalize Pairwise, LFC and PairwiseDESeq2


# grandR 0.2.2

This is the release of grandR including improvements for the submission of the grand-Correct preprint

## Major updates
* Enable simulation with non-constant rates
* 4sU dropout plots and grand-Correct functionality

## Smaller updates
* Fixed a bug that caused warning messages after loading files.


# grandR 0.2.2

This is the release of grandR including improvements for the submission of the grand-Correct preprint

## Major updates
* Enable simulation with non-constant rates
* 4sU dropout plots and grand-Correct functionality

## Smaller updates
* Fixed a bug that caused warning messages after loading files.


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
