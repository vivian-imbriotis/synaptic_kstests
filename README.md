# KS Tests are innappropriate for miniature synaptic events



## Data Generation (data_generation.R)

This file contains procedures for generating paired and unpaired datasets.  
It also contains some simple operations performed on datasets - creating a basic plot of a dataset, and checking whether the dataset is ks-positive (i.e. the KS test will return p<alpha) or lmer-positive (the mixed model comparison will return p<alpha).  
Lastly, it includes functions to generate many datasets and compute the false positive rate (FPR) and power of the KS and mixed model approachs for a given set of data-generating parameters.

## Data Gen Viz (data_gen_vizualization.R)

Production of plots that demonstrate the ways in which data generated as above can vary.

## Heatmaps (figures_heatmaps.R)

Heatmap generation, where FPR and power are encoded as a color while 2 data generating parameters are enconded on the x and y axes

## Line Charts (figures_linecharts.R)

Line chart generation, where FPR and power are on the y-axis

## Mixed model convergence checks

These files are unit tests to confirm that the datasets have the properties that they should and that the linear mixed models can asymptotically extract the appropriate variances and fixed effects.

## Real Data (real_data_vizualization)

This file simply reads in some real-world data, plots it, and uses a mixed model to determine the ratio of the within- and between-neuron variances, so they can be used as defaults for the data generation process.