This folder contains code and data to reproduce analyses from the paper *Molecular Subclusters of Follicular Lymphoma: a report from the UK’s Haematological Malignancy Research Network.

The contents are as follows:

  - `fit_clusters.R`: Script to run the mixture modelling algorithm over the processed genomic data, identifying the AIC and BIC clusters that were reported in the paper. Also plots the heatmaps of enriched mutations and provides the functionality to predict cluster membership for new samples.
  - `plotting_functions.R`: File containing functions used to identify enriched mutations and plot the heatmaps.
  - `data`: Folder containing a csv file
    - `genomic_data.csv`: The processed binary mutation data for the 548 patients of the 31 genetic features used in the clustering.
