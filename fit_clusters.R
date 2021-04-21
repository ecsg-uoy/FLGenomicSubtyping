#######################################################
#
# This script generates the Bernoulli
# mixture model clusters (AIC and ICL)
# that were used in the paper.
#
# Dependencies:
#    - flexmix
#    - tidyverse
#    - gridExtra
#
#######################################################

library(flexmix)
library(tidyverse)
library(gridExtra)

######################################
#   Fit original clusters
######################################
# This first step fits the published clusters using the supplied processed genomic data

# Load library functions used for identifying enriched mutations and plotting heatmaps
source("plotting_functions.R")

# Load data (NB: this is also available in the supplementary material provided with this publication)
muts_all <- read.csv("data/genomic_data.csv")

# There are 548 patients with 32 columns. 
# First column is identifier, remaining 31 are genetic features
dim(muts_all)

# Remove patient id, restricting columns to 117 mutations
muts_df <- muts_all[, -1]
dim(muts_df)

# For reproducibility with the published clusters keep the following 2 lines the same.
# They set the seed to the same used in the paper
RNGversion("4.0.4")
set.seed(17)

# Run the mixture model fitting algorithm through all the possible cluster numbers from 1 to 10,
# running a thousand repetitions of each and selecting the best

ex <- initFlexmix(as.matrix(muts_df) ~ 1,
                  k = 1:10, model = FLXMCmvbinary(),
                  control = list(minprior = 0.10), nrep = 1000, unique=T)

# Visualise number of clusters suggested by AIC and BIC/ICL (3 and 2 respectively)
plot(ex)

# Fit a cluster model with three clusters (according to the AIC).
set.seed(17)
ex_3 <- initFlexmix(as.matrix(muts_df) ~ 1,
                    k = 3, model = FLXMCmvbinary(),
                    control = list(minprior = 0.10), nrep = 1000, unique=T)

# Fit a cluster model with two clusters (according to the IC).
set.seed(17)
ex_2 <- initFlexmix(as.matrix(muts_df) ~ 1,
                    k = 2, model = FLXMCmvbinary(),
                    control = list(minprior = 0.10), nrep = 1000, unique=T)


# Save cluster assignments in main data frame
muts_all$ClusterAIC <- as.factor(paste0('C', ex_3@cluster))
muts_all$ClusterBIC <- as.factor(paste0('C', ex_2@cluster))

summary(muts_all$ClusterAIC)
summary(muts_all$ClusterBIC)

# Plot heatmaps of enriched mutations
genes <- colnames(muts_df)
plt_aic <- heatmap_mutation_extended(muts_all, genes, 'ClusterAIC', y_order = 'group', idcol = 'PID')
grid.arrange(plt_aic)

plt_bic <- heatmap_mutation_extended(muts_all, genes, 'ClusterBIC', y_order = 'group', idcol = 'PID')
grid.arrange(plt_bic)

######################################
#   Assign new samples to clusters
######################################
# This section demonstrates how to input new data into the mixture model to view
# where they would be assigned.

# Function to predict membership for new patients
# Args:
#    - model: A flexmix object.
#    - newdata: Data frame with the same number of columns as muts_df and in the same order.
#
# Returns:
#   - A matrix where each row corresponds to a row of newdata and each column provides the probability
#     of belonging to that cluster. Dimensions NxC where N = nrow(newdata) and C = number of clusters
#     that the model identified


predict_clusters <- function(model, newdata) {
    # Flexmix looks for a data frame with the same name as that used to
    # fit the model in the first place, so we need to overwrite this
    # variable name at the higher scope
    muts_df <<- newdata
    probs <- flexmix::posterior(model, newdata=newdata)
    probs
}

# Predicting cluster for just first 10 samples in dataset
# Having to copy original data frame into new variable as the
# prediction process will overwrite this variable
orig_df <- muts_df
test_data <- orig_df[1:10, ]

# AIC model returns 3 columns
predict_clusters(ex_3, test_data)

# BIC model returns 2 columns
predict_clusters(ex_2, test_data)

# View predicted cluster membership for an imaginary sample.

# This line creates a 1 row data frame with the same columns as muts_df, all set to 0 (no mutation)
new_df <- as.data.frame(setNames(lapply(genes, function(x) 0), genes))
# Individual mutations can be set
new_df$ARID1A <- 1
new_df$IRF8 <- 1
new_df$MYC_S <- 1

# From the heatmap it can be observed that MYD88 is the 5th cluster, which is consistent
# with the predicted 93% assignment for this imaginary sample.
predict_clusters(ex_3, new_df)
