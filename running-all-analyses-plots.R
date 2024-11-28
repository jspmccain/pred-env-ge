#### This script runs all analyses and makes all figures.

## This is a flag to re-run predictive models, which take a lot of time. It's only relevant for 'hackett-full-proteome-prediction.R'.
rerun_all_pred_models <- FALSE
current_date <- format(Sys.Date(), "%Y%m%d")

# Data processing for main analyses
source('scripts/davidi-data-processing.R')#*%(Note this also runs ecoli_fasta_formatting.R)
source('scripts/hackett-data-processing.R')#*%

# Note that the Chubukov data processing is done in the same script that the hierarchical model is run.
source('scripts/chubukov-transcript-flux-hierarchical.R')#*%

# Run hierarchical model analyses
source('scripts/davidi-protein-flux-hierarchical.R')#*%
source('scripts/hackett-protein-flux-hierarchical.R')#*%

# Run single protein to rate analysis (dif functional forms, i.e. pearson/spearman/mutual info)
source('scripts/hackett-davidi-chubukov-single-proteins-additional-metrics.R')#*%

# Plot results from hierarchical model analysis
## Plot of individual coefficients across rate-protein pairs.
## Make dataframe and plots for posterior probability distributions of group level distributions across taxa.
source('scripts/plotting_hierarchical_model_partial_pooling.R')#*%

# Plot two example single protein-rate relationships, and plot the 
# single protein-rate metrics.
source('scripts/plotting-single-protein-rate-relationships.R')#*%

# Run hierarchical model analyses -- 2 level hierarchical model
source('scripts/hackett-protein-flux-hierarchical-2.R') #*%
source('scripts/davidi-protein-flux-hierarchical-2.R') #*%

# Run analysis (LASSO CV) for Hackett full proteome/subsystem rate prediction
source('scripts/hackett-full-proteome-prediction.R') #*%

# Run analysis for Hackett full proteome rate prediction used protein groups as predictors.
source('scripts/hackett-full-proteome-prediction-protein-groups.R')

# Run single protein prediction using a simple linear model 
source('scripts/single-protein-prediction.R') #%

# Plot R^2 comparison from single proteins, pathway-level LASSO, and whole-proteome predictions
source('scripts/plotting-r-squared-distributions.R') #%

# Plot relationship between conserved proteins and slopes across E. coli and yeast.
source('scripts/plotting-yeast-e-coli-homologues.R') ##*%

# Plot a single example (for initial figure)
source('scripts/plotting-single-protein-rate-relationships.R') #%

# Plot the relationship between ridge regression models using whole proteomes and the rate correlation with growth rate
source('scripts/plotting-growth-rate-model-comparison.R') #%

# Plot the aspartate kinase case study
source('scripts/plotting-single-rates-from-whole-proteome-prediction.R')#%

# Plot the glutathione subsystem prediction case study
source('scripts/plotting-glutathione-subsystem.R')#%

# Do the GO term analysis and plot it
source('scripts/go-term-analysis.R')#%

# Make the correlation plot for examining measurement error.
source('scripts/hackett-growth-proteomes-profiles.R')#%

# Plotting the LASSO results
source('scripts/plotting-hackett-non-zero-coef-lasso-whole-proteome.R')
