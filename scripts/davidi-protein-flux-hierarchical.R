######
# This script analyzes data from E. coli protein abundances and corresponding fluxes, specifically fitting a hierarchical Bayesian model to estimate correlation coefficients.
######

### Read in relevant libraries.
library(readxl)
library(magrittr)
library(dplyr)
library(reshape2)
library(randomForestSRC)
library(ggplot2)
library(forcats)
library(stringr)
library(rstan)
library(tidybayes)
library(infotheo)

'%!in%' <- function(x,y)!('%in%'(x,y))

# Estimating the pearson correlation coefficient using a Bayesian hierarchical model -----------------------------------------------

fitness_model_partial_davidi <- stan_model("scripts/flux_protein_hierarchical_partial_pooling.stan")

number_of_proteins_dav <- genes_fluxes_long_df_z_score_cc$gene_as_number %>% unique() %>% length()
protein_amounts_spaced_sampled_dav <- rep(seq(from = -2, to = 2, by = 0.05), number_of_proteins_dav)
protein_ids_spaced_sampled_dav <- rep(c(1:number_of_proteins_dav), each = 81)

davidi_fitness_partial_pool <- sampling(fitness_model_partial_davidi, 
                                     data = list(N = nrow(genes_fluxes_long_df_z_score_cc),
                                                 K = length(unique(genes_fluxes_long_df_z_score_cc$gene_name)),
                                                 id = genes_fluxes_long_df_z_score_cc$gene_as_number,
                                                 flux = genes_fluxes_long_df_z_score_cc$z_score_flux,
                                                 protein_amounts = genes_fluxes_long_df_z_score_cc$z_score_prot,
                                                 protein_amounts_spaced = protein_amounts_spaced_sampled_dav,
                                                 N_spaced = length(protein_amounts_spaced_sampled_dav),
                                                 id_spaced = protein_ids_spaced_sampled_dav), 
                                     iter = 5000, 
                                     chains = 4)

# Examining the model parameters ------------------------------------------

beta_1_estimates_davidi_df <- davidi_fitness_partial_pool %>% 
  spread_draws(beta_1[n]) %>% 
  median_qi() %>% 
  dplyr::rename(gene_as_number = n) %>% 
  inner_join(genes_fluxes_long_df_z_score_cc_gene_as_number_summary,
             by = 'gene_as_number')

beta_1_estimates_davidi_df_summary_states <- davidi_fitness_partial_pool %>% 
  spread_draws(beta_1[n]) %>% 
  median_qi() %>% 
  dplyr::rename(gene_as_number = n) %>% 
  inner_join(genes_fluxes_summary_states_merged,
             by = 'gene_as_number')
