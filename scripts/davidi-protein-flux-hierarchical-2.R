
### Read in relevant libraries.
library(readxl)
library(magrittr)
library(dplyr)
library(reshape2)
library(ggplot2)
library(forcats)
library(stringr)
library(rstan)
library(tidybayes)
library(infotheo)

'%!in%' <- function(x,y)!('%in%'(x,y))

two_level_hier_model_e_coli <- stan_model("scripts/fitness_flux_phi_hierarchical.stan")

## fitting the model in stan
fit_fitness_path_two_hierarch <- sampling(two_level_hier_model_e_coli,
                                          data = list(N = nrow(genes_fluxes_long_df_z_score_cc_2_level),
                                                      K = length(unique(genes_fluxes_long_df_z_score_cc_2_level$gene_name)),
                                                      id = genes_fluxes_long_df_z_score_cc_2_level$gene_as_number,
                                                      flux = genes_fluxes_long_df_z_score_cc_2_level$z_score_flux,
                                                      protein_amounts = genes_fluxes_long_df_z_score_cc_2_level$z_score_prot,
                                                      mean_protein_amount = genes_fluxes_summary_states_merged_cc$log_prot_amount_mean_weighted,
                                                      substrate_connectivity = genes_fluxes_summary_states_merged_cc$geometric_mean,
                                                      delta_g = genes_fluxes_summary_states_merged_cc$delta_g),
                                          iter = 4000,
                                          chains = 4)

gamma0_e_coli <- fit_fitness_path_two_hierarch %>%
  spread_draws(gamma_0) %>%
  ggplot(aes(x = gamma_0)) +
  geom_histogram() +
  geom_vline(xintercept = 0) +
  theme_bw() +
  ylab('Sample Count') +
  xlab(expression(gamma[1])) +
  ggtitle('A. Intercept')

gamma1_e_coli <- fit_fitness_path_two_hierarch %>%
  spread_draws(gamma_1) %>%
  ggplot(aes(x = gamma_1)) +
  geom_histogram() +
  geom_vline(xintercept = 0) +
  theme_bw() +
  ylab('Sample Count') +
  xlab(expression(gamma[1]~(log[e](mmol~Amino~Acid~cGDW^-1)^-1))) +
  theme(text = element_text(size = 16)) +
  ggtitle('D. Protein cost hypothesis');gamma1_e_coli

gamma2_e_coli <- fit_fitness_path_two_hierarch %>%
  spread_draws(gamma_2) %>%
  ggplot(aes(x = gamma_2)) +
  geom_histogram() +
  geom_vline(xintercept = 0) +
  theme_bw() +
  ylab('Sample Count') +
  xlab(expression(gamma[2]~(Substrate~Connectivity~Metric^-1))) +
  theme(text = element_text(size = 16)) +
  ggtitle('E. Substrate connectivity hypothesis')

gamma3_e_coli <- fit_fitness_path_two_hierarch %>%
  spread_draws(gamma_3) %>%
  ggplot(aes(x = gamma_3))  +
  geom_histogram() +
  geom_vline(xintercept = 0) +
  theme_bw() +
  ylab('Sample Count')  +
  xlab(expression(gamma[3]~((kJ~Mol^-1)^-1))) +
  theme(text = element_text(size = 16)) +
  ggtitle(expression(paste('F. ', Delta, 'G hypothesis')))

e_coli_2_level_out <- ggarrange(gamma1_e_coli, gamma2_e_coli, gamma3_e_coli, nrow = 1, align = 'hv')

ggsave(e_coli_2_level_out, 
       filename = 'figures/e_coli_2_level_out.png',
       height = 5, width = 16)

combined_e_coli_yeast_pars <- ggarrange(yeast_2_level_out, 
                                        e_coli_2_level_out, 
                                        nrow = 2)
ggsave(combined_e_coli_yeast_pars, 
       filename = paste0("figures/figsX_", current_date, "_combined_e_coli_yeast_pars.pdf"), 
       height = 8, width = 20)


