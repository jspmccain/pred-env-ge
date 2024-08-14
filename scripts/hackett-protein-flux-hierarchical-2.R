# fitting the two level hierarchical model to examine variation in slope estimates

library(bayesplot)

# Fitting model with two hierarchies and looking at par estimates ------------------------------------------------

# Loading in Stan model
fitness_model_path_yeast <- stan_model("scripts/fitness_flux_phi_hierarchical.stan")

## Fitting the model in stan
fit_fitness_path_two_hierarch_yeast <- sampling(fitness_model_path_yeast, 
                                                data = list(N = nrow(joined_df_flux_prot_mean_no_negative_subset),
                                                            K = length(unique(joined_df_flux_prot_mean_no_negative_subset$Gene)),
                                                            id = joined_df_flux_prot_mean_no_negative_subset$gene_as_factor_number_summary,
                                                            flux = joined_df_flux_prot_mean_no_negative_subset$z_score_unlog_flux,
                                                            protein_amounts = joined_df_flux_prot_mean_no_negative_subset$z_score_prot,
                                                            mean_protein_amount = mean_protein_amount_singles_summaries$log_mean_prot_mf_perc,
                                                            substrate_connectivity = mean_protein_amount_singles_summaries$geometric_mean,
                                                            delta_g = mean_protein_amount_singles_summaries$delta_g), 
                                                            iter = 4000, 
                                                            chains = 4)

## checking the model traces
mcmc_trace(fit_fitness_path_two_hierarch_yeast, 
           pars = c("gamma_0", "gamma_1", "gamma_2", "gamma_3", "tau", "sigma"))

# mcmc_trace(fit_fitness_path_two_hierarch_yeast)

gamma0_yeast <- fit_fitness_path_two_hierarch_yeast %>% 
  spread_draws(gamma_0) %>% 
  ggplot(aes(x = gamma_0)) +
  geom_histogram() +
  geom_vline(xintercept = 0) +
  theme_bw() +
  ylab('Sample Count') +
  xlab(expression(gamma[0])) +
  ggtitle('A) Intercept');gamma0_yeast

gamma1_yeast <- fit_fitness_path_two_hierarch_yeast %>% 
  spread_draws(gamma_1) %>% 
  ggplot(aes(x = gamma_1)) +
  geom_histogram() +
  geom_vline(xintercept = 0) +
  theme_bw() +
  ylab('Sample Count') +
  xlab(expression(gamma[1]~(log[e](Protein~Mass~Fraction*~"%")^-1))) +
  theme(text = element_text(size = 16)) +
  ggtitle('A. Protein cost hypothesis');gamma1_yeast

gamma2_yeast <- fit_fitness_path_two_hierarch_yeast %>% 
  spread_draws(gamma_2) %>% 
  ggplot(aes(x = gamma_2)) +
  geom_histogram() +
  geom_vline(xintercept = 0) +
  theme_bw() +
  ylab('Sample Count') +
  xlab(expression(gamma[2]~(Substrate~Connectivity~Metric^-1))) +
  theme(text = element_text(size = 16)) +
  ggtitle('B. Substrate connectivity hypothesis');gamma2_yeast

gamma3_yeast <- fit_fitness_path_two_hierarch_yeast %>% 
  spread_draws(gamma_3) %>% 
  ggplot(aes(x = gamma_3))  +
  geom_histogram() +
  geom_vline(xintercept = 0) +
  theme_bw() +
  ylab('Sample Count')  +
  xlab(expression(gamma[3]~((kJ~Mol^-1)^-1))) +
  theme(text = element_text(size = 16)) +
  ggtitle(expression(paste('C. ', Delta, 'G hypothesis')));gamma3_yeast

yeast_2_level_out <- ggarrange(gamma1_yeast, 
                               gamma2_yeast, 
                               gamma3_yeast, nrow = 1, align = 'hv')

ggsave(yeast_2_level_out, filename = 'figures/yeast_2_level_out.pdf', 
       height = 5, width = 16)

