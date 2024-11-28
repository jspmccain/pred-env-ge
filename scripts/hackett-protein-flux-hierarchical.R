
# Partial pooling hierarchical model -- fit and par estimates ---------------

fitness_model_partial_yeast <- stan_model("scripts/flux_protein_hierarchical_partial_pooling.stan")

options(mc.cores = 4)
number_of_proteins <- joined_df_flux_prot_mean_no_negative$gene_as_factor_number %>% unique() %>% length()

## Generating fake data for posterior predictions.
protein_amounts_spaced_sampled <- rep(seq(from = -2, to = 2, by = 0.1), number_of_proteins)
protein_ids_spaced_sampled <- rep(c(1:number_of_proteins), 
                                  each = length(seq(from = -2, to = 2, by = 0.1)))

## fitting the model in stan
fit_fitness_partial_pool <- sampling(fitness_model_partial_yeast, 
                                     data = list(N = nrow(joined_df_flux_prot_mean_no_negative),
                                                 K = length(unique(joined_df_flux_prot_mean_no_negative$Gene)),
                                                 id = joined_df_flux_prot_mean_no_negative$gene_as_factor_number,
                                                 flux = joined_df_flux_prot_mean_no_negative$z_score_unlog_flux,
                                                 protein_amounts = joined_df_flux_prot_mean_no_negative$z_score_prot,
                                                 protein_amounts_spaced = protein_amounts_spaced_sampled,
                                                 N_spaced = length(protein_amounts_spaced_sampled),
                                                 id_spaced = protein_ids_spaced_sampled), 
                                     iter = 3000, 
                                     chains = 4)

fit_hackett_out <- rstan::extract(fit_fitness_partial_pool, 
                                  permuted = TRUE)

beta_1_estimates_hackett_df <- fit_fitness_partial_pool %>% 
  spread_draws(beta_1[n]) %>% 
  median_qi() %>% 
  dplyr::rename(gene_as_factor_number = n) %>% 
  inner_join(mean_protein_amount_singles_summaries,
             by = 'gene_as_factor_number')

slope_estimates_yeast <- fit_fitness_partial_pool %>% 
  spread_draws(beta_1[n]) %>% 
  median_qi() %>% 
  ggplot(aes(y = fct_reorder(as.factor(n), beta_1), 
             x = beta_1, 
             xmin = .lower, 
             xmax = .upper)) +
  geom_vline(xintercept = -1, 
             linetype = 2) +
  geom_vline(xintercept = 1, 
             linetype = 2) +
  geom_pointinterval(alpha = 0.2) +
  geom_vline(xintercept = 0) +
  labs(x = 'Slope Estimate between Normalized Enzyme Abundance and Normalized Flux',
       y = 'Enzyme') +
  theme_bw() +
  theme(axis.text.y = element_blank());slope_estimates_yeast

# posterior predictive check from generated quantities --------------------

posterior_pred_hackett_obs <- fit_fitness_partial_pool %>% 
  spread_draws(y_rep[n]) %>% 
  median_qi()

# making a dataframe that has the posterior predictive checks in there
joined_df_flux_prot_mean_no_negative_w_pp <- joined_df_flux_prot_mean_no_negative
joined_df_flux_prot_mean_no_negative_w_pp$y_rep <- posterior_pred_hackett_obs$y_rep
joined_df_flux_prot_mean_no_negative_w_pp$lower_val <- posterior_pred_hackett_obs$.lower
joined_df_flux_prot_mean_no_negative_w_pp$upper_val <- posterior_pred_hackett_obs$.upper

posterior_predictive_check_hackett_single_hierarchy <- joined_df_flux_prot_mean_no_negative_w_pp %>% 
  # filter(Gene == gene_of_interest) %>%
  ggplot(aes(x = z_score_prot, y = y_rep)) +
  geom_ribbon(aes(ymin = lower_val, ymax = upper_val), alpha = 0.3) +
  facet_wrap(~Gene) +
  # geom_line() +
  geom_line() +
  geom_point(data = joined_df_flux_prot_mean_no_negative_w_pp, 
             aes(x = z_score_prot, y = z_score_flux), alpha = 0.5) +
  theme_bw() +
  ylab('Z-Score Normalized Rate') +
  xlab('Z-Score Normalized Protein')

ggsave(filename = "posterior_predictive_check_hackett_single_hierarchy.png", 
       plot = posterior_predictive_check_hackett_single_hierarchy, 
       height = 7.92, 
       width = 11.4)

