#### making some figures that shows the relationship between rate and enzyme using a Bayesian hierarchical model

library(truncnorm)

make_poster_pred <- function(model_output_stan, n_samples = 500000){
  ### function that inputs the partial pooling model estimates and makes a figure that
  ### is a posterior predictive check
  
  # model_output_stan <- davidi_fitness_partial_pool
  tau_vec <- model_output_stan %>% 
    spread_draws(tau)
  tau_vec_1 <- tau_vec$tau
  
  gamma_vec <- model_output_stan %>% 
    spread_draws(gamma_1)
  gamma_vec_1 <- gamma_vec$gamma_1
  
  posterior_draws_out <- c()
  
  for(i in 1:n_samples){
    sampling_index <- sample(x = length(tau_vec_1), size = 1)
    sampled_tau <- tau_vec_1[sampling_index]
    sampled_gamma <- gamma_vec_1[sampling_index]
    sampled_output <- rtruncnorm(n = 1, 
                                mean = sampled_gamma, 
                                sd = sampled_tau, 
                                a = -1, b = 1)
    posterior_draws_out <- c(posterior_draws_out, sampled_output)
  }
  return(data.frame(posterior_draws_out))
}

yeast_col <- 'darkseagreen4'
bsub_col <- 'lightgoldenrod1'
ecoli_col <- 'steelblue4'

## Making figures of individual enzyme rate pairs
beta_1_estimates_davidi <- davidi_fitness_partial_pool %>% 
  spread_draws(beta_1[n]) %>% 
  median_qi() %>% 
  ggplot(aes(y = fct_reorder(as.factor(n), beta_1), 
             x = beta_1, 
             xmin = .lower, 
             xmax = .upper)) +
  geom_vline(xintercept = -1, 
             linetype = 2) +
  ggtitle(expression(paste('B. ', italic(~E.~coli)))) +
  geom_vline(xintercept = 1, 
             linetype = 2) +
  geom_pointinterval(colour = ecoli_col, alpha = 0.7) +
  geom_vline(xintercept = 0) +
  labs(x = 'Pooled Correlation Coefficient between Normalized Enzyme \nAbundance and Normalized Rate',
       y = 'Enzyme') +
  theme_bw() +
  theme(axis.text.y = element_blank());beta_1_estimates_davidi

beta_1_estimates_yeast <- fit_fitness_partial_pool %>% 
  spread_draws(beta_1[n]) %>% 
  median_qi() %>% 
  ggplot(aes(y = fct_reorder(as.factor(n), beta_1), 
             x = beta_1, 
             xmin = .lower, 
             xmax = .upper)) +
  geom_vline(xintercept = -1, 
             linetype = 2) +
  ggtitle(expression(paste('A.', italic(~S.~cerevisiae)))) +
  geom_vline(xintercept = 1, 
             linetype = 2) +
  geom_pointinterval(colour = yeast_col, alpha = 0.7) +
  geom_vline(xintercept = 0) +
  labs(x = 'Pooled Correlation Coefficient between Normalized Enzyme \nAbundance and Normalized Rate',
       y = 'Enzyme') +
  theme_bw() +
  theme(axis.text.y = element_blank());beta_1_estimates_yeast

b_subtilis_out_beta1 <- chu_partial_pool %>% 
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
  geom_pointinterval(colour = bsub_col, alpha = 0.7) +
  geom_vline(xintercept = 0) +
  labs(x = 'Pooled Correlation Coefficient between Normalized\nTranscript Abundance and Normalized Rate',
       y = 'Transcript') +
  ggtitle(expression(paste('C. ', italic(~B.~subtilis)))) +
  theme_bw() +
  theme(axis.text.y = element_blank());b_subtilis_out_beta1

correlation_coefficients_three_species <- ggarrange(beta_1_estimates_yeast, 
                                                    beta_1_estimates_davidi,
                                                    b_subtilis_out_beta1,
          align = 'hv', nrow = 1)

# ggsave(correlation_coefficients_three_species, 
#        filename = 'figures/cor_between_rate_protein_yeast_e_coli_b_sub.pdf', 
#        width = 12, height = 8)

### getting the posterior predictions for correlations

posterior_pred_out_e_coli <- make_poster_pred(model_output_stan = davidi_fitness_partial_pool)
posterior_pred_out_yeast <- make_poster_pred(model_output_stan = fit_fitness_partial_pool)
posterior_pred_out_bsub <- make_poster_pred(model_output_stan = chu_partial_pool)

posterior_pred_e_coli_fig <- posterior_pred_out_e_coli %>% 
  ggplot(aes(x = posterior_draws_out)) +
  geom_histogram(fill = ecoli_col) +
  theme_bw() +
  ggtitle(expression(paste('A. ', italic(E.~coli)))) +
  xlim(-1, 1.00000001) +
  ylab('Density') +
  geom_vline(xintercept = 0) +
  xlab('Distribution of Pooled Correlation Coefficients between \nNormalized Enzyme Abundance and Normalized Rate');posterior_pred_e_coli_fig

posterior_pred_yeast_fig <- posterior_pred_out_yeast %>% 
  ggplot(aes(x = posterior_draws_out)) +
  geom_histogram(fill = yeast_col) +
  theme_bw() +
  xlim(-1, 1.00000001) +
  ggtitle(expression(paste('B. ', italic(S.~cerevisiae)))) +
  ylab('Density') +
  geom_vline(xintercept = 0) +
  xlab('Distribution of Pooled Correlation Coefficients between \nNormalized Enzyme Abundance and Normalized Rate');posterior_pred_yeast_fig

posterior_pred_bsub_fig <- posterior_pred_out_bsub %>% 
  ggplot(aes(x = posterior_draws_out)) +
  geom_histogram(fill = bsub_col) +
  theme_bw() +
  xlim(-1, 1.00000001) +
  ggtitle(expression(paste('C. ', italic(B.~subtilis)))) +
  ylab('Density') +
  geom_vline(xintercept = 0) +
  xlab('Distribution of Pooled Correlation Coefficients between \nNormalized Transcript Abundance and Normalized Rate');posterior_pred_bsub_fig

combined_single_protein_to_rate <- ggarrange(beta_1_estimates_yeast, 
                                             beta_1_estimates_davidi,
                                             b_subtilis_out_beta1,
                                             posterior_pred_yeast_fig,
                                             posterior_pred_e_coli_fig,
                                             posterior_pred_bsub_fig, align = 'hv',
                                             nrow = 2, ncol = 3)
# ggsave(combined_single_protein_to_rate, 
#        filename = 'figures/single_protein_to_rate_par_ppd.pdf',
#        width = 17*0.75, height = 10*0.7)


# replotting group level distributions ------------------------------------
posterior_pred_out_bsub$organism <- rep("B. subtilis", nrow(posterior_pred_out_bsub))
posterior_pred_out_e_coli$organism <- rep("E. coli", nrow(posterior_pred_out_e_coli))
posterior_pred_out_yeast$organism <- rep('S. cerevisiae', nrow(posterior_pred_out_yeast))

posterior_pred_all <- rbind(posterior_pred_out_bsub,
                            posterior_pred_out_e_coli,
                            posterior_pred_out_yeast)
