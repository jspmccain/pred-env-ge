## plotting the number of non-zero coefficients in lasso regression models for Hackett et al full proteome predictions

# This script requires running hackett-full-proteome-prediction.R to get all_yeast_fluxes_predicted_list_summaries

mean_coef_full_proteome_lasso <- all_yeast_fluxes_predicted_list_summaries %>% 
  group_by(reaction_name) %>% 
  summarize(mean_coef_number_per_reaction = mean(coef_number)) %>% 
  ggplot(aes(x = mean_coef_number_per_reaction)) +
  geom_histogram() +
  theme_bw() +
  xlab('Mean Number of Non-zero Coefficients')

ggsave(mean_coef_full_proteome_lasso, 
       filename = paste0('figures/figsX_', current_date, '_mean_coef_full_proteome_lasso.pdf'), 
       height = 5.6, width = 7)
