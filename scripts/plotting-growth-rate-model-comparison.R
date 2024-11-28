### looking at growth rate correlations between the different fluxes

### read in Hackett supp data with dilution rates
hackett_dilution_rates <- read_excel(path = 'data/aaf2786-hackett-sm-table-s9.xlsx', 
                                     sheet = 2) %>% 
  dplyr::select(DR_Actual, ChemostatID) %>% 
  dplyr::rename(exp_condition = ChemostatID)

## join this growth rate data with the individual fluxes
hackett_fluxes_with_dilution_rates <- hackett_fluxes_subset %>% 
  inner_join(hackett_dilution_rates, by = "exp_condition")

## calculate a reaction rate-specific correlation with growth rate
rate_cor_with_reaction_names <- hackett_fluxes_with_dilution_rates %>% 
  group_by(reaction_name) %>% 
  summarize(rate_growth_cor = cor(flux, DR_Actual))

## what does that distribution look like?
rate_cor_hist <- rate_cor_with_reaction_names %>% 
  ggplot(aes(x = rate_growth_cor)) +
  geom_histogram() +
  theme_bw() +
  xlab("Rate correlation with Dilution Rate")
  
## now compare the full proteome prediction models with the growth rate correlations
r_squared_vs_growth_rate_cor <- all_data_subsys_r_squared_ridge %>% 
  inner_join(rate_cor_with_reaction_names, by = "reaction_name") %>% 
  ggplot(aes(x = r_sq, y = rate_growth_cor)) +
  geom_point(alpha = 0.7) +
  theme_bw() +
  xlab(expression(Model~Performance~(Cross~Validated~R^2))) +
  ylab("Rate correlation with Dilution Rate")

succinate_transport_predictions <- all_yeast_fluxes_predicted_list_ridge[[1]] %>% 
  dplyr::filter(reaction_name == 'r_2057') %>% 
  group_by(flux_index_val, actual_flux) %>% 
  summarize(mean_pred = mean(pred_flux)) %>% 
  ggplot(aes(x = actual_flux*1e6, y = mean_pred*1e6)) +
  geom_point() +
  theme_bw() +
  coord_fixed(ratio = 1) +
  scale_y_continuous(limits = c(-85, 0)) +
  scale_x_continuous(limits = c(-85, 0)) +
  labs(x = 'Observed Rate of Succinate Transport<br>(nmol hour<sup>\u22121</sup> mL<sup>\u22121</sup>)',
       y = 'Predicted Rate of Succinate Transport<br>(nmol hour<sup>\u22121</sup> mL<sup>\u22121</sup>)') +
  theme(
    axis.title = ggtext::element_markdown()
  )

comparing_growth_rate_with_flux <- ggarrange(rate_cor_hist + ggtitle('A.') +
                                               theme(text = element_text(size = 13)), 
                                             r_squared_vs_growth_rate_cor + 
                                               ggtitle('B.') +
                                               theme(text = element_text(size = 13)), 
                                             succinate_transport_predictions + 
                                               ggtitle('C. ') +
                                               theme(text = element_text(size = 13)),
                                             nrow = 2, ncol = 2, align = 'hv')

ggsave(comparing_growth_rate_with_flux, 
       filename = paste0("figures/figsX_", current_date, "_growth_rate_correlation_with_whole_proteome_ridge_models.pdf"),
       height = 9.5, 
       width = 9.5)

hackett_fluxes_with_dilution_rates %>% 
  dplyr::filter(reaction_name == 'r_2057') %>% 
  dplyr::mutate(growth_media = substr(exp_condition, start=1, stop=1)) %>% 
  ggplot(aes(x = DR_Actual, y = flux)) + 
  geom_point(aes(shape = growth_media), size = 3) +
  xlab('Growth rate') +
  ylab('Succinate transport rate')

