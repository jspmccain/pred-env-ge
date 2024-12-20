### make a figure of glutathione metabolism predictions

## Read in the standard names that map to the systematic names for the GSH subsystem
glutathione_standard_names <- read.csv('data/glutathione_standard_names_yeast.csv')

## Read in the interpretable names that were scraped from the yeast genome fasta file.
source('scripts/hackett-protein-name-formatting.R')

## retreiving the R squared from this single protein prediction
glutathione_r_squ_linear_reg <- singles_only_r_squared_yeast_rmsd_gsh2 %>% 
  dplyr::filter(reaction_rate_id == 'r_0485')
glutathione_r_squ_number_single_prot <- make_rsq_pretty(glutathione_r_squ_linear_reg$r_sq)

## plotting prediction for single protein linear regression
single_protein_prediction_gsh2 <- gsh2_z_score_transformed_model %>% 
  dplyr::filter(reaction_rate_id == 'r_0485') %>%
  group_by(actual_fluxes) %>% 
  summarize(predicted_flux_mean = mean(predicted_fluxes)) %>% 
  ggplot(aes(y = predicted_flux_mean*1e6, 
             x = actual_fluxes*1e6)) +
  geom_point() +
  theme_bw() +
  labs(x = 'Observed Rate of GSH Synthesis<br>(nmol hour<sup>\u22121</sup> mL<sup>\u22121</sup>)',
       y = 'Predicted Rate of GSH Synthesis<br>(nmol hour<sup>\u22121</sup> mL<sup>\u22121</sup>)') +
  theme(
    axis.title = ggtext::element_markdown()
  ) +
  xlim(0, 1) +
  ylim(0, 1) +
  geom_abline(xintercept = 0, slope = 1) +
  coord_fixed(ratio = 1) +
  annotate(geom = 'text', x = 5e-7*1e6, y = 2.5e-7*1e6, 
           label = glutathione_r_squ_number_single_prot, 
           parse = TRUE, size = 4.5)

## getting r squared for single protein transformed
ag_metrics_across_glut <- ag_metrics_across[ag_metrics_across$gene_name == 'YOL049W',]$spearman_cor %>% round(3)
ag_metrics_across_glut_r2 <- ag_metrics_across[ag_metrics_across$gene_name == 'YOL049W',]$spearman_cor^2 %>% round(3)

## plotting single protein relationship
gsh2_plot <- joined_df_flux_prot_mean_no_negative %>% 
  filter(Gene == 'YOL049W') %>%
  ggplot(aes(x = z_score_prot, y = z_score_unlog_flux)) +
  ggtitle('B. Glutathione synthetase (GSH2)') +
  # ggtitle('c) Monomeric glyoxalase I\n   (Glycolysis/Glutathione Biosynthesis)') +
  geom_point(fill = yeast_col, pch = 21, colour = 'black', size = 3) +
  theme_bw() +
  xlim(-2.5, 2.5) +
  ylim(-2.5, 2.1) +
  annotate(geom = 'text', x = -0.65, 2, 
           label = bquote(rho == .(ag_metrics_across_glut)),
           size = 4.5) +
  annotate(geom = 'text', x = -0.65, 1.7,
           label = bquote(rho^2 == .(ag_metrics_across_glut_r2)),
           size = 4.5) +
  # geom_smooth(method = 'lm', colour = 'darkblue') +
  theme(title = element_text(size = 13),
        strip.text.x = element_text(size = 13),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 13)) +
  xlab('log(Protein abundance)\n(Z-score transformed)') +
  ylab('Reaction rate\n(Z-score transformed)') +
  coord_equal(ratio = 1);gsh2_plot

gamma_glut_predictions <- all_yeast_fluxes_predicted_list_subsystem_ridge_z[[1]] %>% 
  ## first need to summarize the mean and sd of predictions when a given observation was left out
  group_by(reaction_name, flux_index_val) %>% 
  summarize(mean_pred_flux = mean(pred_flux),
            sd_pred_flux = sd(pred_flux),
            mean_actual_flux = mean(actual_flux, na.rm = TRUE), # these are all identical from a single index
            sd_actual_flux = sd(actual_flux), # this is just a check, it should be zero
            n_pre_flux = n())

## retreiving the R squared from this subsystem prediction
glutathione_r_squ <- subsys_only_r_squared_ridge_z %>% 
  dplyr::filter(reaction_name == 'r_0485')
glutathione_r_squ_number <- make_rsq_pretty(glutathione_r_squ$r_sq)

glutathione_metabolism_prediction <- gamma_glut_predictions %>% 
  dplyr::filter(reaction_name == 'r_0485') %>% 
  ggplot(aes(x = mean_actual_flux*1e6, y = mean_pred_flux*1e6)) +
  geom_point(size = 3, colour = "grey30") +
  theme_bw() +
  ggtitle('C. Subsystem protein predictions') +
  coord_fixed(ratio = 1) +
  scale_y_continuous(limits = c(0, 9e-7*1e6)) +
  scale_x_continuous(limits = c(0, 9e-7*1e6)) +
  annotate(geom = 'text', x = 5e-7*1e6, y = 2.5e-7*1e6, 
           label = glutathione_r_squ_number, 
           parse = TRUE, size = 4.5) +
  theme(title = element_text(size = 13),
        strip.text.x = element_text(size = 13),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 13)) +
  labs(x = 'Observed Rate of GSH Synthesis<br>(\u03bcmol hour<sup>\u22121</sup> mL<sup>\u22121</sup>)',
       y = 'Predicted Rate of GSH Synthesis<br>(\u03bcmol hour<sup>\u22121</sup> mL<sup>\u22121</sup>)') +
  theme(
    axis.title = ggtext::element_markdown()
  ) +
  geom_abline(intercept = 0, slope = 1);glutathione_metabolism_prediction

hackett_name_mappings_with_int <- rbind(hackett_name_mappings,
                                        data.frame(hackett_observed_prot_names = '(Intercept)',
                                                   long_names = '(Intercept)'))

gamma_glut_coefs <- all_yeast_fluxes_predicted_list_subsystem_ridge_z[[3]] %>% 
  dplyr::filter(reaction_name == 'r_0485') %>% 
  left_join(hackett_name_mappings_with_int %>% 
               dplyr::rename(coefficient_names = hackett_observed_prot_names), 
             by = 'coefficient_names') %>% 
  left_join(glutathione_standard_names %>% 
              dplyr::rename(coefficient_names = systematic_name),
              by = "coefficient_names")

coefficient_gamma_glut <- gamma_glut_coefs %>% 
  dplyr::filter(long_names != '(Intercept)') %>% 
  group_by(standard_name) %>% 
  summarize(mean_coef_value = mean(coefficient_values),
            sd_coef_value = sd(coefficient_values)) %>% 
  # dplyr::mutate(new_name = paste(long_names, ' (', coefficient_names, ')', sep = '')) %>% 
  ggplot(aes(x = mean_coef_value*1e6, 
             y = fct_reorder(standard_name, mean_coef_value, .desc = FALSE))) +
  geom_point(aes(colour = standard_name), size = 3) +
  ylab('Protein') +
  ggtitle('D. Ridge regression coefficients') +
  scale_x_continuous(label=scientific_10) +
  labs(x = expression(paste('Coefficient Value (\u03bcmol', ' ', hour^-1, mL^-1, ')'))) +
  geom_errorbarh(aes(xmin = mean_coef_value*1e6 - sd_coef_value*1e6, 
                     xmax = mean_coef_value*1e6 + sd_coef_value*1e6,
                     colour = standard_name), height = 0.7, lwd = 1.1) +
  theme_bw() +
  theme(plot.title.position = "plot", #NEW parameter. Apply for subtitle too.
        plot.caption.position =  "plot",
        plot.margin = unit(c(0,0.2,0,1), 'lines'),
        legend.position = 'none',
        title = element_text(size = 13),
        strip.text.x = element_text(size = 13),
        axis.title = element_text(size = 13),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 10)) +
  scale_colour_manual(values = c(rep('black', 7), 'coral3', rep('black', 8))) +
  # theme(axis.text.x = element_text(hjust = 1, angle = 60)) +
  geom_vline(xintercept = 0);coefficient_gamma_glut

blank_plot <- ggplot() + ggtitle('A. Glutathione metabolism subsystem') + theme_bw() + 
  theme(rect = element_blank()) + 
  theme(plot.title.position = "plot", plot.caption.position =  "plot", 
        title = element_text(size = 13))


glutathione_out <- ggarrange(ggarrange(blank_plot, 
                                       gsh2_plot,
                                       glutathione_metabolism_prediction,
                                       nrow = 1, 
                                       align = 'hv'),
          coefficient_gamma_glut, nrow = 2, heights = c(2, 1.07))

## saving main glutathione predictions figure.
ggsave(glutathione_out, 
       filename = "figures/glutathione_out.svg",
       height = 9.57*0.8*0.9, width = 16*0.84)

## saving supplememntary figure of single protein predictions
ggsave(single_protein_prediction_gsh2,
       filename = paste0('figures/figsX_', current_date, '_gsh2_linear_model_predictions.svg'),
       height = 7, 
       width = 8)

# comparing ridge regression coefficient terms ----------------------------

getting_gsh_synth <- get_proteome_flux_dataframe(reaction_name_choice = "r_0485", 
                            gene_flux_data_frame = hackett_fluxes_subset, 
                            wide_proteome_data_frame = hackett_prot_transformed, org = "yeast")

scaled_prot_values <- scale(getting_gsh_synth %>% dplyr::select(-flux)) %>% 
  as.data.frame()

gsh_synth_scaled <- cbind(data.frame(flux = getting_gsh_synth$flux),
                          scaled_prot_values)

gsh_synth_scaled$exp_condition <- hackett_prot_transformed$exp_condition

gsh_synth_scaled$growth_rate <- as.numeric(str_sub(gsh_synth_scaled$exp_condition, end = 5, start = 3))
gsh_synth_scaled$nutrient <- str_sub(gsh_synth_scaled$exp_condition, start = 1, end = 1)

gsh_synth_melt <- gsh_synth_scaled %>% 
  melt(id.vars = c('exp_condition', 'nutrient', 'growth_rate', 'flux'), 
       value.name = 'z_prot', variable.name = 'protein')

gsh_synth_melt %>% 
  inner_join(glutathione_standard_names %>% 
               dplyr::rename(protein = systematic_name), 
             by = 'protein') %>% 
  ggplot(aes(x = z_prot, y = flux)) +
  geom_point(aes(colour = growth_rate, shape = nutrient),
             size = 3) +
  facet_wrap(~standard_name) +
  theme_bw() +
  xlab('Protein Abundance (z score transformed)')
