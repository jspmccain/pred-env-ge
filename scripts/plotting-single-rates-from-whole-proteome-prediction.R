### make a plot of single rate prediction and whole proteome rate prediction

## key data output needed is:

### all_yeast_fluxes_predicted_list_ridge

# Summarizing / formatting model output -----------------------------------

### summarizing across cross validation instances
summarized_df_for_ridge_whole_proteome <- all_yeast_fluxes_predicted_all_data_ridge %>% 
  ## first need to summarize the mean and sd of predictions when a given observation was left out
  group_by(reaction_name, flux_index_val) %>% 
  summarize(mean_pred_flux = mean(pred_flux),
            sd_pred_flux = sd(pred_flux),
            mean_actual_flux = mean(actual_flux), # these are all identical from a single index
            sd_actual_flux = sd(actual_flux), # this is just a check, it should be zero
            n_pre_flux = n())

### summarizing across cross validation instances
summarized_df_for_ridge_whole_proteome_block_cv <- all_yeast_fluxes_predicted_list_ridge_group_cv[[1]] %>% 
  ## first need to summarize the mean and sd of predictions when a given observation was left out
  group_by(reaction_name, flux_index_val) %>% 
  summarize(mean_pred_flux = mean(pred_flux),
            sd_pred_flux = sd(pred_flux),
            mean_actual_flux = mean(actual_flux, na.rm = TRUE), # these are all identical from a single index
            sd_actual_flux = sd(actual_flux), # this is just a check, it should be zero
            n_pre_flux = n())

summarized_df_for_ridge_pathway_leave_two <- all_yeast_fluxes_predicted_list_subsystem_ridge[[1]] %>% 
  ## first need to summarize the mean and sd of predictions when a given observation was left out
  group_by(reaction_name, flux_index_val) %>% 
  summarize(mean_pred_flux = mean(pred_flux),
            sd_pred_flux = sd(pred_flux),
            mean_actual_flux = mean(actual_flux, na.rm = TRUE), # these are all identical from a single index
            sd_actual_flux = sd(actual_flux), # this is just a check, it should be zero
            n_pre_flux = n())

# getting model metrics to annotate figures -------------------------------

## aspartate kinase as a case study YER052C
aspartate_kinase_rxn_number = 'r_0215'

### retreiving metric to put on the figures
leave_2_out_r_0215_rmsd <- all_data_subsys_r_squared_ridge[all_data_subsys_r_squared_ridge$reaction_name == aspartate_kinase_rxn_number, ]$rmsd
leave_2_out_r_0215_rsq <- all_data_subsys_r_squared_ridge[all_data_subsys_r_squared_ridge$reaction_name == aspartate_kinase_rxn_number, ]$r_sq

leave_condition_out_r_0215_rmsd <- all_data_subsys_r_squared_ridge_group_cv[all_data_subsys_r_squared_ridge_group_cv$reaction_name == aspartate_kinase_rxn_number, ]$rmsd
leave_condition_out_r_0215_rsq <- all_data_subsys_r_squared_ridge_group_cv[all_data_subsys_r_squared_ridge_group_cv$reaction_name == aspartate_kinase_rxn_number, ]$r_sq

### retrieving metrics from the single protein
# leave_2_out_r_0215_rmsd <- all_data_subsys_r_squared_ridge[all_data_subsys_r_squared_ridge$reaction_name == aspartate_kinase_rxn_number, ]$rmsd
# leave_2_out_r_0215_rsq <- all_data_subsys_r_squared_ridge[all_data_subsys_r_squared_ridge$reaction_name == aspartate_kinase_rxn_number, ]$r_sq
# 
# leave_condition_out_r_0215_rmsd <- all_data_subsys_r_squared_ridge_group_cv[all_data_subsys_r_squared_ridge_group_cv$reaction_name == aspartate_kinase_rxn_number, ]$rmsd
# leave_condition_out_r_0215_rsq <- all_data_subsys_r_squared_ridge_group_cv[all_data_subsys_r_squared_ridge_group_cv$reaction_name == aspartate_kinase_rxn_number, ]$r_sq

### retrieving metrics from the pathway
leave_2_out_r_0215_rsq_pathway <- subsys_only_r_squared_ridge[subsys_only_r_squared_ridge$reaction_name == aspartate_kinase_rxn_number, ]$r_sq
leave_condition_out_r_0215_rsq_pathway <- subsys_only_r_squared_ridge_group_cv[subsys_only_r_squared_ridge_group_cv$reaction_name == aspartate_kinase_rxn_number, ]$r_sq


### make function for making the number pretty
make_rsq_pretty <- function(number_input){
  rounded_number <- round(number_input, 3)
  return(paste("R^2 == ", rounded_number))
}
make_rmsd_pretty <- function(number_input){
  rounded_number <- round(number_input, 8)
  return(paste("RMSD == ", rounded_number))
}

# plotting ----------------------------------------------------------------

### getting the summary statistics
rsq_label_formatted_two <- make_rsq_pretty(leave_2_out_r_0215_rsq)
rmsd_label_formatted_two <- make_rmsd_pretty(leave_2_out_r_0215_rmsd)

### getting the summary statistics for pathway level
rsq_label_formatted_two_pathway <- make_rsq_pretty(leave_2_out_r_0215_rsq_pathway)

##### plotting the 2-left out cross validation
aspartate_kinase_predictions_two <- summarized_df_for_ridge_whole_proteome %>% 
  dplyr::filter(reaction_name == aspartate_kinase_rxn_number) %>% 
  ## multiply by 1e6 to convert from mmol to nmol
  ggplot(aes(x = mean_actual_flux*1e+6, y = mean_pred_flux*1e+6)) +
  geom_point(size = 3, colour = "grey30") +
  theme_bw() +
  annotate(geom = 'text', 
           x = 2e-5*1e6, 
           y = 5e-6*1e6, 
           label = rsq_label_formatted_two,
           parse = TRUE, size = 6) +
  ggtitle('B.') +
  # scale_y_continuous(label=scientific_10, limits = c(0, 3e-5*1e6)) +
  # scale_x_continuous(label=scientific_10, limits = c(0, 3e-5*1e6)) +
  scale_y_continuous(limits = c(0, 3e-5*1e6)) +
  scale_x_continuous(limits = c(0, 3e-5*1e6)) +
  labs(x = 'Observed Rate of Aspartate Kinase Activity<br>(nmol hour<sup>\u22121</sup> mL<sup>\u22121</sup>)',
       y = 'Predicted Rate of Aspartate Kinase Activity<br>(nmol hour<sup>\u22121</sup> mL<sup>\u22121</sup>)') +
  theme(
    axis.title = ggtext::element_markdown(size = 13),
    plot.title = element_text(size = 13)) +
  # xlab(expression(paste('Observed Rate of Aspartate Kinase Activity (nmol ', hour^-1, mL^-1, ')'))) +
  # ylab(expression(paste('Predicted Rate of Aspartate Kinase Activity (nmol ', hour^-1, mL^-1, ')'))) +
  # theme(axis.title = element_text(size = 14), plot.title = element_text(size = 16)) +
  geom_abline(intercept = 0, slope = 1);aspartate_kinase_predictions_two

##### plotting the culture condition left out cross validation

rsq_label_formatted_group <- make_rsq_pretty(leave_condition_out_r_0215_rsq)
rmsd_label_formatted_group <- make_rmsd_pretty(leave_condition_out_r_0215_rmsd)

aspartate_kinase_predictions_group <- summarized_df_for_ridge_whole_proteome_block_cv %>% 
  dplyr::filter(reaction_name == 'r_0215') %>% 
  ggplot(aes(x = mean_actual_flux, y = mean_pred_flux)) +
  geom_point(size = 3, colour = "grey30") +
  theme_bw() +
  # geom_label_repel(aes(label = round(mean_actual_flux, 8))) +
  annotate(geom = 'text', 
           x = 2e-5, 
           y = 5e-6, 
           label = rsq_label_formatted_group,
           parse = TRUE, size = 6) +
  annotate(geom = 'text', 
           x = 2e-5, 
           y = 6.5e-6, 
           label = rmsd_label_formatted_group,
           parse = TRUE, size = 6) +
  ggtitle('B.') +
  xlab('Observed Rate of Aspartate Kinase Activity\n(mmol per hour per mL)') +
  ylab('Predicted Rate of Aspartate Kinase Activity\n(mmol per hour per mL)') +
  theme(axis.title = element_text(size = 14), plot.title = element_text(size = 16)) +
  geom_abline(intercept = 0, slope = 1);aspartate_kinase_predictions_group

# plotting single protein versus rate -------------------------------------

single_aspartate_kinase_df <- get_proteome_flux_dataframe(reaction_name_choice = aspartate_kinase_rxn_number, 
                                                          gene_flux_data_frame = hackett_fluxes_subset, 
                                                          wide_proteome_data_frame = hackett_prot_transformed, 
                                                          org = "yeast", 
                                                          filter_for_subsystem = 'single_prot')

aspartate_kinase_single_protein_rate <- single_aspartate_kinase_df %>% 
  ggplot(aes(x = YER052C, y = flux*1e6)) +
  geom_point(size = 3, colour = "grey30") +
  ggtitle('A.') +
  # ylab(expression(paste('Observed Rate of Aspartate Kinase Activity (nmol ', hour^-1, mL^-1, ')'))) +
  # ylab(expression(paste('Predicted Rate of GSH Synthesis (nmol ', hour^-1, mL^-1, ')'))) +
  theme_bw() +
  labs(y = 'Observed Rate of Aspartate Kinase Activity<br>(nmol hour<sup>\u22121</sup> mL<sup>\u22121</sup>)',
       x = 'Aspartate Kinase Relative Abundance (log<sub>\0032</sub>-transformed)') +
  theme(
    axis.title = ggtext::element_markdown(size = 13),
    plot.title = element_text(size = 13)
  ) +
  scale_y_continuous(limits = c(0, 3e-5*1e6));aspartate_kinase_single_protein_rate
  # scale_y_continuous(label=scientific_10, limits = c(0, 3e-5*1e6)) +
  # geom_label_repel(aes(label = round(flux, 8))) +
  # theme(axis.title = element_text(size = 14), plot.title = element_text(size = 16))


# arranging all plots -----------------------------------------------------

arranged_main_fig_asparate_predictions <- ggarrange(aspartate_kinase_single_protein_rate + 
                                                      ggtitle('A. Single protein-to-reaction rate relationship') +
                                                      theme(axis.text = element_text(size = 13)), 
          aspartate_kinase_predictions_two + 
            ggtitle('B. Proteome prediction') + 
            coord_fixed(ratio = 1) +
            theme(axis.text = element_text(size = 13)), 
          nrow = 1, align = 'hv')

ggsave(arranged_main_fig_asparate_predictions, 
       filename = 'figures/arranged_main_fig_asparate_predictions.pdf', 
       height = 6*0.7, width = 13.7*0.7)

## separately save them
ggsave(aspartate_kinase_single_protein_rate + ggtitle('A. Single protein-to-reaction rate relationship'),
       filename = 'figures/aspartate_kinase_single_protein_rate.pdf',
       height = 6, width = 13.7/2)

ggsave(aspartate_kinase_predictions_two + ggtitle('B. Proteome prediction'),
       filename = 'figures/aspartate_kinase_predictions_two.pdf',
       height = 6, width = 13.7/2)


# plotting aspartate kinase prediction using pathway level only -----------

pathway_prediction_plot_aspartate <- summarized_df_for_ridge_pathway_leave_two %>% 
  dplyr::filter(reaction_name == aspartate_kinase_rxn_number) %>% 
  ggplot(aes(x = mean_actual_flux*1e6, y = mean_pred_flux*1e6)) +
  geom_point(size = 3, colour = "grey30") +
  theme_bw() +
  # annotate(geom = 'text', 
  #          x = 2e-5, 
  #          y = 5e-6, 
  #          label = rsq_label_formatted_two,
  #          parse = TRUE, size = 6) +
  xlab('Observed Rate of Aspartate Kinase Activity\n(mmol per hour per mL)') +
  ylab('Predicted Rate of Aspartate Kinase Activity\n(mmol per hour per mL)') +
  theme(axis.title = element_text(size = 16), plot.title = element_text(size = 20)) +
  # scale_y_continuous(label=scientific_10, limits = c(-6e-6*1e6, 3e-5*1e6)) +
  scale_y_continuous(limits = c(-6e-6*1e6, 3e-5*1e6)) +
  scale_x_continuous(limits = c(-6e-6*1e6, 3e-5*1e6)) +
  # geom_errorbar(aes(ymin = mean_pred_flux - sd_pred_flux, ymax = mean_pred_flux + sd_pred_flux)) +
  geom_abline(intercept = 0, slope = 1) +
  coord_fixed(ratio = 1);pathway_prediction_plot_aspartate

arranged_main_fig_asparate_predictions_w_pathway <- ggarrange(aspartate_kinase_single_protein_rate, 
                                                              pathway_prediction_plot_aspartate + coord_fixed(ratio = 1),
                                                    aspartate_kinase_predictions_two + coord_fixed(ratio = 1), 
                                                    nrow = 1, align = 'hv')

ggsave(arranged_main_fig_asparate_predictions_w_pathway, 
       filename = 'figures/arranged_main_fig_asparate_predictions_w_pathway.pdf', 
       width = 22, height = 6.5)

# 