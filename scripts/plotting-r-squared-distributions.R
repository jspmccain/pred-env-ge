### plotting r-squared distributions and rmsd distributions and comparisons across scales

### best to run with the entire "running-all-analysis-plots.R", because the scripts required are:

## Data processing for main analyses
# source('scripts/davidi-data-processing.R')
# source('scripts/hackett-data-processing.R')
# source('scripts/chubukov-transcript-flux-hierarchical.R')
# source('scripts/davidi-protein-flux-hierarchical.R')
# source('scripts/hackett-protein-flux-hierarchical.R')
# source('scripts/hackett-full-proteome-prediction.R')
# source('scripts/single-protein-prediction.R')

library(ggforce)

# Calculate R-squared distributions for single protein-to-rate  --------

# Source the function for calculating R^2.
source("scripts/calculating_r_2_from_hierarchical_model.R")

# ecoli_col <- 'darkblue'
# yeast_col <- 'darkseagreen4'
# bsub_col <- 'yellow3'

label <- "Model Performance ("

# Calculating/plotting individual Bayesian R-squared distributions for each dataset --------

bayes_label <- expression(R["B"]^2)

hackett_r_squared <- bayes_r2_from_hierarchical(fit_fitness_partial_pool,
                           gene_as_factor_indices = joined_df_flux_prot_mean_no_negative$gene_as_factor_number)

davidi_r_squared <- bayes_r2_from_hierarchical(davidi_fitness_partial_pool,
                                               gene_as_factor_indices = genes_fluxes_long_df_z_score_cc$gene_as_number)

chubukov_r_squared <- bayes_r2_from_hierarchical(chu_partial_pool,
                                                 gene_as_factor_indices = chu_transformed_cc$flux_name_as_number)

## Making these distributions into a stacked figure

hackett_r_squared_single_hist <- hackett_r_squared %>% 
  ggplot(aes(x = r_sq)) +
  geom_histogram(fill = yeast_col) +
  theme_bw() +
  xlab(bayes_label) +
  xlim(0, 1) +
  ylab('Count') +
  ggtitle(expression(paste('A. Single proteins, ', italic(S.~cerevisiae))));hackett_r_squared_single_hist

davidi_r_squared_single_hist <- davidi_r_squared %>% 
  ggplot(aes(x = r_sq)) +
  geom_histogram(fill = ecoli_col) +
  theme_bw() +
  xlab(bayes_label)+
  ylab('Count') +
  ggtitle(expression(paste('B. Single proteins, ', italic(E.~coli)))) +
  xlim(0, 1);davidi_r_squared_single_hist

chubukov_r_squared_single_hist <- chubukov_r_squared %>% 
  ggplot(aes(x = r_sq)) +
  geom_histogram(fill = bsub_col) +
  theme_bw() +
  xlab(bayes_label) +
  ylab('Count') +
  ggtitle(expression(paste('C. Single transcripts, ', italic(B.~subtilis)))) +
  xlim(0, 1);chubukov_r_squared_single_hist

r_squared_hier_singles <- ggarrange(hackett_r_squared_single_hist,
                                    davidi_r_squared_single_hist,
                                    chubukov_r_squared_single_hist, 
                                    nrow = 1)

ggsave(r_squared_hier_singles, 
       filename = 'figures/r2_1_level_hier_model_protein_rates.png',
       width = 15, height = 6)

# Plotting individual R-squared distributions from the simple linear regressions (CV'd) --------

bsub_linear_reg_rsq <- singles_only_r_squared_b_sub_rmsd %>% 
  ggplot(aes(x = r_sq)) +
  geom_histogram(fill = bsub_col) +
  theme_bw() +
  xlab(expression(Model~Performance~(Cross~Validated~R^2))) +
  geom_vline(xintercept = 0, lty = 2) +
  ylab('Count') +
  ggtitle(expression(paste('B. Single transcripts, ', italic(B.~subtilis)))) +
  xlim(-1, 1);bsub_linear_reg_rsq

yeast_linear_reg_rsq <- singles_only_r_squared_yeast_rmsd %>% 
  ggplot(aes(x = r_sq)) +
  geom_histogram(fill = yeast_col) +
  theme_bw() +
  # xlab(expression(Model~Performance~(Cross~Validated~R^2))) +
  geom_vline(xintercept = 0, lty = 2) +
  xlim(-1, 1) +
  ylab('Count') +
  annotate(geom = 'text', x = -0.4, y  = 7.5, label = paste0('n = ', nrow(singles_only_r_squared_yeast_rmsd))) +
  theme(axis.title.x=element_blank()) +
  ggtitle(expression(paste('A. Single proteins, ', italic(S.~cerevisiae))));yeast_linear_reg_rsq

ecoli_linear_reg_rsq <- singles_only_r_squared_e_coli_rmsd %>% 
  ggplot(aes(x = r_sq)) +
  geom_histogram(fill = ecoli_col) +
  theme_bw() +
  xlab(expression(Model~Performance~(Cross~Validated~R^2)))+
  geom_vline(xintercept = 0, lty = 2) +
  ylab('Count') +
  ggtitle(expression(paste('A. Single proteins, ', italic(E.~coli)))) +
  xlim(-1, 1);ecoli_linear_reg_rsq

# Plotting individual RMSD distributions from the simple linear regressions (CV'd) --------

bsub_linear_reg_rmsd <- singles_only_r_squared_b_sub_rmsd %>% 
  ggplot(aes(x = rmsd)) +
  geom_histogram(fill = bsub_col) +
  theme_bw() +
  xlab("RMSD") +
  ylab('Count') +
  scale_x_log10() +
  ggtitle(expression(paste('C. Single proteins, ', italic(B.~subtilis))));bsub_linear_reg_rmsd

yeast_linear_reg_rmsd <- singles_only_r_squared_yeast_rmsd %>% 
  ggplot(aes(x = rmsd)) +
  geom_histogram(fill = yeast_col) +
  theme_bw() +
  xlab("RMSD") +
  ylab('Count') +
  scale_x_log10() +
  ggtitle(expression(paste('A. Single proteins, ', italic(S.~cerevisiae))));yeast_linear_reg_rmsd

ecoli_linear_reg_rmsd <- singles_only_r_squared_e_coli_rmsd %>% 
  ggplot(aes(x = rmsd)) +
  geom_histogram(fill = ecoli_col) +
  theme_bw() +
  xlab("RMSD") +
  ylab('Count') +
  scale_x_log10() +
  ggtitle(expression(paste('B. Single proteins, ', italic(E.~coli))));ecoli_linear_reg_rmsd

# LASSO Subsystem predictions only ----------------------------------------------

## make a histogram of the subsystem only predictions
subsys_only_r_squared_hist <- subsys_only_r_squared %>% 
  # ggplot(aes(x = cor_val)) +
  ggplot(aes(x = r_sq)) +
  geom_histogram(fill = yeast_col) +
  theme_bw() +
  xlab(expression(Model~Performance~(Cross~Validated~R^2))) +
  xlim(-1, 1) +
  ggtitle(expression(paste('D. Subsystem proteins only, ', italic(S.~cerevisiae)))) +
  ylab('Count');subsys_only_r_squared_hist

subsys_only_rmsd_hist <- subsys_only_r_squared %>% 
  # ggplot(aes(x = cor_val)) +
  ggplot(aes(x = rmsd)) +
  geom_histogram(fill = yeast_col) +
  theme_bw() +
  xlab("RMSD") +
  # xlim(1e-9, 1e-3) +
  # scale_x_log10() +
  scale_x_continuous(trans = "log", 
                     breaks = c(1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3), 
                     limits = c(1e-9, 1e-3)) +
  ggtitle(expression(paste('D. Subsystem proteins only, ', italic(S.~cerevisiae)))) +
  ylab('Count');subsys_only_rmsd_hist


# Ridge subsystem predictions only ----------------------------------------

subsys_only_rmsd_hist_ridge <- subsys_only_r_squared_ridge %>% 
  # ggplot(aes(x = cor_val)) +
  ggplot(aes(x = rmsd)) +
  geom_histogram(fill = yeast_col) +
  theme_bw() +
  xlab("RMSD") +
  # xlim(1e-9, 1e-3) +
  # scale_x_log10() +
  scale_x_continuous(trans = "log", 
                     breaks = c(1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3), 
                     limits = c(1e-9, 1e-3)) +
  ggtitle(expression(paste('D. Subsystem proteins only, ', italic(S.~cerevisiae)))) +
  ylab('Count');subsys_only_rmsd_hist_ridge

subsys_only_r_sq_hist_ridge <- subsys_only_r_squared_ridge %>% 
  # ggplot(aes(x = cor_val)) +
  ggplot(aes(x = r_sq)) +
  geom_histogram(fill = yeast_col) +
  theme_bw() +
  xlab(expression(Model~Performance~(Cross~Validated~R^2))) +
  geom_vline(xintercept = 0, lty = 2) +
  xlim(-1, 1) +
  annotate(geom = 'text', x = -0.4, y = 22, label = paste0('n = ', nrow(subsys_only_r_squared_ridge))) +
  theme(axis.title.x=element_blank()) +
  ggtitle(expression(paste('D. Subsystem proteins only, ', italic(S.~cerevisiae)))) +
  ylab('Count');subsys_only_r_sq_hist_ridge

subsys_only_r_sq_hist_ridge_single_subset <- subsys_only_r_squared_ridge %>% 
  ## filter for only those reactions found in the single protein comparison
  dplyr::filter(reaction_name %in% singles_only_r_squared_yeast_rmsd$reaction_rate_id) %>% 
  # ggplot(aes(x = cor_val)) +
  ggplot(aes(x = r_sq)) +
  geom_histogram(fill = yeast_col) +
  theme_bw() +
  xlab(expression(Model~Performance~(Cross~Validated~R^2))) +
  geom_vline(xintercept = 0, lty = 2) +
  xlim(-1, 1) +
  ggtitle(expression(paste('A. Subsystem proteins only, filtered subset, ', italic(S.~cerevisiae)))) +
  ylab('Count');subsys_only_r_sq_hist_ridge_single_subset

# LASSO Whole proteome predictions  ----------------------------------------------

all_proteome_data_r_squared_hist <- all_data_subsys_r_squared %>% 
  # ggplot(aes(x = cor_val)) +
  ggplot(aes(x = r_sq)) +
  geom_histogram(fill = yeast_col) +
  theme_bw() +
  ggtitle(expression(paste('E. Proteome, ', italic(S.~cerevisiae)))) +
  xlab(label) +
  geom_vline(xintercept = 0, lty = 2) +
  xlim(-1, 1) +
  ylab('Count');all_proteome_data_r_squared_hist

all_proteome_data_rmsd_hist <- all_data_subsys_r_squared %>% 
  # ggplot(aes(x = cor_val)) +
  dplyr::filter(reaction_name %in% single_protein_pairs_rxns) %>% 
  ggplot(aes(x = rmsd)) +
  geom_histogram(fill = yeast_col) +
  theme_bw() +
  scale_x_continuous(trans = "log", 
                     breaks = c(1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3), 
                     limits = c(1e-9, 1e-3)) +
  ggtitle(expression(paste('E. Proteome, ', italic(S.~cerevisiae)))) +
  xlab("RMSD") +
  ylab('Count');all_proteome_data_rmsd_hist

# Ridge Whole proteome predictions  ----------------------------------------------

all_proteome_data_r_squared_hist_ridge <- all_data_subsys_r_squared_ridge %>% 
  # ggplot(aes(x = cor_val)) +
  ggplot(aes(x = r_sq)) +
  geom_histogram(fill = yeast_col) +
  theme_bw() +
  ggtitle(expression(paste('E. Proteome, ', italic(S.~cerevisiae)))) +
  xlab(expression(Model~Performance~(Cross~Validated~R^2))) +
  geom_vline(xintercept = 0, lty = 2) +
  annotate(geom = 'text', x = -0.4, y = 60, label = paste0('n = ', nrow(all_data_subsys_r_squared_ridge))) +
  xlim(-1, 1) +
  ylab('Count');all_proteome_data_r_squared_hist_ridge

all_proteome_data_rmsd_hist_ridge <- all_data_subsys_r_squared_ridge %>% 
  # ggplot(aes(x = cor_val)) +
  dplyr::filter(reaction_name %in% single_protein_pairs_rxns) %>% 
  ggplot(aes(x = rmsd)) +
  geom_histogram(fill = yeast_col) +
  theme_bw() +
  scale_x_continuous(trans = "log", 
                     breaks = c(1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3), 
                     limits = c(1e-9, 1e-3)) +
  ggtitle(expression(paste('E. Proteome, ', italic(S.~cerevisiae)))) +
  xlab("RMSD") +
  ylab('Count');all_proteome_data_rmsd_hist_ridge

all_proteome_data_r_squared_hist_ridge_single_subset <- all_data_subsys_r_squared_ridge %>% 
  dplyr::filter(reaction_name %in% singles_only_r_squared_yeast_rmsd$reaction_rate_id) %>% 
  # ggplot(aes(x = cor_val)) +
  ggplot(aes(x = r_sq)) +
  geom_histogram(fill = yeast_col) +
  theme_bw() +
  ggtitle(expression(paste('B. Proteome, filtered subset, ', italic(S.~cerevisiae)))) +
  xlab(expression(Model~Performance~(Cross~Validated~R^2))) +
  geom_vline(xintercept = 0, lty = 2) +
  xlim(-1, 1) +
  ylab('Count');all_proteome_data_r_squared_hist_ridge_single_subset

# Comparing ridge regression and lasso predictions ------------------------

## aggregating model summaries from the lasso and ridge -- whole proteomes
all_data_r_squared_comparison <- all_data_subsys_r_squared_ridge %>% 
  dplyr::rename(r_sq_ridge = r_sq, rmsd_ridge = rmsd) %>% 
  inner_join(all_data_subsys_r_squared %>% ## appending the LASSO model
               dplyr::rename(r_sq_lasso = r_sq, rmsd_lasso = rmsd), 
             by = "reaction_name")

## aggregating model summaries from the lasso and ridge -- subsystem predictions
subsys_only_r_squared_comparison <- subsys_only_r_squared_ridge %>% 
  dplyr::rename(r_sq_ridge = r_sq, rmsd_ridge = rmsd) %>% 
  inner_join(subsys_only_r_squared %>%  ## appending the LASSO model output
               dplyr::rename(r_sq_lasso = r_sq, rmsd_lasso = rmsd), 
             by = "reaction_name")

## making an axis label
ridge_label <- expression(Model~Performance~(Cross~Validated~R^2)~Ridge)
lasso_label <- expression(Model~Performance~(Cross~Validated~R^2)~LASSO)

### plotting r squared comparisons
whole_prot_model_comparison <- all_data_r_squared_comparison %>% 
  ggplot(aes(x = r_sq_ridge, y = r_sq_lasso)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0) +
  theme_bw() +
  ggtitle("A. Proteome models") +
  xlab(ridge_label) +
  ylab(lasso_label)

subsys_prot_model_comparison <- subsys_only_r_squared_comparison %>% 
  ggplot(aes(x = r_sq_ridge, y = r_sq_lasso)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0) +
  theme_bw() +
  ggtitle("B. Subsystem proteins only") +
  xlab(ridge_label) +
  ylab(lasso_label)

whole_prot_model_comparison_rmsd <- all_data_r_squared_comparison %>% 
  ggplot(aes(x = rmsd_ridge, y = rmsd_lasso)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0) +
  theme_bw() +
  ggtitle("C. Proteome models") +
  xlab("RMSD Ridge") +
  ylab('RMSD LASSO') +
  scale_x_log10() +
  scale_y_log10()

subsys_prot_model_comparison_rmsd <- subsys_only_r_squared_comparison %>% 
  ggplot(aes(x = rmsd_ridge, y = rmsd_lasso)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0) +
  theme_bw() +
  ggtitle("D. Subsystem proteins only") +
  xlab("RMSD Ridge") +
  ylab('RMSD LASSO') +
  scale_y_log10() +
  scale_x_log10()

# Comparing RMSD from single proteins to whole proteome estimates ---------

combined_yeast_rmsd_across_scales <- singles_only_r_squared_yeast_rmsd %>% 
  dplyr::rename(reaction_name = reaction_rate_id,
                single_rmsd = rmsd) %>% 
  ## Note that there are only 45 in this final set because of this step. The subsystem prediction
  ## has a cutoff of requiring 5 proteins per pathway in order to make the prediction.
  inner_join(subsys_only_r_squared_ridge %>%
               dplyr::rename(subsys_rmsd = rmsd), by = 'reaction_name') %>%
  inner_join(all_data_subsys_r_squared_ridge %>% 
               dplyr::rename(full_rmsd = rmsd), by = 'reaction_name')

single_to_full_rmsd <- combined_yeast_rmsd_across_scales %>% 
  ggplot(aes(single_rmsd - full_rmsd)) +
  geom_histogram(fill = yeast_col) +
  scale_x_log10() +
  theme_bw() +
  xlab('RMSD Difference between Single Protein Prediction and Full Proteome Prediction') +
  ggtitle('A. Single Protein Prediction RMSD - Full Proteome RMSD')

subsystem_to_full_rmsd <- combined_yeast_rmsd_across_scales %>% 
  ggplot(aes(subsys_rmsd - full_rmsd)) +
  geom_histogram(fill = yeast_col) +
  scale_x_log10() +
  theme_bw() +
  xlab('RMSD Difference between Subsystem Prediction and Full Proteome Prediction') +
  ggtitle('B. Subsystem Prediction RMSD - Full Proteome RMSD')

# quantile(combined_yeast_rmsd_across_scales$subsys_rmsd - combined_yeast_rmsd_across_scales$full_rmsd)
# quantile(combined_yeast_rmsd_across_scales$single_rmsd - combined_yeast_rmsd_across_scales$full_rmsd)

# Plotting blocked CV for pathways and whole proteomes using ridge --------------------------------

## Pathway level predictions blocked cross validation

block_cross_validation_subsystem_ridge <- subsys_only_r_squared_ridge_group_cv %>% 
  ggplot(aes(x = r_sq)) +
  geom_histogram(fill = yeast_col, bins = 100) +
  theme_bw() +
  ggtitle(expression(paste('A. Pathway proteins, cross validated by leaving out media types, ', italic(S.~cerevisiae)))) +
  xlab(expression(Model~Performance~(Cross~Validated~R^2))) +
  geom_vline(xintercept = 0, lty = 2) +
  # xlim(-1, 1) +
  ylab('Count') +
  facet_zoom(r_sq > -1 & r_sq < 1, xlim = c(-1, 1), zoom.size = 2);block_cross_validation_subsystem_ridge

block_growth_cross_validation_subsystem_ridge <- subsys_only_r_squared_ridge_group_cv_growth %>% 
  ggplot(aes(x = r_sq)) +
  geom_histogram(fill = yeast_col, bins = 100) +
  theme_bw() +
  ggtitle(expression(paste('A. Pathway proteins, cross validated by leaving out growth rates, ', italic(S.~cerevisiae)))) +
  xlab(expression(Model~Performance~(Cross~Validated~R^2))) +
  geom_vline(xintercept = 0, lty = 2) +
  xlim(-1, 1) +
  ylab('Count')

## Whole proteome blocked cross validation

block_cross_validation_whole_proteome_ridge <- all_data_subsys_r_squared_ridge_group_cv %>% 
  ggplot(aes(x = r_sq)) +
  geom_histogram(fill = yeast_col, bins = 100) +
  theme_bw() +
  ggtitle(expression(paste('B. Proteome, cross validated by leaving out media types, ', italic(S.~cerevisiae)))) +
  xlab(expression(Model~Performance~(Cross~Validated~R^2))) +
  geom_vline(xintercept = 0, lty = 2) +
  xlim(-1, 1) +
  ylab('Count')

block_growth_cross_validation_whole_proteome_ridge <- all_data_subsys_r_squared_ridge_group_cv_growth %>% 
  ggplot(aes(x = r_sq)) +
  geom_histogram(fill = yeast_col, bins = 100) +
  theme_bw() +
  ggtitle(expression(paste('B. Proteome, cross validated by leaving out growth rates, ', italic(S.~cerevisiae)))) +
  xlab(expression(Model~Performance~(Cross~Validated~R^2))) +
  geom_vline(xintercept = 0, lty = 2) +
  xlim(-1, 1) +
  ylab('Count')

# Plotting leave 5 out cross validation for pathways and whole pro --------

leave_5_out_cross_validation_subsystem_ridge <- subsys_only_r_squared_ridge_leave_5 %>% 
  ggplot(aes(x = r_sq)) +
  geom_histogram(fill = yeast_col, bins = 100) +
  theme_bw() +
  ggtitle(expression(paste('A. Pathway proteins, leave-5-out cross validation, ', italic(S.~cerevisiae)))) +
  xlab(expression(Model~Performance~(Cross~Validated~R^2))) +
  geom_vline(xintercept = 0, lty = 2) +
  xlim(-1, 1) +
  ylab('Count')

leave_5_out_cross_validation_whole_proteome_ridge <- all_data_subsys_r_squared_ridge_leave_5 %>% 
  ggplot(aes(x = r_sq)) +
  geom_histogram(fill = yeast_col, bins = 100) +
  theme_bw() +
  ggtitle(expression(paste('B. Proteome, leave-5-out cross validation, ', italic(S.~cerevisiae)))) +
  xlab(expression(Model~Performance~(Cross~Validated~R^2))) +
  geom_vline(xintercept = 0, lty = 2) +
  xlim(-1, 1) +
  ylab('Count')

# Arranging all figures ---------------------------------------------------

## Putting together main text figure with R^2 distributions from yeast
yeast_only_r_squared <- ggarrange(yeast_linear_reg_rsq + 
                                    ggtitle('B. Single proteins') + 
                                    xlim(-0.5, 1) + 
                                    theme(axis.title = element_text(size = 13),
                                          axis.text = element_text(size = 13),
                                          plot.margin = unit(c(0,0.2,0,1), 'lines')) +
                                    scale_y_continuous(breaks= pretty_breaks()), 
          subsys_only_r_sq_hist_ridge + 
            ggtitle('C. Subsystem proteins only') + 
            xlim(-0.5, 1) + 
            theme(axis.title = element_text(size = 13),
                  axis.text = element_text(size = 13), 
                  plot.margin = unit(c(0,0.2,0,1), 'lines')),
          all_proteome_data_r_squared_hist_ridge + 
            ggtitle('D. Proteome') + 
            xlim(-0.5, 1) + 
            xlab(expression(Model~Performance~(Cross~Validated~R^2))) + 
            theme(axis.title = element_text(size = 13),
                  axis.text = element_text(size = 13)), 
          nrow = 3, align = 'hv')

ggsave(yeast_only_r_squared, 
       filename = 'figures/r_sq_combined_yeast_only.pdf', height = 8, width = 7)

### Putting togetherfigure with R^2 distributions

r_sq_from_e_coli_b_sub <- ggarrange(ecoli_linear_reg_rsq, bsub_linear_reg_rsq, nrow = 2, align = 'hv')

ggsave(r_sq_from_e_coli_b_sub,
       filename = paste0('figures/figSX_', current_date, '_r_sq_distributions_e_coli_b_sub.pdf'), 
       height = 6.8, width = 8.5)

### Putting together supp text figure with R^2 distributions but only the filtered subset of proteins that are found in the single protein prediction susbet

r_sq_combined_distributions_single_subset <- ggarrange(subsys_only_r_sq_hist_ridge_single_subset, 
          all_proteome_data_r_squared_hist_ridge_single_subset,
          nrow = 2, align = 'hv')

ggsave(r_sq_combined_distributions_single_subset, 
       filename = paste0('figures/figSX_', current_date, '_r_sq_combined_distributions_single_subset.pdf'), 
       height = 6.8, width = 8.5)

### Putting together LASSO and Ridge Regression for pathways and full proteomes

lasso_versus_ridge_rsq <- ggarrange(whole_prot_model_comparison +
                                      scale_y_continuous(limits = c(-0.75, 1)) +
                                      scale_x_continuous(limits = c(-0.75, 1)) +
                                      coord_fixed(ratio = 1), 
                                    subsys_prot_model_comparison +
                                      scale_y_continuous(limits = c(-0.75, 1)) +
                                      scale_x_continuous(limits = c(-0.75, 1)) +
                                      coord_fixed(ratio = 1), 
                                    align = 'hv')

ggsave(lasso_versus_ridge_rsq,
       filename = paste0('figures/figsX_', current_date, '_lasso_versus_ridge_rsq.pdf'),
       width = 13, height = 9)

### Putting together a supp text figure of the block cross validation distribution of R2 for the subsystem to full proteome predictions

block_substrate_cv_comparison <- ggarrange(block_cross_validation_subsystem_ridge, 
                                           block_cross_validation_whole_proteome_ridge, nrow = 2)

ggsave(block_substrate_cv_comparison,
       filename = paste0('figures/figsX_', 
                         current_date, 
                         '_r_squared_difs_cross_validation_substrate.pdf'),
       height = 7, 
       width = 8)

### Putting together a supp text figure of the block cross validation distribution of R2 for the subsystem to full proteome predictions. This time the block cross validation is done across different growth rates.

block_growth_cv_comparison <- ggarrange(block_growth_cross_validation_subsystem_ridge,
                                 block_growth_cross_validation_whole_proteome_ridge, nrow = 2)

ggsave(block_growth_cv_comparison,
       filename = paste0('figures/figsX_', current_date, '_r_squared_difs_cross_validation_growth.pdf'),
       height = 7, 
       width = 8)

## Putting together a supp text figure of the leave 5 out cross validation of the distribution of R2 for the subsystem to full proteome predictions.

leave_5_out_cv_comparison <- ggarrange(leave_5_out_cross_validation_subsystem_ridge, 
                                       leave_5_out_cross_validation_whole_proteome_ridge, 
                                       nrow = 2)

ggsave(leave_5_out_cv_comparison,
       filename = paste0('figures/figsX_', current_date, '_r_squared_leave_5_out_ridge.pdf'),
       height = 7, 
       width = 8)

