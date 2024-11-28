## plotting single examples of enzyme-rate relationships
library(scales)

gene_of_interest1 <- 'YNR050C'
gene_of_interest2 <- 'YML004C'

yeast_col <- 'darkseagreen4'
bsub_col <- 'lightgoldenrod1'
ecoli_col <- 'steelblue4'

## get Spearman correlation coefficients
gene_of_interest1_cor <- round(ag_metrics_across[ag_metrics_across$gene_name == gene_of_interest1, ]$spearman_cor, 2)
gene_of_interest2_cor <- round(ag_metrics_across[ag_metrics_across$gene_name == gene_of_interest2, ]$spearman_cor, 2)

## get raw values used to calculate pearson correlation coefficients & fit Bayesian model.
gene_number1 <- joined_df_flux_prot_mean_no_negative[joined_df_flux_prot_mean_no_negative$Gene == gene_of_interest1, ]$gene_as_factor_number[1]
gene_number2 <- joined_df_flux_prot_mean_no_negative[joined_df_flux_prot_mean_no_negative$Gene == gene_of_interest2, ]$gene_as_factor_number[1]

## sample from full posterior -- just to show model fit.
set.seed(10001)
sample_out_from_posterior <- sample(c(1:nrow(fit_hackett_out$beta_1)), size = 1000, replace = FALSE)

# ## Plotting the data with the hierarchical model fit. -------------------

p_positive_out <- joined_df_flux_prot_mean_no_negative %>% 
  filter(Gene == gene_of_interest1) %>%
  ggplot(aes(x = z_score_prot, y = y_rep)) +
  ggtitle('A. Saccharopine dehydrogenase (LYS9)') +
  geom_abline(intercept = median(fit_hackett_out$beta_0[,gene_number1]),
              slope = median(fit_hackett_out$beta_1[,gene_number1]),
              colour = yeast_col, lwd = 2) +
  geom_abline(data = data.frame(slope = fit_hackett_out$beta_1[sample_out_from_posterior,gene_number1],
                                intercept = fit_hackett_out$beta_0[sample_out_from_posterior,gene_number1]),
              aes(intercept = intercept, slope = slope),
              alpha = 0.05, colour = 'grey') +
  geom_point(data = joined_df_flux_prot_mean_no_negative_w_pp %>% 
               filter(Gene == gene_of_interest1), 
             aes(x = z_score_prot, y = z_score_flux), 
             fill = yeast_col, pch = 21, colour = 'black', size = 4) +
  theme_bw() +
  theme(title = element_text(size = 15)) +
  xlab('log2(Protein abundance) (Z-score transformed)') +
  ylab('Reaction rate (Z-score transformed)')

p_middle_out <- joined_df_flux_prot_mean_no_negative_w_pp %>% 
  filter(Gene == gene_of_interest2) %>%
  ggplot(aes(x = z_score_prot, y = y_rep)) +
  ggtitle('B. Monomeric glyoxalase I (GLO1)') +
  geom_abline(intercept = median(fit_hackett_out$beta_0[,gene_number2]),
              slope = median(fit_hackett_out$beta_1[,gene_number2]),
              colour = yeast_col, lwd = 2) +
  geom_abline(data = data.frame(slope = fit_hackett_out$beta_1[sample_out_from_posterior,gene_number2],
                                intercept = fit_hackett_out$beta_0[sample_out_from_posterior,gene_number2]),
              aes(intercept = intercept, slope = slope),
              alpha = 0.05, colour = 'grey') +
  geom_point(data = joined_df_flux_prot_mean_no_negative_w_pp %>% 
               filter(Gene == gene_of_interest2), 
             aes(x = z_score_prot, y = z_score_flux), 
             fill = yeast_col, pch = 21, colour = 'black', size = 4) +
  theme_bw() +
  theme(title = element_text(size = 15)) +
  xlab('log2(Protein abundance) (Z-score transformed)') +
  ylab('Reaction rate (Z-score transformed)')

singlearranged_out <- ggarrange(p_positive_out,
                                p_middle_out, 
                                nrow = 1, 
                                ncol = 2)

ggsave(singlearranged_out, 
       filename = 'figures/single-enzyme-rate-relationships-bayesian.pdf', 
       height = 4.7, width = 11)


# Plotting the raw data ---------------------------------------------------

p_positive_out_lls <- joined_df_flux_prot_mean_no_negative %>% 
  filter(Gene == gene_of_interest1) %>%
  ggplot(aes(x = z_score_prot, y = z_score_unlog_flux)) +
  ggtitle('A. Saccharopine dehydrogenase') +
  labs(subtitle = '(LYS9)') +
  geom_point(fill = yeast_col, pch = 21, colour = 'black', size = 4) +
  theme_bw() +
  # geom_smooth(method = 'lm', colour = 'darkblue') +
  theme(title = element_text(size = 13),
        strip.text.x = element_text(size = 13),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 11)) +
  xlim(-2.5, 2.5) +
  ylim(-2.5, 2.1) +
  annotate(geom = 'text', x = -1.5, 2, label = bquote(rho == .(gene_of_interest1_cor))) +
  xlab('log(Protein abundance)\n(Z-score transformed)') +
  ylab('Reaction rate\n(Z-score transformed)');p_positive_out_lls

p_middle_out_lls <- joined_df_flux_prot_mean_no_negative %>% 
  filter(Gene == gene_of_interest2) %>%
  ggplot(aes(x = z_score_prot, y = z_score_unlog_flux)) +
  ggtitle('B. Monomeric glyoxalase I') +
  labs(subtitle = '(GLO1)') +
  # ggtitle('c) Monomeric glyoxalase I\n   (Glycolysis/Glutathione Biosynthesis)') +
  geom_point(fill = yeast_col, pch = 21, colour = 'black', size = 4) +
  theme_bw() +
  xlim(-2.5, 2.5) +
  ylim(-2.5, 2.1) +
  annotate(geom = 'text', x = -1.5, 2, label = bquote(rho == .(gene_of_interest2_cor))) +
  # geom_smooth(method = 'lm', colour = 'darkblue') +
  theme(title = element_text(size = 13),
        strip.text.x = element_text(size = 13),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 11)) +
  xlab('log(Protein abundance)\n(Z-score transformed)') +
  ylab('Reaction rate\n(Z-score transformed)');p_middle_out_lls

singlearranged_out_lls <- ggarrange(p_positive_out_lls + coord_fixed(ratio = 1), 
                                    p_middle_out_lls + coord_fixed(ratio = 1), nrow = 2, ncol = 1)

# ggsave(singlearranged_out_lls, 
#        filename = 'figures/single-enzyme-rate-relationships-lls-no-line.pdf', 
#        height = 4.7, width = 9.8)

# plotting different single protein to rate metrics -----------------------

pearson_fig <- ag_metrics_across %>% 
  ggplot(aes(x = pearson_cor, fill = dataset)) +
  geom_histogram() +
  theme_bw() +
  facet_wrap(~dataset, nrow = 3) +
  xlab('Pearson Correlation Coefficient') + 
  theme(strip.text = element_text(face = "italic"),
        legend.position = 'none',
        strip.text.x = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 11)) +
  scale_fill_manual(values = c(bsub_col, ecoli_col, yeast_col))

pearson_log_abundance_fig <- ag_metrics_across %>% 
  ggplot(aes(x = pearson_cor_prot_log, fill = dataset)) +
  geom_histogram() +
  theme_bw() +
  facet_wrap(~dataset, nrow = 3) +
  xlab('Pearson Correlation Coefficient\n(log(Abundance))') + 
  theme(strip.text = element_text(face = "italic"),
        legend.position = 'none',
        strip.text.x = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 11)) +
  scale_fill_manual(values = c(bsub_col, ecoli_col, yeast_col))

mut_info_fig <- ag_metrics_across %>% 
  ggplot(aes(x = mut_info_bits, fill = dataset)) +
  geom_histogram() +
  theme_bw() +
  facet_wrap(~dataset, nrow = 3) +
  xlab('Mutual Information (bits)') + 
  theme(strip.text = element_text(face = "italic"),
        legend.position = 'none',
        strip.text.x = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 11)) +
  scale_fill_manual(values = c(bsub_col, ecoli_col, yeast_col))

### need to run plotting_hierarchical_model_partial_pooling.R to get these objects
posterior_pred_all_plot <- posterior_pred_all %>% 
  ggplot(aes(x = posterior_draws_out)) +
  geom_histogram(aes(fill = organism)) +
  theme_bw() +
  xlim(-1.1, 1.1) +
  facet_wrap(~organism, nrow = 3) +
  labs(y = expression("Posterior Probability Sample Count (x" ~ 10^3~')')) +
  geom_vline(xintercept = 0) +
  theme(strip.text = element_text(face = "italic"),
        legend.position = 'none',
        # title = element_text(size = 13),
        strip.text.x = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 11)) +
  scale_fill_manual(values = c(bsub_col, ecoli_col, yeast_col)) +
  ## This is here to convert the y axis into units of **thousands** of samples
  scale_y_continuous(labels = function(x) x/1000) +
  ggtitle('D. Pooled correlation coefficients') +
  xlab('Pooled Correlated Coefficient\nbetween Abundance and Rate');posterior_pred_all_plot

## want to annotate the two specific examples shown in the first two panels, need to make an annotation dataframe.
ann_text <- data.frame(dataset_f = rep('S. cerevisiae', 2), 
                       dataset = rep('S. cerevisiae', 2), 
                       spearman_cor = c(gene_of_interest1_cor, gene_of_interest2_cor), 
                       lab = "text", gene_name = c('LYS9', 'GLO1'))

## change the ordering of the levels
ag_metrics_across$dataset_f <- factor(ag_metrics_across$dataset, levels = c('S. cerevisiae', 'E. coli', 'B. subtilis'))
ann_text$dataset_f <- factor(ann_text$dataset, levels = c('S. cerevisiae', 'E. coli', 'B. subtilis'))

## calculate the median correlation coefficient for each dataset
ag_metrics_across_mean_spearman <- ag_metrics_across %>% 
  group_by(dataset_f) %>% 
  summarize(mean_spearman_cor = mean(spearman_cor))

## add in the numbers of observations for each
ag_metrics_across_mean_numbers_each_dataset <- ag_metrics_across %>% 
  group_by(dataset_f, dataset) %>% 
  summarize(number_pairs = n()) %>% 
  dplyr::mutate(number_pairs_text = paste0('n = ', number_pairs))

spearman_fig_appended <- ag_metrics_across %>% 
  dplyr::filter(number_vals > 2) %>% ## filter out the one rate/transcript pair that has 2 observations,
  ## the pearson correlation is 1
  ggplot(aes(x = spearman_cor, fill = dataset, colour = dataset)) +
  geom_histogram() +
  theme_bw() +
  facet_wrap(~dataset_f, nrow = 3) +
  geom_text(data = ann_text, aes(label = gene_name, x = spearman_cor), y = 7) +
  geom_text(data = ann_text, aes(label = gene_name, x = spearman_cor), y = 7) +
  geom_label(data = ag_metrics_across_mean_numbers_each_dataset, 
             x = 0.9,
             y = 11, aes(label = number_pairs_text),
             colour = 'black', fill = 'white') + 
  theme(strip.text = element_text(face = "italic"),
        legend.position = 'none',
        title = element_text(size = 13),
        strip.text.x = element_text(size = 13),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 11)) +
  geom_segment(data = ann_text,
               aes(x = spearman_cor, 
                   y = 6, xend = spearman_cor, yend = 3),
               arrow = arrow(length = unit(0.3, "cm"), 
                             type = 'closed')) +
  geom_vline(data = ag_metrics_across_mean_spearman, 
             aes(xintercept = mean_spearman_cor), colour = 'grey40', lty = 2) +
  scale_y_continuous(breaks= pretty_breaks()) +
  xlim(-1.05, 1.05) +
  geom_vline(xintercept = 0) +
  scale_fill_manual(values = c(bsub_col, ecoli_col, yeast_col)) +
  scale_colour_manual(values = c(bsub_col, ecoli_col, yeast_col)) +
  ylab('Protein-Reaction Relationship Count') +
  ggtitle('C.') +
  xlab('Spearman Correlation Coefficient between\nlog(Abundance) and Reaction Rate');spearman_fig_appended

summary_metrics_pearson_bayes <- ggarrange(spearman_fig_appended, 
                                           posterior_pred_all_plot, nrow = 1, align = 'hv')

summary_metrics_pearson_bayes2 <- ggarrange(singlearranged_out_lls, 
          summary_metrics_pearson_bayes, 
          nrow = 1, align = 'hv', widths = c(0.39, 1))

summary_metrics_pearson_bayes3 <- ggarrange(singlearranged_out_lls, 
                                            spearman_fig_appended, 
                                            nrow = 1, align = 'hv', widths = c(0.39, 0.5))

all_single_fig_metrics <- ggarrange(pearson_fig + ggtitle('A.') + scale_y_continuous(breaks= pretty_breaks()), 
                                    pearson_log_abundance_fig + ggtitle('B.') + scale_y_continuous(breaks= pretty_breaks()),
                                    mut_info_fig + ggtitle('C.') + scale_y_continuous(breaks= pretty_breaks()),
                                    posterior_pred_all_plot + ggtitle('D.'),
                                    align = 'hv', 
                                    nrow = 1)
## supplementary figure
ggsave(all_single_fig_metrics, 
       filename = paste0('figures/figsX_', current_date, '_all_single_fig_metrics.pdf'),
       width = 15, height = 11.5*0.6)

current_date <- format(Sys.Date(), "%Y%m%d")

ggsave(summary_metrics_pearson_bayes3, 
       filename = paste0('figures/fig2_', current_date, '_summary_metrics_pearson.pdf'),
       width = 14.5*1.05*0.7, height = 9.6*0.8*1.03)
