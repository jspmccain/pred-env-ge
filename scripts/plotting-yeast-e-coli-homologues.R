### make a figure of the coefficient estimates across yeast and e coli

library(ggplot2)
library(ggrepel)

### first need to run data processing scripts, and also ecoli_fasta_formatting.R, and the hierarchical model scripts

### Key objects
# yeast_gene_to_ecoli_gene_connection: the maps the names of the yeast genes to the e coli ones
# beta_1_estimates_hackett_df: coefficient estimates with the yeast names for the one level Bayesian hierarchical model
# beta_1_estimates_davidi_df: coefficient estimates with the yeast names for the one level Bayesian hierarchical model

### First change the names of the mapping dataframe, so that the coefficient estimate dataframe can be smoothly mapped.
yeast_gene_to_ecoli_gene_connection_renamed <- yeast_gene_to_ecoli_gene_connection %>% 
  dplyr::rename(Gene = sseqid)

coef_combined_df_yeast_e_coli <- beta_1_estimates_hackett_df %>% 
  dplyr::rename(beta_1_yeast = beta_1,
                lower_yeast = .lower,
                upper_yeast = .upper) %>% 
## join the names of the e coli and yeast homologues
  inner_join(yeast_gene_to_ecoli_gene_connection_renamed, by = "Gene") %>% 
## join the beta estimates, but specifying the names of the estimates to have e coli in the name
  inner_join(beta_1_estimates_davidi_df %>% 
               dplyr::rename(beta_1_ecoli = beta_1,
                             lower_ecoli = .lower,
                             upper_ecoli = .upper), 
             by = "gene_name")
  
homologous_protein_comparison_coefficient <- coef_combined_df_yeast_e_coli %>% 
  ggplot(aes(x = beta_1_yeast, y = beta_1_ecoli)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower_ecoli, ymax = upper_ecoli), alpha = 0.3) +
  geom_errorbarh(aes(xmin = lower_yeast, xmax = upper_yeast), alpha = 0.3) +
  theme_bw() +
  xlim(-1, 1) + 
  ylim(-1, 1) +
  xlab(expression(italic(S.~cerevisiae)~'Pooled Correlation Coefficient between Normalised Enzyme Abundance and Rate')) +
  ylab(expression(italic(E.~coli)~'Pooled Correlation Coefficient between Normalised Enzyme Abundance and Rate')) +
  coord_fixed(ratio = 1) +
  theme(axis.text = element_text(size = 13),
        axis.title = element_text(size = 12));homologous_protein_comparison_coefficient

ggsave(homologous_protein_comparison_coefficient,
       filename = paste0('figures/figsX_', current_date, '_homologous_protein_comparison_coefficient.pdf'), 
       height = 10, width = 10)
