## This script examines proteomic profiles from Hackett et al, and the correlation between profiles that differ in their growth rates. The goal for 
## this analysis is to see what an estimated amount of measurement error is in the data. Note that some of the differences across treatments is
## because of biological variation and not measurement error.

single_proteins_observed <- unique(joined_df_flux_prot_mean_no_negative$Gene)

pairwise_growth_correlations <- function(condition_1, condition_2, 
                                         condition_group = 'P', 
                                         hackett_prot_df = hackett_prot,
                                         return_cor_val = TRUE){
  # condition_1 <- '05'
  # condition_2 <- '11'
  # condition_group <- 'P'
  # hackett_prot_df <- hackett_prot
  
  ### Get the subset of the proteome with only the single reaction rate mapped proteins
  hackett_prot_df_singles <- hackett_prot_df %>% 
    dplyr::filter(Gene %in% single_proteins_observed)
  
  ### Subset only those columns with a specific condition type
  hackett_prot_condition <- hackett_prot_df_singles[grepl(x = names(hackett_prot_df_singles),
                                               pattern = condition_group, 
                                               ignore.case = TRUE)]
  
  ### subset only those columns with a specific growth rate
  condition_1_col <- hackett_prot_condition[grepl(x = names(hackett_prot_condition), 
                                        pattern = condition_1)]
  condition_2_col <- hackett_prot_condition[grepl(x = names(hackett_prot_condition), 
                                        pattern = condition_2)]
  ## get the correlation coefficient
  cor_out <- cor(condition_1_col, condition_2_col, method = 'pearson')[1,1]
  
  ## make a conditional statement about what value to return
  if(return_cor_val){
    return(cor_out)
  }
  if(!return_cor_val){
    return(data.frame(condition_1 = condition_1_col, 
                      condition_2 = condition_2_col))
  }
}

all_pairwise_growth <- function(condition_group, condition_list = c('05', '11', '16', '22', '30')){
  cor_list_out <- c()
  for(i in 2:length(condition_list)){
    cor_i <- pairwise_growth_correlations(condition_1 = condition_list[i - 1], 
                                 condition_2 = condition_list[i], 
                                 condition_group = condition_group)
    cor_list_out <- c(cor_i, cor_list_out)
  }
  return(cor_list_out)
}

all_p_cors <- all_pairwise_growth(condition_group = 'P')
all_u_cors <- all_pairwise_growth(condition_group = 'U')
all_c_cors <- all_pairwise_growth(condition_group = 'C')
all_n_cors <- all_pairwise_growth(condition_group = 'N')
all_l_cors <- all_pairwise_growth(condition_group = 'L')

df_cor_out_conditions <- data.frame(cor_val = c(all_p_cors, all_u_cors, all_c_cors, all_n_cors, all_l_cors),
           condition_type = c(rep(c('Phosphate', 'Uracil', 'Carbon', 'Nitrate', 'Leucine'), each = 4)))

p1_out <- df_cor_out_conditions %>% 
  ggplot(aes(x = cor_val)) +
  geom_histogram() +
  facet_wrap(~condition_type) +
  theme_bw() + 
  xlab('Pearson Correlation') +
  theme(text = element_text(size = 16))

# plotting growth condition specific --------------------------------------

pair_of_p_out <- pairwise_growth_correlations(condition_1 = '11', 
                                           condition_2 = '16', 
                                           condition_group = 'L', 
                                           hackett_prot_df = hackett_prot,
                                           return_cor_val = FALSE)

pearson_cor_val <- round(cor(pair_of_p_out$L0.11, pair_of_p_out$L0.16, method = "pearson"), 3)
spearman_cor_val <- round(cor(pair_of_p_out$L0.11, pair_of_p_out$L0.16, method = "spearman"), 3)


comparison_of_singles_prots_leucine_016_011 <- pair_of_p_out %>% 
  ggplot(aes(x = L0.16, y = L0.11)) +
  geom_point(size = 4, colour = 'grey30') +
  geom_abline(intercept = 0, slope = 1) +
  theme_bw() +
  xlim(-2.75, 2.5) +
  ylim(-2.75, 2.5) +
  annotate(geom = 'text', x = -1, y = 2, label = paste('Pearson correlation coefficient = ', pearson_cor_val),
           size = 4) +
  xlab('Relative protein abundance, Leucine Limited Growth\n(Growth rate = 0.16 per hour)') +
  ylab('Relative protein abundance, Leucine Limited Growth\n(Growth rate = 0.11 per hour)') +
  coord_fixed(ratio = 1) +
  theme(text = element_text(size = 13));comparison_of_singles_prots_leucine_016_011

ggsave(filename = paste0('figures/figsX_', current_date, '_comparison_of_singles_prots_leucine_016_011.pdf'), 
       plot = comparison_of_singles_prots_leucine_016_011, 
       height = 9.5*0.8, width = 12*0.8)

