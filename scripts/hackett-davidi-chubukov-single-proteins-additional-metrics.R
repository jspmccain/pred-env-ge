### getting metrics for comparing single protein to rate relationships

source('scripts/mutual_info_b_spline.R')

## dataframe from Davidi
davidi_aggregated_metrics <- genes_fluxes_long_df_z_score_cc %>% 
  group_by(gene_name) %>% 
  summarize(pearson_cor = cor(flux, protein_amount_mean, method = 'pearson'),
            pearson_cor_prot_log = cor(flux, log2(protein_amount_mean), method = 'pearson'),
            pearson_z_score = cor(z_score_flux, z_score_prot),
            spearman_cor = cor(flux, protein_amount_mean, method = 'spearman'),
            number_vals = n(),
            mut_info = mutual.information2(flux, protein_amount_mean, bins = number_vals^(1/3)),
            mut_info_bits = natstobits(mut_info))

## dataframe from Hackett (note that protein amount is log2 transformed)
hackett_aggregated_metrics <- joined_df_flux_prot_mean_no_negative %>%
  ## need to filter out this protein because it doesn't carry any flux through these conditions
  dplyr::filter(Gene != 'YHR128W') %>% 
  dplyr::mutate(protein_amount_unlogged = 2^protein_amount) %>% 
  group_by(Gene) %>% 
  summarize(pearson_cor = cor(flux, protein_amount_unlogged, method = 'pearson'),
            pearson_cor_prot_log = cor(flux, log2(protein_amount_unlogged), method = 'pearson'),
            pearson_z_score = cor(z_score_flux, z_score_prot),
            spearman_cor = cor(flux, protein_amount_unlogged, method = 'spearman'),
            number_vals = n(),
            mut_info = mutual.information2(flux, protein_amount_unlogged, bins = number_vals^(1/3)),
            mut_info_bits = natstobits(mut_info))

## dataframe from chubukov
chubukov_aggregated_metrics <- chu_transformed_cc %>% 
  group_by(enzyme_bsu_numbers) %>% 
  summarize(pearson_cor = cor(flux, mean_expression_val, method = 'pearson'),
            pearson_cor_prot_log = cor(flux, log2(mean_expression_val), method = 'pearson'),
            pearson_z_score = cor(z_score_flux, z_score_exp),
            spearman_cor = cor(flux, mean_expression_val, method = 'spearman'),
            number_vals = n(),
            mut_info = mutual.information2(flux, mean_expression_val, bins = number_vals^(1/3)),
            mut_info_bits = natstobits(mut_info))

## make these have common names
davidi_aggregated_metrics_2 <- davidi_aggregated_metrics %>% 
  mutate(dataset = rep('E. coli'))

hackett_aggregated_metrics_2 <- hackett_aggregated_metrics %>% 
  mutate(dataset = rep('S. cerevisiae')) %>% 
  dplyr::rename(gene_name = Gene)

chubukov_aggregated_metrics_2 <- chubukov_aggregated_metrics %>% 
  mutate(dataset = rep('B. subtilis')) %>% 
  dplyr::rename(gene_name = enzyme_bsu_numbers)

## aggregate into a single dataframe
ag_metrics_across <- rbind(davidi_aggregated_metrics_2, hackett_aggregated_metrics_2, chubukov_aggregated_metrics_2)

