# single protein prediction using cross validation similar to whole proteome prediction

# source('scripts/hackett-full-proteome-prediction.R')
source('scripts/functions_for_lm_cv.R')
source('scripts/functions_for_lasso_cv.R')

set.seed(12345)

# # Cross validation and one-protein linear Regression model for yeast dataset --------

# unique reaction rate id's that have single protein pairs
single_protein_pairs_rxns <- unique(joined_df_flux_prot_mean_no_negative$reaction_id)

subset_yeast_fluxes_singles <- loop_through_and_cv_fluxes_single(gene_flux_data_frame = hackett_fluxes_subset, 
                                                                 wide_proteome_data_frame = hackett_prot_transformed, 
                                                                 reaction_name_vector = single_protein_pairs_rxns, 
                                                                 number_cv_samples = 100, 
                                                                 leave_n_out = 2)

singles_only_r_squared_yeast_rmsd <- convert_model_summaries(subset_yeast_fluxes_singles)

# assessing single protein prediction for GSH2, but z-score transforming this protein abundance so it is directly comparable with the z-score transformed
# model.

## z score the protein abundance
hackett_prot_transformed_z_score_gsh2 <- scale(hackett_prot_transformed$YOL049W) %>% as.numeric()
## make a new dataframe to rewrite the abundance
hackett_prot_transformed_single_gsh2 <- hackett_prot_transformed
## changed abundance to the z-score normalized
hackett_prot_transformed_single_gsh2$YOL049W <- hackett_prot_transformed_z_score_gsh2

gsh2_z_score_transformed_model <- loop_through_and_cv_fluxes_single(gene_flux_data_frame = hackett_fluxes_subset, 
                                                                 wide_proteome_data_frame = hackett_prot_transformed_single_gsh2, 
                                                                 reaction_name_vector = "r_0485", 
                                                                 number_cv_samples = 100, 
                                                                 leave_n_out = 2)

singles_only_r_squared_yeast_rmsd_gsh2 <- convert_model_summaries(gsh2_z_score_transformed_model)

# # Cross validation and one-protein linear Regression model for E coli dataset --------

single_protein_pairs_rxns_e_coli <- genes_fluxes_long_df_z_score$reaction_name %>% unique()

genes_fluxes_long_df_z_score_cc_subset <- genes_fluxes_long_df_z_score %>% 
  dplyr::select(reaction_name, flux, protein_amount_mean, exp_condition) %>% 
  na.omit()

subset_ecoli_fluxes_singles <- loop_through_and_cv_fluxes_single(gene_flux_data_frame = genes_fluxes_long_df_z_score_cc_subset, 
                                                                 wide_proteome_data_frame = data.frame(nothing = c()), 
                                                                 reaction_name_vector = single_protein_pairs_rxns_e_coli, 
                                                                 number_cv_samples = 100, org_choice = "ecoli",
                                                                 leave_n_out = 2)

singles_only_r_squared_e_coli_rmsd <- convert_model_summaries(subset_ecoli_fluxes_singles)

# ## wondering how many observations per regression
# genes_fluxes_long_df_z_score_cc_subset %>% 
#   group_by(reaction_name) %>% 
#   summarize(number_vals = n()) %>% 
#   ggplot(aes(x = number_vals)) +
#   geom_histogram()

# cross validation and one protein linear regression model for B. subtilis --------

single_protein_pairs_rxns_b_sub <- chu_transformed_cc$flux_name %>% unique()

chu_number_obs <- chu_transformed_cc %>% 
  group_by(flux_name) %>% 
  summarize(number_obs_per_flux = n())

#subsetting only those fluxes with at least five observations (2 fold cross validation, need at least three points to fit a Linear regression)
chu_number_obs_rxns <- chu_number_obs[chu_number_obs$number_obs_per_flux >= 5, ]$flux_name

## subsetting the input dataframe for only these reactions
chu_transformed_cc_b_sub_subset <- chu_transformed_cc %>% 
  dplyr::filter(flux_name %in% chu_number_obs_rxns)

subset_bsub_fluxes_singles <- loop_through_and_cv_fluxes_single(gene_flux_data_frame = chu_transformed_cc_b_sub_subset, 
                                                                 wide_proteome_data_frame = data.frame(nothing = c()), 
                                                                 reaction_name_vector = chu_number_obs_rxns, 
                                                                 number_cv_samples = 100, 
                                                                 org_choice = "bsubtilis",
                                                                 leave_n_out = 1)

singles_only_r_squared_b_sub_rmsd <- convert_model_summaries(subset_bsub_fluxes_singles)

