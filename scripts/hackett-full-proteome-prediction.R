
library(reshape2)

# # Hackett et al data proteome prediction of rates

source('scripts/functions_for_lasso_cv.R')

# # Additional data processing for the hackett input data: --------------

hackett_fluxes_subset <- hackett_fluxes %>% 
  ## Remove the flux variability columns
  dplyr::select(Reaction, Model_Reaction_ID, 
                contains('QP'),
                -contains('FVA')) %>% 
  ## make it into a long dataframe format
  melt(variable.name = 'exp_condition_qp', 
       value.name = 'flux') %>% 
  dplyr::rename(reaction_name = Model_Reaction_ID)

## spot checking that this worked fine:
# hackett_fluxes_subset %>% 
#   dplyr::filter(exp_condition_qp == 'P0.05_QP', reaction_name == 'r_0007')

## Remove the '_QP' from the exp_condition string identifier.
hackett_fluxes_subset$exp_condition <- gsub(pattern = "_QP", 
                                            replacement = "", 
                                            x = hackett_fluxes_subset$exp_condition_qp)

## hackett_prot has Gene name and condition name as columns. For the following modelling work,
## we want there to be proteins for column names, and conditions for rows. Note that one of the conditions
## has a lower case p, so we need to transform that.
hackett_prot_transformed <- hackett_prot %>% 
  melt(variable.name = 'exp_condition', 
       value.name = 'protein_level') %>% 
  dcast(exp_condition ~ Gene) %>% 
  mutate(exp_condition = fct_recode(exp_condition, "P0.22" = "p0.22"))

# tester <- subset_yeast_fluxes_singles %>% 
  # dplyr::filter(reaction_rate_id == "r_0989")

# # Cross Validation & LASSO for yeast dataset using whole proteome ------------------------------

unique_hackett_fluxes <- hackett_fluxes_subset$reaction_name %>% unique()

if(rerun_all_pred_models){

## Whole proteome using LASSO regression.
all_yeast_fluxes_predicted_list <- loop_through_and_cv_fluxes(gene_flux_data_frame = hackett_fluxes_subset, 
                                                              wide_proteome_data_frame = hackett_prot_transformed, 
                                                              reaction_name_vector = unique_hackett_fluxes, 
                                                              leave_n_out = 2,
                                                              number_cv_samples = 100, 
                                                              org = 'yeast')

## Whole proteome using ridge regression.
all_yeast_fluxes_predicted_list_ridge <- loop_through_and_cv_fluxes(gene_flux_data_frame = hackett_fluxes_subset, 
                                                              wide_proteome_data_frame = hackett_prot_transformed, 
                                                              reaction_name_vector = unique_hackett_fluxes, 
                                                              leave_n_out = 2,
                                                              number_cv_samples = 100, 
                                                              org = 'yeast', 
                                                              alpha = 0)

## Whole proteome using ridge regression, with n = 5 leave out. With this level of subsampling, sometimes all 5 of the non zero fluxes are removed, making the model
## unable to fit in this framework. So these three reactions are left out with this cross validation.
all_yeast_fluxes_predicted_list_ridge_leave_5_out <- loop_through_and_cv_fluxes(gene_flux_data_frame = hackett_fluxes_subset, 
                                                                    wide_proteome_data_frame = hackett_prot_transformed, 
                                                                    reaction_name_vector = unique_hackett_fluxes[unique_hackett_fluxes %!in% c("r_1074", 
                                                                                                                                               "r_1211",
                                                                                                                                               "r_1272")], 
                                                                    leave_n_out = 5,
                                                                    number_cv_samples = 100, 
                                                                    org = 'yeast', 
                                                                    alpha = 0)

## Whole proteome using ridge regression with block treatment cross validation. Some reactions were removed because these had zero flux in all growth conditions
## except one, making the cross validation approach not valid when using a block design.
all_yeast_fluxes_predicted_list_ridge_group_cv <- loop_through_and_cv_fluxes(gene_flux_data_frame = hackett_fluxes_subset, 
                                                                    wide_proteome_data_frame = hackett_prot_transformed, 
                                                                    reaction_name_vector = unique_hackett_fluxes[unique_hackett_fluxes %!in% c("r_1074", 
                                                                                                                                               "r_1211",
                                                                                                                                               "r_1272")], 
                                                                    leave_n_out = 2,
                                                                    number_cv_samples = 5, 
                                                                    org = 'yeast', 
                                                                    alpha = 0, 
                                                                    lno_cv = FALSE,
                                                                    blocking_structure = 'by_condition')

## Whole proteome using ridge regression with block treatment cross validation. The blocks are specified by growth rate instead of by culture condition.
all_yeast_fluxes_predicted_list_ridge_group_cv_growth <- loop_through_and_cv_fluxes(gene_flux_data_frame = hackett_fluxes_subset, 
                                                                             wide_proteome_data_frame = hackett_prot_transformed, 
                                                                             reaction_name_vector = unique_hackett_fluxes, 
                                                                             leave_n_out = 2,
                                                                             number_cv_samples = 5, 
                                                                             org = 'yeast', 
                                                                             alpha = 0, 
                                                                             lno_cv = FALSE,
                                                                             blocking_structure = 'by_growth_rate')

## Save these model outputs because each takes a while to run.
saveRDS(all_yeast_fluxes_predicted_list, file = 'data/intermediate_data/all_yeast_fluxes_predicted_list.rds')
saveRDS(all_yeast_fluxes_predicted_list_ridge, file = 'data/intermediate_data/all_yeast_fluxes_predicted_list_ridge.rds')
saveRDS(all_yeast_fluxes_predicted_list_ridge_group_cv, file = 'data/intermediate_data/all_yeast_fluxes_predicted_list_ridge_group_cv.rds')
saveRDS(all_yeast_fluxes_predicted_list_ridge_leave_5_out, file = 'data/intermediate_data/all_yeast_fluxes_predicted_list_ridge_leave_5_out.rds')
saveRDS(all_yeast_fluxes_predicted_list_ridge_group_cv_growth, file = 'data/intermediate_data/all_yeast_fluxes_predicted_list_ridge_group_cv_growth.rds')

}

## If necessary, read in the output files:
all_yeast_fluxes_predicted_list <- readRDS(file = 'data/intermediate_data/all_yeast_fluxes_predicted_list.rds')
all_yeast_fluxes_predicted_list_ridge <- readRDS(file = 'data/intermediate_data/all_yeast_fluxes_predicted_list_ridge.rds')
all_yeast_fluxes_predicted_list_ridge_group_cv <- readRDS(file = 'data/intermediate_data/all_yeast_fluxes_predicted_list_ridge_group_cv.rds')
all_yeast_fluxes_predicted_list_ridge_leave_5_out <- readRDS(file = 'data/intermediate_data/all_yeast_fluxes_predicted_list_ridge_leave_5_out.rds')
all_yeast_fluxes_predicted_list_ridge_group_cv_growth <- readRDS(file = 'data/intermediate_data/all_yeast_fluxes_predicted_list_ridge_group_cv_growth.rds')

# # Summarizing/plotting all predictions for yeast -----------------------------------

all_yeast_fluxes_predicted_all_data <- all_yeast_fluxes_predicted_list[[1]]
all_yeast_fluxes_predicted_list_summaries <- all_yeast_fluxes_predicted_list[[2]]
all_yeast_fluxes_predicted_coefs <- all_yeast_fluxes_predicted_list[[3]]

## ridge regression summaries
all_yeast_fluxes_predicted_all_data_ridge <- all_yeast_fluxes_predicted_list_ridge[[1]]
all_yeast_fluxes_predicted_list_summaries_ridge <- all_yeast_fluxes_predicted_list_ridge[[2]]
all_yeast_fluxes_predicted_coefs_ridge <- all_yeast_fluxes_predicted_list_ridge[[3]]

label <- expression(R["CV"]^2)

## subsetting the r squared for the all protein predictions -- LASSO
all_data_subsys_r_squared <- convert_model_summaries(all_yeast_fluxes_predicted_all_data)

## subsetting the r squared for the all protein predictions -- Ridge, leave 2 out CV
all_data_subsys_r_squared_ridge <- convert_model_summaries(all_yeast_fluxes_predicted_all_data_ridge)

## subsetting the r squared for the all protein predictions -- Ridge, block CV
all_data_subsys_r_squared_ridge_group_cv <- convert_model_summaries(all_yeast_fluxes_predicted_list_ridge_group_cv[[1]])

## subsetting the r squared for the all protein predictions -- Ridge, leave 5 out CV
all_data_subsys_r_squared_ridge_leave_5 <- convert_model_summaries(all_yeast_fluxes_predicted_list_ridge_leave_5_out[[1]]) %>% 
  ## this step is because the original models were run inclusive of those three single condition only reactions.
  dplyr::filter(reaction_name %!in% c("r_1074", 
                                      "r_1211",
                                      "r_1272"))

## subsetting the r squared for the all protein predictions -- ridge, block cv using growth rate stratification
all_data_subsys_r_squared_ridge_group_cv_growth <- convert_model_summaries(all_yeast_fluxes_predicted_list_ridge_group_cv_growth[[1]])


# # Subsystems of proteins/reactions to estimate reaction rates--------

# First get a list of all subsystems from the yeast GEM
# yeast_gem is the read-in dataframe from hackett-data-processing.R
all_yeast_subsystems <- unique(yeast_gem$SUBSYSTEM)

### get a list of reaction names that have at least 5 proteins per subsystem
proteins_subsystems_df <- get_proteins_per_subsystem(yeast_gem_df = yeast_gem, 
                                                     vec_of_subsystem = unique(yeast_gem$SUBSYSTEM), 
                                                     observed_proteins_vec = paste(names(hackett_prot_transformed),
                                                                                   collapse = ' '))

## spot checking -- Biotin metabolism has a value of 1 in the above df, but there are 4 proteins in that subsystem.
# get_all_protein_names_subsystem(yeast_gem_df = yeast_gem, subsystem_name = "Biotin metabolism")
at_least_5_proteins_per_subsystem <- proteins_subsystems_df[proteins_subsystems_df$proteins_per_sub > 5, ]$vec_of_subsystem
# at_least_2_proteins_per_subsystem <- proteins_subsystems_df[proteins_subsystems_df$proteins_per_sub > 1, ]$vec_of_subsystem

## Now getting all the estimated fluxes from Hackett, and choosing only those that are in subsystems with >5 proteins
hackett_fluxes_with_subsystems <- hackett_fluxes_subset %>% 
  inner_join(yeast_gem %>% 
               dplyr::select(ID, SUBSYSTEM) %>% 
              dplyr::rename(reaction_name = ID)) %>% 
  dplyr::filter(SUBSYSTEM %in% at_least_5_proteins_per_subsystem)

## getting the vector of fluxes that have more than 5 proteins
unique_hackett_fluxes_subsystem_more_than_5 <- hackett_fluxes_with_subsystems$reaction_name %>% unique()

## remove fluxes that are zero across all nutrient conditions except 1
all_zero_except_1 <- c("r_1074", "r_1211", "r_1272")

## fluxes that are non-zero
unique_hackett_fluxes_subsystem_more_than_5_nz <- unique_hackett_fluxes_subsystem_more_than_5[unique_hackett_fluxes_subsystem_more_than_5 %!in% all_zero_except_1]

if(rerun_all_pred_models){

## same cross validation approach as above
all_yeast_fluxes_predicted_list_subsystem <- loop_through_and_cv_fluxes(gene_flux_data_frame = hackett_fluxes_subset, 
                                                              wide_proteome_data_frame = hackett_prot_transformed, 
                                                              reaction_name_vector = unique_hackett_fluxes_subsystem_more_than_5,
                                                              number_cv_samples = 100, 
                                                              leave_n_out = 2,
                                                              org = 'yeast', 
                                                              subsystem_only = "subsystem")

all_yeast_fluxes_predicted_list_subsystem_ridge <- loop_through_and_cv_fluxes(gene_flux_data_frame = hackett_fluxes_subset, 
                                                                        wide_proteome_data_frame = hackett_prot_transformed, 
                                                                        reaction_name_vector = unique_hackett_fluxes_subsystem_more_than_5,
                                                                        number_cv_samples = 100, 
                                                                        leave_n_out = 2,
                                                                        org = 'yeast', 
                                                                        subsystem_only = "subsystem", 
                                                                        alpha_val = 0)

all_yeast_fluxes_predicted_list_subsystem_ridge_leave_5 <- loop_through_and_cv_fluxes(gene_flux_data_frame = hackett_fluxes_subset, 
                                                                              wide_proteome_data_frame = hackett_prot_transformed, 
                                                                              reaction_name_vector = unique_hackett_fluxes_subsystem_more_than_5_nz,
                                                                              number_cv_samples = 100, 
                                                                              leave_n_out = 5,
                                                                              org = 'yeast', 
                                                                              subsystem_only = "subsystem", 
                                                                              alpha_val = 0)

all_yeast_fluxes_predicted_list_subsystem_ridge_group_cv <- loop_through_and_cv_fluxes(gene_flux_data_frame = hackett_fluxes_subset, 
                                                                              wide_proteome_data_frame = hackett_prot_transformed, 
                                                                              reaction_name_vector = unique_hackett_fluxes_subsystem_more_than_5_nz,
                                                                              number_cv_samples = 5, ### this number is required because it is only possible to do 5-fold cv
                                                                              leave_n_out = 2,
                                                                              org = 'yeast', 
                                                                              subsystem_only = "subsystem", 
                                                                              alpha_val = 0, 
                                                                              lno_cv = FALSE,
                                                                              blocking_structure = 'by_condition')

all_yeast_fluxes_predicted_list_subsystem_ridge_z <- loop_through_and_cv_fluxes(gene_flux_data_frame = hackett_fluxes_subset, 
                                                                              wide_proteome_data_frame = hackett_prot_transformed, 
                                                                              reaction_name_vector = unique_hackett_fluxes_subsystem_more_than_5,
                                                                              number_cv_samples = 100, 
                                                                              leave_n_out = 2,
                                                                              org = 'yeast', 
                                                                              subsystem_only = "subsystem", 
                                                                              alpha_val = 0, 
                                                                              z_score_normalize = TRUE)

all_yeast_fluxes_predicted_list_subsystem_ridge_group_cv_growth <- loop_through_and_cv_fluxes(gene_flux_data_frame = hackett_fluxes_subset, 
                                                                                       wide_proteome_data_frame = hackett_prot_transformed, 
                                                                                       reaction_name_vector = unique_hackett_fluxes_subsystem_more_than_5,
                                                                                       number_cv_samples = 5, # this number is required because it is only possible to do 5-fold cv
                                                                                       leave_n_out = 2,
                                                                                       org = 'yeast', 
                                                                                       subsystem_only = "subsystem", 
                                                                                       alpha_val = 0, 
                                                                                       lno_cv = FALSE,
                                                                                       blocking_structure = 'by_growth_rate')

## Save these model outputs because each takes a while to run.
saveRDS(all_yeast_fluxes_predicted_list_subsystem, file = 'data/intermediate_data/all_yeast_fluxes_predicted_list_subsystem.rds')
saveRDS(all_yeast_fluxes_predicted_list_subsystem_ridge, file = 'data/intermediate_data/all_yeast_fluxes_predicted_list_subsystem_ridge.rds')
saveRDS(all_yeast_fluxes_predicted_list_subsystem_ridge_group_cv, file = 'data/intermediate_data/all_yeast_fluxes_predicted_list_subsystem_ridge_group_cv.rds')
saveRDS(all_yeast_fluxes_predicted_list_subsystem_ridge_z, file = 'data/intermediate_data/all_yeast_fluxes_predicted_list_subsystem_ridge_z.rds')
saveRDS(all_yeast_fluxes_predicted_list_subsystem_ridge_group_cv_growth, 
        file = 'data/intermediate_data/all_yeast_fluxes_predicted_list_subsystem_ridge_group_cv_growth.rds')
saveRDS(all_yeast_fluxes_predicted_list_subsystem_ridge_leave_5, file = 'data/intermediate_data/all_yeast_fluxes_predicted_list_subsystem_ridge_leave_5.rds')

}

## Reading in model outputs (previously saved).
all_yeast_fluxes_predicted_list_subsystem <- readRDS(file = 'data/intermediate_data/all_yeast_fluxes_predicted_list_subsystem.rds')
all_yeast_fluxes_predicted_list_subsystem_ridge <- readRDS(file = 'data/intermediate_data/all_yeast_fluxes_predicted_list_subsystem_ridge.rds')
all_yeast_fluxes_predicted_list_subsystem_ridge_group_cv <- readRDS(file = 'data/intermediate_data/all_yeast_fluxes_predicted_list_subsystem_ridge_group_cv.rds')
all_yeast_fluxes_predicted_list_subsystem_ridge_z <- readRDS(file = 'data/intermediate_data/all_yeast_fluxes_predicted_list_subsystem_ridge_z.rds')
all_yeast_fluxes_predicted_list_subsystem_ridge_group_cv_growth <- readRDS(file = 'data/intermediate_data/all_yeast_fluxes_predicted_list_subsystem_ridge_group_cv_growth.rds')
all_yeast_fluxes_predicted_list_subsystem_ridge_leave_5 <- readRDS(file = 'data/intermediate_data/all_yeast_fluxes_predicted_list_subsystem_ridge_leave_5.rds')

## calculating the rmsd and r squared for the pathway-level predictions
subsys_only_r_squared <- convert_model_summaries(all_yeast_fluxes_predicted_list_subsystem[[1]])

## calculating the rmsd and r squared for the pathway-level predictions, but using the ridge regression
subsys_only_r_squared_ridge <- convert_model_summaries(all_yeast_fluxes_predicted_list_subsystem_ridge[[1]])

## calculating the rmsd and r squared for the pathway-level predictions, but using the ridge regression w/ normalized proteins
subsys_only_r_squared_ridge_z <- convert_model_summaries(all_yeast_fluxes_predicted_list_subsystem_ridge_z[[1]])

## calculating the rmsd and r squared for the pathway-level predictions, but using the ridge regression w/ batch cross validation
subsys_only_r_squared_ridge_group_cv <- convert_model_summaries(all_yeast_fluxes_predicted_list_subsystem_ridge_group_cv[[1]])

## calculating the rmsd and r squared for the pathway-level predictions, but using the ridge regression w/ batch cross validation (growth rates)
subsys_only_r_squared_ridge_group_cv_growth <- convert_model_summaries(all_yeast_fluxes_predicted_list_subsystem_ridge_group_cv_growth[[1]])

## calculating the rmsd and r squared for the pathway-level predictions but using the ridge regression w/ leave 5 out cross validation
subsys_only_r_squared_ridge_leave_5 <- convert_model_summaries(all_yeast_fluxes_predicted_list_subsystem_ridge_leave_5[[1]])


