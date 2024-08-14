##### functions for LASSO regression analysis
library(glmnet)
library(scales)

get_all_protein_names_subsystem <- function(yeast_gem_df, 
                                            subsystem_name){
  # Function to get a list of all protein names that are in that subsystem
  
  # yeast_gem_df <- yeast_gem
  # subsystem_name <- "Glutathione metabolism"
  ## Subset just the subsystem
  yeast_gem_subsystem <- yeast_gem_df %>% 
    dplyr::filter(SUBSYSTEM == subsystem_name)
  ## get all the protein names associated with that subsystem
  vector_of_protein_names_no_format <- yeast_gem_subsystem$`GENE ASSOCIATION`
  ## The following code chunk formats the protein names
  prot_names_collapsed <- paste(vector_of_protein_names_no_format, collapse = " ")
  #### Remove the brackets
  prot_names_collapsed_no_brack <- gsub(x = prot_names_collapsed, pattern = "(", replacement = "", fixed = TRUE)
  prot_names_collapsed_no_brack2 <- gsub(x = prot_names_collapsed_no_brack, pattern = ")", replacement = "", fixed = TRUE)
  #### Get every string that begins with a Y or a Q (all genes)
  list_of_prots <- strsplit(x = prot_names_collapsed_no_brack2, split = " ")[[1]]
  all_subsystem_prot_names_y <- list_of_prots[grepl(pattern = "^Y", x = list_of_prots)]
  all_subsystem_prot_names_q <- list_of_prots[grepl(pattern = "^Q", x = list_of_prots)]
  
  ### Get just the unique protein names
  all_subsystem_prots <- unique(c(all_subsystem_prot_names_y, all_subsystem_prot_names_q))
  return(all_subsystem_prots)
}

## spot checking get all protein names from a subsystem
get_all_protein_names_subsystem(yeast_gem_df = yeast_gem,
                                subsystem_name = "Glutathione metabolism")

## spot checking get all protein names manually
# yeast_gem %>% dplyr::filter(SUBSYSTEM == 'Riboflavin metabolism') %>% 
#   dplyr::select(`GENE ASSOCIATION`) %>% as.vector() %>% unique()

get_proteins_per_subsystem <- function(yeast_gem_df, 
                                       vec_of_subsystem, 
                                       observed_proteins_vec){
  ### get a dataframe of protein numbers per subsystem, 
  ### to eventually subset only the subsystems that have at least 5 proteins
  ### loop through all subsystems
  proteins_per_sub <- c()

  for(i in 1:length(vec_of_subsystem)){
    subsystem_i <- vec_of_subsystem[i]
    get_all_proteins_in_sub_i <- get_all_protein_names_subsystem(yeast_gem_df = yeast_gem_df, 
                                                                 subsystem_name = subsystem_i)
    
    if(length(get_all_proteins_in_sub_i) == 0){
      get_all_proteins_in_sub_i <- "no proteins"
    }
    
    ## this will list all proteins, but many of these proteins were not detected.
    # loop through all these proteins, and if there is a string match found in the vector of names of proteins from the wide proteome
    # dataframe, then append it:
    get_all_proteins_in_sub_i_and_in_observations <- c()

    for(j in 1:length(get_all_proteins_in_sub_i)){
      ## If there is a string match
      if(grepl(pattern = get_all_proteins_in_sub_i[j], 
               x = observed_proteins_vec, 
               fixed = TRUE)){
        get_all_proteins_in_sub_i_and_in_observations <- c(get_all_proteins_in_sub_i_and_in_observations, 
                                                           get_all_proteins_in_sub_i[j])
      }
    }
    # print(get_all_proteins_in_sub_i_and_in_observations)
    proteins_per_subsystem <- length(get_all_proteins_in_sub_i_and_in_observations)
    
    proteins_per_sub <- c(proteins_per_sub, proteins_per_subsystem)
  }
  return(data.frame(vec_of_subsystem, proteins_per_sub))
}

get_single_protein_name <- function(reaction_name_choice, 
                                    hackett_data_processed = joined_df_flux_prot_mean_no_negative){
  # hackett_data_processed <- joined_df_flux_prot_mean_no_negative
  # reaction_name_choice <- "r_0061"
  gene_name_out <- hackett_data_processed[hackett_data_processed$reaction_id == reaction_name_choice, ]$Gene %>% unique()
  return(gene_name_out)
}

get_proteome_flux_dataframe <- function(reaction_name_choice, 
                                        gene_flux_data_frame, 
                                        wide_proteome_data_frame,
                                        org = 'yeast',
                                        filter_for_subsystem = "subsystem",
                                        yeast_gem_df = yeast_gem){
  
  ## Check the filter for subsystem argument
  if(filter_for_subsystem %!in% c('subsystem', 'single_prot')){
    print('Whole proteome used as predictors.')
  }
  
  ### This function gets the flux associated with a given reaction, 
  ### matches that with the full proteome profile, 
  ### and then outputs a dataframe.
  
  # if(org == 'e_coli'){
  #   # Filter for the flux profiles.
  #   df_temp <- gene_flux_data_frame %>% 
  #     dplyr::filter(reaction_name == reaction_name_choice) %>%
  #     dplyr::select(-protein_amount_mean) %>% 
  #     ## These conditions are removed
  #     dplyr::filter(exp_condition %!in% c('GAL_BATCH_mu=0.26_S', 'GAM_BATCH_mu=0.46_S', 
  #                                         'GLC+SALTS_BATCH_mu=0.55_S', 'XYL_BATCH_mu=0.55_S',
  #                                         'GLC_CHEM_mu=0.26_P', 'GLC_CHEM_mu=0.46_P'))
  # }
  if(org == 'yeast'){
    # Filter for the flux profiles.
    df_temp <- gene_flux_data_frame %>% 
      dplyr::filter(reaction_name == reaction_name_choice)
  }
  
  # Join with the full proteome.
  joined_df_protein_profile <- df_temp %>% 
    inner_join(wide_proteome_data_frame, 
    # inner_join(all_raw_renamed_trans_dcast, 
               by = "exp_condition")

  # if(org == 'e_coli'){
  #   # Remove all the additional information in this dataframe -- you want one
  #   joined_df_protein_profile_rf_format <- joined_df_protein_profile %>%
  #     dplyr::select(-reaction_name, -gene_name, -protein_name,
  #                   -formula_left_side, -exp_condition)
  # }
  if(org == 'yeast'){
    joined_df_protein_profile_rf_format <- joined_df_protein_profile %>% 
      dplyr::select(-reaction_name, -Reaction, -exp_condition_qp, -exp_condition)
  }
  
  ## if doing the yeast subsystem only dataframe, need to further filter.
  if(filter_for_subsystem == "subsystem"){
    ## get reaction name subsystem
    reaction_specific_subsystem <- yeast_gem_df[yeast_gem_df$ID == reaction_name_choice, ]$SUBSYSTEM
    
    ## get corresponding proteins for this subsystem
    all_proteins_from_subsys <- get_all_protein_names_subsystem(yeast_gem_df = yeast_gem_df, 
                                                                subsystem_name = reaction_specific_subsystem)
    
    # Select just the flux column and the proteins in that reaction subsystem
    joined_df_protein_profile_rf_format <- joined_df_protein_profile_rf_format[names(joined_df_protein_profile_rf_format) %in% c('flux', all_proteins_from_subsys)] 
    
  }
  
  ### if doing the yeast subsystem defined by KEGG, need to further filter.
  # if(filter_for_subsystem == 'kegg_subsystem'){
  #   ## get the protein names associated with this reaction rate
  #   
  #   ## determine the kegg pathway attributed to these 
  #   
  #   ## get all protein names corresponding to this subsystem
  #   
  #   ## select the protein names (be careful about the protein names from Hackett dataset that have a slash in them!)
  #   
  # }
  
  ## if doing the yeast single protein analysis
  if(filter_for_subsystem == "single_prot"){

    ## get corresponding single protein
    single_protein_name <- get_single_protein_name(reaction_name_choice = reaction_name_choice)
    
    # Select just the flux column and the proteins in that reaction subsystem
    joined_df_protein_profile_rf_format <- joined_df_protein_profile_rf_format[names(joined_df_protein_profile_rf_format) %in% c('flux', single_protein_name)] 
  }

  return(joined_df_protein_profile_rf_format)
}
# gene_assoc_test <- yeast_gem[yeast_gem$ID == 'r_0005', ]$`GENE ASSOCIATION`

## check all gene pathways

# grepl(pattern = yeast_pathway_mappings$gene_name[1], x = gene_assoc_test)
# grepl(pattern = "YPR165W", x = gene_assoc_test)
# 
# pathway_names_associated <- c()
# for(i in 1:length(yeast_pathway_mappings$gene_name)){
#   # i <- 1
#   pathway_within_map <- grepl(pattern = yeast_pathway_mappings$gene_name[i], 
#                                 x = gene_assoc_test)
#   if(pathway_within_map){
#     pathway_names_associated <- c(pathway_names_associated, 
#                                   yeast_pathway_mappings$kegg_pathway[i])
#     print(i)
#   }
# }

### testing hackett_fluxes_subset, 
# tester_out <- get_proteome_flux_dataframe(reaction_name_choice = "r_0005", 
#                             gene_flux_data_frame = hackett_fluxes_subset, 
#                             wide_proteome_data_frame = hackett_prot_transformed)
# scale(tester_out, center = FALSE, scale = TRUE)

get_two_layer_cv_output <- function(formatted_proteome_flux_df, 
                                    reaction_name_choice, 
                                    number_cv_samples = 100,
                                    leave_n_out = 5,
                                    alpha_val = 1,
                                    z_score_normalize = FALSE,
                                    lno_cv = TRUE,
                                    blocking_structure = 'by_condition'){
  
  # formatted_proteome_flux_df <- formatted_proteome_flux_df_tester
  # reaction_name_choice <- "r_0439"
  
  ## remove the flux column
  full_df_no_flux <- formatted_proteome_flux_df %>% 
    dplyr::select(-flux)
  
  ## subset the fluxes as a vector
  flux_vector <- formatted_proteome_flux_df$flux
  
  ### if you need to normalize the proteins for ridge regression, normalize predictors by standard deviation
  if(z_score_normalize){
    full_df_no_flux <- scale(full_df_no_flux) %>% as.data.frame()
    flux_vector <- scale(flux_vector)[,1]
  }
  
  ## setting up empty vectors for aggregating CV samples
  y_predicted_vec <- c()
  test_fluxes_vec <- c()
  mse_vec <- c()
  sample_i_vec <- c()
  test_index_vec <- c()
  coef_number_vec <- c()
  
  ## setting up empty dataframe
  coef_output_df_full <- data.frame(coefficient_values = numeric(),
                                    coefficient_names = numeric())
  
  # loop through all the Cross Validation iterations
  for(i in 1:number_cv_samples){
    
    # total number of observations (from the flux side)
    num_total <- length(flux_vector)
    
    ## Get training data indices and test data indices. This is ** leave x out cross validation **.
    if(lno_cv){
      train_data_indices <- sample(x = c(1:num_total), 
                                   size = num_total - leave_n_out, 
                                   replace = FALSE)
    }
    ## This is group indexed cross validation
    if(!lno_cv){
      ## test that number_cv_samples == 5
      testthat::expect_equal(number_cv_samples, 5)
      
      ## vector of group indices 
      if(blocking_structure == 'by_condition'){
        group_indices <- rep(c(1:5), 
                             each = 5)
      }
      ## group indicies are by growth rate, or by substrate limiting growth
      if(blocking_structure == 'by_growth_rate'){
        group_indices <- rep(c(1:5), 
                             times = 5)
      }

      row_vec <- c(1:25)
      
      ## subsample row vector in chunks of 5 
      ## (5 experimental conditions in the Hackett data, 5 growth rates each)
      train_data_indices <- row_vec[group_indices != i]
    }

    ## The test data indices are those that are not in the train data indices
    test_data_indices <- c(1:num_total)[c(1:num_total) %!in% train_data_indices]
    
    ## Getting the training and testing subsets for fluxes.
    train_fluxes <- flux_vector[train_data_indices]
    test_fluxes <- flux_vector[test_data_indices]
    
    ## Getting the training and testing subsets for proteomes.
    train_prots <- full_df_no_flux[train_data_indices, ]
    test_prots <- full_df_no_flux[test_data_indices, ]
    
    ## Perform k-fold cross-validation to find optimal lambda value
    cv_model <- cv.glmnet(y = train_fluxes, 
                          x = makeX(train_prots), 
                          nfolds = 10, 
                          alpha = alpha_val) ## when alpha = 1, LASSO regression | when alpha = 0, ridge regression

    # Find optimal lambda value that minimizes test MSE
    best_lambda <- cv_model$lambda.min

    # Use fitted best model to make predictions
    y_predicted <- predict(cv_model, 
                           s = best_lambda, 
                           newx = makeX(test_prots))

    # Find mean squared error
    mse_val_i <- sum((y_predicted - test_fluxes)^2)/length(y_predicted)
    
    # getting the model coefficients
    coef_test <- coef(cv_model, 
                      s = best_lambda)
    
    ## Number of non-zero coefficients
    non_zero_coefficients <- coef_test@x %>% length()
    
    ## Making a dataframe of non-zero coefficients
    coefficient_values <- coef_test@x
    coefficient_indices <- coef_test@i
    
    ## Need to add one to this because the indices are zerod.
    coefficient_names <- coef_test@Dimnames[[1]][coefficient_indices + 1]
    
    ## Making a coefficient value to dataframe mapping
    coef_output_df_i <- data.frame(coefficient_values,
                                   coefficient_names)
    
    # Aggregating summaries.
    y_predicted_vec <- c(y_predicted_vec, y_predicted)
    test_fluxes_vec <- c(test_fluxes_vec, test_fluxes)
    mse_vec <- c(mse_vec, mse_val_i)
    sample_i_vec <- c(sample_i_vec, rep(i, length(y_predicted)))
    test_index_vec <- c(test_index_vec, test_data_indices)
    coef_number_vec <- c(coef_number_vec, non_zero_coefficients)
    coef_output_df_full <- rbind(coef_output_df_full, coef_output_df_i)
  }
  
  ## Aggregating all this into one dataframe for the chosen flux value.
  df_out <- data.frame(pred_flux = y_predicted_vec,
                       actual_flux = test_fluxes_vec,
                       cv_sample_id = sample_i_vec,
                       flux_index_val = test_index_vec)
  
  ## Appending a column that identifies the reaction name.
  df_out$reaction_name <- rep(reaction_name_choice, nrow(df_out))
  
  ## Making a dataframe that summaries each model result, rather than the individual predictions
  df_out_summaries <- data.frame(coef_number = coef_number_vec,
                                 mse = mse_vec)
  
  # Add a column that has reaction name
  df_out_summaries$reaction_name <- rep(reaction_name_choice, 
                                        nrow(df_out_summaries))
  
  ## Add a column that has reaction name
  coef_output_df_full$reaction_name <- rep(reaction_name_choice, 
                                           nrow(coef_output_df_full))
  
  return(list(df_out, 
              df_out_summaries, 
              coef_output_df_full))
}

loop_through_and_cv_fluxes <- function(gene_flux_data_frame, 
                                       wide_proteome_data_frame, 
                                       reaction_name_vector, 
                                       number_cv_samples,
                                       org = 'yeast',
                                       leave_n_out = 5,
                                       subsystem_only = "whole-proteome",
                                       alpha_val = 1,
                                       z_score_normalize = FALSE,
                                       lno_cv = TRUE,
                                       blocking_structure = 'by_condition'){
  
  ### Make a dataframe to add values to for predicted and actual flux values
  df_supreme_pred_actual <- data.frame(pred_flux = numeric(),
                                       actual_flux = numeric(),
                                       cv_sample_id = numeric(),
                                       flux_index_val = numeric(),
                                       reaction_name = character())
  
  ## Make a dataframe to add values for the MSE and coef number per CV model
  df_supreme_coef_mse <- data.frame(coef_number = integer(),
                                    mse = numeric())
  
  ## Make a dataframe to add values for the coefficients
  df_supreme_coef_vals <- data.frame(coefficient_values = numeric(),
                                     coefficient_names = character())
  
  ### Go through all fluxes and get the CV model output
  for(unique_reaction_i in 1:length(reaction_name_vector)){
    
    print(paste0('Getting CV model for reaction: ', reaction_name_vector[unique_reaction_i]))
    
    ## First format the dataframe
    formatted_proteome_flux_df <- get_proteome_flux_dataframe(reaction_name_choice = reaction_name_vector[unique_reaction_i], 
                                                              gene_flux_data_frame = gene_flux_data_frame,
                                                              wide_proteome_data_frame = wide_proteome_data_frame,
                                                              org = org, 
                                                              filter_for_subsystem = subsystem_only)
    
    ## Get the dataframe of predictions versus actual fluxes
    cv_reaction_i_list <- get_two_layer_cv_output(formatted_proteome_flux_df = formatted_proteome_flux_df, 
                                                  reaction_name_choice = reaction_name_vector[unique_reaction_i], 
                                                  number_cv_samples = number_cv_samples,
                                                  leave_n_out = leave_n_out, 
                                                  alpha_val = alpha_val,
                                                  z_score_normalize = z_score_normalize, 
                                                  lno_cv = lno_cv,
                                                  blocking_structure = blocking_structure)
    
    ## Append this output to the CV dataframe
    df_supreme_pred_actual <- rbind(df_supreme_pred_actual, cv_reaction_i_list[[1]])
    
    ## Append this output to the model summary dataframe
    df_supreme_coef_mse <- rbind(df_supreme_coef_mse, cv_reaction_i_list[[2]])
    
    ## Appednd this output to the coefficient summary dataframe
    df_supreme_coef_vals <- rbind(df_supreme_coef_vals, cv_reaction_i_list[[3]])
  }
  
  return(list(df_supreme_pred_actual, 
              df_supreme_coef_mse, 
              df_supreme_coef_vals))
  
}

convert_model_summaries <- function(model_output){
  
  # model_output <- subset_yeast_fluxes_singles
  ### this calculates model summaries across reaction names
  if(ncol(model_output) == 4){
    testthat::expect_equal(object = names(model_output), expected = c("predicted_fluxes", 
                                                                      "actual_fluxes", 
                                                                      "cv_sample", 
                                                                      "reaction_rate_id"))
    model_summary_out <- model_output %>% 
      group_by(reaction_rate_id) %>% 
      summarize(mean_actual = mean(actual_fluxes),
                ss_total = sum((actual_fluxes - mean_actual)^2),
                ss_residuals = sum((predicted_fluxes - actual_fluxes)^2),
                number_obs = n(),
                rmsd = sqrt(sum(ss_residuals)/number_obs)) %>% 
      mutate(r_sq = 1 - ss_residuals/ss_total)
  }
  if(ncol(model_output) == 5){
    testthat::expect_equal(object = names(model_output), expected = c("pred_flux", "actual_flux", 
                                                                      "cv_sample_id", "flux_index_val", 
                                                                      "reaction_name"))
    model_summary_out <- model_output %>% 
      group_by(reaction_name) %>% 
      summarize(mean_actual = mean(actual_flux),
                ss_total = sum((actual_flux - mean_actual)^2),
                ss_residuals = sum((pred_flux - actual_flux)^2),
                number_obs = n(),
                rmsd = sqrt(sum(ss_residuals)/number_obs)) %>% 
      mutate(r_sq = 1 - ss_residuals/ss_total) %>% 
      inner_join(yeast_gem %>% 
                   dplyr::select(ID, SUBSYSTEM) %>% 
                   dplyr::rename(reaction_name = ID)) 
  }
  
  return(model_summary_out)
}

# plotting functions ------------------------------------------------------

plot_subsystem_specific_reactions <- function(subsystem_name, prediction_df){
  ## function to plot all the specific predictions from cross validation, subset by subsystem
  
  # join the prediction dataframe with the subsystem id using the reaction name
  prediction_df_w_subsystems <- prediction_df %>% 
    inner_join(yeast_gem %>% 
                 dplyr::select(ID, 
                               SUBSYSTEM) %>% 
                 dplyr::rename(reaction_name = ID))
  
  # make the figure
  p_out <- prediction_df_w_subsystems %>%
    dplyr::filter(SUBSYSTEM == subsystem_name) %>% 
    group_by(flux_index_val, reaction_name) %>%
    summarize(mean_pred_flux = mean(pred_flux),
              sd_pred_flux = sd(pred_flux),
              mean_actual_flux = mean(actual_flux)) %>%
    ggplot(aes(y = mean_pred_flux, x = mean_actual_flux)) +
    # geom_point() +
    facet_wrap(~reaction_name) +
    # geom_errorbar(aes(ymin = mean_pred_flux - sd_pred_flux/2, ymax = mean_pred_flux + sd_pred_flux/2)) +
    xlab('Observed Rate') +
    ylab('LASSO-Predicted Rate') +
    theme_bw()
  
  return(p_out)
}

## need a function that summarizes the ggplot axes 
### Credit to: https://stackoverflow.com/questions/10762287/how-can-i-format-axis-labels-with-exponents-with-ggplot2-and-scales
scientific_10 <- function(x) {
  ifelse(
    x==0, "0",
    parse(text = sub("e[+]?", " %*% 10^", scientific_format()(x)))
  )
} 
