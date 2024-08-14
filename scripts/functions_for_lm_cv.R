# functions for lm cv

## write functions for linear model with only focal protein
get_cv_linear_model_one_coef <- function(formatted_proteome_flux_df,
                                         reaction_name_choice, 
                                         number_cv_samples = 100,
                                         leave_n_out = 2){
  ### this function does a cross validation approach as above, on a linear model with only 
  ### one term that is fit using OLS
  
  # formatted_proteome_flux_df <- joined_df_flux_prot_mean_no_negative %>% 
  #   dplyr::select(reaction_id, flux, mean_protein_amount) %>% 
  #   mutate(protein_amount = 2^mean_protein_amount) %>% 
  #   dplyr::filter(reaction_id == "r_0892")
  # formatted_proteome_flux_df <- formatted_test
  
  # formatted_proteome_flux_df <- df_temp
  
  ## remove the flux column
  protein_values <- formatted_proteome_flux_df %>% 
    dplyr::select(-flux)
  
  ## subset the fluxes as a vector
  flux_vector <- formatted_proteome_flux_df$flux
  
  ## setting up empty vectors for aggregating CV samples
  y_predicted_vec <- c()
  test_fluxes_vec <- c()
  cv_sample_id <- c()
  
  for(i in 1:number_cv_samples){
    
    # total number of observations (from the flux side)
    num_total <- length(flux_vector)
    
    ## Get training data indices and test data indices. This is ** leave 2 out cross validation **.
    train_data_indices <- sample(x = c(1:num_total), 
                                 size = num_total - leave_n_out, 
                                 replace = FALSE)
    ## The test data indices are those that are not in the train data indices
    test_data_indices <- c(1:num_total)[c(1:num_total) %!in% train_data_indices]
    
    ## Getting the training and testing subsets for fluxes.
    train_fluxes <- flux_vector[train_data_indices]
    test_fluxes <- flux_vector[test_data_indices]
    
    ## Getting the training and testing subsets for proteomes.
    train_prots <- protein_values[train_data_indices, ]
    test_prots <- protein_values[test_data_indices, ]
    
    ## fit the linear model
    lm_out <- lm(train_fluxes ~ train_prots)
    
    ## predict on new data
    predictions_i <- predict(object = lm_out, newdata = data.frame(train_prots = test_prots), type = "response")
    
    ## aggregate predictions and actual values
    y_predicted_vec <- c(y_predicted_vec, predictions_i)
    test_fluxes_vec <- c(test_fluxes_vec, test_fluxes)
    cv_sample_id <- c(cv_sample_id, rep(i, length(test_fluxes)))
  }
  
  df_out <- data.frame(predicted_fluxes = y_predicted_vec,
                       actual_fluxes = test_fluxes_vec,
                       cv_sample = cv_sample_id)
  
  return(df_out)
}

get_proteome_flux_dataframe_e_coli <- function(reaction_name_choice, 
                                               protein_gene_df){

  # protein_gene_df <- genes_fluxes_long_df_z_score
  # reaction_name_choice <- "ADSS"
  
  ### This function gets the flux associated with a given reaction, 
  ### matches that with the single protein 
  ### and then outputs a dataframe.
  
  df_temp <- protein_gene_df %>% 
    dplyr::filter(reaction_name == reaction_name_choice) %>%
    ## These conditions are removed
    dplyr::filter(exp_condition %!in% c('GAL_BATCH_mu=0.26_S', 'GAM_BATCH_mu=0.46_S', 
                                        'GLC+SALTS_BATCH_mu=0.55_S', 'XYL_BATCH_mu=0.55_S',
                                        'GLC_CHEM_mu=0.26_P', 'GLC_CHEM_mu=0.46_P')) %>% 
    dplyr::select(flux, protein_amount_mean) %>% 
    dplyr::mutate(log_protein_amount_mean = log2(protein_amount_mean)) %>% 
    dplyr::select(-protein_amount_mean) %>% as.data.frame()

  return(df_temp)
}

get_proteome_flux_dataframe_b_sub <- function(reaction_name_choice, 
                                               protein_gene_df){
  
  # protein_gene_df <- chu_transformed_cc
  # reaction_name_choice <- "pfk"
  
  ### This function gets the flux associated with a given reaction, 
  ### matches that with the single protein 
  ### and then outputs a dataframe.
  
  df_temp <- protein_gene_df %>% 
    dplyr::filter(flux_name == reaction_name_choice) %>%
    dplyr::select(flux, mean_expression_val) %>% 
    dplyr::mutate(log_mean_expression_val = log2(mean_expression_val)) %>% 
    dplyr::select(-mean_expression_val) %>% as.data.frame()
  
  return(df_temp)
}

# get_proteome_flux_dataframe_b_sub(reaction_name_choice = "gapa",
                                  # protein_gene_df = chu_transformed_cc)
# get_cv_linear_model_one_coef(formatted_proteome_flux_df = df_temp, reaction_name_choice = "ADSS", number_cv_samples = 10, leave_n_out = 2)


loop_through_and_cv_fluxes_single <- function(gene_flux_data_frame = hackett_fluxes_subset, 
                                              wide_proteome_data_frame = hackett_prot_transformed, 
                                              org_choice = "yeast",
                                              reaction_name_vector, 
                                              number_cv_samples, 
                                              leave_n_out = 5){
  ### this function loops through a vector of reation rate names, and then does a cross validation using 
  ### a *single protein* linear regression model.
  
  ## instantiate a dataframe to aggregate outputs
  supreme_df <- data.frame(predicted_fluxes = numeric(),
                           actual_fluxes = numeric(),
                           cv_sample = numeric(),
                           reaction_rate_id = numeric())
  
  ## loop through these rates
  for(i in 1:length(reaction_name_vector)){
    
    ## subset reaction rate i 
    reaction_rate_i <- reaction_name_vector[i]
    
    if(org_choice == "yeast"){
      ## first format the dataframe for a single protein-rate pair
      single_protein_dataframe_pair <- get_proteome_flux_dataframe(reaction_name_choice = reaction_rate_i, 
                                                                   gene_flux_data_frame = gene_flux_data_frame, 
                                                                   wide_proteome_data_frame = wide_proteome_data_frame,
                                                                   org = 'yeast',
                                                                   filter_for_subsystem = "single_prot",
                                                                   yeast_gem_df = yeast_gem)
    }
    if(org_choice == "ecoli"){
      single_protein_dataframe_pair <- get_proteome_flux_dataframe_e_coli(reaction_name_choice = reaction_rate_i, 
                                                                          protein_gene_df = gene_flux_data_frame)#genes_fluxes_long_df_z_score)
    }
    if(org_choice == "bsubtilis"){
      single_protein_dataframe_pair <- get_proteome_flux_dataframe_b_sub(reaction_name_choice = reaction_rate_i,
                                                                         protein_gene_df = gene_flux_data_frame)#chu_transformed_cc
    }

    
    ## do cross validation using just a simple linear model
    cv_lm_test <- get_cv_linear_model_one_coef(formatted_proteome_flux_df = single_protein_dataframe_pair, 
                                               reaction_name_choice = reacion_rate_i, 
                                               number_cv_samples = number_cv_samples, 
                                               leave_n_out = leave_n_out)
    
    ## add a new column to identify the rxn
    cv_lm_test$reaction_rate_id <- rep(reaction_rate_i, nrow(cv_lm_test))
    
    ## append to supreme df
    supreme_df <- rbind(supreme_df, cv_lm_test)
  }
  
  return(supreme_df)
}
