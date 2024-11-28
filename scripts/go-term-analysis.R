## go analysis

library(forcats)

# # Get a list of all proteins identified, to then search for GO t --------

all_hackett_proteins_w_extra <- names(hackett_prot_transformed)
all_hacket_proteins <- all_hackett_proteins_w_extra[all_hackett_proteins_w_extra != 'exp_condition']

write.table(x = all_hacket_proteins, file = "data/list_of_yeast_genes.txt", 
            sep = ' ', quote = FALSE, row.names = FALSE, col.names = FALSE)

# Parse the GO term tab output --------------------------------------------

go_process <- read.table(file = 'data/yeast_go_process_mapping.txt', 
                         sep = '\t', header = TRUE)

## change the names
names(go_process) <- tolower(names(go_process))

expand_go_ids_into_df <- function(go_process_df, row_number){
  
  ## need a function that parses this into a reasonable dataframe that's tidy
  ### This one takes in a go process/term dataframe, and a row number, and expands that dataframe
  # go_process_df <- go_component
  annotated_term <- go_process_df[row_number, ]$term
  annotated_genes_i <- go_process_df[row_number, ]$annotated_genes
  
  ## parse these, first remove spaces between genes
  removed_spaces_genes <- str_split(annotated_genes_i, pattern = ' ')[[1]]
  ## then remove the commas
  removed_commas_genes <- gsub(pattern = ',', replacement = '', x = removed_spaces_genes)
  
  ## make these into a tidy df
  df_temp <- data.frame(go_term = rep(annotated_term, length(removed_commas_genes)),
                        gene = removed_commas_genes)
  return(df_temp)
}

# expand_go_ids_into_df(go_process_df = go_component, row_number = 1)
go_through_go_expand_all <- function(go_process_df){
  ### loop through all the go process df and expand each row
  temp_df_out <- data.frame(go_term = character(),
                            gene = character())
  for(row_i in 1:nrow(go_process_df)){
    temp_df_i <- expand_go_ids_into_df(go_process_df = go_process_df, 
                                       row_number = row_i)
    temp_df_out <- rbind(temp_df_out, 
                         temp_df_i)
  }
  return(temp_df_out)
}

go_process_tidy <- go_through_go_expand_all(go_process_df = go_process)

## make a dataframe that has all of the unique go terms, alongside the type of go term it is
unique_process_terms <- unique(go_process_tidy$go_term)

## aggregate these into a dataframe
go_big_terms <- data.frame(go_term = c(unique_process_terms),
                           go_category = c(rep('process', length(unique_process_terms))))

# formatting data from lasso models ---------------------------------------

# all_yeast_fluxes_predicted_coefs is the key dataframe, from the whole proteome predictions
# Need a function to get the number of terms associated with a specific GO term in the LASSO model (N1)

## first need to format this output. But, we want to filter for only the predictions that had an Cross validated R^2 of greater
## than 0.5. all_data_subsys_r_squared is the dataframe with this corresponding information.

## First filter for the R^2 above 0.5.
all_data_subsys_r_squared_r_sq_above_0.5 <- all_data_subsys_r_squared %>% 
  dplyr::filter(r_sq > 0.5)

## Then only choose those coefficients in models that were predictive above a 0.5 threshold.
all_yeast_fluxes_predicted_coefs_formatted <- all_yeast_fluxes_predicted_coefs %>% 
  dplyr::filter(reaction_name %in% all_data_subsys_r_squared_r_sq_above_0.5$reaction_name)

## Formatting the gene names (in the predictive models they have apostrophes)
all_yeast_fluxes_predicted_coefs_formatted$gene <- gsub("`", '', all_yeast_fluxes_predicted_coefs_formatted$coefficient_names)

## Remove the intercept term.
all_yeast_fluxes_predicted_coefs_formatted2 <- all_yeast_fluxes_predicted_coefs_formatted %>% 
  dplyr::filter(gene != '(Intercept)')

## Make a vector of all the genes only. 
all_yeast_fluxes_predicted_coefs_formatted2_genes_only <- all_yeast_fluxes_predicted_coefs_formatted2$gene

### now need to format these -- duplicate those values that are called a single value 'GENEX/GENEY' becomes two.

## for this task we have a function.
format_hackett_proteins_with_slashes <- function(all_hacket_proteins){
  ### this function takes in the list of proteins and then spits out the protein list, but those
  ### proteins that were in one line (because of slashes) are now in two lines)
  empty_protein_list <- c()
  
  for(i in 1:length(all_hacket_proteins)){
    print(i)
    # get the ith protein
    protein_i <- all_hacket_proteins[i]
    
    # if it contains a back slash
    if(grepl(pattern = '/', x = protein_i) == TRUE){
      # split that up
      protein_sub <- str_split(string = protein_i, pattern = '/')[[1]]
    }
    else {
      protein_sub <- protein_i
    }
    ## aggregate it all
    empty_protein_list <- c(empty_protein_list, protein_sub)
  }
  return(empty_protein_list)
}

## Apply the function.
formatted_predictors <- format_hackett_proteins_with_slashes(all_hacket_proteins = all_yeast_fluxes_predicted_coefs_formatted2_genes_only)

## Apply this formatting to all the observed proteins. This is the baseline from which we compare.
all_hackett_proteins_formatted <- format_hackett_proteins_with_slashes(all_hacket_proteins)

## Need a function that assigns GO terms to the list of proteins.
target_go_class_df <- function(all_hacket_proteins = all_hackett_proteins_formatted, go_class){
  
  hackett_proteins_df <- data.frame(gene = all_hacket_proteins) 
  
  if(go_class == 'process'){
    df_out <- hackett_proteins_df %>% 
      inner_join(go_process_tidy, by = 'gene')
  }
  if(go_class == 'function'){
    df_out <- hackett_proteins_df %>% 
      inner_join(go_function_tidy, by = 'gene')
  }
  if(go_class == 'component'){
    df_out <- hackett_proteins_df %>% 
      inner_join(go_component_tidy, by = 'gene')
  }
  
  return(df_out)
}

## Get total terms after making a GO-term joined dataframe. This is for the denominator in the proportion.
get_total_terms <- function(all_hacket_proteins, go_class){
  
  target_df <- target_go_class_df(all_hacket_proteins = all_hacket_proteins, 
                                  go_class = go_class)
  
  return(nrow(target_df))
  
}

## Get a single term terms after making a GO-term joined dataframe. This is for the denominator in the proportion.
get_single_term_from_proteome <- function(all_hacket_proteins, go_class, single_term){
  ## this function finds the number of instances of a specific term in the group
  target_df <- target_go_class_df(all_hacket_proteins = all_hacket_proteins, 
                                  go_class = go_class)
  ## filter for the specified term.
  filtered_df <- target_df %>% 
    dplyr::filter(go_term == single_term)
  
  if(nrow(filtered_df) == 0){
    stop('the term you supplied isnt valid')
  } 
  return(nrow(filtered_df))
}


extract_data_for_binomial <- function(go_class, single_term, 
                                      all_hacket_proteins = all_hackett_proteins_formatted,
                                      predictor_coefs = formatted_predictors){
  ## this function aggregates the above functions to find
  ## 1) N1, the number of a single term from the observed proteome.
  N1_out <- get_single_term_from_proteome(all_hacket_proteins = all_hacket_proteins, 
                                go_class = go_class, 
                                single_term = single_term)
  ## 2) N2, the number of a single term from the predictive models.
  N2_out <- get_single_term_from_proteome(all_hacket_proteins = predictor_coefs, 
                                go_class = go_class, 
                                single_term = single_term)
  ## 3) Y1, the total number of terms in a specific GO class (e.g., function, process, etc.) for the observed proteome.
  Y1_out <- get_total_terms(all_hacket_proteins = all_hackett_proteins_formatted, 
                  go_class = go_class)
  ## 4) Y2, the total number of terms in a specific GO class, for the predictive models.
  Y2_out <- get_total_terms(all_hacket_proteins = predictor_coefs, 
                  go_class = go_class)
  
  ## return in a format that can be read into stan.
  return(list(N1 = N1_out, Y1 = Y1_out, N2 = N2_out, Y2 = Y2_out))
}

# output_test <- extract_data_for_binomial(go_class = 'component', single_term = 'nucleus')

## need a function to format the stan output
get_rhats_for_each_par <- function(stan_model_out){
  return(summary(stan_model_out)$summary %>% as.data.frame())
}

run_stan_model_extract_stats <- function(go_class, single_term){
  
  ## extract data according to the chosen go class and the single term examined
  output_from_subset <- extract_data_for_binomial(go_class = go_class, 
                                                  single_term = single_term)
  
  ## fitting the model in stan
  model_fit_out <- sampling(binomial_model, 
                            data = output_from_subset, 
                            iter = 3000, chains = 4)
  
  ## get the model summary
  model_summary_out <- get_rhats_for_each_par(model_fit_out)
  
  ## append the model summary with the name of the go class and the single term examined
  model_summary_out$go_class <- rep(go_class, nrow(model_summary_out))
  model_summary_out$go_term <- rep(single_term, nrow(model_summary_out))
  model_summary_out$par_out <- rownames(model_summary_out)
  
  return(model_summary_out)
}

# tester_out <- run_stan_model_extract_stats(go_class = 'component', single_term = 'mitochondrion')

## this function loops through all combinations of go categories, and specific go terms, then fits a model, and
## returns the formatted model.
loop_through_go_classes_and_terms <- function(go_big_terms){
  group_df <- data.frame(mean = numeric(),
                         se_mean = numeric(),
                         sd = numeric(),
                         `2.5%` = numeric(),
                         `25%` = numeric(),
                         `50%` = numeric(),
                         `75%` = numeric(),
                         `97.5%` = numeric(),
                         n_eff = numeric(),
                         Rhat = numeric(),
                         go_class = character(),
                         go_term = character(),
                         par_out = character())
  for(i in 1:nrow(go_big_terms)){
    
    print(i)
    
    go_class_i <- go_big_terms[i, ]$go_category
    go_term_i <- go_big_terms[i, ]$go_term
    fit_out_i <- run_stan_model_extract_stats(go_class = go_class_i, 
                                              single_term = go_term_i)
    group_df <- rbind(group_df, 
                      fit_out_i)
  }
  return(group_df)
}

## read in the stan model.
binomial_model <- stan_model("scripts/binomial-difs.stan")

## this fits all the stan models. Takes around 30 mins.
go_big_or_go_home <- loop_through_go_classes_and_terms(go_big_terms = go_big_terms)

# examining model outputs and plotting -------------------------------------------------
  
go_enrichment_process <- go_big_or_go_home %>% 
  dplyr::filter(par_out == 'theta_difference') %>% 
  dplyr::filter(go_class == 'process') %>% 
  ggplot(aes(x = `50%`, y = fct_reorder(go_term, `50%`))) +
  geom_point() +
  geom_linerange(aes(xmin = `25%`, xmax = `75%`), lwd = 2, alpha = 0.4) +
  geom_linerange(aes(xmin = `2.5%`, xmax = `97.5%`), lwd = 1, alpha = 0.2) +
  geom_vline(xintercept = 0) +
  theme_bw() +
  ylab('GO Process Term') +
  theme(text = element_text(size = 16)) +
  xlab('Enrichment for GO Term\n(Posterior Probability of Proportion in LASSO Predictions - \n Proportion in Whole Proteome)')

ggsave(go_enrichment_process, filename = paste0('figures/figsX_', current_date, '_go_enrichment_process.pdf'),
       height = 13, width = 11.6)

