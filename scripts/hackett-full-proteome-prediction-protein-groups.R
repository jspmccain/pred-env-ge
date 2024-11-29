## looking at if protein groups predict reaction rates 

## first need to summarize protein groups by subsystem

## Need a dataframe that maps the protein names to the subsystem
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
# get_all_protein_names_subsystem(yeast_gem_df = yeast_gem,
#                                 subsystem_name = "Glutathione metabolism")

## get all the protein names
protein_names_whole_proteome <- hackett_prot$Gene

## loop through every protein name to retrieve the subsystem
retrieve_subsystem <- function(protein_name_in, yeast_gem_df = yeast_gem){
  ## filter the specific target protein
  sub_df_filtered <- yeast_gem_df %>% 
    dplyr::filter(grepl(pattern = protein_name_in, x = `GENE ASSOCIATION`))
  ## get the subsystem of interest
  subsystem_out <- sub_df_filtered$SUBSYSTEM
  ## return subsystem
  return(subsystem_out)
}

make_df_of_protein_subsystem_map <- function(list_of_protein_names, yeast_gem_df = yeast_gem){
  ## loop through all protein names
  subsystem_vector_out <- c()
  for(i in 1:length(list_of_protein_names)){
    ## get the unique subsystem output
    subsystem_i <- unique(retrieve_subsystem(protein_name_in = list_of_protein_names[i]))
    
    ## check if there is an associated subsystem
    if(length(subsystem_i) == 0){
      subsystem_i <- 'no-assigned-subsystem'
    }
    ## check if there are multiple assigned subsystems (57 in total from preliminary checks)
    if(length(subsystem_i) > 1){
      subsystem_i <- 'multi-subsystem-protein'
    }
  
    ## append the subsystem assignment
    subsystem_vector_out <- c(subsystem_vector_out, subsystem_i)
  }
  return(data.frame(protein_name = list_of_protein_names, subsystem = subsystem_vector_out))
}

## protein to subsystem dataframe
sub_prot_df <- make_df_of_protein_subsystem_map(list_of_protein_names = protein_names_whole_proteome)

## subset the protein to subsystem dataframe to only those proteins that have a single unique subsystem
sub_prot_df_unique <- sub_prot_df %>% 
  dplyr::filter(subsystem %!in% c('no-assigned-subsystem', 'multi-subsystem-protein'))

## make a dataframe that works with the cross validation scheme
hackett_prot_transformed_subsystem_ag <- hackett_prot %>% 
  ## join the subsystem groupings
  inner_join(sub_prot_df_unique %>% 
               dplyr::rename(Gene = protein_name), 
             by = "Gene") %>% 
  melt(variable.name = 'exp_condition', 
       value.name = 'protein_level') %>% 
  dplyr::select(-Gene) %>% 
  group_by(subsystem, exp_condition) %>% 
  summarize(mean_protein_level = mean(protein_level)) %>% 
  dcast(exp_condition ~ subsystem) %>% 
  mutate(exp_condition = fct_recode(exp_condition, "P0.22" = "p0.22"))


if(rerun_all_pred_models){
  ## run cross validation
  all_yeast_fluxes_predicted_list_subsystem_groups <- loop_through_and_cv_fluxes(gene_flux_data_frame = hackett_fluxes_subset, 
                                                                                 wide_proteome_data_frame = hackett_prot_transformed_subsystem_ag, 
                                                                                 reaction_name_vector = unique_hackett_fluxes, 
                                                                                 leave_n_out = 2,
                                                                                 number_cv_samples = 100, 
                                                                                 org = 'yeast')
  
  saveRDS(all_yeast_fluxes_predicted_list_subsystem_groups, file = 'data/intermediate_data/all_yeast_fluxes_predicted_list_subsystem_groups.rds')
}

## read the written results file in 
all_yeast_fluxes_predicted_list_subsystem_groups <- readRDS(file = 'data/intermediate_data/all_yeast_fluxes_predicted_list_subsystem_groups.rds')

## get the model summmaries
model_summaries_sub_groups <- convert_model_summaries(all_yeast_fluxes_predicted_list_subsystem_groups[[1]])

## make a figure of the distribution of cross validated r^2
subsystem_groupings_plot <- model_summaries_sub_groups %>% 
  ggplot(aes(x = r_sq)) +
  geom_histogram(fill = yeast_col) +
  theme_bw() +
  xlab(expression(Model~Performance~(Cross~Validated~R^2))) +
  xlim(-1, 1) +
  ggtitle('A. Subsystem-groupings as predictors, leave-2-out cross validation') +
  ylab('Count')

## save figure
ggsave(subsystem_groupings_plot, 
       filename = paste0("figures/figsX_", current_date, "_subsystem_groupings_as_predictors.pdf"),
       height = 9.5, 
       width = 9.5)

