library(ggplot2)
library(rstan)
library(dplyr)
library(magrittr)
library(readxl)
library(reshape2)
library(tidybayes)
library(tidyr)
library(ggpubr)
library(forcats)
library(stringr)
library(GGally)

'%!in%' <- function(x,y)!('%in%'(x,y))

format_yeast_reactants <- function(yeast_gem_formula){
  ### Function goal: output a vector of reactants, including both substrates and products for 
  ### reversible reactions for each reaction.
  
  ### Include both reactants and products if the reaction is reversible
  if(grepl("<=>", yeast_gem_formula)){
    left_side_only_yeast <- yeast_gem_formula
  } else {
    ### Include only reactants if the reaction is irreversible
    right_and_left_sides_yeast <- str_split(string = yeast_gem_formula, 
                                            pattern = ">")
    ## Then subset only the left side.
    left_side_only_yeast <- sapply(right_and_left_sides_yeast, 
                                   "[[", 
                                   1)
  }
  # print(yeast_gem_formula)
  ## Including a padded zero on the left side so that reactants can be searched using strings, and won't 
  ## get a sub-string. For example, 'moco[c]' would be counted as well as 'bmoco[c]'.
  return(paste0(" ", 
                left_side_only_yeast))
}

format_yeast_reactants_gem <- function(yeast_gem_formula){
  ### Function goal: output a vector of reactants, including both substrates and products for 
  ### reversible reactions for each reaction.
  
  vector_of_reactants <- c()
  for(i in 1:length(yeast_gem_formula)){
    ### Include both reactants and products if the reaction is reversible
    if(grepl("<=>", yeast_gem_formula[i])){
      left_side_only_yeast <- yeast_gem_formula[i]
    } else {
      ### Include only reactants if the reaction is irreversible
      right_and_left_sides_yeast <- str_split(string = yeast_gem_formula[i], 
                                              pattern = ">")
      ## Then subset only the left side.
      left_side_only_yeast <- sapply(right_and_left_sides_yeast, 
                                     "[[", 
                                     1)
    }
    vector_of_reactants <- c(vector_of_reactants, left_side_only_yeast)
  }
  
  single_vector_of_all_reactants <- paste(vector_of_reactants, collapse = " ")
  ## Including a padded zero on the left side so that reactants can be searched using strings, and won't 
  ## get a sub-string. For example, 'moco[c]' would be counted as well as 'bmoco[c]'.
  return(paste0(" ", 
                single_vector_of_all_reactants))
}

# reading in data ---------------------------------------------------------

### read in yeast GEM data to get rate ID to gene name mapping
yeast_gem <- read_excel(path = "data/yeast-GEM.xlsx", 
                        sheet = 1)

### read in yeast metabolites
all_yeast_metabolites <- read_excel(path = "data/yeast-GEM.xlsx", 
                                    sheet = 2)

## read in yeast gene name to name mappings
yeast_gene_names <- read_excel(path = "data/yeast-GEM.xlsx", 
                               sheet = 4)

if(rerun_all_pred_models){
  ## This is to get the yeast_genes_annots.csv data
  source('scripts/format_yeast_gene_names_from_fasta.R')
}

## the above yeast gene names are only those that are represented in the yeast gem, but we want 'em all
yeast_gene_names_from_fasta <- read.csv(file = "data/yeast_genes_annots.csv")

### read in data from Hackett et al 2016, Science
# Flux data
hackett_fluxes <- read_excel(path = 'data/aaf2786-hackett-sm-table-s9.xlsx', 
                             sheet = 5)
# Protein levels
hackett_prot <- read_excel(path = 'data/aaf2786-hackett-sm-table-s9.xlsx', 
                           sheet = 3)

### read in data for rxn delta G, from Yeast GEM
delta_g <- read.csv(file = 'data/model_rxnDeltaG.csv')

## read in protein length data
joined_df_flux_prot_mean_no_negative_mean_proteins_manual <- read.csv(file = "data/joined_df_flux_prot_mean_no_negative_mean_proteins-manual-edit.csv")

yeast_genes_protein_lengths <- joined_df_flux_prot_mean_no_negative_mean_proteins_manual %>%
  dplyr::select(Gene, 
                protein_length)

## read in protein abundance data
ho_proteins <- read_excel(path = "data/1-s2.0-S240547121730546X-mmc5.xlsx", sheet = 1, skip = 2) %>% 
  dplyr::rename(Gene = `Systematic Name`,
                mean_mol_cell = `Mean molecules per cell`,
                med_mol_cell = `Median molecules per cell`)

## match the Hackett et al data with the Yeast GEM data
hackett_fluxes_w_info <- hackett_fluxes %>% 
  inner_join(yeast_gem %>% 
               dplyr::rename(Model_Reaction_ID = ID),
             by = "Model_Reaction_ID")

get_rate_id_gene_name_mapping <- function(yeast_gem, 
                                          hackett_prot){
  ########
  # Function goal: loop through all the proteome and assign reaction IDs to each row
  ########
  supreme_df <- data.frame(reaction_id = character(),
                           gene_name = character(),
                           formula = character())
  
  # go through the dataframe of protein names to get the appropriate gene name
  # then find the associated reaction name
  for(i in 1:nrow(hackett_prot)){
    # i <- 4
    # print(i)
    protein_i <- hackett_prot[i, ]$Gene # extract the gene name
    row_containing_reaction <- dplyr::filter(yeast_gem, grepl(protein_i, `GENE ASSOCIATION`))# extract the reaction name
    
    # what happens if there is no corresponding reaction? append a "no-corresponding-reaction'
    if(nrow(row_containing_reaction) == 0){
      # print(c(protein_i, 'has no corresponding reaction id'))
      reaction_id_i <- "no-corresponding-reaction"
      # reactants_only_i <- "no-corresponding-reactants"
    }
    else {
      reaction_id_i <- row_containing_reaction$ID
      reaction_formula_i <- row_containing_reaction$EQUATION
      # reactants_only_i <- format_yeast_reactants(yeast_gem_formula = reaction_formula_i)
    }
    
    temp_df <- data.frame(reaction_id = reaction_id_i,
                          gene_name = protein_i,
                          formula = reaction_formula_i) # aggregate together for a row
    supreme_df <- rbind(supreme_df, 
                        temp_df) # append row to df
  }
  return(supreme_df)
}

get_rate_id_gene_name_mapping_gem_input <- function(yeast_gem){
  
  ########
  # Function goal: loop through all the proteome and assign reaction IDs to each row
  ########
  
  supreme_df <- data.frame(reaction_id = character(),
                           gene_name = character(),
                           formula = character(),
                           reactants = character())
  
  # unique_gene_associations <- unique(yeast_gem$`GENE ASSOCIATION`)
  ## remove the single NA value from this vector
  # unique_gene_associations_no_na <- unique_gene_associations[!is.na(unique_gene_associations)]
  # go through the dataframe of protein names to get the appropriate gene name
  # then find the associated reaction name
  for(i in 1:nrow(yeast_gem)){
    # i <- 3
    protein_i <- yeast_gem$`GENE ASSOCIATION`[i] # extract the gene name
    reaction_id_i <- yeast_gem$ID[i]
    reaction_reactants_i <- format_yeast_reactants_gem(yeast_gem$EQUATION[i])
    
    # row_containing_reaction <- dplyr::filter(yeast_gem, `GENE ASSOCIATION` == protein_i)# extract the reaction formula
    # reaction_reactants_i <- format_yeast_reactants_gem(row_containing_reaction$EQUATION)
    # reaction_name_i <- row_containing_reaction$ID
    
    ## checking if there are multiple proteins -- there will be an 'or' statement if that's true.
    if(grepl(x = protein_i, 
             pattern = "or", 
             fixed = TRUE)){
      
      ## separating all proteins into a vector of proteins
      full_sub_string_protein <- strsplit(protein_i, split = " or ") %>% 
        unlist()
      
      ### Making a vector of the reactants.
      reactant_vec <- rep(reaction_reactants_i, 
                          length(full_sub_string_protein))
      
      ### Making a vector of the reaction name.
      reaction_name_vec <- rep(reaction_id_i,
                               length(full_sub_string_protein))
      
      ### Putting these all together into a dataframe.
      df_sub_temp <- data.frame(reaction_name = reaction_name_vec, 
                                protein_name = full_sub_string_protein,
                                reactants = reactant_vec)
    }
    ## if its just one gene per reaction, then use the variables above this if statement.
    else {
      df_sub_temp <- data.frame(reaction_name = reaction_id_i, 
                                protein_name = protein_i,
                                reactants = reaction_reactants_i)
    }
    supreme_df <- rbind(supreme_df, 
                        df_sub_temp) # append row to df
  }
  return(supreme_df)
}

rate_to_id <- get_rate_id_gene_name_mapping(yeast_gem = yeast_gem, 
                                            hackett_prot = hackett_prot)

## Note that there are more rates than proteins -- some proteins mediate multiple rates in the Yeast GEM!
### which is why: nrow(rate_to_id) > nrow(hackett_prot)
### One example is r_0310 and r_0815 are both mediated by YAL012W.

## this function is to get the rate to protein mapping for all proteins, not just those observed in Hackett et al.
rate_to_id_from_gem <- get_rate_id_gene_name_mapping_gem_input(yeast_gem = yeast_gem)

### map the reaction rate ids to the protein relative abundances
hackett_fluxes_with_prots <- hackett_prot %>% ## first join the rate_to_id with the Hackett protein relative abundances
  inner_join(rate_to_id %>% 
               dplyr::rename(Gene = gene_name), by = "Gene") %>%
  inner_join(hackett_fluxes %>% 
               dplyr::rename(reaction_id = Model_Reaction_ID), 
             by = "reaction_id") # then join this with the hackett fluxes

## Note that the following dimensional changes occurred above. hackett_prot has 1187 proteins. Once that is joined with rate_to_id, 
## there is an increase in the number of rows to 2592 This is because there are multiple reactions per protein, and the inner_join
## function **duplicates** those rows. Finally, the hackett_fluxes df is joined, which filters out a lot of the pairs above, because there are only 233 reactions (hackett_fluxes$Model_Reaction_ID %>% unique() %>% length() == 233). However, because there are multiple reactions per protein, this joined dataframe becomes 284 observations. For example, hackett_fluxes_with_prots %>% dplyr::filter(Gene %in% "YPL134C") %>% as.data.frame() shows that there are identical rows for all the protein values, but the flux values are different.

# checking the fluxes that have more than one enzyme mapped
more_than_one_enzyme <- hackett_fluxes_with_prots$Gene %>% 
  table() %>% 
  as.data.frame()

more_than_one_rate <- hackett_fluxes_with_prots$reaction_id %>% 
  table() %>% 
  as.data.frame()

# checking the fluxes that have more than one enzyme mapped from the yeast gem (includes proteins that are not measured)
more_than_one_enzyme_gem <- rate_to_id_from_gem$protein_name %>% 
  table() %>% 
  as.data.frame()

more_than_one_rate_gem <- rate_to_id_from_gem$reaction_name %>% 
  table() %>% 
  as.data.frame()

## join the number of enzyme repeats in a df for manual examination
rate_to_id_from_gem_joined <- rate_to_id_from_gem %>% 
  left_join(dplyr::rename(.data = more_than_one_enzyme_gem, 
                          protein_name = .), 
                                   by = 'protein_name') %>% 
  left_join(dplyr::rename(.data = more_than_one_rate_gem,
                          reaction_name = .),
            by = 'reaction_name') %>% 
  mutate(frequencies = Freq.x + Freq.y)

## filter the genes that have one flux per protein, and one protein per flux
singles_only_rates_proteins <- rate_to_id_from_gem_joined %>% 
  # getting the rate/protein pairs that each have a frequency of 1, i.e., they sum to 2
  dplyr::filter(frequencies == 2)

# manual examination
# write.csv(singles_only_rates_proteins, 'tester5.csv')

## There are some Genes that have multiple fluxes per gene.
## here are the examples. (also described above). Identifying the proteins that have a single flux per protein.
singles_only_vec <- singles_only_rates_proteins$protein_name

reactions_w_single_genes_only <- hackett_fluxes_with_prots %>% 
  dplyr::select(reaction_id, 
                Gene) %>% 
  unique() %>% 
  dplyr::filter(Gene %in% singles_only_vec)

# Here is an example. transketolase, which mediates r_1050 and r_1049. There are two enzymes, and two reactions. They are TK1 and TK2. But the protein levels are identical for each -- they are not distinguished in the proteomic data but are in the flux data.
# http://bigg.ucsd.edu/models/iMM904/reactions/TKT2 http://bigg.ucsd.edu/models/iMM904/reactions/TKT1
# both are in the pentose phosphate pathway, but the assignment of each gene to each specific flux is uncertain.

# see:
hackett_fluxes_with_prots %>% 
  filter(reaction_id %in% c('r_1050', 'r_1049')) %>% 
  dplyr::select(P0.11_QP)

hackett_fluxes %>% 
  filter(Model_Reaction_ID %in% c('r_1050', 'r_1049'))

# Formatting data to feed into Stan model
col_names <- names(hackett_fluxes_with_prots)
col_names_proteins <- col_names[c(2:26)]
col_names_fluxes <- col_names[c(30:54)]

# Converting the data into a long dataframe not a wide dataframe:
melted_proteins_conditions <- melt(data = hackett_fluxes_with_prots %>% 
                                     dplyr::select(all_of(c("Gene", "reaction_id", col_names_proteins))), 
                                   variable.name = 'condition-prot', 
                                   value.name = 'protein_amount')
melted_fluxes_conditions <- melt(data = hackett_fluxes_with_prots %>% 
                                   dplyr::select(all_of(c("Gene", "reaction_id", col_names_fluxes))), 
                                 variable.name = 'condition-flux', 
                                 value.name = 'flux')

# This makes a tidy dataframe of the fluxes, conditions, protein names, and reaction ids
tidy_df_flux_prot <- melted_proteins_conditions %>% 
  # one of the letters is not capitalized in the phosphate limited growth condition
  mutate(`condition-prot` = fct_recode(`condition-prot`, "P0.22" = "p0.22")) %>% 
  dplyr::rename(condition_protein = `condition-prot`) %>% 
  inner_join(melted_fluxes_conditions %>% 
               dplyr::rename(condition_protein_qp = `condition-flux`) %>% 
               dplyr::mutate(condition_protein = gsub(pattern = '_QP', 
                                                      replacement = '', 
                                                      x = condition_protein_qp)), 
             by = c('Gene', 'reaction_id', 'condition_protein')) %>% 
  mutate(log_flux = log2(flux)) %>% 
  # remove the condition with the '_QP' appended at the end
  dplyr::select(-condition_protein_qp)

# Getting mean protein level for each protein, then appending that to the entire dataframe.
mean_protein_amount_df <- tidy_df_flux_prot %>% 
  group_by(Gene) %>% 
  mutate(unlog_protein_amount = 2^protein_amount) %>% 
  ## filter out the negative fluxes, those observations cannot be log transformed
  # dplyr::filter(flux > 0) %>% 
  ## Filter out the conditions with knockouts, specifically U0 and L0
  dplyr::filter(!grepl(pattern = 'U0|L0', condition_protein)) %>% 
  summarize(mean_protein_amount = mean(protein_amount, na.rm = TRUE),
            mean_log_flux = mean(log_flux, na.rm = TRUE),
            mean_unlog_flux = mean(flux, na.rm = TRUE),
            sd_log_flux = sd(log_flux, na.rm = TRUE),
            sd_unlog_flux = sd(flux, na.rm = TRUE),
            sd_protein_amount = sd(protein_amount, na.rm = TRUE),
            mean_protein_amount_unlog = mean(unlog_protein_amount, na.rm = TRUE)) %>% 
  inner_join(yeast_genes_protein_lengths, by = 'Gene') %>% 
  rowwise() %>% 
  mutate(log_weighted_protein_mean = log2(mean_protein_amount_unlog*protein_length))

# Now there are 183 unique proteins, which is also the dimension of the dataframe.
# Note that this has decreased the number of unique proteins in tidy_df_flux_prot from 231, then to 222 after filtering
# out negative fluxes. The dataframe yeast_genes_protein_lengths is what brings it down to 183. This is okay, because we then filter out for only those proteins that have unique flux-to-protein mapping.
mean_protein_amount_df %>% dim()
mean_protein_amount_df$Gene %>% unique() %>% length()

mean_protein_amount_singles <- mean_protein_amount_df %>% 
  dplyr::filter(Gene %in% singles_only_vec)

# remove the ambiguous gene to flux mapping genes
tidy_df_flux_prot_norm_singles_only <- tidy_df_flux_prot %>% 
  dplyr::filter(Gene %in% singles_only_vec)


## there are 52 unique proteins in 'tidy_df_flux_prot_norm_singles_only':
## tidy_df_flux_prot_norm_singles_only$Gene %>% unique() %>% length()

### appending the summary statistic
joined_df_flux_prot_mean <- tidy_df_flux_prot_norm_singles_only %>% 
  ## Remove the conditions with the mutatnts
  dplyr::filter(!grepl(pattern = 'U0|L0', 
                       condition_protein)) %>%
  inner_join(mean_protein_amount_singles, 
             by = 'Gene')

## there are 47 unique proteins in 'joined_df_flux_prot_mean':
## joined_df_flux_prot_mean$Gene %>% unique() %>% length()

## The reason for the jump from 81 to 47 is because some of those proteins only mediate negative flux values.

# Getting substrate connectivity metric -----------------------------------
yeast_metabolites_vec <- all_yeast_metabolites$ID

get_unique_prots_per_metabolite_yeast <- function(yeast_metabolites_vec, 
                                                  rate_to_id_from_gem){
  
  ### make an empty vector for unique proteins per metabolite
  unique_proteins_per_metabolite_vec_out <- c()
  ### loop through every metabolite
  for(i in 1:length(yeast_metabolites_vec)){
    # i <- 2757
    ### For every metabolite, filter the above dataframe using grepl and then count the number of unique proteins (or protein complexes) per metabolite
    metabolite_i <- paste0(" ", 
                           yeast_metabolites_vec[i])
    metabolite_i_containing_df_sub <- rate_to_id_from_gem %>% 
      dplyr::filter(grepl(pattern = metabolite_i, 
                          x = reactants, 
                          fixed = TRUE)) ## You need this fixed statement
    if(nrow(metabolite_i_containing_df_sub) == 0){
      unique_proteins_per_metabolite <- 0
    } else {
      all_proteins_subset <- unique(metabolite_i_containing_df_sub$protein_name)
      all_proteins_subset_no_na <- all_proteins_subset[!is.na(all_proteins_subset)]
      unique_proteins_per_metabolite <- length(all_proteins_subset_no_na)
      ## sometimes there is only an NA, in which case assign the value to be zero
      if(length(unique_proteins_per_metabolite) == 0){
        unique_proteins_per_metabolite <- 0
      }
    }
    unique_proteins_per_metabolite_vec_out <- c(unique_proteins_per_metabolite_vec_out, 
                                                unique_proteins_per_metabolite)
    
  }
  return(unique_proteins_per_metabolite_vec_out)
}

unique_proteins_per_metabolite_vec_yeast <- get_unique_prots_per_metabolite_yeast(yeast_metabolites_vec = yeast_metabolites_vec, 
                                                                                  rate_to_id_from_gem = rate_to_id_from_gem)

### arrive at a metabolite connectivity score for each protein. 
metabolite_connect_df_yeast <- data.frame(proteins_per_metabolite = unique_proteins_per_metabolite_vec_yeast,
                                          metabolite_name = yeast_metabolites_vec) %>% 
  # remove water, assuming that water does not limit reaction rates.
  dplyr::filter(!grepl('H2O', metabolite_name, fixed = TRUE))

### manual check of this dataframe.
# write.csv(metabolite_connect_df_yeast, 
#           file = 'metabolite_connect_df_yeast.csv')

multi_substrate_averaging_yeast <- function(rate_to_id_from_gem, 
                                            metabolite_connect_df_yeast){
  ### Function for getting a metric of substrate connectivity when enzymes have multiple substrate
  
  ### this is a vector as long as the vector of unique proteins in 'rate_to_id_from_gem'
  protein_score_vector_arithmetic <- c()
  protein_score_vector_geometric <- c()
  all_unique_proteins <- unique(rate_to_id_from_gem$protein_name)
  
  ### loop through every protein
  for(i in 1:length(all_unique_proteins)){
    # i <- 1
    unique_protein_i <- all_unique_proteins[i]
    
    # subset the rate_to_id_from_gem dataframe for reactions that contain a given protein.
    subset_reactant_df_i <- rate_to_id_from_gem %>% 
      dplyr::filter(grepl(unique_protein_i, 
                          protein_name))
    
    ## get all metabolite reactants associated with protein i. must be pasted if the protein is using multiple reactants.
    all_metabolite_reactants_protein_i <- paste(subset_reactant_df_i$reactants, 
                                                collapse = " ")
    
    ## vector of all the substrate scores
    protein_score_of_all_metabolites <- c()
    
    ### loop through every metabolite from above df
    for(k in 1:nrow(metabolite_connect_df_yeast)){
      ### need to add a space to the front, because we don't want metabolite names that are substrings of other metabolites.
      metabolite_k <- paste0(" ", 
                             metabolite_connect_df_yeast$metabolite_name[k])
      if(grepl(pattern = metabolite_k, 
               x = all_metabolite_reactants_protein_i, 
               fixed = TRUE)){
        protein_score_of_all_metabolites <- c(protein_score_of_all_metabolites, 
                                              metabolite_connect_df_yeast$proteins_per_metabolite[k])
      }
    }
    if(is.null(protein_score_of_all_metabolites)){
      protein_score_of_all_metabolites <- 0
    }
    protein_score_vector_arithmetic <- c(protein_score_vector_arithmetic, 
                                         mean(protein_score_of_all_metabolites, 
                                              na.rm = TRUE))
    protein_score_vector_geometric <- c(protein_score_vector_geometric,
                                        exp(mean(log(protein_score_of_all_metabolites))))
  }
  
  df_protein_level_substrate_connectivity <- data.frame(protein_name = all_unique_proteins,
                                                        arithmetic_mean = protein_score_vector_arithmetic,
                                                        geometric_mean = protein_score_vector_geometric)
  
  return(df_protein_level_substrate_connectivity)
}

df_protein_level_substrate_connectivity_yeast <- multi_substrate_averaging_yeast(rate_to_id_from_gem = rate_to_id_from_gem, 
                                                                                 metabolite_connect_df_yeast = metabolite_connect_df_yeast)

# Filtering out negative fluxes -------------------------------------------

joined_df_flux_prot_mean_no_negative <- joined_df_flux_prot_mean %>% 
  # dplyr::filter(flux > 0) %>% ## two values are removed. ### Changed this
  rowwise() %>% 
  ## joining with data on protein length
  inner_join(yeast_genes_protein_lengths, by = 'Gene') %>% 
  mutate(z_score_unlog_flux = (flux - mean_unlog_flux)/sd_unlog_flux,
         z_score_flux = (log_flux - mean_log_flux)/sd_log_flux,
         z_score_prot = (protein_amount - mean_protein_amount)/sd_protein_amount) %>% 
  rowwise() %>% 
  mutate(gene_as_factor = as.factor(Gene)) %>% 
  ## This dataframe also has YHR128W, which has zero flux over all conditions. It also has YDR300C, which has a paralog,
  ## and one of our filtering conditions is to remove all proteins that don't have a single unambiguous mapping to a rate.
  dplyr::filter(Gene != 'YHR128W', Gene != 'YDR300C') %>% 
  droplevels()

# Making another column that has genes as a numeric factor
joined_df_flux_prot_mean_no_negative$gene_as_factor_number <- joined_df_flux_prot_mean_no_negative$gene_as_factor %>% as.numeric()

# Make a dataframe that has the gene name to gene as factor number mapping
joined_df_flux_prot_mean_no_negative_gene_as_number <- joined_df_flux_prot_mean_no_negative %>% 
  dplyr::select(Gene, gene_as_factor_number) %>% 
  unique()

# Which genes contain a negative flux value?
joined_df_flux_prot_mean_only_negative <- joined_df_flux_prot_mean %>% 
  dplyr::filter(flux < 0)

joined_df_flux_prot_mean_no_negative %>% 
  dplyr::filter(flux < 0)

# appending reaction thermodynamics to the summary dataframe --------------

high_confidence_fluxes_hackett <- hackett_fluxes$Model_Reaction_ID

# Make a dataframe that has both the rate ID's and the Delta G values.
rate_to_id_w_delta_g <- rate_to_id %>% 
  dplyr::filter(gene_name %in% singles_only_vec,
                reaction_id %in% high_confidence_fluxes_hackett) %>%
  inner_join(delta_g %>% 
               dplyr::rename(reaction_id = Var1,
                             delta_g = Var2), 
             by = "reaction_id")

# Filter the very high 1e7 delta_g to manually explore which reactions they are.
missing_delta_g_yeast <- rate_to_id_w_delta_g %>% 
  ### These should be NA values -- no delta G assigned. 
  dplyr::filter(delta_g == 1e7) %>% 
  inner_join(yeast_gem %>% 
               dplyr::rename(reaction_id = ID))

# It's only r_0361
# missing_delta_g_yeast$reaction_id %>% unique()

# Making an explanatory variable dataframe --------------------------------

# Taking the mean protein dataframe
mean_protein_amount_singles_w_gene_factor <- mean_protein_amount_singles %>% 
  inner_join(joined_df_flux_prot_mean_no_negative_gene_as_number, 
             by = 'Gene') %>% 
  # Join this dataframe with the protein abundance from Ho et al 2018, Cell Systems
  inner_join(ho_proteins %>% dplyr::select(Gene, med_mol_cell), by = "Gene") %>% 
  mutate(aa_per_cell = log2(med_mol_cell*protein_length))

## Need to read in this script to calculate the mean protein abundance.
source('scripts/hackett-mean-protein-abundances.R')

# Need to join the mean protein amount dataframe with the substrate connectivity metric, deltaG, and also need to weight the protein amount
mean_protein_amount_singles_summaries <- mean_protein_amount_singles_w_gene_factor %>% 
  inner_join(df_protein_level_substrate_connectivity_yeast %>% 
               rename(Gene = protein_name), 
             by = 'Gene') %>% 
  inner_join(rate_to_id_w_delta_g %>% 
               rename(Gene = gene_name), 
             by = 'Gene') %>%
  inner_join(summarized_abundances_across_conditions %>% 
               dplyr::rename(Gene = protein_name) %>% 
               dplyr::select(Gene, log_mean_prot_mf_perc)) %>% 
  dplyr::filter(delta_g != 10000000.00)

# After making this summary dataframe, I need to make a long-form dataframe that only has this subset of proteins
summary_set_of_proteins <- mean_protein_amount_singles_summaries$Gene
joined_df_flux_prot_mean_no_negative_subset <- joined_df_flux_prot_mean_no_negative %>% 
  dplyr::filter(Gene %in% summary_set_of_proteins)

# Need to make a new gene as factor index, because the filtered one was level 45, and we need the number of total indices to be 46
joined_df_flux_prot_mean_no_negative_subset$gene_as_factor_number_summary <- as.numeric(as.factor(joined_df_flux_prot_mean_no_negative_subset$Gene))

