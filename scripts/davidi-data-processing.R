######
# This script formats data from E. coli protein abundances and corresponding fluxes.
######

### Read in relevant libraries.
library(readxl)
library(magrittr)
library(dplyr)
library(reshape2)
library(randomForestSRC)
library(ggplot2)
library(forcats)
library(stringr)
library(rstan)
library(tidybayes)
library(infotheo)

'%!in%' <- function(x,y)!('%in%'(x,y))

# Read in data ------------------------------------------------------------

## E. coli Genome Scale Model has gene name to reaction mappings.
e_coli_gene_id_to_function <- read_excel(path = 'data/msb201165-sup-0002.xls',
                                         sheet = 8) %>% 
  dplyr::rename(bnumber = `Gene ID`)

## Reading in the metabolite names
all_e_coli_metabolites <- read_excel(path = "data/msb201165-sup-0002.xls", 
                                     sheet = 4)

### data for mapping gene name to reaction names for fluxes. Note that these are from the 
### genome mapped for MG1655, not the strain with proteomic data.
fluxes_to_prots <- read_excel(path = 'data/msb201165-sup-0002.xls', sheet = 3)

### data for fluxes and protein amounts aggregated from Davidi et al PNAS
fluxes_dav <- read_excel(path = "data/pnas.1514240113.sd01.xlsx", 
                         sheet = 7, 
                         skip = 1)
prots_dav <- read_excel(path = 'data/pnas.1514240113.sd01.xlsx', 
                        sheet = 8, 
                        skip = 1)

### converting the prots_dav df into a mapping dataframe
e_coli_b_number_gene_name_mappings <- prots_dav %>% 
  dplyr::select(bnumber, 
                `gene name`) %>% 
  rename(gene_name_lett = `gene name`,
         gene_name = bnumber) %>% 
  unique()

### data for reaction thermodynamics
rxn_delta_g <- read_excel(path = 'data/pcbi.1006492.s005.xlsx', 
                          sheet = 2) %>% 
  dplyr::rename(delta_g = `Delta G0`)

# formatting reaction name
rxn_delta_g$reaction_name <- substr(rxn_delta_g$Reaction, 
                                    3, 
                                    nchar(rxn_delta_g$Reaction))

# Process data to make a matrix of fluxes to protein amounts --------------
format_e_coli_reactants <- function(e_coli_gem_df){
  ### Function goal: output a vector of reactants for each reaction.
  ### For the e coli data, first split the reaction *Formula* column into two
  e_coli_gem_formulas <- e_coli_gem_df$Formula
  ### Include both reactants and products if the reaction is reversible
  if(grepl("<=>", e_coli_gem_formulas)){
    left_side_only_e_coli <- e_coli_gem_formulas
  } else {
    ### Include only reactants if the reaction is irreversible
    right_and_left_sides_e_coli <- str_split(string = e_coli_gem_formulas, 
                                             pattern = ">")
    ## Then subset only the left side.
    left_side_only_e_coli <- sapply(right_and_left_sides_e_coli, 
                                    "[[", 
                                    1)
  }
  ## Including a padded zero on the left side so that reactants can be searched using strings, and won't 
  ## get a sub-string. For example, 'moco[c]' would be counted as well as 'bmoco[c]'.
  return(paste0(" ", 
                left_side_only_e_coli))
}

## function for looping through all fluxes and then getting paired gene names
get_df_matching_genes_rxns <- function(input_df){
  
  ## This function returns a dataframe of matching genes for every flux. Multiple genes per flux are 
  ## coded as multiple rows in the dataframe. 
  ### Last code audit: Dec 15, 2023
  
  df_supreme <- data.frame(reaction_name = character(),
                           gene_name = character(),
                           protein_name = character(),
                           formula_left_side = character())
  
  # loop through all reactions in the spreadsheet
  for(i in 1:nrow(input_df)){
    ## get the reaction abbreviation
    reaction_abbr_i <- input_df$`Reaction Abbreviation`[i]
    ## get the protein name
    protein_i <- input_df$`Protein-Reaction Association`[i]
    ## get the gene name
    gene_i <- input_df$`Gene-Reaction Association`[i]
    # get the formula left hand side for the reaction for row i
    left_side_i <- format_e_coli_reactants(e_coli_gem_df = input_df[i, ])
    ## Check to see if there are multiple proteins associated with a single flux, that can be mediated by **separate** reactions
    if(grepl(x = protein_i, 
             pattern = "or", 
             fixed = TRUE)){
      ## separating all the genes into individual instances
      full_sub_string_protein <- strsplit(protein_i, split = " or ") %>% 
        unlist()
      # removing the brackets from all strings
      sub_string_no_brackets <- gsub(pattern = "(", 
                                     replacement = "", 
                                     x = full_sub_string_protein,
                                     fixed = TRUE)
      sub_string_no_brackets_2 <- gsub(pattern = ")", 
                                       replacement = "", 
                                       x = sub_string_no_brackets,
                                       fixed = TRUE)
      removed_white_space_sub_string <- str_trim(string = sub_string_no_brackets_2)
      
      ## Formatting the gene names
      full_sub_string_gene <- strsplit(gene_i, split = " or ") %>% 
        unlist()
      gene_sub_string_no_brackets <- gsub(pattern = "(", 
                                          replacement = "", 
                                          x = full_sub_string_gene,
                                          fixed = TRUE)
      gene_sub_string_no_brackets2 <- gsub(pattern = ")", 
                                           replacement = "", 
                                           x = gene_sub_string_no_brackets,
                                           fixed = TRUE)
      
      ## Making a vector of the reaction names and proteins.
      flux_vec_sub <- rep(reaction_abbr_i, 
                          length(removed_white_space_sub_string))
      reactant_vec <- rep(left_side_i, 
                          length(removed_white_space_sub_string))
      
      df_sub_temp <- data.frame(reaction_name = flux_vec_sub, 
                                gene_name = gene_sub_string_no_brackets2,
                                protein_name = removed_white_space_sub_string,
                                formula_left_side = reactant_vec)
    }
    ## if its just one gene per reaction, then use the variables above this if statement.
    else {
      df_sub_temp <- data.frame(reaction_name = reaction_abbr_i, 
                                gene_name = gene_i,
                                protein_name = protein_i,
                                formula_left_side = left_side_i)
    }
    
    df_supreme <- rbind(df_supreme, 
                        df_sub_temp)
  }
  return(df_supreme)
}

genes_matching_rxns <- get_df_matching_genes_rxns(input_df = fluxes_to_prots)

##### calculate substrate connectivity metric per gene name

# Data formatting for getting single reaction to protein mappings ---------

### Now get a dataframe that has only the reactions with a single, unambiguous mapping of gene to rxn

## First calculating the number of genes per reaction
genes_per_rxn <- genes_matching_rxns %>% 
  group_by(reaction_name) %>% 
  summarize(genes_per_rxn_number = n())

## Calculate the number of reactions per protein
rxns_per_protein <- genes_matching_rxns %>% 
  group_by(protein_name) %>% 
  summarize(rxns_per_protein_number = n())

# what does the distribution of genes per rxn look like
# genes_per_rxn$genes_per_rxn_number %>% hist()
# genes_per_rxn$genes_per_rxn_number %>% quantile()

## subset the rxns with only one gene
rxns_with_only_one_gene <- genes_per_rxn %>% 
  dplyr::filter(genes_per_rxn_number == 1)

## get this vector of rxn names
rxns_with_only_one_gene_vec <- rxns_with_only_one_gene$reaction_name

## subset the proteins with only one reaction
genes_with_only_one_rxn <- rxns_per_protein %>% 
  dplyr::filter(rxns_per_protein_number == 1)
genes_with_only_one_rxn_vec <- genes_with_only_one_rxn$protein_name

# get those genes with only one rxn
genes_matching_rxns_one_rxn_genes <- genes_matching_rxns %>% 
  dplyr::filter(reaction_name %in% rxns_with_only_one_gene_vec) %>% 
  dplyr::filter(protein_name %in% genes_with_only_one_rxn_vec)

### Now we have the unambiguously mapped proteins to reactions, we want to join the flux matrix and then the 
### protein matrix to this.
fluxes_genes_gene_names <- genes_matching_rxns_one_rxn_genes %>% 
  ### join the flux dataframe that also has paired reaction names
  inner_join(fluxes_dav %>% 
               dplyr::rename(reaction_name = `reaction (model name)`), 
             by = "reaction_name") %>%
  ## convert it to a long format
  melt(variable.name = 'exp_condition', value.name = 'flux')

### there should be equal number of unique reaction names and unique gene names
fluxes_genes_gene_names$reaction_name %>% unique() %>% length()
fluxes_genes_gene_names$gene_name %>% unique() %>% length()

# converting the protein df to a long df not a wide df
melted_prots_dav <- prots_dav %>% 
  dplyr::select(-`gene name`) %>% 
  dplyr::rename(gene_name = bnumber) %>% 
  melt(variable.name = 'exp_condition', value.name = 'protein_amount') %>% 
  group_by(gene_name, exp_condition) %>% 
  ### this summarize statement is weirdly necessary. In the original data sheet there are duplicate rows
  ### of protein levels, e.g. b0221 (fadE) is present 7 times in the spreadsheet. All those values are identical though.
  ### so summarizing it using the mean won't change the replicate value. To verify this, calculate the sd as an additional 
  ### variable in the statement below, and the values are all 0 or NA.
  summarize(protein_amount_mean = mean(protein_amount, na.rm = TRUE))

# Joining the melted (long format) fluxes and protein levels together.
genes_fluxes_long_df <- fluxes_genes_gene_names %>% 
  inner_join(melted_prots_dav, 
             by = c('gene_name', 'exp_condition'))
# This ends up with a dataframe that is 3875 rows. This is because there are genes 
# in the 'melted_prots_dav' that are not uniquely mapping to fluxes.

## These lengths should all be identical (125 different protein-flux pairs)
genes_fluxes_long_df$reaction_name %>% unique() %>% length()
genes_fluxes_long_df$gene_name %>% unique() %>% length()
genes_fluxes_long_df$protein_name %>% unique() %>% length()

## goal is to filter out values where there are multiple fluxes per gene, or multiple genes per flux.
## NOTE: this part is redundant, above a similar thing was done, so there are still 125 flux-protein pairs.
genes_with_multiple_rxns <- genes_fluxes_long_df %>%
  dplyr::select(gene_name, reaction_name) %>%
  unique() %>%
  dplyr::select(gene_name) %>%
  table()

genes_with_multiple_rxns_df <- genes_with_multiple_rxns %>% as.data.frame()
genes_with_multiple_rxns_df_1 <- genes_with_multiple_rxns_df %>%
  filter(Freq == 1)
genes_with_single_rxns_only <- genes_with_multiple_rxns_df_1$gene_name

### e_coli_header_names_df comes from ecoli_fasta_formatting.R
source('scripts/ecoli_fasta_formatting.R')

e_coli_protein_length_gene_names <- e_coli_header_names_df %>% 
  rename(gene_name_lett = gene_name) %>% 
  ## joining with the bnumbers from Davidi et al mapped to gene names
  inner_join(e_coli_b_number_gene_name_mappings) %>% 
  dplyr::select(gene_name_lett, 
                protein_length, 
                gene_name) %>% 
  group_by(gene_name_lett, gene_name) %>% 
  summarize(protein_length_mean = mean(protein_length))
# This mean aggregates only a few protein lengths (eventually for use in the 2-tier hierarchical model), for example bioD1. This is because 
# it is not exactly clear which protein in the E. coli strain this protein is referring to.

# Transforming the fluxes and prots into z scores. First I need to make a summary dataframe that has
# the mean and standard deviations of fluxes and protein levels.
# This results in 125 reactions.
genes_fluxes_long_df_summary_stats <- genes_fluxes_long_df %>% 
  ## This line doesn't matter, it is just to ensure that only proteins with single 
  ## reactions are included. But the df is identical either way.
  dplyr::filter(gene_name %in% genes_with_single_rxns_only) %>% 
  inner_join(e_coli_protein_length_gene_names, 
             by = 'gene_name') %>% 
  group_by(gene_name, 
           reaction_name) %>% 
  ### these conditions are removed because the protein amounts to not correspond with the 
  ### original publications.
  dplyr::filter(exp_condition %!in% c('GAL_BATCH_mu=0.26_S', 'GAM_BATCH_mu=0.46_S', 
                                      'GLC+SALTS_BATCH_mu=0.55_S', 'XYL_BATCH_mu=0.55_S',
                                      'GLC_CHEM_mu=0.26_P', 'GLC_CHEM_mu=0.46_P')) %>% 
  summarize(log_flux_mean = mean(log2(flux), na.rm = TRUE),
            flux_mean = mean(flux, na.rm = TRUE),
            log_prot_amount_mean = mean(log2(protein_amount_mean), na.rm = TRUE),
            log_prot_amount_mean_weighted = log(mean(protein_amount_mean*protein_length_mean, na.rm = TRUE)),
            tester = mean(log2(protein_amount_mean*protein_length_mean), na.rm = TRUE),
            log_flux_sd = sd(log2(flux), na.rm = TRUE),
            flux_sd = sd(flux, na.rm = TRUE),
            log_prot_amount_sd = sd(log2(protein_amount_mean), na.rm = TRUE)) 

# Calculating substrate connectivity metric -------------------------------

e_coli_metabolites_vec <- all_e_coli_metabolites$`Metabolite Abbreviation`

get_unique_prots_per_metabolite <- function(e_coli_metabolites_vec, 
                                            genes_matching_rxns){
  ### make an empty vector for unique proteins per metabolite
  unique_proteins_per_metabolite_vec_out <- c()
  ### loop through every metabolite
  for(i in 1:length(e_coli_metabolites_vec)){
    ### For every metabolite, filter the above dataframe using grepl and then count the number of unique proteins (or protein complexes) per metabolite
    metabolite_i <- paste0(" ", e_coli_metabolites_vec[i])
    metabolite_i_containing_df_sub <- genes_matching_rxns %>% 
      dplyr::filter(grepl(pattern = metabolite_i, 
                          x = formula_left_side, 
                          fixed = TRUE)) ## You need this fixed statement
    if(nrow(metabolite_i_containing_df_sub) == 0){
      unique_proteins_per_metabolite <- 0
    } else {
      all_proteins_subset <- unique(metabolite_i_containing_df_sub$protein_name)
      unique_proteins_per_metabolite <- length(all_proteins_subset)
    }
    unique_proteins_per_metabolite_vec_out <- c(unique_proteins_per_metabolite_vec_out, 
                                                unique_proteins_per_metabolite)
    
  }
  return(unique_proteins_per_metabolite_vec_out)
}

unique_proteins_per_metabolite_vec <- get_unique_prots_per_metabolite(e_coli_metabolites_vec = e_coli_metabolites_vec, 
                                                                      genes_matching_rxns = genes_matching_rxns)

## arrive at a metabolite connectivity score for each protein. 
metabolite_connect_df <- data.frame(proteins_per_metabolite = unique_proteins_per_metabolite_vec,
                                    metabolite_name = e_coli_metabolites_vec) %>% 
  # remove water, assuming that water does not limit reaction rates.
  dplyr::filter(!grepl('h2o', metabolite_name, fixed = TRUE))

## manual check of this dataframe.
write.csv(metabolite_connect_df, 
          file = 'data/metabolite_connect_df_e_coli.csv')

multi_substrate_averaging <- function(genes_matching_rxns, metabolite_connect_df){
  ### Function for getting a metric of substrate connectivity when enzymes have multiple substrate
  
  ### this is a vector as long as the vector of unique proteins in 'genes_matching_rxns'
  protein_score_vector_arithmetic <- c()
  protein_score_vector_geometric <- c()
  all_unique_proteins <- unique(genes_matching_rxns$protein_name)
  
  ### loop through every protein
  for(i in 1:length(unique(genes_matching_rxns$protein_name))){
    # i <- 1
    unique_protein_i <- all_unique_proteins[i]
    # subset the genes_matching_rxns dataframe for reactions that contain a given protein.
    subset_reactant_df_i <- genes_matching_rxns %>% 
      dplyr::filter(grepl(unique_protein_i, 
                          protein_name))
    ## get all metabolite reactants associated with protein i. must be pasted if the protein is using multiple reactants.
    all_metabolite_reactants_protein_i <- paste(subset_reactant_df_i$formula_left_side, 
                                                collapse = " ")
    ## set the score to 0 at first
    protein_score_of_all_metabolites <- c()
    
    ### loop through every metabolite from above df
    for(k in 1:nrow(metabolite_connect_df)){
      ### need to add a space to the front
      metabolite_k <- paste0(" ", 
                             metabolite_connect_df$metabolite_name[k])
      if(grepl(pattern = metabolite_k, 
               x = all_metabolite_reactants_protein_i, 
               fixed = TRUE)){
        protein_score_of_all_metabolites <- c(protein_score_of_all_metabolites, metabolite_connect_df$proteins_per_metabolite[k])
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
  
  df_protein_level_substrate_connectivity <- data.frame(protein_name = unique(genes_matching_rxns$protein_name),
                                                        arithmetic_mean = protein_score_vector_arithmetic,
                                                        geometric_mean = protein_score_vector_geometric) %>% 
    inner_join(genes_matching_rxns, by = 'protein_name')
  
  return(df_protein_level_substrate_connectivity)
}

df_protein_level_substrate_connectivity <- multi_substrate_averaging(genes_matching_rxns = genes_matching_rxns,
                                                                     metabolite_connect_df = metabolite_connect_df)

## spot checking
write.csv(x = df_protein_level_substrate_connectivity, 
          file = "data/df_protein_level_substrate_connectivity.csv")

# Formatting delta G values for reaction thermodynamics -------------------

## Need to add the delta G of reaction to this summary dataframe.
setdiff(genes_fluxes_long_df_summary_stats$reaction_name, rxn_delta_g$reaction_name) %>% length()

# This shows that there at 61 reactions missing from the rxn_delta_g dataframe read in above.

## Which ones are missing?
delta_g_missing <- genes_fluxes_long_df_summary_stats$reaction_name[genes_fluxes_long_df_summary_stats$reaction_name %!in% rxn_delta_g$reaction_name]

## Get a dataframe that has interpretable names of the missing delta G values so they can be 
## manually examined.
fluxes_to_prots_no_delta_g <- fluxes_to_prots %>% 
  dplyr::filter(grepl(pattern = paste(delta_g_missing, collapse='|'), 
                      `Reaction Abbreviation`))

# Write this as a CSV.
write.csv(x = fluxes_to_prots_no_delta_g, 
          file = 'data/fluxes_to_prots_no_delta_g.csv', 
          row.names = FALSE)

## some of the delta G values are not represented in the spreadsheet rxn_delta_g
unique(genes_fluxes_long_df_summary_stats$reaction_name)
write.csv(data.frame(gene_name = unique(genes_fluxes_long_df_summary_stats$gene_name)), 
          file = "data/e_coli_gene_names_in_set.csv", row.names = FALSE)


# Merging the summary statistics into one dataframe -------------------

genes_fluxes_summary_states_merged <- genes_fluxes_long_df_summary_stats %>% 
  left_join(rxn_delta_g, ## some delta G values are not present, so left join enables a NA to be added
            by = 'reaction_name') %>% 
  ## merging the substrate connectivity metric
  left_join(df_protein_level_substrate_connectivity, 
            by = c('reaction_name', 'gene_name'))

## looking at this distribution
genes_fluxes_summary_states_merged$geometric_mean %>% quantile()

### getting the vector of gene names in the second level df
gene_names_in_merged_df <- genes_fluxes_summary_states_merged$gene_name %>% unique()
reaction_names_in_merged_df <- genes_fluxes_summary_states_merged$reaction_name %>% unique()

## merge these summary statistics for calculating the transformed values
genes_fluxes_long_df_z_score <- genes_fluxes_long_df %>% 
  ### these conditions are removed because unsure if the corresponding fluxes match.
  dplyr::filter(exp_condition %!in% c('GAL_BATCH_mu=0.26_S', 'GAM_BATCH_mu=0.46_S',
                                      'GLC+SALTS_BATCH_mu=0.55_S', 'XYL_BATCH_mu=0.55_S',
                                      'GLC_CHEM_mu=0.26_P', 'GLC_CHEM_mu=0.46_P')) %>%
  left_join(genes_fluxes_summary_states_merged, by = c('gene_name', 'reaction_name')) %>% 
  rowwise() %>% 
  mutate(z_score_prot = (log2(protein_amount_mean) - log_prot_amount_mean)/log_prot_amount_sd,
         z_score_flux = (flux - flux_mean)/flux_sd)

## get only the rows that have no NA values, only considering the flux and protein levels. This removes specific cases
## where there are no protein levels, for example.
genes_fluxes_long_df_z_score_cc <- genes_fluxes_long_df_z_score[complete.cases(genes_fluxes_long_df_z_score %>% 
                                                                                 dplyr::select(gene_name,
                                                                                               z_score_flux,
                                                                                               z_score_prot,
                                                                                               exp_condition)), ]

## do complete cases on the entire df, which includes the explanatory variables (so this would exclude rows that have all
# required fluxes and protein levels, but are also missing for example, delta G values)
genes_fluxes_long_df_z_score_cc_2_level <- genes_fluxes_long_df_z_score[complete.cases(genes_fluxes_long_df_z_score), ]
genes_fluxes_summary_states_merged_cc <- genes_fluxes_summary_states_merged[complete.cases(genes_fluxes_summary_states_merged),]

# Making the gene label as a number
genes_fluxes_long_df_z_score_cc$gene_as_number <- genes_fluxes_long_df_z_score_cc$gene_name %>% 
  as.factor() %>% 
  as.numeric()
genes_fluxes_long_df_z_score_cc_gene_as_number_summary <- genes_fluxes_long_df_z_score_cc %>% 
  dplyr::select(gene_name, gene_as_number) %>% 
  unique()
genes_fluxes_summary_states_merged$gene_as_number <- genes_fluxes_summary_states_merged$gene_name %>% 
  as.factor() %>% 
  as.numeric()
genes_fluxes_long_df_z_score_cc_2_level$gene_as_number <- genes_fluxes_long_df_z_score_cc_2_level$gene_name %>% 
  as.factor() %>% 
  as.numeric()
genes_fluxes_summary_states_merged_cc$gene_as_number <- genes_fluxes_summary_states_merged_cc$gene_name %>% 
  as.factor() %>% 
  as.numeric()

## merging the summary dataframe to add formulas to it
genes_fluxes_summary_states_merged_w_formulas <- genes_fluxes_summary_states_merged %>% 
  left_join(fluxes_to_prots %>% 
              dplyr::rename(reaction_name = `Reaction Abbreviation`),
            by = 'reaction_name')

## converting function for metabolite abbreviation into metabolite full name for searching deltaG values
get_formula_full_names <- function(gene_summary_df, all_e_coli_metabolites){
  #### This function parses the reaction formula and converts it to full names
  formula_vector_to_parse <- paste0(" ", gene_summary_df$Formula)
  ## go through each reaction formula
  for(i in 1:length(formula_vector_to_parse)){
    # i <- 1
    ## for each reaction formula, go through every metabolite abbreviation.
    formula_i_from_vector <- formula_vector_to_parse[i]
    for(k in 1:nrow(all_e_coli_metabolites)){
      # k <- 307
      ## need to add a space so you don't get substring matches
      ## if that metabolite abbreviation is present, replace the text with the full metabolite name
      if(grepl(pattern = paste0(" ", all_e_coli_metabolites$`Metabolite Abbreviation`[k]), 
               x = formula_i_from_vector, 
               fixed = TRUE)){
        formula_vector_to_parse[i] <- gsub(pattern = paste0(" ", all_e_coli_metabolites$`Metabolite Abbreviation`[k]), 
                                           replacement = paste0(" ", all_e_coli_metabolites$`Metabolite Name`[k]), 
                                           x = formula_vector_to_parse[i], fixed = TRUE)
      }
    }
  }
  return(formula_vector_to_parse)
}

reaction_parsed_e_coli <- get_formula_full_names(gene_summary_df = genes_fluxes_summary_states_merged_w_formulas, 
                                                 all_e_coli_metabolites = all_e_coli_metabolites)

genes_fluxes_summary_states_merged_w_formulas$Formula_Full_Name <- reaction_parsed_e_coli

write.csv(genes_fluxes_summary_states_merged_w_formulas, 
          file = "data/e_coli_genes_fluxes_summary_states_merged_full_name.csv", 
          row.names = FALSE)


#### Overall summary of key dataframes:

# genes_fluxes_long_df_z_score_cc: This is a long-format df of paired fluxes with protein abundances.
# genes_fluxes_long_df_z_score_cc_2_level: This is a long-format df of paired fluxes with protein abundances, but filtered for those proteins that have paired information on weighted abundance (all), reaction deltaG, and substrate connectivity.
# genes_fluxes_summary_states_merged: This is a summary of the long-format df, that is used for inference on individual protein characteristics in the 2 level hierarchical model.

