library(magrittr)
library(dplyr)
library(ggplot2)

yeast_lengths <- read.csv("data/yeast_S288C_lengths.csv")

load('data/20130313ProtPepMatrices.Rdata')
# load('data/20130313ProtPepMatrices.Rdata', ex <- new.env())
# ls.str(ex)

# heavyIC is the heavy standard
# lightIC is the experimental 
# PepMatrix is the log2 transformed ratio of the light/heavy

## We are interested in the experimental treatments, and getting a protein-level estimate

## First, we need to normalize the data, so that variation in protein loaded on the column is standardized.
## This is not the case in the raw data:

# lightIC %>% colSums(na.rm = TRUE) %>% hist() ## These values are not all constant.

## Normalize for the sum of all peptides ionized.
## First calculate normalization factors:

normalization_factors_col_wise_light <- colSums(lightIC, na.rm = TRUE)

## Now divide every column by this vector:

norm_light <- sweep(lightIC, 2, 
                    normalization_factors_col_wise_light, 
                    FUN = '/')

## These should now all sum to one:

# colSums(norm_light, na.rm = TRUE)

## There is one column that indicates that these samples were collected on an orbitrap, not a QTOF. This is removed here to keep the MS modality consistent.
norm_light_no_orbi <- norm_light[, -which(colnames(norm_light) == "P0.30_ORBI")]

## Now make a function that inputs a specific protein and experimental condition. With that protein, we then want to subset all the 
## corresponding peptides. We want to ensure that there are at least x (x = 5, perhaps) corresponding peptides. Get the peptide
## abundances for that protein, and then choose the top three. Then take the mean, which gives the absolute abundance estimate.
## Then we want to apply this to every experimental condition, to get a cross-experiment measurement of protein abundance.

# input_protein <- "YBR011C"
# pepprot_mapping <- ProtPepMatrix


# subset_target_protein_peptide %>% names()

get_corresponding_peptides <- function(input_protein, pepprot_mapping = ProtPepMatrix){
  ## get the column that corresponds
  target_protein_peptide_mapped <- pepprot_mapping[, input_protein]
  ## subset only the entries that have a 1 not a 0
  subset_target_protein_peptide <- target_protein_peptide_mapped[which(target_protein_peptide_mapped == 1)]
  ## return the names of this, which are the peptide sequences. (Note, these have been gut checked with SGD)
  return(names(subset_target_protein_peptide))
}

# trx1 <- get_corresponding_peptides(input_protein = "YLR043C")
# trx2 <- get_corresponding_peptides(input_protein = "YGR209C")
# 
# intersect(trx1, trx2)
# 
# trr1 <- get_corresponding_peptides(input_protein = "YDR353W")
# trr2 <- get_corresponding_peptides(input_protein = "YHR106W")
# 
# intersect(trr1, trr2)
# 
# get_corresponding_peptides(input_protein = "YDR300C")
# get_corresponding_peptides(input_protein = "YHR033W")
# 
# get_corresponding_peptides(input_protein = "YJL101C")

get_abundances_for_a_condition <- function(peptide_set, condition_chosen, normalized_peptide_abundance_matrix = norm_light_no_orbi){
  ## subset the peptide abundances using the determined peptide set and the condition chosen
  subsetted_abundances <- normalized_peptide_abundance_matrix[peptide_set, condition_chosen]
  ## now get all the abundances that are non NA and positive
  number_of_non_na_or_non_zero_obs <- subsetted_abundances[!is.na(subsetted_abundances) & subsetted_abundances > 0]
  ## return these abundances
  return(number_of_non_na_or_non_zero_obs)
}

# abundances_tester <- get_abundances_for_a_condition(peptide_set = peptides_tester, condition_chosen = condition_test)

top3_for_condition <- function(input_protein, condition_chosen, number_peps_cutoff = 4){
  
  ## first get the peptides
  peptides_i <- get_corresponding_peptides(input_protein = input_protein)
  
  ## then get the abundances for the corresponding peptides
  peptide_abundances_i <- get_abundances_for_a_condition(peptide_set = peptides_i, 
                                                         condition_chosen = condition_chosen)
  
  ## if there are very few peptides per protein, based on the cutoff, then return an NA
  if(length(peptide_abundances_i) < number_peps_cutoff){
    return(NA)
  }
  else {
    ## get the top three peptide abundances
    top_three_peps <- peptide_abundances_i[rank(peptide_abundances_i) < 4]
    return(mean(top_three_peps))
  }
}

mass_fraction_for_condition <- function(input_protein, condition_chosen){
  ## first get the peptides
  peptides_i <- get_corresponding_peptides(input_protein = input_protein)
  ## then get the abundances for the corresponding peptides
  peptide_abundances_i <- get_abundances_for_a_condition(peptide_set = peptides_i, 
                                                         condition_chosen = condition_chosen)
  return(sum(peptide_abundances_i))
}

# top3_for_condition(input_protein = input_protein, condition_chosen = condition_test)

# get_corresponding_peptides(input_protein = "YBR011C")
vector_conditions_test <- colnames(norm_light_no_orbi)

get_protein_abundances_across_conditions <- function(input_protein, 
                                                     vector_of_conditions, 
                                                     number_peps_cutoff = 5){
  ## collect protein estimates across conditions into a vector
  condition_specific_protein_estimates <- c()
  
  ## collect protein estimates across conditions into a vector
  condition_specific_mass_fraction <- c()
  
  ## loop through and get a protein abundance estimate across all conditions and replicates
  for(i in 1:length(vector_of_conditions)){
    ## get the single top3 estimate
    top3_estimate_condition_i <- top3_for_condition(input_protein = input_protein, 
                                                     condition_chosen = vector_of_conditions[i], 
                                                     number_peps_cutoff = number_peps_cutoff)
    mf_estimate_condition_i <- mass_fraction_for_condition(input_protein = input_protein,
                                                          condition_chosen = vector_of_conditions[i])
    # aggregate these
    condition_specific_protein_estimates <- c(condition_specific_protein_estimates, top3_estimate_condition_i)
    condition_specific_mass_fraction <- c(condition_specific_mass_fraction, mf_estimate_condition_i)
  }
  ## return that! woop!
  return(list(condition_specific_protein_estimates, condition_specific_mass_fraction))
}

get_protein_abundances_across_conditions_and_proteins <- function(input_protein_vector, 
                                                                  vector_of_conditions, 
                                                                  number_peps_cutoff = 5){
  
  ## go through all the input proteins and get the vector of protein abundances
  protein_abundances_df <- data.frame(top3_abund = numeric(),
                                      mf_abund = numeric(),
                                      condition = character(),
                                      protein_name = character())
  
  ## loop through all the proteins in the input protein vector
  for(i in 1:length(input_protein_vector)){
    protein_i <- input_protein_vector[i]
    protein_abundance_vecs_i <- get_protein_abundances_across_conditions(input_protein = protein_i, 
                                                                        vector_of_conditions = vector_of_conditions)
    temp_df_i <- data.frame(top3_abund = protein_abundance_vecs_i[[1]],
                            mass_frac = protein_abundance_vecs_i[[2]],
                            condition = vector_of_conditions, 
                            protein_name = rep(protein_i, length(vector_of_conditions)))
    protein_abundances_df <- rbind(protein_abundances_df, 
                                   temp_df_i)
  }
  return(protein_abundances_df)
}


## this dataframe is from hackett-data-processing.R
vector_of_proteins_from_subset <- mean_protein_amount_singles_w_gene_factor$Gene %>% unique()

## now run this function across all proteins that are in the subset, and across all conditions
abundances_out_conditions <- get_protein_abundances_across_conditions_and_proteins(input_protein_vector = vector_of_proteins_from_subset, 
                                                                                   vector_of_conditions = vector_conditions_test)

## summarize these to get a single mean estimate across conditions
summarized_abundances_across_conditions <- abundances_out_conditions %>% 
  group_by(protein_name) %>% 
  summarize(mean_prot_across_conditions = mean(top3_abund, na.rm = TRUE),
            sd_prot_across_conditions = sd(top3_abund, na.rm = TRUE),
            mean_prot_mf = mean(mass_frac, na.rm = TRUE),
            sd_prot_mf = sd(mass_frac, na.rm = TRUE),
            number_non_nas_across_conditions = sum(!is.na(top3_abund))) %>% 
  dplyr::mutate(mean_prot_mf_percentage = mean_prot_mf*100,
                log_mean_prot_mf = log(mean_prot_mf),
                log_mean_prot_mf_perc = log(mean_prot_mf_percentage)) 
  
