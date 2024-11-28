#### examining relationship between fluxes and transcript abundances from Chubukov et al

library(ggplot2)
library(dplyr)
library(magrittr)
library(readxl)

# read in data and format fluxes and transcript values ------------------------------------------------------------

## This data contains the transcript abundance
chubukov_transcripts <- read_excel('data/b_subtilis_data/TableS2.xlsx')

## This data contains the growth condition names
chubukov_cond <- read_excel('data/b_subtilis_data/TableS1.xlsx')
## Making the naming easier to work with
names(chubukov_cond)[1] <- 'chip_n' 
names(chubukov_cond)[2] <- 'chip_id'
names(chubukov_cond)[3] <- 'condition'

## Cnverting the chip number to numeric
chubukov_cond$chip_n <- as.numeric(chubukov_cond$chip_n)

## We want to only subset the chubukov data treatments
chubukov_cond_subset <- chubukov_cond %>% 
  dplyr::filter(chip_n > 84, chip_n < 109)

## There are eight treatments in this study, with three replicates each. Therefore there should be 24
## rows in this subset
nrow(chubukov_cond_subset) == 24

## We want the chip id numbers from this chubukov_cond_subset dataframe because those will tell us which 
## columns to subset from the chubukov_transcripts dataframe
chip_id_to_subset <- chubukov_cond_subset$chip_id

# this is a unique experimental identifier, so it should equal 24
unique(chip_id_to_subset) %>% length() == 24

## Now filter the transcript dataframe if it contains these strings.
chubukov_transcripts_sub <- chubukov_transcripts %>% 
  dplyr::select(contains(c('Name', 
                           'Locus_tag', 
                           chip_id_to_subset))) %>% 
  ## Then convert this to a long-form dataframe
  melt(id.vars = c('Name', 'Locus_tag'), 
       variable.name = 'condition_chip_id', 
       value.name = 'expression_val')

all_transcript_locus_tags <- chubukov_transcripts_sub$Locus_tag %>% unique()
all_transcript_gene_names <- chubukov_transcripts_sub$Name %>% unique()

## Read in the flux data
chubukov_flux <- read_excel('data/b_subtilis_data/msb201366-sup-0002.xlsx', sheet = 4, skip = 3)
names(chubukov_flux)[1] <- 'enzyme_bsu_numbers'
names(chubukov_flux)[2] <- 'flux_name'

## The gene names have capital letters, here I convert the capitals to lower case.
chubukov_flux$flux_name <- tolower(chubukov_flux$flux_name)

names(chubukov_flux) <- c('enzyme_bsu_numbers', 'flux_name', 'flux', 'Glu', 'Fru', 'Glucon', 'G_S', 'Gly', 'Mal', 'M_G', 'Pyr',
                          'standard_dev', 'Glu_sd', 'Fru_sd', 'Glucon_sd', 'G_S_sd', 'Gly_sd', 'Mal_sd', 'M_G_sd', 'Pyr_sd')

## Getting the fluxes that have single locus tag matches. This excludes several fluxes, either
## that are associated with multiple gene loci e.g. BSU37050;BSU23550;BSU29880;BSU29220, or that are not
## associated with any gene loci.
chubukov_flux_sub <- chubukov_flux %>% 
  dplyr::filter(enzyme_bsu_numbers %in% all_transcript_locus_tags) %>% 
  dplyr::select(-flux, -standard_dev, -contains('sd')) %>% 
  melt(id.vars = c('enzyme_bsu_numbers', 'flux_name'),
       variable.name = 'condition', 
       value.name = 'flux')

## Now getting the standard deviations associated with those measured fluxes.
chubukov_flux_sub_sd <- chubukov_flux %>% 
  dplyr::filter(enzyme_bsu_numbers %in% all_transcript_locus_tags) %>% 
  dplyr::select(-flux, -standard_dev, -Glu, -Fru, -Glucon, -G_S, -Gly, -Mal, -M_G, -Pyr,
                enzyme_bsu_numbers, 
                flux_name, 
                contains('sd')) %>% 
  melt(id.vars = c('enzyme_bsu_numbers', 'flux_name'),
       variable.name = 'condition_sd', 
       value.name = 'flux_sd')

## Making a joining column that has the condition name, instead of condition_sd, as a column
chubukov_flux_sub_sd$condition <- gsub(pattern = "_sd", 
                                       replacement = "", 
                                       x = chubukov_flux_sub_sd$condition_sd)

## The dimensions of the Flux df and the SD Flux df should be the same
nrow(chubukov_flux_sub_sd) == nrow(chubukov_flux_sub)

## The order of enzyme bsu numbers should be the same
chubukov_flux_sub_sd$enzyme_bsu_numbers == chubukov_flux_sub$enzyme_bsu_numbers

## Merging the long dataframes to include both flux and flux SD
chubukov_flux_long <- chubukov_flux_sub %>% 
  inner_join(chubukov_flux_sub_sd, by = c('enzyme_bsu_numbers', 'flux_name', 'condition'))

## Now I need to merge the formatted flux dataframe with the transcript data, which is in this dataframe: chubukov_transcripts_sub.

## First I need to format the transcript condition column so it matches the flux condition column
chubukov_transcripts_sub$condition1 <- substr(as.character(chubukov_transcripts_sub$condition_chip_id), 1, 
                                             nchar(as.character(chubukov_transcripts_sub$condition_chip_id)) - 12)

## Divide the condition string into two, the first will be the description of the carbon source, the second will be the replicate.
string_split_chubukov <- str_split(string = chubukov_transcripts_sub$condition1, 
                                        pattern = "_")
## Then subset only the condition name.
condition_name_chubukov <- sapply(string_split_chubukov, 
                               "[[", 
                               1)
## Then subset only the condition name.
replicate_number_chubukov <- sapply(string_split_chubukov, 
                                  "[[", 
                                  2)
## Removing the + in the condition description
chubukov_transcripts_sub$condition <- gsub(pattern = "+", 
                                           replacement = "_", 
                                           x = condition_name_chubukov, fixed = TRUE)

## Adding the replicate number
chubukov_transcripts_sub$replicate_number <- as.numeric(replicate_number_chubukov)

## Averaging the transcript abundance across replicates
chubukov_transcripts_aggregated <- chubukov_transcripts_sub %>% 
  rename(enzyme_bsu_numbers = Locus_tag) %>% 
  group_by(condition, enzyme_bsu_numbers) %>% 
  summarize(mean_expression_val = mean(as.numeric(expression_val), na.rm = TRUE))

## Add the aggregated transcript dataframe to the flux dataframe.
chubukov_transcripts_sub_only_in_flux <- chubukov_flux_long %>% 
  inner_join(chubukov_transcripts_aggregated, 
             by = c('enzyme_bsu_numbers', 'condition'))

# transform flux values and transcript values ----------------------------

## Filter out the fluxes that are zero. The mean of the **logged** values is taken.
chu_trans_exp_flux_summaries <- chubukov_transcripts_sub_only_in_flux %>% 
  # dplyr::filter(flux > 0) %>%
  dplyr::filter(enzyme_bsu_numbers != "BSU17890") %>%  ## this enzyme is taken out because it maps to *both* Tkt1 and Tkt2
  group_by(flux_name) %>% 
  summarize(mean_flux_across_cond = mean(flux, na.rm = TRUE),
            mean_log_exp_across_cond = mean(log2(mean_expression_val), na.rm = TRUE),
            sd_flux_across_cond = sd(flux, na.rm = TRUE),
            sd_log_exp_across_cond = sd(log2(mean_expression_val), na.rm = TRUE))

## Z-score transform the fluxes and transcript abundances
chu_transformed <- chubukov_transcripts_sub_only_in_flux %>% 
  inner_join(chu_trans_exp_flux_summaries, 
             by = 'flux_name') %>%
  rowwise() %>%
  mutate(z_score_exp = (log2(mean_expression_val) - mean_log_exp_across_cond)/sd_log_exp_across_cond,
         z_score_flux = (flux - mean_flux_across_cond)/sd_flux_across_cond)

## remove infinity and NA values
chu_transformed_cc <- chu_transformed[!is.na(chu_transformed$sd_flux_across_cond) & 
                                        is.finite(chu_transformed$z_score_flux) &
                                        is.finite(chu_transformed$z_score_exp), ]

chu_transformed_cc$flux_name_as_number <- as.numeric(as.factor(chu_transformed_cc$flux_name))

## plotting raw data output
chu_transformed_cc %>% 
  ggplot(aes(x = z_score_exp, 
             group = flux_name, 
             y = z_score_flux)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  facet_wrap(~flux_name) +
  ylim(-2.5, 2.5)

# fitting the model in stan -----------------------------------------------

fitness_model_partial_bsub <- stan_model("scripts/flux_protein_hierarchical_partial_pooling.stan")

options(mc.cores = 4)

number_of_proteins_chu <- chu_transformed_cc$flux_name_as_number %>% unique() %>% length()
protein_amounts_spaced_sampled_chu <- rep(seq(from = -2, to = 2, by = 0.05), number_of_proteins_chu)
protein_ids_spaced_sampled_chu <- rep(c(1:number_of_proteins_chu), each = 81)

## fitting the model in stan
chu_partial_pool <- sampling(fitness_model_partial_bsub, 
                                     data = list(N = nrow(chu_transformed_cc),
                                                 K = length(unique(chu_transformed_cc$flux_name)),
                                                 id = chu_transformed_cc$flux_name_as_number,
                                                 flux = chu_transformed_cc$z_score_flux,
                                                 protein_amounts = chu_transformed_cc$z_score_exp,
                                                 protein_amounts_spaced = protein_amounts_spaced_sampled_chu,
                                                 N_spaced = length(protein_amounts_spaced_sampled_chu),
                                                 id_spaced = protein_ids_spaced_sampled_chu), 
                                     iter = 10000, 
                             chains = 10)

beta_1_estimates_chu_df <- chu_partial_pool %>% 
  spread_draws(beta_1[n]) %>% 
  median_qi()

## posterior predictive check
posterior_pred_chu_obs <- chu_partial_pool %>% 
  spread_draws(y_rep[n]) %>% 
  median_qi()

# making a dataframe that has the posterior predictive checks in there
chu_transformed_cc_w_pp <- chu_transformed_cc
chu_transformed_cc_w_pp$y_rep <- posterior_pred_chu_obs$y_rep
chu_transformed_cc_w_pp$lower_val <- posterior_pred_chu_obs$.lower
chu_transformed_cc_w_pp$upper_val <- posterior_pred_chu_obs$.upper

## plotting posterior predictions for Chubukov et al dataset
posterior_predictions_check_chubukov <- chu_transformed_cc_w_pp %>% 
  ggplot(aes(x = z_score_exp, y = y_rep)) +
  geom_ribbon(aes(ymin = lower_val, ymax = upper_val), alpha = 0.3) +
  facet_wrap(~flux_name) +
  geom_line() +
  theme_bw() +
  geom_point(data = chu_transformed_cc_w_pp, 
             aes(x = z_score_exp, y = z_score_flux), alpha = 0.5) 

ggsave(posterior_predictions_check_chubukov, 
       filename = 'figures/posterior_predictions_check_chubukov.png', height = 11.1, width = 11.3)
