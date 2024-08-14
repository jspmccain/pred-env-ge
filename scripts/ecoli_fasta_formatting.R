library(seqinr)
library(magrittr)
'%!in%' <- function(x,y)!('%in%'(x,y))

blastp_out <- read.table(file = 'data/yeast_v_e_coli.blastp', 
                         col.names = c("qseqid", "sseqid", "pident", "length", 
                                       "mismatch", "gapopen", "qstart", "qend", 
                                       "sstart", "send", "evalue", "bitscore"))
blastp_out_argc_e_coli <- read.table(file = 'data/e_coli_argc_against_BW25113.blastp',
                                     col.names = c("qseqid", "sseqid", "pident", "length", 
                                                               "mismatch", "gapopen", "qstart", "qend", 
                                                               "sstart", "send", "evalue", "bitscore"))
blastp_out_psd_e_coli <- read.table(file = 'data/e_coli_psd_against_BW25113.blastp',
                                     col.names = c("qseqid", "sseqid", "pident", "length", 
                                                   "mismatch", "gapopen", "qstart", "qend", 
                                                   "sstart", "send", "evalue", "bitscore"))
e_coli_fasta <- read.fasta(file = 'data/GCF_004355105.2_ASM435510v2_translated_cds.faa')

fasta_header_to_gene_name_e_coli <- function(e_coli_fast_file){

  name_to_header_df <- data.frame(fasta_header = character(),
                                  gene_name = character(),
                                  protein_length = numeric())
  
  # e_coli_fast_file <- e_coli_fasta
  
  for(i in 1:length(e_coli_fast_file)){
    
    ## subsetting the fasta characteristics from the ith entry
    element_i_val <- attributes(e_coli_fast_file[[i]])
    features_i <- element_i_val["Annot"]$Annot
    fasta_header_i <- element_i_val["name"]$name
    protein_length_i <- e_coli_fast_file[[i]] %>% length()
    
    ## formatting the fasta string to get out the gene name value
    
    ### Split up the features by a space
    parsed_features <- str_split(string = features_i, pattern = " ")
    
    ### Get the split up substring that designates the gene string
    if(sum(grepl(pattern = "[gene=", 
                 x = parsed_features[[1]], fixed = TRUE)) > 0){
      gene_name_i <- parsed_features[[1]][grepl(pattern = "[gene=", 
                                           x = parsed_features[[1]], fixed = TRUE)]
    } else {
      gene_name_i <- sapply(parsed_features, "[[", 2)
    }

    ### Format the gene name
    gene_name_i_form1 <- gsub(pattern = "[gene=", replacement = "", 
                              x = gene_name_i, fixed = TRUE)
    gene_name_i_form2 <- gsub(pattern = "]", replacement = "", 
                              x = gene_name_i_form1, fixed = TRUE)
    ### Make a temporary dataframe
    temp_df <- data.frame(fasta_header = fasta_header_i,
                          gene_name = gene_name_i_form2,
                          protein_length = protein_length_i)
    
    ### Append that dataframe
    name_to_header_df <- rbind(name_to_header_df, 
                               temp_df)
  }
  return(name_to_header_df)
}

e_coli_header_names_df <- fasta_header_to_gene_name_e_coli(e_coli_fast_file = e_coli_fasta)

#### now need to use data from davidi-data-processing.R specifically to get the b number - to - protein letter name mapping

### there are some issues with using the gene names in letters, and specifically, we do not have equivalent letters in the genome file.
unmapped_gene_names <- e_coli_b_number_gene_name_mappings$gene_name_lett[e_coli_b_number_gene_name_mappings$gene_name_lett %!in% e_coli_header_names_df$gene_name]

# from the list of b numbers to gene name mappings, which ones are not mapped in the e coli header names df
e_coli_b_number_gene_name_mappings[e_coli_b_number_gene_name_mappings$gene_name_lett %!in% e_coli_header_names_df$gene_name,]

# e_coli_header_names_df %>% 
#   # filter(gene_name == 'trpGD')
#   filter(grepl('ps', gene_name))
# e_coli_header_names_df %>% 
#   filter(protein_length == 322)

### manually adjusting the names from the fasta file to match the names from the Davidi et al 2016 sheet. Double checked these on EcoCyc by looking at the bnumber and the amino acid number, and specifically identifying the name in the fasta file by looking at the synonym component of EcoCyc.

## e.g. in the fasta file, the gene ispC but in Davidi et al it's called dxr
e_coli_header_names_df[e_coli_header_names_df$gene_name == 'ispC',]$gene_name <- 'dxr' 
e_coli_header_names_df[e_coli_header_names_df$gene_name == 'bioD',]$gene_name <- c('bioD1', 'bioD1')
e_coli_header_names_df[e_coli_header_names_df$gene_name == 'trpCF',]$gene_name <- 'trpC'
e_coli_header_names_df[e_coli_header_names_df$gene_name == 'trpD',]$gene_name <- 'trpGD'
e_coli_header_names_df[e_coli_header_names_df$gene_name == 'ribE',]$gene_name <- c('ribC', 'ribE')
e_coli_header_names_df[e_coli_header_names_df$gene_name == 'hisIE',]$gene_name <- 'hisI'
e_coli_header_names_df[e_coli_header_names_df$gene_name == 'gndA',]$gene_name <- 'gnd'
e_coli_header_names_df[e_coli_header_names_df$gene_name == 'iscU',]$gene_name <- 'nifU'
e_coli_header_names_df[e_coli_header_names_df$gene_name == 'nadK',]$gene_name <- 'ppnK' # note this is a bit confusing. In Schmidt et al, they interchange using ppnK and nadK, and this cascades to Davidi et al, who use ppnK.

#### Unclear if E. coli this strain has an argC protein -- this seems to not be the case based on BLASTP the argC protein from E. coli MG1655 sequence to the E. coli strain used in Schmidt et al. 
#### Unclear if E. coli this strain has psd protein as well. A blastp does not reveal any homologues using this genome.

yeast_gene_to_ecoli_gene_connection <- blastp_out %>% 
  dplyr::rename(fasta_header = qseqid) %>% 
  inner_join(e_coli_header_names_df, by = "fasta_header") %>% 
  dplyr::filter(evalue < 1e-30) %>% 
  dplyr::select(sseqid, gene_name) %>% 
  unique() %>% 
  dplyr::rename(gene_name_lett = gene_name) %>% 
  inner_join(e_coli_b_number_gene_name_mappings, by = "gene_name_lett")


