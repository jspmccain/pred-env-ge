#### making a dataframe of yeast protein names

library(seqinr)

## read in the fasta file
yeast_fasta <- read.fasta('data/orf_trans.fasta')

## parse the annotations from the fasta header
annot_yeast <- getAnnot(yeast_fasta)

### need a function to go through all the protein annotations, take the interpretable protein name
### and then take the rest of the stuff on the first column

get_yeast_gene_names_and_descriptions <- function(annot_from_fasta_file){
  ## make instances of vectors to aggregate
  left_side_containing_name <- c()
  gene_description <- c()
  short_gene_description <- c()
  short_name <- c()
  
  ## loop through all the annotations
  for(i in 1:length(annot_from_fasta_file)){
    # i <- 1
    left_side_i <- strsplit(annot_from_fasta_file[[i]], "\"")[[1]][1]
    name_only_i <- strsplit(x = left_side_i, split = " ")[[1]][1]
    name_only_i_form <- gsub(pattern = ">", 
                             replacement = "", 
                             x = name_only_i)
    
    description_i <- strsplit(annot_from_fasta_file[[i]], "\"")[[1]][2]
    short_description_i <- strsplit(description_i, ";")[[1]][1]
    
    ## append these to the vectors
    short_name <- c(short_name, name_only_i_form)
    left_side_containing_name <- c(left_side_containing_name, 
                                   left_side_i)
    gene_description <- c(gene_description, 
                          description_i)
    short_gene_description <- c(short_gene_description, 
                                short_description_i)
    
  }
  
  df_of_names_and_desc <- data.frame(left_side_containing_name, 
                                     gene_description,
                                     short_desc = short_gene_description,
                                     gene_name = short_name)
  return(df_of_names_and_desc)
}

df_yeast_annots <- get_yeast_gene_names_and_descriptions(annot_from_fasta_file = annot_yeast)

write.csv(x = df_yeast_annots, file = "data/yeast_genes_annots.csv", row.names = FALSE)
