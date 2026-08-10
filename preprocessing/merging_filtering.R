library(Seurat)
library(magrittr)
library(tidyverse)
library(data.table)
library(pbmcapply)


save_name <- 'TIERS'
#tumor_type <- 'SCC'

getwd()

#### load data and determine variable genes ############################
dat <- readRDS('./LungCAIRE_merged_annotation_fine.Rds')

variable_genes <- VariableFeatures(dat, nfeatures=5000) 
all_genes <- dat %>% Features() %>% unique()


########################################################################
############ prepare cell meta data #####################################
cell_meta_data <- dat@meta.data
cell_meta_data$id = rownames(cell_meta_data)

View(head(cell_meta_data))

########################################################################
############ prepare cell rna data #####################################
array_all <- dat@assays$RNA


########################################################################
#################prepare file names and functions for loops
file_names <- names(array_all@layers)
file_names_counts <- file_names[startsWith(file_names, 'counts')]
file_names_normalized <- file_names[startsWith(file_names, 'data')]


get_data <- function(current_name, gene_names = NULL) {
  cell_identifier <- data.frame('id' = colnames(array_all[current_name]))
  dense_matrix <- as.matrix(GetAssayData(array_all, slot = current_name)) %>% t()
  filtered_dense_matrix <- dense_matrix[,colnames(dense_matrix) %in% gene_names]
  print(filtered_dense_matrix %>% dim())
  specific_cell_meta_data <- cell_identifier %>% left_join(cell_meta_data) %>% 
    #dplyr::select(id, Pseudo = orig.ident, Tier_1, Tier_2, Tier_3, Tier_4)
    dplyr::select(id, Pseudo = orig.ident)
  
  cbind(specific_cell_meta_data, filtered_dense_matrix)
  }

############################################################################################
#prepare dataframes with rna data and cell meta data

filtered_res_normalized <- rbindlist(pbmclapply(file_names_normalized, function(x) get_data(x,variable_genes), mc.cores=20L))



library(arrow)
#write_parquet(filtered_res_normalized,paste0('./data/MERGED_normalized_5000genes',save_name,'.parquet'))

