library(data.table)
library(dplyr)
library(magrittr)
library(ggplot2)
library(cluster)
library(tidyverse)
library(Rtsne)
library(nnet)
library(Hmisc)
library(glmnet)
library(arrow)
library(survival)
library(survminer)
library(survcomp)
library(caret)
library(coxme)
library(Seurat)
library(ggrepel)
library(pbmcapply)
library(ggsankey)
library(randomForestSRC)
library(slingshot)
library(dyno)
library(viridis)
library(cowplot)


getwd()
setwd('plots_statistics')

cancer_type <- 'AC'

patient_meta_data <- read.csv('../data/Full_cohort_with_clinical_data.csv') 

clean_cell_annotation <- read.csv('./use_data/VICREG_clustering_highinputdropout_5000epochs_rev_annotation.csv') %>% 
  mutate(Tier_2 = ifelse(Tier_3 == 'Tumor', 'Tumor', Tier_2)) %>% 
  filter(Tier_3 != 'Mixed') 

cluster_result <- read_parquet('../scGPT_embeddings/clusterings/VICREG_clustering_highinputdropout_5000epochs.parquet') %>%
  mutate(cell_position = seq(dim(.)[1])) %>% 
  dplyr::select(-c(Tier_1, Tier_2, Tier_3, Tier_4)) %>%  
  left_join(clean_cell_annotation, by = c('cluster' = 'Cluster_ID')) %>%                
  filter(Tier_2 == 'Tumor')                                

tumor_cell_coefficients_annotated <- read.csv(paste0('./use_data/important_clusters_', cancer_type, '.csv')) %>% 
  mutate(cluster = as.integer(gsub('X', '', cluster))) 


tumor_cell_coefficients_annotated_MP <- read.csv(paste0('./use_data/important_clusters_', cancer_type, '_MP.csv')) %>% 
  mutate(cluster = as.integer(gsub('X', '', cluster))) 

cell_expression <- read_parquet('../data/MERGED_normalized_5000genesTIERS.parquet')

embeddings <- read_parquet('../scGPT_embeddings/embeddings/VICREG_embedding_highinputdropout_5000epochs.parquet')

malignant_ranking <- read.csv( './use_data/latent_metastatic_cells_candidates.csv') %>% 
  mutate(cluster = as.integer(gsub('X', '', cluster))) %>% 
  dplyr::select(-X)

tropism_scores <- read.csv('./use_data/Cell_clusters_tropism_scores_AC.csv') %>% 
  dplyr::select(-X)


pseudotimes <- fread('./use_data/tumor_dpt.csv') %>%  
  dplyr::select(-Tier_3, -V1) %>%
  filter(dpt_pseudotime !=0) %>%
  left_join(patient_meta_data %>% dplyr::select(Pseudo, Entity, UICC8_edition2)) %>%
  inner_join(cluster_result %>% dplyr::select(id, Pseudo, cluster), by =c('cell_id' = 'id', 'Pseudo', 'cluster')) %>% 
  left_join(malignant_ranking) %>%
  mutate(UICC8_edition3 = ifelse(UICC8_edition2 %in% c('IA', 'IB'), 'I', UICC8_edition2)) %>% 
  filter(!is.na(malignancy_rank)) %>%
  filter(Pseudo!= 'LC395') 


within_patient_time_vs_malignancy <- pseudotimes %>% 
  mutate(rounded_dpt_pseudotime = round(dpt_pseudotime,1)) %>% 
  group_by(rounded_dpt_pseudotime, Pseudo, Entity, UICC8_edition3) %>% 
  dplyr::summarize(malignancy_rank = mean(malignancy_rank)) %>%
  group_by(Pseudo, Entity, UICC8_edition3) %>% 
  mutate(N=n()) %>% filter(N>4) %>%
  dplyr::summarize(r = (rcorr(rounded_dpt_pseudotime, malignancy_rank, type = 'spearman')$r[1,2]))


within_patient_time_vs_malignancy %>% 
  group_by(Entity, UICC8_edition3) %>% 
  #group_by(Entity) %>% 
  dplyr::summarize(p = wilcox.test(r)$p.value)

ggplot(within_patient_time_vs_malignancy, aes(x = UICC8_edition3, y = r)) +
  geom_boxplot() +
    geom_point() +
  facet_wrap(~Entity) +
  geom_hline(yintercept = 0) +
  theme_minimal()

#result for manuscript
within_patient_time_vs_malignancy %>% 
  group_by(Entity) %>% 
  dplyr::summarize(mean_r = mean(r), P = wilcox.test(r)$p.value) %>% 
  arrange(Entity)


#################################################
#calculate key marker genes for malignancy scores
#################################################
cluster_result_with_malignancy <- cluster_result %>% 
  left_join(patient_meta_data %>% dplyr::select(Pseudo, Entity)) %>%
  left_join(malignant_ranking)

tumor_cell_expression <- cell_expression %>% 
  filter(Tier_3 == 'Tumor') %>%  
  dplyr::select(-c('Tier_1', 'Tier_2', 'Tier_3', 'Tier_4'))

cell_expression_and_malignancy <- cluster_result_with_malignancy %>% 
  inner_join(tumor_cell_expression)

cell_expression_and_malignancy_AC <- cell_expression_and_malignancy %>% filter(Entity == 'AC', !is.na(r)) #%>% .[1:100000,]
lasso_model_AC <- glmnet(cell_expression_and_malignancy_AC[,-c(1:11)], cell_expression_and_malignancy_AC$r, alpha = 1, lambda = 1e-3) 
coefs_AC <- coef(lasso_model_AC) %>% as.data.frame() %>% 
  filter(s0!=0)
cell_expression_AC_filtered <- cell_expression_and_malignancy_AC %>% 
  dplyr::select(colnames(.)[colnames(.) %in% rownames(coefs_AC)])
lasso_model_AC_clean <- glmnet(cell_expression_AC_filtered, cell_expression_and_malignancy_AC$r, alpha = 0, lambda = 0) 
clean_coefs_AC <- data.frame(coef(lasso_model_AC_clean)) %>% 
  mutate(gene = rownames(.), cancer_type = 'AC') %>% 
  arrange((s0)) %>% 
  mutate(gene = factor(gene, levels = gene)) %>% 
  filter(gene != '(Intercept)')

coefficient_plot_AC <- ggplot(clean_coefs_AC, aes(x = s0, y = gene)) +
  geom_segment(aes(xend = 0, yend=gene)) +
  geom_point(size=3) +
  geom_vline(xintercept=0) +
  theme_minimal() +
  theme(text=element_text(size=20), 
    axis.title = element_blank())


cell_expression_and_malignancy_SCC <- cell_expression_and_malignancy %>% filter(Entity == 'SCC', !is.na(r))
lasso_model_SCC <- glmnet(cell_expression_and_malignancy_SCC[,-c(1:11)], cell_expression_and_malignancy_SCC$r, alpha = 1, lambda = 2e-3) 
coefs_SCC <- coef(lasso_model_SCC) %>% as.data.frame() %>% 
  filter(s0!=0)
cell_expression_SCC_filtered <- cell_expression_and_malignancy_SCC %>% 
  dplyr::select(colnames(.)[colnames(.) %in% rownames(coefs_SCC)])
lasso_model_SCC_clean <- glmnet(cell_expression_SCC_filtered, cell_expression_and_malignancy_SCC$r, alpha = 0, lambda = 0) 
clean_coefs_SCC <- data.frame(coef(lasso_model_SCC_clean)) %>% 
  mutate(gene = rownames(.), cancer_type = 'SCC') %>% 
  arrange((s0)) %>% 
  mutate(gene = factor(gene, levels = gene)) %>% 
  filter(gene != '(Intercept)')

coefficient_plot_SCC <- ggplot(clean_coefs_SCC, aes(x = s0, y = gene)) +
  geom_segment(aes(xend = 0, yend=gene)) +
  geom_point(size=3) +
  geom_vline(xintercept=0) +
  theme_minimal() +
  theme(text=element_text(size=20), 
    axis.title = element_blank())

plot_grid(coefficient_plot_AC, coefficient_plot_SCC)


gene_expression_dat <- tumor_cell_expression %>% 
  dplyr::select(cell_id = id, Pseudo, gene_expression = LAMC2)
gene_expression_dat


###############
#plot it
################


all_patients <- pseudotimes$Pseudo %>% unique()
one_patient <- pseudotimes %>% 
  filter(Pseudo %in%all_patients[1:5]) %>% #296:300
  left_join(gene_expression_dat)

dim(one_patient)

one_patient_plot <- rbind(one_patient %>% mutate(type = 'time', value = dpt_pseudotime), 
                        one_patient %>% mutate(type = 'malignancy', value = malignancy_rank), 
                        one_patient %>% mutate(type = 'malignancy', value = log(1+gene_expression)))

time_plot <- ggplot(one_patient, aes(x = UMAP1, y = UMAP2, color = dpt_pseudotime)) + 
  geom_point(size=1) +
  facet_wrap(~Pseudo, ncol=1) +
  scale_color_viridis() +
  theme_minimal()

malignancy_plot <- ggplot(one_patient, aes(x = UMAP1, y = UMAP2, color = malignancy_rank)) + 
  geom_point(size=1) +
  facet_wrap(~Pseudo, ncol=1) +
  #scale_color_gradient(low = 'gray80', mid = 'darkred')+
  scale_color_gradientn(colors=c('gray80', mid = 'darkred', 'red'))+
  theme_minimal()

gene_expression_plot <- ggplot(one_patient, aes(x = UMAP1, y = UMAP2, color = log(1+gene_expression))) + 
  geom_point(size=1, alpha=0.4) +
  facet_wrap(~Pseudo, ncol=1) +
  #scale_color_gradient(low = 'gray80', mid = 'darkred')+
  scale_color_gradientn(colors=c('gray', mid = 'darkblue', '#1818df'))+
  theme_minimal()

png('./figures/dpt_and_malignancy.png', width=4500, height=3000, res=200)
plot_grid(time_plot, malignancy_plot, gene_expression_plot, nrow=1)
dev.off()

pseudotimes_quantiled <- pseudotimes %>% 
  #mutate(coarse_dpt = floor(dpt_pseudotime*10)/10) %>% 
  mutate(quantile_pseudotime = ntile(dpt_pseudotime,10)) %>% 
  dplyr::mutate(UICC8_edition3 = ifelse(UICC8_edition2 %in% c('IA', 'IB'), 'I', UICC8_edition2)) %>%
  group_by(quantile_pseudotime, Entity, Pseudo, UICC8_edition3) %>% 
  dplyr::summarize(average_malignancy_rank = mean(malignancy_rank))  %>% 
  mutate(cancer_type = ifelse(Entity == 'AC', 'LUAD', 'LUSC'))

#for manuscript
pseudotimes_quantiled_correlations <- pseudotimes_quantiled %>%  
    group_by(Entity, Pseudo, UICC8_edition3) %>%
    mutate(N=n()) %>% filter(N>4) %>%
   dplyr::summarize(r = (rcorr(quantile_pseudotime, average_malignancy_rank, type = 'spearman')$r[1,2])) 

pseudotimes_quantiled_correlations %>%    
  #group_by(Entity) %>% 
  group_by(Entity, UICC8_edition3) %>% 
   dplyr::summarize(average_r = mean(r))



pseudotimes_quantiled_correlations %>% 
  #group_by(Entity, UICC8_edition3) %>% 
  group_by(Entity) %>% 
  dplyr::summarize(p = wilcox.test(r)$p.value)



png('./figures/dpt_and_malignancy_line_plot.png', height=1500, width=4000, res=200)
ggplot(pseudotimes_quantiled %>% filter(!is.na(UICC8_edition3)), aes(x = as.factor(quantile_pseudotime), y = average_malignancy_rank, fill = cancer_type)) +
  geom_boxplot(width=0.5, alpha=0.8) +
  #facet_wrap(~UICC8_edition3, nrow=2)+
  facet_grid(cancer_type~UICC8_edition3)+
  scale_fill_manual(values = c('#ff7700', '#1414ad')) +
  theme_minimal() +
  theme(text=element_text(size=25), 
    legend.title = element_blank()) +
  xlab('Pseudotime (quantiles)') +
  ylab('Estimated metastatic potential')
dev.off()

pdf(paste0('./figures/dpt_vs_malignancy_boxplots.pdf'), width = 15, height = 15, onefile = FALSE)
ggplot(pseudotimes_quantiled %>% filter(!is.na(UICC8_edition3)), aes(x = as.factor(quantile_pseudotime), y = average_malignancy_rank, fill = cancer_type)) +
  geom_boxplot(width=0.5, alpha=0.8, show.legend=F) +
  #facet_wrap(~UICC8_edition3, nrow=2)+
  facet_grid(cancer_type~UICC8_edition3)+
  scale_fill_manual(values = c('#ff7700', '#1414ad')) +
  theme_minimal() +
  theme(text=element_text(size=25), 
    legend.title = element_blank()) +
    scale_fill_manual(values = c("#11B5AE", "#4046C9")) + 
  xlab('Pseudotime (quantiles)') +
  ylab('Average metastatic potential') 
  dev.off()

gradient_colors <-  c("#11B5AE", "#4046C9", "#F68512", "#DE3C82") # plasma(5)[1:4]
