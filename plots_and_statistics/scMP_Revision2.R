library(data.table)
library(dplyr)
library(magrittr)
library(ggplot2)
library(tidyverse)
library(Hmisc)
library(glmnet)
library(arrow)
library(survival)
library(survminer)
library(survcomp)
library(Seurat)
library(ggsankey)
library(readxl)
library(pROC)
library(viridis)
library(cowplot)

getwd()
setwd('plots_statistics')
###CAVE biggest change to revision 1: using k=4

match_ids <- read.csv('/mnt/ssd/shared/LungCAIRE_SingleCell/data/LungCAIRE_snRNAseq_cohort.csv') %>% 
  dplyr::select(Match_ID, Pseudo) %>% mutate(ENR = Match_ID) #%>% 
  #filter(Pseudo != 'LC395')

patient_meta_data <- read.csv('/mnt/ssd/shared/LungCAIRE_SingleCell/data/Metadata_LungCAIRE_retrospective_Berlin_Cologne_20250303.csv') %>% 
  left_join(match_ids) %>%
  mutate(UICC8_edition3 = ifelse(UICC8_EDITION2 %in% c('IA', 'IB'), 'I', UICC8_EDITION2)) %>% 
  mutate(Entity = ENTITY, Event = EVENT, OS_m = OS_M, Grade = GRADE, Adj_therapy = ADJ_THERAPY) %>% 
  #filter(Pseudo != 'LC395') %>% 
  filter(ENTITY != 'ASC') 

#######################
#load clean annotation and replace old annotation
clean_cell_annotation <- read.csv('./use_data/VICREG_clustering_highinputdropout_5000epochs_rev_annotation.csv') %>% 
  dplyr::select(cluster = Cluster_ID, Tier_1, Tier_2, Tier_3)

cluster_result <- read_parquet('../scGPT_embeddings/clusterings/VICREG_clustering_highinputdropout_5000epochs.parquet') %>%
  mutate(cell_position = seq(dim(.)[1]))   %>% 
  dplyr::select(-c(Tier_1, Tier_2, Tier_3, Tier_4)) %>% 
  left_join(clean_cell_annotation)
########################

annotated_clusters_FINAL <- read.csv('./use_data/VICREG_clustering_highinputdropout_5000epochs_rev_annotation.csv') %>%
  mutate(Tier_2 = ifelse(Tier_3 == 'Tumor', 'Tumor', Tier_2))
clusters_not_mixed <- annotated_clusters_FINAL  %>% filter(Tier_3 != 'Mixed')

cell_types_fine_and_coarse <- annotated_clusters_FINAL %>% 
  mutate(max_fine_Tier = Tier_3, max_coarse_Tier = Tier_2) %>%
  mutate(cluster = paste0('X', Cluster_ID)) %>% 
  dplyr::select(-Cluster_ID)

patient_vector_long <- cluster_result %>% 
  group_by(cluster, Pseudo) %>%
  dplyr::summarize(N=n(), cluster = first(cluster))  %>% 
  group_by(Pseudo) %>% 
  mutate(proportion = N/sum(N)) %>% 
  dplyr::select(cluster, proportion, Pseudo, N) %>% 
  ungroup() %>% 
  filter(cluster %in% clusters_not_mixed$Cluster_ID)
  
patient_vector <- patient_vector_long %>%
  dplyr::select(Pseudo, cluster, proportion) %>%
  group_by(cluster) %>% 
  pivot_wider(names_from = cluster, values_from = proportion)

patient_vector[is.na(patient_vector)] <- 0

patients_with_snRNA <- patient_vector_long$Pseudo %>% unique()

filtered_patient_meta_data <- patient_meta_data %>% 
  filter(Pseudo %in% patients_with_snRNA)


##########################################################
##presence of any superkiller tumor cells?
###########################################################
extract_largest_number <- Vectorize(function(input_string) {
  # Extract all numbers (integers or decimals) from the string
  numbers <- as.numeric(unlist(regmatches(input_string, gregexpr("\\d+\\.?\\d*", input_string))))
  
  # Check if any numbers were found
  if (length(numbers) == 0) {
    return(NA)  # Return NA if no numbers are found
  }
  # Return the largest number  #NOW AVERAGE IS USED!
  return(mean(numbers, na.rm = TRUE))
})


tumor_diameter <- filtered_patient_meta_data %>%
  dplyr::select(Pseudo, DIAMETER_MM) %>% 
  mutate(DIAMETER_MM = as.numeric(extract_largest_number(DIAMETER_MM))) %>%
  filter(!is.na(DIAMETER_MM)) %>%
  group_by(Pseudo) %>%
  dplyr::summarize(DIAMETER_MM = mean(DIAMETER_MM), N=n()) %>% 
  dplyr::mutate(volume = 4/3 * pi * (DIAMETER_MM /2)**3) %>%
  mutate(centered_diameter = volume) %>%
  dplyr::select(Pseudo, centered_diameter)

metastasized_yes_no <- filtered_patient_meta_data %>% 
  dplyr::select(ENR, PTNM_N = PTNM_N_SUBGROUPS, PTNM_M) %>% unique() %>%
  inner_join(filtered_patient_meta_data %>% dplyr::select(Pseudo, ENR = Match_ID) %>% unique()) %>%
  dplyr::select(Pseudo, PTNM_N, PTNM_M) %>% 
  mutate(numeric_N = as.numeric(factor(PTNM_N, levels = c('0', '1a','1b', '2a1', '2a2', '2b', '3'))), 
          numeric_M = as.numeric(factor(PTNM_M, levels = c('0', '1a','1b', '1c'))))%>% 
  filter(!is.na(PTNM_N), !is.na(PTNM_M)) %>%
  dplyr::mutate(PTNM_N= PTNM_N!=0, PTNM_M = PTNM_M !=0) %>%
  group_by(Pseudo) %>% 
  dplyr::summarize(PTNM_N = mean(PTNM_N), PTNM_M = mean(PTNM_M),
                numeric_N = mean(numeric_N), numeric_M = mean(numeric_M)) %>% 
  mutate(malign = (PTNM_N==1 | PTNM_M==1)*1)


relapse <-patient_meta_data %>% 
  dplyr::select(ENR, TIME_TO_RELAPS_M, RELAPS_TYPE, DIA_DATE, LAST_CONTACT, R_STATUS, R_TYPE) %>% unique() %>%
  mutate(DIA_DATE = as.Date(DIA_DATE, format = "%d.%m.%Y"), LAST_CONTACT = as.Date(LAST_CONTACT, format = "%d.%m.%Y")) %>%
  mutate(TIME_DIAGNOSIS2LASTCONTACT = as.numeric(interval(DIA_DATE, LAST_CONTACT) / months(1))) %>%
  inner_join(filtered_patient_meta_data %>% dplyr::select(Pseudo, ENR = Match_ID, OS_m, Event) %>% unique()) %>%
  dplyr::select(Pseudo, TIME_TO_RELAPS_M, RELAPS_TYPE, TIME_DIAGNOSIS2LASTCONTACT, R_STATUS, R_TYPE) %>% 
  unique()

relapse %>% dplyr::select(R_STATUS, R_TYPE) %>%group_by(R_STATUS, R_TYPE) %>% dplyr::summarize(N=n())


######
##use SCM to determine influence of cell clusters on metastasis
######
get_coefs_logistic <- function(outcome, a, b, g) {
  df <- data.frame(outcome,a,b)
  logistic_model <- glm(outcome ~ a + b, data = df, family = "binomial")
  
  log_reg_res <- summary(logistic_model)
  data.frame('r' = log_reg_res$coefficients[2,1] )$r
  }


UICC_correlating_cell_clusters_preparation <- patient_vector %>%
  pivot_longer(!Pseudo, names_to = 'cluster', values_to = 'proportion') %>% 
  left_join(metastasized_yes_no) %>% 
  left_join(patient_meta_data %>% dplyr::select(Pseudo, Entity)) %>% 
  inner_join(tumor_diameter) %>%
  left_join(annotated_clusters_FINAL %>% dplyr::select(cluster = Cluster_ID, Tier_3, Tier_2) %>% mutate(cluster = paste0('', cluster))) %>% 
  filter(Tier_2 == 'Tumor') %>% 
  group_by(Pseudo, Entity, malign) %>% 
  mutate(proportion_within_tumorcells = proportion / sum(proportion)) 

potential_latent_metastatic_cells <- UICC_correlating_cell_clusters_preparation %>% 
  group_by(cluster) %>%
  filter(!is.na(malign)) %>%
  dplyr::summarize(r = get_coefs_logistic(malign, ntile(proportion_within_tumorcells,100), centered_diameter))


malignant_ranking <- potential_latent_metastatic_cells %>% 
  dplyr::select(r, cluster) %>% 
  unique() %>%
  #filter(r>0) %>% 
  #mutate(r = ifelse(r>0, r, 0)) %>% 
  #dplyr::mutate(malignancy_rank = rank(r, ties.method = 'max')) %>% 
  dplyr::mutate(malignancy_rank = rank(r)) %>% 
  ungroup() %>%
  mutate(malignancy_rank = malignancy_rank / max(malignancy_rank) * 100 ) %>%
  mutate(cluster = paste0('X', cluster)) 

#write.csv(malignant_ranking, './use_data/latent_metastatic_cells_candidates.csv')



get_kth_max <- function(x){
  if (length(x)==0) {return(0)}
  x = sort(x, decreasing=T)
  k <- 4
  res <- ifelse(length(x)>=k, x[k], min(x))
}



patients_with_malignancy_ranking <- patient_vector_long %>% 
  mutate(cluster = paste0('X', cluster)) %>%
  left_join(cell_types_fine_and_coarse) %>% 
  filter(Tier_2 == 'Tumor') %>% 
  group_by(Pseudo) %>%
  mutate(proportion_within_tumorcells = proportion/sum(proportion)) %>%
  left_join(patient_meta_data %>% dplyr::select(Pseudo, UICC8_edition3, OS_m, Event, Entity)) %>%
  inner_join(malignant_ranking) %>%
  #filter(proportion_within_tumorcells>=0.005) %>% 
  filter(proportion_within_tumorcells>0.002, N>=2) %>% #0.001
  group_by(Pseudo, UICC8_edition3, Event, OS_m) %>% 
  mutate(kth_max = get_kth_max(malignancy_rank)) %>%
  ungroup() %>%
  filter(malignancy_rank >= kth_max) %>%
  dplyr::group_by(Pseudo, UICC8_edition3, Event, OS_m) %>%   #calculate malignancy rank for AC and SCC together
  dplyr::summarize(malignancy_rank = mean(malignancy_rank)) %>% #weigh malignancy rank by proportion?
  ungroup() %>%
  left_join(patient_meta_data %>% dplyr::select(Pseudo, Grade, Entity, Adj_therapy,  
                                                packyears = PACK.YEARS, 
                                                NEOADJ_THERAPY, 
                                                RAD_THERAPY, 
                                                R_STATUS, R_TYPE, 
                                                ALK, MET, P63, pdl1 = `PD.L1....`, ki67 = `ki67....`)) %>%
  mutate(grade_and_MPR = Grade *malignancy_rank) %>% 
    left_join(metastasized_yes_no %>% dplyr::select(Pseudo, PTNM_N, PTNM_M, malign)) %>% 
  dplyr::mutate(KM_risk_group = ifelse(Entity == 'AC', 
                                    case_when(#malignancy_rank>=99.1 ~'r4', 
                                      malignancy_rank<210 & malignancy_rank>=93.38 ~'r2', 
                                      #malignancy_rank<=85.7 & malignancy_rank>72.5~'r2', 
                                      malignancy_rank<93.38 ~'r1'), 
                                      case_when(#malignancy_rank>=99.1 ~'r4', 
                                      malignancy_rank<210 & malignancy_rank>=93.93~'r2', 
                                      #malignancy_rank<=95.6 & malignancy_rank>92.0~'r2', 
                                      malignancy_rank<93.93 ~'r1'))) 


patients_with_malignancy_ranking %>% 
  group_by(Entity, KM_risk_group) %>% 
  dplyr::summarize(N= n())

colnames(patient_meta_data)  

#write.csv(patients_with_malignancy_ranking, './use_data/patient_metastatic_potential.csv')




#correlation between grading and malignancy?
patients_with_malignancy_ranking %>% 
  group_by(Entity) %>% 
  filter(!is.na(Grade)) %>%
  dplyr::summarize(r_tograde = rcorr(malignancy_rank, Grade, type = 'spearman')$r[1,2], 
  P = rcorr(malignancy_rank, Grade, type = 'spearman')$P[1,2]) 


get_c_index <- function(pred, time, event) {
    concordance.index(
    x = pred,  # The predictor variable
    surv.time = time,  # The survival time
    surv.event = event, 
    method = 'noether'  # Concordance method (optional, "noether" or "unbiased")
    )$c.index
  }

get_c_index_p <- function(pred, time, event) {
    res <- concordance.index(
    x = pred,  # The predictor variable
    surv.time = time,  # The survival time
    surv.event = event, 
    method = 'noether'  # Concordance method (optional, "noether" or "unbiased")
    )
    res$p.value
  }


patients_with_malignancy_ranking %>% 
  ungroup() %>%   
  group_by(Entity, UICC8_edition3) %>% 
  dplyr::summarize(c_index = get_c_index(malignancy_rank, OS_m, Event),
                  P = get_c_index_p(malignancy_rank, OS_m, Event),
                  N=n()) 

patients_with_malignancy_ranking %>% 
  ungroup() %>%   
  group_by(UICC8_edition3) %>% 
  dplyr::summarize(c_index = get_c_index(malignancy_rank, OS_m, Event),
                  P = get_c_index_p(malignancy_rank, OS_m, Event),
                  N=n()) 

patients_with_malignancy_ranking$Adj_therapy %>% unique()



################
#show distribution of first k scMP per patient
#################
get_some_max <- function(x){
  if (length(x)==0) {return(0)}
  x = sort(x, decreasing=T)
  k <- 4
  res <- ifelse(length(x)>=k, x[k], min(x))
}

highest_scMP_per_patient <- patient_vector_long %>% 
  mutate(cluster = paste0('X', cluster)) %>%
  left_join(cell_types_fine_and_coarse) %>% 
  filter(Tier_2 == 'Tumor') %>% 
  group_by(Pseudo) %>%
  mutate(proportion_within_tumorcells = proportion/sum(proportion)) %>%
  left_join(patient_meta_data %>% dplyr::select(Pseudo, UICC8_edition3, OS_m, Event, Entity)) %>%
  left_join(malignant_ranking) %>%
  filter(proportion_within_tumorcells>=0.005) %>% 
  group_by(Pseudo, UICC8_edition3, Event, OS_m) %>% 
  mutate(kth_max = get_some_max(malignancy_rank)) %>%
  ungroup() %>%
  filter(malignancy_rank >= kth_max) %>% 
  group_by(Pseudo) %>% 
  mutate(rankl = ((rank(desc(malignancy_rank))))) %>% 
  filter(rankl!=1.5, rankl != 3.5)

pseudo_order <- highest_scMP_per_patient %>% 
  filter(rankl == 1) %>% 
  arrange(malignancy_rank)

highest_scMP_per_patient$Pseudo <- factor(highest_scMP_per_patient$Pseudo, levels = pseudo_order$Pseudo)
#ggplot(highest_scMP_per_patient, aes(x = Pseudo, y = malignancy_rank, color = rankl, group = rankl, size=N)) +
#  geom_point() +
#  scale_color_viridis(option = 'C', begin = 0.0, end = 0.9) +
#  coord_flip()
####################


#GRADING!
patients_with_malignancy_ranking %>%   
  filter(!is.na(Grade)) %>%
  ungroup() %>%
  group_by(Entity, UICC8_edition3) %>% 
  #group_by(UICC8_edition3) %>% 
  dplyr::summarize(c_index = get_c_index(Grade, OS_m, Event),
                  P = get_c_index_p(Grade, OS_m, Event), 
                c_index_withMPR = get_c_index(grade_and_MPR, OS_m, Event),
                  P_withMPR = get_c_index_p(grade_and_MPR, OS_m, Event), 
                  N=n()) 

#KI67!
patients_with_malignancy_ranking %>% 
  mutate(ki67 = as.numeric(ki67)) %>%  
  filter(!is.na(ki67)) %>%
  ungroup() %>%
  group_by(Entity, UICC8_edition3) %>% 
  #group_by(UICC8_edition3) %>% 
  dplyr::summarize(c_index = get_c_index(ki67, OS_m, Event),
                  P = get_c_index_p(ki67, OS_m, Event), 
                  N=n()) 


patients_with_malignancy_ranking %>% 
  mutate(ki67 = as.numeric(ki67)) %>% 
  filter(!is.na(ki67)) %>% 
  ungroup() %>% 
  dplyr::summarize(N=n())


nquantiles <- 4
malignancy_ranking_for_KM <- patients_with_malignancy_ranking %>% 
  ungroup() %>%
  #dplyr::mutate(KM_risk_group = ntile(malignancy_rank,nquantiles)) %>%
    mutate(Entity_adj = paste0(Entity, Adj_therapy ))

gradient_colors <- scales::gradient_n_pal(c("#049404",'#ffc400', "red"))(seq(0, 1, length.out = nquantiles))

png('./figures/Survival_internal_k4.png', width=2000, height=1500, res=250)
ggsurvplot_facet(
  surv_fit_malignancy_ranking <- survfit(Surv(OS_m, Event) ~ KM_risk_group, data = malignancy_ranking_for_KM),
  data = malignancy_ranking_for_KM %>% filter(UICC8_edition3 == 'I'),
  #facet.by = c('UICC8_edition3'),
  facet.by = c('UICC8_edition3', 'Entity'),
  pval = F,           # Add p-value
  pval.method = TRUE,    # Show test method in the p-value plot
  pval.size = 4,         # Set the p-value font size
  conf.int = F,       # Add confidence intervals
  censor = TRUE,         # Show censored data points
  legend.title = "Risk quantile",
  title = 'Survival Internal Dataset',
  xlab = "Time",
  ylab = "Survival Probability",
  palette = gradient_colors,
  linewidth=0.1,
  xlim = c(0,60),
  ggtheme=theme_minimal(),
  text=element_text(size=20))
dev.off()

survival_internal_all_stages_plot <- ggsurvplot_facet(
  surv_fit_malignancy_ranking <- survfit(Surv(OS_m, Event) ~ KM_risk_group, data = malignancy_ranking_for_KM),
  #data = malignancy_ranking_for_KM %>% filter(UICC8_edition3 == 'I'),
  data = malignancy_ranking_for_KM %>% mutate(UICC = UICC8_edition3, 
                                                KM_risk_group = case_when( KM_risk_group == 'r1'~'RGL', 
                                                KM_risk_group == 'r2'~'RGH', 
                                                KM_risk_group == 'r3'~'NOTEXISTENT', 
                                                ), Entity = ifelse(Entity == 'AC', 'LUAD', 'LUSC')),
  #facet.by = c('UICC8_edition3'),
  facet.by = c('UICC', 'Entity'),
  pval = F,           # Add p-value
  pval.method = TRUE,    # Show test method in the p-value plot
  pval.size = 4,         # Set the p-value font size
  conf.int = F,       # Add confidence intervals
  censor = TRUE,         # Show censored data points
  legend.title = "Risk group",
  title = 'Overall Survival',
  xlab = "Time (Months)",
  ylab = "Survival Probability",
  #palette = gradient_colors,
  palette = c("#13B5AE",'#3F46C9', "#F68512"),
  linewidth=0.1,
  xlim = c(0,60),
  ggtheme=theme_minimal(base_size=15))

png('./figures/Survival_internal_all_stages_k4.png', width=2000, height=2500, res=250)
survival_internal_all_stages_plot
dev.off()

pdf('./figures/Survival_internal_all_stages_k4.pdf', width=8, height=10,onefile = F)
survival_internal_all_stages_plot
dev.off()

patients_with_malignancy_ranking_and_relapse <- patients_with_malignancy_ranking %>% 
  left_join(relapse) %>% 
  mutate(inverse_TIME_TO_RELAPS_M = as.numeric(ifelse(TIME_TO_RELAPS_M == 'N/APP', 0, 1/as.numeric(TIME_TO_RELAPS_M)))) %>% 
  mutate(numeric_TIME_TO_RELAPS_M = as.numeric(ifelse(TIME_TO_RELAPS_M == 'N/APP', Inf, as.numeric(TIME_TO_RELAPS_M)))) %>% 
  mutate(censored_TIME_TO_RELAPS_M = as.numeric(ifelse(TIME_TO_RELAPS_M == 'N/APP', TIME_DIAGNOSIS2LASTCONTACT, as.numeric(TIME_TO_RELAPS_M)))) %>% 
  filter(! (R_STATUS==1 & RELAPS_TYPE == 'locoregional') ) %>%
  mutate(relapse_occurred_within_2_years = (numeric_TIME_TO_RELAPS_M<=24))  %>% 
  mutate(relapse_event = ifelse(TIME_TO_RELAPS_M == 'N/APP',0,1)) %>%
  mutate(NM_status = ifelse(PTNM_M == 0, ifelse(PTNM_N == 0, 'N0M0', 'N+M0'), 'M+')) %>% 
  mutate(NM_status = factor(NM_status, levels = c('N0M0', 'N+M0', 'M+'))) %>%
  filter(!is.na(Entity))   #why is this necessary?! 


#write.csv(patients_with_malignancy_ranking_and_relapse %>% dplyr::select(Pseudo, censored_TIME_TO_RELAPS_M, relapse_event), './use_data/time2relapse.csv')

ggplot(patients_with_malignancy_ranking_and_relapse, aes(x = KM_risk_group,  fill = Adj_therapy)) +
  #geom_bar(position  = 'fill') +
  geom_bar() +
  facet_grid(UICC8_edition3 ~ Entity)

patients_with_malignancy_ranking_and_relapse_AC<- patients_with_malignancy_ranking_and_relapse %>% filter(Entity=='AC', UICC8_edition3 == 'I')
adj_therapy_contingency_table <- table(patients_with_malignancy_ranking_and_relapse_AC$KM_risk_group, patients_with_malignancy_ranking_and_relapse_AC$Adj_therapy)
adj_therapy_contingency_table
fisher.test(adj_therapy_contingency_table)

#####################
#risk groups vs sub staging
#####################
predata_for_sankey_from_substaging <- patients_with_malignancy_ranking_and_relapse %>% 
  filter(Entity == 'AC', UICC8_edition3 == 'I') %>% 
  left_join(patient_meta_data %>% dplyr::select(Pseudo, UICC8_EDITION2)) %>% 
  dplyr::select(Pseudo, Treatment_adj = Adj_therapy, Substaging = UICC8_EDITION2, , KM_risk_group)
  
data_for_sankey_from_substaging <- predata_for_sankey_from_substaging%>% 
  make_long(!c(Pseudo))

ggplot(data_for_sankey_from_substaging, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = node)) +
  geom_sankey(flow.alpha = 0.6, node.color = "black") +
  theme_sankey(base_size = 14) +
  labs(title = "Sankey Diagram Example", x = "Stage", y = "Count") +
  theme_minimal()

###########
#calculate cutoffs for risk groups
############
cutoff_for_no_relapse_within_5_years <- patients_with_malignancy_ranking_and_relapse %>% 
  #filter(TIME_TO_RELAPS_M<=60) %>% 
  filter(relapse_event==1, UICC8_edition3 %in% c('I')) %>%
  #filter(Entity == 'AC') %>%
  .$malignancy_rank %>% min()
cutoff_for_no_relapse_within_5_years

cutoff_for_middle_groups <- patients_with_malignancy_ranking_and_relapse %>%   
  #filter(malignancy_rank>cutoff_for_no_relapse_within_5_years) %>%
  dplyr::summarize(middle_value = (100+min(malignancy_rank))/2)

cutoff_for_middle_groups

cutoff_for_no_relapse_within_5_years_SCC <- patients_with_malignancy_ranking_and_relapse_all %>% 
  filter(relapse_event==1, UICC8_edition3 %in% c('I')) %>%
  filter(Entity == 'SCC') %>%
  .$malignancy_rank %>% min()
cutoff_for_no_relapse_within_5_years

cutoff_for_middle_groups_SCC <- patients_with_malignancy_ranking_and_relapse_all %>%   
  filter(Entity == 'SCC') %>%
  #filter(malignancy_rank>cutoff_for_no_relapse_within_5_years_SCC) %>%
  dplyr::summarize(middle_value = (100+min(malignancy_rank))/2)
cutoff_for_middle_groups_SCC

###mean and median cutoff
patients_with_malignancy_ranking_and_relapse %>% 
  filter(UICC8_edition3 == 'I') %>% 
  group_by(Entity) %>% 
  mutate(maxv = quantile(malignancy_rank, 0.98), minv = quantile(malignancy_rank, 0.02)) %>% 
  filter(malignancy_rank < maxv, malignancy_rank> minv) %>% 
  dplyr::summarize(meanv = mean(malignancy_rank), medianv = median(malignancy_rank), 
    q1 = quantile(malignancy_rank, 0.25), q2 = quantile(malignancy_rank, 0.75))

patients_with_malignancy_ranking_and_relapse %>% 
  filter(UICC8_edition3 == 'I') %>% 
  group_by(Entity) %>% 
  dplyr::summarize(maxv = quantile(malignancy_rank, 0.98), minv = quantile(malignancy_rank, 0.02)) %>% 
    mutate(q1 = minv + 0.33 * (maxv-minv), q2 = minv + 0.66 * (maxv-minv), middlev = minv + 0.5 * (maxv-minv))

ggplot(patients_with_malignancy_ranking_and_relapse  %>%  filter(UICC8_edition3 == 'I'), aes(x = malignancy_rank)) +
  geom_histogram() +
  facet_wrap(~Entity, ncol=1)

library(survminer)
surv_cutpoint(patients_with_malignancy_ranking_and_relapse %>% filter(Entity == 'AC', UICC8_edition3 == 'I'), time = 'OS_m', event = 'Event', variables = c('malignancy_rank'), minprop = 0.25) 
surv_cutpoint(patients_with_malignancy_ranking_and_relapse %>% filter(Entity == 'SCC', UICC8_edition3 == 'I'), time = 'OS_m', event = 'Event', variables = c('malignancy_rank'), minprop=0.25) 



#now do KM to relapse
a <- patients_with_malignancy_ranking_and_relapse %>% 
  dplyr::select(Pseudo, OS_m, Event, censored_TIME_TO_RELAPS_M, relapse_event) %>% 
  mutate(OS_m >= censored_TIME_TO_RELAPS_M)

patients_with_malignancy_ranking_and_relapse$KM_risk_group






# Create truncated survival variables: any OS_m > 60 becomes 60, and the event is set to 0 (censored)
malignancy_ranking_for_KM_60 <- malignancy_ranking_for_KM %>%
  mutate(OS_m_trunc = pmin(OS_m, 120),
         Event_trunc = if_else(OS_m > 120, 0, Event))

malignancy_ranking_for_KM_2 = malignancy_ranking_for_KM_60 %>% filter(UICC8_edition3 == 'I')
malignancy_ranking_for_KM_2_AC = malignancy_ranking_for_KM_2 %>% filter(Entity == 'AC')
malignancy_ranking_for_KM_2_SCC = malignancy_ranking_for_KM_2 %>% filter(Entity == 'SCC')


#gradient_colors <- scales::gradient_n_pal(c("#13B5AE",'#3F46C9', "#DE3C82"))(seq(0, 1, length.out = nquantiles))

pdf('./figures/Survival_internal_AC_k4.pdf', width=5, height=6,onefile = F)
plotfile <- ggsurvplot(
  surv_fit_malignancy_ranking <- survfit(Surv(OS_m_trunc, Event) ~ KM_risk_group, data = malignancy_ranking_for_KM_2_AC),
  data = malignancy_ranking_for_KM_2_AC,
  #facet.by = c('UICC8_edition3'),
  #facet.by = c('Entity'),
  pval = T,           # Add p-value
  pval.method = TRUE,    # Show test method in the p-value plot
  pval.size = 4,         # Set the p-value font size
  conf.int = F,       # Add confidence intervals
  censor = TRUE,         # Show censored data points
  legend.title = "Risk quantile",
  title = 'Survival Internal Dataset',
  xlab = "Time",
  ylab = "Survival Probability",
  break.time.by = 12,
  risk.table = T,
  #palette = c("#13B5AE",'#3F46C9', "#F68512"),
  palette = c('#3F46C9', "#F68512"),
  linewidth=0.1,
  #xlim = c(0,60),
  ggtheme=theme_minimal(),
  text=element_text(size=20))

plotfile$table <- plotfile$table + 
  scale_y_discrete(labels = c("RG3", "RG2", "RG1")) # Customize y-axis labels

plotfile
dev.off()

pdf('./figures/Survival_internal_SCC_k4.pdf', width=5, height=6,onefile = F)
plotfile <- ggsurvplot(
  surv_fit_malignancy_ranking <- survfit(Surv(OS_m_trunc, Event) ~ KM_risk_group, data = malignancy_ranking_for_KM_2_SCC),
  data = malignancy_ranking_for_KM_2_SCC,
  #facet.by = c('UICC8_edition3'),
  #facet.by = c('Entity'),
  pval = T,           # Add p-value
  pval.method = TRUE,    # Show test method in the p-value plot
  pval.size = 4,         # Set the p-value font size
  conf.int = F,       # Add confidence intervals
  censor = TRUE,         # Show censored data points
  legend.title = "Risk quantile",
  title = 'Survival Internal Dataset',
  xlab = "Time",
  ylab = "Survival Probability",
  break.time.by = 12,
  risk.table = T,
  #palette = c("#13B5AE",'#3F46C9', "#F68512"),
  palette = c('#3F46C9', "#F68512"),
  linewidth=0.1,
  #xlim = c(0,60),
  ggtheme=theme_minimal(),
  text=element_text(size=20))

plotfile$table <- plotfile$table + 
  scale_y_discrete(labels = c("RG3", "RG2", "RG1")) # Customize y-axis labels

plotfile
dev.off()



patients_with_malignancy_ranking_and_relapse_all <- patients_with_malignancy_ranking %>% 
  left_join(relapse) %>% 
  mutate(inverse_TIME_TO_RELAPS_M = as.numeric(ifelse(TIME_TO_RELAPS_M == 'N/APP', 0, 1/as.numeric(TIME_TO_RELAPS_M)))) %>% 
  mutate(numeric_TIME_TO_RELAPS_M = as.numeric(ifelse(TIME_TO_RELAPS_M == 'N/APP', Inf, as.numeric(TIME_TO_RELAPS_M)))) %>% 
  mutate(censored_TIME_TO_RELAPS_M = as.numeric(ifelse(TIME_TO_RELAPS_M == 'N/APP', TIME_DIAGNOSIS2LASTCONTACT, as.numeric(TIME_TO_RELAPS_M)))) %>% 
  filter(! (R_STATUS==1 & RELAPS_TYPE == 'locoregional') ) %>%
  mutate(relapse_occurred_within_2_years = (numeric_TIME_TO_RELAPS_M<=24))  %>% 
  mutate(relapse_event = ifelse(TIME_TO_RELAPS_M == 'N/APP',0,1)) %>%
  mutate(NM_status = ifelse(PTNM_M == 0, ifelse(PTNM_N == 0, 'N0M0', 'N+M0'), 'M+')) %>% 
  mutate(NM_status = factor(NM_status, levels = c('N0M0', 'N+M0', 'M+'))) %>%
  filter(!is.na(Entity))   #why is this necessary?! 





# Create truncated survival variables: any OS_m > 60 becomes 60, and the event is set to 0 (censored)
patients_with_malignancy_ranking_and_relapse_60 <- patients_with_malignancy_ranking_and_relapse_all %>%
  mutate(censored_TIME_TO_RELAPS_M = pmin(censored_TIME_TO_RELAPS_M, 60),
         relapse_event = if_else(censored_TIME_TO_RELAPS_M > 60, 0, relapse_event))

patients_with_malignancy_ranking_and_relapse_60 = patients_with_malignancy_ranking_and_relapse_60 %>% filter(UICC8_edition3 == 'I')
patients_with_malignancy_ranking_and_relapse_60_AC = patients_with_malignancy_ranking_and_relapse_60 %>% filter(Entity == 'AC')
patients_with_malignancy_ranking_and_relapse_60_SCC = patients_with_malignancy_ranking_and_relapse_60 %>% filter(Entity == 'SCC')


#gradient_colors <- scales::gradient_n_pal(c("#13B5AE",'#3F46C9', "#DE3C82"))(seq(0, 1, length.out = nquantiles))

pdf('./figures/Relapse_internal_AC_k4.pdf', width=5, height=6,onefile = F)
plotfile <- ggsurvplot(
  surv_fit_malignancy_ranking <- survfit(Surv(censored_TIME_TO_RELAPS_M, relapse_event) ~ KM_risk_group, data = patients_with_malignancy_ranking_and_relapse_60_AC),
  data = patients_with_malignancy_ranking_and_relapse_60_AC,
  #facet.by = c('UICC8_edition3'),
  #facet.by = c('Entity'),
  pval = T,           # Add p-value
  pval.method = TRUE,    # Show test method in the p-value plot
  pval.size = 4,         # Set the p-value font size
  conf.int = F,       # Add confidence intervals
  censor = TRUE,         # Show censored data points
  legend.title = "Risk quantile",
  title = 'Survival Internal Dataset',
  xlab = "Time",
  ylab = "Survival Probability",
  break.time.by = 12,
  risk.table = T,
  palette = c("#13B5AE",'#3F46C9', "#F68512"),
  linewidth=0.1,
  #xlim = c(0,60),
  ggtheme=theme_minimal(),
  text=element_text(size=20))

plotfile$table <- plotfile$table + 
  scale_y_discrete(labels = c("RG3", "RG2", "RG1")) # Customize y-axis labels

plotfile
dev.off()


summary_5yr_internal <- summary(surv_fit_malignancy_ranking, times = 60)
summary_5yr_internal


pdf('./figures/Relapse_internal_SCC_k4.pdf', width=5, height=6,onefile = F)
plotfile <- ggsurvplot(
  surv_fit_malignancy_ranking <- survfit(Surv(censored_TIME_TO_RELAPS_M, relapse_event) ~ KM_risk_group, data = patients_with_malignancy_ranking_and_relapse_60_SCC),
  data = patients_with_malignancy_ranking_and_relapse_60_SCC,
  #facet.by = c('UICC8_edition3'),
  #facet.by = c('Entity'),
  pval = T,           # Add p-value
  pval.method = TRUE,    # Show test method in the p-value plot
  pval.size = 4,         # Set the p-value font size
  conf.int = F,       # Add confidence intervals
  censor = TRUE,         # Show censored data points
  legend.title = "Risk quantile",
  title = 'Survival Internal Dataset',
  xlab = "Time",
  ylab = "Survival Probability",
  break.time.by = 12,
  risk.table = T,
  palette = c("#13B5AE",'#3F46C9', "#F68512"),
  linewidth=0.1,
  #xlim = c(0,60),
  ggtheme=theme_minimal(),
  text=element_text(size=20))

plotfile$table <- plotfile$table + 
  scale_y_discrete(labels = c("RG3", "RG2", "RG1")) # Customize y-axis labels

plotfile
dev.off()


gradient_colors_relapse <- (scales::gradient_n_pal(c("#049404",'#ffc400', "red"))(seq(0, 1, length.out =patients_with_malignancy_ranking_and_relapse$KM_risk_group %>% unique() %>% length() )))


relapse_internal_all_patients_plot <- ggsurvplot_facet(
  surv_fit_malignancy_ranking <- survfit(Surv(censored_TIME_TO_RELAPS_M, relapse_event) ~ KM_risk_group, data = patients_with_malignancy_ranking_and_relapse),
  #data = patients_with_malignancy_ranking_and_relapse %>% filter(UICC8_edition3 %in% c('I')),
  data = patients_with_malignancy_ranking_and_relapse_all %>% mutate(UICC = UICC8_edition3, 
                                                KM_risk_group = case_when( KM_risk_group == 'r1'~'RGL', 
                                                KM_risk_group == 'r2'~'RGH', 
                                                KM_risk_group == 'r3'~'NOTEXISTENT', 
                                                ), Entity = ifelse(Entity == 'AC', 'LUAD', 'LUSC')),
  #facet.by = c('UICC8_edition3'),
  #facet.by = c('Entity'),
  facet.by = c('UICC', 'Entity'),
  pval = F,           # Add p-value
  pval.method = TRUE,    # Show test method in the p-value plot
  pval.size = 4,         # Set the p-value font size
  conf.int = F,       # Add confidence intervals
  censor = TRUE,         # Show censored data points
  legend.title = "Risk group",
  title = 'Time-to-relapse',
  xlab = "Time (Months)",
  ylab = "Relapse-free Probability",
  #palette = gradient_colors,
  palette = c("#13B5AE",'#3F46C9', "#F68512"),
  linewidth=0.1,
  xlim = c(0,60),
  ggtheme=theme_minimal(base_size=15))



png('./figures/Relapse_internal_all_stages_k4.png', width=2000, height=2500, res=250)
relapse_internal_all_patients_plot
dev.off()

pdf('./figures/Relapse_internal_all_stages_k4.pdf', width=8, height=10,onefile = F)
relapse_internal_all_patients_plot
dev.off()


ggsurvplot(
  surv_fit_malignancy_ranking <- survfit(Surv(censored_TIME_TO_RELAPS_M, relapse_event) ~ KM_risk_group, data = patients_with_malignancy_ranking_and_relapse),
  data = patients_with_malignancy_ranking_and_relapse,,# %>% filter(UICC8_edition3 %in% c('I')),
  #facet.by = c('UICC8_edition3'),
  facet.by = c( 'UICC8_edition3', 'Entity'),
  pval = F,           # Add p-value
  pval.method = TRUE,    # Show test method in the p-value plot
  pval.size = 4,         # Set the p-value font size
  conf.int = F,       # Add confidence intervals
  censor = TRUE,         # Show censored data points
  legend.title = "Risk quantile",
  title = 'Time to relapse',
  xlab = "Time",
  ylab = "Probability",
  palette = gradient_colors_relapse,
  linewidth=0.1,
  xlim = c(0,60),
  risk.table=T,
  ncol=2,
  ggtheme=theme_minimal(),
  text=element_text(size=20))

#how many patients relapse in RG1-RG3?
patients_with_malignancy_ranking_and_relapse %>% 
  filter(UICC8_edition3 == 'I') %>%
  #filter(censored_TIME_TO_RELAPS_M<=60) %>%
  group_by(Entity, KM_risk_group) %>% 
    dplyr::summarize(N_relapses = sum(relapse_event), N_overall = n()) %>% mutate(percentage = N_relapses/N_overall)

#hazard ratios
hazard_ratio_data_AC <- patients_with_malignancy_ranking_and_relapse %>% filter(Entity == 'AC', UICC8_edition3=='I') 
hazard_ratios_cox_model_relapse_AC <- coxph(Surv(censored_TIME_TO_RELAPS_M, relapse_event) ~ KM_risk_group, data = hazard_ratio_data_AC)
hazard_ratios_cox_model_survival_AC <- coxph(Surv(OS_m, Event) ~ KM_risk_group, data = hazard_ratio_data_AC)
summary(hazard_ratios_cox_model_relapse_AC)
summary(hazard_ratios_cox_model_survival_AC)

hazard_ratio_data_SCC <- patients_with_malignancy_ranking_and_relapse %>% filter(Entity == 'SCC', UICC8_edition3=='I') 
hazard_ratios_cox_model_relapse_SCC <- coxph(Surv(censored_TIME_TO_RELAPS_M, relapse_event) ~ KM_risk_group, data = hazard_ratio_data_SCC)
hazard_ratios_cox_model_survival_SCC <- coxph(Surv(OS_m, Event) ~ KM_risk_group, data = hazard_ratio_data_SCC)
summary(hazard_ratios_cox_model_relapse_SCC)
summary(hazard_ratios_cox_model_survival_SCC)

#ggforest(hazard_ratios_cox_model, data = hazard_ratio_data, ref='r2')

patients_with_malignancy_ranking_and_relapse %>% group_by(Entity, KM_risk_group, UICC8_edition3) %>% 
  dplyr::summarize(N=n())


#hazard ratios multivariate across all stages
for_hazard_ratios_multivariate_data <- patients_with_malignancy_ranking_and_relapse %>% 
  dplyr::mutate(ki67 = as.numeric(ki67), Grade = as.numeric(Grade)) %>% 
  #filter(!is.na(ki67), !is.na(Grade)) %>% 
  filter(!is.na(ki67)) %>% 
  filter(Entity == 'AC') %>% 
  filter(UICC8_edition3 == 'I') %>% 
  mutate(numeric_adj_therapy = ifelse(Adj_therapy == 'yes', 1, 0))

dim(for_hazard_ratios_multivariate_data)
hazard_ratios_cox_model_multivariate <- coxph(Surv(censored_TIME_TO_RELAPS_M, relapse_event) ~ malignancy_rank + ki67, data = for_hazard_ratios_multivariate_data) #%>% filter(KM_risk_group != 'r1'))
summary(hazard_ratios_cox_model_multivariate)


rcorr(for_hazard_ratios_multivariate_data$ki67, for_hazard_ratios_multivariate_data$malignancy_rank)

patients_with_malignancy_ranking_and_relapse_all %>% 
  ungroup() %>%   
  group_by(Entity, UICC8_edition3) %>% 
  dplyr::summarize(c_index = get_c_index(malignancy_rank, censored_TIME_TO_RELAPS_M, relapse_event),
                  P = get_c_index_p(malignancy_rank, censored_TIME_TO_RELAPS_M, relapse_event),
                  N=n()) 

patients_with_malignancy_ranking_and_relapse %>% 
  ungroup() %>%   
  group_by(Entity) %>% 
  dplyr::summarize(c_index = get_c_index(malignancy_rank, censored_TIME_TO_RELAPS_M, relapse_event),
                  P = get_c_index_p(malignancy_rank, censored_TIME_TO_RELAPS_M, relapse_event),
                  N=n()) 

#cindex ki67
patients_with_malignancy_ranking_and_relapse %>% 
  mutate(ki67 = as.numeric(ki67)) %>%
  filter(!is.na(ki67)) %>% 
  ungroup() %>%   
  group_by(Entity, UICC8_edition3) %>% 
  dplyr::summarize(c_index = get_c_index(ki67, censored_TIME_TO_RELAPS_M, relapse_event),
                  P = get_c_index_p(ki67, censored_TIME_TO_RELAPS_M, relapse_event),
                  N=n(), 
                  n_events = sum(relapse_event)) 



##############
#of how malignant cells is each patient's tumor comprised
###############
malignant_ranking_factor <- malignant_ranking %>% arrange(malignancy_rank) %>% 
  mutate(malignancy_rank_ordered = factor(malignancy_rank, levels = (malignancy_rank))) %>% 
  mutate(ordered_cluster = factor(cluster, levels = rev(c(cluster))))

patients_with_malignancy_ranking <- patients_with_malignancy_ranking %>%
  distinct(Pseudo, .keep_all = TRUE)

patient_meta_data <- patient_meta_data %>%
  distinct(Pseudo, .keep_all = TRUE)

patients_with_UICC <- patient_meta_data %>% dplyr::select(Pseudo, UICC8_edition3) %>% 
  left_join(patients_with_malignancy_ranking %>% dplyr::select(Pseudo, malignancy_rank) %>% ungroup()) %>%
  #arrange(UICC8_edition3) %>% 
  arrange(malignancy_rank) %>% 
  mutate(ordered_pseudo = factor(Pseudo, levels = Pseudo)) %>% 
  dplyr::select(Pseudo, ordered_pseudo)


malignancy_composition <- patient_vector_long %>% 
  mutate(cluster = paste0('X', cluster)) %>%
  left_join(cell_types_fine_and_coarse) %>% 
  filter(Tier_2 == 'Tumor') %>% 
  group_by(Pseudo) %>%
  mutate(proportion_within_tumorcells = proportion/sum(proportion)) %>%
  left_join(patient_meta_data %>% dplyr::select(Pseudo, UICC8_edition3, OS_m, Event, Entity)) %>%
  left_join(malignant_ranking_factor) %>% 
  ungroup() %>%
  left_join(patients_with_UICC) %>%
  ungroup() %>% 
  mutate(scMP = malignancy_rank, UICC = as.factor(UICC8_edition3)) 

UICC_legend <- ggplot(malignancy_composition  %>% dplyr::select(ordered_pseudo, UICC) %>% unique(), aes(x = ordered_pseudo,y=1, fill = UICC)) +
  geom_bar(stat='identity', width=1) +
  theme_void() +
  theme(legend.text=element_text(size=15), 
  legend.title=element_text(size=15), 
  legend.position='bottom') +
  scale_fill_manual(values = c("#11B5AE", "#4046C9", "#F68512", "#DE3C82"))
 
 malignancy_rank_per_patient_plot <-plot_grid(ggplot(malignancy_composition, aes(x = ordered_pseudo, y  = proportion_within_tumorcells, group= ordered_cluster, fill = scMP)) +
    geom_bar(stat = 'identity') +
    scale_fill_viridis_c(option='D') +
      theme_classic() +
      ylab('Proportion of tumor cells') + 
      xlab(element_blank())+
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank(), 
            text=element_text(size=15)), 
  UICC_legend, ncol=1, rel_heights = c(10,1), align='v')

#pdf('./figures/malignancy_ranks_per_patient.pdf', width=12, height=7, onefile = F)
malignancy_rank_per_patient_plot
#dev.off()

#png('./figures/malignancy_ranks_per_patient.png', width=2000, height=1200, res=200)
malignancy_rank_per_patient_plot
#dev.off()

###Manuskript deskriptive scMP Ergebnisse
scMP_descriptive_statistics <- malignancy_composition %>% 
  group_by(Pseudo) %>% 
  dplyr::summarize(N=n(), iqr = IQR(malignancy_rank), median_malignancy_rank = median(malignancy_rank), min_malignancy_rank = min(malignancy_rank))

rcorr(scMP_descriptive_statistics$iqr, scMP_descriptive_statistics$median_malignancy_rank, type = 'spearman')

#descriptive by entity - PMP
PMP_descriptive_statistics <- patients_with_malignancy_ranking %>% 
  group_by(Entity, UICC8_edition3) %>% 
  dplyr::summarize(mean_PMP = mean(malignancy_rank), median_PMP = median(malignancy_rank), std_PMP = sd(malignancy_rank), IQR_PMP = IQR(malignancy_rank))

PMP_in_III_IV <- PMP_descriptive_statistics <- patients_with_malignancy_ranking %>% 
  filter(UICC8_edition3 %in% c('III', 'IV')) %>%
  group_by(Entity) %>% 
  dplyr::summarize(mean_PMP = mean(malignancy_rank), median_PMP = median(malignancy_rank), std_PMP = sd(malignancy_rank), IQR_PMP = IQR(malignancy_rank))

################
#Reviewer question: compare scMP with Stage
################

correlation_scMP_Stage <- malignancy_composition %>% 
  mutate(UICC8_numeric = as.numeric(factor(UICC8_edition3, levels = c('I', 'II', 'III', 'IV')))) %>% 
  group_by(Entity,UICC8_edition3,UICC8_numeric,  Pseudo) %>% 
  dplyr::summarize(weighted_scMP = sum(proportion_within_tumorcells * scMP)) %>% 
  group_by(Entity) %>% 
  dplyr::summarize(corr_v = rcorr(UICC8_numeric, weighted_scMP, type = 'spearman')$r[1,2])


Distribution_scMP_Stage <- cluster_result %>% 
  filter(Tier_2 == 'Tumor') %>% 
  mutate(cluster = paste0('X', cluster)) %>% 
  left_join(malignant_ranking) %>% 
  left_join(patient_meta_data %>% dplyr::select(Pseudo, UICC8_edition3, Entity)) %>% 
    mutate(UICC8_numeric = as.numeric(factor(UICC8_edition3, levels = c('I', 'II', 'III', 'IV'))))  %>% 
    mutate(Entity = ifelse(Entity == 'AC', 'LUAD', 'LUSC')) %>% 
    mutate(UICC = UICC8_edition3)


correlations_Distribution_scMP_Stage <- Distribution_scMP_Stage %>% 
    group_by(Entity) %>% 
  dplyr::summarize(corr_v = rcorr(UICC8_numeric, malignancy_rank, type = 'spearman')$r[1,2], 
  p_v = rcorr(UICC8_numeric, malignancy_rank, type = 'spearman')$P[1,2])

correlations_Distribution_scMP_Stage$p_v


#png('./figures/scMP_vs_Staging_Revision.png', width=2000, height=1500, res=200)
ggplot(Distribution_scMP_Stage, aes(x = malignancy_rank, fill = UICC, color = UICC)) +
  facet_wrap(~Entity, ncol=1) +
  geom_density(alpha = 0.4, position = "identity") +
  xlab('scMP') +
  theme_minimal() +
  theme(text=element_text(size=15))
#dev.off()

View(head(correlation_scMP_Stage))
##################
#external validation for malignancy rank
##################
#matched_cluster_Frame <- data.frame(id = names(matched_clusters$orig.ident), Pseudo = matched_clusters$orig.ident, cluster = matched_clusters$predicted.id)
validation_meta_data <- read.csv('../data/Met_Primary.csv')
validation_meta_data_any_metastasis <- read.csv('./use_data/Matched_LC_Data_no_name.csv') %>% 
  filter(metastasis_cat ==1) %>% 
  dplyr::select(Pseudo)

aligned_clusters <- read_parquet('../scGPT_embeddings/clusterings/external_data_clustered_aligned.parquet')
matched_cluster_Frame <- aligned_clusters


validation_meta_data_malignancy_rank <- read.csv('./use_data/Matched_LC_Data_no_name.csv') %>% 
  dplyr::select(Pseudo, Entity = Histology)

validation_data_malignancy_rank <- matched_cluster_Frame %>% 
  left_join(validation_meta_data_malignancy_rank) %>%
  #filter(Pseudo %in% validation_meta_data$Pseudo) %>%
  mutate(cluster = paste0('X', cluster)) %>%
  group_by(cluster, Pseudo, Entity) %>% 
  dplyr::summarize(countv = n()) %>%
  left_join(cell_types_fine_and_coarse) %>% 
  filter(Tier_2 == 'Tumor') %>% 
  #filter(Entity == 'AC') %>%
  left_join(malignant_ranking) %>% 
  group_by(Pseudo, Entity) %>% 
  mutate(proportion_within_tumorcells = countv / sum(countv)) %>%
  #filter(proportion_within_tumorcells>=0.005, countv>=0) %>%
  filter(proportion_within_tumorcells>=0.002, countv>=2) %>%
  mutate(kth_max = get_kth_max(malignancy_rank)) %>% 
  filter(malignancy_rank >= kth_max) %>%
  #filter(proportion_within_tumorcells>0.001) %>%
  dplyr::summarize(#average_malignancy = sum(malignancy_rank*proportion_within_tumorcells), 
                                            malignancy_rank = mean(malignancy_rank), 
                                            min_malignancy = quantile(malignancy_rank, 0.1, na.rm=T)) %>% 
  #dplyr::summarize(average_malignancy = mean(malignancy_rank), malignancy_rank = max(malignancy_rank)) %>% 
  mutate(is_brain_metastasis = ifelse(Pseudo %in% validation_meta_data$Pseudo,1,0)) %>% 
  mutate(develops_any_metastasis = ifelse(Pseudo %in% validation_meta_data_any_metastasis$Pseudo,1,0))   


#write.csv(validation_data_malignancy_rank, './use_data/external_data_average_malignancy.csv')
#write.csv(validation_data_malignancy_rank %>% mutate(PMP = malignancy_rank), './use_data/external_data_patient_metastatic_potential.csv')

#against survival

external_surviva_data <- read.csv('./use_data/Matched_LC_Data_no_name.csv') %>% 
  dplyr::select(Pseudo, Entity = Histology, UICC8, UICC8_alt, date_of_diagnosis, last_contact, status, dob) %>%
  mutate(Time = 1/365 * as.numeric(difftime(as.Date(last_contact, format="%d.%m.%Y"), as.Date(date_of_diagnosis, format= "%d.%m.%Y"), units='days'))) %>% 
  mutate(dob_date = as.Date(dob, format="%d.%m.%Y"), date_of_diagnosis_date = as.Date(date_of_diagnosis, format= "%d.%m.%Y")) %>%
  mutate(age = 1/365 * as.numeric(difftime(date_of_diagnosis_date, dob_date, units='days'))) %>%
  left_join(validation_data_malignancy_rank) %>% 
  filter(!is.na(last_contact), last_contact != "") %>% 
  #filter(Entity == 'AC') %>%
  mutate(UICC8 = case_when(
    UICC8 %in% c('IA1', 'IA2', 'IA3', 'IB') ~ 'I',
    UICC8 %in% c('IIA', 'IIB') ~ 'II',
    UICC8 %in% c('IIIA', 'IIIB') ~ 'III',
    UICC8 %in% c('IVA', 'IVB') ~ 'IV',
    T ~ UICC8
  )) %>% 
  group_by(UICC8, Entity) %>%
  filter(!is.na(malignancy_rank)) %>% 
  filter(!is.na(Time)) %>% 
  filter(!(Pseudo %in% validation_meta_data$Met)) %>% 
  dplyr::mutate(KM_risk_group =  case_when(#malignancy_rank>=99.1 ~'r4', 
                                      malignancy_rank<110 & malignancy_rank>=94.6~'r2', 
                                      #malignancy_rank<92.8 & malignancy_rank>=84.5~'r2', 
                                      malignancy_rank<94.6 ~'r1')) 

                        

    #mutate(KM_risk_group = ntile(malignancy_rank,4)) 

gradient_colors <- scales::gradient_n_pal(c("#049404",'#ffc400', "red"))(seq(0, 1, length.out = 3))

external_surviva_data$Time_m <- external_surviva_data$Time

external_surviva_data$Time <- external_surviva_data$Time_m*12

# Create truncated survival variables: any OS_m > 60 becomes 60, and the event is set to 0 (censored)
external_surviva_data <- external_surviva_data %>%
  mutate(Time = pmin(Time, 60),
         status = if_else(Time > 60, 0, status))

external_surviva_data = external_surviva_data %>% filter(UICC8 == 'I')
external_surviva_data = external_surviva_data %>% filter(Entity == 'AC')


#gradient_colors <- scales::gradient_n_pal(c("#13B5AE",'#3F46C9', "#DE3C82"))(seq(0, 1, length.out = nquantiles))

#pdf('./figures/Survival_external_AC.pdf', width=5, height=6,onefile = F)
plotfile <- ggsurvplot(
  surv_fit_malignancy_ranking <- survfit(Surv(Time, status) ~ KM_risk_group, data = external_surviva_data),
  data = external_surviva_data,
  #facet.by = c('UICC8_edition3'),
  #facet.by = c('Entity'),
  pval = T,           # Add p-value
  pval.method = TRUE,    # Show test method in the p-value plot
  pval.size = 4,         # Set the p-value font size
  conf.int = F,       # Add confidence intervals
  censor = TRUE,         # Show censored data points
  legend.title = "Risk quantile",
  title = 'Survival External Dataset',
  xlab = "Time",
  ylab = "Survival Probability",
  break.time.by = 12,
  risk.table = T,
  palette = c("#13B5AE",'#3F46C9', "#F68512"),
  linewidth=0.1,
  #xlim = c(0,60),
  ggtheme=theme_minimal(),
  text=element_text(size=20))

plotfile$table <- plotfile$table + 
  scale_y_discrete(labels = c("RG3", "RG2", "RG1")) # Customize y-axis labels

plotfile
#dev.off()




#png('./figures/validation_malignancy_rank_survival.png', width=2000, height=1500, res=200)
ggsurvplot_facet(
  surv_fit_external <- survfit(Surv(Time, status) ~ KM_risk_group, data = external_surviva_data),
  data = external_surviva_data %>% filter(UICC8=='I'),
  facet.by = c('UICC8', 'Entity'),
  #facet.by = 1,
  #facet.by = c('Entity'),
  pval = F,           # Add p-value
  pval.method = TRUE,    # Show test method in the p-value plot
  pval.size = 4,         # Set the p-value font size
  conf.int = F,       # Add confidence intervals
  censor = TRUE,         # Show censored data points
  legend.title = "Risk quantile",
  title = 'Overall survival',
  xlab = "Time",
  ylab = "Survival Probability",
  palette = gradient_colors,
  linewidth=0.1,
  #xlim = c(0,2000),
  xlim = c(0,10),
  ggtheme=theme_minimal()
)
#dev.off()

external_surviva_data %>% 
  group_by(UICC8)

external_surviva_data %>% 
  group_by(Entity, UICC8) %>% 
  dplyr::summarize(c_index = get_c_index(malignancy_rank, Time, status),
                  P = get_c_index_p(malignancy_rank, Time, status), 
                  N=n()) 

external_surviva_data %>% 
  group_by(Entity) %>% 
  dplyr::summarize(c_index = get_c_index(malignancy_rank, Time, status),
                  P = get_c_index_p(malignancy_rank, Time, status), 
                  N=n()) 


external_cox_model <- coxph(Surv(Time, status) ~ UICC8  +  malignancy_rank, data = external_surviva_data)

# View the summary of the model
summary(external_cox_model)

########################
#validation relapse
#########################
external_surviva_data_for_relapse <- read.csv('./use_data/Matched_LC_Data_no_name.csv') %>% 
  dplyr::select(Pseudo, Entity = Histology, UICC8, UICC8_alt, date_of_diagnosis, last_contact, status, dob, recurrence_cat, recurrrence_date1) %>%
  mutate(Time2relapse = 1/365 * as.numeric(difftime(as.Date(recurrrence_date1, format="%d.%m.%y"), as.Date(date_of_diagnosis, format= "%d.%m.%Y"), units='days'))) %>% 
  mutate(dob_date = as.Date(dob, format="%d.%m.%Y"), date_of_diagnosis_date = as.Date(date_of_diagnosis, format= "%d.%m.%Y")) %>%
  mutate(age = 1/365 * as.numeric(difftime(date_of_diagnosis_date, dob_date, units='days'))) %>%
  mutate(last_contact_date = as.Date(last_contact, format= "%d.%m.%Y")) %>%
  mutate(TIME2LASTCONTACT = as.numeric(interval(date_of_diagnosis_date, last_contact_date) / years(1))) %>%
  mutate(censored_time_to_relapse = ifelse(recurrence_cat == 1, Time2relapse, TIME2LASTCONTACT)) %>%
  left_join(validation_data_malignancy_rank) %>%
  filter(!is.na(last_contact), last_contact != "") %>% 
  #filter(Entity == 'AC') %>%
  mutate(UICC8 = case_when(
    UICC8 %in% c('IA1', 'IA2', 'IA3', 'IB') ~ 'I',
    UICC8 %in% c('IIA', 'IIB') ~ 'II',
    UICC8 %in% c('IIIA', 'IIIB') ~ 'III',
    UICC8 %in% c('IVA', 'IVB') ~ 'IV',
    T ~ UICC8
  )) %>% 
  group_by(UICC8, Entity) %>%
  #mutate(KM_risk_group = ntile(malignancy_rank,4)) %>% 
  filter(!is.na(malignancy_rank)) %>%  
  filter(!is.na(recurrence_cat)) %>% 
  filter(!(Pseudo %in% validation_meta_data$Met)) %>% 
  mutate(has_relapse_within_2_years = (recurrence_cat==1 & Time2relapse<=2)) %>% 
  dplyr::mutate(KM_risk_group =  case_when(#malignancy_rank>=99.1 ~'r4', 
                                      malignancy_rank>=94.36~'r2',    #92.8, 84.5
                                      #malignancy_rank<92.5 & malignancy_rank>=84.8~'r2', 
                                      malignancy_rank<94.36 ~'r1'))



ggplot(external_surviva_data_for_relapse, aes(x = has_relapse_within_2_years, y = malignancy_rank)) +
  geom_violin() +
    #geom_point() +
  geom_jitter() +
  facet_grid(Entity~UICC8)






external_surviva_data_for_relapse$Time_m <- external_surviva_data_for_relapse$censored_time_to_relapse

external_surviva_data_for_relapse$Time <- external_surviva_data_for_relapse$Time_m*12

# Create truncated survival variables: any OS_m > 60 becomes 60, and the event is set to 0 (censored)
external_surviva_data_for_relapse <- external_surviva_data_for_relapse %>%
  mutate(Time = pmin(Time, 60),
         recurrence_cat = if_else(Time > 60, 0, recurrence_cat))

external_surviva_data_for_relapse = external_surviva_data_for_relapse %>% filter(UICC8 == 'I')
external_surviva_data_for_relapse = external_surviva_data_for_relapse %>% filter(Entity == 'AC')


#gradient_colors <- scales::gradient_n_pal(c("#13B5AE",'#3F46C9', "#DE3C82"))(seq(0, 1, length.out = nquantiles))

#pdf('./figures/Relapse_external_AC.pdf', width=5, height=6,onefile = F)
plotfile <- ggsurvplot(
  surv_fit_malignancy_ranking_external <- survfit(Surv(Time, recurrence_cat) ~ KM_risk_group, data = external_surviva_data_for_relapse),
  data = external_surviva_data_for_relapse,
  #facet.by = c('UICC8'),
  #facet.by = c('Entity'),
  pval = T,           # Add p-value
  pval.method = TRUE,    # Show test method in the p-value plot
  pval.size = 4,         # Set the p-value font size
  conf.int = F,       # Add confidence intervals
  censor = TRUE,         # Show censored data points
  legend.title = "Risk quantile",
  title = 'Relapse External Dataset',
  xlab = "Time",
  ylab = "Survival Probability",
  break.time.by = 12,
  risk.table = T,
  palette = c("#13B5AE",'#3F46C9', "#F68512"),
  linewidth=0.1,
  #xlim = c(0,60),
  ggtheme=theme_minimal(),
  text=element_text(size=20))

plotfile$table <- plotfile$table + 
  scale_y_discrete(labels = c("RG3", "RG2", "RG1")) # Customize y-axis labels

plotfile
#dev.off()

summary_5yr_external <- summary(surv_fit_malignancy_ranking_external, times = 60)
summary_5yr_external



#validation KM curve
#png('./figures/validation_malignancy_rank_relapse.png', width=2000, height=1500, res=200)
ggsurvplot_facet(
  surv_fit_malignancy_ranking <- survfit(Surv(censored_time_to_relapse, recurrence_cat) ~ KM_risk_group, data = external_surviva_data_for_relapse),
  data = external_surviva_data_for_relapse,# %>% filter(UICC8=='I'),
  #facet.by = c('UICC8_edition3'),
  facet.by = c('UICC8','Entity'),
  pval = F,           # Add p-value
  pval.method = TRUE,    # Show test method in the p-value plot
  pval.size = 4,         # Set the p-value font size
  conf.int = F,       # Add confidence intervals
  censor = TRUE,         # Show censored data points
  legend.title = "Risk quantile",
  title = 'Time to relapse',
  xlab = "Time",
  ylab = "Probability",
  palette = gradient_colors_relapse,
  linewidth=0.1,
  xlim = c(0,5),
  ncol=2,
  ggtheme=theme_minimal(),
  text=element_text(size=20))
#dev.off()


external_surviva_data_for_relapse %>% 
  ungroup() %>%   
  group_by(Entity, UICC8) %>% 
  dplyr::summarize(c_index = get_c_index(malignancy_rank, censored_time_to_relapse, recurrence_cat),
                  P = (get_c_index_p(malignancy_rank, censored_time_to_relapse, recurrence_cat)),
                  N=n())  #%>% .$P


#########hazard ratios external


hazard_ratio_data_AC_external_relapse <- external_surviva_data_for_relapse %>% filter(Entity == 'AC', UICC8=='I') 
hazard_ratios_cox_model_relapse_AC_external <- coxph(Surv(censored_time_to_relapse, recurrence_cat) ~ KM_risk_group, data = hazard_ratio_data_AC_external_relapse)
hazard_ratio_data_AC_external_os <- external_surviva_data %>% filter(Entity == 'AC', UICC8=='I') 
hazard_ratios_cox_model_os_AC_external <- coxph(Surv(Time, status) ~ KM_risk_group, data = hazard_ratio_data_AC_external_os)

hazard_ratios_cox_model_survival_AC <- coxph(Surv(OS_m, Event) ~ KM_risk_group, data = hazard_ratio_data_AC)
summary(hazard_ratios_cox_model_relapse_AC_external)
summary(hazard_ratios_cox_model_os_AC_external)



###################################
#compare primary tumors and brain metastases, from NatMed paper
###################################
library(ggridges)
scMP_scores_PT_BM <- read.csv('/mnt/ssd/shared/LungCAIRE_SingleCell/plots_statistics/use_data/ingest_genebasis_external_clusters_PTvsBM_with_scMP.csv') %>% 
  mutate(group = ifelse(group == 'PT','Primary Tumor', 'Brain Metastasis'))

average <- scMP_scores_PT_BM %>% group_by(group) %>% 
  dplyr::summarize(mean_scMP = mean(malignancy_rank), 
                  median_scMP = median(malignancy_rank), 
                  std_scMP = sd(malignancy_rank), 
                  IQR_scMP = IQR(malignancy_rank))
                  
average_per_patient <- scMP_scores_PT_BM %>% group_by(Pseudo, group) %>% 
  dplyr::summarize(mean_scMP = mean(malignancy_rank), 
                  median_scMP = median(malignancy_rank), 
                  std_scMP = sd(malignancy_rank), 
                  IQR_scMP = IQR(malignancy_rank), 
                  N=n())

  



pdf('./figures/External_validation_PT_BM.pdf', width=12, height=10,onefile = F)
ggplot(scMP_scores_PT_BM, aes(x = Pseudo, y = malignancy_rank, fill = group)) + 
  geom_boxplot(outlier.size=0.2, alpha=0.5) + 
  geom_hline(data = average, aes(yintercept = mean_scMP, color = group), linetype='dashed',  linewidth=3) +
  scale_fill_manual(values = c( "#F68512", '#3F46C9')) +
  scale_color_manual(values = c( "#F68512", '#3F46C9')) +
  #geom_violin() +
    theme_minimal()+ 
    ylab('scMP') +
    xlab('') + 
  theme(axis.text.x = element_text(angle=90), 
    axis.text = element_text(size=15),
    legend.text = element_text(size=15),
    axis.title = element_text(size=17),
    legend.title = element_blank())
dev.off()

png('./figures/External_validation_PT_BM.png', width=2000, height=1600, res=200)
ggplot(scMP_scores_PT_BM, aes(x = Pseudo, y = malignancy_rank, fill = group)) + 
  geom_boxplot(outlier.size=0.2, alpha=0.5) + 
  geom_hline(data = average, aes(yintercept = mean_scMP, color = group), linetype='dashed',  linewidth=3) +
  scale_fill_manual(values = c( "#F68512", '#3F46C9')) +
  scale_color_manual(values = c( "#F68512", '#3F46C9')) +
  #geom_violin() +
    theme_minimal()+ 
    ylab('scMP') +
    xlab('') + 
  theme(axis.text.x = element_text(angle=90), 
    axis.text = element_text(size=15),
    legend.text = element_text(size=15),
    axis.title = element_text(size=17),
    legend.title = element_blank())
dev.off()

ggplot(average_per_patient, aes(x = group, y = mean_scMP, fill = group)) + 
  geom_boxplot() + 
  geom_point(aes(x = group, y = mean_scMP), alpha = 0.5) + 
  #geom_hline(data = average, aes(x = group, yintercept = mean_scMP, color = group), linetype='dashed',  linewidth=3) +
  #geom_violin() +
  theme_minimal()


wilcox.test(malignancy_rank ~ group, data = scMP_scores_PT_BM)
average



scMP_scores_PT_BM %>% group_by(group) %>% 
  dplyr::select(Pseudo, group) %>%
  unique() %>% 
  group_by(group) %>% 
  dplyr::summarize(N=n())


write.csv(scMP_scores_PT_BM %>% dplyr::select(Pseudo, malignancy_rank, group), './figures/figure_data/F3G_PTBM.csv')




#####calculate PMP for external data

get_3rd_max <- function(x){
  if (length(x)==0) {return(0)}
  x = sort(x, decreasing=T)
  k <- 3
  res <- ifelse(length(x)>=k, x[k], min(x))
}



patients_with_malignancy_ranking_external <- scMP_scores_PT_BM %>% 
  mutate(cluster = paste0('X', cluster)) %>%
  left_join(cell_types_fine_and_coarse) %>% 
  filter(Tier_2 == 'Tumor') %>% 
  group_by(Pseudo, group, cluster, malignancy_rank) %>%
  dplyr::summarize(sumv = n(), N=n()) %>% 
  group_by(Pseudo, group, malignancy_rank) %>% 
  mutate(proportion_within_tumorcells = sumv/sum(sumv)) %>%
  #inner_join(malignant_ranking) %>%
  #filter(proportion_within_tumorcells>=0.005) %>% 
  filter(proportion_within_tumorcells>0.002, N>=2) %>% #0.001
  #filter(proportion_within_tumorcells>0.002) %>% #0.001
  group_by(Pseudo) %>% 
  mutate(kth_max = get_3rd_max(malignancy_rank)) %>%
  ungroup() %>%
  filter(malignancy_rank >= kth_max) %>%
  dplyr::group_by(Pseudo) %>%   #calculate malignancy rank for AC and SCC together
  dplyr::summarize(malignancy_rank = mean(malignancy_rank)) %>% #weigh malignancy rank by proportion?
  ungroup() %>% 
  mutate(dataset = 'External_PT_BM') %>% 
  filter(startsWith(Pseudo, 'PT')) %>% 
  rbind(patients_with_malignancy_ranking %>% filter(Entity == 'AC', UICC8_edition3 == 'IV') %>% dplyr::select(Pseudo, malignancy_rank) %>% mutate(dataset = 'Main_Cohort')) %>% 
  ungroup() %>%
  mutate(rankl = rank(malignancy_rank)) %>% 
  mutate(rankl = rankl/max(rankl))




ggplot(patients_with_malignancy_ranking_external %>% filter(dataset == 'External_PT_BM'), aes(y = Pseudo, x = malignancy_rank, color = dataset)) +
  #geom_boxplot() +
  geom_jitter() +
  theme_minimal() +
  theme(legend.position = 'none') +
  ylab('scMP') +
  xlab('Dataset') +
  theme(text=element_text(size=15))
