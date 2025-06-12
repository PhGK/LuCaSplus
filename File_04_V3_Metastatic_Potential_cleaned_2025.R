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
  dplyr::mutate(malignancy_rank = rank(r)) %>% 
  ungroup() %>%
  mutate(malignancy_rank = malignancy_rank / max(malignancy_rank) * 100 ) %>%
  mutate(cluster = paste0('X', cluster)) 

write.csv(malignant_ranking, './use_data/latent_metastatic_cells_candidates.csv')
  


get_kth_max <- function(x){
  if (length(x)==0) {return(0)}
  x = sort(x, decreasing=T)
  k <- 3
  res <- ifelse(length(x)>=k, x[k], min(x))
}

patients_with_malignancy_ranking <- patient_vector_long %>% 
  mutate(cluster = paste0('X', cluster)) %>%
  left_join(cell_types_fine_and_coarse) %>% 
  filter(Tier_2 == 'Tumor') %>% 
  group_by(Pseudo) %>%
  mutate(proportion_within_tumorcells = proportion/sum(proportion)) %>%
  left_join(patient_meta_data %>% dplyr::select(Pseudo, UICC8_edition3, OS_m, Event, Entity)) %>%
  left_join(malignant_ranking) %>%
  filter(proportion_within_tumorcells>0.001) %>% 
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
                                      malignancy_rank<210 & malignancy_rank>=92.5~'r3', 
                                      malignancy_rank<92.5 & malignancy_rank>=84.8~'r2', 
                                      malignancy_rank<84.8 ~'r1'), 
                                      case_when(#malignancy_rank>=99.1 ~'r4', 
                                      malignancy_rank<210 & malignancy_rank>=92.5~'r3', 
                                      malignancy_rank<92.5 & malignancy_rank>=84.8~'r2', 
                                      malignancy_rank<84.8 ~'r1'))) 


patients_with_malignancy_ranking %>% 
  group_by(Entity, KM_risk_group) %>% 
  dplyr::summarize(N= n())

colnames(patient_meta_data)  

write.csv(patients_with_malignancy_ranking, './use_data/patient_metastatic_potential.csv')

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

png('./figures/Survival_internal.png', width=2000, height=1500, res=250)
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

View(patients_with_malignancy_ranking_and_relapse)

write.csv(patients_with_malignancy_ranking_and_relapse %>% dplyr::select(Pseudo, censored_TIME_TO_RELAPS_M, relapse_event), './use_data/time2relapse.csv')

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
  filter(malignancy_rank>cutoff_for_no_relapse_within_5_years) %>%
  dplyr::summarize(middle_value = (100+min(malignancy_rank))/2)

cutoff_for_middle_groups

cutoff_for_no_relapse_within_5_years_SCC <- patients_with_malignancy_ranking_and_relapse %>% 
  filter(relapse_event==1, UICC8_edition3 %in% c('I')) %>%
  filter(Entity == 'SCC') %>%
  .$malignancy_rank %>% min()
cutoff_for_no_relapse_within_5_years

cutoff_for_middle_groups_SCC <- patients_with_malignancy_ranking_and_relapse %>%   
  filter(Entity == 'SCC') %>%
  filter(malignancy_rank>cutoff_for_no_relapse_within_5_years_SCC) %>%
  dplyr::summarize(middle_value = (100+min(malignancy_rank))/2)
cutoff_for_middle_groups_SCC


#now do KM to relapse
a <- patients_with_malignancy_ranking_and_relapse %>% 
  dplyr::select(Pseudo, OS_m, Event, censored_TIME_TO_RELAPS_M, relapse_event) %>% 
  mutate(OS_m >= censored_TIME_TO_RELAPS_M)

patients_with_malignancy_ranking_and_relapse$KM_risk_group






# Create truncated survival variables: any OS_m > 60 becomes 60, and the event is set to 0 (censored)
malignancy_ranking_for_KM <- malignancy_ranking_for_KM %>%
  mutate(OS_m_trunc = pmin(OS_m, 60),
         Event_trunc = if_else(OS_m > 60, 0, Event))

malignancy_ranking_for_KM_2 = malignancy_ranking_for_KM %>% filter(UICC8_edition3 == 'I')
malignancy_ranking_for_KM_2 = malignancy_ranking_for_KM_2 %>% filter(Entity == 'AC')


#gradient_colors <- scales::gradient_n_pal(c("#13B5AE",'#3F46C9', "#DE3C82"))(seq(0, 1, length.out = nquantiles))

pdf('./figures/Survival_internal_AC.pdf', width=5, height=6,onefile = F)
plotfile <- ggsurvplot(
  surv_fit_malignancy_ranking <- survfit(Surv(OS_m_trunc, Event) ~ KM_risk_group, data = malignancy_ranking_for_KM_2),
  data = malignancy_ranking_for_KM_2,
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





# Create truncated survival variables: any OS_m > 60 becomes 60, and the event is set to 0 (censored)
patients_with_malignancy_ranking_and_relapse <- patients_with_malignancy_ranking_and_relapse %>%
  mutate(censored_TIME_TO_RELAPS_M = pmin(censored_TIME_TO_RELAPS_M, 60),
         relapse_event = if_else(censored_TIME_TO_RELAPS_M > 60, 0, relapse_event))

patients_with_malignancy_ranking_and_relapse = patients_with_malignancy_ranking_and_relapse %>% filter(UICC8_edition3 == 'I')
patients_with_malignancy_ranking_and_relapse = patients_with_malignancy_ranking_and_relapse %>% filter(Entity == 'AC')


#gradient_colors <- scales::gradient_n_pal(c("#13B5AE",'#3F46C9', "#DE3C82"))(seq(0, 1, length.out = nquantiles))

pdf('./figures/Relapse_internal_AC.pdf', width=5, height=6,onefile = F)
plotfile <- ggsurvplot(
  surv_fit_malignancy_ranking <- survfit(Surv(censored_TIME_TO_RELAPS_M, relapse_event) ~ KM_risk_group, data = patients_with_malignancy_ranking_and_relapse),
  data = patients_with_malignancy_ranking_and_relapse %>% filter(UICC8_edition3 %in% c('I')),
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



gradient_colors_relapse <- (scales::gradient_n_pal(c("#049404",'#ffc400', "red"))(seq(0, 1, length.out =patients_with_malignancy_ranking_and_relapse$KM_risk_group %>% unique() %>% length() )))

png('./figures/Relapse_internal.png', width=2000, height=1500, res=250)
ggsurvplot_facet(
  surv_fit_malignancy_ranking <- survfit(Surv(censored_TIME_TO_RELAPS_M, relapse_event) ~ KM_risk_group, data = patients_with_malignancy_ranking_and_relapse),
  data = patients_with_malignancy_ranking_and_relapse %>% filter(UICC8_edition3 %in% c('I')),
  #facet.by = c('UICC8_edition3'),
  facet.by = c('Entity'),
  #facet.by = c( 'UICC8_edition3', 'Entity'),
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

#hazard ratios across all stages
hazard_ratio_data <- patients_with_malignancy_ranking_and_relapse %>% filter(Entity == 'AC') %>% filter(KM_risk_group != 'r2')
hazard_ratios_cox_model <- coxph(Surv(censored_TIME_TO_RELAPS_M, relapse_event) ~ KM_risk_group, data = hazard_ratio_data)
summary(hazard_ratios_cox_model)

#ggforest(hazard_ratios_cox_model, data = hazard_ratio_data, ref='r2')

patients_with_malignancy_ranking_and_relapse %>% group_by(Entity, KM_risk_group, UICC8_edition3) %>% 
  dplyr::summarize(N=n())


#hazard ratios multivariate across all stages
for_hazard_ratios_multivariate_data <- patients_with_malignancy_ranking_and_relapse %>% 
  dplyr::mutate(ki67 = as.numeric(ki67), Grade = as.numeric(Grade)) %>% 
  #filter(!is.na(ki67), !is.na(Grade)) %>% 
  filter(!is.na(ki67)) %>% 
  filter(Entity == 'AC') %>% 
  filter(UICC8_edition3 == 'I')

dim(for_hazard_ratios_multivariate_data)
hazard_ratios_cox_model_multivariate <- coxph(Surv(censored_TIME_TO_RELAPS_M, relapse_event) ~ malignancy_rank + ki67, data = for_hazard_ratios_multivariate_data %>% filter(KM_risk_group != 'r1'))
summary(hazard_ratios_cox_model_multivariate)

rcorr(for_hazard_ratios_multivariate_data$ki67, for_hazard_ratios_multivariate_data$malignancy_rank)

patients_with_malignancy_ranking_and_relapse %>% 
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
  ungroup() 

UICC_legend <- ggplot(malignancy_composition  %>% dplyr::select(ordered_pseudo, UICC8_edition3) %>% unique(), aes(x = ordered_pseudo,y=1, fill = as.factor(UICC8_edition3))) +
  geom_bar(stat='identity', width=1) +
  theme_void() +
  scale_fill_manual(values = c("#11B5AE", "#4046C9", "#F68512", "#DE3C82"))

pdf('./figures/malignancy_ranks_per_patient.pdf', width=12, height=4, onefile = F)
  plot_grid(ggplot(malignancy_composition, aes(x = ordered_pseudo, y  = proportion_within_tumorcells, group= ordered_cluster, fill = (malignancy_rank))) +
    geom_bar(stat = 'identity') +
    scale_fill_viridis_c(option='D') +
      theme_classic() +
      xlab(element_blank())+
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank()), 
  UICC_legend, ncol=1, rel_heights = c(10,0.5), align='v')
dev.off()


###Manuskript deskriptive scMP Ergebnisse
scMP_descriptive_statistics <- malignancy_composition %>% 
  group_by(Pseudo) %>% 
  dplyr::summarize(N=n(), iqr = IQR(malignancy_rank), median_malignancy_rank = median(malignancy_rank), min_malignancy_rank = min(malignancy_rank))

rcorr(scMP_descriptive_statistics$iqr, scMP_descriptive_statistics$median_malignancy_rank, type = 'spearman')



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
  filter(proportion_within_tumorcells>0.001, countv>=0) %>%
  mutate(kth_max = get_kth_max(malignancy_rank)) %>% 
  filter(malignancy_rank >= kth_max) %>%
  #filter(proportion_within_tumorcells>0.001) %>%
  dplyr::summarize(average_malignancy = sum(malignancy_rank*proportion_within_tumorcells), 
                                            malignancy_rank = mean(malignancy_rank), 
                                            min_malignancy = quantile(malignancy_rank, 0.1, na.rm=T)) %>% 
  #dplyr::summarize(average_malignancy = mean(malignancy_rank), malignancy_rank = max(malignancy_rank)) %>% 
  mutate(is_brain_metastasis = ifelse(Pseudo %in% validation_meta_data$Pseudo,1,0)) %>% 
  mutate(develops_any_metastasis = ifelse(Pseudo %in% validation_meta_data_any_metastasis$Pseudo,1,0))   


write.csv(validation_data_malignancy_rank, './use_data/external_data_average_malignancy.csv')

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
                                      malignancy_rank<110 & malignancy_rank>=92.8~'r3', 
                                      malignancy_rank<92.8 & malignancy_rank>=84.5~'r2', 
                                      malignancy_rank<84.55 ~'r1')) 

                        

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

pdf('./figures/Survival_external_AC.pdf', width=5, height=6,onefile = F)
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
dev.off()




png('./figures/validation_malignancy_rank_survival.png', width=2000, height=1500, res=200)
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
dev.off()

external_surviva_data %>% 
  group_by(UICC)

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
                                      malignancy_rank>=92.5~'r3',    #92.8, 84.5
                                      malignancy_rank<92.5 & malignancy_rank>=84.8~'r2', 
                                      malignancy_rank<84.8 ~'r1'))



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

pdf('./figures/Relapse_external_AC.pdf', width=5, height=6,onefile = F)
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
dev.off()

summary_5yr_external <- summary(surv_fit_malignancy_ranking_external, times = 60)
summary_5yr_external



#validation KM curve
png('./figures/validation_malignancy_rank_relapse.png', width=2000, height=1500, res=200)
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
dev.off()


external_surviva_data_for_relapse %>% 
  ungroup() %>%   
  group_by(Entity, UICC8) %>% 
  dplyr::summarize(c_index = get_c_index(malignancy_rank, censored_time_to_relapse, recurrence_cat),
                  P = (get_c_index_p(malignancy_rank, censored_time_to_relapse, recurrence_cat)),
                  N=n())  #%>% .$P



##########################
#validation of malignancy rank between primary tumor and metastasis
##########################
aligned_clusters_metprim <- read_parquet('../scGPT_embeddings/clusterings/external_data_clustered_aligned_MetPrim.parquet')

validation_meta_data_metastasis <- read.csv('./use_data/Matched_LC_Data_no_name.csv') %>% 
  dplyr::select(dob,Pseudo)

aligned_clusters_metprim$Pseudo %>% unique()

validation_data_metastasis <- aligned_clusters_metprim %>% 
  left_join(validation_meta_data_malignancy_rank) %>%
  mutate(cluster = paste0('X', cluster)) %>%
  group_by(cluster, Pseudo, Entity) %>% 
  dplyr::summarize(countv = n()) %>%
  left_join(cell_types_fine_and_coarse) %>% 
  filter(Tier_2 == 'Tumor') %>% 
  left_join(malignant_ranking) %>% 
  group_by(Pseudo, Entity) %>% 
  mutate(proportion_within_tumorcells = countv / sum(countv)) %>%
  filter(proportion_within_tumorcells>0.001) %>%
  mutate(is_brain_metastasis = ifelse(Pseudo %in% validation_meta_data$Met,1,0)) %>% 
  mutate(is_primary = ifelse(Pseudo %in% validation_meta_data$Primary,1,0)) %>% 
  #mutate(is_primary = Pseudo %in% c('LCVAL094', 'LCVAL095', 'LCVAL096')) %>% 
  filter(is_brain_metastasis | is_primary) %>% 
  mutate(location = ifelse(is_primary, 'primary', 'metastasis')) %>%
  left_join(validation_meta_data_metastasis) %>% 
  group_by(dob) %>% 
  mutate(has_match = length(unique(Pseudo))==2) %>% 
  filter(has_match)

  validation_data_metastasis_wide <- validation_data_metastasis %>% 
    ungroup() %>%
    dplyr::select(cluster, dob, malignancy_rank, r, proportion_within_tumorcells, location) %>% 
    pivot_wider(names_from = location, values_from = proportion_within_tumorcells) %>% 
    mutate(metastasis = ifelse(is.na(metastasis), 0, metastasis), 
          primary = ifelse(is.na(primary), 0, primary)) 

validation_data_metastasis_for_plot <- validation_data_metastasis_wide %>% 
  pivot_longer(c(primary, metastasis), values_to = 'proportion_within_tumorcells', names_to = 'location') %>% 
  ungroup() %>%
  mutate(malignancy_quantile = ntile(malignancy_rank,10)) %>% 
  group_by(cluster) %>% 
  mutate(primary_is_zero = mean(proportion_within_tumorcells[location == 'primary']) == 0) %>% 
  left_join(validation_data_metastasis %>% filter(location == 'primary') %>% dplyr::select(Primary_ID = Pseudo, dob) %>% unique())

validation_data_metastasis_for_plot_averages <- validation_data_metastasis_for_plot %>% 
  group_by(Primary_ID, location) %>% 
  dplyr::summarize(malignancy_rank = sum(proportion_within_tumorcells* malignancy_rank))

validation_data_metastasis_for_plot_quantiles <- validation_data_metastasis_for_plot %>% 
  group_by(malignancy_quantile, dob, location) %>% 
  dplyr::summarize(proportion_within_tumorcells = sum(proportion_within_tumorcells))

write.csv(validation_data_metastasis_for_plot, './use_data/comparison_primary_vs_metastasis_malignancy_rank.csv')

cor(validation_data_metastasis_wide$malignancy_rank, validation_data_metastasis_wide$metastasis, method = 'spearman')
cor(validation_data_metastasis_wide$malignancy_rank, validation_data_metastasis_wide$primary, method = 'spearman')

ggplot(validation_data_metastasis_for_plot, aes(x = malignancy_rank, y = proportion_within_tumorcells, color = location)) +
  #geom_point() +
  geom_smooth(method = 'lm') + 
  facet_wrap(~dob, ncol=1) 

ggplot(validation_data_metastasis_for_plot, aes(x = malignancy_rank, y = proportion_within_tumorcells, color = location)) +
  geom_point() +
  facet_wrap(~dob, ncol=1) 

ggplot(validation_data_metastasis_for_plot, aes(x = malignancy_rank, y = proportion_within_tumorcells, color = location)) +
  geom_point() +
  facet_grid(primary_is_zero~dob) 
 


primary_vs_met_plot <- ggplot() +
  #geom_point(data = validation_data_metastasis_for_plot %>% filter(location == 'metastasis'), aes(x = (malignancy_rank), y = proportion_within_tumorcells), color = 'darkgreen', size=2, alpha=1.0) +
  #geom_line(data = validation_data_metastasis_for_plot %>% filter(location == 'metastasis'), aes(x = (malignancy_rank), y = proportion_within_tumorcells), color = 'darkgreen', alpha=1.0) +
  geom_bar(data = validation_data_metastasis_for_plot %>% filter(location == 'primary'), aes(x = (malignancy_rank), y = proportion_within_tumorcells), stat = 'identity', fill ='#2525bf' , alpha=0.5) +
  geom_bar(data = validation_data_metastasis_for_plot %>% filter(location == 'metastasis'), aes(x = (malignancy_rank), y = proportion_within_tumorcells), stat = 'identity', fill = '#d22020', alpha=0.5) +
  geom_vline(data = validation_data_metastasis_for_plot_averages, aes(xintercept = malignancy_rank, color = location), linewidth=1, linetype=2) +
  facet_wrap(~Primary_ID, ncol=1, scales='free') +
  scale_color_manual(values = c('#d22020', '#2525bf')) +
  theme_minimal() +
  theme(text=element_text(size=20), 
        legend.title = element_blank(), 
        legend.position = 'bottom') +
  xlab('Metastatic potential') +
  ylab('Proportion')

primary_vs_met_plot

gradient_colors
##look only at differences between metastasis and primary

validation_data_metastasis_differences <- validation_data_metastasis_for_plot %>% 
  dplyr::select(-malignancy_quantile) %>%
  pivot_wider(names_from = location, values_from = proportion_within_tumorcells) %>% 
  mutate(diffv = metastasis - primary)

diff_plot <- ggplot(validation_data_metastasis_differences, aes(x = malignancy_rank, y = diffv)) +
  geom_bar(stat = 'identity', alpha=1.0) + 
  #geom_bar(data = validation_data_metastasis_for_plot %>% filter(location == 'primary'), aes(x = (malignancy_rank), y = proportion_within_tumorcells), stat = 'identity', color ='#8a8ad1' , alpha=0.2) +
  #geom_bar(data = validation_data_metastasis_for_plot %>% filter(location == 'metastasis'), aes(x = (malignancy_rank), y = proportion_within_tumorcells), stat = 'identity', color = '#ed9d9d', alpha=0.2) +
  facet_wrap(~Primary_ID, ncol=1, scales='free') +
  #scale_color_manual(values = c('#d22020', '#2525bf')) +
  scale_color_manual(values = gradient_colors[c(2,3)]) +
  theme_minimal() +
  theme(text=element_text(size=20), 
        legend.title = element_blank(), 
        legend.position = 'bottom') +
  xlab('Metastatic potential') +
  ylab('Proportion difference')

png('./figures/comparison_primary_vs_metastasis_malignancy.png', height=1500, width =2000, res=150)
plot_grid(primary_vs_met_plot, diff_plot)
dev.off()

pvalues_differences <- validation_data_metastasis_differences %>% 
  group_by(Primary_ID) %>% 
  dplyr::summarize(P = rcorr(diffv, malignancy_rank, type = 'spearman')$P[1,2], 
                    average_primary = sum(primary * malignancy_rank), 
                    average_metastasis = sum(metastasis  *malignancy_rank)) %>% 
  mutate(P_sci = formatC(P, format = "e", digits = 2))

pvalues_differences



############################
#validation on visium
###########################
visium_filenames <- c('VISIUM_filtered_expressionAK15_2563.csv',
                      'VISIUM_filtered_expressionAK15_6865.csv', 
                      'VISIUM_filtered_expressionAK16_4161.csv', 
                      'VISIUM_filtered_expressionAK17_1184.csv', 
                      'VISIUM_filtered_expressionAK17_15128.csv')

visium_data <- fread(paste0('./use_data/', visium_filenames[4])) %>%
  dplyr::select(V1, cluster_id, x_coord, y_coord) %>% 
  #filter(cluster_prob>0.3) %>%
  dplyr::mutate(cluster = paste0('X', cluster_id)) %>% 
  left_join(malignant_ranking)
  

png('./figures/visium_plot_malignancy_rank.png', width=4000, height=4000, res=100)
ggplot(visium_data, aes(x = x_coord, y = y_coord, color = as.numeric(malignancy_rank))) + 
  geom_point(size =0.1, alpha=1.0) +
  scale_color_viridis_c(option = 'C') +
  theme(
    panel.background = element_rect(fill = "black"),
    plot.background = element_rect(fill = "black"))
dev.off()

#############################
#Validation of protector rank
#############################






################################################
#compare malignancy rank to immune cell presence
#################################################
malignant_ranking

immune_cells_and_malignant_ranking <- patient_vector_long %>% 
  mutate(cluster = paste0('X', cluster)) %>%
  left_join(cell_types_fine_and_coarse) %>% 
  mutate(istumor = Tier_3 == 'Tumor') %>% 
  group_by(istumor, Pseudo) %>% 
  mutate(relative_cell_proportion = N/sum(N)) %>%
  #filter(Tier_3 == 'Tumor' | Tier_1 == 'Immune cells') %>%
  filter(Tier_3 == 'Tumor' | Tier_1 == 'Immune cells') %>%

  left_join(patient_meta_data %>% dplyr::select(Pseudo, Entity, UICC8_edition3)) %>%
  left_join(malignant_ranking) %>% 
  ungroup() 


immune_tumor_cancer_type <- 'AC'
immune_cell_proportions <- immune_cells_and_malignant_ranking %>% 
  #filter(Entity == immune_tumor_cancer_type) %>%
  #filter(UICC8_edition3=='IV') %>%
  filter(Tier_1 == 'Immune cells') %>% 
  dplyr::select(Pseudo, Entity, cluster, relative_cell_proportion) %>% 
  pivot_wider(names_from = cluster, values_from = relative_cell_proportion) 

patients_with_malignancy_ranking_correlation_analysis <- immune_cell_proportions %>% 
  dplyr::select(Pseudo) %>%
  ungroup() %>%
  left_join(patients_with_malignancy_ranking %>% ungroup() %>% 
  left_join(metastasized_yes_no %>% dplyr::select(Pseudo, malign)) %>%
  dplyr::select(Pseudo, malign) %>% ungroup())

immune_cell_proportions[is.na(immune_cell_proportions)] <- 0

correlation_malignancy_rank_immune <- immune_cell_proportions %>% 
  pivot_longer(!c(Pseudo, Entity), names_to = 'cluster', values_to = 'proportion') %>%
  left_join(patients_with_malignancy_ranking_correlation_analysis) %>%
  group_by(cluster) %>% 
  dplyr::summarize(correlation_with_malignancy_rank = rcorr(malign, proportion, type='spearman')$r[1,2], 
                  correlation_with_malignancy_rank_pvalue = rcorr(malign, proportion, type='spearman')$P[1,2]) %>%
  mutate(immune_clusters = cluster) %>% 
  left_join(cell_types_fine_and_coarse %>% dplyr::select(cluster, Tier_2, Tier_3), by = c('immune_clusters'  = 'cluster')) 

#View(correlation_malignancy_rank_immune)

ranking_immune_cells <- correlation_malignancy_rank_immune %>% 
  group_by(Tier_3) %>% 
  dplyr::summarize(average_correlation_with_malignancy_rank = mean(correlation_with_malignancy_rank))

correlation_malignancy_rank_immune_ordered <- correlation_malignancy_rank_immune %>% 
  left_join(ranking_immune_cells) %>% 
  arrange(average_correlation_with_malignancy_rank) %>%
  mutate(Tier_3 = factor(Tier_3, levels = unique(Tier_3)))

png(paste0('./figures/correlation_metastasis_immune_cells_', immune_tumor_cancer_type, '.png'), width=2000, height=2500, res=150)
ggplot(correlation_malignancy_rank_immune_ordered, aes(x = correlation_with_malignancy_rank, y = Tier_3, color = correlation_with_malignancy_rank_pvalue<0.05)) +
  geom_vline(xintercept = 0, linewidth=1, alpha=0.8)+
  geom_point(size=4) +
  #geom_smooth(se=F, method = 'lm') +
  #facet_wrap(~UICC8_edition3, ncol=1) +
  geom_hline(yintercept = 0, linewidth=1)+
  theme_minimal() +
  theme(legend.position = 'bottom', 
  text = element_text(size=20),
  axis.text.x = element_text(angle=0)) +
  scale_color_manual(values = c('#2e2ed8','#b12929'))
dev.off()

###############
#make protective scoring system
###############

protective_scores <- correlation_malignancy_rank_immune %>% 
  mutate(protective_value = ifelse(correlation_with_malignancy_rank>=0, 0,-1)) %>% 
  #mutate(protective_value =correlation_with_malignancy_rank) %>%
  filter(protective_value!=0) %>% 
  #filter(correlation_with_malignancy_rank_pvalue<0.05) %>% 
  dplyr::select(cluster, protective_value)


patients_with_protector_ranking <- patient_vector_long %>% 
  mutate(cluster = paste0('X', cluster)) %>%
  left_join(cell_types_fine_and_coarse) %>% 
  filter(Tier_1 == 'Immune cells' | Tier_3 == 'Tumor') %>% 
  group_by(Pseudo) %>% 
  mutate(immune_proportion = proportion / sum(proportion)) %>%
  inner_join(protective_scores) %>% 
  mutate(weighed_protective_values = immune_proportion * protective_value) %>%   #should the proportion be calculated only relative to tumor tissue or other?
  left_join(patient_meta_data %>% dplyr::select(Pseudo, UICC8_edition3, OS_m, Event, Entity)) %>%
  #dplyr::group_by(Pseudo, UICC8_edition3, Event, OS_m, Entity) %>%
  dplyr::group_by(Pseudo, UICC8_edition3, Event, OS_m) %>%
  dplyr::summarize(protector_rank = sum(weighed_protective_values)) %>% #weigh malignancy rank by proportion?
  left_join(patient_meta_data %>% dplyr::select(Pseudo, Grade)) %>% 
  left_join(patients_with_malignancy_ranking) %>% 
  mutate(low_malignancy = malignancy_rank < 160)


patients_with_protector_ranking %>% 
  group_by(Entity, UICC8_edition3) %>% 
  #group_by(UICC8_edition3) %>% 
  #group_by(Entity) %>% 
  dplyr::summarize(N = n(), 
                    c_index_protector = get_c_index(protector_rank, OS_m, Event), 
                    p_protector = get_c_index_p(protector_rank, OS_m, Event), 
                ) 





################
#correlation between immune cells:
################
dat_for_corr <- patient_vector_long %>% 
  mutate(cluster = paste0('X', cluster)) %>%
  left_join(cell_types_fine_and_coarse) %>% 
  filter(Tier_1 == 'Immune cells') %>% 
  mutate(annotated_cluster = paste0(Tier_3, '_', cluster)) %>%
  group_by(Tier_3, Pseudo) %>%
  dplyr::summarize(proportion = sum(proportion)) %>% 
  mutate(annotated_cluster = Tier_3) %>%
  ungroup() %>%
  dplyr::select(annotated_cluster, Pseudo, proportion) %>% 
  pivot_wider(names_from = annotated_cluster, values_from = proportion)

dat_for_corr[is.na(dat_for_corr)]<- 0

immune_corr <- cor(dat_for_corr[,-1]) %>% 
  as.data.frame() 
  
immune_corr_long <- immune_corr %>%
  mutate(annotated_cluster = rownames(.)) %>%
  pivot_longer(!annotated_cluster, names_to = 'annotated_cluster2', values_to = 'immune_r')

Heatmap(immune_corr %>% as.matrix())

################################################################
#mixed effects model to see if dendritic cells influence T cells  -FAIL
################################################################
  cell_types_fine_and_coarse$Tier_3 %>% unique()

dat_for_mixed_T_dendritic <- patient_vector_long %>% 
  mutate(cluster = paste0('X', cluster)) %>%
  left_join(cell_types_fine_and_coarse) %>% 
  filter(Tier_1 == 'Immune cells') %>% 
  mutate(annotated_cluster = Tier_3) %>%
  filter(annotated_cluster != 'T cells, cytotoxic, exhausted') %>%
    mutate(annotated_cluster = ifelse(startsWith(annotated_cluster, 'T cells, helper'), 'T cell', annotated_cluster)) %>% 
      mutate(annotated_cluster = ifelse(startsWith(annotated_cluster, 'Mast'), 'Dendritic', annotated_cluster)) %>% 
      group_by(annotated_cluster, Pseudo) %>%
  dplyr::summarize(proportion = sum(proportion)) %>% 
  ungroup() %>%
  dplyr::select(annotated_cluster, Pseudo, proportion) %>% 
  filter(startsWith(annotated_cluster, 'T cell') | startsWith(annotated_cluster, 'Dendritic')) %>% 
  filter((annotated_cluster == 'T cell') | (annotated_cluster == 'Dendritic')) %>% 

  pivot_wider(names_from = annotated_cluster, values_from = proportion)

dat_for_mixed_T_dendritic[is.na(dat_for_mixed_T_dendritic)] <- 0


dat_for_mixed_T_dendritic_with_survival <- dat_for_mixed_T_dendritic %>% 
  left_join(patient_meta_data %>% dplyr::select(Pseudo, UICC8_edition3, OS_m, Event, Entity)) #%>% 
  #filter(Entity == 'AC')

colnames(dat_for_mixed_T_dendritic)

Tcell_dendritic_cox_model <- coxph(Surv(OS_m, Event) ~ UICC8_edition3 + `T cell` + Dendritic , data = dat_for_mixed_T_dendritic_with_survival)

summary(Tcell_dendritic_cox_model)

####################################
#compare Immune cells to staging? --probably FAIL
######################################

immune_cells_per_stage <- patient_vector %>%
  pivot_longer(!Pseudo,, names_to = 'cluster', values_to = 'proportion') %>%
  mutate(cluster = paste0('X', cluster)) %>%
  left_join(cell_types_fine_and_coarse) %>% 
  #filter(Tier_1 == 'Mesenchymal cells') %>% 
  group_by(Tier_3, Pseudo) %>%
  dplyr::summarize(proportion = sum(proportion)) %>% 
  ungroup() %>%
  dplyr::select(Tier_3, Pseudo, proportion) %>% 
  left_join(patient_meta_data %>% dplyr::select(Pseudo, UICC8_edition3, OS_m, Event, Entity)) %>% 
  left_join(metastasized_yes_no %>% dplyr::select(Pseudo, malign))

immune_cells_per_stage_pvalue <- immune_cells_per_stage %>% 
  group_by(Tier_3) %>% 
  dplyr::summarize(r=rcorr(as.numeric(as.factor(UICC8_edition3)), proportion, type = 'spearman')$r[1,2], 
    P = rcorr(as.numeric(as.factor(UICC8_edition3)), proportion)$P[1,2], type = 'spearman') %>% 
  ungroup() %>% 
  mutate(p_adj = p.adjust(P, method = 'fdr')) %>% 
  filter(p_adj<0.01) %>% 
  dplyr::select(Tier_3) %>% 
  unique()


immune_cells_per_stage_summarized <- immune_cells_per_stage %>% 
  group_by(Tier_3, UICC8_edition3, Entity) %>% 
  dplyr::summarize(proportion = mean(proportion)) %>% 
  inner_join(immune_cells_per_stage_pvalue)

immune_cells_per_stage_significant <-immune_cells_per_stage %>% 
  inner_join(immune_cells_per_stage_pvalue)

some_colors <- plasma(5)

ggplot(immune_cells_per_stage_summarized, aes(x = proportion, y = Tier_3, color = UICC8_edition3)) +
  facet_wrap(~Entity) +
  geom_point(size=5, alpha=0.6) +
  theme_bw() +
  scale_x_continuous(trans = 'log10') +
  scale_color_manual(values = some_colors[1:4]) +
  theme(text = element_text(size=15))


ggplot(immune_cells_per_stage_significant, aes(x = proportion, y = Tier_3, fill = UICC8_edition3)) +
  facet_wrap(~Entity) +
  geom_boxplot() +
  theme_bw() +
  scale_x_continuous(trans = 'log10') +
  scale_fill_manual(values = some_colors[1:4]) +
  theme(text = element_text(size=15))




