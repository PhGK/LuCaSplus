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
library(readxl)
library(compareC)


setwd('/mnt/ssd/shared/LungCAIRE_SingleCell/plots_statistics')
cancer_type <- 'SCC'

####
## load data
####

fill_left_zeros <- function(dummy_matrix) {
  # Apply the transformation to each row
  filled_matrix <- t(apply(dummy_matrix, 1, function(row) {
    for (i in seq_along(row)) {
      if (row[i] == 1) {
        row[1:i] <- 1
      }
    }
    return(row)
  }))
  
  return(as.data.frame(filled_matrix))
}

match_ids <- read.csv('/mnt/ssd/shared/LungCAIRE_SingleCell/data/LungCAIRE_snRNAseq_cohort.csv') %>% 
  dplyr::select(Match_ID, Pseudo) %>% mutate(ENR = Match_ID)

patient_meta_data <- read.csv('/mnt/ssd/shared/LungCAIRE_SingleCell/data/Metadata_LungCAIRE_retrospective_Berlin_Cologne_20250303.csv') %>% 
  mutate(UICC8_edition2 = UICC8_EDITION2, OS_m = OS_M, Event = EVENT,
            age = DIA_AGE_Y, Sex = SEX, Entity = ENTITY, 
            Adj_therapy = ADJ_THERAPY, Grade = GRADE) %>% 
  left_join(match_ids) %>% 
  filter(Pseudo != 'LC125')

clean_cell_annotation <- read.csv('./use_data/VICREG_clustering_highinputdropout_5000epochs_rev_annotation.csv') %>% 
  dplyr::select(cluster = Cluster_ID, Tier_1, Tier_2, Tier_3)


relapse_data <- read.csv('./use_data/time2relapse.csv')

##########
# for manuscript
########
patient_malignancy_ranks <- read.csv('./use_data/patient_metastatic_potential.csv') %>% 
    dplyr::select(Pseudo, malignancy_rank)


staging_and_survival <- patient_meta_data%>% 
    unique() %>% 
    mutate(UICC8_edition3 = ifelse(UICC8_edition2 %in% c('IA', 'IB'), 'I', UICC8_edition2)) %>%
    mutate(UICC8_edition3 = factor(UICC8_edition3, levels = c('I', 'II', 'III', 'IV'), ordered=T)) %>% 
    dplyr::select(-ADJ_CHEMO_DRUGS) %>%
    left_join(patient_malignancy_ranks) 


###############
#prepare clinical data
###############

stage_dummies <- model.matrix(~ UICC8_edition3 - 1, data = staging_and_survival) %>% as.data.frame() %>%
        fill_left_zeros()

entity_dummies <- staging_and_survival %>% dplyr::select(Entity_ac = Entity) %>% mutate(Entity_ac = ifelse(Entity_ac == 'AC',1,0))
adj_dummies <- staging_and_survival %>% dplyr::select(adj_therapy = Adj_therapy) %>% mutate(adj_therapy = ifelse(adj_therapy == 'yes',1,0)) %>% 
  mutate(adj_therapy = ifelse(is.na(adj_therapy),0,adj_therapy))


all_clinical_data <- cbind(stage_dummies, entity_dummies) %>% #, age_dummies) %>% 
    mutate(Pseudo = staging_and_survival$Pseudo)

all_clinical_data_with_survival <- all_clinical_data%>%
        mutate(OS_m =  staging_and_survival$OS_m, Event = staging_and_survival$Event) 



##################
#prepare single cell profiles
##################

annotated_clusters_FINAL <- read.csv('./use_data/VICREG_clustering_highinputdropout_5000epochs_rev_annotation.csv') %>%
  mutate(Tier_2 = ifelse(Tier_3 == 'Tumor', 'Tumor', Tier_2))
clusters_not_mixed <- annotated_clusters_FINAL  %>% filter(Tier_3 != 'Mixed') 

cluster_result <- read_parquet('../scGPT_embeddings/clusterings/VICREG_clustering_highinputdropout_5000epochs.parquet') %>%
  mutate(cell_position = seq(dim(.)[1])) %>% 
  dplyr::select(-c(Tier_1, Tier_2, Tier_3, Tier_4)) %>% 
  left_join(clean_cell_annotation)

patient_vector_long <- cluster_result %>% 
  group_by(cluster, Pseudo) %>%
  dplyr::summarize(N=n(), cluster = first(cluster))  %>% 
  group_by(Pseudo) %>% 
  mutate(proportion = N/sum(N)) %>% 
  dplyr::select(cluster, proportion, Pseudo) %>% 
  ungroup() %>% 
  filter(cluster %in% clusters_not_mixed$Cluster_ID) 
  
 patient_vector <- patient_vector_long %>%
  dplyr::select(Pseudo, cluster, proportion) %>%
  group_by(cluster) %>% 
  pivot_wider(names_from = cluster, values_from = proportion)

patient_vector[is.na(patient_vector)] <- 0

patients_with_snRNA <- patient_vector_long$Pseudo %>% unique()


#########################
#### filter clinical data for patients with existing snRNA
##########################
filtered_all_clinical_data <- all_clinical_data %>% 
  filter(Pseudo %in% patients_with_snRNA) 

filtered_staging_and_survival <- staging_and_survival %>% 
  filter(Pseudo %in% patients_with_snRNA) 

all_clinical_data_with_tumor <- filtered_all_clinical_data %>%
  left_join(patient_vector) %>% 
  dplyr::select(-Pseudo) 

#####################

only_tumor_data <- all_clinical_data_with_tumor %>% select(matches("^\\d+$"))
only_heterogeneity <- filtered_all_clinical_data %>% 
  dplyr::select(Pseudo) %>% 
  left_join(cluster_per_patient) %>% 
  dplyr::select(-Pseudo)

dim(all_clinical_data_with_tumor)
dim(only_tumor_data)

#determine hyperparameter
general_alpha <- 0.1 #was 0.1

test_coef_all <- function(X,surv_time, surv_event, indices, mode = 'coefs') {
  alpha=general_alpha

  y <- Surv(surv_time, surv_event)
  best_lambda1 <- cv.glmnet(X %>% as.matrix(), y, family = 'cox', alpha = alpha)$lambda.min
  best_lambda2 <- cv.glmnet(X %>% as.matrix(), y, family = 'cox', alpha = alpha)$lambda.min
  best_lambda3 <- cv.glmnet(X %>% as.matrix(), y, family = 'cox', alpha = alpha)$lambda.min

  best_lambda <- mean(best_lambda1,best_lambda2,best_lambda3)

  model <- glmnet(X %>% as.matrix(), y, family = 'cox', alpha=alpha, lambda=best_lambda)

  coefs <- data.frame(x = coef(model) %>% as.matrix) %>% t()
  data.frame(coefs)
  }

test_predictions<- function(X,patient_ids, surv_time, surv_event, indices) {
  alpha=general_alpha

  X_train <- X[-indices,]
  X_test <- X[indices,]

  surv_time_train <- surv_time[-indices]
  surv_time_test <- surv_time[indices]

  surv_event_train <- surv_event[-indices]
  surv_event_test <- surv_event[indices]

  y_train <- Surv(surv_time_train, surv_event_train)
  
  #train_data <- cbind(surv_time_train, surv_event_train, X_train)

  best_lambda1 <- cv.glmnet(X_train %>% as.matrix(), y_train, family = 'cox', alpha = alpha)$lambda.min
  best_lambda2 <- cv.glmnet(X_train %>% as.matrix(), y_train, family = 'cox', alpha = alpha)$lambda.min
  best_lambda3 <- cv.glmnet(X_train %>% as.matrix(), y_train, family = 'cox', alpha = alpha)$lambda.min

  best_lambda <- mean(best_lambda1,best_lambda2,best_lambda3)
  print(best_lambda)

  model <- glmnet(X_train %>% as.matrix(), y_train, family = 'cox', alpha=alpha, lambda=best_lambda)
  lp <- predict(model, newx = X_test %>% as.matrix(), s = best_lambda, type = 'link') %>% as.vector()

  data.frame('pred' = lp, surv_time_test, surv_event_test, patient_ids = patient_ids[indices])
  }


##
#last filtering of train and test data

adj_mask <- filtered_staging_and_survival$Adj_therapy != 'yes'
adj_mask[is.na(adj_mask)] <- T

entity_mask <- (filtered_all_clinical_data$Entity_ac == ifelse(cancer_type == 'SCC', 0,1))
mask <- entity_mask & adj_mask 

#generate seeds for LOO-CV
set.seed(123)
folds <- caret::createFolds(filtered_staging_and_survival$UICC8_edition3[mask], k=sum(mask), returnTrain = F)
folds

############
#Kaplan Meier for PJ
############

calculate_median_survival <- function(survival_time, death_event) {
  # Create a survival object
  surv_obj <- Surv(survival_time, death_event)
  # Fit the Kaplan-Meier estimator
  fit <- survfit(surv_obj ~ 1)
  # Extract and return the median survival time
  return(summary(fit)$table["median"])
}



#execute LOO-CV for clinical, snRNA, and clinical # snRNA
all_clinical_predictions <- rbindlist(lapply(folds, function(x) test_predictions(filtered_all_clinical_data[mask,]%>% dplyr::select(-Pseudo),filtered_staging_and_survival$Pseudo[mask], filtered_staging_and_survival$OS_m[mask], filtered_staging_and_survival$Event[mask], x)))
tumor_cell_predictions_no_clinical <- rbindlist(pbmclapply(folds, function(x) test_predictions(only_tumor_data[mask,],filtered_staging_and_survival$Pseudo[mask], filtered_staging_and_survival$OS_m[mask], filtered_staging_and_survival$Event[mask], x), mc.cores=20))
tumor_cell_predictions <- rbindlist(pbmclapply(folds, function(x) test_predictions(all_clinical_data_with_tumor[mask,],filtered_staging_and_survival$Pseudo[mask], filtered_staging_and_survival$OS_m[mask], filtered_staging_and_survival$Event[mask], x), mc.cores=20))


tnm_predictions <- filtered_staging_and_survival[mask,] %>% 
  mutate(surv_time_test = filtered_staging_and_survival$OS_m[mask], surv_event_test = filtered_staging_and_survival$Event[mask]) %>% 
    mutate(q = as.numeric(UICC8_edition3)) %>% 
    mutate(patient_ids = filtered_staging_and_survival$Pseudo[mask])


#combine all predictions and UICC staging in one dataframe
all_predictions <- tumor_cell_predictions_no_clinical %>% 
    dplyr::select(Pseudo = patient_ids, only_tumor_pred = pred, survival_time = surv_time_test, survival_event = surv_event_test) %>% 
  left_join(tumor_cell_predictions %>% dplyr::select(Pseudo =patient_ids, full_model_pred = pred)) %>% 
  left_join(all_clinical_predictions %>% dplyr::select(Pseudo =patient_ids, clinical_model_pred = pred)) %>% 
  left_join(tnm_predictions %>% dplyr::select(Pseudo =patient_ids, UICC = UICC8_edition3, surv_time_UICC = surv_time_test, surv_event_UICC = surv_event_test)) %>%
  mutate(UICC = as.numeric(UICC)) 
  
#not sure if this is used anywhere
UICC_survival_times <- all_predictions %>% group_by(UICC) %>% 
  dplyr::summarize(UICC_median_survival = calculate_median_survival(survival_time, survival_event))


###calculate nquantile risk groups with similar risk groups for KM visualization and future clinical application###
nquantiles<-4

min_max_prediction_values <- all_predictions %>% 
  dplyr::summarize(minv = quantile(full_model_pred, 0.02), maxv = quantile(full_model_pred, 0.98), diffv = maxv-minv) 
equal_sized_q <- data.frame(data.frame(cutoffs = (seq(nquantiles)-1)) * (min_max_prediction_values$diffv/nquantiles)+ min_max_prediction_values$minv) %>% 
  t() 

colnames(equal_sized_q) = c('A', 'B', 'C', 'D')
write.csv(equal_sized_q, paste0('./use_data/q_cutoffs_', cancer_type, '.csv'))

all_predictions_q <- all_predictions %>% 
          cbind(equal_sized_q) %>% 
          mutate(q = (full_model_pred>= (-Inf)) + (full_model_pred>=B) + (full_model_pred>=C) + (full_model_pred>=D))


direct_comparison_internal <- all_predictions_q %>% 
              dplyr::summarize(Pvalue = compareC(survival_time, survival_event, full_model_pred, UICC)$pval , 
                               Pvalue_compareC_alternative = cindex.comp( concordance.index(x = full_model_pred, surv.time = survival_time, surv.event = survival_event, outx=F, method = 'noether'), 
                                                                          concordance.index(UICC, surv.time = survival_time, surv.event = survival_event, outx = F, method = 'noether'))$p.value,)
direct_comparison_internal

all_predictions_long <- all_predictions %>% pivot_longer(!c(Pseudo, survival_time, survival_event), names_to = 'model', values_to = 'prediction') 
all_predictions_long_by_UICC <- all_predictions %>% pivot_longer(!c(Pseudo, survival_time, survival_event, UICC), names_to = 'model', values_to = 'prediction') 

all_predictions_long %>% 
  group_by(model) %>% 
  dplyr::summarize(c_index = concordance.index(x = prediction, surv.time = survival_time, surv.event = survival_event, outx=F)$c.index, 
                  p = concordance.index(x = prediction, surv.time = survival_time, surv.event = survival_event, outx=F, method = 'noether')$p.value, 
                  N=n()) 

all_predictions_long_by_UICC_c_indices <- all_predictions_long_by_UICC %>% 
  group_by(model, UICC) %>% 
  dplyr::summarize(c_index = concordance.index(x = prediction, surv.time = survival_time, surv.event = survival_event, outx=F)$c.index, N=n())




###############################
#Visualize prediction results in KM curves
###############################


####
# 💩 code by PJ #
####

#Prepare data for UICC and Full model
# Cap survival times at 120 months
tnm_predictions <- tnm_predictions %>%
  mutate(survival_time_capped = pmin(surv_time_test, 120))
surv_fit_TNM <- survfit(Surv(survival_time_capped, surv_event_test) ~ UICC8_edition3, data = tnm_predictions)

# Cap survival times at 120 months
all_predictions_q <- all_predictions_q %>%
  mutate(survival_time_capped = pmin(survival_time, 120))
surv_fit_tumor_cell <- survfit(Surv(survival_time_capped, survival_event) ~ q, data = all_predictions_q)


# Generate the plot
tnm <- ggsurvplot(
  surv_fit_TNM,
  data = tnm_predictions,
  pval = T,
  pval.method = TRUE,
  pval.size = 4,
  conf.int = FALSE,
  censor = TRUE,
  legend.title = "UICC stage",
  title = "UICC8",
  xlab = "Time (Months)",
  ylab = "Survival Probability",
  palette = c("#11B5AE", "#4046C9", "#F68512", "#DE3C82"),
  linewidth = 0.1,
  ggtheme = theme_minimal(),
  lwd = 0.5,
  risk.table = TRUE,       # Add the risk table below the plot
  risk.table.title = "Number at Risk", # Title for the risk table
  risk.table.height = 0.25 # Adjust height of the risk table
)

tnm$table <- tnm$table + 
  scale_y_discrete(labels = c("IV", "III", "II", "I")) # Customize y-axis labels


# Generate the plot
new <- ggsurvplot(
  surv_fit_tumor_cell,
  data = all_predictions_q,
  pval = T,
  pval.method = TRUE,
  pval.size = 4,
  conf.int = FALSE,
  censor = TRUE,
  legend.title = "Risk quantile",
  title = "Full model",
  xlab = "Time (Months)",
  ylab = "Survival Probability",
  palette = c("#11B5AE", "#4046C9", "#F68512", "#DE3C82"),
  linewidth = 0.1,
  ggtheme = theme_minimal(),
  lwd = 0.5,
  risk.table = TRUE,       # Add the risk table below the plot
  risk.table.title = "Number at Risk", # Title for the risk table
  risk.table.height = 0.25 # Adjust height of the risk table
)

new$table <- new$table + 
  scale_y_discrete(labels = rev(c("1", "2", "3", "4"))) # Customize y-axis labels



pdf(paste0('./figures/KM_',cancer_type,'_tnm_with_risk_table.pdf'), width = 5, height = 6, onefile = FALSE)
tnm
dev.off()

pdf(paste0('./figures/KM_',cancer_type,'_snRNAseq_with_risk_table.pdf'), width = 5, height = 6, onefile = FALSE)
new
dev.off()






gradient_colors <- c("#11B5AE", "#4046C9", "#F68512", "#DE3C82", "#11B5AE", "#4046C9", "#F68512", "#DE3C82")

########
#new survival stratification by STAGE
########

tumor_cell_predictions_with_staging <- tumor_cell_predictions %>% 
  left_join(filtered_staging_and_survival %>% dplyr::select(patient_ids = Pseudo, UICC8_edition3))

#surv_fit_snRNA_with_staging <- survfit(Surv(surv_time_test, surv_event_test) ~ q, data = tumor_cell_predictions_with_staging)
surv_fit_snRNA_with_staging <- survfit(Surv(survival_time, survival_event) ~ q, data = all_predictions_q)

png(paste0('./figures/KM_by_UICC_',cancer_type,'_MP.png'), width=1500, height=1500, res=250)
ggsurvplot_facet(
  surv_fit_tumor_cell,
  data = all_predictions_q,
  facet.by = c('UICC'),
  #facet.by = c('Entity'),
  pval = F,           # Add p-value
  pval.method = TRUE,    # Show test method in the p-value plot
  pval.size = 4,         # Set the p-value font size
  conf.int = F,       # Add confidence intervals
  censor = TRUE,         # Show censored data points
  legend.title = "Risk quantile",
  title = cancer_type,
  xlab = "Time (Months)",
  ylab = "Survival Probability",
  palette = gradient_colors,
  linewidth=0.1,
  xlim = c(0,60),
  ggtheme=theme_minimal()
)
dev.off()



#############
#sankey plots
############
calculate_5yr_survival <- function(survival_time, death_event) {
  # Convert inputs to numeric
  survival_time <- as.numeric(survival_time)
  death_event <- as.numeric(death_event)
  
  # Create a survival object
  surv_obj <- Surv(survival_time, death_event)
  
  # Fit the Kaplan-Meier estimator
  fit <- survfit(surv_obj ~ 1)
  
  # Determine the proper time point:
  time_point <- ifelse(max(survival_time, na.rm = TRUE) < 10, 5, 60)
  
  # Get the survival summary at the chosen time point
  surv_summary <- summary(fit, times = time_point)
  
  # If no survival estimate is returned because all patients died at 5yr, return 1
  surv_5yr <- if (length(surv_summary$surv) == 0) 0 else surv_summary$surv
  
  return(surv_5yr)
}




median_survival_for_sankey <- all_predictions_q %>% 
  dplyr::select(Pseudo, survival_time, survival_event, UICC, q) %>% 
  group_by(q) %>% 
  dplyr::mutate(new_median_survival = as.character(round(calculate_5yr_survival(survival_time, survival_event),3))) %>% 
  group_by(UICC) %>% 
  dplyr::mutate(old_median_survival = as.character(round(calculate_5yr_survival(survival_time, survival_event),3))) %>% 
  mutate(old_median_survival = ifelse(is.na(old_median_survival), 'Not reached', old_median_survival)) %>%
  mutate(new_median_survival = ifelse(is.na(new_median_survival), 'Not reached', new_median_survival)) %>%
  mutate(q = paste0(q, ' (', new_median_survival, ')'), UICC = paste0(UICC, ' (', old_median_survival, ')')) %>% 
  ungroup() 

restratification_numbers <- median_survival_for_sankey %>% 
  group_by(UICC, q) %>%
  dplyr::summarize(N = n()) %>% 
  group_by(UICC) %>% mutate(UICC_numbers = sum(N)) %>%
  group_by(q) %>% mutate(q_numbers = sum(N))

dat_for_sankey_wide <- median_survival_for_sankey %>% 
  dplyr::select(UICC, snRNA = q, Pseudo) 
  
dat_for_sankey <- dat_for_sankey_wide%>% 
  make_long(!Pseudo)

write.csv(dat_for_sankey_wide, paste0('./use_data/sankey_data_',cancer_type, '.csv'))


png(paste0('./figures/sankey_',cancer_type,'_MP.png'), width=1800, height=1500, res=220)
ggplot(dat_for_sankey, aes(x = x,
               next_x = next_x,
               node = node,
               next_node = next_node,
               fill = factor(node),
                label = node)) +
    geom_sankey(flow.alpha = 0.5,          #This Creates the transparency of your node
    node.color = "black",     # This is your node color
    show.legend = F) +
    theme_minimal() +
    theme(axis.text.y = element_blank(),
        axis.text.x = element_text(size=15),
        axis.title = element_blank()) +
        geom_sankey_label(size = 6, fill = 'white', alpha=0.7) +
    #scale_fill_manual(values = c('#0d0887','#fad824',   '#f1814d',    '#d45270', '#ac2694'))
    scale_fill_manual(values = c(gradient_colors,gradient_colors))
dev.off()

pdf(paste0('./figures/sankey_',cancer_type,'_MP.pdf'), width=12, height=10, onefile = FALSE)
ggplot(dat_for_sankey, aes(x = x,
               next_x = next_x,
               node = node,
               next_node = next_node,
               fill = factor(node),
                label = node)) +
    geom_sankey(flow.alpha = 0.5,          #This Creates the transparency of your node
    node.color = "black",     # This is your node color
    show.legend = F) +
    theme_minimal() +
    theme(axis.text.y = element_blank(),
        axis.text.x = element_text(size=15),
        axis.title = element_blank()) +
        geom_sankey_label(size = 6, fill = 'white', alpha=0.7) +
    #scale_fill_manual(values = c('#0d0887','#fad824',   '#f1814d',    '#d45270', '#ac2694'))
    scale_fill_manual(values = c("#11B5AE", "#4046C9", "#F68512", "#DE3C82", "#11B5AE", "#4046C9", "#F68512", "#DE3C82"))
dev.off()

################
#input_relevances
################

all_clinical_coefficients <- test_coef_all(filtered_all_clinical_data[mask,]%>% dplyr::select(-Pseudo), filtered_staging_and_survival$OS_m[mask], filtered_staging_and_survival$Event[mask]) %>% 
  mutate(index = rownames(.)) %>% 
  pivot_longer(!index, names_to = 'cluster', values_to = 'coefficient') %>% 
  filter(coefficient!=0)

tumor_cell_coefficients <- test_coef_all(all_clinical_data_with_tumor[mask,], filtered_staging_and_survival$OS_m[mask], filtered_staging_and_survival$Event[mask]) %>%
    mutate(index = rownames(.)) %>% 
  pivot_longer(!index, names_to = 'cluster', values_to = 'coefficient') %>% 
    filter(coefficient!=0)

#############
#determine cluster cell types
############

cell_types_fine_and_coarse <- annotated_clusters_FINAL %>% 
  mutate(max_fine_Tier = Tier_3, max_coarse_Tier = Tier_2) %>%
  mutate(cluster = paste0('X', Cluster_ID)) %>% 
  dplyr::select(-Cluster_ID)

how_many_cell_types <- cell_types_fine_and_coarse %>% group_by(Tier_1) %>% 
  dplyr::summarize(N=n())
how_many_cell_types

tumor_cell_coefficients_annotated <-  tumor_cell_coefficients %>% 
  left_join(cell_types_fine_and_coarse) %>% 
  mutate(annotated_cluster_coarse = ifelse(is.na(max_coarse_Tier), cluster, max_coarse_Tier)) %>% 
  group_by(annotated_cluster_coarse) %>% 
  mutate(cluster_rank = as.character(rank(cluster)))


#selected vs possible cluster numbers
N_all_cell_types <- cell_types_fine_and_coarse %>% 
  dplyr::select(Tier_2, Tier_3) %>% unique() %>% mutate(type_ = 'all')

N_selected_cell_types <- tumor_cell_coefficients_annotated %>% 
  ungroup() %>%
  filter(startsWith(cluster, 'X')) %>% 
  dplyr::select(Tier_2, Tier_3) %>% unique() %>%
  mutate(type = 'selected') %>% 
  right_join(N_all_cell_types)

  
write.csv(tumor_cell_coefficients_annotated %>% mutate(cancer_type = cancer_type), paste0('./use_data/important_clusters_', cancer_type, '_MP.csv'))
write.csv(N_selected_cell_types %>% mutate(cancer_type = cancer_type), paste0('./use_data/important_celltypes_', cancer_type, '_MP.csv'))

write.csv(tumor_cell_predictions, paste0('./use_data/patient_outcome_predictions_full_model', cancer_type, '_MP.csv'))



###########################
#plot coefficients
###########################
# Load your data
important_clusters_AC <- read.csv("/mnt/ssd/shared/LungCAIRE_SingleCell/plots_statistics/use_data/important_clusters_AC_MP_fixed.csv")
important_clusters_SCC <- read.csv("/mnt/ssd/shared/LungCAIRE_SingleCell/plots_statistics/use_data/important_clusters_SCC_MP_fixed.csv")
annotation_data <- read.csv("/mnt/ssd/shared/LungCAIRE_SingleCell/plots_statistics/use_data/VICREG_clustering_highinputdropout_5000epochs_rev_annotation.csv") %>% 
  dplyr::select(-c(Freq))

# Remove the 'X' prefix from cluster columns in both AC and SCC datasets
important_clusters_AC$cluster <- gsub("X", "", important_clusters_AC$cluster)
important_clusters_SCC$cluster <- gsub("X", "", important_clusters_SCC$cluster)

# Convert the cluster column to numeric to allow proper merging
important_clusters_AC$cluster <- as.numeric(important_clusters_AC$cluster)
important_clusters_SCC$cluster <- as.numeric(important_clusters_SCC$cluster)

# Merge the important_clusters with annotation data
merged_AC <- merge(important_clusters_AC, annotation_data, by.x = "cluster", by.y = "Cluster_ID", all.x = TRUE)
merged_SCC <- merge(important_clusters_SCC, annotation_data, by.x = "cluster", by.y = "Cluster_ID", all.x = TRUE)

# Add a column to indicate tumor type
merged_AC$tumor_type <- "AC"
merged_SCC$tumor_type <- "SCC"

# Combine both datasets
combined_data <- bind_rows(merged_AC, merged_SCC)

combined_data$coefficient_asinh <- asinh(combined_data$coefficient)

filtered_data <- combined_data %>%
  filter(!is.na(Tier_1.y) & Tier_2.y != "Normal epithelial cells" & Tier_1.y != "Mixed" & Tier_3.y != "Skeletal muscle cells")


# Now, the ggplot code
#pdf("Dotplot_prognosis_stromal_immune_cells.pdf", height = 8, width =8)
pdf("./figures/Dotplot_prognosis_stromal_immune_cells.pdf", width = 8, height = 8)
ggplot(filtered_data, aes(x = coefficient_asinh, y = Tier_3.x, color = coefficient_asinh, size = Freq)) +
  geom_point(alpha = 1) +  # Use geom_point() without jitter
  facet_grid(Tier_2.y ~ tumor_type, scales = "free_y", space = "free_y") +  # Allow free y-space for better fitting
  geom_vline(xintercept = 0, color = "#0A3659", linetype = "solid") +  # Add vertical line at x = 0
  scale_color_gradientn(colors = c("#15607A", "#C7CDD1", "#A63717"), limits = c(-6, 6), oob = scales::squish) +  # Custom 3-color palette
  scale_y_discrete() +  # Now that Tier_3_rev.y is a factor, no need to specify limits again
  #scale_x_continuous(limits = c(-2, 1)) +  # Adjust the x-axis range
  theme_minimal() +
  labs(x = "Asinh(Coefficient)",
       y = "Tier 3 Annotation",
       color = "Coefficient") +  # Color legend label
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.spacing.x = unit(3.5, "lines"))  # Adjust the spacing between panels
dev.off()


######################################
#Certain patients with interesting risk restratification
#####################################
calc_wilcox <- function(y,x) {
  if (length(unique(x))!=2) {return(NA)}

  return(wilcox.test(y ~ x)$p.value)
  }

tumor_cell_coefficients_annotated

relevant_switchers <- all_predictions_q %>% 
  dplyr::select(Pseudo, UICC, q) 

switch4_1 <- relevant_switchers %>% 
  filter(UICC==4,q!=2,  q!=3) %>%
  mutate(risk_strat = ifelse(q==4, '4_4', '4_1')) %>% 
  left_join(patient_vector_long) %>% 
  mutate(cluster = paste0('X', cluster)) %>% 
  inner_join(tumor_cell_coefficients_annotated) %>% 
  mutate(contribution = proportion * coefficient) %>% 
  mutate(cluster_Tier_3 = paste0(Tier_2, cluster)) %>% 
  mutate(Tier_3_tumor = ifelse(Tier_2 == 'Tumor', cluster_Tier_3, Tier_3)) %>%
  group_by(Tier_2, Tier_3, Tier_3_tumor, Pseudo, risk_strat) %>% 
  dplyr::summarize(contribution = sum(contribution)) %>% 
  group_by(Tier_3_tumor) %>% 
  #mutate(P_wilcox = wilcox.test(contribution ~ risk_strat)$p.value) %>% 
  mutate(P_wilcox = calc_wilcox(contribution,risk_strat)) %>% 
  ungroup() %>% 
  filter(Tier_2 != 'Normal epithelial cells')

patient_meta_data %>% dplyr::select(Pseudo, OS_M, Event) %>% 
  filter(Pseudo %in% (switch4_1 %>% filter(risk_strat =='4_1') %>% .$Pseudo))

switch4_1_p <- switch4_1 %>% 
  group_by(Tier_2, Tier_3_tumor, P_wilcox) %>% 
  unique() %>% 
  mutate(p_adj = signif(p.adjust(P_wilcox, method = 'fdr'),3), p_str = ifelse(p_adj<0.05, as.character(p_adj), ''), 
  p_star = ifelse(p_adj<0.001, '***', ifelse(p_adj<0.01, '**', '*'))) %>% 
  mutate(Restratification = ifelse(risk_strat == '4_1', 'Stage IV to FM1', 'Stage IV to FM4'))



png(paste0('./figures/switch_4_1_',cancer_type,'.png'), width=2000, height=2000, res=200)
ggplot(switch4_1_p %>% filter(p_adj < 0.05), aes(y = Tier_3, x= asinh(contribution), fill = risk_strat)) +
  geom_boxplot() +
  geom_text(data = switch4_1_p %>% filter(p_adj<0.05), aes(y = Tier_3, x =0.8, label = p_str)) +
  facet_grid(Tier_2~1, scales='free_y', space='free_y') +
  theme_minimal() +
  scale_fill_manual(values = c("#15607A", "#A63717"))
dev.off()

pdf(paste0('./figures/switch_4_1_',cancer_type,'.pdf'), width=10, height=6, onefile=T)
ggplot(switch4_1_p %>% filter(p_adj < 0.05), aes(y = Tier_3, x= asinh(contribution), fill = Restratification)) +
  geom_point(aes(color = Restratification), position = position_dodge(0.8), size=1) +
  geom_boxplot(alpha=0.1, coef=Inf) +
  geom_text(data = switch4_1_p %>% filter(p_adj<0.05), aes(y = Tier_3, x =0.8, label = p_star)) +
  facet_grid(Tier_2~'', scales='free_y', space='free_y') +
  theme_minimal() +
  scale_fill_manual(values = c("#15607A", "#A63717")) +
  scale_color_manual(values = c("#15607A", "#A63717")) +
  theme(axis.title.y = element_blank(), 
        panel.spacing = unit(1,'cm'), 
        strip.text.y = element_text(angle=0))
dev.off()

switch1_34 <- relevant_switchers %>% 
  filter(UICC==1, q!=2) %>%
  mutate(risk_strat = ifelse(q>=3, '1_34', '1_1')) %>% 
  left_join(patient_vector_long) %>% 
  mutate(cluster = paste0('X', cluster)) %>% 
  inner_join(tumor_cell_coefficients_annotated) %>% 
  mutate(contribution = proportion * coefficient) %>% 
  mutate(cluster_Tier_3 = paste0(Tier_2, cluster)) %>% 
  mutate(Tier_3_tumor = ifelse(Tier_2 == 'Tumor', cluster_Tier_3, Tier_3)) %>%
  #mutate(Tier_3_tumor = ifelse(Tier_3 == 'Tumor', cluster_Tier_3, cluster_Tier_3)) %>%
  group_by(Tier_2, Tier_3_tumor,Tier_3,  Pseudo, risk_strat) %>% 
  dplyr::summarize(contribution = sum(contribution)) %>% 
  group_by(Tier_3_tumor) %>% 
  #mutate(P_wilcox = wilcox.test(contribution ~ risk_strat)$p.value) %>% 
   mutate(P_wilcox = calc_wilcox(contribution,risk_strat)) %>% 
  ungroup() %>% 
    filter(Tier_2 != 'Normal epithelial cells')


switch1_34_p <- switch1_34 %>% 
  group_by(Tier_2, Tier_3_tumor, P_wilcox) %>% 
  unique() %>% 
  mutate(p_adj = signif(p.adjust(P_wilcox, method = 'fdr'),3), p_str = ifelse(p_adj<0.05, as.character(p_adj), ''), 
        p_star = ifelse(p_adj<0.001, '***', ifelse(p_adj<0.01, '**', '*'))) %>% 
    mutate(Restratification = ifelse(risk_strat == '1_34', 'Stage I to FM3/4', 'Stage I to FM1'))



png(paste0('./figures/switch_1_34_',cancer_type,'.png'), width=2000, height=2000, res=200)
ggplot(switch1_34_p %>% filter(p_adj < 0.05), aes(y = Tier_3, x= asinh(contribution), fill = risk_strat)) +
  geom_boxplot() +
  geom_text(data = switch1_34_p %>% filter(p_adj<0.05), aes(y = Tier_3, x =0.8, label = p_star)) +
  facet_grid(Tier_2~'', scales='free_y', space='free_y') +
  theme_minimal() +
  scale_fill_manual(values = c("#15607A", "#A63717")) +
  theme(axis.title.y = element_blank(), 
        panel.spacing = unit(1,'cm'), 
        strip.text.y = element_text(angle=0))
dev.off()


pdf(paste0('./figures/switch_1_34_',cancer_type,'.pdf'), width=10, height=12, onefile=T)
ggplot(switch1_34_p %>% filter(p_adj < 0.05), aes(y = Tier_3, x= asinh(contribution), fill = Restratification)) +
  geom_point(aes(color = Restratification), position = position_dodge(0.8), size=1) +
  geom_boxplot(alpha=0.1, coef=Inf) +
  geom_text(data = switch1_34_p %>% filter(p_adj<0.05), aes(y = Tier_3, x =0.8, label = p_star)) +
  facet_grid(Tier_2~'', scales='free_y', space='free_y') +
  theme_minimal() +
  scale_fill_manual(values = c("#15607A", "#A63717")) +
  scale_color_manual(values = c("#15607A", "#A63717")) +
  theme(axis.title.y = element_blank(), 
        panel.spacing = unit(1,'cm'), 
        strip.text.y = element_text(angle=0))
dev.off()

ggplot(switch1_34 %>% filter(!startsWith(Tier_3_tumor,'Tumor')), aes(y = Tier_3_tumor, x= contribution, fill = risk_strat)) +
  geom_boxplot() +
  facet_grid(Tier_2~1, scales='free_y', space='free_y') +
  theme_minimal()


###################
#PDL1 clusters by PJ - overrepresentation in switchers?
####################
pdl1_clusters <- paste0('X', c(253,259,254,363,278,255,279,277,479,281,260,258,261,431,491))

switchers_pdl1 <-relevant_switchers %>% 
  #filter(UICC==1, q!=2) %>% mutate(risk_strat = ifelse(q>=3, '1_34', '1_1')) %>% 
  filter(UICC>=3) %>% mutate(risk_strat = ifelse(q>=3, '34_34', '34_12')) %>%
  left_join(patient_vector_long) %>% 
  mutate(cluster = paste0('X', cluster)) %>% 
  #left_join(cell_types_fine_and_coarse) %>%
  inner_join(tumor_cell_coefficients_annotated) %>% 
  #mutate(contribution = proportion * coefficient) %>% 
  mutate(cluster_Tier_3 = paste0(Tier_3, cluster)) %>%
  #filter(Tier_3 == 'Tumor') %>% 
  mutate(ingroup = cluster %in% pdl1_clusters)

png(paste0('./figures/switchers_pdl1_PJ',cancer_type,'.png'), width=2000, height=2000, res=200)
ggplot(switchers_pdl1, aes(x = cluster, y = proportion, fill = risk_strat)) +
  geom_boxplot() +
  facet_wrap(~ingroup, scales = 'free')
dev.off()


###################
#compare coefficients to scMP
###################



scMP <- read.csv('./use_data/latent_metastatic_cells_candidates.csv')
FP <- read.csv('./use_data/cluster_fingerprint.csv') %>% mutate(cluster = paste0('X', cluster))


coefficients_and_scMP <- tumor_cell_coefficients_annotated %>% 
  dplyr::select(cluster, coefficient, Tier_2, Tier_3) %>% 
  left_join(scMP) %>% 
  left_join(FP) %>%
  filter(Tier_2 == 'Tumor')

#Manuscript
rcorr(coefficients_and_scMP$malignancy_rank, coefficients_and_scMP$coefficient)


contributions_and_scMP1_34 <- switch1_34 %>%
  #mutate() %>%
  inner_join(scMP %>% left_join(cell_types_fine_and_coarse) %>% mutate(Tier_3_tumor = paste0('Tumor', cluster))) %>% 
  left_join(FP) %>%
  filter(Tier_2 == 'Tumor') %>% 
  group_by(Tier_3_tumor, risk_strat) %>% 
  dplyr::summarize(contribution = mean(contribution), malignancy_rank = mean(malignancy_rank), fingerprint = mean(fingerprint))

contributions_and_scMP4_1 <- switch4_1 %>%
  #mutate() %>%
  inner_join(scMP %>% left_join(cell_types_fine_and_coarse) %>% mutate(Tier_3_tumor = paste0('Tumor', cluster))) %>% 
  left_join(FP) %>%
  filter(Tier_2 == 'Tumor') %>% 
  group_by(Tier_3_tumor, risk_strat) %>% 
  dplyr::summarize(contribution = mean(contribution), malignancy_rank = mean(malignancy_rank), fingerprint = mean(fingerprint))


ggplot(coefficients_and_scMP, aes(x = coefficient, y =malignancy_rank, label = cluster)) +
  geom_point() +
  geom_text_repel() +
  ylim(0,100)

ggplot(contributions_and_scMP1_34, aes(x = contribution, y =malignancy_rank, label = Tier_3_tumor)) +
  geom_point() +
  facet_wrap(~risk_strat, ncol=1) +
  geom_text_repel() +
  ylim(0,100)

ggplot(coefficients_and_scMP, aes(x = coefficient, y =fingerprint, label = cluster)) +
  geom_point() +
  geom_text_repel() 

ggplot(contributions_and_scMP1_34, aes(x = contribution, y =fingerprint, label = Tier_3_tumor)) +
  geom_point() +
  facet_wrap(~risk_strat, ncol=1) +
  geom_text_repel() 

ggplot(contributions_and_scMP4_1, aes(x = contribution, y =fingerprint, label = Tier_3_tumor)) +
  geom_point() +
  facet_wrap(~risk_strat, ncol=1) +
  geom_text_repel() 


rcorr(contributions_and_scMP1_34$contribution, contributions_and_scMP1_34$malignancy_rank, type = 'pearson')



##################################
#Validation of survival prediction
###################################
library(compareC)

cancer_type_validation <- 'SCC'

q_cutoffs_validation <- read.csv(paste0('./use_data/q_cutoffs_', cancer_type_validation, '.csv'))

validation_entity_mask <- (filtered_all_clinical_data$Entity_ac == ifelse(cancer_type_validation == 'SCC', 0,1))
validation_mask <- adj_mask & validation_entity_mask

internal_X <- all_clinical_data_with_tumor[validation_mask,]
internal_OS <- filtered_staging_and_survival$OS_m[validation_mask]
internal_event <- filtered_staging_and_survival$Event[validation_mask]

clusters = data.frame(cluster = colnames(internal_X)) %>% 
  filter(!(cluster %in% colnames(all_clinical_data))) #%>% 
#  filter(!(cluster %in% colnames(malignancy_x_staging_dummies)))
#clusters

validation_cells <- read_parquet('../scGPT_embeddings/clusterings/external_data_clustered_aligned.parquet') %>% 
  mutate(cluster= as.character(cluster)) %>%
  group_by(Pseudo, cluster) %>% 
  dplyr::summarize(N = n()) %>% 
  group_by(Pseudo) %>% 
  mutate(proportion = N/sum(N)) %>% 
  dplyr::select(-N) %>%
  arrange(cluster) %>%
  left_join(clusters,.) %>%
  pivot_wider(names_from = cluster, values_from = proportion) %>% 
  filter(!is.na(Pseudo))

patients_with_external_data <- validation_cells$Pseudo
patients_with_external_data

validation_cells[is.na(validation_cells)]<- 0
metastasis_data <- read.csv('../data/Met_Primary.csv')

validation_malignancy_ranks <- read.csv('./use_data/external_data_average_malignancy.csv') %>%  
  dplyr::select(Pseudo, malignancy_rank) %>% 
  mutate(malignancy_rank = malignancy_rank/max(malignancy_rank, na.rm=T))

validation_clinical_data_raw <- read.csv('./use_data/Matched_LC_Data_no_name.csv') %>% 
  dplyr::select(Pseudo, Entity = Histology, UICC8, UICC8_alt, date_of_diagnosis, last_contact, status, dob, sex, adj_therapy = adjuvant_cat) %>%
  mutate(Time = as.numeric(difftime(as.Date(last_contact, format="%d.%m.%Y"), as.Date(date_of_diagnosis, format= "%d.%m.%Y"), units='days'))) %>% 
  mutate(dob_date = as.Date(dob, format="%d.%m.%Y"), date_of_diagnosis_date = as.Date(date_of_diagnosis, format= "%d.%m.%Y")) %>%
  mutate(age = 1/365 * as.numeric(difftime(date_of_diagnosis_date, dob_date, units='days'))) %>%
  filter(!is.na(last_contact), last_contact != "") %>% 
  mutate(UICC8_edition3 = case_when(
    UICC8 %in% c('IA1', 'IA2', 'IA3', 'IB') ~ 'I',
    UICC8 %in% c('IIA', 'IIB') ~ 'II',
    UICC8 %in% c('IIIA', 'IIIB') ~ 'III',
    UICC8 %in% c('IVA', 'IVB') ~ 'IV',
    T ~ UICC8
  )) %>% 
  filter(!is.na(Time)) %>%
  filter(!(Pseudo %in% metastasis_data$Primary), !is.na(UICC8_alt)) %>% 
  mutate(UICC8_edition3 = factor(UICC8_edition3, levels = c('I', 'II', 'III', 'IV'))) %>%
  filter(Entity == cancer_type_validation) %>%
  mutate(adj_therapy = ifelse(is.na(adj_therapy),0,adj_therapy)) %>% 
  filter(Pseudo %in% patients_with_external_data) %>% 
  left_join(validation_malignancy_ranks) %>%
  filter(adj_therapy==0)
  #filter(!(UICC8 %in% c('III','IVA', 'IVB')))


stage_dummies_validation <- model.matrix(~ UICC8_edition3 - 1, data = validation_clinical_data_raw) %>% as.data.frame() %>%
        fill_left_zeros()
stage_dummies_validation_numeric <- validation_clinical_data_raw$UICC8_edition3 %>% as.factor() %>% as.numeric()

entity_dummies_validation <- validation_clinical_data_raw %>% dplyr::select(Entity_ac = Entity) %>% mutate(Entity_ac = ifelse(Entity_ac == 'AC',1,0))
adj_dummies_validation = validation_clinical_data_raw %>% dplyr::select(adj_therapy)

#validation_stage_x_malignancy_dummies <- (model.matrix(~ UICC8_edition3 - 1, data = validation_clinical_data_raw) %>% as.data.frame()) * validation_clinical_data_raw$malignancy_rank
#colnames(validation_stage_x_malignancy_dummies) <- paste0('malignancyx',colnames(validation_stage_x_malignancy_dummies) )
#validation_stage_x_malignancy_dummies$Pseudo <- validation_clinical_data_raw$Pseudo

clinical_data_validation <- cbind(Pseudo = validation_clinical_data_raw$Pseudo, stage_dummies_validation, entity_dummies_validation)
dim(clinical_data_validation)

external_X <- clinical_data_validation %>% left_join(validation_cells) #%>% left_join(validation_stage_x_malignancy_dummies)
external_OS <- validation_clinical_data_raw$Time
external_event <- validation_clinical_data_raw$status
UICC_baseline_extern <- validation_clinical_data_raw$UICC8_alt

all(colnames(internal_X) == colnames(external_X)[-1])
dim(internal_X)
dim(external_X)

dim(validation_cells)
dim(patient_vector)


###train model
internal_y <- Surv(internal_OS, internal_event)

best_lambda_validation1 <- cv.glmnet(internal_X %>% as.matrix(), internal_y, family = 'cox', alpha = 0.1, nfolds=20)$lambda.min
best_lambda_validation <- best_lambda_validation1
set.seed(123)
internal_model <- glmnet(internal_X %>% as.matrix(), internal_y, family = 'cox', alpha= 0.1, lambda=best_lambda_validation)

external_prediction <- predict(internal_model, external_X %>% dplyr::select(-Pseudo) %>% as.matrix(), s = best_lambda_validation, type = 'link')
external_prediction

coefs_validation <- coef(internal_model) %>%  data.frame() %>% filter(s0!=0)

res_combined_validation  <- data.frame(prediction = as.vector(external_prediction), external_OS, external_event, UICC_baseline_extern) %>% 
  cbind(q_cutoffs_validation) %>% 
  mutate(q_validation = (prediction>= (-Inf)) + (prediction>=B) + (prediction>=C) + (prediction>=D))

res_combined_validation %>% 
  group_by(UICC_baseline_extern) %>%
  dplyr::summarize(C = concordance.index(prediction, 
                                          surv.time = external_OS, 
                                          surv.event = external_event, outx = F)$c.index, 
                  p_solo = concordance.index(x = prediction, surv.time = external_OS, surv.event = external_event, outx=F, method = 'noether')$p.value,
                  C_UICC = concordance.index(UICC_baseline_extern, 
                                          surv.time = external_OS, 
                                          surv.event = external_event, outx = F)$c.index, 
                                          N=n(), 
                                          nevents = sum(external_event),
                  p_UICC_solo = concordance.index(x = UICC_baseline_extern, surv.time = external_OS, surv.event = external_event, outx=F, method = 'noether')$p.value,
                    Pvalue = compareC(external_OS, external_event, prediction, UICC_baseline_extern)$pval)


res_combined_validation %>% 
  dplyr::summarize(C = concordance.index(prediction, 
                                          surv.time = external_OS, 
                                          surv.event = external_event, outx = F, method = 'noether')$c.index, 
                  p_solo = concordance.index(x = prediction, surv.time = external_OS, surv.event = external_event, outx=F, method = 'noether')$p.value,

                  #C_RF = concordance.index(prediction_rf, 
                  #                        surv.time = external_OS, 
                  #                        surv.event = external_event, outx = F)$c.index, 
                  C_UICC = concordance.index(UICC_baseline_extern, 
                                          surv.time = external_OS, 
                                          surv.event = external_event, outx = F, method = 'noether')$c.index, 
                  p_UICC_solo = concordance.index(x = UICC_baseline_extern, surv.time = external_OS, surv.event = external_event, outx=F, method = 'noether')$p.value,
                    Pvalue = compareC(external_OS, external_event, 1-prediction, 1-UICC_baseline_extern)$pval,
                    Pvalue_compareC_alternative = cindex.comp( concordance.index(x = prediction, surv.time = external_OS, surv.event = external_event, outx=F, method = 'noether'), 
                                                                concordance.index(x=UICC_baseline_extern, surv.time = external_OS, surv.event = external_event, outx = F, method = 'noether'))$p.value,
                                          N=n(), 
                                          nevents = sum(external_event))



validation_prediction <- ggsurvplot(
  survfit(Surv(external_OS, external_event) ~ q_validation, data = res_combined_validation),
  data = res_combined_validation,
  #facet.by = c('Entity'),
  pval = F,           # Add p-value
  pval.method = TRUE,    # Show test method in the p-value plot
  pval.size = 4,         # Set the p-value font size
  conf.int = F,       # Add confidence intervals
  censor = TRUE,         # Show censored data points
  legend.title = "Risk quantile",
  title = 'snRNA',
  xlab = "Time (Days)",
  ylab = "Survival Probability",
  palette = gradient_colors,
  linewidth=0.1,
  xlim = c(0,365*5),
  ggtheme=theme_minimal())

validation_UICC <- ggsurvplot(
  survfit(Surv(external_OS, external_event) ~ UICC_baseline_extern,data = res_combined_validation),
  data = res_combined_validation,
  #facet.by = c('Entity'),
  pval = F,           # Add p-value
  pval.method = TRUE,    # Show test method in the p-value plot
  pval.size = 4,         # Set the p-value font size
  conf.int = F,       # Add confidence intervals
  censor = TRUE,         # Show censored data points
  legend.title = "Risk quantile",
  title = 'UICC',
  xlab = "Time (Days)",
  ylab = "Survival Probability",
  palette = gradient_colors,
  linewidth=0.1,
  xlim = c(0,365*5),
  ggtheme=theme_minimal())


arrange_ggsurvplots(list(validation_UICC, validation_prediction), ncol=2, print=T)

#####
combined_predictions_internal_external <- all_predictions_q %>% 
  dplyr::select(pred = full_model_pred, UICC, survival_time, survival_event) %>% 
  rbind(res_combined_validation %>% 
                        mutate(external_OS = external_OS/30) %>% 
                        dplyr::select(pred = prediction, UICC = UICC_baseline_extern, survival_time = external_OS, survival_event = external_event))

combined_predictions_internal_external %>%   
  dplyr::summarize(c_index = concordance.index(x = pred, surv.time = survival_time, surv.event = survival_event, outx=F, method = 'noether')$c.index, 
                  UICC_c_index = concordance.index(x = UICC, surv.time = survival_time, surv.event = survival_event, outx=F, method = 'noether')$c.index, 
                  p = concordance.index(x = pred, surv.time = survival_time, surv.event = survival_event, outx=F, method = 'noether')$p.value, 
                  Pvalue_compareC = compareC(survival_time, survival_event, pred, UICC)$pval,
                  Pvalue_compareC_alternative = cindex.comp(concordance.index(x = pred, surv.time = survival_time, surv.event = survival_event, outx=F, method = 'noether'), 
                                                          concordance.index(x = UICC, surv.time = survival_time, surv.event = survival_event, outx=F, method = 'noether'))$p.value,

                  N=n()) 


#########################################
#combine internal and external data ()
#########################################
both_X <- rbind(internal_X, external_X %>% dplyr::select(-Pseudo))

both_OS <- c(internal_OS, external_OS/30)
both_event <- c(internal_event, external_event)
both_Pseudo <- c(filtered_staging_and_survival$Pseudo[validation_mask], external_X$Pseudo)
both_UICC_baseline <- data.frame(patient_ids = both_Pseudo, 
                                UICC8 = c(as.numeric(filtered_staging_and_survival$UICC8_edition3[validation_mask]), UICC_baseline_extern))

set.seed(123)
both_validation_folds <- caret::createFolds(both_UICC_baseline$UICC8, k=length(both_UICC_baseline$UICC8), returnTrain = F)
both_validation_folds

#comment in for analysis
#all_clinical_predictions_both_datasets <- rbindlist(lapply(both_validation_folds, function(x) test_predictions(both_X,both_Pseudo, both_OS, both_event, x))) 

all_clinical_predictions_both_datasets %>% 
  left_join(both_UICC_baseline) %>%
  mutate(UICC8_numeric = as.numeric(as.factor(UICC8))) %>%
  #group_by(UICC8_numeric) %>%
  dplyr::summarize(c_index = concordance.index(x = pred, surv.time = surv_time_test, surv.event = surv_event_test, outx=F)$c.index, 
                  UICC_c_index = concordance.index(x = UICC8_numeric, surv.time = surv_time_test, surv.event = surv_event_test, outx=F)$c.index, 
                  p = concordance.index(x = pred, surv.time = surv_time_test, surv.event = surv_event_test, outx=F, method = 'noether')$p.value, 
                  Pvalue_compareC = compareC(surv_time_test, surv_event_test, pred, UICC8_numeric)$pval,
                  Pvalue_compareC_alternative = cindex.comp(concordance.index(x = pred, surv.time = surv_time_test, surv.event = surv_event_test, outx=F, method = 'noether'), 
                                                          concordance.index(x = UICC8_numeric, surv.time = surv_time_test, surv.event = surv_event_test, outx=F, method = 'noether'))$p.value,
                  N=n()) 


coefs_both <- test_coef_all(both_X,both_OS, both_event, both_Pseudo) %>% 
  mutate(index = rownames(.)) %>% 
  pivot_longer(!index, names_to = 'cluster', values_to = 'coefficient') %>% 
  filter(coefficient!=0) %>% 
  left_join(cell_types_fine_and_coarse)

ggplot(coefs_both, aes(x = asinh(coefficient), y = Tier_3, color = coefficient, size = Freq)) +
  geom_point(alpha = 1) +  # Use geom_point() without jitter
  facet_grid(Tier_2 ~ '', scales = "free_y", space = "free_y") +  # Allow free y-space for better fitting
  geom_vline(xintercept = 0, color = "#0A3659", linetype = "solid") +  # Add vertical line at x = 0
  scale_color_gradientn(colors = c("#15607A", "#C7CDD1", "#A63717"), limits = c(-6, 6), oob = scales::squish) +  # Custom 3-color palette
  scale_y_discrete() +  # Now that Tier_3_rev.y is a factor, no need to specify limits again
  #scale_x_continuous(limits = c(-2, 1)) +  # Adjust the x-axis range
  theme_minimal() +
  labs(x = "Asinh(Coefficient)",
       y = "Tier 3 Annotation",
       color = "Coefficient") +  # Color legend label
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.spacing.x = unit(3.5, "lines"))  # Adjust the spacing between panels

write.csv(coefs_both, './use_data/important_clusters_SCC_combined_datasets.csv')






min_max_prediction_values_both_datasets <- all_clinical_predictions_both_datasets %>% 
  dplyr::summarize(minv = quantile(pred, 0.05), maxv = quantile(pred, 0.95), diffv = maxv-minv) 
  #dplyr::summarize(minv = min(full_model_pred), maxv = max(full_model_pred), diffv = maxv-minv) 

equal_sized_q_both_datasets <- data.frame(data.frame(cutoffs = (seq(nquantiles)-1)) * (min_max_prediction_values_both_datasets$diffv/nquantiles)+ min_max_prediction_values_both_datasets$minv) %>% 
  t() 
colnames(equal_sized_q_both_datasets) = c('A', 'B', 'C', 'D')
write.csv(equal_sized_q_both_datasets, paste0('./use_data/q_cutoffs_both_datasets_', cancer_type, '.csv'))

all_predictions_q_both_datasets <- all_clinical_predictions_both_datasets %>% 
    left_join(both_UICC_baseline) %>%
  mutate(UICC8_numeric = as.numeric(as.factor(UICC8))) %>%
  #mutate(q = case_when(full_model_pred <=q_fav$q_fav_cutoff ~ 1, 
  #                    full_model_pred >q_fav$q_fav_cutoff & full_model_pred <= q_unfav$q_unfav_cutoff ~ 2, 
  #                    full_model_pred >q_unfav$q_unfav_cutoff ~3)) %>% 
          cbind(equal_sized_q_both_datasets) %>% 
          mutate(q = (pred>= (-Inf)) + (pred>=B) + (pred>=C) + (pred>=D)) %>% 
          mutate(surv_time_test_capped = pmin(surv_time_test, 120))




surv_fit_tumor_cell_both_datasets <- survfit(Surv(surv_time_test_capped, surv_event_test) ~ q, data = all_predictions_q_both_datasets)
surv_fit_TNM_both_datasets <- survfit(Surv(surv_time_test_capped, surv_event_test) ~ UICC8, data = all_predictions_q_both_datasets)


new_both_datasets <- ggsurvplot(
  surv_fit_tumor_cell_both_datasets,
  data = all_predictions_q_both_datasets,
  pval = T,
  pval.method = TRUE,
  pval.size = 4,
  conf.int = FALSE,
  censor = TRUE,
  legend.title = "Risk quantile",
  title = "Full model",
  xlab = "Time (Months)",
  ylab = "Survival Probability",
  palette = c("#11B5AE", "#4046C9", "#F68512", "#DE3C82"),
  linewidth = 0.1,
  ggtheme = theme_minimal(),
  lwd = 0.5,
  risk.table = TRUE,       # Add the risk table below the plot
  risk.table.title = "Number at Risk", # Title for the risk table
  risk.table.height = 0.25 # Adjust height of the risk table
)

#new_both_datasets$table <- new$table + 
#  scale_y_discrete(labels = rev(c("1", "2", "3", "4"))) # Customize y-axis labels

new_both_datasets

tnm_both_datasets <- ggsurvplot(
  surv_fit_TNM_both_datasets,
  data = all_predictions_q_both_datasets,
  pval = T,
  pval.method = TRUE,
  pval.size = 4,
  conf.int = FALSE,
  censor = TRUE,
  legend.title = "UICC stage",
  title = "UICC8",
  xlab = "Time (Months)",
  ylab = "Survival Probability",
  palette = c("#11B5AE", "#4046C9", "#F68512", "#DE3C82"),
  linewidth = 0.1,
  ggtheme = theme_minimal(),
  lwd = 0.5,
  risk.table = TRUE,       # Add the risk table below the plot
  risk.table.title = "Number at Risk", # Title for the risk table
  risk.table.height = 0.25 # Adjust height of the risk table
)

tnm_both_datasets

arrange_ggsurvplots(list(tnm_both_datasets, new_both_datasets), ncol=2, print=T)


all_predictions_q_both_datasets %>% dim()




