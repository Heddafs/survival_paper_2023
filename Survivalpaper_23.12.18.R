## Survival paper

## ---- Setup----------------------------------------------------------
setwd("~/Documents/Forskerlinje/Bioinformatics_mammatumor/Survival_paper")
library(tidyverse)
library(DESeq2)
library(splitstackshape)
library(survival)
library(glmnet)
library(c060)
library(caret)
source("helper_coxph_all.R") 
# helper_time_dependant.R for predictions using cv.glmnet
#use helper_coxph_all.R for the helper-script for predictions from  coxph models instead of cv.glmnet models. 
library(VIM)
library(faux)
load("data/tumor_samples_with_subtype.Rdata")
load("data/counts_tumors.Rdata")

## ------Formatting annotation file--------------------------------------
# A bit of formatting: making sure variables that we will use later are the correct data type 
# Also imputing missing values with KNN and dropping dogs with no Survival days information. 
# Here, values are imputed for 10 dogs (With 1-2 missing values each dog). 11 dogs were excluded due to NAs for Survival_days
clin_vars = c("Age", "ER_status", "Lymphatic_invasion")

clinical_tbl <- tumor_samples_with_subtype %>% 
  mutate(subtype_MPAM50 = as.factor(subtype_MPAM50)) %>% 
  kNN(variable = clin_vars, k=5) %>% 
  drop_na("Survival_days") %>% 
  mutate(Age = as.vector(scale(Age))) %>%
  select(c(Sample_ID, all_of(clin_vars), Status_survival, Survival_days)) %>%  
  # Only keeping variables we will use
  as_tibble()

# Making treatment contrast groups (there is probably a smarter way to do this 
# that does not involve using columnnames...)
clinical_onehot <- clinical_tbl %>%
  mutate(ER_positive = ifelse(ER_status == "N", 0, 1)) %>%
  mutate(Lymphatic_invasion_present = ifelse(Lymphatic_invasion == "Absent", 0, 1)) %>%
  select(-ER_status, -Lymphatic_invasion)

#skimr::skim(clinical_onehot) 
# All variables are numeric data type (Except for sample_ID)
# Mean survival days is 327 --> makes sense to set cut-off at 1 year if we do time dependant discrimination metric.

## -------Formatting RNAseq count file ---------------------------------
# Making sure counts_tumors (Integer matrix with number of counts per gene in each sample) columns are the same as clinical_tbl Sample ID column - this is essential since DESeq doesn't make any guesses###
counts_tumors <- counts_tumors %>% 
  select(clinical_tbl$Sample_ID)
all(colnames(counts_tumors) %in% clinical_tbl$Sample_ID) 
all(colnames(counts_tumors) == clinical_tbl$Sample_ID) #true --> same order of dogIDs

#Preprocessing: normalizing, variance filtering, removing highly correlated genes, 
# log2 transforming, scaling and univariate cox filtering. 
gene_data <- preprocess_genes(x = counts_tumors, y=clinical_tbl, variance_cutoff = 1000, corr_cutoff = 0.95, fdr_cutoff = 0.1)

###### Checking correlations ######## 

clin_vars <- c("Age", "ER_positive", "Lymphatic_invasion_present")
corr_matrix <- matrix(nrow = length(clin_vars), ncol = ncol(gene_data), dimnames = list(clin_vars, colnames(gene_data)))
for (i in colnames(gene_data)) {
  for (j in clin_vars) {
  corr <- stats::cor(gene_data[,i], clinical_onehot[,j])
  corr_matrix[j,i] <- corr
}
}

corrplot(corr_matrix)

### Maybe change to cancor::CCA #####
install.packages("candisc")

correlations <- stats::cancor(gene_data, clinical_onehot[,c(2,5,6)])
correlations$xcoef

prcomp()

# Looks like the correlation is no more than 0.46 between each clinical variable and gene

## -------Model fitting ---------------------------------
### dataset splitting, alpha tuning (both and gene model), model fitting, 
# model prediction and C- index calculation (with a self made function model_evaluation() 
# (see helper.R for code)
Results_model_eval <- model_evaluation(x = clinical_onehot, 
                            y = gene_data, 
                            nrep = 20,                   # number of iterations
                            split.percentage = 0.8, 
                            method =  '1se', 
                            nfolds = 3,                 # number of folds
                            type.measure = "C", 
                            lambda = "lambda.min", 
                            pf = c(rep(0, length(clin_vars)), # Clinical variables are not penalized
                                rep(1, ncol(gene_data)))) 

save(Results_model_eval, file = "results_modelevaluation_23.12.16_tester.Rda") # saving results table 

# Plotting histogram of the alpha value chosen for both models for each of the iterations
alphas_plot <- Results_model_eval %>%   
  filter(Model_id != "Cox.clinical") %>% 
  group_by(Iteration, Model_id) %>% 
  summarise(alpha = mean(Best_alpha)) %>% 
  ggplot(aes(x=alpha))+
  geom_histogram(aes(y = after_stat(count / (sum(count)/2))), bins = 11)+
  facet_wrap("Model_id")+
  ggtitle("Tuned alpha value used for cv.glmnet()")
ggsave(alphas_plot, file = "alphas_plot_23.12.16_tester.jpeg", height = 6, width = 8)

#plotting obtained C-indexes for the three models.
performance_plot <- ggplot(Results_model_eval, aes(x = Model_id, y = value))+ 
  geom_boxplot()+
  facet_wrap(~Metric_id)+
  theme_bw()+
  ylim(0.2,0.9)
ggsave(performance_plot,file="performance_plot_23.12.16_tester.jpeg", height = 6, width = 8)


