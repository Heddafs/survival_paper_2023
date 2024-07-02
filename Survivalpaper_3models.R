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
source("helper.R") 
library(VIM)
library(faux)
library(corrplot)
library(riskRegression)
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
preprocess_res <- preprocess_genes(x = counts_tumors, y=clinical_tbl, variance_cutoff = 1000, corr_cutoff = 0.95, fdr_cutoff = 0.1)
results_univariate_genes <- as.data.frame(preprocess_res$univariate_result) %>% rownames_to_column(var="Gene names") # results univariate tests genes - fdr adj p-value. 
writexl::write_xlsx(results_univariate_genes, "plots/results_univariate_genes.xlsx")

str(preprocess_res$genes_uncorrelated) # 999 genes were tested in univariate --> only one gene removed due to high correlation.

gene_data <- preprocess_res$univariate_significant_genes


###### Checking correlations ######## 

clin_vars_corr <- c("Age" ,"ER_positive", "Lymphatic_invasion_present")
variables_both <- clinical_onehot %>% 
  column_to_rownames("Sample_ID") %>% 
  select(all_of(clin_vars_corr)) %>% 
  merge(gene_data, by=0) %>%
  column_to_rownames("Row.names")

corr_matrix <- matrix(nrow = ncol(variables_both), ncol = ncol(variables_both), dimnames = list(colnames(variables_both), colnames(variables_both)))
for (i in colnames(variables_both)) {
  for (j in colnames(variables_both)) {
  corr <- stats::cor(variables_both[,i], variables_both[,j])
  corr_matrix[j,i] <- corr
}
}

corr_plot <- corrplot(corr_matrix, order = "hclust", tl.col = "red4", tl.cex = 0.75, tl.srt = 45, type = "lower", diag = F)
# Looks like the correlation is no more than 0.46 between each clinical variable and gene
colnames(corr_matrix)[1:3] <- ""
rownames(corr_matrix)[1:3] <- ""
corr_plot <- corrplot(corr_matrix, order = "hclust", tl.col = "black", tl.cex = 0.75, tl.srt = 45, type = "lower", diag = F, add = T)


## -------Model fitting ---------------------------------
### dataset splitting, alpha tuning (both and gene model), model fitting, 
# model prediction and C- index calculation (with a self made function model_evaluation() 
# (see helper.R for code)
Results_model_eval <- model_evaluation(x = clinical_onehot, 
                            y = gene_data, 
                            nrep = 100,                   # number of iterations
                            split.percentage = 0.8, 
                            method =  '1se', 
                            nfolds = 3,                 # number of folds
                            type.measure = "C", 
                            lambda = "lambda.min", 
                            pf = c(rep(0, length(clin_vars)), # Clinical variables are not penalized
                                rep(1, ncol(gene_data)))) 

save(Results_model_eval, file = "results_modelevaluation_2024_06_25.Rda") # saving results table 

# Plotting histogram of the alpha value chosen for both models for each of the iterations
alphas_plot <- Results_model_eval %>%   
  filter(Model_id != "Cox.clinical") %>% 
  group_by(Iteration, Model_id) %>% 
  summarise(alpha = mean(Best_alpha)) %>% 
  ggplot(aes(x=alpha))+
  geom_histogram(aes(y = after_stat(count / (sum(count)/2))), bins = 11)+
  facet_wrap("Model_id")+
  ggtitle("Tuned alpha value used for cv.glmnet()")+
  theme_bw()
ggsave(alphas_plot, file = "alphas_plot_2024_06_25.jpeg", height = 6, width = 8)

# Median and range of number of terms included in the combined model
Results_model_eval %>% filter(Metric_id == "Uno's C-index", Model_id == "Cox.both") %>% pull(nfeatures_gene_model) %>% median
Results_model_eval %>% filter(Metric_id == "Uno's C-index", Model_id == "Cox.both") %>% pull(nfeatures_both_model) %>% range

# Median and range of number of terms included in gene only model
Results_model_eval %>% filter(Metric_id == "Uno's C-index", Model_id == "Cox.genes") %>% pull(nfeatures_gene_model) %>% median
Results_model_eval %>% filter(Metric_id == "Uno's C-index", Model_id == "Cox.genes") %>% pull(nfeatures_gene_model) %>% range


#plotting obtained C-indexes for the three models.
my_comparisons <- list( c("Cox.both", "Cox.clinical"), c("Cox.both", "Cox.genes"), c("Cox.clinical", "Cox.genes") )

# Overall C-indexes:
C_index_plot <- Results_model_eval %>% 
  filter(Metric_id == "Uno's C-index") %>% 
  mutate(across(Model_id, 
                ~factor(., levels=c("Cox.clinical", "Cox.genes", "Cox.both")))) %>% 
  ggplot(., aes(x = Model_id, y = value))+ 
  geom_boxplot(fill=c("red4", "dodgerblue4", "green4"))+
  ylim(0,1.05)+
  ylab("Uno's C-index")+
  xlab("")+
  theme_bw(base_size = 15)+ 
  scale_x_discrete(labels = c("1) Clinical","2) Genes" , "3) Combined"))+
  ggpubr::stat_compare_means(comparisons = my_comparisons, label.y = c(0.89, 0.95, 0.99))

#calculate means for Uno's C-indexes
Results_model_eval %>% filter(Model_id=="Cox.clinical", Metric_id == "Uno's C-index") %>% pull(value) %>% mean()
Results_model_eval %>% filter(Model_id=="Cox.genes", Metric_id == "Uno's C-index") %>% pull(value) %>% mean()
Results_model_eval %>% filter(Model_id=="Cox.both", Metric_id == "Uno's C-index") %>% pull(value) %>% mean()


# Time dependant C-index = ROC-AUC(t)
timepoints <- c("Time_dependant_AUC_6m" = "6 months",
                "Time_dependant_AUC_1y" = "1 year",
                "Time_dependant_AUC_2y" = "2 years")


ROC_AUC_plot <- Results_model_eval %>% 
  filter(Metric_id %in% c("Time_dependant_AUC_6m", "Time_dependant_AUC_1y", "Time_dependant_AUC_2y")) %>% 
  mutate(across(Metric_id, 
                ~factor(., levels=c("Time_dependant_AUC_6m", "Time_dependant_AUC_1y", "Time_dependant_AUC_2y")))) %>%
  mutate(across(Model_id, 
                ~factor(., levels=c("Cox.clinical", "Cox.genes", "Cox.both")))) %>% 
  ggplot(., aes(x = Model_id, y = value))+ 
  geom_boxplot(fill=rep(c("red4", "dodgerblue4", "green4"), 3))+
  facet_wrap(~Metric_id, labeller = as_labeller(timepoints))+
  theme_bw(base_size = 13)+ 
  ylim(0,1.15)+
  ylab("ROC-AUC(t)")+
  xlab("")+
  scale_x_discrete(labels = c("1) Clinical","2) Genes" , "3) Combined"))+
  ggpubr::stat_compare_means(comparisons = my_comparisons, label.y = c(0.99, 1.05, 1.1))

Results_model_eval %>% filter(Model_id=="Cox.clinical", Metric_id == "Time_dependant_AUC_6m") %>% pull(value) %>% mean()


ggsave(C_index_plot,file="plots/C_index_plot_2024_06_25.jpeg", height = 5, width = 4)
ggsave(ROC_AUC_plot,file="plots/ROC_AUC_plot_2024_06_25.jpeg", height = 6, width = 9.5)




