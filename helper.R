#################### Count matrix preprocessing function #####################################
#' @description
#' Function that inputs a count matrix and
#' - Normalizes counts using DESEq2
#' - Filters top variable genes
#' - Removes highly correlated genes
#' - log2 transformes
#' - scales
#' - Filters genes with univariate cox filtering
#' 
#' 
#' @param x count matrix with genes as rownames and sample IDs as column names
#' @param y tibble with sample IDs in a column called Sample_ID, clinical variables and survival information (in colums calles Survival_days and Status_survival). Make sure SampleIDs are in the same order in object x and y
#' @param variance_cutoff number of top variable genes to keep
#' @param corr_cutoff remove one gene in each gene pair with a variance above cutoff
#' @param fdr_cutoff Keeps genes with an fdr below cutoff
#' 
#' @returns Returns a gene expression matrix. Sample IDs are rownames and gene_names are column names

preprocess_genes <- function(x, y, variance_cutoff = 1000, corr_cutoff = 0.95, fdr_cutoff = 0.1){
  # Changing gene names to avoid problems later on. 
  rownames(x) <- make.names(rownames(x))
  
  # Normalizing count data with DESeq2. 
  # Used clinical_tbl(which has the original variables) and not the clinical_onehot. 
  dds_tumor <- DESeqDataSetFromMatrix(countData = x, 
                                      colData = y, 
                                      design = ~ 1) 
  dds_tumor <- estimateSizeFactors(dds_tumor)
  normalized_counts <- counts(dds_tumor, normalized=TRUE)
  normalized_counts <- normalized_counts[,clinical_tbl$Sample_ID] #making sure 
  # normalized counts are in the correct order. 
  
  ## Feature selection with M3C::featurefilter() = Keeps top 1000 most variable genes
  ngenes = nrow(normalized_counts)
  ngenes_tokeep = variance_cutoff
  per = (ngenes_tokeep/ngenes)*100
  filtered = M3C::featurefilter(normalized_counts, percentile = per, method = 'MAD')
  filtered_normcounts = filtered$filtered_data
  genes_topVar <- t(filtered_normcounts)
  
  # Removing highly correlated variables: removing one gene in each gene pair with over 95% correlation
  genes_topVar_cor <- cor(genes_topVar)
  idx <- findCorrelation(genes_topVar_cor, cutoff=corr_cutoff) 
  idx <- sort(idx)
  length(idx) # number of genes removed = 72
  gene_data <- genes_topVar[,-c(idx)]
  
  
  # log2 transformation
  gene_data <- log(gene_data+1)
  
  # Standardizing gene expression
  gene_data <- scale(gene_data, center = T, scale = T)
  ncol(gene_data) # 928 genes are included up until this point
  
  # Univariate cox filtering
  gene_df <- cbind(gene_data, "Status_survival" = clinical_tbl$Status_survival, 
                   "Survival_days" = clinical_tbl$Survival_days) %>% as.data.frame()
  p.values <- vapply(colnames(gene_data), function(gene){
    print(gene)
    model <- coxph(as.formula(sprintf("Surv(Survival_days, Status_survival) ~ %s", gene)), data = gene_df)
    summary(model)$coefficients[,"Pr(>|z|)"]
  }, numeric(1))
  
  p.values.adj <- p.adjust(p.values, method = "fdr")
  
  genes_univariate_keep <- names(p.values.adj[p.values.adj<fdr_cutoff])
  length(genes_univariate_keep) # 47 genes are kept
  
  gene_data[,genes_univariate_keep]
  
}





#################### tune_alpha function #####################################
#' @description
#' A function that can be used for tuning alphas for a cv_glmnet() model.
#' The function outputs a list with a tibble of the C-indexes or deviances for
#' the tested alphas ($tune_alpha_tbl) and the best alpha after tuning (e.g. the
#' one with the highest C-index ($best_alpha)).

#' @param alphas sequence of alphas to tune
#' @param x train data matrix
#' @param y Survival object (time and event) for train data matrix
#' @param foldid vector of length nrow(x) identyfying the fold group for each sample
#' @param pf Penealty factor, which penalty factor to use for each variable (column in x). 0 = do not penalize, 1 = can be penalized.
#' @param type.measure discrimination metric used for tuning lambda. "C" for C- index, "deviance"
#' @param method How to selected the best alpha. `best` selects the alpha with the
#' best score (depends on the `type.measure`, can be the alpha with the highest
#' mean C-index, or the one with the lowest mean Partial Likelihood Deviance).
#' `1se` selects the higher alpha that results in a mean score which is within
#' 1 standard-error from the best alpha (thus favoring closer to Lasso-based models).
#'
#' @returns a list with a table of the tested alphas and measures (C-indexes or
#' deviance), a graph of alphas vs the chosen metric and the final selected alpha
#' based on the method

tune_alpha <- function(alphas = seq(0, 1, 0.1), x, y, foldid, type.measure = "C", pf=rep(1, ncol(x)), method = "best") {
  atbl = future.apply::future_lapply(1:length(alphas), function(i) {
    alpha = alphas[i]
    coxnet <- cv.glmnet(x = x,
      y = y, family = "cox", foldid = foldid, standardize = FALSE, penalty.factor = pf,
      alpha = alpha, type.measure = type.measure, maxit = 10000)

    lambda_min_index <- coxnet$index["min",]
    nfolds = length(unique(foldid))

    tibble(
      alpha = alpha,
      low  = coxnet$cvlo[lambda_min_index],
      mean = coxnet$cvm[lambda_min_index],
      up   = coxnet$cvup[lambda_min_index],
      sem  = coxnet$cvsd[lambda_min_index]/sqrt(nfolds), # standard error of the mean
      model = list(coxnet)
    )
  }, future.seed = TRUE) %>% dplyr::bind_rows()

  # select best alpha based on mean CV score
  alpha_seq <- atbl %>%
    arrange(desc(mean)) %>%
    pull(alpha)
  alpha_cv = ifelse(type.measure == 'C', alpha_seq[1], alpha_seq[length(alphas)])

  if (method == "best") {
    # alpha with best mean CV score (C-index or deviance)
    sel_alpha = alpha_cv
  } else {
    # alpha within 1se of the best mean CV score
    cvm = atbl %>% filter(alpha == alpha_cv) %>% pull(mean)
    sem = atbl %>% filter(alpha == alpha_cv) %>% pull(sem)

    if (type.measure == 'C') { # C-index, higher is better
      sel_alpha = atbl %>%
        filter(mean <= cvm, mean >= cvm - sem) %>% # C-index - SEM
        slice_tail() %>% # choose higher alpha
        pull(alpha)
    } else { # deviance, lower is better
      sel_alpha = atbl %>%
        filter(mean >= cvm, mean <= cvm + sem) %>% # PLDeviance + SEM
        slice_tail() %>% # choose higher alpha
        pull(alpha)
    }
  }

  # get trained model for the selected alpha
  best_model = atbl %>%
    filter(alpha == sel_alpha) %>%
    pull(model) %>%
    `[[`(1) # unlist it

  plot <- ggplot(atbl, aes(x=alpha)) +
    geom_ribbon(aes(ymin = low, ymax = up), fill="grey80") +
    geom_ribbon(aes(ymin = mean - sem, ymax = mean + sem), fill="grey30") +
    geom_line(aes(y=mean)) +
    # add selected alpha to the plot
    geom_vline(xintercept = sel_alpha, col = 'red', linetype = 'dashed') +
    ylab(ifelse(type.measure == "C", "C-index", "Partial Likelihood Deviance"))

  list(tune_alpha_tbl = atbl, best_alpha = sel_alpha,
       best_model = best_model, plot = plot)
}

#' @param x one hot encoded clinical variables, including outcome columns (Status_survival and Survival_days) and a column specifying sample id (Sample_ID)
#' @param y matrix with normalized, log2 transformed and scaled gene expression values for samples(rows). Genenames in columns.
#' @param split.percentage percentage of samples that should be included in train dataset. Default is 0.67
#' @param nrep number of times to repeat model fitting and prediction
#' @param nfolds number of folds to use for crossvalidation. Default is 5
#' @param type.measure Discrimination metric type to use in tune alpha function
#' @param lambda which lambda will be used for prediction (`lambda.min` or `lambda.1se`)
#' @param pf penalty factor to use for the shared model. 

model_evaluation <- function(x, y, split.percentage = 0.67, nrep, nfolds = 3,
  type.measure = "C", method = "best", lambda = "lambda.min", pf =  c(rep(0, (ncol(train_x_clin))), rep(1, ncol(train_x_both)-ncol(train_x_clin)))) {

  results_list = lapply(1:nrep, function(i) {
    print(i)

    ## ---------Splitting data into train and test------------------------
    train <- stratified(x,'Status_survival', split.percentage) %>%
      column_to_rownames("Sample_ID")
    test <- x %>% dplyr::filter(!Sample_ID %in% rownames(train)) %>%
      column_to_rownames("Sample_ID")

    ## Making test/train datasets for **predictors** (clinical variables, genes and a shared dataset)
    train_x_clin <- train %>%
      select(-c(Survival_days, Status_survival))
    test_x_clin <- test %>%
      select(-c(Survival_days, Status_survival))

    train_x_genes <- y[rownames(train),]
    test_x_genes <- y[rownames(test),]

    train_x_both <- train_x_clin %>%
      merge(train_x_genes, by=0) %>%
      column_to_rownames("Row.names") %>%
      as.matrix()

    test_x_both <- test_x_clin %>%
      merge(test_x_genes, by=0) %>%
      column_to_rownames("Row.names") %>%
      as.matrix()

    # Making the survival object for train and test datasets = **Outcome variables**
    train_y <- Surv(train$Survival_days,train$Status_survival)
    test_y <- Surv(test$Survival_days,test$Status_survival)

    ## ----------Fitting Cox.clinical ------------------------------------
    clin_formula <- as.formula(paste0("Surv(Survival_days, Status_survival) ~ ", paste0(colnames(train_x_clin), collapse = " + ")))
    Cox.clinical <- coxph(clin_formula, data = train)

    ## ----------Fitting CoxNet.genes and optimizing alpha-------------------
    all(rownames(train_x_genes) == train$Sample_ID)
    utils::capture.output({
      fold <- c060::balancedFolds(class.column.factor = train$Status_survival, cross.outer = nfolds) #making random folds (stratified by censoring) for crossvalidation. same folds are used for fitting CoxNet.both model
    })

    tune_alpha_res_genes <- tune_alpha(x = train_x_genes, y = train_y, foldid = fold, type.measure = type.measure, method = method) #tuning alpha

    CoxNet.genes <- tune_alpha_res_genes$best_model # Get best CoxNet.gene model
    
    vec_coef_genes <- coef(CoxNet.genes, s= CoxNet.genes[[lambda]])[,1]
    vec_coef_genes <- vec_coef_genes[vec_coef_genes!= 0]
    gene_formula <- as.formula(paste0("Surv(Survival_days, Status_survival) ~ ", paste0(names(vec_coef_genes), collapse = " + ")))
    train_both <- train %>%
      merge(train_x_genes, by=0) %>%
      column_to_rownames("Row.names")
    train_x_genes <- train_x_genes %>% as.data.frame()
    test_x_genes <- test_x_genes %>%  as.data.frame()
    
    Cox.genes <- coxph(gene_formula, data = train_both)
    
    ## -------------Fitting CoxNet.both (clinical and genes)----------------
    pf <- pf
    tune_alpha_res_clingen <- tune_alpha(x = train_x_both, y = train_y, foldid = fold, pf = pf, type.measure = type.measure, method = method) #tuning alpha

    CoxNet.both <- tune_alpha_res_clingen$best_model # get best model
    vec_coef <- coef(CoxNet.both, s= CoxNet.both[[lambda]])[,1]
    vec_coef <- vec_coef[vec_coef!= 0]
    
    train_x_both <- train_x_both %>% as.data.frame()
    test_x_both <- test_x_both %>%  as.data.frame()
    train_both <- train %>%
      merge(train_x_genes, by=0) %>%
      column_to_rownames("Row.names")
   
    both_formula <- as.formula(paste0("Surv(Survival_days, Status_survival) ~ ", paste0(names(vec_coef), collapse = " + ")))
    Cox.both <- coxph(both_formula, data = train_both)
    

    ## ------------- Predicting and evaluating performance------------------------
    ## Cox.clinical
    predict_clinical <- predict(Cox.clinical, newdata = test_x_clin, type = "lp")
    Cindex_Harrell_model1 <- glmnet::Cindex(pred = predict_clinical, y = test_y)
    Cindex_Uno_model1 <- survAUC::UnoC(train_y, test_y, predict_clinical)

    ## CoxNet.genes
    predict_genes <- predict(Cox.genes, newdata = test_x_genes, type = "lp")
    Cindex_Harrell_model2 <- glmnet::Cindex(pred = predict_genes, y = test_y)
    Cindex_Uno_model2 <- survAUC::UnoC(Surv.rsp = train_y, Surv.rsp.new = test_y, lpnew = predict_genes)
    
    ## CoxNet.both
    predict_both <- predict(Cox.both, newdata = test_x_both, type = "lp")
    Cindex_Harrell_model3 <- glmnet::Cindex(pred = predict_both, y = test_y)
    Cindex_Uno_model3 <- survAUC::UnoC(Surv.rsp = train_y, Surv.rsp.new = test_y, lpnew = predict_both)
    #browser()
    
    ### Time dependant discrimination metric, at 3 timepoints
    Predicting_tbl <- tibble(Clinical_lp = predict_clinical, 
                             Both_lp = predict_both, 
                             Genes_lp = predict_genes, 
                             time=test$Survival_days, 
                             status=test$Status_survival)
    
    Time_dependant_clinical <- Score(object = list(Predicting_tbl$Clinical_lp),
                                     formula = Surv(time, status) ~ 1,
                                     data = Predicting_tbl,
                                     times = c(183, 365, 730),
                                     metrics = 'auc')
    ROC_AUC_clinical_6m <- Time_dependant_clinical$AUC$score$AUC[1]
    ROC_AUC_clinical_1y <- Time_dependant_clinical$AUC$score$AUC[2]
    ROC_AUC_clinical_2y <- Time_dependant_clinical$AUC$score$AUC[3]
    
    
    Time_dependant_genes <- Score(object = list(Predicting_tbl$Genes_lp),
                                  formula = Surv(time, status) ~ 1,
                                  data = Predicting_tbl,
                                  times = c(183, 365, 730),
                                  metrics = 'auc')
    ROC_AUC_genes_6m <- Time_dependant_genes$AUC$score$AUC[1]  
    ROC_AUC_genes_1y <- Time_dependant_genes$AUC$score$AUC[2]  
    ROC_AUC_genes_2y <- Time_dependant_genes$AUC$score$AUC[3]  
    
    
    Time_dependant_both <- Score(object = list(Predicting_tbl$Both_lp),
                                 formula = Surv(time, status) ~ 1,
                                 data = Predicting_tbl,
                                 times = c(183, 365, 730),
                                 metrics = 'auc')
    ROC_AUC_both_6m <- Time_dependant_both$AUC$score$AUC[1]
    ROC_AUC_both_1y <- Time_dependant_both$AUC$score$AUC[2]
    ROC_AUC_both_2y <- Time_dependant_both$AUC$score$AUC[3]
    
    
    ## Storing results
    tibble("Iteration" = i,
           "Cox.clinical" = c(Cindex_Harrell_model1, Cindex_Uno_model1, 
                              ROC_AUC_clinical_6m, ROC_AUC_clinical_1y, ROC_AUC_clinical_2y),
           "Cox.genes" = c(Cindex_Harrell_model2, Cindex_Uno_model2, 
                           ROC_AUC_genes_6m, ROC_AUC_genes_1y, ROC_AUC_genes_2y),
           "Cox.both" = c(Cindex_Harrell_model3, Cindex_Uno_model3, 
                          ROC_AUC_both_6m, ROC_AUC_both_1y, ROC_AUC_both_2y),
           "Metric_id" = c("Harrel's C-index", "Uno's C-index", "Time_dependant_AUC_6m", "Time_dependant_AUC_1y", "Time_dependant_AUC_2y"), 
           "nfeatures_both_model" = length(vec_coef), 
           "features_both_model" = list(names(vec_coef)), 
           "nfeatures_gene_model" = length(vec_coef_genes), 
           "features_gene_model" = list(names(vec_coef_genes))
      ) %>%
      pivot_longer(cols = c(Cox.clinical, Cox.genes, Cox.both),
        names_to = "Model_id") %>%
      mutate("Best_alpha" = ifelse(Model_id == "Cox.clinical", NA, ifelse(Model_id == "Cox.genes", tune_alpha_res_genes$best_alpha, tune_alpha_res_clingen$best_alpha)))
  })

  dplyr::bind_rows(results_list)
}
