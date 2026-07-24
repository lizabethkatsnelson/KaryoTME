library(dplyr)
library(purrr)
library(stringr)




parse_cytoband_info <- function(bands) {
  # strip the "Arm" prefix (if present) so arm-level names like "Arm1p"
  # are parsed the same way as cytoband-level names like "1p36.33"
  clean_bands <- str_remove(bands, "^Arm")
  
  tibble(Cytoband = bands) %>%
    mutate(
      Chromosome = str_extract(clean_bands, "^(\\d+|X|Y)"),
      Arm = str_extract(clean_bands, "[pq]"),
      ChrArm = paste0(Chromosome, Arm)
    )
}

fit_one_band <- function(y, cnv, band_name, model, covariates = NULL) {
  
  if(is.null(model)){model = "ALL"}
  
  # base data
  dat <- data.frame(
    y = y,
    CNV = cnv
  )
  
  # if covariates are provided, bind them in
  if (!is.null(covariates)) {
    covariates <- as.data.frame(covariates)
    dat <- cbind(dat, covariates)
  }
  
  # drop missing values
  dat <- dat[complete.cases(dat), ]
  
  # CNV status
  dat$cnv_status <- cut(
    dat$CNV,
    breaks = c(-Inf, -0.2, 0.2, Inf),
    labels = c("Loss", "Neu", "Gain")
  )
  if(model != "ALL"){
    dat <- dat[dat$cnv_status == model | dat$cnv_status == "Neu", ]
  }
  
  if (nrow(dat) < 10) {
    return(tibble(
      Cytoband = band_name,
      Coefficients = NA,
      T_val = NA,
      P_val = NA,
      cnv_type = "con",
      Event = model))
  }
  
  # automatically build the formula
  cov_names <- setdiff(colnames(dat), c("y", "CNV", "cnv_status"))
  
  if (length(cov_names) == 0) {
    form <- as.formula("y ~ CNV")
  } else {
    form <- as.formula(
      paste("y ~ CNV +", paste(cov_names, collapse = " + "))
    )
  }
  
  fit <- tryCatch(
    lm(form, data = dat),
    error = function(e) NULL
  )
  
  if (is.null(fit)) {
    return(tibble(
      Cytoband = band_name,
      Coefficients = NA,
      T_val = NA,
      P_val = NA,
      cnv_type = "con",
      Event = model
    ))
  }
  
  sm <- summary(fit)
  
  if (!"CNV" %in% rownames(sm$coefficients)) {
    return(tibble(
      Cytoband = band_name,
      Coefficients = NA,
      T_val = NA,
      P_val = NA,
      cnv_type = "con",
      Event = model
    ))
  }
  
  res = tibble(
    Cytoband = band_name,
    Coefficients = sm$coefficients["CNV", "Estimate"],
    T_val = sm$coefficients["CNV", "t value"],
    P_val = sm$coefficients["CNV", "Pr(>|t|)"],
    cnv_type = "con",
    Event = model
  )
  
  fit = fit
  return(list("fit" = fit,
              "res" = res))
}




eval_one_band <- function(fit, band, y, cnv, covariates, test_idx) {
  
  if (is.null(fit)) {
    return(tibble(
      Cytoband = band,
      R2 = NA,
      MSE = NA,
      Cor = NA,
      n_test = 0
    ))
  }
  
  # build test data
  dat_test <- data.frame(
    y = y,
    CNV = cnv
  )
  
  if (!is.null(covariates)) {
    dat_test <- cbind(dat_test, covariates)
  }
  
  dat_test <- dat_test[complete.cases(dat_test), ]
  
  if (nrow(dat_test) < 10) {
    return(tibble(
      Cytoband = band,
      R2 = NA,
      MSE = NA,
      Cor = NA,
      n_test = nrow(dat_test)
    ))
  }
  
  # predict
  y_pred <- tryCatch(
    predict(fit, newdata = dat_test),
    error = function(e) rep(NA, nrow(dat_test))
  )
  
  y_true <- dat_test$y
  
  # metrics
  mse <- mean((y_true - y_pred)^2)
  r2 <- cor(y_true, y_pred)^2
  cor_val <- cor(y_true, y_pred)
  
  tibble(
    Cytoband = band,
    R2 = r2,
    MSE = mse,
    Cor = cor_val,
    n_test = nrow(dat_test)
  )
}




linear_function = function(thread_cutoff = NULL,  ###using top and bottom #thread_cutoff precentage patients as
                           model = NULL, ###sep, NULL = ALL patients
                           covariates = NULL, ###cov df
                           y, cnv
){
  
  patient_name = rownames(cnv)
  thread_cutoff = thread_cutoff
  avail_idx = y > quantile(y,1-thread_cutoff) |  y <= quantile(y,thread_cutoff)
  if(is.null(model)){
    list_res = map(cytobands, function(band) {
      fit_one_band(y = y[avail_idx],
                   cnv = cnv[avail_idx, band],
                   band_name = band,
                   covariates = covariates[avail_idx,,drop = F],
                   model = NULL)})
    summary_df = lapply(list_res,function(x){ x[["res"]]}) %>% bind_rows %>%
      left_join(cyto_info, by = "Cytoband")
    list_fit = lapply(list_res,function(x){x[["fit"]]})
    res = summary_df
  }
  
  # ======================
  # regression
  # ======================
  if(!is.null(model)){
    ## gain
    list_res = map(cytobands, function(band) {
      fit_one_band(y = y[avail_idx],
                   cnv = cnv[avail_idx, band],
                   band_name = band,
                   covariates = covariates[avail_idx,,drop = F],
                   model = "Gain")})
    summary_df_gain = lapply(list_res,function(x){ x[["res"]]}) %>% bind_rows %>%
      left_join(cyto_info, by = "Cytoband")
    list_fit = lapply(list_res,function(x){x[["fit"]]})
    
    ## loss
    list_res = map(cytobands, function(band) {
      fit_one_band(y = y[avail_idx],
                   cnv = cnv[avail_idx, band],
                   band_name = band,
                   covariates = covariates[avail_idx,,drop = F],
                   model = "Loss")})
    summary_df_loss = lapply(list_res,function(x){ x[["res"]]}) %>% bind_rows %>%
      left_join(cyto_info, by = "Cytoband")
    list_fit = lapply(list_res,function(x){x[["fit"]]})
    
    res =  rbind(summary_df_gain, summary_df_loss)
  }
  
  return(res)
}




















linear_cv_function = function(kfold = NULL, ##split patients into k folds
                              seed = 12345, ##seeds used for splitting
                              nfold_iter = NULL, ### number of folds in each iteration
                              thread_cutoff = NULL,  ###using top and bottom #thread_cutoff precentage patients as
                              model = NULL, ###sep, NULL = ALL patients
                              covariates = NULL, ###cov df
                              y, cnv
){
  
  patient_name = rownames(cnv)
  n = nrow(cnv)
  kfold =kfold
  set.seed(seed)
  fold_id <- sample(rep(1:kfold, length.out = n))
  names(fold_id) <- patient_name
  combn_list = combn(kfold,nfold_iter,simplify = FALSE)
  niter = length(combn_list)
  
  fold_results <- list()
  thread_cutoff = thread_cutoff
  
  for (iter in 1:niter) {
    cat("iteration:", iter, "\n")
    fold_iter =  combn_list[[iter]]
    train_idx = grep(paste0(fold_iter, collapse="|"),fold_id )
    test_idx = grep(paste0(fold_iter, collapse="|"),fold_id,invert = T )
    train_avail_idx = y[train_idx] > quantile(y[train_idx],1-thread_cutoff) |  y[train_idx] <= quantile(y[train_idx],thread_cutoff)
    test_avail_idx = y[test_idx] > quantile(y[test_idx],1-thread_cutoff) |  y[test_idx] <= quantile(y[test_idx],thread_cutoff)
    
    if(is.null(model)){
      list_res = map(cytobands, function(band) {
        fit_one_band(y = y[train_idx][train_avail_idx],
                     cnv = cnv[train_idx, band][train_avail_idx],
                     band_name = band,
                     covariates = covariates[train_idx,,drop = F][train_avail_idx,,drop = F],
                     model = NULL)})
      summary_df = lapply(list_res,function(x){ x[["res"]]}) %>% bind_rows %>%
        left_join(cyto_info, by = "Cytoband") %>%
        mutate(Fold = iter)
      list_fit = lapply(list_res,function(x){x[["fit"]]})
      
      
      test_df <- purrr::map_dfr(cytobands, function(band) {
        eval_one_band(
          fit = list_fit[[band]],
          band = band,
          y = y[test_idx][test_avail_idx],
          cnv = cnv[test_idx,band][test_avail_idx],
          covariates = covariates[test_idx,,drop=F][test_avail_idx,,drop=F])})
      test_df = left_join(test_df,cyto_info, by = "Cytoband") %>%
        mutate(Fold = iter)
      
      res = list("summary_df" = summary_df,
                 "test_df" = test_df)
      
    }
    
    # ======================
    # regression
    # ======================
    if(!is.null(model)){
      ## gain
      list_res = map(cytobands, function(band) {
        fit_one_band(y = y[train_idx][train_avail_idx],
                     cnv = cnv[train_idx, band][train_avail_idx],
                     band_name = band,
                     covariates = covariates[train_idx,,drop = F][train_avail_idx,,drop = F],
                     model = "Gain")})
      summary_df_gain = lapply(list_res,function(x){ x[["res"]]}) %>% bind_rows %>%
        left_join(cyto_info, by = "Cytoband") %>%
        mutate(Fold = iter)
      list_fit = lapply(list_res,function(x){x[["fit"]]})
      
      
      test_df_gain <- purrr::map_dfr(cytobands, function(band) {
        eval_one_band(
          fit = list_fit[[band]],
          band = band,
          y = y[test_idx][test_avail_idx ],
          cnv = cnv[test_idx,band][test_avail_idx],
          covariates = covariates[test_idx,,drop=F][test_avail_idx,,drop=F])})
      test_df_gain = left_join(test_df_gain,cyto_info, by = "Cytoband") %>%
        mutate(Fold = iter)
      test_df_gain$Event = "Gain"
      
      
      ## loss
      list_res = map(cytobands, function(band) {
        fit_one_band(y = y[train_idx][train_avail_idx],
                     cnv = cnv[train_idx, band][train_avail_idx],
                     band_name = band,
                     covariates = covariates[train_idx,,drop = F][train_avail_idx,,drop = F],
                     model = "Loss")})
      summary_df_loss = lapply(list_res,function(x){ x[["res"]]}) %>% bind_rows %>%
        left_join(cyto_info, by = "Cytoband") %>%
        mutate(Fold = iter)
      list_fit = lapply(list_res,function(x){x[["fit"]]})
      
      
      test_df_loss <- purrr::map_dfr(cytobands, function(band) {
        eval_one_band(
          fit = list_fit[[band]],
          band = band,
          y = y[test_idx][test_avail_idx],
          cnv = cnv[test_idx,band][test_avail_idx],
          covariates = covariates[test_idx,,drop=F][test_avail_idx,,drop=F] )})
      test_df_loss$Event = "Loss"
      test_df_loss = left_join(test_df_loss,cyto_info, by = "Cytoband") %>%
        mutate(Fold = iter)
      
      
      res = list( "summary_df" = rbind(summary_df_gain, summary_df_loss),
                  "test_df" = rbind(test_df_gain, test_df_loss) )
      
      #res_gain = map_dfr(cytobands, function(band) {
      #   fit_one_band(y = y[train_idx][avail_idx],
      #                cnv = cnv[train_idx, band][avail_idx],
      #                band_name = band,
      #                covariates = covariates[train_idx,,drop = F][avail_idx,,drop = F],
      #                model = "Gain") }) %>%
      #  left_join(cyto_info, by = "Cytoband") %>%
      #  mutate(Fold = iter)
      
      #res_loss = map_dfr(cytobands, function(band) {
      #   fit_one_band(y = y[train_idx][avail_idx],
      #                cnv = cnv[train_idx, band][avail_idx],
      #                band_name = band,
      #                covariates = covariates[train_idx,,drop = F][avail_idx,,drop = F],
      #                model = "Loss") }) %>%
      #  left_join(cyto_info, by = "Cytoband") %>%
      #  mutate(Fold = iter)
      
      #res = rbind(res_gain, res_loss)
      
    }
    
    fold_results[[iter]] <- res
  }
  return(fold_results)
}















logisticfit_one_band <- function(y, cnv, band_name, model,is_thresh = 0.2, covariates = NULL) {
  
  if(is.null(model)){model = "ALL"}
  is_thresh = 0.3
  # base data
  dat <- data.frame(
    y = y,
    CNV = cnv
  )
  
  dat$y = ifelse(dat$y>= unname(quantile(dat$y, probs = 1-is_thresh)), 1,ifelse(dat$y <= unname(quantile(dat$y, probs = is_thresh)), 0, NA) )
  
  # if covariates are provided, bind them in
  if (!is.null(covariates)) {
    covariates <- as.data.frame(covariates)
    dat <- cbind(dat, covariates)
  }
  
  # drop missing values
  dat <- dat[complete.cases(dat), ]
  
  # CNV status
  dat$cnv_status <- cut(
    dat$CNV,
    breaks = c(-Inf, -0.2, 0.2, Inf),
    labels = c("Loss", "Neu", "Gain")
  )
  if(model != "ALL"){
    dat <- dat[dat$cnv_status == model | dat$cnv_status == "Neu", ]
    
  }
  
  if (nrow(dat) < 10) {
    return(tibble(
      Cytoband = band_name,
      Coefficients = NA,
      Z_val = NA,
      P_val = NA,
      cnv_type = "con",
      Event = model))
  }
  
  # automatically build the formula
  cov_names <- setdiff(colnames(dat), c("y", "CNV", "cnv_status"))
  
  if (length(cov_names) == 0) {
    form <- as.formula("y ~ CNV")
  } else {
    form <- as.formula(
      paste("y ~ CNV +", paste(cov_names, collapse = " + "))
    )
  }
  
  
  
  fit <- tryCatch(
    glm(form, data = dat, family = binomial),
    error = function(e) NULL
  )
  
  if (is.null(fit)) {
    return(tibble(
      Cytoband = band_name,
      Coefficients = NA,
      Z_val = NA,
      P_val = NA,
      cnv_type = "con",
      Event = model
    ))
  }
  
  sm <- summary(fit)
  
  if (!"CNV" %in% rownames(sm$coefficients)) {
    return(tibble(
      Cytoband = band_name,
      Coefficients = NA,
      Z_val = NA,
      P_val = NA,
      cnv_type = "con",
      Event = model
    ))
  }
  
  res = tibble(
    Cytoband = band_name,
    Coefficients = sm$coefficients["CNV", "Estimate"],
    Z_val = sm$coefficients["CNV", "z value"],
    P_val = sm$coefficients["CNV", "Pr(>|z|)"],
    cnv_type = "con",
    Event = model
  )
  
  fit = fit
  return(list("fit" = fit,
              "res" = res))
}




logisticeval_one_band <- function(fit, band, y, cnv, covariates, is_thresh = 0.2,test_idx) {
  
  if (is.null(fit)) {
    return(tibble(
      Cytoband = band,
      R2 = NA,
      MSE = NA,
      Cor = NA,
      n_test = 0
    ))
  }
  
  # build test data
  dat_test <- data.frame(
    y = y,
    CNV = cnv
  )
  dat_test$y = ifelse(dat_test$y>= unname(quantile(dat_test$y, probs = 1-is_thresh)), 1,ifelse(dat_test$y <= unname(quantile(dat_test$y, probs = is_thresh)), 0, NA) )
  
  if (!is.null(covariates)) {
    dat_test <- cbind(dat_test, covariates)
  }
  
  dat_test <- dat_test[complete.cases(dat_test), ]
  
  if (nrow(dat_test) < 10) {
    return(tibble(
      Cytoband = band,
      R2 = NA,
      MSE = NA,
      Cor = NA,
      n_test = nrow(dat_test)
    ))
  }
  
  # predict
  y_pred <- tryCatch(
    predict(fit, newdata = dat_test),
    error = function(e) rep(NA, nrow(dat_test))
  )
  
  y_true <- dat_test$y
  
  # metrics
  mse <- mean((y_true - y_pred)^2)
  r2 <- cor(y_true, y_pred)^2
  cor_val <- cor(y_true, y_pred)
  
  tibble(
    Cytoband = band,
    R2 = r2,
    MSE = mse,
    Cor = cor_val,
    n_test = nrow(dat_test)
  )
}






logistic_cv_function = function(kfold = NULL, ##split patients into k folds
                                seed = 12345, ##seeds used for splitting
                                nfold_iter = NULL, ### number of folds in each iteration
                                thread_cutoff = NULL,  ###using top and bottom #thread_cutoff precentage patients as
                                model = NULL, ###sep, NULL = ALL patients
                                covariates = NULL, ###cov df
                                y, cnv
){
  
  patient_name = rownames(cnv)
  n = nrow(cnv)
  kfold =kfold
  set.seed(seed)
  fold_id <- sample(rep(1:kfold, length.out = n))
  names(fold_id) <- patient_name
  combn_list = combn(kfold,nfold_iter,simplify = FALSE)
  niter = length(combn_list)
  
  fold_results <- list()
  thread_cutoff = thread_cutoff
  
  for (iter in 1:niter) {
    cat("iteration:", iter, "\n")
    fold_iter =  combn_list[[iter]]
    train_idx = grep(paste0(fold_iter, collapse="|"),fold_id )
    test_idx = grep(paste0(fold_iter, collapse="|"),fold_id,invert = T )
    train_avail_idx = y[train_idx] > quantile(y[train_idx],1-thread_cutoff) |  y[train_idx] <= quantile(y[train_idx],thread_cutoff)
    test_avail_idx = y[test_idx] > quantile(y[test_idx],1-thread_cutoff) |  y[test_idx] <= quantile(y[test_idx],thread_cutoff)
    
    if(!is.null(model)){
      ## gain
      list_res = map(cytobands, function(band) {
        logisticfit_one_band(y = y[train_idx][train_avail_idx],
                             cnv = cnv[train_idx, band][train_avail_idx],
                             band_name = band,
                             covariates = covariates[train_idx,,drop = F][train_avail_idx,,drop = F],
                             model = "Gain")})
      summary_df_gain = lapply(list_res,function(x){ x[["res"]]}) %>% bind_rows %>%
        left_join(cyto_info, by = "Cytoband") %>%
        mutate(Fold = iter)
      list_fit = lapply(list_res,function(x){x[["fit"]]})
      
      
      test_df_gain <- purrr::map_dfr(cytobands, function(band) {
        logisticeval_one_band(
          fit = list_fit[[band]],
          band = band,
          y = y[test_idx][test_avail_idx ],
          cnv = cnv[test_idx,band][test_avail_idx],
          covariates = covariates[test_idx,,drop=F][test_avail_idx,,drop=F])})
      test_df_gain = left_join(test_df_gain,cyto_info, by = "Cytoband") %>%
        mutate(Fold = iter)
      test_df_gain$Event = "Gain"
      
      
      ## loss
      list_res = map(cytobands, function(band) {
        logisticfit_one_band(y = y[train_idx][train_avail_idx],
                             cnv = cnv[train_idx, band][train_avail_idx],
                             band_name = band,
                             covariates = covariates[train_idx,,drop = F][train_avail_idx,,drop = F],
                             model = "Loss")})
      summary_df_loss = lapply(list_res,function(x){ x[["res"]]}) %>% bind_rows %>%
        left_join(cyto_info, by = "Cytoband") %>%
        mutate(Fold = iter)
      list_fit = lapply(list_res,function(x){x[["fit"]]})
      
      
      test_df_loss <- purrr::map_dfr(cytobands, function(band) {
        logisticeval_one_band(
          fit = list_fit[[band]],
          band = band,
          y = y[test_idx][test_avail_idx],
          cnv = cnv[test_idx,band][test_avail_idx],
          covariates = covariates[test_idx,,drop=F][test_avail_idx,,drop=F] )})
      test_df_loss$Event = "Loss"
      test_df_loss = left_join(test_df_loss,cyto_info, by = "Cytoband") %>%
        mutate(Fold = iter)
      
      
      res = list( "summary_df" = rbind(summary_df_gain, summary_df_loss),
                  "test_df" = rbind(test_df_gain, test_df_loss) )
      
      #res_gain = map_dfr(cytobands, function(band) {
      #   fit_one_band(y = y[train_idx][avail_idx],
      #                cnv = cnv[train_idx, band][avail_idx],
      #                band_name = band,
      #                covariates = covariates[train_idx,,drop = F][avail_idx,,drop = F],
      #                model = "Gain") }) %>%
      #  left_join(cyto_info, by = "Cytoband") %>%
      #  mutate(Fold = iter)
      
      #res_loss = map_dfr(cytobands, function(band) {
      #   fit_one_band(y = y[train_idx][avail_idx],
      #                cnv = cnv[train_idx, band][avail_idx],
      #                band_name = band,
      #                covariates = covariates[train_idx,,drop = F][avail_idx,,drop = F],
      #                model = "Loss") }) %>%
      #  left_join(cyto_info, by = "Cytoband") %>%
      #  mutate(Fold = iter)
      
      #res = rbind(res_gain, res_loss)
      
    }
    
    fold_results[[iter]] <- res
  }
  return(fold_results)
}





















#########sep model#####

### Functino ot plot IS CNV regression model output
plot_cnv_IS_regression <- function(df, cnvlevel=c("Gene", "Cytoband", "Arm"), scorename="IS",
                                   stat=NULL, cnvtype=c("bi", "con"),
                                   signif=2,sec_signif =1, plottitle="", savepath=NULL) {
  
  ### process table
  df <- df[!is.na(df$Chromosome),] # remove rows with missing chrom info
  df$Chromosome <- factor(gsub("Chr", "", df$Chromosome), levels = c(1:22, "X", "Y"))
  df$ChrArm <- factor(gsub("Chr", "", df$ChrArm), levels = na.omit(unique(gsub("Chr", "", df$ChrArm))))
  
  if (cnvlevel == "Gene") {
    df <- df[order(df$Chromosome, df$Start),] # order by chr and start site
  }
  if (cnvlevel == "Cytoband") {
    # df <- df[order(df$Cytoband),] # order by cytoband
    df$Cytoband <- factor(df$Cytoband, levels = unique(df$Cytoband))
  }
  if (cnvlevel == "Arm") {
    df <- df[order(df$ChrArm),] # order by arm
  }
  
  # plotting with p val and binary cnv
  if (stat == "P_val" & cnvtype == "bi") {
    df$signedlogP <- ifelse(df$Coefficients > 0, -log10(df$P_val), -log10(df$P_val)*-1 ) # if coef pos, logP is pos, if neg logP is neg
    df$signedlogP[is.na(df$signedlogP)] <- 0 # change NA to 0
    stat = "signedlogP" # change which var to use for plot
  }
  
  # plotting with pval and continuous cnv
  if (stat == "P_val" & cnvtype == "con") {
    df$signedlogP <- -log10(df$P_val)
    df$signedlogP <- ifelse(df$Event == "Gain" & df$Coefficients < 0, df$signedlogP *-1, # if gain and (-) weight, logP (-)
                            ifelse(df$Event == "Loss" & df$Coefficients > 0, df$signedlogP *-1, df$signedlogP)) # if loss and (+) weight, logP (-)
    df$signedlogP[is.na(df$signedlogP)] <- 0 # change NA to 0
    stat = "signedlogP" # change which var to use for plot
  }
  
  # plotting with z val (log) or t val (lin) and binary cnv - don't need to change anything, statistic already correct direction
  
  # plotting with z or t(lin) val and continuous cnv
  if (stat %in% c("Z_val", "T_val") & cnvtype == "con") { ### because stat already signed, just need to change direction for losses
    df[,stat] <- ifelse(df$Event == "Loss", df[,stat]*-1, df[,stat]) # if loss, change dir of stat (neg becomes pos, pos becomes neg)
  }
  
  ### get x axis cutoffs for genome
  chr.cut = na.omit(df %>% group_by(Ref = Chromosome) %>% summarise(vals = last(get(cnvlevel)))) # end of chr
  
  arm.cut <- df %>% group_by(Ref = ChrArm) %>% summarise(vals = last(get(cnvlevel))) # end of arms
  arm.cut$vals <- ifelse(grepl("q", arm.cut$Ref), NA, as.character(arm.cut$vals)) # remove q arm (don't add line to end of q)
  arm.cut$Ref <- gsub("p|q", "", arm.cut$Ref)
  arm.cut <- arm.cut[!duplicated(arm.cut$Ref),]
  
  chr.mid <- df %>% group_by(Ref = Chromosome) %>% filter(row_number()==ceiling(n()/2)) %>% dplyr::select(Ref, vals=cnvlevel) # middle of each chr
  
  y = sym(stat) ## turn string into symbol
  x = sym(cnvlevel)
  
  ymax = ceiling(max(df[,stat]))
  ymin = floor(min(df[,stat]))
  
  if( length(grep("R2|MSE|Cor",stat)) == 1){ }
  
  ymax = max(abs(ymin), ymax)
  #if (ymax < 4) { ymax = 4 }
  if( length(grep("R2|MSE|Cor",stat)) == 1){ymax = ymax * 1.1
  } else {if (ymax < 4) { ymax = 4 } }
  
  
  if (cnvlevel %in% c("Gene", "Cytoband")) {
    p <- ggplot(df) +
      geom_area(aes(x = !!x, y = !!y, group = Event, fill = Event), position = "identity") +
      scale_fill_manual(values = alpha(c("#ad0224","#020dad"),0.6)) +
      labs(title = plottitle) +
      theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
            title = element_text(size=10), axis.text.y = element_text(size=9)) +
      scale_x_discrete(limits = unique(df[,cnvlevel])) +
      ylim(-ymax, ymax) +
      geom_vline(xintercept = arm.cut$vals, linetype="dotted", color = "#8F8F8F", linewidth=0.15) +
      geom_vline(xintercept = chr.cut$vals, linetype="longdash", color = "#8F8F8F", linewidth=0.3) +
      geom_hline(yintercept = c(-signif, signif), linetype="dashed", color = "black", linewidth=0.3) +
      geom_hline(yintercept = c(-sec_signif, sec_signif ), linetype="dashed", color = "black", linewidth=0.15) +
      geom_hline(yintercept = 0, color = "black", linewidth=0.15) +
      geom_text(data = chr.mid, mapping = aes(x = vals, y = -ymax, label = Ref, hjust = 0, vjust = 0), angle=90, size=2)
  }
  
  if (cnvlevel == "Arm") {
    p <- ggplot() +
      geom_bar(data = df[df$Event == "Loss",], aes(x = !!x, y = !!y, fill = Event), stat = "identity", width = 1) +
      geom_bar(data = df[df$Event == "Gain",], aes(x = !!x, y = !!y, fill = Event), stat = "identity", width = 1) +
      scale_fill_manual(values = alpha(c("#ad0224","#020dad"),0.6)) +
      labs(title = plottitle) + ylim(-ymax, ymax) + theme_bw()
  }
  
  ### save
  if (is.null(savepath)) { return(p) } else { # if no path given to save pdf, return the plot as output for function
    pdf(savepath, width = 10, height = 5)
    print(p)
    junk <- dev.off()
  }
}




######full model#####

plot_cnv_cv_signedlogp <- function(df,
                                   cnvlevel = c("Cytoband", "Arm"),
                                   stat_col = "signedlogP",
                                   signif = 2,
                                   sec_signif = 1,
                                   plottitle = "",
                                   cnv_df = NULL  ####used to determine the sign
) {
  
  cnvlevel <- match.arg(cnvlevel)
  
  # drop missing values
 # df <- df[!is.na(df$Chromosome), ]
  df <- df[!is.na(df[[stat_col]]), ]
  
  # chromosome order
  df$Chromosome <- factor(gsub("Chr", "", df$Chromosome),
                          levels = c(as.character(1:22), "X", "Y"))
  
  df$ChrArm <- factor(gsub("Chr", "", df$ChrArm),
                      levels = na.omit(unique(gsub("Chr", "", df$ChrArm))))
  
  # Cytoband order
  if (cnvlevel == "Cytoband") {
    df$Cytoband <- factor(df$Cytoband, levels = unique(df$Cytoband))
  }
  
  # Arm order
  if (cnvlevel == "Arm") {
    df <- df[order(df$ChrArm), ]
  }
  
  # -----------------------------
  # chr / arm boundary positions
  # -----------------------------
  chr.cut <- na.omit(
    df %>%
      group_by(Ref = Chromosome) %>%
      summarise(vals = last(get(cnvlevel)), .groups = "drop")
  )
  
  arm.cut <- df %>%
    group_by(Ref = ChrArm) %>%
    summarise(vals = last(get(cnvlevel)), .groups = "drop")
  
  # keep only the p-arm boundary lines (consistent with the original code)
  arm.cut$vals <- ifelse(grepl("q", arm.cut$Ref), NA, as.character(arm.cut$vals))
  arm.cut$Ref <- gsub("p|q", "", arm.cut$Ref)
  arm.cut <- arm.cut[!duplicated(arm.cut$Ref), ]
  
  # chr label position
  chr.mid <- df %>%
    group_by(Ref = Chromosome) %>%
    filter(row_number() == ceiling(n() / 2)) %>%
    dplyr::select(Ref, vals = all_of(cnvlevel))
  
  # -----------------------------
  # y-axis range
  # -----------------------------
  ymax <- ceiling(max(df[[stat_col]], na.rm = TRUE))
  ymin <- floor(min(df[[stat_col]], na.rm = TRUE))
  ymax <- max(abs(ymin), ymax)
  
  if( length(grep("R2|MSE|Cor",stat_col)) == 1){ymax = ymax * 1.1
  } else {if (ymax < 4) { ymax = 4 } }
  
  
  #if (ymax < 4) ymax <- 4
  
  x <- rlang::sym(cnvlevel)
  y <- rlang::sym(stat_col)
  
  df$cnv_status = "not defined"
  color_set = c()
  if(!is.null(cnv_df)){
    cnv_status = sign(apply(cnv_df,2,mean))
    df[[stat_col]] = cnv_status * df[[stat_col]]
    df$cnv_status = ifelse(cnv_status > 0, "Gain","Loss")
    df_padding = df
    df_padding[[stat_col]] = 0
    df_padding$cnv_status = ifelse(df_padding$cnv_status =="Gain","Loss","Gain")
    df = rbind(df,df_padding)
    color_set = c("#ad0224","#020dad")
  }else{
    color_set = "#4C78A8"
  }
  
  # -----------------------------
  # Cytoband plot
  # -----------------------------
  if (cnvlevel == "Cytoband") {
    
    p <- ggplot(df) +
      geom_area(aes(x = !!x, y = !!y, group = cnv_status,fill = cnv_status)) +
      scale_fill_manual(values = alpha(color_set,0.6)) +
      labs(title = plottitle, y = stat_col) +
      
      theme_classic() +
      theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 10),
        axis.text.y = element_text(size = 9)
      ) +
      
      #      scale_x_discrete(limits = unique(df[, cnvlevel])) +
      ylim(-ymax, ymax) +
      
      # arm boundary lines
      geom_vline(xintercept = arm.cut$vals,
                 linetype = "dotted",
                 color = "#8F8F8F",
                 linewidth = 0.15) +
      
      # chromosome boundary lines
      geom_vline(xintercept = chr.cut$vals,
                 linetype = "longdash",
                 color = "#8F8F8F",
                 linewidth = 0.3) +
      
      # significance lines
      geom_hline(yintercept = c(-signif, signif),
                 linetype = "dashed",
                 color = "black",
                 linewidth = 0.3) +
      
      # ±1 reference lines
      geom_hline(yintercept = c(-sec_signif, sec_signif),
                 linetype = "dashed",
                 color = "black",
                 linewidth = 0.15) +
      
      # zero line
      geom_hline(yintercept = 0,
                 color = "black",
                 linewidth = 0.15) +
      
      # chromosome labels
      geom_text(
        data = chr.mid,
        mapping = aes(x = vals, y = -ymax, label = Ref, hjust = 0, vjust = 0),
        angle = 90,
        size = 2,
        inherit.aes = FALSE
      )
  }
  
  # -----------------------------
  # Arm-level plot
  # -----------------------------
  if (cnvlevel == "Arm") {
    
    p <- ggplot(df, aes(x = !!x, y = !!y)) +
      geom_col(fill = "#4C78A8", width = 1) +
      labs(title = plottitle, y = "Signed -log10(P)") +
      ylim(-ymax, ymax) +
      theme_bw() +
      geom_hline(yintercept = c(-signif, signif),
                 linetype = "dashed",
                 color = "black",
                 linewidth = 0.3) +
      geom_hline(yintercept = c(-sec_signif , sec_signif ),
                 linetype = "dashed",
                 color = "black",
                 linewidth = 0.15) +
      geom_hline(yintercept = 0,
                 color = "black",
                 linewidth = 0.15)
  }
  
  return(p)
}




lr_plot_function=function(df,
                          cnv_df = cnv_df,
                          signif = 2,
                          sec_signif = 1,
                          stat_col =NULL,
                          model = c("sep","full"),
                          plottitle = plottitle, # "wgdcov linear IS ~ Cytoband_CNV(con0.2) + AS"
                          savepath =savepath ####paste0("../output/wgdcov_linearregression_", fold, "_is_Cytobandcnv_continuous0.2_plot_logP.pdf")
){
  
  if(model == "sep"){
    plot_cnv_IS_regression(df = df,
                           cnvlevel = "Cytoband",
                           scorename = "IS",
                           stat = stat_col,
                           cnvtype = "con",
                           signif = signif,
                           sec_signif = sec_signif,
                           plottitle = plottitle,
                           savepath = savepath )
    p = NULL
  }
  
  if(model == "full"){
    p = plot_cnv_cv_signedlogp(
      df = df,
      cnvlevel = "Cytoband",
      stat_col = stat_col,
      signif = signif,
      cnv_df = cnv_df,
      sec_signif = sec_signif,
      plottitle = plottitle)
  }
  
  return(p)
}

















variable_compute_function = function(model = c("sep", "full"),
                                     list_df = fold_results,
                                     stat_col = c("signedlogP")) {
  
  model <- match.arg(model)
  fold_results_df <- bind_rows(list_df)
  
  ## mean + var expressions
  mean_exprs <- setNames(
    lapply(stat_col, function(col) {
      rlang::expr(mean(.data[[!!col]], na.rm = TRUE))
    }),
    paste0("mean_", stat_col)
  )
  
  var_exprs <- setNames(
    lapply(stat_col, function(col) {
      rlang::expr(var(.data[[!!col]], na.rm = TRUE))
    }),
    paste0("var_", stat_col)
  )
  
  summary_exprs <- c(mean_exprs, var_exprs)
  
  if (model == "sep") {
    
    summary_df <- fold_results_df %>%
      group_by(Cytoband, Chromosome, Arm, ChrArm, Event) %>%
      summarise(!!!summary_exprs, .groups = "drop") %>%
      mutate(
        Cytoband = factor(Cytoband, levels = cyto_levels),
        Chromosome = factor(Chromosome, levels = c(as.character(1:22), "X", "Y"))
      )
    
    ## kept for compatibility with existing plotting logic (e.g. when signedlogP is included)
    
    fold_1 <- fold_results_df[fold_results_df$Fold == 1, ]
    summary_df <- as.data.frame(summary_df)
    
    summary_df_gain <- summary_df[summary_df$Event == "Gain", ]
    summary_df_loss <- summary_df[summary_df$Event == "Loss", ]
    
    rownames(summary_df_gain) <- summary_df_gain$Cytoband
    rownames(summary_df_loss) <- summary_df_loss$Cytoband
    
    summary_df_gain <- summary_df_gain[fold_1$Cytoband, ]
    summary_df_loss <- summary_df_loss[fold_1$Cytoband, ]
    
    summary_df <- rbind(summary_df_gain, summary_df_loss)
    
    return(summary_df)
  }
  
  if (model == "full") {
    
    summary_df <- fold_results_df %>%
      group_by(Cytoband, Chromosome, Arm, ChrArm) %>%
      summarise(!!!summary_exprs, .groups = "drop") %>%
      mutate(
        Cytoband = factor(Cytoband, levels = cyto_levels),
        Chromosome = factor(Chromosome, levels = c(as.character(1:22), "X", "Y"))
      )
    
    
    fold_1 <- fold_results_df[fold_results_df$Fold == 1, ]
    summary_df <- as.data.frame(summary_df)
    rownames(summary_df) <- summary_df$Cytoband
    summary_df <- summary_df[fold_1$Cytoband, ]
    
    return(summary_df)
  }
}






###

compute_geneset_score <- function(counts, geneset) {
  
  library(DESeq2)
  
  dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = data.frame(row.names = colnames(counts)),
    design = ~ 1
  )
  
  dds <- estimateSizeFactors(dds)
  norm_counts <- counts(dds, normalized = TRUE)
  
  log_counts <- log2(norm_counts + 1)
  
  z_mat <- t(scale(t(log_counts)))
  
  geneset <- intersect(geneset, rownames(z_mat))
  
  score <- colMeans(z_mat[geneset, , drop = FALSE])
  
  return(score)
}