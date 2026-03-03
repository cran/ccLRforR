ps4.logistic <- function(gene = c("BRCA1", "BRCA2", "PALB2", "CHEK2", "ATM", "TP53", "custom"),
                      genotypes,
                      geno_notation = c("n", "n/n"),
                      phenotype,
                      custom_rules = NULL,
                      outdir = NULL,
                      output = "PS4",
                      stratifyby = NULL,
                      agefilter = c(0, 80),
                      exportcsv = FALSE,
                      progress = FALSE) {
  
  if (is.null(outdir)) {
    data.folder.location <- tempdir()
  } else {
    data.folder.location <- outdir
  }
  
  
  if (!(gene %in% c("BRCA1", "BRCA2", "PALB2", "CHEK2", "ATM", "TP53", "custom"))) {
    stop("The gene type (stated via the gene argument) should be one of: BRCA1 or BRCA2 or PALB2 or CHEK2 or ATM or TP53 or custom.")
  }
  
  
  if (!(geno_notation %in% c("n", "n/n"))) {
    stop("The genotype notation format should be one of: n or n/n.")
  }
  
  
  if (!is.null(stratifyby) && !(stratifyby %in% c("country", "ethnicity", "study"))) {
    stop("The stratification must be by country, ethnicity, study, or NULL (no stratification). Revise your stratifyby choice.")
  }
  
  
  if (
    is.null(agefilter) ||
    !is.numeric(agefilter) ||
    length(agefilter) != 2 ||
    any(is.na(agefilter)) ||
    agefilter[1] > agefilter[2]
  ) {
    stop(
      "agefilter must be a numeric vector of length 2 specifying ",
      "the minimum and maximum ages (e.g., c(0, 80) or c(21, 80))."
    )
  }
  if (agefilter[1] < 0) {
    stop("agefilter minimum age cannot be negative.")
  }
  

  
  gen_expected_column <- c("sample_ids")
  if (!(colnames(genotypes)[1]) == gen_expected_column) {
    stop("The genotype data do not match the required first column name which is: sample_ids.")
  }
  
  
  phenotype <- as.data.frame(phenotype)
  phen_ncols <- ncol(phenotype)
  if (identical(stratifyby, "country")) {
    if (phen_ncols != 5) stop("The number of columns in the phenotype data frame is not 5 as required.")
    phen_expected_columns <- c("sample_ids", "status", "ageInt", "AgeDiagIndex", "StudyCountry")
    if (!all(sort(colnames(phenotype)) == sort(phen_expected_columns))) {
      stop("The phenotype data frame does not match the required column names which are: \nsample_ids, status, ageInt, AgeDiagIndex, StudyCountry.")
    }
    phenotype <- phenotype[ , c("sample_ids", "status", "ageInt", "AgeDiagIndex", "StudyCountry")]
    phenotype[,3] <- floor(phenotype[,3])
    phenotype[,4] <- floor(phenotype[,4])
  } else {
    if (identical(stratifyby, "ethnicity")) {
      if (phen_ncols != 5) stop("The number of columns in the phenotype data frame is not 5 as required.")
      phen_expected_columns <- c("sample_ids", "status", "ageInt", "AgeDiagIndex", "ethnicityClass")
      if (!all(sort(colnames(phenotype)) == sort(phen_expected_columns))) {
        stop("The phenotype data frame does not match the required column names which are: \nsample_ids, status, ageInt, AgeDiagIndex, ethnicityClass")
      }
      phenotype <- phenotype[ , c("sample_ids", "status", "ageInt", "AgeDiagIndex", "ethnicityClass")]
      phenotype[,3] <- floor(phenotype[,3])
      phenotype[,4] <- floor(phenotype[,4])
    } else {
      if (identical(stratifyby, "study")) {
        if (phen_ncols != 5) stop("The number of columns in the phenotype data frame is not 5 as required.")
        phen_expected_columns <- c("sample_ids", "status", "ageInt", "AgeDiagIndex", "study")
        if (!all(sort(colnames(phenotype)) == sort(phen_expected_columns))) {
          stop("The phenotype data frame does not match the required column names which are: \nsample_ids, status, ageInt, AgeDiagIndex, study")
        }
        phenotype <- phenotype[ , c("sample_ids", "status", "ageInt", "AgeDiagIndex", "study")]
        phenotype[,3] <- floor(phenotype[,3])
        phenotype[,4] <- floor(phenotype[,4])
      } else {
        if (is.null(stratifyby)) {
          if (ncol(phenotype) != 4) stop("The number of columns in the phenotype data frame is not 4 as required when no stratification is used.")
          phen_expected_columns <- c("sample_ids", "status", "ageInt", "AgeDiagIndex")
          if (!all(sort(colnames(phenotype)) == sort(phen_expected_columns))) {
            stop("The phenotype data frame must have columns: sample_ids, status, ageInt, AgeDiagIndex.")
          }
          phenotype <- phenotype[, phen_expected_columns]
          phenotype[,3] <- floor(phenotype[,3])
          phenotype[,4] <- floor(phenotype[,4])
        }
      }
    }
  }
  
  df <- merge(genotypes, phenotype, by="sample_ids")
  
  df$AgeDiagIndex[is.na(df$AgeDiagIndex)] <- 888
  df$ageInt[is.na(df$ageInt)] <- 888
  
  list <- list(colnames(df))
  
  colnum <- ncol(df) - 4
  total_variants <- colnum - 1
  
  results <- data.frame()

  df$age_pen <- assignAgePen(df$ageInt, df$AgeDiagIndex, df$status)
  df <- df[df$age_pen != 888, ]
  
  
  min_age <- agefilter[1]
  max_age <- agefilter[2]
  
  out_of_range <- df$age_pen < min_age | df$age_pen > max_age
  
  if (any(out_of_range, na.rm = TRUE)) {
    ndimdf <- nrow(df)
    df <- df[!out_of_range, ]
    del.df <- ndimdf - nrow(df)
    
    warning_message <- paste0(
      "This analysis is only applicable for samples diagnosed or interviewed ",
      "between the ages of ", min_age, " and ", max_age, ".\n",
      del.df, " sample(s) not following these criteria were excluded automatically."
    )
    
    warning(warning_message)
  }
  
  
  
  for (i in 2:colnum) {
    if (identical(stratifyby, "country")) df1 <- df[, c("sample_ids", list[[1]][i], "status", "ageInt", "AgeDiagIndex", "StudyCountry", "age_pen")]
    if (identical(stratifyby, "ethnicity")) df1 <- df[, c("sample_ids", list[[1]][i], "status", "ageInt", "AgeDiagIndex", "ethnicityClass", "age_pen")]
    if (identical(stratifyby, "study")) df1 <- df[, c("sample_ids", list[[1]][i], "status", "ageInt", "AgeDiagIndex", "study", "age_pen")]
    if (is.null(stratifyby)) df1 <- df[, c("sample_ids", list[[1]][i], "status", "ageInt", "AgeDiagIndex", "age_pen")]
    
    df1 <- df1 %>% drop_na(list[[1]][i])
    
    if (geno_notation == "n") {
      df1 <- df1[df1[,2]!="-1",]
    } else {
      if (geno_notation == "n/n") {
        df1 <- df1[df1[,2]!="./.",]
        df1[,2] <- ifelse(df1[,2] == "0/1", 1, 
                          ifelse(df1[,2] == "0/0", 0, 
                                 ifelse(df1[,2] == "1/1", 1, NA)))# make homos be hets
      }
    }
    
    
    variant <- paste(list[[1]][i])
  
    
    if(nrow(df1)>0){
      
      allcases <- df1[df1$status == 1 , ]
      allcontrols <- df1[df1$status == 0 , ]
      
      MAF_cases <- ( (nrow(df1[df1$status == 1 & df1[,2] == 1,])) + (2*nrow(df1[df1$status == 1 & df1[,2] == 2,])) ) / (2*nrow(allcases))
      MAF_controls <- ( (nrow(df1[df1$status == 0 & df1[,2] == 1,])) + (2*nrow(df1[df1$status == 0 & df1[,2] == 2,])) ) / (2*nrow(allcontrols))
      
      carriers <- df1[df1[[2]] == 1 & !is.na(df1[[2]]), ]
      
      carriers_cases <- carriers[carriers$status == 1 , ]
      carriers_controls <- carriers[carriers$status == 0 , ]
      
      if(nrow(carriers_cases)>0 & nrow(carriers_controls>0)){
        
        stats_case <- round( calculateStatistics(carriers_cases$age_pen), digits=2 )
        stats_control <- round( calculateStatistics(carriers_controls$age_pen), digits=2 )
        
        if (identical(stratifyby, "country")) {
          df1$StratifiedFactor <- as.factor(df1$StudyCountry)
        } else if (identical(stratifyby, "ethnicity")) {
          df1$StratifiedFactor <- as.factor(df1$ethnicityClass)
        } else if (identical(stratifyby, "study")) {
          df1$StratifiedFactor <- as.factor(df1$study)
        } else if (is.null(stratifyby)) {
          df1$StratifiedFactor <- factor("all")
        }
        
        df1$genotype <- df1[[ list[[1]][i] ]]
        
        # df1 <- df1 %>%
        #   dplyr::select(.data$status, .data$age_pen, .data$StratifiedFactor, .data$genotype) %>%
        #   tidyr::drop_na()
        
        df1 <- df1[, c("status", "age_pen", "StratifiedFactor", "genotype")]
        df1 <- df1[stats::complete.cases(df1), ]
        
        
        if (length(unique(df1$StratifiedFactor)) < 2) {
          model_full <- glm(status ~ genotype + age_pen, data = df1, family = binomial)
          model_reduced <- glm(status ~ age_pen, data = df1, family = binomial)
        } else {
          model_full <- glm(status ~ genotype + age_pen + StratifiedFactor, data = df1, family = binomial)
          model_reduced <- glm(status ~ age_pen + StratifiedFactor, data = df1, family = binomial)
        }
        
        LRT <- anova(model_reduced, model_full, test="LRT")
        
        pval_chisq <- LRT$"Pr(>Chi)"[2]
        
        log_odds_ratio <- as.numeric(model_full$coefficients[2])
        SE <- as.numeric(sqrt(diag(vcov(model_full)))[2])
        Odds_Ratio <- exp(log_odds_ratio)
        # CI <- suppressMessages(exp(confint(model_full)["genotype", ]))
        CI <- exp(log_odds_ratio + c(-1, 1) * 1.96 * SE)
        LCI <- CI[1]
        UCI <- CI[2]
        
        results[variant,"Variant ID"] <- variant
        results[variant,"Stratified By"] <- ifelse(is.null(stratifyby), "none", stratifyby)
        
        results[variant,"Total Carriers"] <- nrow(carriers)
        
        results[variant,"Case Carriers, N (%)"] <- paste(nrow(carriers_cases),"/",nrow(allcases), " (",formatC(MAF_cases, format = "e", digits = 2),")",sep="")
        
        if (nrow(carriers_cases) > 1) {
          results[variant,"Case Carriers Age"] <- paste(stats_case[1],"\u00b1",stats_case[2],
                                                        " (",stats_case[3],"-",stats_case[4],")",sep="")
        } else {
          if (nrow(carriers_cases) == 1) {
            results[variant,"Case Carriers Age"] <- stats_case[1]
          } else {
            results[variant,"Case Carriers Age"] <- "-"
          }
        } 
        
        
        results[variant,"Control Carriers, N (%)"] <- paste(nrow(carriers_controls),"/",nrow(allcontrols), " (",formatC(MAF_controls, format = "e", digits = 2),")",sep="")
        
        
        if (nrow(carriers_controls) > 1) {
          results[variant,"Control Carriers Age"] <- paste(stats_control[1],"\u00b1",stats_control[2],
                                                           " (",stats_control[3],"-",stats_control[4],")",sep="")
        } else {
          if (nrow(carriers_controls) == 1) {
            results[variant,"Control Carriers Age"] <- stats_control[1]
          } else { 
            results[variant,"Control Carriers Age"] <- "-"
          }
        }
        
        
        results[variant,"Logistic Regression Coefficient (Log Odds)"] <- log_odds_ratio
        results[variant,"Logistic Regression SE"] <- SE
        results[variant,"Logistic Regression OR"] <- Odds_Ratio
        results[variant,"Logistic Regression Lower Limit of the 95% CI"] <- LCI
        results[variant,"Logistic Regression Upper Limit of the 95% CI"] <- UCI
        results[variant,"Logistic Regression p-value"] <- pval_chisq
        
        results[variant, "PS4 criterion"] <- ps4_flag_one(
          gene = gene,
          OR   = Odds_Ratio,
          LCI  = LCI,
          UCI  = UCI,
          pval = pval_chisq,
          custom_rules = custom_rules
        )
        
        
        if (progress == TRUE) {
          cat("Variant", i-1, "out of", total_variants, "completed:", variant, "\n")
        }
        
      
      } else {
        if (progress == TRUE) {
          cat("Variant", i-1, "out of", colnum-1, "skipped (no carriers in cases or in controls):", variant, "\n")
        }
      } 
    } else {
      if (progress == TRUE) {
        cat("Variant", i-1, "out of", colnum-1, "skipped (no usable genotypes):", variant, "\n")
      }
    }
    
}
  
  rownames(results) <- NULL
  
  results$`Case Carriers Age` <- enc2utf8(as.character(results$`Case Carriers Age`)) # it forces underlying data to be encoded with 'utf-8'   
  results$`Control Carriers Age` <- enc2utf8(as.character(results$`Control Carriers Age`)) # it forces underlying data to be encoded with 'utf-8'   
  
  results$`Variant ID` <- gsub("chr","",results$`Variant ID`)
  
  print(results, quote = FALSE)
  
  
  if (exportcsv == TRUE) {
    write.csv(results, paste(data.folder.location, "/", output, ".csv", sep=""), row.names = FALSE)
  }
  
  results
  
}
