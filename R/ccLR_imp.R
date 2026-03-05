ccLR.imp <- function(cancer = c("breast", "ovarian", "custom"),
                     gene = c("BRCA1", "BRCA2", "PALB2", "CHEK2", "ATM", "TP53", "custom"),
                     genotypes,
                     phenotype,
                     penetrance = c("Dorling", "Kuchenbaecker", "Antoniou", "Fortuno", "Li", "Hall", "Yang", "Momozawa", "custom"),
                     custom_penetrance = NULL, 
                     incidence_rate = c("England", "USA", "Japan", "Finland", "custom"),
                     custom_incidence = NULL,
                     outdir = NULL,
                     output = "ccLR",
                     stratifyby = NULL,
                     agefilter = c(0, 80),
                     exportcsv = FALSE,
                     progress = FALSE) {
  
  if (is.null(outdir)) {
    data.folder.location <- tempdir()
  } else {
    data.folder.location <- outdir
  }
  
  logical_args <- list(
    exportcsv = exportcsv,
    progress = progress
  )
  for (arg_name in names(logical_args)) {
    arg_value <- logical_args[[arg_name]]
    if (!is.logical(arg_value) || length(arg_value) != 1) {
      stop(sprintf("Argument '%s' must be a single logical value (TRUE or FALSE).", arg_name), call. = FALSE)
    }
  }
  
  
  if (!(cancer %in% c("breast", "ovarian", "custom"))) {
    stop("The cancer type (stated via the cancer argument) should be one of: breast or ovarian or custom.")
  }
  
  if ((cancer %in% "breast") &  !(gene %in% c("BRCA1", "BRCA2", "PALB2", "CHEK2", "ATM", "TP53", "custom"))) {
    stop("For breast cancer, the gene type (stated via the gene argument) should be one of: BRCA1 or BRCA2 or PALB2 or CHEK2 or ATM or TP53 or custom.")
  }
  
  if ((cancer %in% "ovarian") &  !(gene %in% c("BRCA1", "BRCA2", "PALB2", "ATM", "custom"))) {
    stop("For ovarian cancer, the gene type (stated via the gene argument) should be one of: BRCA1 or BRCA2 or PALB2 or ATM or custom.")
  }
  
  if ((cancer %in% "breast") & (gene %in% c("BRCA1", "BRCA2")) & !(penetrance %in% c("Dorling", "Antoniou", "Kuchenbaecker", "Momozawa", "custom")) ) {
    stop("The penetrance for your chosen cancer (stated via the cancer argument) and gene type (stated via the gene argument) is not available. Choose Dorling or Antoniou or Kuchenbaecker or Momozawa or provide custom rates.")
  }
  
  if ((cancer %in% "breast") & (gene %in% c("PALB2")) & !(penetrance %in% c("Dorling", "Antoniou", "custom")) ) {
    stop("The penetrance for your chosen cancer (stated via the cancer argument) and gene type (stated via the gene argument) is not available. Choose Dorling or Antoniou or provide custom penetrance.")
  }
  
  if ((cancer %in% "breast") & (gene %in% c("TP53")) & !(penetrance %in% c("Fortuno", "Li", "custom")) ) {
    stop("The penetrance for your chosen cancer (stated via the cancer argument) and gene type (stated via the gene argument) is not available. Choose Fortuno or Li or provide custom penetrance.")
  }
  
  if ((cancer %in% "breast") & (gene %in% c("CHEK2", "ATM")) & !(penetrance %in% c("Dorling", "custom")) ) {
    stop("The penetrance for your chosen gene type (stated via the gene argument) is not available. Choose Dorling or provide custom penetrance.")
  }
  
  if ((cancer %in% "ovarian") & (gene %in% c("BRCA1", "BRCA2")) & !(penetrance %in% c("Kuchenbaecker", "Momozawa", "custom")) ) {
    stop("The penetrance for your chosen cancer (stated via the cancer argument) and gene type (stated via the gene argument) is not available. Choose Kuchenbaecker or Momozawa or provide custom penetrance.")
  }
  
  if ((cancer %in% "ovarian") & (gene %in% c("PALB2")) & !(penetrance %in% c("Yang", "custom")) ) {
    stop("The penetrance for your chosen cancer (stated via the cancer argument) and gene type (stated via the gene argument) is not available. Choose Yang or provide custom penetrance.")
  }
  
  if ((cancer %in% "ovarian") & (gene %in% c("ATM")) & !(penetrance %in% c("Hall", "custom")) ) {
    stop("The penetrance for your chosen cancer (stated via the cancer argument) and gene type (stated via the gene argument) is not available. Choose Hall or provide custom penetrance.")
  }
  
  
  if (!(penetrance %in% c("Dorling", "Kuchenbaecker", "Antoniou", "Fortuno", "Li", "Hall", "Yang", "Momozawa", "custom"))) {
    stop("The penetrance (stated via the penetrance argument) should be one of: Dorling or Kuchenbaecker or Antoniou or Fortuno or Li or Hall or Yang or Momozawa or custom.")
  }
  
  if (!(incidence_rate %in% c("England", "USA", "Japan", "Finland", "custom"))) {
    stop("The incidence rates (stated via the incidence_rate argument) should be one of: England or USA or Japan or Finland or custom.")
  }
  
  if (!is.null(stratifyby) && !(stratifyby %in% c("country", "ethnicity", "study"))) {
    stop("The stratification must be by country, ethnicity, study, or NULL (no stratification). Revise your stratifyby choice.")
  }
  
  if ( identical(cancer, "custom") & !(identical(penetrance, "custom")) ) {
    stop("The cancer type (stated via the cancer argument) is set to custom, hence, the penetrance (stated via the penetrance argument) must be set as custom. \n Note that, a custom_penetrance dataset must be set via the custom_penetrance argument.")
  }
  
  if ( identical(gene, "custom") & !(identical(penetrance, "custom")) ) {
    stop("The gene type (stated via the gene argument) is set to custom, hence, the penetrance (stated via the penetrance argument) must be set as custom. \n Note that, a custom_penetrance dataset must be set via the custom_penetrance argument.")
  }
  
  
  if ( identical(penetrance, "custom") & is.null(custom_penetrance) ) {
    stop("The penetrance (stated via the penetrance argument) is set to custom, but custom rates data were not provided in the custom_penetrance argument.")
  }
  
  if ( identical(incidence_rate, "custom") & is.null(custom_incidence) ) {
    stop("The incidence rates (stated via the incidence_rate argument) is set to custom, but custom incidence rates data were not provided in the custom_incidence argument.")
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
  
  
  if (identical(penetrance, "custom")) {
    
    breast_genes  <- c("BRCA1", "BRCA2", "ATM", "PALB2", "CHEK2", "TP53", "custom")
    ovarian_genes <- c("BRCA1", "BRCA2", "PALB2", "ATM", "custom")
    prefix <- ifelse(cancer == "breast", "BC", ifelse(cancer == "ovarian", "OC", "Other"))
    
    expected_cols <- if (gene == "custom") {
      c("Age", "Penetrance_Carriers")
    } else if (cancer == "custom") {
      c("Age", paste0("Penetrance_Carriers_", gene))
      
    } else {
      c("Age", paste0(prefix, "_Penetrance_Carriers_", gene))
    }
    
    if (!all(sort(colnames(custom_penetrance)) == sort(expected_cols))) {
      stop(
        "Custom penetrance data have incorrect column names.\n",
        "Expected columns: ", paste(expected_cols, collapse = ", "), "\n",
        "Provided columns: ", paste(colnames(custom_penetrance), collapse = ", ")
      )
    }
  }
  
  if (identical(incidence_rate, "custom")) {
    expected_cols <- c("Age", "Incidence_rates")
    if (!all(sort(colnames(custom_incidence)) == sort(expected_cols))) {
      stop(
        "Custom incidence rates data have incorrect column names.\n",
        "Expected columns: ", paste(expected_cols, collapse = ", "), "\n",
        "Provided columns: ", paste(colnames(custom_incidence), collapse = ", ")
      )
    }
  }
  
  
  
  if (identical(penetrance, "custom")) {
    gene_rates <- custom_penetrance
    colnames(gene_rates) <- c("Age", "RR")
  } else {
    prefix <- ifelse(cancer == "breast", "BC", "OC")
    penetrance_obj_name <- paste0(prefix, "_Penetrance_Carriers_", gene)
    penetrance <- get(penetrance)
    gene_rates <- penetrance[, c("Age", penetrance_obj_name)]
    # gene_rates <- data.frame(penetrance$Age, penetrance[[penetrance_obj_name]])
    colnames(gene_rates) <- c("Age", "RR")
  }
  
  
  if (identical(incidence_rate, "custom")) {
    inc_rates <- custom_incidence
    colnames(inc_rates) <- c("Age", "IR")
  } else {
    prefix <- ifelse(cancer == "breast", "BC", ifelse(cancer == "ovarian", "OC", "Other"))
    incedent_obj_name <- paste0(prefix, "_", incidence_rate)
    
    inc_rates <- data.frame(incidence_data$Age, incidence_data[[incedent_obj_name]])
    colnames(inc_rates) <- c("Age", "IR")
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
  stratified_result <- data.frame()
  
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
    
    variant <- paste(list[[1]][i])
    
    if(nrow(df1)>0){
      
      ## for imputed: column is dosage (0..2)
      orig_vals <- df1[,2]
      df1[,2] <- suppressWarnings(as.numeric(df1[,2]))
      na_conv <- sum(is.na(df1[,2]) & !is.na(orig_vals))
      if (na_conv > 0) {
        warning(
          na_conv,
          " genotype value(s) were non-numeric and converted to NA. ",
          "These samples will be excluded from the ccLR calculation."
        )
      }
      
      df1 <- df1[!is.na(df1[,2]), ]   
      
      bad_dosage <- df1[,2] < 0 | df1[,2] > 2
      
      if (any(bad_dosage, na.rm = TRUE)) {
        warning(
          sum(bad_dosage, na.rm = TRUE),
          " genotype dosage value(s) were outside the expected range [0,2] ",
          "and were excluded from the analysis."
        )
        df1 <- df1[!bad_dosage, ]
      }
      
      allcases <- df1[df1$status == 1 , ]
      allcontrols <- df1[df1$status == 0 , ]
      
      ## imputed MAF from dosage
      MAF_cases <- mean(allcases[,2], na.rm = TRUE) / 2
      MAF_controls <- mean(allcontrols[,2], na.rm = TRUE) / 2
      
      ## expected carrier count (sum of dosages)
      K_all <- sum(df1[,2], na.rm = TRUE)
      
      if(K_all > 0){
        
        ## expected counts (sum of dosages) by case/control
        K_cases <- sum(allcases[,2], na.rm = TRUE)
        K_controls <- sum(allcontrols[,2], na.rm = TRUE)
        
        ## weighted age summaries (weights = dosage)
        stats_case <- round(wstats(allcases$age_pen, allcases[,2]), digits=2)
        stats_control <- round(wstats(allcontrols$age_pen, allcontrols[,2]), digits=2)
        
        if (identical(stratifyby, "country")) df1$StratifiedFactor <- as.factor(df1$StudyCountry)
        if (identical(stratifyby, "ethnicity")) df1$StratifiedFactor <- as.factor(df1$ethnicityClass)
        if (identical(stratifyby, "study")) df1$StratifiedFactor <- as.factor(df1$study)
        if (is.null(stratifyby)) df1$StratifiedFactor <- factor("all")  
        
        strat <- levels(df1$StratifiedFactor)
        
        for (h in seq_along(strat)) {
          
          df2 <- df1[which(df1$StratifiedFactor == strat[h]), ]
          N <- nrow(df2)
          
          m2 <- df2
          rownames(m2) <- NULL
          
          hazard <- inc_rates$IR            
          cumhazard <- cumsum(hazard)             
          h_carrier <- inc_rates$IR * gene_rates$RR   
          cumhazard_carrier <- cumsum(h_carrier)
          
          m2$RR <- gene_rates$RR[ match(m2$age_pen, gene_rates$Age) ]
          m2$S0 <- exp( - cumhazard[m2$age_pen] )
          m2$S1 <- exp( - cumhazard_carrier[m2$age_pen] )
          
          ## imputed ccLR
          d_i <- as.numeric(m2[,2])              # dosages
          K_star <- sum(d_i, na.rm = TRUE)       # expected carriers (sum dosages)
          
          m2$Likelihood <- calculateLikelihood(m2$S0, m2$S1, m2$RR, m2$status)
          
          if (K_star > 0) {
            log_num <- sum(d_i * log(m2$Likelihood))
            log_denom <- K_star * log(mean(m2$Likelihood))
            log_LR <- log_num - log_denom
            LR <- exp(log_LR)
          } else {
            LR <- 1
          }
          
          stratified_result[strat[h],"LR_stratum"] <- LR 
        }
        
        results[variant,"Variant ID"] <- variant
        results[variant,"Stratified By"] <- ifelse(is.null(stratifyby), "none", stratifyby)
        
        ## keep same column name; now it's expected carriers (sum dosages)
        results[variant,"Total Carriers"] <- K_all
        
        results[variant,"Case Carriers, N (%)"] <- paste(K_cases,"/",nrow(allcases), " (",formatC(MAF_cases, format = "e", digits = 2),")",sep="")
        
        if (is.finite(stats_case[1]) && is.finite(stats_case[2])) {
          results[variant,"Case Carriers Age"] <- paste(stats_case[1],"\u00b1",stats_case[2],
                                                        " (",stats_case[3],"-",stats_case[4],")",sep="")
        } else {
          results[variant,"Case Carriers Age"] <- "-"
        }
        
        results[variant,"Control Carriers, N (%)"] <- paste(K_controls,"/",nrow(allcontrols), " (",formatC(MAF_controls, format = "e", digits = 2),")",sep="")
        
        if (is.finite(stats_control[1]) && is.finite(stats_control[2])) {
          results[variant,"Control Carriers Age"] <- paste(stats_control[1],"\u00b1",stats_control[2],
                                                           " (",stats_control[3],"-",stats_control[4],")",sep="")
        } else { 
          results[variant,"Control Carriers Age"] <- "-"
        }
        
        finalLR <- prod(stratified_result$LR_stratum)
        
        if (identical(finalLR, "NaN")) {
          results[variant,"LR"] <- "-"
          results[variant,"ACMG/AMP Evidence"] <- "-"
        } else {
          if( is.na(finalLR) ) {
            results[variant,"LR"] <- "-"
            results[variant,"ACMG/AMP Evidence"] <- "-"
          } else {
            results[variant,"LR"] <- finalLR
            if (finalLR >= 350) results[variant,"ACMG/AMP Evidence"] <- "Pathogenic Very Strong"
            if (finalLR >= 18.7 & finalLR < 350) results[variant,"ACMG/AMP Evidence"] <- "Pathogenic Strong"
            if (finalLR >= 4.33 & finalLR < 18.7) results[variant,"ACMG/AMP Evidence"] <- "Pathogenic Moderate"
            if (finalLR >= 2.08 & finalLR < 4.33) results[variant,"ACMG/AMP Evidence"] <- "Pathogenic Supporting"
            if (finalLR <= 0.0029) results[variant,"ACMG/AMP Evidence"] <- "Bening Very Strong"
            if (finalLR > 0.0029 & finalLR <= 0.053) results[variant,"ACMG/AMP Evidence"] <- "Bening Strong"
            if (finalLR > 0.053 & finalLR <= 0.231) results[variant,"ACMG/AMP Evidence"] <- "Bening Moderate"
            if (finalLR > 0.231 & finalLR <= 0.48) results[variant,"ACMG/AMP Evidence"] <- "Bening Supporting"
            if (finalLR < 2.08 & finalLR > 0.48) results[variant,"ACMG/AMP Evidence"] <- "No Evidence"
          }
        }
        
        if (progress == TRUE) {
          cat("Variant", i-1, "out of", total_variants, "completed:", variant, "\n")
        }
        
      } else {
        if (progress == TRUE) {
          cat("Variant", i-1, "out of", colnum-1, "skipped (no dosage > 0):", variant, "\n")
        }
      } 
    } else {
      if (progress == TRUE) {
        cat("Variant", i-1, "out of", colnum-1, "skipped (no usable genotypes):", variant, "\n")
      }
    }
    
  }
  
  rownames(results) <- NULL
  
  results$`Case Carriers Age` <- enc2utf8(as.character(results$`Case Carriers Age`))
  results$`Control Carriers Age` <- enc2utf8(as.character(results$`Control Carriers Age`))
  
  results$`Variant ID` <- gsub("chr","",results$`Variant ID`)
  
  print(results, quote = FALSE)
  
  if (exportcsv == TRUE) {
    write.csv(results, paste(data.folder.location, "/", output, ".csv", sep=""), row.names = FALSE)
  }
  
  results
}





