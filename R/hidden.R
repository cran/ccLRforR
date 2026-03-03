ps4_flag_one <- function(gene, OR, LCI, UCI, pval,
                         custom_rules = NULL) {
  g <- toupper(gsub("\\s+", "", gene %||% ""))  
  if (is.na(OR) || is.na(LCI) || is.na(UCI) || is.na(pval)) return("No")
  if (!(pval <= 0.05)) return("No")
  ci_excludes <- function(value, lci, uci) !(lci <= value && uci >= value)
  
  # Default rules
  default_rules <- list(
    BRCA1 = function() ifelse(OR >= 4 && ci_excludes(2.0, LCI, UCI), "Yes", "No"),
    BRCA2 = function() ifelse(OR >= 4 && ci_excludes(2.0, LCI, UCI), "Yes", "No"),
    ATM   = function() ifelse(OR >= 2 || LCI >= 1.5, "Yes", "No"),
    CHEK2 = function() ifelse(OR >= 2 || LCI >= 1.5, "Yes", "No"),
    PALB2 = function() ifelse(OR >= 3 || LCI >= 1.5, "Yes", "No"),
    TP53  = function() ifelse(OR > 5 && ci_excludes(1.0, LCI, UCI), "Yes", "No")
  )
  
  # Merge in user-provided rules (overrides defaults)
  if (!is.null(custom_rules)) {
    default_rules[names(custom_rules)] <- custom_rules
  }
  
  # Apply relevant rule
  if (g %in% names(default_rules)) {
    return(default_rules[[g]]())
  } else {
    return("No")
  }
}


utils::globalVariables(c("incidence_data"))




