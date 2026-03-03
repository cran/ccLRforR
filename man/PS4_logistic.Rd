\name{ps4.logistic}

\alias{ps4.logistic}

\title{Logistic Regression PS4 Criterion}

\description{
  This function performs logistic regression, calculates the likelihood ratio test, the odds ratio, and confidence intervals around it to compare against gene-specific PS4 criteria. The function performs based on input genotype and phenotype data. The factors assessed in the model are the ages and an optional stratification factor which is either country or ethnicity.
}

\usage{
  ps4.logistic(
    gene = c("BRCA1", "BRCA2", "PALB2", "CHEK2", "ATM", "TP53", "custom"),
    genotypes,
    geno_notation = c("n", "n/n"),
    phenotype,
    custom_rules = NULL,
    outdir = NULL,
    output = "PS4",
    stratifyby = NULL,
    agefilter = c(0, 80),
    exportcsv = FALSE,
    progress = FALSE
  )
}

\arguments{
  \item{gene}{
    A character string specifying the gene of interest. Options are \code{"BRCA1"}, \code{"BRCA2"}, \code{"PALB2"}, \code{"CHEK2"}, \code{"ATM"}, \code{"TP53"} and \code{custom}.
  }
  \item{genotypes}{
    A data frame containing genotype data with the first column named \code{"sample_ids"} and subsequent columns for genotype information.
  }
    \item{geno_notation}{
    A character string specifying the format of the genotypes notation. Options are  \code{"n"}, or \code{"n/n"} only. In context, if variants take entries 0 (homozygous reference), 1 (heterozygous), 2 (homozygous alternate), and -1 (missing) then choose \code{geno_notation}=\code{"n"}. Alternatively, if variants take entries 0/0 (homozygous reference), 0/1 (heterozygous), 1/1 (homozygous alternate), and ./. (missing) then choose \code{geno_notation}=\code{"n/n"}. For other formats, please tranform your dataset to one of the accepted/implemented formats. 
  }
  \item{phenotype}{
    A data frame containing phenotype data. The required columns depend on the \code{stratifyby} parameter. If single strata is considered, i.e., if \code{stratifyby}=NULL, the data frame must include columns \code{"sample_ids"}, \code{"status"}, \code{"ageInt"}, \code{"AgeDiagIndex"}. If stratification is considered, the data frame must have an additional stratification column (\code{"StudyCountry"}, \code{"ethnicityClass"}, or \code{"study"}) depending on the stratification variable.
  }
\item{custom_rules}{
    Optional. A named list of functions that define user-specified PS4 decision rules for one or more genes. Each function must return `"Yes"` or `"No"` when evaluated, and will be passed the arguments `OR`, `LCI`, `UCI`, and `pval`. By default, hard-coded thresholds for BRCA1, BRCA2, ATM, CHEK2, PALB2, and TP53 are applied (see Details). Supplying a `custom_rules` list allows users to: (a) Override the default criteria for one or more of these genes, and (b) Define thresholds for `"custom"` genes. Check the Examples section for an example. 
  }
  \item{outdir}{
    Optional. A character string specifying the output directory. The default is set to NULL and in this case the output file containing the results is stored to a temporary file. To specify a permanent location this argument needs be specified.
  }
  \item{output}{
    Optional. A character string specifying the output file name. Defaults to \code{"PS4"}.
  }
  \item{stratifyby}{
    A character string specifying the stratification variable. Options are \code{"country"}, \code{"ethnicity"}, or \code{"study"}, or NULL for single strata. The default entry is NULL. 
  }
\item{agefilter}{
    A numeric vector of length 2 specifying the age range to include in the analysis. Defaults to ages 0 to 80. 
  }
    \item{exportcsv}{
    Optional. A logical value indicating whether to export the results as a CSV file (on top of printing the results in R). Defaults to \code{FALSE}.
  }
    \item{progress}{
      Optional. If \code{TRUE}, it returns the progress of the variants analysed. The default entry is FALSE.  
    }
}

\details{
The function implements the case-control likelihood ratio methodology for different genetic variants and stratifies results by the specified variable. It validates inputs, applies the calculations based on the chosen method, and generates a summary of the results. Only samples diagnosed or interviewed between the ages of 21 and 80 are included in the analysis.

The function implements ClinGen-specified, gene-specific criteria for applying the ACMG/AMP rule PS4 (case–control evidence of pathogenicity). It evaluates each variant using the odds ratio (OR), relative risk (RR), Wald confidence interval (CI), and association test p-value from logistic regression, and then applies thresholds that differ by gene. For BRCA1/2, PS4 is assigned when p <= 0.05, OR >= 4, and the 95% CI excludes 2.0. For ATM and CHEK2, PS4 is assigned when p <= 0.05 and either OR >= 2 or the lower 95% CI bound is >= 1.5. For PALB2, PS4 is assigned when p-value <= 0.05 and either OR >= 3 or the lower 95% CI bound is >= 1.5. For TP53, PS4 is assigned when p <= 0.05, the OR is > 5.0, and the 95% CI does not include 1.0. Variants not meeting these conditions, or lacking sufficient statistical information, are classified as not fulfilling PS4.
}

\value{
  A data frame containing the results of the logistic regression and likelihood ratio test analysis, evaluated against the PS4 criteria. If \code{exportcsv = TRUE}, the results are saved as a CSV file.
}

\examples{
  
  ## Example 1:
  ## Define simulated inputs - genotypes and phenotype
  
  genotypes <- data.frame(
    sample_ids = 1:100,
    variant1 = rbinom(100, 2, 0.3),
    variant2 = rbinom(100, 2, 0.2)
  )
  
  phenotype <- data.frame(
    sample_ids = 1:100,
    status = rbinom(100, 1, 0.5),
    ageInt = floor(runif(100, 21, 80)),
    AgeDiagIndex = floor(runif(100, 21, 80)),
    StudyCountry = sample(c("USA", "UK", "Canada"), 100, replace = TRUE)
  )
  
  # Run the function

  ps4.logistic(
    gene = "CHEK2",
    genotypes = genotypes,
    geno_notation="n",
    phenotype = phenotype,
    stratifyby = "country",
    exportcsv = TRUE, 
    progress = FALSE
  )


  ## Example 2:
  ## Define simulated inputs - genotypes and phenotype
  
  genotypes <- data.frame(
    sample_ids = 1:100,
    variantX = rbinom(100, 2, 0.1)
  )
  
  phenotype <- data.frame(
    sample_ids = 1:100,
    status = rbinom(100, 1, 0.5),
    ageInt = floor(runif(100, 21, 80)),
    AgeDiagIndex = floor(runif(100, 21, 80)),
    ethnicityClass = sample(c("European", "Asian", "African"), 100, replace = TRUE)
  )
  
  ## Define a custom rule for a "custom" gene:
  ### Flag "Yes" if OR >= 2.5 and CI lower bound >= 1.2

  custom_rules <- list(
    CUSTOM = function() ifelse(OR >= 2.5 && LCI >= 1.2, "Yes", "No")
  )
  
  ## Run the function
  ps4.logistic(
    gene = "custom",
    genotypes = genotypes,
    geno_notation = "n",
    phenotype = phenotype,
    custom_rules = custom_rules,
    stratifyby = "ethnicity",
    exportcsv = FALSE,
    progress = TRUE,
  )


}


\author{
  Damianos Michaelides \email{damianosm@cing.ac.cy}, Maria Zanti, Christian Carrizosa, Theodora Nearchou, Kyriaki Michailidou
}

\references{
Parsons, M. T. et al. Evidence-based recommendations for gene-specific ACMG/AMP variant classification from the ClinGen ENIGMA BRCA1 and BRCA2 Variant Curation Expert Panel. Am J Hum Genet (2024).

Zanti, M. et al. (2023). A likelihood ratio approach for utilizing case-control data in the clinical classification of rare sequence variants: application to BRCA1 and BRCA2. Hum Mutat.

Zanti M et al. (2025). Analysis of more than 400,000 women provides case-control evidence for BRCA1 and BRCA2 variant classification. Nature Communications.
}

\keyword{case-control}
\keyword{likelihood ratio}
\keyword{case-control likelihood ratio}
\keyword{Breast Cancer research}
\keyword{genetic variants}

