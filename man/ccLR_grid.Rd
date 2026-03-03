\name{ccLR.grid}

\alias{ccLR.grid}

\title{Refined Case-Control Likelihood Ratio (ccLR) Analysis for performing Grid Search by Scaling the Relative Risk}

\description{
  This function performs a grid-search case-control likelihood ratio analysis based on input genotype and phenotype data, optionally stratifying the results by country, ethnicity, or study. The function supports predefined or custom gene risk rates and allows for a grid scaling of the penetrance to investigate what magnitude of relative risk best fits the data. 
}

\usage{
  ccLR.grid(cancer = c("breast", "ovarian", "custom"),
            gene = c("BRCA1", "BRCA2", "PALB2", "CHEK2", "ATM", "TP53", "custom"),
            genotypes,
            geno_notation = c("n", "n/n"),
            phenotype,
	          grid = seq(0.5, 2, by=0.5),
            penetrance = c("Dorling", "Kuchenbaecker", "Antoniou", "Fortuno", 
                           "Li", "Hall", "Yang", "Momozawa", "custom"),
            custom_penetrance = NULL, 
            incidence_rate = c("England", "USA", "Japan", "Finland", "custom"),
            custom_incidence = NULL,
            outdir = NULL,
            output = "ccLR",
            stratifyby = NULL,
	          agefilter = c(0, 80),
            exportcsv = FALSE,
            progress = FALSE
  )
}

\arguments{
  \item{cancer}{
    A character string specifying the cancer type under investigation. Options are \code{"breast"}, \code{"ovarian"} or \code{"custom"} only.
  }
  \item{gene}{
    A character string specifying the gene of interest. Options are \code{"BRCA1"}, \code{"BRCA2"}, \code{"PALB2"}, \code{"CHEK2"}, \code{"ATM"}, \code{"TP53"} or \code{"custom"} only.
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
  \item{grid}{
    Optional. A vector of grid/scaling parameters that is applied to the age-specific relative risk curve. It represents how much more (or less) penetrant a specific variant may be compared to the average pathogenic variant for the same gene. For example: at 1 it assumes average gene-level pathogenicity, at 2 it assumes double the risk, and at 0.5 it assumes half the risk. Defaults to a sequence from 0.5 to 2 by 0.5 increments. 
  }
  \item{penetrance}{
    A character string specifying the penetrance method. Options are \code{"Dorling"}, \code{"Kuchenbaecker"}, \code{"Antoniou"}, \code{"Fortuno"}, \code{"Li"}, \code{"Hall"}, \code{"Yang"} or \code{"custom"}. Dorling contains breast rates for genes BRCA1, BRCA2, PALB2, CHEK2, and ATM. Kuchenbaecker contains breast and ovarian cancer rates for BRCA1 and BRCA2. Antoniou contains breast cancer rates for BRCA1, BRCA2, and PALB2. Fortuno and Li contain breast cancer rates for TP53. Hall contains ovarian cancer rates for ATM. Yang contains ovarian cancer rates for PALB2. If penetrance is set to \code{"custom"} the next argument \code{"custom_penetrance"} must be specified. 
  }
\item{custom_penetrance}{
  A data frame containing user-specified age-specific penetrance rates for variant carriers.
  Defaults to \code{NULL} but must be specified if \code{penetrance = "custom"}.
  
  The required column structure depends on the values of \code{cancer} and \code{gene}:
  
  \itemize{
    \item If \code{gene = "custom"}, the data frame must contain exactly two columns:
    \code{"Age"} and \code{"Penetrance_Carriers"}.
    
    \item If \code{cancer = "custom"} and \code{gene} is not \code{"custom"}, the data frame must contain exactly two columns:
    \code{"Age"} and \code{"Penetrance_Carriers_<gene>"}.
    
    \item If \code{cancer} is \code{"breast"} or \code{"ovarian"} and \code{gene} is not \code{"custom"}, the data frame must contain exactly two columns:
    \code{"Age"} and \code{"BC_Penetrance_Carriers_<gene>"} (for breast cancer) or
    \code{"OC_Penetrance_Carriers_<gene>"} (for ovarian cancer).
  }
  
  Column names are case-sensitive and no additional columns are permitted.
}
  \item{incidence_rate}{
   A character string specifying the population incidence rates to be used in the analysis. 
   Supported options are: \code{"England"}, \code{"USA"}, \code{"Japan"}, \code{"Finland"}, or \code{"custom"}.
   If incidence_rate is set to \code{"custom"} the next argument \code{"custom_incidence"} must be specified. 
   }
  \item{custom_incidence}{
  A data frame containing user-specified age-specific incidence rates.
  Defaults to \code{NULL} but must be specified if \code{incidence_rate = "custom"}.
  
  The data frame must contain exactly two columns:
  
  \itemize{
    \item \code{"Age"}: Age (in years).
    \item \code{"Incidence_rates"}: Population incidence rate at the corresponding age.
  }
  
  Column names are case-sensitive and no additional columns are permitted.
  }
  \item{outdir}{
Optional. A character string specifying the output directory. The default is set to NULL and in this case the output file containing the results is stored to a temporary file. To specify a permanent location this argument needs be specified.
  }
  \item{output}{
    Optional. A character string specifying the output file name. Defaults to \code{"ccLR"}.
  }
  \item{stratifyby}{
    Optional. A character string specifying the stratification variable. Options are \code{"country"}, \code{"ethnicity"}, or \code{"study"}, or NULL for single strata. The default entry is NULL. 
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
  The function implements a grid-search case-control likelihood ratio methodology for different genetic variants and optionally stratifies results by the specified variable. The grid search ccLR approach makes use of the grid scaling values (which are subject to choice) to scale the gene-specific average relative risk and identify what risk best fits the case-control data. Only samples diagnosed or interviewed between the ages of 21 and 80 are included in the analysis. The likelihood ratios derived are evaluated against the ACMG/AMP thresholds.
}

\value{
  A data frame containing the results of the case-control likelihood ratio analysis. If \code{exportcsv = TRUE}, the results are saved as a CSV file.
}

\examples{
  
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
  ccLR.grid(
    cancer = "breast",
    gene = "PALB2",
    genotypes = genotypes,
    geno_notation="n",
    phenotype = phenotype,
    penetrance = "Antoniou",
    incidence_rate = "England",
    stratifyby = "country",
    exportcsv = TRUE,
    progress = FALSE
  )
}


\author{
  Damianos Michaelides \email{damianosm@cing.ac.cy}, Maria Zanti, Christian Carrizosa, Theodora Nearchou, Kyriaki Michailidou
}

\references{
Antoniou, A., Pharoah, P. D. P., Narod, S., Risch, H. A., Eyfjord, J. E., Hopper, J. L., et al. (2003). Average risks of breast and ovarian cancer associated with BRCA1 or BRCA2 mutations detected in case Series unselected for family history: a combined analysis of 22 studies. Am. J. Hum. Genet. 72, 1117–1130.

Antoniou, A. C., Casadei, S., Heikkinen, T., Barrowdale, D., Pylkas, K., Roberts, J., ... and Tischkowitz, M. (2014). Breast-cancer risk in families with mutations in PALB2. New England Journal of Medicine, 371(6), 497-506.

Dorling, L. et al. (2021). Breast Cancer Risk Genes - Association Analysis in More than 113,000 Women. N Engl J Med 384, 428-439.

Fortuno, C., Feng, B. J., Carroll, C., Innella, G., Kohlmann, W., Lázaro, C., ..., and Spurdle, A. B. (2024). Cancer risks associated with TP53 pathogenic variants: Maximum likelihood analysis of extended pedigrees for diagnosis of first cancers beyond the Li-Fraumeni syndrome spectrum. JCO Precision Oncology, 8, e2300453.

Hall, M. J., Bernhisel, R., Hughes, E., Larson, K., Rosenthal, E. T., Singh, N. A., ... & Kurian, A. W. (2021). Germline pathogenic variants in the ataxia telangiectasia mutated (ATM) gene are associated with high and moderate risks for multiple cancers. Cancer Prevention Research, 14(4), 433-440.

Kuchenbaecker, K. B. J. L. Hopper, D. R. Barnes et al. (2017). Risks of breast, ovarian, and contralateral breast cancer for BRCA1 and BRCA2 mutation carriers. JAMA, vol. 317, no. 23, pp. 2402–2416.

Li, S., MacInnis, R. J., Lee, A., Nguyen-Dumont, T., Dorling, L., Carvalho, S., ..., and Antoniou, A. C. (2022). Segregation analysis of 17,425 population-based breast cancer families: evidence for genetic susceptibility and risk prediction. The American Journal of Human Genetics, 109(10), 1777-1788.

Yang, X., Leslie, G., Doroszuk, A., Schneider, S., Allen, J., Decker, B., ... & Tischkowitz, M. (2020). Cancer risks associated with germline PALB2 pathogenic variants: an international study of 524 families. Journal of clinical oncology, 38(7), 674-685.

Zanti, M. et al. (2023). A likelihood ratio approach for utilizing case-control data in the clinical classification of rare sequence variants: application to BRCA1 and BRCA2. Hum Mutat.

Zanti M et al. (2025). Analysis of more than 400,000 women provides case-control evidence for BRCA1 and BRCA2 variant classification. Nature Communications.

Momozawa, Y., Sasai, R., Usui, Y., Shiraishi, K., Iwasaki, Y., Taniyama, Y., ... & Kubo, M. (2022). Expansion of cancer risk profile for BRCA1 and BRCA2 pathogenic variants. JAMA oncology, 8(6), 871-878.
}

\keyword{case-control}
\keyword{likelihood ratio}
\keyword{case-control likelihood ratio}
\keyword{Breast Cancer research}
\keyword{genetic variants}

