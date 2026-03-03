\name{ccLRforR}

\alias{ccLRforR-package}
\alias{ccLRforR}

\docType{package}

\title{
\packageTitle{ccLRforR}
}

\description{
This package provides tools for implementing the case-control likelihood ratio (ccLR) analyses and logistic regression applying the PS4 criteria for genetic data, supporting optional stratification by country, ethnicity, or study. It includes functionality for built-in or custom cancer gene risk rates.

The package is designed for researchers analysing genetic breast, ovarian, or custom cancer data using the standard and grid-search ccLR approaches or a logistic regression applying the PS4 criteria. 
  
The ccLR method compares the likelihood of the distribution of the variant of interest among cases and controls, under the hypothesis that the variant is associated with similar risks of the disease in question, as the "average" pathogenic variant, compared to the likelihood under the  hypothesis that it is a benign variant not associated with increased risk.

The grid search ccLR approach makes use of the grid scaling values (which are subject to choice) to scale the gene-specific average relative risk and identify what risk best fits the case-control data.
    
The package includes functions for:
    
    \itemize{
      \item ccLR analysis for built-in gene-specific and age-specific penetrances (\code{\link{ps4.ccLR}}).
      \item Grid search ccLR analysis that uses grid parameters to scale the penetrance levels and find which best fits the data (\code{\link{ccLR.grid}}).
      \item Likelihood Ratio Test via Logistic Regression to evaluate variants' pathogenicity against the PS4 criteria (\code{\link{ps4.logistic}}). 
    }
  
  Key features:
    \itemize{
      \item Optional stratification by country, ethnicity, or study.
      \item Built-in penetrances by Dorling et. al., 2021, Kuchenbaecker et. al., 2017, Antoniou et. al., 2003, Fortuno et. al., 2024, Li et. al., 2022, Hall et. al., 2021, and Yang et. al., 2020.
      \item Support for custom cancer types, genes penetrances, and incidence rates.
      \item Exclusion of samples outside the age range 21-80.
    }
  
}

\author{
\packageAuthor{ccLRforR}

Maintainer: \packageMaintainer{ccLRforR}
}


\references{
Antoniou, A., Pharoah, P. D. P., Narod, S., Risch, H. A., Eyfjord, J. E., Hopper, J. L., et al. (2003). Average risks of breast and ovarian cancer associated with BRCA1 or BRCA2 mutations detected in case Series unselected for family history: a combined analysis of 22 studies. Am. J. Hum. Genet. 72, 1117–1130.

Antoniou, A. C., Casadei, S., Heikkinen, T., Barrowdale, D., Pylkas, K., Roberts, J., ... and Tischkowitz, M. (2014). Breast-cancer risk in families with mutations in PALB2. New England Journal of Medicine, 371(6), 497-506.

Dorling, L. et al. (2021). Breast Cancer Risk Genes - Association Analysis in More than 113,000 Women. N Engl J Med 384, 428-439.

Fortuno, C., Feng, B. J., Carroll, C., Innella, G., Kohlmann, W., Lázaro, C., ..., and Spurdle, A. B. (2024). Cancer risks associated with TP53 pathogenic variants: Maximum likelihood analysis of extended pedigrees for diagnosis of first cancers beyond the Li-Fraumeni syndrome spectrum. JCO Precision Oncology, 8, e2300453.

Kuchenbaecker, K. B. J. L. Hopper, D. R. Barnes et al. (2017). Risks of breast, ovarian, and contralateral breast cancer for BRCA1 andBRCA2 mutation carriers. JAMA, vol. 317, no. 23, pp. 2402–2416.

Li, S., MacInnis, R. J., Lee, A., Nguyen-Dumont, T., Dorling, L., Carvalho, S., ..., and Antoniou, A. C. (2022). Segregation analysis of 17,425 population-based breast cancer families: evidence for genetic susceptibility and risk prediction. The American Journal of Human Genetics, 109(10), 1777-1788.

Parsons, M. T. et al. (2024). Evidence-based recommendations for gene-specific ACMG/AMP variant classification from the ClinGen ENIGMA BRCA1 and BRCA2 Variant Curation Expert Panel. Am J Hum Genet.

Zanti, M. et al. (2023). A likelihood ratio approach for utilizing case-control data in the clinical classification of rare sequence variants: application to BRCA1 and BRCA2. Hum Mutat.

Zanti M et al. (2025). Analysis of more than 400,000 women provides case-control evidence for BRCA1 and BRCA2 variant classification. Nature Communications.
}

\seealso{
  \code{\link{ps4.ccLR}}, \code{\link{ps4.logistic}}, \code{\link{ccLR.grid}}
}

\keyword{case-control}
\keyword{likelihood ratio}
\keyword{case-control likelihood ratio}
\keyword{Breast Cancer research}
\keyword{genetic analysis}
