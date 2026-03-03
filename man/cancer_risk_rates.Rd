\name{risk_rates}

\alias{risk_rates}

\alias{Dorling}
\alias{Kuchenbaecker}
\alias{Antoniou}
\alias{Fortuno}
\alias{Li}
\alias{Hall}
\alias{Yang}
\alias{Momozawa}


\docType{data}

\title{
  Breast and Ovarian Cancer Risk Rates: Dorling et al. (2021), Kuchenbaecker et al. (2017), Antoniou et al. (2003), Fortuno et al. (2024), Li et al. (2022), Hall et al. (2021), Yang et al. (2020), and Momozawa et al. (2022)
}

\description{
  These datasets provide age-specific disease penetrances (relative risks) for breast and ovarian cancers. The datasets are derived from Dorling et al. (2021), Kuchenbaecker et al. (2017), Antoniou et al. (2003), Antoniou et al. (2014), Fortuno et al. (2024), Li et al. (2022), Hall et al. (2021), Yang et al. (2020), and Momozawa et al. (2022). The datasets are used for the calculation of the case-control likelihood ratio (ccLR) analyses of the BRCA1, BRCA2, PALB2, CHEK2, ATM, and TP53 genetic variants. Dorling contains breast cancer rates for genes BRCA1, BRCA2, PALB2, CHEK2, and ATM. Kuchenbaecker contains breast and ovarian cancer rates for BRCA1 and BRCA2. Antoniou contains breast cancer rates for BRCA1, BRCA2, and PALB2. Fortuno and Li contain breast cancer rates for TP53. Hall contains ovarian cancer rates for ATM. Yang contains ovarian cancer rates for PALB2. Momozawa contains breast and ovarian cancer rates for BRCA1 and BRCA2.
}

\usage{
  Dorling
  Kuchenbaecker
  Antoniou
  Fortuno
  Li
  Hall 
  Yang
  Momozawa
}

\details{
  Each dataset is a data frame that contains:
    
    \describe{
      \item{\code{Age}}{Numeric. The age range or specific ages.}
      \item{\code{Relative_risk}}{Numeric. The relative risk for carriers of the pre-mentioned gene variants.}
    }
}

\references{
Antoniou, A., Pharoah, P. D. P., Narod, S., Risch, H. A., Eyfjord, J. E., Hopper, J. L., et al. (2003). Average risks of breast and ovarian cancer associated with BRCA1 or BRCA2 mutations detected in case Series unselected for family history: a combined analysis of 22 studies. Am. J. Hum. Genet. 72, 1117–1130.

Antoniou, A. C., Casadei, S., Heikkinen, T., Barrowdale, D., Pylkas, K., Roberts, J., ... and Tischkowitz, M. (2014). Breast-cancer risk in families with mutations in PALB2. New England Journal of Medicine, 371(6), 497-506.

Dorling, L. et al. (2021). Breast Cancer Risk Genes - Association Analysis in More than 113,000 Women. N Engl J Med 384, 428-439.

Fortuno, C., Feng, B. J., Carroll, C., Innella, G., Kohlmann, W., Lázaro, C., ..., and Spurdle, A. B. (2024). Cancer risks associated with TP53 pathogenic variants: Maximum likelihood analysis of extended pedigrees for diagnosis of first cancers beyond the Li-Fraumeni syndrome spectrum. JCO Precision Oncology, 8, e2300453.

Hall, M. J., Bernhisel, R., Hughes, E., Larson, K., Rosenthal, E. T., Singh, N. A., ... & Kurian, A. W. (2021). Germline pathogenic variants in the ataxia telangiectasia mutated (ATM) gene are associated with high and moderate risks for multiple cancers. Cancer Prevention Research, 14(4), 433-440.

Kuchenbaecker, K. B. J. L. Hopper, D. R. Barnes et al. (2017). Risks of breast, ovarian, and contralateral breast cancer for BRCA1 andBRCA2 mutation carriers. JAMA, vol. 317, no. 23, pp. 2402–2416.

Li, S., MacInnis, R. J., Lee, A., Nguyen-Dumont, T., Dorling, L., Carvalho, S., ..., and Antoniou, A. C. (2022). Segregation analysis of 17,425 population-based breast cancer families: evidence for genetic susceptibility and risk prediction. The American Journal of Human Genetics, 109(10), 1777-1788.

Yang, X., Leslie, G., Doroszuk, A., Schneider, S., Allen, J., Decker, B., ... & Tischkowitz, M. (2020). Cancer risks associated with germline PALB2 pathogenic variants: an international study of 524 families. Journal of clinical oncology, 38(7), 674-685.

Momozawa, Y., Sasai, R., Usui, Y., Shiraishi, K., Iwasaki, Y., Taniyama, Y., ... & Kubo, M. (2022). Expansion of cancer risk profile for BRCA1 and BRCA2 pathogenic variants. JAMA oncology, 8(6), 871-878.
}




\examples{
  ## Load the Dorling dataset
  data(Dorling)
  head(Dorling)
  
  ## Load the Kuchenbaecker dataset
  data(Kuchenbaecker)
  head(Kuchenbaecker)
  
  ## Load the Antoniou dataset
  data(Antoniou)
  head(Antoniou)

  ## Load the Fortuno dataset
  data(Fortuno)
  head(Fortuno)
    
  ## Load the Li dataset
  data(Li)
  head(Li) 

  ## Load the Hall dataset
  data(Hall)
  head(Hall) 

  ## Load the Yang dataset
  data(Yang)
  head(Yang) 

  ## Load the Momozawa dataset
  data(Momozawa)
  head(Momozawa) 
}


\keyword{breast cancer research}
\keyword{ovarian cancer research}

