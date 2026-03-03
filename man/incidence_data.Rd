\name{incidence_data}

\alias{incidence_data}


\docType{data}

\title{Age-specific population incidence rates for breast and ovarian cancer}


\description{
  Age-specific population incidence rates used to compute baseline hazards in the likelihood ratio model.
  
  The incidence rates data are stored with one row per age (from 0 to 80 years) and separate columns for each cancer type and population.
  
  The first column is:
  \itemize{
    \item \code{"Age"}: Age in years.
  }
  
  The remaining columns follow the naming convention:
  \code{"<Cancer>_<Population>"},
  where:
  
  \itemize{
    \item \code{<Cancer>} is \code{"BC"} (breast cancer) or \code{"OC"} (ovarian cancer)
    \item \code{<Population>} is one of \code{"England"}, \code{"USA"}, \code{"Finland"}, or \code{"Japan"}
  } 

  Each cell contains the annual incidence rate for the corresponding age, cancer type,
  and population. These rates are used to construct cumulative baseline hazards
  in the likelihood calculations.

}

\usage{
incidence_data
}


\format{
A data frame with 81 rows and 9 columns containing age-specific annual
population incidence rates.

\describe{
  \item{Age}{Age in years (integer, ranging from 0 to 80).}
  \item{BC_England}{Breast cancer incidence rate for England at the given age.}
  \item{OC_England}{Ovarian cancer incidence rate for England at the given age.}
  \item{BC_Finland}{Breast cancer incidence rate for Finland at the given age.}
  \item{OC_Finland}{Ovarian cancer incidence rate for Finland at the given age.}
  \item{BC_USA}{Breast cancer incidence rate for the USA at the given age.}
  \item{OC_USA}{Ovarian cancer incidence rate for the USA at the given age.}
  \item{BC_Japan}{Breast cancer incidence rate for Japan at the given age.}
  \item{OC_Japan}{Ovarian cancer incidence rate for Japan at the given age.}
}
}



\examples{
  ## Load the incidence data
  data(incidence_data)
  head(incidence_data)
  
 }


\keyword{breast cancer research}
\keyword{ovarian cancer research}

