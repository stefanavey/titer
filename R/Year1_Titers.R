#' Year 1 titers.
#'
#' Antibody titers to 3 strains of influenza in a cohort of young and older adults from Yale during the 2010-2011 flu season.
#' 
#' @docType data
#' @format A data frame with 42 rows and 11 variables:
#' \describe{
#'   \item{YaleID}{subject identifier, unique}
#'   \item{AgeGroup}{age of subject. 20-35 (Young), >= 65 (Older)}
#'   \item{...}{Other columns folow the format <type>_<strain> where <type> is either Day 0 ("d0"), Day 28 ("d28"), or fold change ("fc").}
#' }
#' @references Thakar J, et al. (2015) Aging-dependent alterations in gene expression and a mitochondrial signature of responsiveness to human influenza vaccination. Aging (Albany NY) 7(1):38–52. \url{https://www.ncbi.nlm.nih.gov/pubmed/25596819}
"Year1_Titers"

