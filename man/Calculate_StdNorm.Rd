% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Calculate_padjMFC.R
\name{Calculate_StdNorm}
\alias{Calculate_StdNorm}
\title{Calculate Normalized Titers}
\usage{
Calculate_StdNorm(
  dat,
  type,
  fcToOne = FALSE,
  idCol = "SubjectID",
  cols = grep(paste0(type, "_[AB]"), colnames(dat), value = TRUE)
)
}
\arguments{
\item{dat}{Data frame containing \code{fcStdCols}}

\item{type}{What should be standarized. Either "d0", or "fc".}

\item{fcToOne}{Logical. Are titer fold changes allowed to be less than 1
or should these be changed to 1 before standardization?
Default is FALSE and no changes will be made. Only relevant
when \code{type == "fc"}}

\item{idCol}{Name of column containing subject IDs}

\item{cols}{column names containing the titer measurements
for each strain. Length must be > 1}
}
\value{
A data frame like \code{dat} but with standarized columns added
}
\description{
\code{Calculate_StdNorm} calculates the standardized d0 or fc titers
}
\details{
This must be run on only 1 cohort at a time because titers will be normalized
across all subjects. The median is used but unlike the original reference,
the standard deviation is calculated rather than the maximum absolute deviation.
}
\examples{
## First Example

}
\references{
Tsang JS, et al. (2014) Global analyses of human immune variation reveal baseline predictors of postvaccination responses. Cell 157(2):499-513.
}
\author{
Stefan Avey
}
