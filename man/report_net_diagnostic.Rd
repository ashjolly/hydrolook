% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/report_net_diagnostic.R
\name{report_net_diagnostic}
\alias{report_net_diagnostic}
\title{Generate diagnostic report}
\usage{
report_net_diagnostic(output_type = "pdf", PROV_TERR_STATE_LOC = "BC")
}
\arguments{
\item{output_type}{the type of file to be outputted. Currently html and pdf are supported. defaults to pdf}

\item{PROV_TERR_STATE_LOC}{Province to be surveyed. Defaults to BC.}
}
\description{
run this command to render the net_diag report. The reports are then outputted to the report folder
}
\examples{
\dontrun{
report_net_diagnostic(output_type = "pdf", PROV_TERR_STATE_LOC = "PE")
}
}
\seealso{
Other report_generators: \code{\link{report_regional_streamflow}},
  \code{\link{report_station}}
}
