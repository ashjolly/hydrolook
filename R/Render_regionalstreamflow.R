# clear workspace
rm(list = ls())

drive = "\\\\DRAIN.dmz\\Shared"
drive_G = "\\\\Backhoe\\s63101\\Watershare\\rfc"
drive_Q = "\\\\wwwt.env.gov.bc.ca\\envwwwt\\rfc"

source("C:/Users/AJOLLYMO/RProjects/hydrolook/R/droughtstats_function.R")

year <- year(Sys.Date())
start.time1 <- Sys.time()

devtools::install_github('rstudio/rmarkdown')

# Run the Rmarkdown for the daily statistics files
#rmarkdown::render(paste0("C:/Users/AJOLLYMO/RProjects/hydrolook/inst/templates/regional_streamflow.Rmd"),
#                  output_file = paste0("C:/Users/AJOLLYMO/RProjects/hydrolook/inst/templates/WestCoastStats_",
#                                       Sys.Date(), ".html"),
#                  params = list(region = c("West Coast Natural Resource Region")))

natural_resource_regions = unique(nr_regions()$ORG_UNIT_NAME)

for(i in 1:length(natural_resource_regions)){
  rmarkdown::render(paste0("C:/Users/AJOLLYMO/RProjects/hydrolook/inst/templates/regional_streamflow.Rmd"),
                    output_file = paste0("C:/Users/AJOLLYMO/RProjects/hydrolook/inst/templates/",
                                         natural_resource_regions[i], "_", Sys.Date(), ".html"),
                    params = list(region = natural_resource_regions[i]))
}
