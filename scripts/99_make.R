# PARAMETERS
# start date of REFNET
# NB this parameter is set to define the time period during which the temporal trends are assessed
refnet_start <- 2013

# if whole continental France, set selected_depts <- NULL
# selected_depts <- c(22, 29, 35, 56)
selected_depts <- NULL

# PACKAGES
# Install necessary libraries ----
## CRAN packages
packages_list <- c("rmarkdown",
                   "remotes",
                   "tidyverse",
                   "flextable",
                   "VennDiagram",
                   "trend",
                   "sf")

installed_packages <- packages_list %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages_list[!installed_packages])
}

## GitHub packages if not already installed
# remotes::install_github("pascalirz/aspe")
# remotes::install_github("ropensci/rnaturalearthhires")
# remotes::install_github("pascalirz/tod")


# make background map ----
source(file = "scripts/10_make_background_maps.R")


# select data ----
rmarkdown::render(
  input = 'scripts/20_select_data.Rmd',
  params = list(refnet_start = refnet_start,
                depts = selected_depts)
)

# calculate all indicators at the survey ----
rmarkdown::render(
  input = 'scripts/30_calculate_fbi_and_pop_indicators.Rmd',
  params = list(refnet_start = refnet_start)
)


# calculate temporal trends ----
rmarkdown::render(
  input = 'scripts/40_assess_temporal_trends.Rmd'
)

# assemble results ----
rmarkdown::render(
  input = 'scripts/70_output_main.Rmd',
  params = list(refnet_start = refnet_start)
)

# export figures in HD ----
rmarkdown::render(
  input = 'scripts/90_export_figures.Rmd',
  output_file = 'scripts/90_export_figures.html'
)


# export figures in HD ----
rmarkdown::render(
  input = 'scripts/95_appendix.Rmd',
  output_file = '../manuscript/fbi_trends_supplementary_material.html'
)
