if (!requireNamespace("here", quietly = TRUE)) {
  install.packages("here", repos = "https://cloud.r-project.org")
}

library(here)

source(here::here("requirements.R"))
source(here::here("R/functions.R"))
source(here::here("analysis/01_preparing.R"))

######## Modeling scenarios (execution loop) ########
rhats <- c(1:3)
rolling_cases <- c(FALSE, TRUE)

for(rolling in rolling_cases) {
  for(rhat in rhats){
    print(glue::glue("\n\nWorking in scenario rhat: {rhat} and rolling case: {rolling}\n"))
    source(here::here("analysis/02_models.R"))
  }
}
#####################################################

source(here::here("analysis/03_figures_tables.R"))


#------- Optional (see README.md file for instructions) --------#
#
# rmarkdown::render(here::here("analysis/04_outputs.Rmd"),
#                   output_file = "outputs.html",
#                   output_dir = here::here("outputs"),
#                   quiet = TRUE)
#
#---------------------------------------------------------------#





