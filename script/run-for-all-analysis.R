################################################################################
##                                                                            ##
## Master Script to Generate All Results                                      ##
## Note: This script sequentially calls all R scripts for the project,        ##
##       executing the analyses and producing the final outputs.              ##
################################################################################


# starting time
message(
  paste0("Running the master file on ", Sys.time(), "\n")
)

# libraries
library(here)


# sourcing functions
message("Sourcing functions ...")
source(here("script", "00_functions.R"))

# sourcing R scripts
source(here("script", "01_data_pre_processing.R"))
source(here("script", "02_iptw_adjustment_for_nonresp.R")) 
source(here("script", "03_descriptives.R"))
source(here("script", "04_apc_bounding_apci_analysis.R"))
source(here("script", "05_apc_bounding_apci_analysis_logit.R"))



# print time and session info
message(
  paste0("\n done on ", Sys.time())
)

### END OF CODE ###