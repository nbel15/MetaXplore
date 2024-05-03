
shinyServer(function(input, output, session) {
  
  #######################################################
  #Analysis pipeline function
  #######################################################
  
  source("server_functions/Import_analysis_files.R", local=TRUE)
  source("server_functions/Alpha_div_functions.R", local=TRUE)
  source("server_functions/Beta_div_functions.R", local=TRUE)
  source("server_functions/R_abundance_functions.R", local=TRUE)
  source("server_functions/Diff_abund_functions.R", local=TRUE)
  source("server_functions/core_functions.R", local=TRUE)
  source("server_functions/Biomarker_functions.R", local=TRUE)
  
  
  #######################################################
  #Guide
  #######################################################
  source("server_functions/Guide_functions.R", local=TRUE)
  
})