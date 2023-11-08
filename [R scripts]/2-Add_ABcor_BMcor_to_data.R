#### PACKAGE LOADING ####
library(tidyverse)

#### CALCULATION OF AB_cor and BM_cor VARIABLES ####
# Function that adds a column AB_cor and BM_cor to report Abundance and Biomass to 1 m².
Create_T_Brut_Final <- function(data_VDT_input, table_name) {
  
  print(paste("Processing:", table_name))
  
  if (any(is.na(data_VDT_input$Nbr_VDT))) {
    message = paste("The Nbr_VDT column has empty box(es) in the following rows:",
                    paste(which(is.na(data_VDT_input$Nbr_VDT)), collapse = ", "),
                    ".\n\n Please fill in the relevant blanks."
    )
    stop(message)
  }
  
  data_VDT_input <- data_VDT_input %>%
    rowwise() %>%
    mutate(
      AB_cor = case_when(
        Code_Methode %in% c("F", "M", "FHS") ~ Nbr_VDT,
        Code_Methode %in% c("HS_4","HSAITC_4") ~ 4 * Nbr_VDT,
        Code_Methode == "HSAITC_7.95775385" ~ 7.95775385 * Nbr_VDT,
        Code_Methode == "AITC_11.1111111111111" ~ 11.1111111111111 * Nbr_VDT,
        Code_Methode == "TB" ~ 25 * Nbr_VDT,
        Code_Methode == "HSAITC_6.25" ~ 6.25 * Nbr_VDT,
        Code_Methode %in% c("AITCTM", "HSAITC_16","HS_16") ~ 16 * Nbr_VDT,
        Code_Methode == "TM" & (is.null(Cadre) | is.na(Cadre)) ~ 16 * Nbr_VDT,
        Code_Methode == "TM" & Cadre %in% c("IR", "R") ~ 8 * Nbr_VDT,
        TRUE ~ as.numeric(NA)
      )
    )
  
  
  if (!is.null(data_VDT_input$Pds)) {
    data_VDT_input <- data_VDT_input %>%
      rowwise() %>%
      mutate(
        BM_cor = case_when(
          Code_Methode %in% c("F", "M", "FHS") ~ Pds,
          Code_Methode %in% c("HS_4","HSAITC_4") ~ 4 * Pds,
          Code_Methode == "HSAITC_7.95775385" ~ 7.95775385 * Pds,
          Code_Methode == "TB" ~ 25 * Pds,
          Code_Methode == "HSAITC_6.25" ~ 6.25 * Pds,
          Code_Methode %in% c("AITCTM", "HSAITC_16","HS_16") ~ 16 * Pds,
          Code_Methode == "TM" & (is.null(Cadre) | is.na(Cadre)) ~ 16 * Pds,
          Code_Methode == "TM" & Cadre %in% c("IR", "R") ~ 8 * Pds,
          TRUE ~ as.numeric(NA)
        )
      )
  }
  
  if (is.null(data_VDT_input$Pds)) {
    message("The \"AB_cor\" variable has been created.")
  } else {
    message("The variables \"AB_cor\" and \"BM_cor\" have been created.")
  }
  
  return(data_VDT_input)
}

# Get the names of data frames in the global environment that end with "_rd"
data_frames_names <- ls(pattern = "_rd$", envir = .GlobalEnv)

# Loop through each data frame and apply the function
for (table_name in data_frames_names) {
  modified_data_frame <- Create_T_Brut_Final(get(table_name), table_name)  # Apply the function to the data frame
  
  # Check if the result is not empty before assigning it back to the global environment
  if (!is.null(modified_data_frame) && nrow(modified_data_frame) > 0) {
    assign(table_name, modified_data_frame, envir = .GlobalEnv)  # Replace the modified data frame in the global environment
  } 
}
