library(dplyr)

# Function to calculate mean by "Site" or "Plot" based on an identifier 
mean_by_identifier <- function(data) {
  # List of identification columns
  id_cols <- c("Programme", "Protocole", "Code_Methode", "ID_Site","Code_Parcelle", "Annee", "Modalite", "Bloc","Date_prelevement")
  
  # Filter the identification columns that are present in the data
  present_id_cols <- id_cols[id_cols %in% colnames(data)]
  
  # Check if any identification columns are present
  if (length(present_id_cols) > 0) {
    # Calculate the mean for numeric columns, ignoring NA values, and group by identifier columns
    data %>%
      dplyr::filter(if_any(where(is.numeric), ~!is.na(.))) %>%
      dplyr::group_by(across(all_of(present_id_cols))) %>%
      dplyr::summarise(across(where(is.numeric), ~mean(., na.rm = TRUE)),
                       across(where(is.factor), ~first(.)),
                       .groups = "drop") %>%
      dplyr::select(-Repetition) -> new_data
    
    # Create the new table name with "_site"
    new_table_name <- paste(deparse(substitute(data)), "_site", sep = "")
    
    # Assign the new data frame to the parent frame environment
    assign(new_table_name, new_data, envir = parent.frame())
    
    # Return the new data frame
    return(new_data)
  } else {
    # If no identification columns are present, stop execution and throw an error
    stop("No identification columns found in the data.")
  }
}

# Get the names of data frames in the global environment that end with "_rd"
data_frames_names <- ls(pattern = "_rd_rep$", envir = .GlobalEnv)

# Loop through each data frame to apply the function
for (table_name in data_frames_names) {
  # Define the name for the new data frame with "_site" appended
  new_table_name <- paste0(table_name, "_site")
  
  # Use tryCatch to handle any errors that might occur during the function call
  tryCatch({
    # Print the name of the current table being processed
    cat("Processing table:", table_name, "\n")
    
    # Call the mean_by_identifier function
    output_table <- mean_by_identifier(get(table_name))
    
    # Check if the function returns a valid output
    if (!is.null(output_table)) {
      cat("Averaged by plot calculated for the table", table_name, "\n")
      # Assign the output table to the global environment with the new name
      assign(new_table_name, output_table, envir = .GlobalEnv)
    } else {
      cat("Error in creating the table", new_table_name, "\n")
    }
  }, error = function(e) {
    # If an error occurs, print the error message
    cat("Error while processing the table", table_name, ":", conditionMessage(e), "\n")
  })
}
