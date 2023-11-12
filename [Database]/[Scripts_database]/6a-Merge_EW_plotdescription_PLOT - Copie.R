# Load required package
library(dplyr)

# Get list of all objects in the environment
all_objects <- ls()

# Filter tables ending with "_description" and "_rd_rep_site"
meta_tables <- grep("_description$", all_objects, value = TRUE)
rd_rep_site_tables <- grep("_rd_rep_site$", all_objects, value = TRUE)

# Create a list to store the merged tables
fusioned_tables <- list()

# Function to check for the presence of identifier columns in a given table
check_id_columns <- function(data) {
  id_columns <- c("Programme", "ID_Site", "Annee", "Modalite", "Bloc")
  missing_columns <- setdiff(id_columns, colnames(data))
  if (length(missing_columns) > 0) {
    cat("The following identifier columns are missing in the table:\n")
    cat(paste(missing_columns, collapse = ", "), "\n")
    return(FALSE)
  }
  return(TRUE)
}

# Function to perform the join between a "_rd_rep_site" table and a "_description" table
perform_join <- function(meta_table_name, rd_rep_site_table_name) {
  # Load the tables from the environment
  meta_data <- get(meta_table_name)
  rd_rep_site_data <- get(rd_rep_site_table_name)
  
  # Check for the presence of identifier columns in both tables
  if (!check_id_columns(meta_data) || !check_id_columns(rd_rep_site_data)) {
    message(paste("Some identifier columns are missing in", meta_table_name, "or", rd_rep_site_table_name))
    return(NULL)
  }
  
  # Create the identifier as a string while accounting for NAs
  meta_data$ID <- with(meta_data, as.character(paste(Programme, ID_Site, Annee, Modalite, Bloc, sep = "_")))
  rd_rep_site_data$ID <- with(rd_rep_site_data, as.character(paste(Programme, ID_Site, Annee, Modalite, Bloc, sep = "_")))
  rd_rep_site_data <- rd_rep_site_data %>% select(-Programme, -ID_Site, -Annee, -Modalite, -Bloc)
  
  # Perform the join using left_join
  joined_data <- right_join(meta_data, rd_rep_site_data, by = "ID", keep = FALSE)
  
  # Rename the resulting table by appending "_fusion"
  new_table_name <- gsub("_rd_rep_site$", "_fusion_site", rd_rep_site_table_name)
  assign(new_table_name, joined_data, envir = .GlobalEnv)
  
  return(new_table_name)
}


# Iterate through the pairs of "_description" and "_rd_rep_site" tables and perform the joins
for (meta_table_name in meta_tables) {
  rd_rep_site_table_name <- sub("_description$", "_rd_rep_site", meta_table_name)
  if (rd_rep_site_table_name %in% rd_rep_site_tables) {
    new_table_name <- perform_join(meta_table_name, rd_rep_site_table_name)
    if (!is.null(new_table_name)) {
      fusioned_tables[[new_table_name]] <- get(new_table_name)
    }
  }
}

# Perform a vertical join using bind_rows()
LandWorm_dataset_site <- bind_rows(fusioned_tables)

