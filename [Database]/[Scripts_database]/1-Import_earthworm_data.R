#### PACKAGE LOADING ####
library(tidyverse)

#### FILE IMPORT #### 
# Specify the directory path containing the .csv files
directory_path <- "[Database]/[Raw datasets]/EW_datasets"

# List the .csv files in the specified directory
# Files ending by "_rd" are raw data sets.
csv_files <- list.files(directory_path, pattern = "*_rd.csv", full.names = TRUE)

# Load the .csv files into the environment
for (csv_file in csv_files) {
  table_name <- tools::file_path_sans_ext(basename(csv_file))  # Table name based on the file name
  
  # Read the .csv file with manually set colClasses
  assign(
    table_name,
    read.csv(
      csv_file,
      stringsAsFactors = FALSE,
      colClasses = "character",  # Specify all columns as character initially
      na.strings = c("", "NA", "N/A", "null"),  # Specify NA values that can be in the data
      dec = ",", # Specify that the decimal separator is a comma
      row.names = NULL  # Avoid reading row names as a column
    )
  )
}

#### DATA HOMOGENISATION ####
# useful for processing data and aggregating tables afterwards
# not all columns are necessarily useful and can be deleted later
# list of desired column names :
desired_headers <- c(
  "Programme", "Protocole","Code_Methode" ,"Annee", "Date_Prelevement", "ID_Site",
  "Code_Parcelle", "Site", "Parcelle", "Modalite", "Cadre", "Sous.cadre", "Bloc",
  "Repetition", "CE.Deter.Terrain", "Deter.Code.Taxon", "Code_Taxon", "GF", "Stade",
  "Incomplet", "Nbr_VDT", "Pds", "Pds_GF", "Commentaires"
)

# Function to homogenize column names and data types of a given data frame
homogenize_dataframe <- function(df, desired_headers) {
  # Check if the column names of the data frame match the desired headers
  if (!identical(colnames(df), desired_headers)) {
    # Filter and rearrange existing headers
    existing_headers <- intersect(desired_headers, colnames(df))
    df <- df[, existing_headers]
    
    # Add missing columns with NA values
    missing_headers <- setdiff(desired_headers, colnames(df))
    for (header in missing_headers) {
      df[[header]] <- NA
    }
    
    # Rearrange columns in the desired order
    df <- df[, desired_headers]
  }
  
  # Convert all columns (except for specific ones) to character type
  columns_to_convert <- setdiff(colnames(df), c("Nbr_VDT", "Pds", "Pds_GF"))
  for (col in columns_to_convert) {
    df[[col]] <- as.character(df[[col]])
  }
  
  # Convert specific columns to numeric
  df$Nbr_VDT <- as.numeric(as.character(df$Nbr_VDT))
  df$Pds <- as.numeric(as.character(df$Pds))
  df$Pds_GF <- as.numeric(as.character(df$Pds_GF))
  
  # Translate the levels in the "Protocol" column in English
  if ("Protocole" %in% colnames(df)) {
    df$Protocole <- gsub("MTM", "M_HS", df$Protocole)
    df$Protocole <- gsub("FTM", "F_HS", df$Protocole)
    df$Protocole <- gsub("TB", "HS", df$Protocole)
  }
  
  return(df)
}

# Apply the homogenize_dataframe function to all data frames in the environment
for (df_name in ls()) {
  if (is.data.frame(get(df_name))) {
    assign(df_name, homogenize_dataframe(get(df_name), desired_headers))
  }
}

#### IMPORT DATA ALREADY AGGREGATED (1 line = 1 plot) ####
cp_rd_rep_site=read_csv("./[Database]/[Raw datasets]/EW_datasets/cp_rd_rep_site.csv")


