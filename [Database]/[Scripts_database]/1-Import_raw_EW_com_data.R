#### PACKAGE LOADING ####
library(tidyverse)
library(readxl)

#### FILE IMPORT #### 
# Initialiser un vecteur pour stocker les messages concernant les Code_Taxon sans correspondance
unmatched_messages <- character()

# Load the Taxa_list_fev2021 table
taxa_list <- read_excel("[Database]/[Raw Datasets]/T_TAXONOMIE_v2024.02.09.xlsx") %>%
  dplyr::select(Code_Taxon, Taxon)  # Select only Code_Taxon and Taxon columns


# Specify the directory path containing the .csv files
directory_path <- "[Database]/[Raw datasets]/EW_datasets"

# List the .csv files in the specified directory
csv_files <- list.files(directory_path, pattern = "*_rd.csv", full.names = TRUE)

# Load the .csv files
for (csv_file in csv_files) {
  table_name <- tools::file_path_sans_ext(basename(csv_file))  # Determine the table name based on the file name
  
  # Load the .csv file
  data <- read.csv(
    csv_file,
    stringsAsFactors = FALSE,
    colClasses = "character",  # Treat all columns as character initially
    na.strings = c("", "NA", "N/A", "null"),  # Specify NA representations in the data
    dec = ",",  # Set the decimal separator as comma
    row.names = NULL  # Do not treat any column as row names
  )
  
  # Supprimer les "?" dans Code_Taxon
  data$Code_Taxon <- gsub("\\?", "", data$Code_Taxon)
  
  # Merge with taxa_list based on Code_Taxon
  merged_data <- left_join(data, taxa_list, by = "Code_Taxon")

  # Identifier les Code_Taxon sans correspondance dans les deux colonnes
  no_match <- unique(merged_data$Code_Taxon[is.na(merged_data$Taxon) & !merged_data$Code_Taxon %in% taxa_list$Taxon])
  
  # Vérifier si des Code_Taxon sans correspondance ont été identifiés et les ajouter aux messages
  if (length(no_match) > 0) {
    unmatched_msg <- paste("Dans", table_name, "les Code_Taxon non trouvés sont:", paste(no_match, collapse = ", "))
    unmatched_messages <- c(unmatched_messages, unmatched_msg)
  }
  
  # Assigner le résultat à une variable dans l'environnement
  assign(table_name, merged_data)
}

# Vérifier si des messages sur les Code_Taxon sans correspondance ont été générés
if (length(unmatched_messages) == 0) {
  message("Tous les taxons ont trouvé correspondance.")
} else {
  message("Des problèmes ont été identifiés avec certains Code_Taxon :")
  for (msg in unmatched_messages) {
    message(msg)
  }
}


#### DATA HOMOGENISATION ####
# useful for processing data and aggregating tables afterwards
# not all columns are necessarily useful and can be deleted later
# list of desired column names :
desired_headers <- c(
  "Programme", "Protocole","Code_Methode" ,"Annee", "Date_Prelevement", "ID_Site",
  "Code_Parcelle", "Site", "Parcelle", "Modalite", "Cadre", "Sous.cadre", "Bloc",
  "Repetition","Code_Taxon","Taxon", "Stade",
  "Incomplet", "Nbr_VDT", "Pds", "Commentaires"
)

# Function to homogenize column names and data types of a given data frame
homogenize_dataframe <- function(df, desired_headers) {
  # Gérez les dataframes vides
  if (ncol(df) == 0) {
    df <- data.frame(matrix(ncol = length(desired_headers), nrow = 0))
    colnames(df) <- desired_headers
  } else {
    # Assurez-vous que les colonnes désirées existent dans df
    existing_headers <- intersect(desired_headers, colnames(df))
    df <- df[, existing_headers, drop = FALSE]
    
    # Ajoutez les colonnes manquantes
    missing_headers <- setdiff(desired_headers, colnames(df))
    df[missing_headers] <- NA
    
    # Réorganisez les colonnes
    df <- df[, desired_headers, drop = FALSE]
  }
  
  # Convert all columns (except for specific ones) to character type
  columns_to_convert <- setdiff(colnames(df), c("Nbr_VDT", "Pds"))
  for (col in columns_to_convert) {
    df[[col]] <- as.character(df[[col]])
  }
  
  # Convert specific columns to numeric
  df$Nbr_VDT <- as.numeric(as.character(df$Nbr_VDT))
  df$Pds <- as.numeric(as.character(df$Pds))

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
INPN_rd_rep_site=read_csv("./[Database]/[Raw datasets]/EW_datasets/INPN_rd_rep_site.csv")


