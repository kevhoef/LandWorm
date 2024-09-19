#### PACKAGE LOADING ####
library(tidyverse)
library(openxlsx)
library(readxl)

# Function to homogenize dates into common format (JJ_MM_AAAA, MM_AAAA, or AAAA)
homogenize_dates <- function(dates) {
  # Convert character to Date if possible
  parsed_dates <- suppressWarnings(parse_date_time(dates, orders = c("dmy", "my", "y", "ym", "ymd", "ym")))
  
  # Handle failed parsing by replacing NA dates with original dates for manual inspection
  failed_parsing <- is.na(parsed_dates)
  if (any(failed_parsing)) {
    warning("Some dates could not be parsed: ", paste(dates[failed_parsing], collapse = ", "))
    parsed_dates[failed_parsing] <- dates[failed_parsing]
  }
  
  # Format dates into desired format (JJ_MM_AAAA, MM_AAAA, AAAA)
  formatted_dates <- ifelse(!is.na(parsed_dates), format(parsed_dates, "%d_%m_%Y"), NA)
  
  # Handle cases where the original format is MM_AAAA or AAAA
  formatted_dates <- ifelse(
    str_detect(dates, "^\\d{2}_\\d{4}$"), dates, formatted_dates
  )
  formatted_dates <- ifelse(
    str_detect(dates, "^\\d{4}$"), dates, formatted_dates
  )
  
  return(formatted_dates)
}



#### FILE IMPORT #### 
# Load the Taxa_list_fev2021 table
taxa_list <- read_excel("[Database]/[Raw datasets]/Table_Taxonomy_LandWorm_07_03_2024.xlsx") %>%
  dplyr::select(Code_Taxon, Taxon)  # Select only Code_Taxon and Taxon columns

# Spécifier le chemin du répertoire contenant les fichiers .xlsx
directory_path <- "[Database]/[Raw datasets]/EW_datasets"

# Lister les fichiers .xlsx dans le répertoire spécifié
xlsx_files <- list.files(directory_path, pattern = "_rd\\.xlsx$", full.names = TRUE)

#### FILE IMPORT #### 
# Initialiser un vecteur pour stocker les messages relatifs aux Tax_Codes sans correspondance
unmatched_messages <- character()

# Charger les fichiers .xlsx
for (xlsx_file in xlsx_files) {
  table_name <- tools::file_path_sans_ext(basename(xlsx_file))  # Déterminer le nom de la table basé sur le nom du fichier
  
  # Extraire le propriétaire du jeu de données à partir du nom du fichier
  owner <- sub("^(.*?)_.*", "\\1", table_name)
  
  # Charger le fichier .xlsx
  data <- read_excel(
    xlsx_file,
    na = c("", "NA", "N/A", "null"),  # Spécifier les représentations NA dans les données
    col_types = "text"  # Lire toutes les colonnes comme des caractères
  ) %>%
    mutate(across(everything(), as.character)) %>%  # Assurer que toutes les colonnes sont bien des caractères
    mutate(owner = owner)  # Ajouter la colonne owner avec la sous-chaîne extraite

  # Check if the Code_Taxon column exists (Mainly for Daniel Cluzeau's data)
  if ("Code_Taxon" %in% colnames(data)) {
    # Remove "?" in Code_Taxon (It may be debatable, but if the taxon code is positioned, then there is evidence)
    data$Code_Taxon <- gsub("\\?", "", data$Code_Taxon)
    
    # Create a named vector with Code_Taxon as names and Taxon as values
    taxa_vector <- setNames(taxa_list$Taxon, taxa_list$Code_Taxon)
    
    # Update or create the Taxon column based on Code_Taxon
    data$Taxon <- taxa_vector[data$Code_Taxon]
    
    # Identify unmatched Code_Taxon, excluding NAs
    no_match_code_taxon <- data$Code_Taxon[is.na(data$Taxon)]
    
    # Generate messages based on whether there are unmatched Code_Taxon values
    if (length(unique(no_match_code_taxon)) > 0) {
      unmatched_msg_code_taxon <- paste("In", table_name, "unmatched Code_Taxon are:", paste(unique(no_match_code_taxon), collapse = ", "))
      unmatched_messages <- c(unmatched_messages, unmatched_msg_code_taxon)
    } else {
      # Generate a message indicating all Code_Taxon values have a corresponding Taxon in taxa_list
      all_matched_code_taxon_msg <- paste("All Code_Taxon values in", table_name, "have a corresponding Taxon in taxa_list.")
      unmatched_messages <- c(unmatched_messages, all_matched_code_taxon_msg)
    }
  }
  # Check if the Taxon column exists
  if ("Taxon" %in% colnames(data)) {
    # Verify if each Taxon in the dataframe exists in taxa_list
    taxon_exist <- data$Taxon %in% taxa_list$Taxon
    
    # Identify Taxon values that do not exist in taxa_list
    no_match_taxon <- data$Taxon[!taxon_exist]
    
    # Generate messages based on whether there are unmatched Taxon values
    if (length(unique(no_match_taxon)) > 0) {
      unmatched_msg_taxon <- paste("In", table_name, "unmatched Taxon values are:", paste(unique(no_match_taxon), collapse = ", "))
      unmatched_messages <- c(unmatched_messages, unmatched_msg_taxon)
    } else {
      all_matched_taxon_msg <- paste("All Taxon values in", table_name, "have a corresponding entry in taxa_list.")
      unmatched_messages <- c(unmatched_messages, all_matched_taxon_msg)
    }
  }
  # Assign the result to a variable in the environment
  assign(table_name, data, envir = .GlobalEnv)
}

# Print the unmatched messages after the loop
if (length(unmatched_messages) > 0) {
  cat("Messages:\n")
  for (message in unmatched_messages) {
    cat(message, "\n")  # Print each message on a new line
  }
} else {
  cat("All Code_Taxon and Taxon values in all files have corresponding entries in taxa_list.\n")
}




#### DATA HOMOGENISATION ####
# useful for processing data and aggregating tables afterwards
# not all columns are necessarily useful and can be deleted later
# list of desired column names :
desired_headers <- c(
  "owner","Programme", "Protocole","Code_Methode" ,"Annee", "Date_Prelevement", "ID_Site",
  "Code_Parcelle", "Site", "Parcelle", "Modalite", "Cadre", "Sous.cadre", "Bloc",
  "Repetition","Code_Taxon","Taxon", "Stade",
  "Incomplet", "Nbr_VDT", "Pds", "Commentaires"
)

# Function to homogenize column names and data types of a given data frame
homogenize_dataframe <- function(df, desired_headers, df_name) {
  cat("Processing dataframe:", df_name, "\n")
  
  # Manage empty dataframes
  if (ncol(df) == 0) {
    df <- data.frame(matrix(ncol = length(desired_headers), nrow = 0))
    colnames(df) <- desired_headers
  } else {
    # Make sure that the columns exist in df
    existing_headers <- intersect(desired_headers, colnames(df))
    df <- df[, existing_headers, drop = FALSE]
    
    # Add missing columns
    missing_headers <- setdiff(desired_headers, colnames(df))
    df[missing_headers] <- NA
    
    # Rearrange columns
    df <- df[, desired_headers, drop = FALSE]
  }
  
  # Convert all columns (except for specific ones) to character type
  columns_to_convert <- setdiff(colnames(df), c("Nbr_VDT", "Pds"))
  for (col in columns_to_convert) {
    df[[col]] <- as.character(df[[col]])
  }
  
  # Remove values containing "?" in Nbr_VDT and Pds
  df$Nbr_VDT <- ifelse(grepl("\\?", df$Nbr_VDT), NA, df$Nbr_VDT)
  df$Pds <- ifelse(grepl("\\?", df$Pds), NA, df$Pds)
  
  # Replace commas with dots in Nbr_VDT and Pds
  df$Nbr_VDT <- gsub(",", ".", df$Nbr_VDT)
  df$Pds <- gsub(",", ".", df$Pds)
  
  # Convert Nbr_VDT and Pds to numeric
  df$Nbr_VDT <- as.numeric(df$Nbr_VDT)
  df$Pds <- as.numeric(df$Pds)
  
  # Translate the levels in the "Protocole" column in English
  if ("Protocole" %in% colnames(df)) {
    df$Protocole <- gsub("MTM", "M_HS", df$Protocole)
    df$Protocole <- gsub("FTM", "F_HS", df$Protocole)
    df$Protocole <- gsub("TB", "HS", df$Protocole)
  }
  
  return(df)
}


# Apply the homogenize_dataframe function to all data frames in the environment
df_names <- ls()
df_names_to_process <- df_names[grep("_rd$", df_names)] ## Filtrer les noms de dataframe qui se terminent par "_rd"

# Appliquer la boucle uniquement aux dataframes sélectionnés
for (df_name in df_names_to_process) {
  print(df_name)
  tryCatch({
    df <- get(df_name)
    df <- homogenize_dataframe(df, desired_headers, df_name)  # Passer df_name comme argument supplémentaire
    assign(df_name, df, envir = .GlobalEnv)
  }, error = function(e) {
    print(paste("Error in dataframe:", df_name))
    print(e)
  }, warning = function(w) {
    print(paste("Warning in dataframe:", df_name))
    print(w)
  })
}

#### IMPORT DATA ALREADY AGGREGATED (1 line = 1 plot) ####
cp_rd_rep_site <- read_excel("./[Database]/[Raw datasets]/EW_datasets/cp_rd_rep_site.xlsx")
INPN_rd_rep_site <- read_excel("./[Database]/[Raw datasets]/EW_datasets/INPN_rd_rep_site.xlsx")









##################
##################


##################
##################

# Get the names of dataframes ending with _rd
df_names <- ls(pattern = "_rd$")
df_names_to_process <- df_names[grep("_rd$", df_names)]

# Initialize a list to store results
results <- list()

# Loop over each dataframe and check the condition
for (df_name in df_names_to_process) {
  # Get the dataframe object
  df <- get(df_name)
  
  # Check if Pds > 15 when Nbr_VDT = 1
  check <- df %>%
    filter(Nbr_VDT == 1 & Pds > 15)
  
  # Store the result
  results[[df_name]] <- check
}

# Print the results
for (df_name in names(results)) {
  cat("Results for", df_name, ":\n")
  print(results[[df_name]])
  cat("\n")
}

