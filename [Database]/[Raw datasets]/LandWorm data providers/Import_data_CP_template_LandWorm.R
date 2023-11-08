# Load the necessary libraries
library(readxl)
library(tidyverse)

# Specify the path to your Excel file
excel_file_path <- "[Database]/Raw datasets/cp_data_site.xlsx"

# Read data from each sheet
sheet_data_Study <- read_excel(excel_file_path, sheet = "5. data Study") %>%
  rename(GPS_Y = gps_x,GPS_X = gps_y) %>%
  mutate(ID_bloc = na_if(ID_bloc, "na")) %>%
  mutate(ID_modality = na_if(ID_modality, "na"))%>%
  mutate_all(~gsub(" ", "_", .))
sheet_data_Land_use_management <- read_excel(excel_file_path, sheet = "7.data Land_use_management")%>%
  mutate(ID_bloc = na_if(ID_bloc, "na")) %>%
  mutate(ID_modality = na_if(ID_modality, "na"))%>%
  mutate_all(~gsub(" ", "_", .))
sheet_data_Soil <- read_excel(excel_file_path, sheet = "9.data Soil")%>%
  mutate(ID_bloc = na_if(ID_bloc, "na")) %>%
  mutate(ID_modality = na_if(ID_modality, "na"))%>%
  mutate_all(~gsub(" ", "_", .))
sheet_data_Sampling <- read_excel(excel_file_path, sheet = "11. data Sampling")%>%
  mutate(ID_bloc = na_if(ID_bloc, "na")) %>%
  mutate(ID_modality = na_if(ID_modality, "na"))%>%
  mutate_all(~gsub(" ", "_", .))
cp_rd_rep_site <- read_excel(excel_file_path, sheet = "13. Data Earth. com.") %>%
  mutate(ID_bloc = na_if(ID_bloc, "na")) %>%
  mutate(ID_modality = na_if(ID_modality, "na"))

##### EARTHWORM COMMUNITIES
# Extract and replace the terms in EW_species and update stage_dev
cp_rd_rep_site <- cp_rd_rep_site %>%
  mutate(stage_dev = case_when(
    grepl("sub-adulte$", EW_species) ~ "sub-adulte",
    grepl("adulte$", EW_species) ~ "adulte",
    grepl("juvenile$", EW_species) ~ "juvenile",
    grepl("juveniles$", EW_species) ~ "juveniles",
    TRUE ~ NA_character_  # Set to NA if none of the conditions match
  )) %>%
  mutate(EW_species = sub("(adulte|sub-adulte|juvenile|juveniles)$", "", EW_species)) %>%
  mutate(EW_species = as.factor(EW_species))
  

levels(cp_rd_rep_site$EW_species)
cp_rd_rep_site <- cp_rd_rep_site %>%
  mutate(EW_species = case_when(
    EW_species == "A. chlorotica " ~ "AC",
    EW_species == "A. cupulifera " ~ "ACU",
    EW_species == "A. giardi " ~ "NG",
    EW_species == "A. longa " ~ "NLX",
    EW_species == "A. longa/giardi " ~ "NLX/NG",
    EW_species == "A. muldali " ~ "AM",
    EW_species == "A. muldali/rosea " ~ "AM/ARR",
    EW_species == "A. rosea " ~ "ARR",
    EW_species == "Allolobophora chlorotica" ~ "ACX",
    EW_species == "Allolobophora muldali" ~ "AM", 
    EW_species == "Anécique lumbricus " ~ "LX_A",
    EW_species == "Aporrectodea caliginosa" ~ "NCX",
    EW_species == "Aporrectodea giardi" ~ "NG",
    EW_species == "Aporrectodea icterica" ~ "AI",
    EW_species == "Aporrectodea longa" ~ "NLX",
    EW_species == "Aporrectodea rosea" ~ "ARR",
    EW_species == "D. mammalis " ~ "DM",
    EW_species == "D. pygmea " ~ "DPC",
    EW_species == "D. subrubicunda " ~ "DSS",
    EW_species == "E. tetraedra " ~ "ET",
    EW_species == "Eisenia " ~ "EX",
    EW_species == "endoge indetermine " ~ "END",
    EW_species == "L. castaneus " ~ "LC",
    EW_species == "L. centralis " ~ "LCE",
    EW_species == "L. friendi " ~ "LFR", 
    EW_species == "L. friendi/centralis " ~ "LFR/LCE",
    EW_species == "L. rubellus friendoides " ~ "LRFR",
    EW_species == "L. rubellus rubellus " ~ "LRR",
    EW_species == "Lumbricus castaneus" ~ "LC",
    EW_species == "Lumbricus festivus" ~ "LFE",
    EW_species == "Lumbricus terrestris" ~ "LT",
    EW_species == "M. dubius " ~ "MD",
    EW_species == "N. caliginosa " ~ "NCX",
    EW_species == "N. caliginosus meridionalis " ~ "NCM",
    EW_species == "N. terrestris terrestris " ~ "NG",
    EW_species == "non identifie epige" ~ "EPI",
    EW_species == "O. cyaneum " ~ "OC",
    EW_species == "O. lacteum " ~ "OL",
    EW_species == "Octolasium " ~ "OX",
    EW_species == "Octolasium cyaneum" ~ "OC",
    EW_species == "P. fragilis " ~ "PF",
    TRUE ~ EW_species  # Si aucune condition n'est satisfaite, conservez la valeur d'origine
  )) %>%
  mutate(EW_species = as.factor(EW_species))

cp_rd_rep_site <- cp_rd_rep_site %>%
  group_by(
    project_acronym,
    sampling_method_acronym,
    ID_site,
    ID_plot,
    date_plot_collection,
    ID_modality,
    ID_bloc,
    EW_species
  ) %>%
  summarize(
    AB = sum(abundance, na.rm = TRUE),
    BM = sum(biomass, na.rm = TRUE)
  ) %>%
  ungroup()

cp_rd_rep_site <- cp_rd_rep_site %>%
  pivot_wider(
    id_cols = c(
      project_acronym,
      sampling_method_acronym,
      ID_site,
      ID_plot,
      date_plot_collection,
      ID_modality,
      ID_bloc
    ),
    names_from = EW_species,
    values_from = c(AB, BM),
    values_fill = NA,
    names_sep = "_"
  )


cp_rd_rep_site <- cp_rd_rep_site %>%
  mutate(AB_tot = rowSums(select(., starts_with("AB_")), na.rm = TRUE)) %>%
  mutate(BM_tot = rowSums(select(., starts_with("BM_")), na.rm = TRUE)) %>%
  mutate(Richesse = rowSums(select(., starts_with("AB_")) > 0, na.rm = TRUE))



#### RENOMMER LES VARIABLES
cp_rd_rep_site <- cp_rd_rep_site %>%
  rename(
    Programme = project_acronym,
    Protocole = sampling_method_acronym,
    Parcelle = ID_plot,
    Annee= date_plot_collection,
    Modalite = ID_modality,
    Bloc = ID_bloc,
    ID_Site = ID_site
    
    # Ajoutez autant de renommages que nécessaire
  )


#### HOMOGENEISATION DES COLONNES ####
# List of desired column names
desired_headers <- c(
  "Programme", "Protocole", "Code_Methode", "Annee", "Date_Prelevement", "ID_Site",
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
  
  return(df)
}

cp_rd_rep_site=homogenize_dataframe(cp_rd_rep_site, desired_headers)
  
cp_rd_rep_site <- cp_rd_rep_site %>%
  mutate_at(vars(Programme, ID_Site, Annee, Modalite, Bloc), as.character)

write.csv(cp_rd_rep_site, "[Database]/Raw datasets/cp_rd_rep_site.csv", row.names = FALSE)










##### META

# Créez l'identifiant commun dans chaque tableau en utilisant paste() pour le concaténer avec "_"
sheet_data_Study <- sheet_data_Study %>%
  mutate(ID_common = as.character(paste(ID_site, ID_plot, date_plot_collection, ID_modality, sep = "_")))

sheet_data_Land_use_management <- sheet_data_Land_use_management %>%
  mutate(ID_common = as.character(paste(ID_site, ID_plot, date_plot_collection, ID_modality, sep = "_")))

sheet_data_Soil <- sheet_data_Soil %>%
  mutate(ID_common = as.character(paste(ID_site, ID_plot, date_plot_collection, ID_modality, sep = "_")))

sheet_data_Sampling <- sheet_data_Sampling %>%
  mutate(ID_common = as.character(paste(ID_site, ID_plot, date_plot_collection, ID_modality, sep = "_")))


# Supprimez les colonnes individuelles d'identifiants dans les tableaux supplémentaires
sheet_data_Land_use_management <- sheet_data_Land_use_management %>%
  select(-c("project_acronym","sampling_method_acronym","ID_site", "ID_plot", "date_plot_collection", "ID_modality","ID_bloc"))

sheet_data_Soil <- sheet_data_Soil %>%
  select(-c("project_acronym","sampling_method_acronym","ID_site", "ID_plot", "date_plot_collection", "ID_modality","ID_bloc"))

sheet_data_Sampling <- sheet_data_Sampling %>%
  select(-c("project_acronym","sampling_method_acronym","ID_site", "ID_plot", "date_plot_collection", "ID_modality","ID_bloc"))


# Effectuez les jointures en utilisant l'identifiant commun
cp_meta <- sheet_data_Study %>%
  left_join(sheet_data_Land_use_management, by = "ID_common") %>%
  left_join(sheet_data_Soil, by = "ID_common") %>%
  left_join(sheet_data_Sampling, by = "ID_common")


# Supprimez la colonne ID_common si vous n'en avez plus besoin
cp_meta <- cp_meta %>%
  select(-ID_common)


# Récupérez les noms des colonnes
noms_colonnes <- colnames(cp_meta)

# Affichez les noms des colonnes
print(noms_colonnes)


cp_meta <- cp_meta %>%
  rename(
    Programme = project_acronym,
    Protocole = sampling_method_acronym,
    Parcelle = ID_plot,
    Annee= date_plot_collection,
    Modalite = ID_modality,
    Bloc = ID_bloc,
    ID_Site=ID_site

    # Ajoutez autant de renommages que nécessaire
  )

cp_meta <- cp_meta %>%
  mutate_at(vars(Programme, ID_Site, Annee, Modalite, Bloc), as.character)


write.csv(cp_meta, "[Database]/Raw datasets/descriptive_data/cp_meta.csv", row.names = FALSE)
