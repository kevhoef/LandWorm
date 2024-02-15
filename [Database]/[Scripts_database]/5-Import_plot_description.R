library(tidyverse)
library(readxl)


print_duplicate_columns <- function(data) {
  duplicated_cols <- names(data)[duplicated(names(data))]
  if (length(duplicated_cols) > 0) {
    message("Colonnes en double dans '", table_name, "': ", paste(duplicated_cols, collapse = ", "))
  }
}
# Chemin vers le fichier Excel de correspondance
correspondance_path <- "[Database]/[Raw datasets]/Correspondance_CLC_LandWorm.xlsx"# Charger les feuilles de correspondance pour clcm_lvl2 et clcm_lvl3
correspondance_colonnes <- read_excel(correspondance_path, sheet = "Colonnes")
correspondance_clcm_lvl1 <- read_excel(correspondance_path, sheet = "clcm_lvl1")
correspondance_clcm_lvl2 <- read_excel(correspondance_path, sheet = "clcm_lvl2")
correspondance_clcm_lvl3 <- read_excel(correspondance_path, sheet = "clcm_lvl3", col_types ="text")
correspondance_clcm_lvl4 <- read_excel(correspondance_path, sheet = "clcm_lvl4")

# Specify the directory path containing the .csv files
directory_path <- "[Database]/[Raw datasets]/Plot_description_datasets"

# List the .csv files in the specified directory
csv_files <- list.files(directory_path, pattern = "*_description.csv", full.names = TRUE)



for (csv_file in csv_files) {
  table_name <- tools::file_path_sans_ext(basename(csv_file))
  
  # Charger le fichier CSV
  data <- read.csv(
    csv_file,
    sep = ",",
    colClasses = "character",
    stringsAsFactors = FALSE,
    na.strings = c("", "NA", "N/A", "null")
  )
  
  # Renommer les colonnes selon la feuille "Colonnes"
  for (i in 1:nrow(correspondance_colonnes)) {
    old_name <- correspondance_colonnes$Ancien[i]
    new_name <- correspondance_colonnes$Nouveau[i]
    if (old_name %in% colnames(data)) {
      colnames(data)[colnames(data) == old_name] <- new_name
    }
  }


  #########
  
  # Assurer que clcm_lvl2 est au format correct et non NA
  if ("clcm_lvl2" %in% colnames(data)) {
    # Filtrer correspondance_clcm_lvl4 pour exclure les lignes où Old_names_lvl3 est NA
    filtered_correspondance <- correspondance_clcm_lvl2 %>% 
      select(clcm_lvl2, code_clcm_lvl2)    # Sélectionner uniquement les colonnes nécessaires
    
    # Jointure avec la table de correspondance filtrée
    data <- left_join(data, filtered_correspondance, by = c("clcm_lvl2" = "clcm_lvl2"))
    
    # S'assurer que les valeurs NA dans clcm_lvl3 n'entraînent pas de changements indésirables
    data$clcm_lvl2 <- ifelse(is.na(data$clcm_lvl2), NA, data$clcm_lvl2)
  }
  
  #######
  
  # Assurer que clcm_lvl3 est au format correct et non NA
  if ("clcm_lvl3" %in% colnames(data)) {
    data$clcm_lvl3 <- gsub(" ", "_", data$clcm_lvl3)
    
    # Filtrer correspondance_clcm_lvl4 pour exclure les lignes où Old_names_lvl3 est NA
    filtered_correspondance <- correspondance_clcm_lvl4 %>% 
      filter(!is.na(Old_names_lvl3))%>%
      select(Old_names_lvl3, clcm_lvl4, code_clcm_lvl4)    # Sélectionner uniquement les colonnes nécessaires
    
    # Jointure avec la table de correspondance filtrée
    data <- left_join(data, filtered_correspondance, by = c("clcm_lvl3" = "Old_names_lvl3"))
    
    # S'assurer que les valeurs NA dans clcm_lvl3 n'entraînent pas de changements indésirables
    data$clcm_lvl4 <- ifelse(is.na(data$clcm_lvl3), NA, data$clcm_lvl4)
  }
  
  # Initialize or check code_clcm_lvl columns
  lvl_cols <- paste0("code_clcm_lvl", 1:4)
  for (lvl_col in lvl_cols) {
    if (!lvl_col %in% colnames(data)) {
      data[[lvl_col]] <- NA_character_
    } else {
      data <- data %>%
        mutate(!!lvl_col := if_else(!is.na(.data[[lvl_col]]), as.character(.data[[lvl_col]]), NA_character_))
    }
  }
  # Mettre à jour les colonnes code_clcm_lvl1
  data <- data %>%
    mutate(
      code_clcm_lvl1 = if_else(!is.na(code_clcm_lvl2) & nchar(code_clcm_lvl2) >= 1, substr(code_clcm_lvl2, 1, 1), code_clcm_lvl1),
    )
  
  
  # Mettre à jour les colonnes code_clcm_lvl1, code_clcm_lvl2 et code_clcm_lvl3
  data <- data %>%
    mutate(
      code_clcm_lvl1 = if_else(!is.na(code_clcm_lvl4) & nchar(code_clcm_lvl4) >= 1, substr(code_clcm_lvl4, 1, 1), code_clcm_lvl1),
      code_clcm_lvl2 = if_else(!is.na(code_clcm_lvl4) & nchar(code_clcm_lvl4) >= 2, substr(code_clcm_lvl4, 1, 2), code_clcm_lvl2),
      code_clcm_lvl3 = if_else(!is.na(code_clcm_lvl4) & nchar(code_clcm_lvl4) >= 3, substr(code_clcm_lvl4, 1, 3), code_clcm_lvl3)
    )

  # Écraser clcm_lvl1 si nécessaire
  if (any(!is.na(data$code_clcm_lvl1))) {
    data <- left_join(data, correspondance_clcm_lvl1 %>% select(code_clcm_lvl1, clcm_lvl1), by = "code_clcm_lvl1") %>%
      mutate(clcm_lvl1 = coalesce(clcm_lvl1.y, clcm_lvl1.x)) %>%
      select(-clcm_lvl1.x, -clcm_lvl1.y)
  }
  
  # Écraser clcm_lvl2 si nécessaire
  if (any(!is.na(data$code_clcm_lvl2))) {
    data <- left_join(data, correspondance_clcm_lvl2 %>% select(code_clcm_lvl2, clcm_lvl2), by = "code_clcm_lvl2") %>%
      mutate(clcm_lvl2 = coalesce(clcm_lvl2.y, clcm_lvl2.x)) %>%
      select(-clcm_lvl2.x, -clcm_lvl2.y)
  }
  
  # Écraser clcm_lvl3 si nécessaire
  if (any(!is.na(data$code_clcm_lvl3))) {
    correspondance_clcm_lvl3 <- correspondance_clcm_lvl3 %>% 
      distinct(code_clcm_lvl3, .keep_all = TRUE)
    data <- left_join(data, correspondance_clcm_lvl3 %>% select(code_clcm_lvl3, clcm_lvl3), by = "code_clcm_lvl3") %>%
      mutate(clcm_lvl3 = coalesce(clcm_lvl3.y, clcm_lvl3.x)) %>%
      select(-clcm_lvl3.x, -clcm_lvl3.y)
  }
  
      # Assigner les données modifiées au nom de la table
  assign(table_name, data, envir = .GlobalEnv)
}
  


#############
#############
#############

# Identifiez les noms de variables qui se terminent par "_description"
tableaux_final <- ls(pattern = "_description$")

# Fonction pour homogénéiser chaque dataframe
homogenize_dataframe <- function(df) {
  # Convertir toutes les colonnes en caractères
  df <- df %>%
    dplyr::mutate_all(as.character)

  
  # Colonnes spécifiées
  cols_to_keep <- c(
    "Programme","Annee", "Date_Prelevement", "ID_Site", "Modalite", "Bloc",
    "postal_code",	"gps_x",	"gps_y",	"Altitude",
    "clcm_lvl1",	"clcm_lvl2",	"clcm_lvl3", "clcm_lvl4","land_cover_detail", "code_clcm_lvl1", "code_clcm_lvl2", "code_clcm_lvl3", "code_clcm_lvl4",    "ph_kcl",	"ph_eau",	"c_tot",	"c_org",	"n_tot",	"c/n",	"om",	"cu_tot",	"cu_EDTA",	"soil_temperature",	"soil_humidity",	"fine_sand",	"coarse_sand",	"sand",	"fine_silt",	"coarse_silt",	"silt",	"clay",
    "type_tillage",	"tillage_depth",	"tillage_frequency_intra",	"tillage_frequency_inter",	"tillage_date",	
    "fertilisation",	"ferti_min_product",	"ferti_min_qtty",	"ferti_orga_product",	"ferti_orga_qtty",	"ferti_orga_freq",	
    "fungicide_freq",	"insecticide_freq",	"herbicide_freq",	"molluscicide_freq",	"nematicide_freq",	"tfi_fungicide",	"tfi_insecticide",	"tfi_herbicide",	"tfi_mollucicide",	"tfi_nematicide",	"total_tfi",
    "mecanical_weed_control",	"thermal_weed_control",	"crop_rotation_yr",	"rotation_plant_div",	"intercrop_div",	"rotation_grassland",	"crop_residues_management",	
    "amdmt_orga_freq",	"amdmt_orga_names",	"amdmt_orga_qtty",	"amdmt_calcic",	"amdmt_calcic_names",	"amdmt_calcic_qtty",	
    "herbage_use",	"mowing_frequency_yr",	"herb_age",	"animal_loading",	"trampling_nature", "grassland_type"
  )
  
  # Vérifier les colonnes manquantes et les ajouter avec des valeurs NA
  missing_cols <- setdiff(cols_to_keep, names(df))
  df[missing_cols] <- NA
  
  # Filtrer les colonnes pour ne garder que celles qui existent dans le dataframe
  df <- df %>%
    dplyr::select(all_of(cols_to_keep))
  
  return(df)
}


# Homogénéiser chaque dataframe et mettre à jour les tableaux de l'environnement
for (tableau_name in tableaux_final) {
  # Récupérer l'objet
  obj <- get(tableau_name)
  
  # Homogénéiser le dataframe
  df_homogenized <- homogenize_dataframe(obj)
  
  # Mettre à jour le tableau de l'environnement avec le dataframe homogénéisé
  assign(tableau_name, df_homogenized)
}

