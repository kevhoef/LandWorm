# Load the necessary libraries
library(readxl)
library(tidyverse)

# Specify the path to your Excel file
excel_file_path <- "[Database]/[Raw datasets]/LandWorm data providers/cp_data_site.xlsx"

# Read data from each sheet
sheet_data_Study <- read_excel(excel_file_path, sheet = "5. data Study") %>%
  dplyr::rename(GPS_Y = gps_x,GPS_X = gps_y) %>%
  mutate(ID_bloc = na_if(ID_bloc, "na")) %>%
  mutate(ID_modality = na_if(ID_modality, "na"))%>%
  mutate_all(~gsub(" ", "_", .))
sheet_data_Land_use_management <- read_excel(excel_file_path, sheet = "7.data Land_use_management")%>%
  mutate(ID_bloc = na_if(ID_bloc, "na")) %>%
  mutate(ID_modality = na_if(ID_modality, "na"))%>%
  mutate_all(~gsub(" ", "_", .))%>%
  dplyr::rename(clcm_lvl1 =land_use_lvl1,
         clcm_lvl2 = land_use_lvl2,
         clcm_lvl3 = land_use_lvl3,
         land_cover_detail = land_use_lvl4)
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
    EW_species == "A. chlorotica " ~ "Allolobophora_chlorotica_chlorotica",
    EW_species == "A. cupulifera " ~ "Aporrectodea_cupulifera",
    EW_species == "A. giardi " ~ "Aporrectodea_giardi",
    EW_species == "A. longa " ~ "Aporrectodea_longa",
    EW_species == "A. longa/giardi " ~ "Aporrectodea_longa/giardi",
    EW_species == "A. muldali " ~ "Murchieona_muldali",
    EW_species == "A. muldali/rosea " ~ "A._muldali/rosea",
    EW_species == "A. rosea " ~ "Aporrectodea_rosea",
    EW_species == "Allolobophora chlorotica" ~ "Allolobophora_chlorotica",
    EW_species == "Allolobophora muldali" ~ "Murchieona_muldali", 
    EW_species == "Anécique lumbricus " ~ "Lumbricus_indéterminable_anecic",
    EW_species == "Aporrectodea caliginosa" ~ "Aporrectodea_caliginosa_indéterminable",
    EW_species == "Aporrectodea giardi" ~ "Aporrectodea_giardi",
    EW_species == "Aporrectodea icterica" ~ "Aporrectodea_icterica",
    EW_species == "Aporrectodea longa" ~ "Aporrectodea_longa",
    EW_species == "Aporrectodea rosea" ~ "Aporrectodea_rosea",
    EW_species == "D. mammalis " ~ "Satchellius_mammalis",
    EW_species == "D. pygmea " ~ "Dendrobaena_pygmea",
    EW_species == "D. subrubicunda " ~ "Dendrodrilus_rubidus_subrubicundus",
    EW_species == "E. tetraedra " ~ "Eiseniella_tetraedra",
    EW_species == "Eisenia " ~ "Eisenia_indéterminable",
    EW_species == "endoge indetermine " ~ "indéterminable_endogeic",
    EW_species == "L. castaneus " ~ "Lumbricus_castaneus",
    EW_species == "L. centralis " ~ "Lumbricus_centralis",
    EW_species == "L. friendi " ~ "Lumbricus_friendi", 
    EW_species == "L. friendi/centralis " ~ "Lumbricus_friendi/centralis",
    EW_species == "L. rubellus friendoides " ~ "Lumbricus_rubellus_friendoides",
    EW_species == "L. rubellus rubellus " ~ "Lumbricus_rubellus_rubellus",
    EW_species == "Lumbricus castaneus" ~ "Lumbricus_castaneus",
    EW_species == "Lumbricus festivus" ~ "Lumbricus_festivus",
    EW_species == "Lumbricus terrestris" ~ "Lumbricus_terrestris",
    EW_species == "M. dubius " ~ "Microscolex_dubius",
    EW_species == "N. caliginosa " ~ "Aporrectodea_caliginosa_indéterminable",
    EW_species == "N. caliginosus meridionalis " ~ "Aporrectodea_caliginosa_meridionalis",
    EW_species == "N. terrestris terrestris " ~ "Aporrectodea_giardi",
    EW_species == "non identifie epige" ~ "Indéterminable_epigeic",
    EW_species == "O. cyaneum " ~ "Octolasion_cyaneum",
    EW_species == "O. lacteum " ~ "Octolasion_lacteum_lacteum",
    EW_species == "Octolasium " ~ "Octolasion_indéterminable",
    EW_species == "Octolasium cyaneum" ~ "Octolasion_cyaneum",
    EW_species == "P. fragilis " ~ "Prosellodrilus_fragilis_fragilis",
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
  mutate(BM_tot = rowSums(select(., starts_with("BM_")), na.rm = TRUE))


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


write.csv(cp_rd_rep_site, "[Database]/[Raw datasets]/EW_datasets/cp_rd_rep_site.csv", row.names = FALSE)










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
cp_description <- sheet_data_Study %>%
  left_join(sheet_data_Land_use_management, by = "ID_common") %>%
  left_join(sheet_data_Soil, by = "ID_common") %>%
  left_join(sheet_data_Sampling, by = "ID_common")


# Supprimez la colonne ID_common si vous n'en avez plus besoin
cp_description <- cp_description %>%
  select(-ID_common)


# Récupérez les noms des colonnes
noms_colonnes <- colnames(cp_description)

# Affichez les noms des colonnes
print(noms_colonnes)


cp_description <- cp_description %>%
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

cp_description <- cp_description %>%
  mutate_at(vars(Programme, ID_Site, Annee, Modalite, Bloc), as.character)


write.csv(cp_description, "[Database]/[Raw datasets]/Plot_description_datasets/cp_description.csv", row.names = FALSE)
