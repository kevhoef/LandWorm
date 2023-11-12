# Load the necessary libraries
library(readxl)
library(tidyverse)

# Specify the path to your Excel file
excel_file_path <- "[Database]/[Raw datasets]/LandWorm data providers/Habitat_and_Earthworm_data_and_explanations.xlsx"

# Read data from each sheet
cp_biobio_expla <- read_excel(excel_file_path, sheet = "HabitatType_Explanation") # Pour avoir la signification des codes habitats
cp_biobio_description <- read_excel(excel_file_path, sheet = "Habitat_data_France")%>% # Pour que les colonnes du tableau corespondent aux scripts de concaténation
  dplyr::rename(gps_x = XCoord_WGS1984,gps_y = YCoord_WGS1984, ID_Site = Site)%>%
  filter(str_starts(ID_Site, "FR"))%>%
  mutate(Programme = "biobio",
         Annee = "2010")
cp_biobio_rd <- read_excel(excel_file_path, sheet = "Earthworm_data")%>% # Pour que les colonnes du tableau corespondent aux scripts de concaténation
  dplyr::rename(ID_Site = Site, Repetition = Location, Code_Methode =Method, Nbr_VDT = Abundance, Code_Taxon = Species)%>%
  filter(str_starts(ID_Site, "FR")) %>%
  mutate(Code_Methode = case_when(
    Code_Methode == "X" ~ "AITC_11.1111111111111",
    Code_Methode == "Y" ~ "TM",
    TRUE ~ Code_Methode),
    Programme = "biobio",
    Protocole = "AITCTM",
    Stade = "X",
    Code_Taxon = ifelse(Code_Taxon == "zero.pseudo", NA, Code_Taxon),
    Annee = "2010")
cp_biobio_rd <- cp_biobio_rd %>% # A REVOIR
  mutate(Code_Taxon = case_when(
    Code_Taxon == "Aporrectodea caliginosa subsp caliginosa" ~ "Aporrectodea_caliginosa_caliginosa",
    Code_Taxon == "Aporrectodea caliginosa subsp meridionalis" ~ "Aporrectodea_caliginosa_meridionalis",
    Code_Taxon == "Prosellodrilus fragilis" ~ "Prosellodrilus_fragilis_indéterminable",
    Code_Taxon == "Allolobophora chlorotica" ~ "Allolobophora_chlorotica_indéterminable",
    Code_Taxon == "Aporrectodea rosea" ~ "Aporrectodea_rosea_indéterminable",
    Code_Taxon == "Dendrobaena_mammalis" ~ "Satchellius_mammalis",
    Code_Taxon == "Allolobophora_muldali" ~ "Murchieona_muldali",
    TRUE ~ Code_Taxon  # Si aucune condition n'est satisfaite, conservez la valeur d'origine
  )) %>%
  mutate(Code_Taxon = as.factor(Code_Taxon))%>%
  mutate(Code_Taxon = str_replace_all(Code_Taxon, " ", "_"))

cp_biobio_description <- cp_biobio_description %>% # Pour récupérer le corine land cover à partir de la description des habitats
  left_join(select(cp_biobio_expla, Code, Category, Description), by = c("HabitatType" = "Code")) %>%
  mutate(
    clcm_lvl3 = case_when(
      Description %in% c("Lines of scrub", "Grass strips", "Herbaceous strips") ~ "213_Bande enherbée",
      Description == "Lines of trees" ~ "232_Haies",
      Description == "Private track grass strips" ~ "219_Autre",
      Description == "Woody crops vines" ~ "218_Vignes et autres Cultures pérennes",
      Description %in% c("Not entomophilic and/or bee attracting annual spring crops", "Not entomophilic and/or bee attracting annual winter crops", "Perennials, e.g. forage crops", "Entomophilic and/or bee attracting annual crops") ~ "214_Culture annuelle", #"Perennials, e.g. forage crops" parce que c'est dans "cultivated"
      Category == "Vegetated herbaceous" ~ "120_Végétation clairsemée",
      Description == "Forest phanerophytes winter deciduous" ~ "111_Forêt de feuillus",
      Description == "Forest phanerophytes coniferous" ~ "112_Forêt de conifères",
      Description %in% c("Tall phanerophytes winter deciduous", "Mid phanerophytes winter deciduous", "Mid phanerophytes evergreen","Tall phanerophytes coniferous","Low phanerophytes winter deciduous","Tall phanerophytes evergreen") ~ "115_Autre",
    TRUE ~ NA_character_),
    clcm_lvl2 = case_when(
      Description %in% c("Lines of scrub", "Grass strips", "Herbaceous strips") ~ "21_Agricole ouvert",
      Description == "Lines of trees" ~ "22_Agricole boisé",
      Description == "Private track grass strips" ~ "21_Agricole ouvert",
      Description == "Woody crops vines" ~ "21_Agricole ouvert",
      Description %in% c("Not entomophilic and/or bee attracting annual spring crops", "Not entomophilic and/or bee attracting annual winter crops", "Perennials, e.g. forage crops", "Entomophilic and/or bee attracting annual crops") ~ "21_Agricole ouvert", #"Perennials, e.g. forage crops" parce que c'est dans "cultivated"
      Category == "Vegetated herbaceous" ~ "12_Naturel ouvert",
      Description == "Forest phanerophytes winter deciduous" ~ "11_Naturel fermé",
      Description == "Forest phanerophytes coniferous" ~ "12_Naturel ouvert",
      Description %in% c("Tall phanerophytes winter deciduous", "Mid phanerophytes winter deciduous", "Mid phanerophytes evergreen","Tall phanerophytes coniferous","Low phanerophytes winter deciduous","Tall phanerophytes evergreen") ~ "11_Naturel fermé",
      TRUE ~ NA_character_),
    clcm_lvl1 = case_when(
      Description %in% c("Lines of scrub", "Grass strips", "Herbaceous strips") ~ "2_Agricole",
      Description == "Lines of trees" ~ "2_Agricole",
      Description == "Private track grass strips" ~ "2_Agricole",
      Description == "Woody crops vines" ~ "2_Agricole",
      Description %in% c("Not entomophilic and/or bee attracting annual spring crops", "Not entomophilic and/or bee attracting annual winter crops", "Perennials, e.g. forage crops", "Entomophilic and/or bee attracting annual crops") ~ "2_Agricole", #"Perennials, e.g. forage crops" parce que c'est dans "cultivated"
      Category == "Vegetated herbaceous" ~ "1_Naturel",
      Description == "Forest phanerophytes winter deciduous" ~ "1_Naturel",
      Description == "Forest phanerophytes coniferous" ~ "1_Naturel",
      Description %in% c("Tall phanerophytes winter deciduous", "Mid phanerophytes winter deciduous", "Mid phanerophytes evergreen","Tall phanerophytes coniferous","Low phanerophytes winter deciduous","Tall phanerophytes evergreen") ~ "1_Naturel",
      TRUE ~ NA_character_)) %>%
  dplyr::rename(land_cover_detail = Description)


write.csv(cp_biobio_rd, "[Database]/[Raw datasets]/EW_datasets/cp_biobio_rd.csv", row.names = FALSE)
write.csv(cp_biobio_description, "[Database]/[Raw datasets]/Plot_description_datasets/cp_biobio_description.csv", row.names = FALSE)

