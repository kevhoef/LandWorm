# Load the necessary libraries
library(readxl)
library(tidyverse)

# Specify the path to your Excel file
excel_file_path <- "[Database]/Raw datasets/LandWorm data providers/LandWorm_data_template_V15_05_23-MHEDDE.xlsx"

# Read data from each sheet
sheet_data_Study <- read_excel(excel_file_path, sheet = "5. data Study")%>%
  mutate_all(~gsub("-", "_", .))
sheet_data_Land_use_management <- read_excel(excel_file_path, sheet = "7.data Land_use_management")%>%
  mutate_all(~gsub("-", "_", .))
sheet_data_Soil <- read_excel(excel_file_path, sheet = "9.data Soil")%>%
  mutate_all(~gsub("-", "_", .))
sheet_data_Sampling <- read_excel(excel_file_path, sheet = "11. data Sampling")%>%
  mutate_all(~gsub("-", "_", .))%>%
  mutate(date_plot_collection = as.character(date_plot_collection))
mh_rd <- read_excel(excel_file_path, sheet = "13. Data Earth. com.")%>%
  mutate_all(~gsub("-", "_", .))


##### EARTHWORM COMMUNITIES
#Récupération de la surface d'échantillonnage
#Probl : les identifiants project_acronym	sampling_method_acronym	ID_site	ID_plot... ne corespondent pas entre les feuilles excels
#Les surface d'échantillonnage sont similaires si on comnine project_acronym et date_plot_collection
data_sampling_unique <- sheet_data_Sampling %>%
  group_by(project_acronym, date_plot_collection) %>%
  summarise(sampled_area_for_each_subsample = first(sampled_area_for_each_subsample))
#Jointure
mh_rd <- mh_rd %>%
  left_join(data_sampling_unique, 
            by = c("project_acronym", "date_plot_collection"))


# Extract and replace the terms in EW_species and update stage_dev
mh_rd <- mh_rd %>%
  mutate(stage_dev = recode(stage_dev, 
                            "Subadult" = "SA",
                            "Adult" = "AD",
                            "Juvenile" = "JV",
                            .default = NA_character_)) %>%
  filter(!(stage_dev == "cocon")) %>%
  filter(!(sampling_method_acronym  %in% c("other (qualitative HS)" , "other"))) %>%
  mutate(stage_dev = as.factor(stage_dev))%>%
  mutate(EW_species = as.factor(EW_species)) %>%
  mutate(biomass = as.numeric(biomass))



levels(mh_rd$EW_species)
mh_rd <- mh_rd %>% # A REVOIR
  mutate(EW_species = case_when(
    EW_species == "Allolobophora chlorotica (Savigny 1826)" ~ "AC",
    EW_species == "Allolobophora cupulifera (Tétry 1937)" ~ "ACU",
    EW_species == "Allolobophoridella eiseni (Levinsen 1884)" ~ "AE",
    EW_species == "Aporrectodea" ~ "NX",
    EW_species == "Aporrectodea giardi (Ribaucourt 1910)" ~ "NG",
    EW_species == "Aporrectodea longa (Ude 1885)" ~ "NLX",
    EW_species == "Aporrectodea terrestris (Savigny, 1826)" ~ "ATerrestris",
    EW_species == "Aporrectodea rosea (Savigny 1826)" ~ "ARR",
    EW_species == "Aporrectodea rubra" ~ "ARubra",
    EW_species == "Aporrectodea trapezoides (Dugès, 1828)" ~ "Atrapez", 
    EW_species == "Bimastos eiseni" ~ "Beiseni",
    EW_species == "Aporrectodea caliginosa (Savigny 1826)" ~ "NCX",
    EW_species == "Aporrectodea icterica (Savigny 1826)" ~ "AI",
    EW_species == "Aporrectodea rosea" ~ "ARR",
    EW_species == "Dendrobaena" ~ "DX",
    EW_species == "Dendrobaena alpina (Rosa, 1884)" ~ "DA?",
    EW_species == "Dendrobaena alpina juvenile (Rosa, 1884)" ~ "DA?",
    EW_species == "Dendrobaena mammalis (Savigny 1826)" ~ "DM",
    EW_species == "Dendrobaena octaedra (Savigny 1826)" ~ "DO",
    EW_species == "Dendrodrilus rubidus (Savigny 1826)" ~ "DR?",
    EW_species == "Dendrodrilus subrubicundus (Eisen, 1874)" ~ "DSS?",
    EW_species == "Dendrodrilus rubidus (Savigny 1826)" ~ "DR?",
    EW_species == "Eiseniella tetraedra  (Savigny 1826)" ~ "ET?",
    EW_species == "Lumbricus castaneus (Savigny 1826)" ~ "LC",
    EW_species == "Lumbricus centralis Bouché 1972" ~ "LCE",
    EW_species == "Lumbricus festivus (Savigny 1826)" ~ "LFE", 
    EW_species == "Lumbricus friendi friendi Cognetti 1804" ~ "LFR", # verif
    EW_species == "Lumbricus Linnaeus 1758" ~ "LLinnaeus",
    EW_species == "L. rubellus rubellus " ~ "LRR",
    EW_species == "Lumbricus meliboeus Rosa, 1884" ~ "LM?", #
    EW_species == "Lumbricus rubellus castanoides Bouché 1972" ~ "LRC",
    EW_species == "Lumbricus rubellus rubellus Hoffmeister 1843" ~ "LRR",
    EW_species == "Lumbricus terrestris L. 1758" ~ "LT",
    EW_species == "Microscolex" ~ "MX", #verif
    EW_species == "Microscolex dubius (Fletcher, 1887)" ~ "MD",
    EW_species == "Microscolex phosphoreus (Dugès, 1837)" ~ "MP",
    EW_species == "N. caliginosus meridionalis " ~ "NCM",
    EW_species == "Murchieona minuscula (Rosa 1905/6)" ~ "MM", #vérif
    EW_species == "Octolasion cyaneum (Savigny 1826)" ~ "OC",
    EW_species == "Octolasion lacteum (Orley, 1885)" ~ "OL",
    EW_species == "Octolasion Orley, 1885" ~ "OX", #verif
    EW_species == "Octolasion tyrtaeum (Savigny 1826)" ~ "OT",#verif
    EW_species == "Prosellodrilus Bouché, 1972" ~ "PX",# vérif
    EW_species == "Oligochaeta" ~ "X",# vérif
    EW_species == "Satchellius mammalis (Savigny 1826)" ~ "DM",# vérif
    EW_species == "Scherotheca Bouché, 1972" ~ "SX",# vérif
    EW_species == "Scherotheca (Opothedrilus) chicharia (Bouché, 1967)" ~ "SC",# vérif
    EW_species == "Scherotheca (Scherotheca) mifuga (Bouche, 1972)" ~ "SM",# vérif
    EW_species == "Scherotheca gigas (Duges, 1828)" ~ "SG",# vérif
    EW_species == "Scherotheca rhodana Bouche, 1972" ~ "SR",# vérif
    TRUE ~ EW_species  # Si aucune condition n'est satisfaite, conservez la valeur d'origine
  )) %>%
  mutate(EW_species = as.factor(EW_species))

#### RENOMMER LES VARIABLES
mh_rd <- mh_rd %>% rename(
    Programme = project_acronym,
    Protocole = sampling_method_acronym,
    Parcelle = ID_site,
    Annee= date_plot_collection,
    Modalite = ID_modality,
    Bloc = ID_bloc,
    ID_Site = ID_plot,
    Repetition = ID_subsample,
    Code_Taxon = EW_species,
    Stade = stage_dev,
    Nbr_VDT = abundance,
    Pds = biomass
    ) %>%
  mutate(Protocole = case_when(
    Protocole == "HS" & sampled_area_for_each_subsample == 0.25 ~ "HS_4",
    Protocole == "HS" & sampled_area_for_each_subsample == 0.0625 ~ "HS_16",
    Protocole == "HSAITC" & sampled_area_for_each_subsample == 0.1256636 ~ "HSAITC_7.95775385",
    Protocole == "HSAITC" & sampled_area_for_each_subsample == 0.16 ~ "HSAITC_6.25",
    Protocole == "HSAITC" & sampled_area_for_each_subsample == 0.25 ~ "HSAITC_4",
    Protocole == "HSAITC" & sampled_area_for_each_subsample == 0.0625 ~ "HSAITC_16",
    Protocole == "HS" & sampled_area_for_each_subsample == 0.25 ~ "HS_4",
    TRUE ~ Protocole
  ))%>%
  mutate(Code_Methode = Protocole)
mh_rd = mh_rd%>% mutate(Protocole = as.factor(Protocole))

levels(mh_rd$Protocole)

write.csv(mh_rd, "[Database]/Raw datasets/mh_rd.csv", row.names = FALSE)


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
mh_meta <- sheet_data_Study %>%
  left_join(sheet_data_Land_use_management, by = "ID_common") %>%
  left_join(sheet_data_Soil, by = "ID_common") %>%
  left_join(sheet_data_Sampling, by = "ID_common")


# Supprimez la colonne ID_common si vous n'en avez plus besoin
mh_meta <- mh_meta %>%
  select(-ID_common)


# Récupérez les noms des colonnes
noms_colonnes <- colnames(mh_meta)

# Affichez les noms des colonnes
print(mh_meta)


mh_meta <- mh_meta %>%
  rename(
    Programme = project_acronym,
    Protocole = sampling_method_acronym,
    Parcelle = ID_site,
    Annee= date_plot_collection,
    Modalite = ID_modality,
    Bloc = ID_bloc,
    ID_Site=ID_plot,
    )

mh_meta <- mh_meta %>%
  mutate_at(vars(Programme, ID_Site, Annee, Modalite, Bloc, Parcelle), as.character)


write.csv(mh_meta, "[Database]/Raw datasets/descriptive_data/mh_meta.csv", row.names = FALSE)
