library(dplyr)

#### Nettoyage sbt_tb_rd => A faire à chaque nouveau JDD ####
sbt_tb_rd <- read_excel("[Database]/[Raw datasets]/UR_raw_datasets/SBT_ENI/ENI_Datas_Compilation_TB_2019-2023_V25.10.23.xlsx", sheet = "DATAS ENI")
sbt_codes = read_excel("[Database]/[Raw datasets]/UR_raw_datasets/SBT_ENI/codes_MNHN_SBT-ENI_V09.06.2023.xlsx", sheet = "MNHN -> UR1")
sbt_tb_rd <- left_join(sbt_tb_rd, sbt_codes %>% select(ID_site, id_parcelle), by = c("ID_Site" = "ID_site"))%>%
  dplyr::rename(ID_Site_UR1 = ID_Site) %>%
  dplyr::rename(ID_Site = id_parcelle) %>%
  mutate(Programme ="SBT-ENI-TB")
#  unite(ID_Site, id_parcelle,Annee, sep = "_", remove = FALSE) %>%
#  mutate(Site = id_parcelle) %>%
#  select(-id_parcelle)
write.table(sbt_tb_rd, "./[Database]/[Raw datasets]/EW_datasets/sbt_tb_rd.csv", sep=",", col.names=TRUE, row.names=TRUE, quote=TRUE, na="NA")



#### Chargement des meta avec id_parcelle et annee permettant de créer un ID unique avec les autres tableaux
sbt_tb_meta <- read_csv2("./[Database]/[Raw datasets]/UR_raw_datasets/SBT_ENI/Export_biovigilance_MNHN_15_06_23/obs_comptage.csv") %>%
  filter(protocole == "Vers de terre - BECHE") %>%
  filter(!(id_obs %in% c(26242, 32474,30423,23473,29316,30498))) %>% # doublons détectés, garde la dernière saisie en date
  mutate(Programme = "SBT-ENI") %>%
  dplyr::rename(ID_Site = id_parcelle, Annee = annee)


##### Chargement des meta avec coordonnées GPS 
sbt_parcelle= read_csv2("./[Database]/[Raw datasets]/UR_raw_datasets/SBT_ENI/Export_biovigilance_MNHN_15_06_23/parcelle.csv") %>%
  dplyr::rename(ID_Site = ID_Parcelle)


#####Joiture meta id_parcelle + coordonnées GPS
sbt_tb_meta=left_join(sbt_tb_meta, sbt_parcelle, by = "ID_Site") 

sbt_tb_meta <- sbt_tb_meta %>%
  mutate(Details_Milieu_Niv3 = case_when(# Créer la colonne Details_Milieu_Niv3 en fonction des conditions
    plante_reference %in% c("Blé", "Maïs") ~ "214_Culture annuelle",
    plante_reference %in% "Vigne" ~ "218_Vignes et autres Cultures pérennes",
    plante_reference %in% "Salade" ~ "216_Légume plein champ",
    TRUE ~ NA_character_
  )) %>%
  mutate(Categorie_Milieu_Niv1 = "2_Agricole")%>%# Créer la colonne Categorie_Milieu_Niv1
  
  mutate(SousCategorie_Milieu_Niv2 = "21_Agricole ouvert")# Créer la colonne SousCategorie_Milieu_Niv2



##### Chargement des meta avec analyses de sol
sbt_analyses_sol <- read_excel("[Database]/[Raw datasets]/UR_raw_datasets/SBT_ENI/ResultatsAnalyseSol_500ENI_2022_verif.xlsx", sheet = "Data") %>%
  filter(!(id_parcelle %in% c("Nouvelle"))) %>%
  mutate(id_parcelle = as.numeric(id_parcelle)) %>%
  dplyr::rename(ID_Site = id_parcelle) %>%
  clay
colnames(sbt_analyses_sol)
#####Joiture meta avec coordonnées GPS
sbt_tb_meta=left_join(sbt_tb_meta, sbt_analyses_sol, by = c("ID_Site"= "ID_Site"))

# Renommer les colonnes "coord_x" en "GPS_X" et "coord_y" en "GPS_Y" dans sbt_meta
colnames(sbt_tb_meta)[colnames(sbt_tb_meta) == "coord_x"] <- "GPS_X"
colnames(sbt_tb_meta)[colnames(sbt_tb_meta) == "coord_y"] <- "GPS_Y"

sbt_tb_meta$GPS_X <- as.numeric(gsub(",", ".", sbt_tb_meta$GPS_X))
sbt_tb_meta$GPS_Y <- as.numeric(gsub(",", ".", sbt_tb_meta$GPS_Y))

sbt_tb_meta <- sbt_tb_meta %>%
  dplyr::mutate(
    temp_X = if_else(GPS_X >= 40 & GPS_X <= 50, GPS_Y, GPS_X), #Il y a des inversion X et Y
    temp_Y = if_else(GPS_X >= 40 & GPS_X <= 50, GPS_X, GPS_Y),
    ID = paste(ID_Site, Annee, sep = "_"))%>%
  dplyr::mutate(Programme ="SBT-ENI-TB",
         GPS_X = temp_X,
         GPS_Y = temp_Y) %>%
  select(-temp_X, -temp_Y)

write.csv(sbt_tb_meta, "[Database]/[Raw datasets]/Plot_description_datasets/sbt_tb_meta.csv", row.names = FALSE)
