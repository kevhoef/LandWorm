#### Nettoyage sbt_m_rd => A faire à chaque nouveau JDD ####
#sbt_m_rd = read_csv("./[Database]/Raw datasets/sbt_m_rd.csv")
#sbt_m_rd <- read_excel("[Database]/Raw datasets/National/SBT_ENI/ENI _ Data Compilation Moutarde 2012-2018 V2023.03.17.xlsx", sheet = "Compilation data_Brut")
#sbt_codes = read_excel("[Database]/Raw datasets/National/SBT_ENI/codes_MNHN_SBT-ENI_V09.06.2023.xlsx", sheet = "MNHN -> UR1")
#sbt_m_rd <- left_join(sbt_m_rd, sbt_codes %>% select(ID_site, id_parcelle), by = c("ID_Site" = "ID_site"))%>%
#  rename(ID_Site_UR1 = ID_Site) %>%
#  rename(ID_Site = id_parcelle)
#  unite(ID_Site, id_parcelle,Annee, sep = "_", remove = FALSE) %>%
#  mutate(Site = id_parcelle) %>%
#  select(-id_parcelle)
#write.table(sbt_m_rd, "./[Database]/Raw datasets/sbt_m_rd.csv", sep=",", col.names=TRUE, row.names=TRUE, quote=TRUE, na="NA")



library(dplyr)

#### Nettoyage sbt_m_rd => A faire à chaque nouveau JDD ####
sbt_m_rd <- read_excel("[Database]/[Raw datasets]/UR_raw_datasets/SBT_ENI/ENI _ Data Compilation Moutarde 2012-2018 V2023.03.17.xlsx", sheet = "Compilation data_Brut")
sbt_codes = read_excel("[Database]/[Raw datasets]/UR_raw_datasets/SBT_ENI/codes_MNHN_SBT-ENI_V09.06.2023.xlsx", sheet = "MNHN -> UR1")
sbt_m_rd <- left_join(sbt_m_rd, sbt_codes %>% select(ID_site, id_parcelle), by = c("ID_Site" = "ID_site"))%>%
  dplyr::rename(ID_Site_UR1 = ID_Site) %>%
  dplyr::rename(ID_Site = id_parcelle) %>%
  mutate(Programme = "SBT-ENI-M")
#  unite(ID_Site, id_parcelle,Annee, sep = "_", remove = FALSE) %>%
#  mutate(Site = id_parcelle) %>%
#  select(-id_parcelle)
write.table(sbt_m_rd, "./[Database]/[Raw datasets]/EW_datasets/sbt_m_rd.csv", sep=",", col.names=TRUE, row.names=TRUE, quote=TRUE, na="NA")



#### Chargement des meta avec id_parcelle et annee permettant de créer un ID unique avec les autres tableaux
sbt_m_meta <- read_csv2("./[Database]/[Raw datasets]/UR_raw_datasets/SBT_ENI/Export_biovigilance_MNHN_15_06_23/obs_comptage.csv") %>%
  filter(protocole == "Vers de terre - MOUTARDE") %>%
  filter(!(id_obs %in% c(22379, 4321))) %>%
  mutate(Programme = "SBT-ENI") %>%
  dplyr::rename(ID_Site = id_parcelle, Annee = annee)
  


##### Chargement des meta avec coordonnées GPS 
sbt_parcelle= read_csv2("./[Database]/[Raw datasets]/UR_raw_datasets/SBT_ENI/Export_biovigilance_MNHN_15_06_23/parcelle.csv") %>%
  dplyr::rename(ID_Site = ID_Parcelle)


#####Joiture meta id_parcelle + coordonnées GPS
sbt_m_meta=left_join(sbt_m_meta, sbt_parcelle, by = "ID_Site") 
  
sbt_m_meta <- sbt_m_meta %>%
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
  dplyr::rename(ID_Site = id_parcelle) 

#####Joiture meta avec coordonnées GPS
sbt_m_meta=left_join(sbt_m_meta, sbt_analyses_sol, by = c("ID_Site"= "ID_Site"))

# Renommer les colonnes "coord_x" en "GPS_X" et "coord_y" en "GPS_Y" dans sbt_meta
colnames(sbt_m_meta)[colnames(sbt_m_meta) == "coord_x"] <- "GPS_X"
colnames(sbt_m_meta)[colnames(sbt_m_meta) == "coord_y"] <- "GPS_Y"

sbt_m_meta$GPS_X <- as.numeric(gsub(",", ".", sbt_m_meta$GPS_X))
sbt_m_meta$GPS_Y <- as.numeric(gsub(",", ".", sbt_m_meta$GPS_Y))

sbt_m_meta <- sbt_m_meta %>%
  dplyr::mutate(
    temp_X = if_else(GPS_X >= 40 & GPS_X <= 50, GPS_Y, GPS_X), # Il y a des inversions X et Y
    temp_Y = if_else(GPS_X >= 40 & GPS_X <= 50, GPS_X, GPS_Y),
    ID = paste(ID_Site, Annee, sep = "_")
  ) %>%
  dplyr::mutate(
    Programme = "SBT-ENI-M",
    GPS_X = if_else(temp_X == 1.000000, NA_real_, temp_X), # Remplacer par NA si GPS_X est 1.000000
    GPS_Y = if_else(temp_Y == 1.000000, NA_real_, temp_Y)  # Remplacer par NA si GPS_Y est 1.000000
  ) %>%
  select(-temp_X, -temp_Y)

write.csv(sbt_m_meta, "[Database]/[Raw datasets]/Plot_description_datasets/sbt_m_meta.csv", row.names = FALSE)

# Afficher le dataframe sbt_meta après les modifications
print(sbt_m_meta)







sbt_travail_sol= read_csv2("./[Database]/Raw datasets/National/SBT_ENI/Export_biovigilance_MNHN_15_06_23/pcu_travail_sol.csv")
sbt_travail_sol <- sbt_travail_sol %>%
  separate(date_pratique, into = c("annee", "mois", "jour"), sep = "-", remove = FALSE)
sbt_comptage= read_csv2("./[Database]/Raw datasets/National/SBT_ENI/Export_biovigilance_MNHN_15_06_23/obs_comptage.csv")
sbt_comptage <- sbt_comptage %>%
  filter(protocole == "Vers de terre - MOUTARDE")
sbt_comptage$date_echantillonnage <- paste(sbt_comptage$annee, sbt_comptage$mois, sbt_comptage$jour, sep = "-")


# Supposons que sbt_comptage est votre dataframe
valeurs_a_supprimer <- c(4878, 15603, 4321, 15606, 25955, 3721, 17207, 9826, 8118, 16947)

# Supprimer les lignes ayant les valeurs d'id_obs spécifiques
sbt_comptage <- filter(sbt_comptage, !(id_obs %in% valeurs_a_supprimer))

# Supposons que votre dataframe s'appelle df
doublons <- sbt_comptage[duplicated(sbt_comptage$ID) | duplicated(sbt_comptage$ID, fromLast = TRUE), ]


sbt_comptage$ID <- paste(sbt_comptage$id_parcelle, sbt_comptage$annee, sep = "_")
sbt_travail_sol$ID <- paste(sbt_travail_sol$ID_Parcelle, sbt_travail_sol$annee, sep = "_")





sbt_comptage_subset <- sbt_comptage[, c("ID", "date_echantillonnage")]

# Joindre les deux tableaux en utilisant ID comme clé de jointure
jointure <- merge(sbt_travail_sol, sbt_comptage_subset, by = "ID", all.x = TRUE)

# Utilisez anti_join pour obtenir les lignes de sbt_comptage qui n'ont pas de correspondance dans jointure
non_liees <- anti_join(sbt_travail_sol, jointure, by = "ID")

# Affichez les lignes non liées
print(non_liees)


library(dplyr)

# Convertir les dates en objets Date
jointure$date_echantillonnage <- as.Date(jointure$date_echantillonnage)
jointure$date_pratique <- as.Date(jointure$date_pratique)

# Calcul de la différence en jours entre date_echantillonnage et date_pratique
jointure <- jointure %>%
  mutate(difference_jours = as.numeric(difftime(date_echantillonnage, date_pratique, units = "days")))

# Ajouter une nouvelle colonne basée sur la différence en jours
jointure <- jointure %>%
  mutate(resultat = case_when(
    difference_jours < -365 ~ "sup_1an",     # Antérieur de plus de 1 an
    difference_jours >= -365 & difference_jours < 0 ~ "moins_1an",  # Antérieur de moins de 1 an
    difference_jours >= 0 ~ "futur"        # Supérieur ou égal à 0 jour (futur)
  ))

# Supprimer la colonne difference_jours si elle n'est plus nécessaire
jointure <- jointure %>%
  select(-difference_jours)


jointure_moins_1an <- jointure %>%
  filter(resultat == "moins_1an")


library(dplyr)

# Groupement par "ID" et création de la colonne "type_tillage"
resultat_grouped <- jointure_moins_1an %>%
  group_by(ID) %>%
  summarize(type_tillage = paste(nom_pratique_agricole, collapse = "_"))

# Maintenant, resultat_grouped contient les valeurs de "nom_pratique_agricole" regroupées par "ID" dans "type_tillage"

# Supposons que vous ayez déjà créé les dataframes resultat_grouped et sbt_m_meta

# Utilisez merge pour joindre les dataframes par la colonne "ID" commune
sbt_m_meta <- left_join(sbt_m_meta, resultat_grouped, by = "ID")

# Le dataframe jointure_finale contient maintenant la jointure des deux dataframes par la colonne "ID" commune










#Joiture DATA VDT  + META
# Créez une nouvelle colonne combinant les trois colonnes avec "_"

cp_meta <- cp_meta %>%
  mutate_at(vars(Programme, ID_Site, Annee, Modalite, Bloc), as.character)

cp_rd_rep_site <- cp_rd_rep_site %>%
  mutate_at(vars(Programme, ID_Site, Annee, Modalite, Bloc), as.character)

cp_meta <- cp_meta %>%
  mutate(identifiant = paste(Programme, ID_Site, Annee, Modalite, Bloc, sep = "_"))
cp_rd_rep_site <- cp_rd_rep_site %>%
  mutate(identifiant = paste(Programme, ID_Site, Annee, Modalite, Bloc, sep = "_"))


str(sbt_m_meta)
str(sbt_m_rd_rep_site)

# Effectuez la jointure en utilisant left_join
resultat <- right_join(sbt_m_meta,sbt_m_rd_rep_site, by = "identifiant")

# Supprimez la colonne identifiant si vous le souhaitez
resultat <- resultat %>%
  select(-identifiant)

colnames(resultat)


library(dplyr)

head(cp_meta$identifiant)
head(cp_rd_rep_site$identifiant)
resultat <- left_join(sbt_m_meta, sbt_m_rd_rep_site, by = c("identifiant" = "ID"))
resultat <- left_join(sbt_m_meta, sbt_m_rd_rep_site, by = c("identifiant" = "ID_rep"))

