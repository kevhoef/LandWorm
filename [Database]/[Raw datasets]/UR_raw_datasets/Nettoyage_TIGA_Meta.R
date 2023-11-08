library(dplyr)
tiga_rd=read_csv("./[Database]/Raw datasets/tiga_rd.csv")
# Remplacer "URBA" par "URBAN" dans la colonne ID_Site
tiga_rd$ID_Site <- gsub("\\bURBA\\b", "URBAN", tiga_rd$ID_Site)
tiga_rd$ID_Site <- gsub("URBANNN", "URBAN", tiga_rd$ID_Site)

# Remplacer "RUR" par "RURAL" dans la colonne ID_Site
tiga_rd$ID_Site <- gsub("\\bRUR\\b", "RURAL", tiga_rd$ID_Site)
tiga_rd$ID_Site <- gsub("RURALAL", "RURAL", tiga_rd$ID_Site)

write.csv(tiga_rd, "[Database]/Raw datasets/tiga_rd.csv", row.names = FALSE)












tiga_meta=read_xlsx("./[Database]/Raw datasets/descriptive_data/ti_meta.xlsx")

# If you want to remove the date from the "code_ech" column, you can use sub() again with an empty string.
tiga_meta$code_ech <- str_replace(tiga_meta$code_ech, "_\\d{4}_", "_") 
tiga_meta <- tiga_meta %>%
  mutate(Programme = "TIGA")

tiga_meta = dplyr::rename(tiga_meta, ID_Site = code_ech)

# Remplacer "URBA" par "URBAN" dans la colonne ID_Site
tiga_meta$ID_Site <- gsub("URBA", "URBAN", tiga_meta$ID_Site)


# Remplacer le code parcelle qui déconne 
tiga_meta$ID_Site <- gsub("TIDM_GEN_GG_504", "TIDM_GEN_GC_504", tiga_meta$ID_Site)
tiga_meta$ID_Site <- gsub("TIDM_URBAN_VP_0111", "TIDM_URBAN_V _0020", tiga_meta$ID_Site)

library(sp)

tiga_meta <- tiga_meta %>%
  mutate(
    Details_Milieu_Niv3 = case_when(
      usage1 == "Agriculture_urbaine" ~ "340_Potager",
      usage1 == "Cultures" ~ "214_Culture annuelle",
      usage1 == "Forets" & localisation %in% c("BAC_Jeute", "Rural") ~ "113_Forêt mixte",
      usage1 == "Forets" & localisation == "urbain" ~ "310_Forêt urbaine",
      usage1 == "Loisir" ~ NA_character_,
      usage1 == "Maraichages" ~ "216m_Maraichage",
      usage1 == "Pelouse_calcaire" ~ NA_character_,
      usage1 == "Prairies" ~ "231_Prairie",
      usage1 == "Vignes" ~ "218_Vignes et autres Cultures pérennes",
      usage1 == "Voirie" ~ NA_character_,
      TRUE ~ NA_character_
    ),
    SousCategorie_Milieu_Niv2 = case_when(
      usage1 %in% c("Cultures", "Maraichages", "Prairies", "Vignes") ~ "21_Agricole ouvert",
      usage1 == "Forets" & localisation %in% c("BAC_Jeute", "Rural") ~ "11_Naturel fermé",
      usage1 == "Forets" & localisation == "urbain" ~ "31_Espace vert boisé",
      usage1 == "Loisir" ~ "32_Espace vert ouvert",
      usage1 == "Pelouse_calcaire" ~ "32_Espace vert ouvert",
      usage1 == "Voirie" ~ "35_Bordure de voirie",
      usage1 == "Agriculture_urbaine" ~ "34_Espace cultivé urbain",
      TRUE ~ NA_character_
    ),
    Categorie_Milieu_Niv1 = case_when(
      usage1 %in% c("Cultures", "Maraichages", "Prairies", "Vignes") ~ "2_Agricole",
      usage1 == "Forets" & localisation %in% c("BAC_Jeute", "Rural") ~ "1_Naturel",
      usage1 == "Forets" & localisation == "urbain" ~ "3_Artificialisé",
      usage1 == "Loisir" ~ "3_Artificialisé",
      usage1 == "Pelouse_calcaire" ~ "3_Artificialisé",
      usage1 == "Voirie" ~ "3_Artificialisé",
      usage1 == "Agriculture_urbaine" ~ "3_Artificialisé",
      TRUE ~ NA_character_
    )
  )



# Mise à jour pour les lignes 101 à 140
tiga_meta$Annee[101:140] <- "2021"
# Mise à jour pour les lignes 141 à 170
tiga_meta$Annee[141:170] <- "2022"

# Renommer les colonnes "coord_x" en "GPS_X" et "coord_y" en "GPS_Y" dans sbt_meta
colnames(tiga_meta)[colnames(tiga_meta) == "X_L93"] <- "GPS_X"
colnames(tiga_meta)[colnames(tiga_meta) == "Y_L93"] <- "GPS_Y"

tiga_meta$GPS_X <- as.numeric(gsub(",", ".", tiga_meta$GPS_X))
tiga_meta$GPS_Y <- as.numeric(gsub(",", ".", tiga_meta$GPS_Y))


# Définir le système de coordonnées de départ (Lambert 93) et le système de coordonnées d'arrivée (WGS84)
src_crs <- CRS("EPSG:2154")  # Lambert 93 (EPSG:2154)
dst_crs <- CRS("+init=epsg:4326")  # WGS84 (EPSG:4326)


# Créer un objet "SpatialPoints" à partir des colonnes GPS_X et GPS_Y du tableau rmqs_meta
coords <- SpatialPoints(coords = tiga_meta[, c("GPS_X", "GPS_Y")], proj4string = src_crs)

# Effectuer la transformation de Lambert 93 vers WGS84 (degrés décimaux)
coords_transformed <- spTransform(coords, dst_crs)

# Mettre à jour les colonnes GPS_X et GPS_Y du tableau rmqs_meta avec les nouvelles valeurs transformées
tiga_meta$GPS_X <- coords_transformed@coords[, 1]
tiga_meta$GPS_Y <- coords_transformed@coords[, 2]

write.csv(tiga_meta, "[Database]/Raw datasets/descriptive_data/tiga_meta.csv", row.names = FALSE)
