##############################################################
####### Changing the nature of numerical variables ########
##############################################################
LandWorm_dataset_site_m <- LandWorm_dataset_site %>%
  dplyr::mutate(across(c(gps_x, gps_y, Altitude, ph_kcl, ph_eau, c_tot, c_org, n_tot, om,
                         cu_tot, cu_EDTA, soil_temperature, soil_humidity, fine_sand, coarse_sand, sand,
                         fine_silt, coarse_silt, silt, clay), ~ as.numeric(gsub("[^0-9.]", "", .))))



##############################################################
####### Correction of CLC according to classification ########
##############################################################
#### CLC LVL 1 ####
LandWorm_dataset_site_m <- LandWorm_dataset_site_m %>%
  mutate(clc_lvl3 = case_when(
    clcm_lvl3 %in% c("214_Culture_annuelle","214c_Culture_annuelle") ~ "211_Arable land*",
    clcm_lvl3 %in% c("218_Vignes_et_autres_Cultures_pérennes") ~ "221_Vineyards",
    clcm_lvl3 %in% c("220ir_Agroforesterie", "220r_Agroforesterie") ~ "244_Agro-forestry areas",
    clcm_lvl3 %in% c("Prairie_agricole","130_Prairie_humide","210_Prairie_agricole_permanente", "211_Prairie_agricole_temporaire","210_Prairie_agricole_permanente?","121_Prairie_naturelle_&_Paturage", "131_Mégaphorbiaie","231_Prairies","Prairie agricole") ~ "231_Pastures, meadows and other permanent grasslands under agricultural use",
    clcm_lvl3 %in% c("111_Forêt_de_feuillus") ~ "311_Broad-leaved forest",
    clcm_lvl3 %in% c("113_Forêt_mixte","310_Forêt_urbaine") ~ "313_Mixed forest",
    clcm_lvl3 %in% c("112_Forêt_de_conifères") ~ "312_Coniferous forest",
    clcm_lvl3 %in% c("221_Arboriculture","221_Verger","312_Verger") ~ "222_Fruit trees and berry plantations",
    clcm_lvl3 %in% c("114c_Bois_ou_forêt_contaminé", "114_Bois","311_Bois_urbain") ~ "324_Transitional woodland-shrub",
    clcm_lvl3 %in% c("321_Pelouse","380_Friche_industrielle", "320_Prairie_urbaine","321_Pelouse_urbaine","391_Bassin_de_rétention", "330_Massif_ornemental","351_Massif_ornemental_contraint") ~ "141a_Green urban areas",
    clcm_lvl3 %in% c("340_Potager","340p_Potager-pelouse","340p_Potager-Pelouse") ~ "141b_Urban agriculture*",
    clcm_lvl3 %in% c("212_Jachère","215_Interculture","232_Haies","213_Bande_enherbée","216l_Légume_plein_champ","216_Maraîchage_et_légume_plein_champ","230_Talus_et_bordure_enherbés","217_Serre,_tunnel_fixe","232_Haie","216m_Maraichage") ~ "242_Complex cultivation patterns",
    clcm_lvl3 %in% c("370_Stock_de_terre") ~ "133_Construction sites",
    clcm_lvl3 %in% c("120_Végétation_clairsemée","132_Lande_hygrophile","122_Lande_mésophile_et_broussailles") ~ "322_Moors and heathland",
    clcm_lvl3 %in% c("356_Noue_arborée","355_Noue_enherbée","352_Alignement_d'arbres_en_fosse_continue","350_Bande_enherbée_contrainte","353_Arbre_en_fosse_unitaire","354_Rond-point") ~ "122_Road and rail networks and associated land",
    clcm_lvl3 %in% c("324_Prairie_aéroport") ~ "124_Airports",
    clcm_lvl3 %in% c("322_Terrain_récréatif") ~ "142_Sport and leisure facilities",
    clcm_lvl3 %in% c("124_Plage,_dune_et_sable") ~ "331_Beaches, dunes, sands",
    clcm_lvl3 %in% c("343_Autre","323_Autre","313_Autre","223_Autre","219_Autre","115_Autre","357_Autre") ~ "Other",
    TRUE ~ as.character(clcm_lvl3)  # Pour gérer d'autres valeurs non spécifiées
  ))

#### CLC LVL 2 ####
LandWorm_dataset_site_m <- LandWorm_dataset_site_m %>%
  mutate(clc_lvl2 = case_when(
    clcm_lvl3 %in% c("214_Culture_annuelle","214c_Culture_annuelle") ~ "21_Arable land",
    clcm_lvl3 %in% c("218_Vignes_et_autres_Cultures_pérennes", "221_Arboriculture","221_Verger","312_Verger") ~ "22_Permanent crops",
    clcm_lvl3 %in% c("220ir_Agroforesterie", "220r_Agroforesterie") ~ "24_Heterogeneous agricultural areas",
    clcm_lvl3 %in% c("130_Prairie_humide","Prairie_agricole","210_Prairie_agricole_permanente", "211_Prairie_agricole_temporaire","210_Prairie_agricole_permanente?","121_Prairie_naturelle_&_Paturage", "131_Mégaphorbiaie","231_Prairies","Prairie agricole") ~ "23_Pastures",
    clcm_lvl3 %in% c("111_Forêt_de_feuillus","113_Forêt_mixte","310_Forêt_urbaine","112_Forêt_de_conifères") ~ "31_Forests",
    clcm_lvl3 %in% c("114c_Bois_ou_forêt_contaminé", "114_Bois","311_Bois_urbain","132_Lande_hygrophile","122_Lande_mésophile_et_broussailles","120_Végétation_clairsemée") ~ "32_Scrub and/or herbaceous vegetation associations",
    clcm_lvl3 %in% c("324_Prairie_aéroport","321_Pelouse","380_Friche_industrielle", "320_Prairie_urbaine","321_Pelouse_urbaine","391_Bassin_de_rétention", "330_Massif_ornemental","351_Massif_ornemental_contraint", "356_Noue_arborée","356_Noue_arborée","355_Noue_enherbée","352_Alignement_d'arbres_en_fosse_continue","350_Bande_enherbée_contrainte","353_Arbre_en_fosse_unitaire","354_Rond-point","322_Terrain_récréatif") ~ "14a_Artificial, non-agricultural vegetated areas",
    clcm_lvl3 %in% c("340_Potager","340p_Potager-pelouse","340p_Potager-Pelouse") ~ "14b_Urban agriculture",
    clcm_lvl3 %in% c("212_Jachère","215_Interculture","232_Haies","213_Bande_enherbée","216l_Légume_plein_champ","216_Maraîchage_et_légume_plein_champ","230_Talus_et_bordure_enherbés","217_Serre,_tunnel_fixe","232_Haie","216m_Maraichage") ~ "24_Heterogeneous agricultural areas",
    clcm_lvl3 %in% c("370_Stock_de_terre") ~ "13_Mine, dump and construction sites",
    clcm_lvl3 %in% c("124_Plage,_dune_et_sable") ~ "33_Open spaces with little or no vegetation",
    clcm_lvl3 %in% c("343_Autre","323_Autre","313_Autre","223_Autre","219_Autre","115_Autre","357_Autre") ~ "Other",
    TRUE ~ as.character(clcm_lvl3)  # Pour gérer d'autres valeurs non spécifiées
  ))

#### CLC LVL 3 ####
LandWorm_dataset_site_m <- LandWorm_dataset_site_m %>%
  mutate(clc_lvl1 = case_when(
    clcm_lvl3 %in% c("Prairie_agricole","223_Autre","219_Autre","214_Culture_annuelle","214c_Culture_annuelle","218_Vignes_et_autres_Cultures_pérennes", "221_Arboriculture","221_Verger","312_Verger","220ir_Agroforesterie", "220r_Agroforesterie", "210_Prairie_agricole_permanente", "211_Prairie_agricole_temporaire","210_Prairie_agricole_permanente?","121_Prairie_naturelle_&_Paturage", "131_Mégaphorbiaie","231_Prairies","Prairie agricole","212_Jachère","215_Interculture","232_Haies","213_Bande_enherbée","216l_Légume_plein_champ","216_Maraîchage_et_légume_plein_champ","230_Talus_et_bordure_enherbés","217_Serre,_tunnel_fixe","232_Haie","216m_Maraichage") ~ "2_Agricultural areas",
    clcm_lvl3 %in% c("115_Autre","124_Plage,_dune_et_sable","130_Prairie_humide","111_Forêt_de_feuillus","113_Forêt_mixte","310_Forêt_urbaine","112_Forêt_de_conifères","114c_Bois_ou_forêt_contaminé", "114_Bois","311_Bois_urbain","132_Lande_hygrophile","122_Lande_mésophile_et_broussailles","120_Végétation_clairsemée") ~ "3_Forest and semi natural areas",
    clcm_lvl3 %in% c("324_Prairie_aéroport","321_Pelouse","357_Autre","343_Autre","323_Autre","313_Autre","370_Stock_de_terre","380_Friche_industrielle", "320_Prairie_urbaine","321_Pelouse_urbaine","391_Bassin_de_rétention", "330_Massif_ornemental","351_Massif_ornemental_contraint", "356_Noue_arborée","356_Noue_arborée","355_Noue_enherbée","352_Alignement_d'arbres_en_fosse_continue","350_Bande_enherbée_contrainte","353_Arbre_en_fosse_unitaire","354_Rond-point","322_Terrain_récréatif","340_Potager","340p_Potager-pelouse","340p_Potager-Pelouse") ~ "1_Artificial surfaces",
    TRUE ~ as.character(clcm_lvl3)  # Pour gérer d'autres valeurs non spécifiées
  ))

#########################################
####### Soil texture correction ########
#########################################

#Calculation of total sand and silt from finer fractions
LandWorm_dataset_site_m <- LandWorm_dataset_site_m %>%
  mutate(
    sand = if_else(is.na(sand) & !is.na(fine_sand) & !is.na(coarse_sand),
                   fine_sand + coarse_sand, 
                   sand),
    silt = if_else(is.na(silt) & !is.na(fine_silt) & !is.na(coarse_silt),
                   fine_silt + coarse_silt, 
                   silt))


#Transformation of the unit to get %
LandWorm_dataset_site_m <- LandWorm_dataset_site_m %>%
  # First condition for dividing by 10 if the sum exceeds 900
  dplyr::mutate(total_1 = fine_sand + coarse_sand + fine_silt + coarse_silt + clay,
         across(c(fine_sand, coarse_sand, fine_silt, coarse_silt, clay, silt, sand),
                ~ if_else(total_1 > 900, . / 10, .))) %>%
  # Second condition for dividing by 10 if the sum exceeds 900
  dplyr::mutate(total_2 = sand + silt + clay,
         across(c(clay, silt, sand),
                ~ if_else(total_2 > 900, . / 10, .))) %>%
  # Delete temporary columns
  dplyr::select(-total_1, -total_2)





###############
##### STOP ####
###############

###########################################
#### GRAPHIQUES TEXTURE = A SUPPRIMER #####
# Données de texture du sol
texture_data <- na.omit(LandWorm_dataset_site_m[, c("sand", "silt", "clay")])%>%
  mutate(total = clay + silt + sand)%>%
  mutate(clay = (clay / total) * 100,
         silt = (silt / total) * 100,
         sand = (sand / total) * 100)%>%
  select(-total)%>%
  dplyr::rename(SAND= sand,
         SILT= silt,
         CLAY= clay)

# Tracer le triangle de texture sans les arguments 'xlab', 'ylab', et 'ID'
ggplot(     # 'Regular' ggplot 
  data=texture_data,         
  aes(
    x = SAND,
    y = CLAY,
    z = SILT                 
  )) +
  coord_tern(                # Add z coordinate to ggplot
    L='x',                   # Left pole
    T='y',                   # Top pole
    R='z'                    # Right pole
  )+
  geom_point()


textures_ggplot<-textures_ggplot+
  geom_polygon(
    data=USDA,aes(x=Sand,y=Clay,z=Silt,group=Label),
    fill=NA,size = 0.3,alpha=0.5,
    color = "grey30"
  )
