library(tidyverse)
# Specify the directory path containing the .csv files
directory_path <- "[Database]/[Raw datasets]/Plot_description_datasets"

# List the .csv files in the specified directory
csv_files <- list.files(directory_path, pattern = "*_description.csv", full.names = TRUE)

# Load and rename the .csv files into the environment
for (csv_file in csv_files) {
  table_name <- tools::file_path_sans_ext(basename(csv_file))  # Table name based on the file name
  
  # Read the .csv file without manually setting colClasses
  data <- read.csv(
    csv_file,
    sep = ",",
    stringsAsFactors = FALSE,
    na.strings = c("", "NA", "N/A", "null")  # Specify NA values that can be in the data
  )


  # Rename specific columns
  if ("Protocole" %in% colnames(data)) {
    colnames(data)[colnames(data) == "Protocole"] <- "Protocole_description"
  }
  if ("GPS_X" %in% colnames(data)) {
    colnames(data)[colnames(data) == "GPS_X"] <- "gps_x"
  }
  if ("GPS_Y" %in% colnames(data)) {
    colnames(data)[colnames(data) == "GPS_Y"] <- "gps_y"
  }
  
  if ("Code_Postal" %in% colnames(data)) {
    colnames(data)[colnames(data) == "Code_Postal"] <- "postal_code"
  }
  if ("Climat" %in% colnames(data)) {
    colnames(data)[colnames(data) == "Climat"] <- "Climate Zone"
  }
  if ("altitude" %in% colnames(data)) {
    colnames(data)[colnames(data) == "altitude"] <- "Altitude"
  }
  if ("Categorie_Milieu_Niv1" %in% colnames(data)) {
    colnames(data)[colnames(data) == "Categorie_Milieu_Niv1"] <- "clcm_lvl1"
  }
  if ("SousCategorie_Milieu_Niv2" %in% colnames(data)) {
    colnames(data)[colnames(data) == "SousCategorie_Milieu_Niv2"] <- "clcm_lvl2"
  }
  if ("Details_Milieu_Niv3" %in% colnames(data)) {
    colnames(data)[colnames(data) == "Details_Milieu_Niv3"] <- "clcm_lvl3"
  }
  if ("Specificites_Parcelle_Niv4" %in% colnames(data)) {
    colnames(data)[colnames(data) == "Specificites_Parcelle_Niv4"] <- "land_cover_detail"
  }
  if ("SableF" %in% colnames(data)) {
    colnames(data)[colnames(data) == "SableF"] <- "fine_sand"
  }
  if ("SableG" %in% colnames(data)) {
    colnames(data)[colnames(data) == "SableG"] <- "coarse_sand"
  }
  if ("LimonF" %in% colnames(data)) {
    colnames(data)[colnames(data) == "LimonF"] <- "fine_silt"
  }
  if ("LimonG" %in% colnames(data)) {
    colnames(data)[colnames(data) == "LimonG"] <- "coarse_silt"
  }
  if ("Argile" %in% colnames(data)) {
    colnames(data)[colnames(data) == "Argile"] <- "clay"
  }
  if ("Argiles..0.à.0.002.mm._g.kg" %in% colnames(data)) {
    colnames(data)[colnames(data) == "Argiles..0.à.0.002.mm._g.kg"] <- "clay"
  }
  if ("Argile....2.µm..g.kg" %in% colnames(data)) {
    colnames(data)[colnames(data) == "Argile....2.µm..g.kg"] <- "clay"
  }
  if ("Limon" %in% colnames(data)) {
    colnames(data)[colnames(data) == "Limon"] <- "silt"
  }
  if ("Sable" %in% colnames(data)) {
    colnames(data)[colnames(data) == "Sable"] <- "sand"
  }
  if ("Limons" %in% colnames(data)) {
    colnames(data)[colnames(data) == "Limons"] <- "silt"
  }
  if ("Limons.grossiers..20.à.50.µm..g.kg" %in% colnames(data)) {
    colnames(data)[colnames(data) == "Limons.grossiers..20.à.50.µm..g.kg"] <- "coarse_silt"
  }
  if ("Limons.grossiers..0.02.à.0.05.mm._g.kg" %in% colnames(data)) {
    colnames(data)[colnames(data) == "Limons.grossiers..0.02.à.0.05.mm._g.kg"] <- "coarse_silt"
  }
  if ("Limons.fins..0.002.à.0.02.mm._g.kg" %in% colnames(data)) {
    colnames(data)[colnames(data) == "Limons.fins..0.002.à.0.02.mm._g.kg"] <- "fine_silt"
  }
  if ("Limons.fins..2.à.20.µm..g.kg" %in% colnames(data)) {
    colnames(data)[colnames(data) == "Limons.fins..2.à.20.µm..g.kg"] <- "fine_silt"
  }
  if ("Sables" %in% colnames(data)) {
    colnames(data)[colnames(data) == "Sables"] <- "sand"
  }
  if ("Sables.grossiers..200.à.2000.µm..g.kg" %in% colnames(data)) {
    colnames(data)[colnames(data) == "Sables.grossiers..200.à.2000.µm..g.kg"] <- "coarse_sand"
  }
  if ("Sables.grossiers..0.2.à.2.0.mm._g.kg" %in% colnames(data)) {
    colnames(data)[colnames(data) == "Sables.grossiers..0.2.à.2.0.mm._g.kg"] <- "coarse_sand"
  }
  if ("Sables.fins..50.à.200.µm..g.kg" %in% colnames(data)) {
    colnames(data)[colnames(data) == "Sables.fins..50.à.200.µm..g.kg"] <- "fine_sand"
  }
  if ("Sables.fins..0.05.à.0.2.mm._g.kg" %in% colnames(data)) {
    colnames(data)[colnames(data) == "Sables.fins..0.05.à.0.2.mm._g.kg"] <- "fine_sand"
  }
  if ("pH_eau" %in% colnames(data)) {
    colnames(data)[colnames(data) == "pH_eau"] <- "ph_eau"
  }
  if ("pH.eau.pH.units.No.unit" %in% colnames(data)) {
    colnames(data)[colnames(data) == "pH.eau.pH.units.No.unit"] <- "ph_eau"
  }
  if ("pH" %in% colnames(data)) { ### Pour opvtBZH et TIGA
    colnames(data)[colnames(data) == "pH"] <- "ph_eau"
  }
  if ("pH_KCl" %in% colnames(data)) {
    colnames(data)[colnames(data) == "pH_KCl"] <- "ph_kcl"
  }
  if ("C_tot" %in% colnames(data)) {
    colnames(data)[colnames(data) == "C_tot"] <- "c_tot"
  }
  if ("Carbone..C..total_g.kg" %in% colnames(data)) {
    colnames(data)[colnames(data) == "Carbone..C..total_g.kg"] <- "c_tot"
  }
  if ("C" %in% colnames(data)) { 
    colnames(data)[colnames(data) == "C"] <- "c_tot"
  }
  if ("C_org" %in% colnames(data)) {
    colnames(data)[colnames(data) == "C_org"] <- "c_org"
  }
  if ("Carbone.organique.g.kg" %in% colnames(data)) {
    colnames(data)[colnames(data) == "Carbone.organique.g.kg"] <- "c_org"
  }
  if ("Carbone.organique.total..COT._g.kg" %in% colnames(data)) {
    colnames(data)[colnames(data) == "Carbone.organique.total..COT._g.kg"] <- "c_org"
  }
  if ("Azote..N..total_g.kg" %in% colnames(data)) { 
    colnames(data)[colnames(data) == "Azote..N..total_g.kg"] <- "n_tot"
  }
  if ("Azote.Dumas.g.kg" %in% colnames(data)) {
    colnames(data)[colnames(data) == "Azote.Dumas.g.kg"] <- "n_tot"
  }
  if ("N_tot" %in% colnames(data)) {
    colnames(data)[colnames(data) == "N_tot"] <- "n_tot"
  }
  if ("N" %in% colnames(data)) { ## sicarex meta
    colnames(data)[colnames(data) == "N"] <- "n_tot"
  }
  if ("C.N" %in% colnames(data)) {
    colnames(data)[colnames(data) == "C.N"] <- "c/n"
  }
  if ("C_N" %in% colnames(data)) {
    colnames(data)[colnames(data) == "C_N"] <- "c/n"
  }
  if ("MO_pourc" %in% colnames(data)) {
    colnames(data)[colnames(data) == "MO_pourc"] <- "om_perc"
  }
  if ("MO" %in% colnames(data)) {
    colnames(data)[colnames(data) == "MO"] <- "om"
  }
  if ("Matières.Organiques..Carbone.x.1.73..g.kg" %in% colnames(data)) {
    colnames(data)[colnames(data) == "Matières.Organiques..Carbone.x.1.73..g.kg"] <- "om"
  }
  if ("Matière.organique_g.kg" %in% colnames(data)) {
    colnames(data)[colnames(data) == "Matière.organique_g.kg"] <- "om"
  }
  if ("mo_tot" %in% colnames(data)) {
    colnames(data)[colnames(data) == "mo_tot"] <- "om"
  }
  if ("Cuivre" %in% colnames(data)) {
    colnames(data)[colnames(data) == "Cuivre"] <- "cu_tot"
  }
  if ("cu_tot" %in% colnames(data)) {
    colnames(data)[colnames(data) == "cu_tot"] <- "cu_tot"
  }
  if ("Cuivre..Cu..1" %in% colnames(data)) {
    colnames(data)[colnames(data) == "Cuivre..Cu..1"] <- "cu_tot"
  }
  if ("CU" %in% colnames(data)) {
    colnames(data)[colnames(data) == "Cuivre"] <- "cu_tot"
  }
  if ("Cuivre..Cu..mg.kg" %in% colnames(data)) {
    colnames(data)[colnames(data) == "Cuivre..Cu..mg.kg"] <- "cu_tot"
  }
  if ("Cuivre..Cu." %in% colnames(data)) {
    colnames(data)[colnames(data) == "Cuivre..Cu."] <- "cu_tot"
  }
  if ("Cuivre_total" %in% colnames(data)) {
    colnames(data)[colnames(data) == "Cuivre_total"] <- "cu_tot"
  }
  if ("Cuivre.extractible.EDTA..Cu..mg.kg" %in% colnames(data)) {
    colnames(data)[colnames(data) == "Cuivre.extractible.EDTA..Cu..mg.kg"] <- "cu_EDTA"
  }
  if ("Profondeur_Sol_cm" %in% colnames(data)) {
  colnames(data)[colnames(data) == "Profondeur_Sol_cm"] <- "soil_depth"
  }
  if ("Travail_sol" %in% colnames(data)) {
    colnames(data)[colnames(data) == "Travail_sol"] <- "type_tillage"
  }
  if ("gc_tech_trav_sol_2021" %in% colnames(data)) {
  colnames(data)[colnames(data) == "gc_tech_trav_sol_2021"] <- "type_tillage"
  }
  if ("Type_W_sol" %in% colnames(data)) {
  colnames(data)[colnames(data) == "Type_W_sol"] <- "type_tillage"
  }
  if ("Nbre_Interventions_W_sol" %in% colnames(data)) {
  colnames(data)[colnames(data) == "Nbre_Interventions_W_sol"] <- "tillage_frequency_intra"
  }
  if ("X15.ferti_min_produit...213" %in% colnames(data)) {
  colnames(data)[colnames(data) == "X15.ferti_min_produit...213"] <- "ferti_min_product"
  }
  if ("gc_ferti_min_2021" %in% colnames(data)) {
  colnames(data)[colnames(data) == "gc_ferti_min_2021"] <- "ferti_min_product"
  }
  if ("X15.ferti_min_qtte_kg.ha...214" %in% colnames(data)) {
  colnames(data)[colnames(data) == "X15.ferti_min_qtte_kg.ha...214"] <- "ferti_min_qtty"
  }
  if ("Qte_ferti_min_Unite_N_Efficace" %in% colnames(data)) {
  colnames(data)[colnames(data) == "Qte_ferti_min_Unite_N_Efficace"] <- "ferti_min_qtty"
  }
  if ("X15.ferti_orga_produit...216" %in% colnames(data)) {
  colnames(data)[colnames(data) == "X15.ferti_orga_produit...216"] <- "ferti_orga_product"
  }
  if ("X15.ferti_orga_qtte_kg.ha...217" %in% colnames(data)) {
  colnames(data)[colnames(data) == "X15.ferti_orga_qtte_kg.ha...217"] <- "ferti_orga_qtty"
  }
  if ("Qte_fert_org_Unite_N_Ttotal" %in% colnames(data)) {
  colnames(data)[colnames(data) == "Qte_fert_org_Unite_N_Ttotal"] <- "ferti_orga_qtty"
  }
  if ("X15.ferti_orga_passage.an...218" %in% colnames(data)) {
  colnames(data)[colnames(data) == "X15.ferti_orga_passage.an...218"] <- "ferti_orga_freq"
  }
  if ("gc_ferti_org_2021" %in% colnames(data)) {
  colnames(data)[colnames(data) == "gc_ferti_org_2021"] <- "ferti_orga_product"
  }
  if ("Fertilisation" %in% colnames(data)) {
  colnames(data)[colnames(data) == "Fertilisation"] <- "fertilisation"
  }
  if ("Nbre_total_herbicides" %in% colnames(data)) {
  colnames(data)[colnames(data) == "Nbre_total_herbicides"] <- "herbicide_freq"
  }
  if ("gc_nb_pass_her_2021" %in% colnames(data)) {
  colnames(data)[colnames(data) == "gc_nb_pass_her_2021"] <- "herbicide_freq"
  }
  if ("gc_nb_pass_fon_2021" %in% colnames(data)) {
  colnames(data)[colnames(data) == "gc_nb_pass_fon_2021"] <- "fungicide_freq"
  }
  if ("Nbre_fongicides" %in% colnames(data)) {
  colnames(data)[colnames(data) == "Nbre_fongicides"] <- "fungicide_freq"
  }
  if ("Nbre_insecticides" %in% colnames(data)) {
  colnames(data)[colnames(data) == "Nbre_insecticides"] <- "insecticide_freq"
  }
  if ("gc_nb_pass_ins_2021" %in% colnames(data)) {
  colnames(data)[colnames(data) == "gc_nb_pass_ins_2021"] <- "insecticide_freq"
  }
  if ("Nbre_molluscicides" %in% colnames(data)) {
  colnames(data)[colnames(data) == "Nbre_molluscicides"] <- "molluscicide_freq"
  }
  if ("IFT_Herbicide_12mois_Av_pv" %in% colnames(data)) {
  colnames(data)[colnames(data) == "IFT_Herbicide_12mois_Av_pv"] <- "tfi_herbicide"
  }
  if ("gc_qte_ha_her_2021" %in% colnames(data)) {
  colnames(data)[colnames(data) == "gc_qte_ha_her_2021"] <- "tfi_herbicide"
  }
  if ("gc_qte_ha_fon_2021" %in% colnames(data)) {
  colnames(data)[colnames(data) == "gc_qte_ha_fon_2021"] <- "tfi_fongicide"
  }
  if ("gc_qte_ha_ins_2021" %in% colnames(data)) {
  colnames(data)[colnames(data) == "gc_qte_ha_ins_2021"] <- "tfi_insecticide"
  }
  if ("IFT_total_12mois_Av_pv" %in% colnames(data)) {
  colnames(data)[colnames(data) == "IFT_total_12mois_Av_pv"] <- "total_tfi"
  }
  if ("Nbre_culture_dans_rotation" %in% colnames(data)) {
  colnames(data)[colnames(data) == "Nbre_culture_dans_rotation"] <- "rotation_plant_div"
  }
  if ("gc_residus_2021" %in% colnames(data)) {
  colnames(data)[colnames(data) == "gc_residus_2021"] <- "crop_residues_management"
  }
  if ("Gestion_Herbe" %in% colnames(data)) {
  colnames(data)[colnames(data) == "Gestion_Herbe"] <- "herbage_use"
  }
  if ("X5.age_prairie_ans" %in% colnames(data)) {
  colnames(data)[colnames(data) == "X5.age_prairie_ans"] <- "herb_age"
  }
  if ("Age" %in% colnames(data)) {
  colnames(data)[colnames(data) == "Age"] <- "herb_age"
  }
  if ("X7.type_prairie" %in% colnames(data)) {
  colnames(data)[colnames(data) == "X7.type_prairie"] <- "grassland_type"
  }
  if ("p_usage_parcelle" %in% colnames(data)) {
  colnames(data)[colnames(data) == "p_usage_parcelle"] <- "grassland_type"
  }
  if ("Type_de_prairie" %in% colnames(data)) {
  colnames(data)[colnames(data) == "Type_de_prairie"] <- "grassland_type"
  }
  if ("X10.paturage_animaux" %in% colnames(data)) {
  colnames(data)[colnames(data) == "X10.paturage_animaux"] <- "animal_loading"
  }
  if ("X19.Paturage_animaux" %in% colnames(data)) {
  colnames(data)[colnames(data) == "X19.Paturage_animaux"] <- "trampling_nature"
  }
  if ("X22.paturage_UGB" %in% colnames(data)) {
  colnames(data)[colnames(data) == "X22.paturage_UGB "] <- "animal_loading"
  }
  if ("UGB" %in% colnames(data)) {
  colnames(data)[colnames(data) == "UGB "] <- "animal_loading"
  }
  if ("pp_ugb_ha" %in% colnames(data)) {
  colnames(data)[colnames(data) == "pp_ugb_ha "] <- "animal_loading"
  }
  if ("Gestion_Herbe" %in% colnames(data)) {
  colnames(data)[colnames(data) == "Gestion_Herbe "] <- "herbage_use"
  }
  if ("pp_foin" %in% colnames(data)) {
  colnames(data)[colnames(data) == "pp_foin "] <- "mowing_frequency_yr"
  }



  
  
  # Replace spaces with underscores in clcm_lvl2 column
  if ("clcm_lvl2" %in% colnames(data)) {
    data$clcm_lvl2 <- gsub(" ", "_", data$clcm_lvl2)
    
    # Additional replacements
    data$clcm_lvl2 <- gsub("22_Agricole_Boisé", "22_Agricole_boisé", data$clcm_lvl2)
    data$clcm_lvl2 <- gsub("23_Bord_de_champ", "23_Bord_de_champs", data$clcm_lvl2)
    data$clcm_lvl2 <- gsub("31_Espace_vert_fermé", "31_Espace_vert_boisé", data$clcm_lvl2)
    data$clcm_lvl2 <- gsub("33_Espace_vert_ouverts", "32_Espace_vert_ouvert", data$clcm_lvl2)
    data$clcm_lvl2 <- gsub("32_Espace_vert_ouverts", "32_Espace_vert_ouvert", data$clcm_lvl2)
    data$clcm_lvl2 <- gsub("33_Espace_vert_ouverts", "32_Espace_vert_ouvert", data$clcm_lvl2)
    data$clcm_lvl2 <- gsub("32_Espaces_verts_ouverts", "32_Espace_vert_ouvert", data$clcm_lvl2)
    data$clcm_lvl2 <- gsub("33_Espace_vert_ouvert", "32_Espace_vert_ouvert", data$clcm_lvl2)
    
  }
  
  if ("clcm_lvl3" %in% colnames(data)) {
    data$clcm_lvl3 <- gsub(" ", "_", data$clcm_lvl3)
    
    data$clcm_lvl3 <- gsub('210_Prairie_agricole_permanente?', "210_Prairie_agricole_permanente", data$clcm_lvl3)
    data$clcm_lvl3 <- gsub('216_Maraichage_et_légumes_plein_champs', "216_Maraîchage_et_légume_plein_champ", data$clcm_lvl3)
    data$clcm_lvl3 <- gsub('216l_Légume_plein_champ', "216l_Légume_plein_champ", data$clcm_lvl3)
    data$clcm_lvl3 <- gsub('216_Légume_plein_champ', "216l_Légume_plein_champ", data$clcm_lvl3)
    data$clcm_lvl3 <- gsub("218_Culture_anuelle", "214_Culture_annuelle", data$clcm_lvl3)
    data$clcm_lvl3 <- gsub('218_Culture_perenne', "218_Vignes_et_autres_Cultures_pérennes", data$clcm_lvl3)
    data$clcm_lvl3 <- gsub('218_Culture_perenne_Vignes', "218_Vignes_et_autres_Cultures_pérennes", data$clcm_lvl3)
    data$clcm_lvl3 <- gsub('218_Cultures_pérennes', "218_Vignes_et_autres_Cultures_pérennes", data$clcm_lvl3)
    data$clcm_lvl3 <- gsub('218_Vignes_et_autres_Cultures_pérennes_Vignes', "218_Vignes_et_autres_Cultures_pérennes", data$clcm_lvl3)
    data$clcm_lvl3 <- gsub("Grasslands", "231_Prairies", data$clcm_lvl3)
    data$clcm_lvl3 <- gsub('210_Prairie_agricole_permanente?', "210_Prairie_agricole_permanente", data$clcm_lvl3)
    data$clcm_lvl3 <- gsub('111_Forêt_feuillus', "111_Forêt_de_feuillus", data$clcm_lvl3)
    data$clcm_lvl3 <- gsub('113_Foret_mixte', "113_Forêt_mixte", data$clcm_lvl3)
    data$clcm_lvl3 <- gsub('115_Autres', "115_Autre", data$clcm_lvl3)
    data$clcm_lvl3 <- gsub('121_Prairie_naturelle', "121_Prairie_naturelle_&_Paturage", data$clcm_lvl3)
    data$clcm_lvl3 <- gsub('121_Prairie_naturelle_&_Paturage_&_Paturage', "121_Prairie_naturelle_&_Paturage", data$clcm_lvl3)
    data$clcm_lvl3 <- gsub('320_Pelouse_urbaine', "321_Pelouse_urbaine", data$clcm_lvl3)
    data$clcm_lvl3 <- gsub('322_Pelouse_urbaine', "321_Pelouse_urbaine", data$clcm_lvl3)
    data$clcm_lvl3 <- gsub('320_Pelouse_urbaine', "321_Pelouse_urbaine", data$clcm_lvl3)
    
  }  
  
  # Assign the renamed data to the table_name
  assign(table_name, data)
}



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
    "postal_code",	"gps_x",	"gps_y",	"Altitude",	"Climate Zone",
    "clcm_lvl1",	"clcm_lvl2",	"clcm_lvl3", "land_cover_detail",
    "ph_kcl",	"ph_eau",	"c_tot",	"c_org",	"n_tot",	"c/n",	"om",	"cu_tot",	"cu_EDTA",	"soil_temperature",	"soil_humidity",	"fine_sand",	"coarse_sand",	"sand",	"fine_silt",	"coarse_silt",	"silt",	"clay",
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

