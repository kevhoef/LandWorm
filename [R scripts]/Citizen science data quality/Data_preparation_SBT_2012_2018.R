
#Chargement de l'environnement
#load("D:/Copie PC/Post-doc KH/LANDWORM/DATA/BDD R - NATHAN/Environnements/BDD_ECOBIO_v2022.11.01.Rdata")


#### PACKAGE LOADING ####
library(tidyverse)
library(forcats)
library(readxl)
library(readr)


#### DATASET LOADING ####
source("./[R scripts]/Citizen science data quality/Data_loading_SBT_2012_2018.R", echo=FALSE) 


#### DATASET CURATION/CLEANUP ####

#Data check
str(sbt_m)
unique(sbt_m$Protocole)
unique(sbt_m$Annee)
unique(sbt_m$CE.Deter.Terrain)
summary(sbt_m$CE.Deter.Terrain)


#Dataset cleanup
CS_error <- sbt_m %>%
  mutate(Protocole = fct_drop(Protocole, "MTM")) %>% # Focus on the mustard protocol
  rename(FIEC = CE.Deter.Terrain, Nbr_EW = Nbr_VDT, LIEC = GF) %>% 
  mutate_at(c("FIEC", "LIEC", "Stade", "ID_Site", "Annee", "Code_Taxon"), as.factor) %>%
  mutate(Nbr_EW = as.numeric(Nbr_EW)) %>%
  mutate(FIEC = fct_recode(FIEC, "ANE_BH" = "ATN", "ANE_RH" = "ATR")) %>%
  mutate(LIEC = fct_recode(LIEC, "ANE_BH" = "ATN", "ANE_RH" = "ATR")) %>%
  mutate(Region = substr(ID_Site, 1, 2)) %>%
  unite(FIEC_LIEC, FIEC, LIEC, sep = "_", remove = FALSE) %>%
  mutate_at(c("FIEC_LIEC", "Region"), as.factor) %>%
  filter(!is.na(FIEC)) # deletion of NA lines because no ecological category can be identified in the field

#Check
str(CS_error)
summary(CS_error$FIEC)
any(is.na(CS_error$Nbr_EW))==F # if true = OK
any(is.na(CS_error$LIEC))==F # if true = OK
any(is.na(CS_error$Region))==F # if true = OK
unique(CS_error$Region)
any(is.na(CS_error$FIEC_LIEC))==F # if true = OK

####VIEUX CODES A VERIFIER####
#nrow(CS_error)==20743 # A VERIFIER
#length(unique(paste(CS_error$Annee,CS_error$ID_Site)))==1074 # A VERIFIER
#View(CS_error[duplicated2(paste(CS_error$Annee,CS_error$ID_Site,CS_error$Repetition,CS_error$Code_Methode,CS_error$Code_Taxon,CS_error$Deter.Code.Taxon,CS_error$FIEC,CS_error$Stade)),])
#View(CS_error[duplicated2(CS_error),])
#nrow(CS_error[duplicated2(CS_error),])==36
# ya des tableaux ou on a litt?ralement 1 ligne = 1 vdt (1 seul individu par ligne => Nbr_EW==1 pour toutes les lignes)
#ex :
#SBT-ENI	M	M	2012	BO_09	BOURG-09-VI	NA	NA	NA	NA	2	EPI	LX_A	LT	EPA	JV

# ya aussi des doublons provenants de tableaux avec une difference dans l'ancienne 
#colonne Stockage, que j'ai remplac? par pilu de stockage dans l'uniformisation
#ex :
#SBT-ENI	M	2016	PC	PCB9	PC_12	M	1	EPI	LC	LC	LC	EPI	AD
#SBT-ENI	M	M	2016	PC_12	PCB9	NA	NA	NA	NA	1	EPI	LC	LC	EPI	AD

#lignes avec doublon ET il y avait pas de difference dans l'ancienne colonne "Stockage
#ex :
#SBT-ENI	M	2018	AQ		AQ-034	AQ_34	M	1	EPI	LX-E	LCD	LCD	LCD	EPI	JV
#SBT-ENI	M	M	2018	AQ_28	AQ-028	NA	NA	NA	NA	2	END	ARR	ARR	END	AD

#avec les verifs ci dessus je considere que c'est ok, je passe ? la suite






#### CALCULATION OF UR AND MR MORPHOGROUP INDICES ####

{

###### Concatenation of certain variables of interest ######  
CS_error <- CS_error %>%
  mutate(FIEC_LIEC_Stade = paste(FIEC_LIEC, Stade, sep = "_"),
         ID_Site_Annee = paste(ID_Site, Annee, sep = "_")) %>%
  mutate_at(c("FIEC_LIEC_Stade", "ID_Site_Annee"), as.factor)
  
str(CS_error)

#### Groupement des CE d?ter/terrain ####
CS_error_group = CS_error %>%  
  group_by(FIEC_LIEC,ID_Site_Annee,ID_Site,Annee)  %>% 
  summarize(Nbr_EW=sum(Nbr_EW) ) 

# Transposage de la colonne FIEC_LIEC en plusieurs colonnes
CS_error_group_pivot =CS_error_group %>% 
  pivot_wider(id_cols = c( "ID_Site_Annee", "ID_Site","Annee"),
              names_from = FIEC_LIEC, 
              values_from = Nbr_EW,
              values_fill = list(Nbr_EW = 0))

#Calcul AB TOT
CS_error_AB = CS_error %>%  
  group_by(ID_Site_Annee) %>% 
  summarize(Nbr_EW=sum(Nbr_EW) ) 

CS_error_group_pivot=left_join(CS_error_AB,CS_error_group_pivot, by = "ID_Site_Annee")
any(anti_join(CS_error_AB, CS_error_group_pivot, by = "ID_Site_Annee"))==F



# Calcul stade DEV
CS_error_group_pivot_SD <- CS_error %>%
  group_by(Stade, ID_Site_Annee) %>%
  summarize(Nbr_EW = sum(Nbr_EW)) %>%
  pivot_wider(names_from = Stade, values_from = Nbr_EW, values_fill = 0)

CS_error_group_pivot=left_join(CS_error_group_pivot_SD,CS_error_group_pivot, by = "ID_Site_Annee")
any(anti_join(CS_error_group_pivot_SD, CS_error_group_pivot, by = "ID_Site_Annee"))==F

# Création variable ratio ((AD+SA)/(JV)) x 100
CS_error_group_pivot = mutate(CS_error_group_pivot, ADJV = ((AD+SA)*100)/Nbr_EW)



# Calcul diversit? LIEC
CS_error_group_LIEC = CS_error %>%  
  group_by(LIEC,ID_Site_Annee)  %>% 
  summarize(Nbr_EW=sum(Nbr_EW) ) 
CS_error_group_pivot_LIEC =CS_error_group_LIEC %>% 
  pivot_wider(id_cols = c( "ID_Site_Annee"),
              names_from = LIEC, 
              values_from = Nbr_EW,
              values_fill = list(Nbr_EW = 0))


CS_error_group_pivot_LIEC = mutate(CS_error_group_pivot_LIEC, Div_EPI = ifelse (EPI > 0, 1, 0)) %>% 
  mutate(CS_error_group_pivot_LIEC, Div_EPA = ifelse (ANE_RH > 0, 1, 0)) %>% 
  mutate(CS_error_group_pivot_LIEC, Div_ANS = ifelse (ANE_BH > 0, 1, 0)) %>% 
  mutate(CS_error_group_pivot_LIEC, Div_END = ifelse (END > 0, 1, 0)) %>% 
  mutate(CS_error_group_pivot_LIEC, Div_LIEC = (Div_EPI + Div_EPA + Div_ANS + Div_END)) %>% 
  select(ID_Site_Annee, Div_LIEC)

CS_error_group_pivot=left_join(CS_error_group_pivot_LIEC,CS_error_group_pivot, by ="ID_Site_Annee")
any(anti_join(CS_error_group_pivot_LIEC, CS_error_group_pivot, by = "ID_Site_Annee"))==F


#### MR et UR at national scale : EC & Stade developpement #####

# Ecological categories
CS_error_national = CS_error %>%  
  group_by(FIEC_LIEC)  %>% 
  summarize(Nbr_EW=sum(Nbr_EW) ) 

CS_error_national = pivot_wider(CS_error_national,
                                names_from = FIEC_LIEC, 
                                values_from = Nbr_EW)

CS_error_national = mutate(CS_error_national, MR_EPI = ((EPI_EPI+EPI_ANE_RH+EPI_ANE_BH+EPI_ANE_BH)-EPI_EPI)/(EPI_EPI+EPI_ANE_RH+EPI_ANE_BH+EPI_ANE_BH)*100)
CS_error_national = mutate(CS_error_national, MR_ANE_RH = ((ANE_RH_EPI+ANE_RH_ANE_RH+ANE_RH_ANE_BH+ANE_RH_END)-ANE_RH_ANE_RH)/(ANE_RH_EPI+ANE_RH_ANE_RH+ANE_RH_ANE_BH+ANE_RH_END)*100)
CS_error_national = mutate(CS_error_national, MR_ANE_BH = ((ANE_BH_EPI+ANE_BH_ANE_RH+ANE_BH_ANE_BH+ANE_BH_END)-ANE_BH_ANE_BH)/(ANE_BH_EPI+ANE_BH_ANE_RH+ANE_BH_ANE_BH+ANE_BH_END)*100)
CS_error_national = mutate(CS_error_national, MR_END = ((END_EPI+END_ANE_RH+END_ANE_BH+END_END)-END_END)/(END_EPI+END_ANE_RH+END_ANE_BH+END_END)*100)

CS_error_national = mutate(CS_error_national, UR_EPI = ((EPI_EPI+ANE_RH_EPI+ANE_BH_EPI+END_EPI+GF_X_EPI)-EPI_EPI)/(EPI_EPI+ANE_RH_EPI+ANE_BH_EPI+END_EPI+GF_X_EPI)*100)
CS_error_national = mutate(CS_error_national, UR_ANE_RH = ((EPI_ANE_RH+ANE_RH_ANE_RH+ANE_BH_ANE_RH+END_ANE_RH+GF_X_ANE_RH)-ANE_RH_ANE_RH)/(EPI_ANE_RH+ANE_RH_ANE_RH+ANE_BH_ANE_RH+END_ANE_RH+GF_X_ANE_RH)*100)
CS_error_national = mutate(CS_error_national, UR_ANE_BH = ((EPI_ANE_BH+ANE_RH_ANE_BH+ANE_BH_ANE_BH+END_ANE_BH+GF_X_ANE_BH)-ANE_BH_ANE_BH)/(EPI_ANE_BH+ANE_RH_ANE_BH+ANE_BH_ANE_BH+END_ANE_BH+GF_X_ANE_BH)*100)
CS_error_national = mutate(CS_error_national, UR_END = ((EPI_END+ANE_RH_END+ANE_BH_END+END_END+GF_X_END)-END_END)/(EPI_END+ANE_RH_END+ANE_BH_END+END_END+GF_X_END)*100)

CS_error_national <- CS_error_national %>% dplyr:: select(grep("MR", names(CS_error_national)), grep("UR", names(CS_error_national)))

#write.table(sbt_colerror_national, "D:/Home/khoeffner/Downloads/sbt_colerror_national.csv", sep=",", col.names=TRUE, row.names=TRUE, quote=TRUE, na="NA")
#view(sbt_colerror_national)


# Stage of development


CS_error_national_SD = CS_error %>%  
  group_by(FIEC_LIEC, Stade)  %>% 
  summarize(Nbr_EW=sum(Nbr_EW) ) 

CS_error_national_SD = pivot_wider(CS_error_national_SD,
                                    id_cols = c( "Stade"),
                                    names_from = FIEC_LIEC, 
                                    values_from = Nbr_EW)

CS_error_national_SD = mutate(CS_error_national_SD, MR_EPI = ((EPI_EPI+EPI_ANE_RH+EPI_ANE_BH+EPI_ANE_BH)-EPI_EPI)/(EPI_EPI+EPI_ANE_RH+EPI_ANE_BH+EPI_ANE_BH)*100)
CS_error_national_SD = mutate(CS_error_national_SD, MR_ANE_RH = ((ANE_RH_EPI+ANE_RH_ANE_RH+ANE_RH_ANE_BH+ANE_RH_END)-ANE_RH_ANE_RH)/(ANE_RH_EPI+ANE_RH_ANE_RH+ANE_RH_ANE_BH+ANE_RH_END)*100)
CS_error_national_SD = mutate(CS_error_national_SD, MR_ANE_BH = ((ANE_BH_EPI+ANE_BH_ANE_RH+ANE_BH_ANE_BH+ANE_BH_END)-ANE_BH_ANE_BH)/(ANE_BH_EPI+ANE_BH_ANE_RH+ANE_BH_ANE_BH+ANE_BH_END)*100)
CS_error_national_SD = mutate(CS_error_national_SD, MR_END = ((END_EPI+END_ANE_RH+END_ANE_BH+END_END)-END_END)/(END_EPI+END_ANE_RH+END_ANE_BH+END_END)*100)

CS_error_national_SD = mutate(CS_error_national_SD, UR_EPI = ((EPI_EPI+ANE_RH_EPI+ANE_BH_EPI+END_EPI+GF_X_EPI)-EPI_EPI)/(EPI_EPI+ANE_RH_EPI+ANE_BH_EPI+END_EPI+GF_X_EPI)*100)
CS_error_national_SD = mutate(CS_error_national_SD, UR_ANE_RH = ((EPI_ANE_RH+ANE_RH_ANE_RH+ANE_BH_ANE_RH+END_ANE_RH+GF_X_ANE_RH)-ANE_RH_ANE_RH)/(EPI_ANE_RH+ANE_RH_ANE_RH+ANE_BH_ANE_RH+END_ANE_RH+GF_X_ANE_RH)*100)
CS_error_national_SD = mutate(CS_error_national_SD, UR_ANE_BH = ((EPI_ANE_BH+ANE_RH_ANE_BH+ANE_BH_ANE_BH+END_ANE_BH+GF_X_ANE_BH)-ANE_BH_ANE_BH)/(EPI_ANE_BH+ANE_RH_ANE_BH+ANE_BH_ANE_BH+END_ANE_BH+GF_X_ANE_BH)*100)
CS_error_national_SD = mutate(CS_error_national_SD, UR_END = ((EPI_END+ANE_RH_END+ANE_BH_END+END_END+GF_X_END)-END_END)/(EPI_END+ANE_RH_END+ANE_BH_END+END_END+GF_X_END)*100)

CS_error_national_SD <- CS_error_national_SD %>% dplyr:: select(grep("MR", names(CS_error_national_SD)), grep("UR", names(CS_error_national_SD)), Stade)
#write.table(CS_error_national_SD, "D:/Home/khoeffner/Downloads/CS_error_national_SD.csv", sep=",", col.names=TRUE, row.names=TRUE, quote=TRUE, na="NA")
#view(CS_error_national_SD)




#### MR et UR at national scale : EC & Stade developpement #####
CS_error_plot = mutate(CS_error_group_pivot, MR_EPI = ((EPI_EPI+EPI_ANE_RH+EPI_ANE_BH+EPI_ANE_BH)-EPI_EPI)/(EPI_EPI+EPI_ANE_RH+EPI_ANE_BH+EPI_ANE_BH))
CS_error_plot = mutate(CS_error_plot, MR_ANE_RH = ((ANE_RH_EPI+ANE_RH_ANE_RH+ANE_RH_ANE_BH+ANE_RH_END)-ANE_RH_ANE_RH)/(ANE_RH_EPI+ANE_RH_ANE_RH+ANE_RH_ANE_BH+ANE_RH_END))
CS_error_plot = mutate(CS_error_plot, MR_ANE_BH = ((ANE_BH_EPI+ANE_BH_ANE_RH+ANE_BH_ANE_BH+ANE_BH_END)-ANE_BH_ANE_BH)/(ANE_BH_EPI+ANE_BH_ANE_RH+ANE_BH_ANE_BH+ANE_BH_END))
CS_error_plot = mutate(CS_error_plot, MR_END = ((END_EPI+END_ANE_RH+END_ANE_BH+END_END)-END_END)/(END_EPI+END_ANE_RH+END_ANE_BH+END_END))

CS_error_plot = mutate(CS_error_plot, UR_EPI = ((EPI_EPI+ANE_RH_EPI+ANE_BH_EPI+END_EPI+GF_X_EPI)-EPI_EPI)/(EPI_EPI+ANE_RH_EPI+ANE_BH_EPI+END_EPI+GF_X_EPI))
CS_error_plot = mutate(CS_error_plot, UR_ANE_RH = ((EPI_ANE_RH+ANE_RH_ANE_RH+ANE_BH_ANE_RH+END_ANE_RH+GF_X_ANE_RH)-ANE_RH_ANE_RH)/(EPI_ANE_RH+ANE_RH_ANE_RH+ANE_BH_ANE_RH+END_ANE_RH+GF_X_ANE_RH))
CS_error_plot = mutate(CS_error_plot, UR_ANE_BH = ((EPI_ANE_BH+ANE_RH_ANE_BH+ANE_BH_ANE_BH+END_ANE_BH+GF_X_ANE_BH)-ANE_BH_ANE_BH)/(EPI_ANE_BH+ANE_RH_ANE_BH+ANE_BH_ANE_BH+END_ANE_BH+GF_X_ANE_BH))
CS_error_plot = mutate(CS_error_plot, UR_END = ((EPI_END+ANE_RH_END+ANE_BH_END+END_END+GF_X_END)-END_END)/(EPI_END+ANE_RH_END+ANE_BH_END+END_END+GF_X_END))

}



#### ADDING DESCRIPTIVE VARIABLES ####
Codes_MNHN_UR=read_xlsx("./[Database]/Raw datasets/National/SBT_ENI/codes_MNHN_SBT-ENI_V09.06.2023.xlsx")
str(sbt_descr)

#Fusion des ID_site avec table VDT
CS_error_codes=left_join(CS_error_plot, Codes_MNHN_UR, by = c("ID_Site"= "ID_site"))

#Vérification des parcelles non jointes ?
any(anti_join(CS_error_codes, Codes_MNHN_UR, by = c("ID_Site"= "ID_site")))==F
#write.table(CS_error_codes, "D:/Home/khoeffner/Downloads/CS_error_codes.csv", sep=",", col.names=TRUE, row.names=TRUE, quote=TRUE, na="NA")


#Concatenation de ID_site avec annee pour avoir un ID unique
CS_error_codes_ID <-unite(CS_error_codes, ID_parcelle_Annee, id_parcelle,Annee, sep = "_", remove = FALSE)
#write.table(CS_error_codes_ID, "D:/Home/khoeffner/Downloads/CS_error_codes_ID.csv", sep=",", col.names=TRUE, row.names=TRUE, quote=TRUE, na="NA")



#Chargement des m?ta_compatage 
data_MNHN_earthworm = read_csv2(file ="./[Database]/Raw datasets/National/SBT_ENI/Export_biovigilance_MNHN_15_06_23/obs_comptage.csv")%>%
  subset(protocole == "Vers de terre - MOUTARDE")%>%
  unite(ID_parcelle_Annee, id_parcelle,annee, sep = "_", remove = FALSE)
#write.table(CS_error_eni_comptage_VDT, "D:/Home/khoeffner/Downloads/sbt_comptage.csv", sep=",", col.names=TRUE, row.names=TRUE, quote=TRUE, na="NA")
#View(duplicated(paste(CS_error_eni_comptage_VDT$ID_parcelle_Annee)))

#Fusion des SBT_fusion_error avec sbt_eni_comptage
CS_error_MNHN_earthworm=right_join(CS_error_codes_ID,data_MNHN_earthworm, by ="ID_parcelle_Annee")
NoID =anti_join(CS_error_codes_ID, data_MNHN_earthworm, by = "ID_parcelle_Annee")

# Sélection des colonnes du dataframe de gauche qui ne sont pas présentes dans le dataframe de droite
extra_cols <- setdiff(names(CS_error_codes_ID), names(data_MNHN_earthworm))
NoID <- NoID[, c(extra_cols, "ID_parcelle_Annee")]
#write.table(NoID, "D:/Home/khoeffner/Downloads/NoID.csv", sep=",", col.names=TRUE, row.names=TRUE, quote=TRUE, na="NA")

# Affichage du message
if (nrow(NoID) > 0) {
  message("Il y a des lignes non jointes.")
} else {
  message("Il n'y a pas de lignes non jointes.")
}
"====> Il y a 48 parcelles_annee qui n'ont pas été remplies dans la BDD MNHN"
write.table(NoID, "D:/Home/khoeffner/Downloads/NOID.csv", sep=",", col.names=TRUE, row.names=TRUE, quote=TRUE, na="NA")
write.table(data_MNHN_earthworm, "D:/Home/khoeffner/Downloads/data_MNHN_earthworm.csv", sep=",", col.names=TRUE, row.names=TRUE, quote=TRUE, na="NA")


#SBT_fusion_error_comptage$scds_field = as.numeric(difftime(SBT_fusion_error_comptage$heure_fin,SBT_fusion_error_comptage$heure_debut, units = "secs"))
#SBT_fusion_error_comptage = mutate(SBT_fusion_error_comptage, scds_field=ifelse(scds_field ==0 | scds_field < 0, NA, scds_field))


#Chargement des m?ta_obs
sbt_eni_VDT = read_csv2(file ="D:/Copie PC/Post-doc KH/LANDWORM/DATA/@EcoBioSoil_v25.02.2021 - NATHAN/Data_EcoBioSoil/FR/SBT_ENI/Export MNHN/2022/Export_biovigilance/obs_VDT.csv")
sbt_eni_VDT_aggreg = sbt_eni_VDT %>%  
  group_by(id_obs)  %>% 
  summarize(across(c(abondance_placette_1, abondance_placette_2,abondance_placette_3), sum))%>% 
  mutate(AB_field = abondance_placette_1 +abondance_placette_2+abondance_placette_3)
#write.table(SBT_fusion_error_comptage_obs, "D:/Home/khoeffner/Downloads/SBT_fusion_error_comptage_obs.csv", sep=",", col.names=TRUE, row.names=TRUE, quote=TRUE, na="NA")

SBT_fusion_error_comptage_obs=left_join(SBT_fusion_error_comptage, sbt_eni_VDT_aggreg, by ="id_obs")
SBT_error = SBT_fusion_error_comptage_obs %>%
  select(ID_Site,id_obs,Annee, MR_EPI, MR_ANE_RH, MR_ANE_BH, MR_END, UR_EPI, UR_ANE_RH, UR_ANE_BH, UR_END,Nbr_EW,AB_field,AD,SA,JV,ADJV, Div_LIEC, observateur,nuages,pluie,vent,temperature_air)%>%
  mutate_at(c("ID_Site","nuages","observateur","Annee"), as.factor) %>%
  drop_na("ID_Site","ADJV", "nuages", "Div_LIEC", "Nbr_EW","observateur") %>%
  mutate(observateur_annee = paste(observateur, Annee, sep = "_"))

###### Calcul de l'expérience des observateurs
Train_obs = SBT_error %>%  
  group_by(observateur, Annee) %>% 
  summarize(Nbr_EW= sum(Nbr_EW, na.rm=TRUE))  %>%
  drop_na(Annee)


Train_obs <- Train_obs %>% 
  pivot_wider(id_cols = c("observateur"),
              names_from = Annee, 
              values_from = Nbr_EW)


########

Train_obs <- Train_obs %>%
  mutate(across(starts_with("20"), ~ifelse(is.na(.), 0,.)))

EW_count <- Train_obs %>%
  mutate(EW_2012 = ifelse(`2012`== 0, 0, `2012`))%>%
  mutate(EW_2013 = ifelse(`2013`== 0, 0, `2012` + `2013`)) %>%
  mutate(EW_2014 = ifelse(`2014`== 0, 0, `2012`+ `2013`+ `2014`)) %>%
  mutate(EW_2015 = ifelse(`2015`== 0, 0, `2012`+ `2013`+ `2014` + `2015`)) %>%
  mutate(EW_2016 = ifelse(`2016`== 0, 0, `2012`+ `2013`+ `2014` + `2015` + `2016`)) %>%
  mutate(EW_2017 = ifelse(`2017`== 0, 0, `2012`+ `2013`+ `2014` + `2015` + `2016` + `2017`)) %>%
  mutate(EW_2018 = ifelse(`2018`== 0, 0, `2012`+ `2013`+ `2014` + `2015` + `2016` + `2017` + `2018`))


EW_count_transposed <- EW_count %>%
  pivot_longer(cols = starts_with("EW_"), 
               names_to = "EW_year",
               values_to = "EW",
               values_drop_na = TRUE) %>%
  mutate(Annee = sub("^EW_", "", EW_year)) %>%
  group_by(observateur) %>%
  mutate(observateur_annee = paste(observateur, Annee, sep = "_")) %>%
  ungroup() %>%
  select(observateur_annee, EW)


#########

Train_obs <- Train_obs %>%
  mutate(across(starts_with("20"), ~ifelse(is.na(.), 0, 1)))

view(EW_count)
Train_obs <- Train_obs %>%
  mutate(count_2012 = ifelse(`2012`== 0, NA, `2012`))%>%
  mutate(count_2013 = ifelse(`2013`== 0, NA, `2012` + `2013`)) %>%
  mutate(count_2014 = ifelse(`2014`== 0, NA, `2012`+ `2013`+ `2014`)) %>%
  mutate(count_2015 = ifelse(`2015`== 0, NA, `2012`+ `2013`+ `2014` + `2015`)) %>%
  mutate(count_2016 = ifelse(`2016`== 0, NA, `2012`+ `2013`+ `2014` + `2015` + `2016`)) %>%
  mutate(count_2017 = ifelse(`2017`== 0, NA, `2012`+ `2013`+ `2014` + `2015` + `2016` + `2017`)) %>%
  mutate(count_2018 = ifelse(`2018`== 0, NA, `2012`+ `2013`+ `2014` + `2015` + `2016` + `2017` + `2018`))

Train_obs_transposed <- Train_obs %>%
  pivot_longer(cols = starts_with("count_"), 
               names_to = "Count_year",
               values_to = "Train",
               values_drop_na = TRUE) %>%
  mutate(Annee = sub("^count_", "", Count_year)) %>%
  group_by(observateur) %>%
  mutate(observateur_annee = paste(observateur, Annee, sep = "_")) %>%
  ungroup() %>%
  select(observateur_annee, Train)

  

if (any(duplicated(Train_obs_transposed$observateur_annee))) {
  print("Il y a des doublons dans la colonne.")
} else {
  print("Il n'y a pas de doublons dans la colonne.")
}

############ Calcul du nombre de parcelle observée par observateur
Parcelle_obs <- SBT_error %>%  
  group_by(observateur, Annee) %>% 
  summarize(Parcelle = n_distinct(ID_Site), .groups = "drop")  %>%
  drop_na(Annee) %>%
  mutate(observateur_annee = paste(observateur, Annee, sep = "_"))

Parcelle_obs <- Parcelle_obs %>% 
  pivot_wider(id_cols = c("observateur"),
              names_from = Annee, 
              values_from = Parcelle)


Parcelle_obs <- Parcelle_obs %>%
  mutate(across(starts_with("20"), ~ifelse(is.na(.), 0,.)))

Parcelle_count <- Parcelle_obs %>%
  mutate(Parcelle_2012 = ifelse(`2012`== 0, 0, `2012`))%>%
  mutate(Parcelle_2013 = ifelse(`2013`== 0, 0, `2012` + `2013`)) %>%
  mutate(Parcelle_2014 = ifelse(`2014`== 0, 0, `2012`+ `2013`+ `2014`)) %>%
  mutate(Parcelle_2015 = ifelse(`2015`== 0, 0, `2012`+ `2013`+ `2014` + `2015`)) %>%
  mutate(Parcelle_2016 = ifelse(`2016`== 0, 0, `2012`+ `2013`+ `2014` + `2015` + `2016`)) %>%
  mutate(Parcelle_2017 = ifelse(`2017`== 0, 0, `2012`+ `2013`+ `2014` + `2015` + `2016` + `2017`)) %>%
  mutate(Parcelle_2018 = ifelse(`2018`== 0, 0, `2012`+ `2013`+ `2014` + `2015` + `2016` + `2017` + `2018`))

Parcelle_count_transposed <- Parcelle_count %>%
  pivot_longer(cols = starts_with("Parcelle_"), 
               names_to = "Parcelle_year",
               values_to = "Parcelle",
               values_drop_na = TRUE) %>%
  mutate(Annee = sub("^Parcelle_", "", Parcelle_year)) %>%
  group_by(observateur) %>%
  mutate(observateur_annee = paste(observateur, Annee, sep = "_")) %>%
  ungroup() %>%
  select(observateur_annee, Parcelle)

if (any(duplicated(Parcelle_count_transposed$observateur_annee))) {
  print("Il y a des doublons dans la colonne.")
} else {
  print("Il n'y a pas de doublons dans la colonne.")
}


############ Fichier final
SBT_error_train <- right_join(Train_obs_transposed, SBT_error, by = "observateur_annee")
SBT_error_train <- right_join(Parcelle_count_transposed, SBT_error_train, by = "observateur_annee")
SBT_error_train <- right_join(EW_count_transposed, SBT_error_train, by = "observateur_annee")


SBT_error_train =drop_na(SBT_error_train,"ID_Site","ADJV", "nuages", "Div_LIEC", "Nbr_EW", "observateur_annee", "Train")
non_jointes <- anti_join(Train_obs_transposed, SBT_error, by = "observateur_annee")

# Calculer la moyenne et les quantiles de la colonne numérique
mean_value <- mean(SBT_error_train$Nbr_EW)
lower_quantile <- quantile(SBT_error_train$Nbr_EW, 0.025)
upper_quantile <- quantile(SBT_error_train$Nbr_EW, 0.975)

# Filtrer les lignes basées sur les quantiles
SBT_error_train <- subset(SBT_error_train, Nbr_EW >= lower_quantile & Nbr_EW <= upper_quantile)




############ Regroupement pour une ligne par Observateur_Annee
SBT_error_obs_annee=SBT_error_train %>% 
  group_by(observateur_annee, Train, Parcelle) %>% 
  summarize(MR_EPI=mean(MR_EPI, na.rm= T),
            MR_ANE_RH = mean(MR_ANE_RH, na.rm = TRUE),
            MR_ANE_BH = mean(MR_ANE_BH, na.rm = TRUE),
            MR_END = mean(MR_END, na.rm = TRUE),
            UR_EPI=mean(UR_EPI, na.rm= T),
            UR_ANE_RH = mean(UR_ANE_RH, na.rm = TRUE),
            UR_ANE_BH = mean(UR_ANE_BH, na.rm = TRUE),
            UR_END = mean(UR_END, na.rm = TRUE),
            EW = mean(EW, na.rm = TRUE))%>%
  separate(observateur_annee, into = c("observateur", "annee"), sep = "_")%>%
  mutate_at(c("observateur", "annee"), as.factor) %>%
  mutate(EWTrain= EW/Train) %>%
  mutate(TrainParcelle= Train/Parcelle)




SBT_error_obs_annee=SBT_error_train %>% 
  group_by(observateur_annee, Train, Parcelle) %>% 
  summarize(MR_EPI=sd(MR_EPI, na.rm= T),
            MR_ANE_RH = sd(MR_ANE_RH, na.rm = TRUE),
            MR_ANE_BH = sd(MR_ANE_BH, na.rm = TRUE),
            MR_END = sd(MR_END, na.rm = TRUE),
            UR_EPI=sd(UR_EPI, na.rm= T),
            UR_ANE_RH = sd(UR_ANE_RH, na.rm = TRUE),
            UR_ANE_BH = sd(UR_ANE_BH, na.rm = TRUE),
            UR_END = sd(UR_END, na.rm = TRUE))%>%
  separate(observateur_annee, into = c("observateur", "annee"), sep = "_")

############ Regroupement pour une ligne par Observateur
SBT_error_obs=SBT_error_train %>% 
  group_by(observateur) %>% 
  summarize(MR_EPI=mean(MR_EPI, na.rm= T),
            MR_ANE_RH = mean(MR_ANE_RH, na.rm = TRUE),
            MR_ANE_BH = mean(MR_ANE_BH, na.rm = TRUE),
            MR_END = mean(MR_END, na.rm = TRUE),
            UR_EPI=mean(UR_EPI, na.rm= T),
            UR_ANE_RH = mean(UR_ANE_RH, na.rm = TRUE),
            UR_ANE_BH = mean(UR_ANE_BH, na.rm = TRUE),
            UR_END = mean(UR_END, na.rm = TRUE),
            Parcelle = sum(Parcelle, na.rm = TRUE),
            Train = mean(Train, na.rm = TRUE))%>%
  mutate(TrainParcelle= Train/Parcelle)


SBT_error_obs=SBT_error_train %>% 
  group_by(observateur) %>% 
  summarize(MR_EPI=sd(MR_EPI, na.rm= T),
            MR_ANE_RH = sd(MR_ANE_RH, na.rm = TRUE),
            MR_ANE_BH = sd(MR_ANE_BH, na.rm = TRUE),
            MR_END = sd(MR_END, na.rm = TRUE),
            UR_EPI=sd(UR_EPI, na.rm= T),
            UR_ANE_RH = sd(UR_ANE_RH, na.rm = TRUE),
            UR_ANE_BH = sd(UR_ANE_BH, na.rm = TRUE),
            UR_END = sd(UR_END, na.rm = TRUE),
            Parcelle = sum(Parcelle, na.rm = TRUE),
            Train = sd(Train, na.rm = TRUE))




Train <- SBT_error_train %>% 
  pivot_wider(id_cols = c("observateur"),
              names_from = Annee, 
              values_from = MR_END) %>%
  select(observateur, across('2012':'2018'))







Train_MR_END_transposed<- na.omit(Train_MR_END_transposed[, c("observateur", "Value", "Count")])



##############################################





SBT_error_MREPI = drop_na(SBT_error, "ID_Site","ADJV", "MR_EPI", "nuages", "Div_LIEC", "Nbr_EW")
SBT_error_MRANE_RH = drop_na(SBT_error, "ID_Site","ADJV", "MR_ANE_RH", "nuages", "Div_LIEC", "Nbr_EW")
SBT_error_MRANE_BH = drop_na(SBT_error_train, "ID_Site","ADJV", "MR_ANE_BH", "nuages", "Div_LIEC", "Nbr_EW")
SBT_error_MREND = drop_na(SBT_error, "ID_Site","ADJV", "MR_END", "nuages", "Div_LIEC", "Nbr_EW")
SBT_error_UREPI = drop_na(SBT_error, "ID_Site","ADJV", "UR_EPI", "nuages", "Div_LIEC", "Nbr_EW")
SBT_error_URANE_RH = drop_na(SBT_error, "ID_Site","ADJV", "UR_ANE_RH", "nuages", "Div_LIEC", "Nbr_EW")
SBT_error_URANE_BH = drop_na(SBT_error, "ID_Site","ADJV", "UR_ANE_BH", "nuages", "Div_LIEC", "Nbr_EW")
SBT_error_UREND = drop_na(SBT_error, "ID_Site","ADJV", "UR_END", "nuages", "Div_LIEC", "Nbr_EW")



SBT_error_UR = drop_na(SBT_error, "EC", "ID_Site","ADJV", "UR", "nuages", "Div_LIEC", "Nbr_EW")

# Calculer la moyenne et les quantiles de la colonne numérique
mean_value <- mean(SBT_error_MR$Nbr_EW)
lower_quantile <- quantile(SBT_error_MR$Nbr_EW, 0.025)
upper_quantile <- quantile(SBT_error_MR$Nbr_EW, 0.975)

# Filtrer les lignes basées sur les quantiles
SBT_error_MR <- subset(SBT_error_MR, Nbr_EW >= lower_quantile & Nbr_EW <= upper_quantile)

# Calculer la moyenne et les quantiles de la colonne numérique
mean_value <- mean(SBT_error_UR$Nbr_EW)
lower_quantile <- quantile(SBT_error_MR$Nbr_UR, 0.025)
upper_quantile <- quantile(SBT_error_MR$Nbr_EW, 0.975)

# Filtrer les lignes basées sur les quantiles
SBT_error_MR <- subset(SBT_error_MR, Nbr_EW >= lower_quantile & Nbr_EW <= upper_quantile)



##### TRANSFO POUR RECUPERER LES OBSERVATEUR * ANNEE ####
SBT_error_train = SBT_error %>%  
  group_by(observateur,Annee)  %>% 
  summarize(MR_EPI= mean(MR_EPI, na.rm=TRUE),MR_ANE_RH = mean(MR_ANE_RH, na.rm=TRUE), MR_ANE_BH = mean(MR_ANE_BH, na.rm=TRUE), MR_END=mean(MR_END,na.rm=TRUE),
            UR_EPI= mean(UR_EPI, na.rm=TRUE),UR_ANE_RH = mean(UR_ANE_RH, na.rm=TRUE), UR_ANE_BH = mean(UR_ANE_BH, na.rm=TRUE), UR_END=mean(UR_END,na.rm=TRUE)) %>%
  mutate_at(c("Annee"), as.factor)
  
view(SBT_error_train)

##### TRANSFO TABLEAU OBSERVATEUR * ANNEE ####
Train_MR_END <- SBT_error_train %>% 
  pivot_wider(id_cols = c("observateur"),
              names_from = Annee, 
              values_from = MR_END)
  select(observateur, across('2012':'2018'))

Train_MR_END <- Train_MR_EPI %>% 
  relocate(observateur, `2012`, `2013`, `2014`, `2015`, `2016`, `2017`, `2018`)
library(dplyr)

Train_MR_END_count <- Train_MR_END %>%
  mutate(count_2013 = rowSums(!is.na(.[, -2])),
         count_2014 = rowSums(!is.na(.[, -c(2:3)])),
         count_2015 = rowSums(!is.na(.[, -c(2:4)])),
         count_2016 = rowSums(!is.na(.[, -c(2:5)])),
         count_2017 = rowSums(!is.na(.[, -c(2:6)])),
         count_2018 = rowSums(!is.na(.[, -c(2:7)])))



library(dplyr)

Train_MR_EPI_count <- Train_MR_EPI %>%
  mutate(count_2013 = rowSums(!is.na(select(., starts_with("2012", exclude = "2013")))),
         count_2014 = rowSums(!is.na(select(., starts_with("2012", exclude = c("2013", "2014"))))),
         count_2015 = rowSums(!is.na(select(., starts_with("2012", exclude = c("2013", "2014", "2015"))))),
         count_2016 = rowSums(!is.na(select(., starts_with("2012", exclude = c("2013", "2014", "2015", "2016"))))),
         count_2017 = rowSums(!is.na(select(., starts_with("2012", exclude = c("2013", "2014", "2015", "2016", "2017"))))),
         count_2018 = rowSums(!is.na(select(., starts_with("2012", exclude = c("2013", "2014", "2015", "2016", "2017", "2018"))))))

Train_MR_EPI_count <- Train_MR_EPI %>%
  mutate(
    count_2012 = ifelse(!is.na(`2012`), 1, 0),
    count_2013 = ifelse(!is.na(`2012`) & !is.na(`2013`), 2, ifelse(!is.na(`2013`), 1, 0)),
    count_2014 = ifelse(!is.na(`2012`) & !is.na(`2013`) & !is.na(`2014`), 3, ifelse(!is.na(`2013`) & !is.na(`2014`), 2, ifelse(!is.na(`2014`), 1, 0))),
    count_2015 = ifelse(!is.na(`2012`) & !is.na(`2013`) & !is.na(`2014`) & !is.na(`2015`), 4, ifelse(!is.na(`2013`) & !is.na(`2014`) & !is.na(`2015`), 3, ifelse(!is.na(`2014`) & !is.na(`2015`), 2, ifelse(!is.na(`2015`), 1, 0)))),
    count_2016 = ifelse(!is.na(`2012`) & !is.na(`2013`) & !is.na(`2014`) & !is.na(`2015`) & !is.na(`2016`), 5, ifelse(!is.na(`2013`) & !is.na(`2014`) & !is.na(`2015`) & !is.na(`2016`), 4, ifelse(!is.na(`2014`) & !is.na(`2015`) & !is.na(`2016`), 3, ifelse(!is.na(`2015`) & !is.na(`2016`), 2, ifelse(!is.na(`2016`), 1, 0))))),
    count_2017 = ifelse(!is.na(`2012`) & !is.na(`2013`) & !is.na(`2014`) & !is.na(`2015`) & !is.na(`2016`) & !is.na(`2017`), 6, ifelse(!is.na(`2013`) & !is.na(`2014`) & !is.na(`2015`) & !is.na(`2016`) & !is.na(`2017`), 5, ifelse(!is.na(`2014`) & !is.na(`2015`) & !is.na(`2016`) & !is.na(`2017`), 4, ifelse(!is.na(`2015`) & !is.na(`2016`) & !is.na(`2017`), 3, ifelse(!is.na(`2016`) & !is.na(`2017`), 2, ifelse(!is.na(`2017`), 1, 0)))))),
    count_2018 = ifelse(!is.na(`2012`) & !is.na(`2013`) & !is.na(`2014`) & !is.na(`2015`) & !is.na(`2016`) & !is.na(`2017`) & !is.na(`2018`), 7, ifelse(!is.na(`2013`) & !is.na(`2014`) & !is.na(`2015`) & !is.na(`2016`) & !is.na(`2017`) & !is.na(`2018`), 6, ifelse(!is.na(`2014`) & !is.na(`2015`) & !is.na(`2016`) & !is.na(`2017`) & !is.na(`2018`), 5, ifelse(!is.na(`2015`) & !is.na(`2016`) & !is.na(`2017`) & !is.na(`2018`), 4, ifelse(!is.na(`2016`) & !is.na(`2017`) & !is.na(`2018`), 3, ifelse(!is.na(`2017`) & !is.na(`2018`), 2, ifelse(!is.na(`2018`), 1, 0)))))))
  )




Train_MR_END_transposed <- Train_MR_EPI %>%
  pivot_longer(cols = starts_with(c("2012", "2013", "2014", "2015", "2016", "2017", "2018")), 
               names_to = "MR_EPI", 
               values_to = "Value") %>%
  pivot_longer(cols = starts_with(c("count_2012", "count_2013", "count_2014", "count_2015", "count_2016", "count_2017", "count_2018")), 
               names_to = "MR_EPI_count", 
               values_to = "Count")
Train_MR_END_transposed<- na.omit(Train_MR_END_transposed[, c("observateur", "Value", "Count")])



#MR_EPI ####

##### TRANSFO TABLEAU OBSERVATEUR * ANNEE ####
Train_MR_EPI = SBT_error_train %>% 
  pivot_wider(id_cols = c("observateur"),
              names_from = Annee, 
              values_from = MR_EPI) 
Train_MR_EPI <- Train_MR_EPI %>%
  mutate(Occurrence_MR_EPI = rowSums(across(matches("^20\\d{2}$")) > 0, na.rm = TRUE)) %>%
  rowwise() %>%
  mutate(MR_EPI_mean = mean(c_across(matches("^20\\d{2}$")), na.rm = TRUE) / Occurrence_MR_EPI)

Train_MR_EPI<- na.omit(Train_MR_EPI[, c("observateur", "Occurrence_MR_EPI", "MR_EPI_mean")])


# Créer le graphique avec un nuage de points
theme_set(theme_minimal(base_size = 14))  # Choix d'un thème minimaliste avec une taille de police de base de 14

# Calcul du coefficient de détermination et de la significativité
lm_model <- lm(MR_EPI_mean ~ Occurrence_MR_EPI, data = Train_MR_EPI)
r_squared <- summary(lm_model)$r.squared
p_value <- summary(lm_model)$coefficients[2, 4]
significance <- ifelse(p_value < 0.05, "*", " ")

# Création du graphique avec des personnalisations
ggplot(Train_MR_EPI, aes(x = Occurrence_MR_EPI, y = MR_EPI_mean)) +
  geom_point(color = "#FF7F00", size = 4) +  # Couleur orange pour les points, taille de 4
  geom_smooth(method = "lm", se = FALSE, color = "#1F77B4", linetype = "dashed") +  # Couleur bleue pour la droite de régression, style en pointillé
  labs(x = "Occurrences", y = "MR_EPI") +
  ggtitle("Relation entre Occurences et MR_EPI") +
  geom_text(x = max(Train_MR_EPI$Occurrence_MR_EPI), y = max(Train_MR_EPI$MR_EPI_mean), 
            label = paste("R² =", round(r_squared, 2), significance),
            hjust = 1, vjust = 1, color = "#FF0000", size = 5, fontface = "bold") +  # Couleur rouge, taille de 5, texte en gras pour le R²
  theme(plot.title = element_text(size = 18, face = "bold", margin = margin(b = 20)),  # Personnalisation du titre du graphique
        axis.title = element_text(size = 16, face = "bold"),  # Personnalisation des titres des axes
        axis.text = element_text(size = 14),  # Personnalisation des étiquettes des axes
        legend.position = "none")  # Suppression de la légEPIe


#MR_EPA ####

##### TRANSFO TABLEAU OBSERVATEUR * ANNEE ####
Train_MR_ANE_RH = SBT_error_train %>% 
  pivot_wider(names_from = Annee, 
              values_from = MR_ANE_RH) 
Train_MR_ANE_RH <- Train_MR_ANE_RH %>%
  mutate(Occurrence_MR_ANE_RH = rowSums(across(matches("^20\\d{2}$")) > 0, na.rm = TRUE)) %>%
  rowwise() %>%
  mutate(MR_ANE_RH_mean = mean(c_across(matches("^20\\d{2}$")), na.rm = TRUE) / Occurrence_MR_ANE_RH)

Train_MR_ANE_RH<- na.omit(Train_MR_ANE_RH[, c("observateur", "Occurrence_MR_ANE_RH", "MR_ANE_RH_mean")])


# Créer le graphique avec un nuage de points
theme_set(theme_minimal(base_size = 14))  # Choix d'un thème minimaliste avec une taille de police de base de 14

# Calcul du coefficient de détermination et de la significativité
lm_model <- lm(MR_ANE_RH_mean ~ Occurrence_MR_ANE_RH, data = Train_MR_ANE_RH)
r_squared <- summary(lm_model)$r.squared
p_value <- summary(lm_model)$coefficients[2, 4]
significance <- ifelse(p_value < 0.05, "*", " ")

# Création du graphique avec des personnalisations
ggplot(Train_MR_ANE_RH, aes(x = Occurrence_MR_ANE_RH, y = MR_ANE_RH_mean)) +
  geom_point(color = "#FF7F00", size = 4) +  # Couleur orange pour les points, taille de 4
  geom_smooth(method = "lm", se = FALSE, color = "#1F77B4", linetype = "dashed") +  # Couleur bleue pour la droite de régression, style en pointillé
  labs(x = "Occurrences", y = "MR_ANE_RH") +
  ggtitle("Relation entre Occurences et MR_ANE_RH") +
  geom_text(x = max(Train_MR_ANE_RH$Occurrence_MR_ANE_RH), y = max(Train_MR_ANE_RH$MR_ANE_RH_mean), 
            label = paste("R² =", round(r_squared, 2), significance),
            hjust = 1, vjust = 1, color = "#FF0000", size = 5, fontface = "bold") +  # Couleur rouge, taille de 5, texte en gras pour le R²
  theme(plot.title = element_text(size = 18, face = "bold", margin = margin(b = 20)),  # Personnalisation du titre du graphique
        axis.title = element_text(size = 16, face = "bold"),  # Personnalisation des titres des axes
        axis.text = element_text(size = 14),  # Personnalisation des étiquettes des axes
        legend.position = "none")  # Suppression de la légANE_RHe


#MR_EPA ####

##### TRANSFO TABLEAU OBSERVATEUR * ANNEE ####
Train_MR_ANE_BH = SBT_error_train %>% 
  pivot_wider(id_cols = c("observateur"),
              names_from = Annee, 
              values_from = MR_ANE_BH) 
Train_MR_ANE_BH <- Train_MR_ANE_BH %>%
  mutate(Occurrence_MR_ANE_BH = rowSums(across(matches("^20\\d{2}$")) > 0, na.rm = TRUE)) %>%
  rowwise() %>%
  mutate(MR_ANE_BH_mean = mean(c_across(matches("^20\\d{2}$")), na.rm = TRUE) / Occurrence_MR_ANE_BH)

Train_MR_ANE_BH<- na.omit(Train_MR_ANE_BH[, c("observateur", "Occurrence_MR_ANE_BH", "MR_ANE_BH_mean")])


# Créer le graphique avec un nuage de points
theme_set(theme_minimal(base_size = 14))  # Choix d'un thème minimaliste avec une taille de police de base de 14

# Calcul du coefficient de détermination et de la significativité
lm_model <- lm(MR_ANE_BH_mean ~ Occurrence_MR_ANE_BH, data = Train_MR_ANE_BH)
r_squared <- summary(lm_model)$r.squared
p_value <- summary(lm_model)$coefficients[2, 4]
significance <- ifelse(p_value < 0.05, "*", " ")

# Création du graphique avec des personnalisations
library(ggplot2)

# Personnalisation des options esthétiques
theme_set(theme_minimal(base_size = 14))  # Choix d'un thème minimaliste avec une taille de police de base de 14

# Création du graphique avec des personnalisations
ggplot(Train_MR_ANE_BH, aes(x = Occurrence_MR_ANE_BH, y = MR_ANE_BH_mean)) +
  geom_point(size = 4) +  # Taille de 4 pour les points
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +  # Style en pointillé pour la droite de régression
  labs(x = "Occurrences", y = "MR_END") +
  ggtitle("Relation entre Occurences et MR_END") +
  geom_text(x = max(Train_MR_ANE_BH$Occurrence_MR_ANE_BH), y = max(Train_MR_ANE_BH$MR_ANE_BH_mean), 
            label = paste("R² =", round(r_squared, 2), significance),
            hjust = 1, vjust = 1, color = "#FF0000", size = 5, fontface = "bold") +  # Couleur rouge, taille de 5, texte en gras pour le R²
  theme(plot.title = element_text(size = 18, face = "bold", margin = margin(b = 20)),  # Personnalisation du titre du graphique
        axis.title = element_text(size = 16, face = "bold"),  # Personnalisation des titres des axes
        axis.text = element_text(size = 14),  # Personnalisation des étiquettes des axes
        legend.position = "none")  # Position de la légende à droite


library(ggplot2)

# Personnalisation des options esthétiques
theme_set(theme_minimal(base_size = 14))  # Choix d'un thème minimaliste avec une taille de police de base de 14

# Création du graphique avec des personnalisations
ggplot(Train_MR_ANE_BH, aes(x = Occurrence_MR_ANE_BH, y = MR_ANE_BH_mean, color = observateur, shape = observateur)) +
  geom_point(size = 4) +  # Taille de 4 pour les points
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +  # Style en pointillé pour la droite de régression
  labs(x = "Occurrences", y = "MR_ANE_BH") +
  ggtitle("Relation entre Occurences et MR_ANE_BH") +
  geom_text(x = max(Train_MR_ANE_BH$Occurrence_MR_ANE_BH), y = max(Train_MR_ANE_BH$MR_ANE_BH_mean), 
            label = paste("R² =", round(r_squared, 2), significance),
            hjust = 1, vjust = 1, color = "#FF0000", size = 5, fontface = "bold") +  # Couleur rouge, taille de 5, texte en gras pour le R²
  scale_shape_manual(values = 1:length(unique(Train_MR_ANE_BH$observateur))) +  # Définition des formes manuellement
  theme(plot.title = element_text(size = 18, face = "bold", margin = margin(b = 20)),  # Personnalisation du titre du graphique
        axis.title = element_text(size = 16, face = "bold"),  # Personnalisation des titres des axes
        axis.text = element_text(size = 14),  # Personnalisation des étiquettes des axes
        legend.position = "none")  # Position de la légANE_BHe à droite











#### MR_END ####
Train_MR_END <- Train_MR_END %>%
  mutate(Occurrence_MR_END = rowSums(across(matches("^20\\d{2}$")) > 0, na.rm = TRUE)) %>%
  rowwise() %>%
  mutate(MR_END_mean = mean(c_across(matches("^20\\d{2}$")), na.rm = TRUE) / Occurrence_MR_END)

Train_MR_END<- na.omit(Train_MR_END[, c("observateur", "Occurrence_MR_END", "MR_END_mean")])


# Créer le graphique avec un nuage de points
theme_set(theme_minimal(base_size = 14))  # Choix d'un thème minimaliste avec une taille de police de base de 14

# Personnalisation des options esthétiques
theme_set(theme_minimal(base_size = 14))  # Choix d'un thème minimaliste avec une taille de police de base de 14

# Calcul du coefficient de détermination et de la significativité
lm_model <- lm(MR_END_mean ~ Occurrence_MR_END, data = Train_MR_END)
r_squared <- summary(lm_model)$r.squared
p_value <- summary(lm_model)$coefficients[2, 4]
significance <- ifelse(p_value < 0.05, "*", " ")

# Création du graphique avec des personnalisations
ggplot(Train_MR_END, aes(x = Occurrence_MR_END, y = MR_END_mean)) +
  geom_point(color = "#FF7F00", size = 4) +  # Couleur orange pour les points, taille de 4
  geom_smooth(method = "lm", se = FALSE, color = "#1F77B4", linetype = "dashed") +  # Couleur bleue pour la droite de régression, style en pointillé
  labs(x = "Occurrences", y = "MR_END") +
  ggtitle("Relation entre Occurences et MR_END") +
  geom_text(x = max(Train_MR_END$Occurrence_MR_END), y = max(Train_MR_END$MR_END_mean), 
            label = paste("R² =", round(r_squared, 2), significance),
            hjust = 1, vjust = 1, color = "#FF0000", size = 5, fontface = "bold") +  # Couleur rouge, taille de 5, texte en gras pour le R²
  theme(plot.title = element_text(size = 18, face = "bold", margin = margin(b = 20)),  # Personnalisation du titre du graphique
        axis.title = element_text(size = 16, face = "bold"),  # Personnalisation des titres des axes
        axis.text = element_text(size = 14),  # Personnalisation des étiquettes des axes
        legend.position = "none")  # Suppression de la légende



  
  group_by(observateur,Annee)  %>% 
  summarize(MR_EPI= mean(MR_EPI, na.rm=TRUE),MR_ANE_RH = mean(MR_ANE_RH, na.rm=TRUE), MR_ANE_BH = mean(MR_ANE_BH, na.rm=TRUE), MR_END=mean(MR_END,na.rm=TRUE),
            UR_EPI= mean(UR_EPI, na.rm=TRUE),UR_ANE_RH = mean(UR_ANE_RH, na.rm=TRUE), UR_ANE_BH = mean(UR_ANE_BH, na.rm=TRUE), UR_END=mean(UR_END,na.rm=TRUE)) %>%
  mutate_at(c("Annee"), as.factor)




















##################################
##################################
##################################
##################################
fviz_mca_ind(res.mca)
fviz_mca_var(res.mca)
fviz_mca_biplot(res.mca)
fviz_screeplot(res.mca, addlabels = TRUE, ylim = c(0, 45))

fviz_mca_ind(res.mca, 
             label = "none", # hide individual labels
             habillage = "Vomiting", # color by groups 
             palette = c("#00AFBB", "#E7B800"),
             addEllipses = TRUE, ellipse.type = "confidence",
             ggtheme = theme_minimal()) 
fviz_mca_biplot(res.mca)

print(res.mca)
SBT_TEST = count(SBT_fusion_error_comptage, Nbr_EW)



SBT_fusion_error_comptage$FIEC_LIEC=as.factor(SBT_fusion_error_comptage$FIEC_LIEC)

table(SBT_fusion_error_comptage$FIEC_LIEC, SBT_fusion_error_comptage$Study_Site) 

table(SBT_error$ID_Site, SBT_error$Annee) 

SBT_TEST4= spread(SBT_fusion_error_comptage, key=FIEC_LIEC, value=Nbr_EW, drop = FALSE)


SBT_test <- xtabs( Nbr_EW ~ ID_site_annee +FIEC_LIEC  , data=SBT_fusion_error_comptage,  na.action=na.pass)
SBT_test= data.frame(SBT_test)


SBT_TEST_2 = SBT_fusion_error_comptage %>%  group_by(FIEC_LIEC) %>%  summarize(Nbr_EW= sum(Nbr_EW, na.rm=TRUE))

SBT_TEST2 <- SBT_fusion_error_comptage %>%  spread(FIEC_LIEC)

SBT_test=SBT_fusion_error_comptage %>% group_by(FIEC_LIEC, ID_site_annee) %>% summarise(Nbr_EW = sum())

SBT_test %>% spread(FIEC_LIEC,ID_site_annee) %>% mutate(Total=rowSums(.)) %>% rbind(.,Total=colSums(.))
SBT_test_2= data.frame(SBT_test_2)


SBT_test <- xtabs(Freq ~ ID_site_annee +FIEC_LIEC  , data=SBT_test,  na.action=na.pass)
?spread

SBT_test_3 =xtabs(~FIEC_LIEC+Freq, data=SBT_test_2)
SBT_test_3= data.frame(SBT_test_3)

view(surveys)
SBT_test_2 <- xtabs( Nbr_EW ~ ID_site_annee +FIEC_LIEC  , data=SBT_fusion_error_comptage,  na.action=na.pass)






SBT_fusion_error_comptage %>% spread(FIEC_LIEC,ID_site_annee) %>% mutate(Total=rowSums(.)) %>% rbind(.,Total=colSums(.))


SBT_test <- xtabs( Nbr_EW ~ ID_site_annee +FIEC_LIEC  , data=SBT_fusion_error_comptage,  na.action=na.pass)

SBT_test<-table(SBT_fusion_error_comptage$FIEC_LIEC,SBT_fusion_error_comptage$Nbr_EW)
print(SBT_test)
SBT_test= data.frame( matrix(SBT_test ) )






##### nombre de prelevement par annee et region #####





colID=c("Programme","Protocole","Code_Methode","ID_Site","Code_Parcelle","Annee",
        "Modalite","Bloc")#"Date_Prelevement", dans les metas => pas grave
colID=intersect(colID,colnames(CS_error))

CS_error$cle=as.vector(apply(CS_error[,colID],1,paste,collapse="|"))
length(unique(CS_error$cle))==1074#ok


dataID <- unique(CS_error[,c("cle",colID,"Region")])
nrow(dataID)==1074

table(dataID$Region,dataID$Annee)
#tout comme pr?vu














##### Region bien formee : National #####

nrow(Tab2)==21044
dataID <- unique(CS_error[,c("cle",colID,"Region")])
nrow(dataID)==1109

n2012=table(dataID$Annee)[1]
n2013=table(dataID$Annee)[2]
n2014=table(dataID$Annee)[3]
n2015=table(dataID$Annee)[4]
n2016=table(dataID$Annee)[5]
n2017=table(dataID$Annee)[6]
n2018=table(dataID$Annee)[7]

vectReussiteNAT=c(
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC%in%c("EPI_EPI","EPA_EPA","ANS_ANS","END_END") & Tab2$Annee%in%2013])/sum(Tab2$Nbr_EW[Tab2$Annee%in%2013]),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC%in%c("EPI_EPI","EPA_EPA","ANS_ANS","END_END") & Tab2$Annee%in%2014])/sum(Tab2$Nbr_EW[Tab2$Annee%in%2014]),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC%in%c("EPI_EPI","EPA_EPA","ANS_ANS","END_END") & Tab2$Annee%in%2015])/sum(Tab2$Nbr_EW[Tab2$Annee%in%2015]),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC%in%c("EPI_EPI","EPA_EPA","ANS_ANS","END_END") & Tab2$Annee%in%2016])/sum(Tab2$Nbr_EW[Tab2$Annee%in%2016]),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC%in%c("EPI_EPI","EPA_EPA","ANS_ANS","END_END") & Tab2$Annee%in%2017])/sum(Tab2$Nbr_EW[Tab2$Annee%in%2017]),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC%in%c("EPI_EPI","EPA_EPA","ANS_ANS","END_END") & Tab2$Annee%in%2018])/sum(Tab2$Nbr_EW[Tab2$Annee%in%2018]))


##### taux pour les _EPI #####


# calcule des sommes

som_EPItt=sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPI_EPI"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPA_EPI"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="ANS_EPI"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="END_EPI"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="LIEC_X_EPI"])

som_EPI13=sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPI_EPI" & Tab2$Annee=="2013"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPA_EPI" & Tab2$Annee=="2013"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="ANS_EPI" & Tab2$Annee=="2013"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="END_EPI" & Tab2$Annee=="2013"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="LIEC_X_EPI" & Tab2$Annee=="2013"])

som_EPI14=sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPI_EPI" & Tab2$Annee=="2014"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPA_EPI" & Tab2$Annee=="2014"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="ANS_EPI" & Tab2$Annee=="2014"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="END_EPI" & Tab2$Annee=="2014"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="LIEC_X_EPI" & Tab2$Annee=="2014"])

som_EPI15=sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPI_EPI" & Tab2$Annee=="2015"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPA_EPI" & Tab2$Annee=="2015"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="ANS_EPI" & Tab2$Annee=="2015"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="END_EPI" & Tab2$Annee=="2015"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="LIEC_X_EPI" & Tab2$Annee=="2015"])

som_EPI16=sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPI_EPI" & Tab2$Annee=="2016"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPA_EPI" & Tab2$Annee=="2016"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="ANS_EPI" & Tab2$Annee=="2016"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="END_EPI" & Tab2$Annee=="2016"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="LIEC_X_EPI" & Tab2$Annee=="2016"])

som_EPI17=sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPI_EPI" & Tab2$Annee=="2017"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPA_EPI" & Tab2$Annee=="2017"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="ANS_EPI" & Tab2$Annee=="2017"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="END_EPI" & Tab2$Annee=="2017"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="LIEC_X_EPI" & Tab2$Annee=="2017"])

som_EPI18=sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPI_EPI" & Tab2$Annee=="2018"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPA_EPI" & Tab2$Annee=="2018"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="ANS_EPI" & Tab2$Annee=="2018"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="END_EPI" & Tab2$Annee=="2018"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="LIEC_X_EPI" & Tab2$Annee=="2018"])




# vecteurs pour chaque annee

vectEPI_EPI=c(
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPI_EPI" & Tab2$Annee=="2013"])/
    (som_EPI13),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPI_EPI" & Tab2$Annee=="2014"])/
    (som_EPI14),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPI_EPI" & Tab2$Annee=="2015"])/
    (som_EPI15),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPI_EPI" & Tab2$Annee=="2016"])/
    (som_EPI16),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPI_EPI" & Tab2$Annee=="2017"])/
    (som_EPI17),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPI_EPI" & Tab2$Annee=="2018"])/
    (som_EPI18))

vectEPA_EPI=c(
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPA_EPI" & Tab2$Annee=="2013"])/
    (som_EPI13),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPA_EPI" & Tab2$Annee=="2014"])/
    (som_EPI14),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPA_EPI" & Tab2$Annee=="2015"])/
    (som_EPI15),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPA_EPI" & Tab2$Annee=="2016"])/
    (som_EPI16),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPA_EPI" & Tab2$Annee=="2017"])/
    (som_EPI17),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPA_EPI" & Tab2$Annee=="2018"])/
    (som_EPI18))

vectANS_EPI=c(
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="ANS_EPI" & Tab2$Annee=="2013"])/
    (som_EPI13),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="ANS_EPI" & Tab2$Annee=="2014"])/
    (som_EPI14),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="ANS_EPI" & Tab2$Annee=="2015"])/
    (som_EPI15),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="ANS_EPI" & Tab2$Annee=="2016"])/
    (som_EPI16),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="ANS_EPI" & Tab2$Annee=="2017"])/
    (som_EPI17),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="ANS_EPI" & Tab2$Annee=="2018"])/
    (som_EPI18))

vectEND_EPI=c(
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="END_EPI" & Tab2$Annee=="2013"])/
    (som_EPI13),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="END_EPI" & Tab2$Annee=="2014"])/
    (som_EPI14),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="END_EPI" & Tab2$Annee=="2015"])/
    (som_EPI15),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="END_EPI" & Tab2$Annee=="2016"])/
    (som_EPI16),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="END_EPI" & Tab2$Annee=="2017"])/
    (som_EPI17),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="END_EPI" & Tab2$Annee=="2018"])/
    (som_EPI18))

vectLIEC_X_EPI=c(
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="LIEC_X_EPI" & Tab2$Annee=="2013"])/
    (som_EPI13),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="LIEC_X_EPI" & Tab2$Annee=="2014"])/
    (som_EPI14),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="LIEC_X_EPI" & Tab2$Annee=="2015"])/
    (som_EPI15),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="LIEC_X_EPI" & Tab2$Annee=="2016"])/
    (som_EPI16),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="LIEC_X_EPI" & Tab2$Annee=="2017"])/
    (som_EPI17),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="LIEC_X_EPI" & Tab2$Annee=="2018"])/
    (som_EPI18))




##### taux pour les _EPA #####


# calcule des sommes

som_EPAtt=sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPI_EPA"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPA_EPA"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="ANS_EPA"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="END_EPA"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="LIEC_X_EPA"])

som_EPA13=sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPI_EPA" & Tab2$Annee=="2013"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPA_EPA" & Tab2$Annee=="2013"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="ANS_EPA" & Tab2$Annee=="2013"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="END_EPA" & Tab2$Annee=="2013"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="LIEC_X_EPA" & Tab2$Annee=="2013"])

som_EPA14=sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPI_EPA" & Tab2$Annee=="2014"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPA_EPA" & Tab2$Annee=="2014"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="ANS_EPA" & Tab2$Annee=="2014"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="END_EPA" & Tab2$Annee=="2014"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="LIEC_X_EPA" & Tab2$Annee=="2014"])

som_EPA15=sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPI_EPA" & Tab2$Annee=="2015"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPA_EPA" & Tab2$Annee=="2015"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="ANS_EPA" & Tab2$Annee=="2015"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="END_EPA" & Tab2$Annee=="2015"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="LIEC_X_EPA" & Tab2$Annee=="2015"])

som_EPA16=sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPI_EPA" & Tab2$Annee=="2016"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPA_EPA" & Tab2$Annee=="2016"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="ANS_EPA" & Tab2$Annee=="2016"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="END_EPA" & Tab2$Annee=="2016"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="LIEC_X_EPA" & Tab2$Annee=="2016"])

som_EPA17=sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPI_EPA" & Tab2$Annee=="2017"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPA_EPA" & Tab2$Annee=="2017"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="ANS_EPA" & Tab2$Annee=="2017"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="END_EPA" & Tab2$Annee=="2017"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="LIEC_X_EPA" & Tab2$Annee=="2017"])

som_EPA18=sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPI_EPA" & Tab2$Annee=="2018"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPA_EPA" & Tab2$Annee=="2018"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="ANS_EPA" & Tab2$Annee=="2018"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="END_EPA" & Tab2$Annee=="2018"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="LIEC_X_EPA" & Tab2$Annee=="2018"])






# vecteurs pour chaque annee

vectEPI_EPA=c(
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPI_EPA" & Tab2$Annee=="2013"])/
    (som_EPA13),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPI_EPA" & Tab2$Annee=="2014"])/
    (som_EPA14),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPI_EPA" & Tab2$Annee=="2015"])/
    (som_EPA15),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPI_EPA" & Tab2$Annee=="2016"])/
    (som_EPA16),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPI_EPA" & Tab2$Annee=="2017"])/
    (som_EPA17),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPI_EPA" & Tab2$Annee=="2018"])/
    (som_EPA18))

vectEPA_EPA=c(
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPA_EPA" & Tab2$Annee=="2013"])/
    (som_EPA13),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPA_EPA" & Tab2$Annee=="2014"])/
    (som_EPA14),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPA_EPA" & Tab2$Annee=="2015"])/
    (som_EPA15),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPA_EPA" & Tab2$Annee=="2016"])/
    (som_EPA16),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPA_EPA" & Tab2$Annee=="2017"])/
    (som_EPA17),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPA_EPA" & Tab2$Annee=="2018"])/
    (som_EPA18))

vectANS_EPA=c(
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="ANS_EPA" & Tab2$Annee=="2013"])/
    (som_EPA13),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="ANS_EPA" & Tab2$Annee=="2014"])/
    (som_EPA14),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="ANS_EPA" & Tab2$Annee=="2015"])/
    (som_EPA15),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="ANS_EPA" & Tab2$Annee=="2016"])/
    (som_EPA16),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="ANS_EPA" & Tab2$Annee=="2017"])/
    (som_EPA17),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="ANS_EPA" & Tab2$Annee=="2018"])/
    (som_EPA18))

vectEND_EPA=c(
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="END_EPA" & Tab2$Annee=="2013"])/
    (som_EPA13),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="END_EPA" & Tab2$Annee=="2014"])/
    (som_EPA14),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="END_EPA" & Tab2$Annee=="2015"])/
    (som_EPA15),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="END_EPA" & Tab2$Annee=="2016"])/
    (som_EPA16),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="END_EPA" & Tab2$Annee=="2017"])/
    (som_EPA17),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="END_EPA" & Tab2$Annee=="2018"])/
    (som_EPA18))

vectLIEC_X_EPA=c(
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="LIEC_X_EPA" & Tab2$Annee=="2013"])/
    (som_EPA13),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="LIEC_X_EPA" & Tab2$Annee=="2014"])/
    (som_EPA14),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="LIEC_X_EPA" & Tab2$Annee=="2015"])/
    (som_EPA15),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="LIEC_X_EPA" & Tab2$Annee=="2016"])/
    (som_EPA16),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="LIEC_X_EPA" & Tab2$Annee=="2017"])/
    (som_EPA17),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="LIEC_X_EPA" & Tab2$Annee=="2018"])/
    (som_EPA18))




##### taux pour les _ANS #####


# calcule des sommes

som_ANStt=sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPI_ANS"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPA_ANS"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="ANS_ANS"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="END_ANS"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="LIEC_X_ANS"])

som_ANS13=sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPI_ANS" & Tab2$Annee=="2013"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPA_ANS" & Tab2$Annee=="2013"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="ANS_ANS" & Tab2$Annee=="2013"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="END_ANS" & Tab2$Annee=="2013"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="LIEC_X_ANS" & Tab2$Annee=="2013"])

som_ANS14=sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPI_ANS" & Tab2$Annee=="2014"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPA_ANS" & Tab2$Annee=="2014"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="ANS_ANS" & Tab2$Annee=="2014"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="END_ANS" & Tab2$Annee=="2014"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="LIEC_X_ANS" & Tab2$Annee=="2014"])

som_ANS15=sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPI_ANS" & Tab2$Annee=="2015"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPA_ANS" & Tab2$Annee=="2015"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="ANS_ANS" & Tab2$Annee=="2015"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="END_ANS" & Tab2$Annee=="2015"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="LIEC_X_ANS" & Tab2$Annee=="2015"])

som_ANS16=sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPI_ANS" & Tab2$Annee=="2016"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPA_ANS" & Tab2$Annee=="2016"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="ANS_ANS" & Tab2$Annee=="2016"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="END_ANS" & Tab2$Annee=="2016"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="LIEC_X_ANS" & Tab2$Annee=="2016"])

som_ANS17=sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPI_ANS" & Tab2$Annee=="2017"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPA_ANS" & Tab2$Annee=="2017"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="ANS_ANS" & Tab2$Annee=="2017"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="END_ANS" & Tab2$Annee=="2017"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="LIEC_X_ANS" & Tab2$Annee=="2017"])

som_ANS18=sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPI_ANS" & Tab2$Annee=="2018"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPA_ANS" & Tab2$Annee=="2018"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="ANS_ANS" & Tab2$Annee=="2018"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="END_ANS" & Tab2$Annee=="2018"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="LIEC_X_ANS" & Tab2$Annee=="2018"])






# vecteurs pour chaque annee

vectEPI_ANS=c(
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPI_ANS" & Tab2$Annee=="2013"])/
    (som_ANS13),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPI_ANS" & Tab2$Annee=="2014"])/
    (som_ANS14),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPI_ANS" & Tab2$Annee=="2015"])/
    (som_ANS15),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPI_ANS" & Tab2$Annee=="2016"])/
    (som_ANS16),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPI_ANS" & Tab2$Annee=="2017"])/
    (som_ANS17),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPI_ANS" & Tab2$Annee=="2018"])/
    (som_ANS18))

vectEPA_ANS=c(
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPA_ANS" & Tab2$Annee=="2013"])/
    (som_ANS13),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPA_ANS" & Tab2$Annee=="2014"])/
    (som_ANS14),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPA_ANS" & Tab2$Annee=="2015"])/
    (som_ANS15),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPA_ANS" & Tab2$Annee=="2016"])/
    (som_ANS16),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPA_ANS" & Tab2$Annee=="2017"])/
    (som_ANS17),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPA_ANS" & Tab2$Annee=="2018"])/
    (som_ANS18))

vectANS_ANS=c(
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="ANS_ANS" & Tab2$Annee=="2013"])/
    (som_ANS13),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="ANS_ANS" & Tab2$Annee=="2014"])/
    (som_ANS14),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="ANS_ANS" & Tab2$Annee=="2015"])/
    (som_ANS15),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="ANS_ANS" & Tab2$Annee=="2016"])/
    (som_ANS16),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="ANS_ANS" & Tab2$Annee=="2017"])/
    (som_ANS17),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="ANS_ANS" & Tab2$Annee=="2018"])/
    (som_ANS18))

vectEND_ANS=c(
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="END_ANS" & Tab2$Annee=="2013"])/
    (som_ANS13),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="END_ANS" & Tab2$Annee=="2014"])/
    (som_ANS14),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="END_ANS" & Tab2$Annee=="2015"])/
    (som_ANS15),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="END_ANS" & Tab2$Annee=="2016"])/
    (som_ANS16),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="END_ANS" & Tab2$Annee=="2017"])/
    (som_ANS17),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="END_ANS" & Tab2$Annee=="2018"])/
    (som_ANS18))

vectLIEC_X_ANS=c(
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="LIEC_X_ANS" & Tab2$Annee=="2013"])/
    (som_ANS13),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="LIEC_X_ANS" & Tab2$Annee=="2014"])/
    (som_ANS14),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="LIEC_X_ANS" & Tab2$Annee=="2015"])/
    (som_ANS15),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="LIEC_X_ANS" & Tab2$Annee=="2016"])/
    (som_ANS16),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="LIEC_X_ANS" & Tab2$Annee=="2017"])/
    (som_ANS17),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="LIEC_X_ANS" & Tab2$Annee=="2018"])/
    (som_ANS18))




##### taux pour les _END #####


# calcule des sommes

som_ENDtt=sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPI_END"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPA_END"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="ANS_END"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="END_END"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="LIEC_X_END"])

som_END13=sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPI_END" & Tab2$Annee=="2013"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPA_END" & Tab2$Annee=="2013"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="ANS_END" & Tab2$Annee=="2013"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="END_END" & Tab2$Annee=="2013"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="LIEC_X_END" & Tab2$Annee=="2013"])

som_END14=sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPI_END" & Tab2$Annee=="2014"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPA_END" & Tab2$Annee=="2014"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="ANS_END" & Tab2$Annee=="2014"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="END_END" & Tab2$Annee=="2014"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="LIEC_X_END" & Tab2$Annee=="2014"])

som_END15=sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPI_END" & Tab2$Annee=="2015"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPA_END" & Tab2$Annee=="2015"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="ANS_END" & Tab2$Annee=="2015"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="END_END" & Tab2$Annee=="2015"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="LIEC_X_END" & Tab2$Annee=="2015"])

som_END16=sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPI_END" & Tab2$Annee=="2016"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPA_END" & Tab2$Annee=="2016"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="ANS_END" & Tab2$Annee=="2016"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="END_END" & Tab2$Annee=="2016"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="LIEC_X_END" & Tab2$Annee=="2016"])

som_END17=sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPI_END" & Tab2$Annee=="2017"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPA_END" & Tab2$Annee=="2017"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="ANS_END" & Tab2$Annee=="2017"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="END_END" & Tab2$Annee=="2017"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="LIEC_X_END" & Tab2$Annee=="2017"])

som_END18=sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPI_END" & Tab2$Annee=="2018"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPA_END" & Tab2$Annee=="2018"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="ANS_END" & Tab2$Annee=="2018"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="END_END" & Tab2$Annee=="2018"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="LIEC_X_END" & Tab2$Annee=="2018"])






# vecteurs pour chaque annee

vectEPI_END=c(
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPI_END" & Tab2$Annee=="2013"])/
    (som_END13),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPI_END" & Tab2$Annee=="2014"])/
    (som_END14),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPI_END" & Tab2$Annee=="2015"])/
    (som_END15),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPI_END" & Tab2$Annee=="2016"])/
    (som_END16),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPI_END" & Tab2$Annee=="2017"])/
    (som_END17),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPI_END" & Tab2$Annee=="2018"])/
    (som_END18))

vectEPA_END=c(
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPA_END" & Tab2$Annee=="2013"])/
    (som_END13),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPA_END" & Tab2$Annee=="2014"])/
    (som_END14),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPA_END" & Tab2$Annee=="2015"])/
    (som_END15),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPA_END" & Tab2$Annee=="2016"])/
    (som_END16),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPA_END" & Tab2$Annee=="2017"])/
    (som_END17),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPA_END" & Tab2$Annee=="2018"])/
    (som_END18))

vectANS_END=c(
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="ANS_END" & Tab2$Annee=="2013"])/
    (som_END13),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="ANS_END" & Tab2$Annee=="2014"])/
    (som_END14),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="ANS_END" & Tab2$Annee=="2015"])/
    (som_END15),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="ANS_END" & Tab2$Annee=="2016"])/
    (som_END16),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="ANS_END" & Tab2$Annee=="2017"])/
    (som_END17),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="ANS_END" & Tab2$Annee=="2018"])/
    (som_END18))

vectEND_END=c(
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="END_END" & Tab2$Annee=="2013"])/
    (som_END13),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="END_END" & Tab2$Annee=="2014"])/
    (som_END14),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="END_END" & Tab2$Annee=="2015"])/
    (som_END15),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="END_END" & Tab2$Annee=="2016"])/
    (som_END16),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="END_END" & Tab2$Annee=="2017"])/
    (som_END17),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="END_END" & Tab2$Annee=="2018"])/
    (som_END18))

vectLIEC_X_END=c(
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="LIEC_X_END" & Tab2$Annee=="2013"])/
    (som_END13),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="LIEC_X_END" & Tab2$Annee=="2014"])/
    (som_END14),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="LIEC_X_END" & Tab2$Annee=="2015"])/
    (som_END15),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="LIEC_X_END" & Tab2$Annee=="2016"])/
    (som_END16),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="LIEC_X_END" & Tab2$Annee=="2017"])/
    (som_END17),
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="LIEC_X_END" & Tab2$Annee=="2018"])/
    (som_END18))





##### Graphiques #####

plot(2013:2018,vectEPI_EPI*100,col="#000099",type="l",lwd=4,ylim=c(0,131),
     main="Taux d'identification terrain r?ussi Aquitaine",
     xlab="",
     ylab="",
     yaxt='n')
lines(2013:2018,vectEPA_EPA*100,col="#FFFF00",type="l",lwd=4)
lines(2013:2018,vectANS_ANS*100,col="#009900",type="l",lwd=4)
lines(2013:2018,vectEND_END*100,col="#FF6699",type="l",lwd=4)
lines(2013:2018,vectReussiteNAT*100,col="black",type="l",lwd=5)
axis(side = 2, at=c(0,50,100),las=1,labels=c("0%","50%","100%"))
axis(side = 1, at=2013:2018, pos=-15, tick=F, labels=c(paste("n=",n2013,sep=""),
                                                       paste("n=",n2014,sep=""),paste("n=",n2015,sep=""),
                                                       paste("n=",n2016,sep=""),paste("n=",n2017,sep=""),
                                                       paste("n=",n2018,sep="")))

##### Region bien formee : Aquitaine #####

tabAQ=Tab2[Tab2$Region%in%"AQ",]
nrow(tabAQ)==1959
dataID <- unique(CS_error[,c("cle",colID,"Region")])
nrow(dataID)==1109
dataID=dataID[dataID$Region%in%"AQ",]
nrow(dataID)==129

n2013=table(dataID$Annee)[1]
n2014=table(dataID$Annee)[2]
n2015=table(dataID$Annee)[3]
n2016=table(dataID$Annee)[4]
n2017=table(dataID$Annee)[5]
n2018=table(dataID$Annee)[6]

vectReussiteAQ=c(
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC%in%c("EPI_EPI","EPA_EPA","ANS_ANS","END_END") & tabAQ$Annee%in%2013])/sum(tabAQ$Nbr_EW[tabAQ$Annee%in%2013]),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC%in%c("EPI_EPI","EPA_EPA","ANS_ANS","END_END") & tabAQ$Annee%in%2014])/sum(tabAQ$Nbr_EW[tabAQ$Annee%in%2014]),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC%in%c("EPI_EPI","EPA_EPA","ANS_ANS","END_END") & tabAQ$Annee%in%2015])/sum(tabAQ$Nbr_EW[tabAQ$Annee%in%2015]),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC%in%c("EPI_EPI","EPA_EPA","ANS_ANS","END_END") & tabAQ$Annee%in%2016])/sum(tabAQ$Nbr_EW[tabAQ$Annee%in%2016]),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC%in%c("EPI_EPI","EPA_EPA","ANS_ANS","END_END") & tabAQ$Annee%in%2017])/sum(tabAQ$Nbr_EW[tabAQ$Annee%in%2017]),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC%in%c("EPI_EPI","EPA_EPA","ANS_ANS","END_END") & tabAQ$Annee%in%2018])/sum(tabAQ$Nbr_EW[tabAQ$Annee%in%2018]))


##### taux pour les _EPI #####


# calcule des sommes

som_EPItt=sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPI_EPI"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPA_EPI"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="ANS_EPI"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="END_EPI"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="LIEC_X_EPI"])

som_EPI12=sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPI_EPI" & tabAQ$Annee=="2012"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPA_EPI" & tabAQ$Annee=="2012"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="ANS_EPI" & tabAQ$Annee=="2012"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="END_EPI" & tabAQ$Annee=="2012"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="LIEC_X_EPI" & tabAQ$Annee=="2012"])

som_EPI13=sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPI_EPI" & tabAQ$Annee=="2013"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPA_EPI" & tabAQ$Annee=="2013"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="ANS_EPI" & tabAQ$Annee=="2013"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="END_EPI" & tabAQ$Annee=="2013"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="LIEC_X_EPI" & tabAQ$Annee=="2013"])

som_EPI14=sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPI_EPI" & tabAQ$Annee=="2014"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPA_EPI" & tabAQ$Annee=="2014"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="ANS_EPI" & tabAQ$Annee=="2014"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="END_EPI" & tabAQ$Annee=="2014"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="LIEC_X_EPI" & tabAQ$Annee=="2014"])

som_EPI15=sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPI_EPI" & tabAQ$Annee=="2015"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPA_EPI" & tabAQ$Annee=="2015"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="ANS_EPI" & tabAQ$Annee=="2015"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="END_EPI" & tabAQ$Annee=="2015"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="LIEC_X_EPI" & tabAQ$Annee=="2015"])

som_EPI16=sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPI_EPI" & tabAQ$Annee=="2016"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPA_EPI" & tabAQ$Annee=="2016"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="ANS_EPI" & tabAQ$Annee=="2016"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="END_EPI" & tabAQ$Annee=="2016"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="LIEC_X_EPI" & tabAQ$Annee=="2016"])

som_EPI17=sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPI_EPI" & tabAQ$Annee=="2017"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPA_EPI" & tabAQ$Annee=="2017"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="ANS_EPI" & tabAQ$Annee=="2017"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="END_EPI" & tabAQ$Annee=="2017"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="LIEC_X_EPI" & tabAQ$Annee=="2017"])

som_EPI18=sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPI_EPI" & tabAQ$Annee=="2018"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPA_EPI" & tabAQ$Annee=="2018"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="ANS_EPI" & tabAQ$Annee=="2018"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="END_EPI" & tabAQ$Annee=="2018"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="LIEC_X_EPI" & tabAQ$Annee=="2018"])




# vecteurs pour chaque annee

vectEPI_EPI=c(
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPI_EPI" & tabAQ$Annee=="2012"])/
    (som_EPI12),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPI_EPI" & tabAQ$Annee=="2013"])/
    (som_EPI13),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPI_EPI" & tabAQ$Annee=="2014"])/
    (som_EPI14),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPI_EPI" & tabAQ$Annee=="2015"])/
    (som_EPI15),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPI_EPI" & tabAQ$Annee=="2016"])/
    (som_EPI16),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPI_EPI" & tabAQ$Annee=="2017"])/
    (som_EPI17),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPI_EPI" & tabAQ$Annee=="2018"])/
    (som_EPI18))

vectEPA_EPI=c(
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPA_EPI" & tabAQ$Annee=="2012"])/
    (som_EPI12),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPA_EPI" & tabAQ$Annee=="2013"])/
    (som_EPI13),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPA_EPI" & tabAQ$Annee=="2014"])/
    (som_EPI14),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPA_EPI" & tabAQ$Annee=="2015"])/
    (som_EPI15),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPA_EPI" & tabAQ$Annee=="2016"])/
    (som_EPI16),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPA_EPI" & tabAQ$Annee=="2017"])/
    (som_EPI17),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPA_EPI" & tabAQ$Annee=="2018"])/
    (som_EPI18))

vectANS_EPI=c(
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="ANS_EPI" & tabAQ$Annee=="2012"])/
    (som_EPI12),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="ANS_EPI" & tabAQ$Annee=="2013"])/
    (som_EPI13),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="ANS_EPI" & tabAQ$Annee=="2014"])/
    (som_EPI14),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="ANS_EPI" & tabAQ$Annee=="2015"])/
    (som_EPI15),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="ANS_EPI" & tabAQ$Annee=="2016"])/
    (som_EPI16),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="ANS_EPI" & tabAQ$Annee=="2017"])/
    (som_EPI17),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="ANS_EPI" & tabAQ$Annee=="2018"])/
    (som_EPI18))

vectEND_EPI=c(
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="END_EPI" & tabAQ$Annee=="2012"])/
    (som_EPI12),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="END_EPI" & tabAQ$Annee=="2013"])/
    (som_EPI13),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="END_EPI" & tabAQ$Annee=="2014"])/
    (som_EPI14),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="END_EPI" & tabAQ$Annee=="2015"])/
    (som_EPI15),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="END_EPI" & tabAQ$Annee=="2016"])/
    (som_EPI16),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="END_EPI" & tabAQ$Annee=="2017"])/
    (som_EPI17),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="END_EPI" & tabAQ$Annee=="2018"])/
    (som_EPI18))

vectLIEC_X_EPI=c(
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="LIEC_X_EPI" & tabAQ$Annee=="2012"])/
    (som_EPI12),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="LIEC_X_EPI" & tabAQ$Annee=="2013"])/
    (som_EPI13),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="LIEC_X_EPI" & tabAQ$Annee=="2014"])/
    (som_EPI14),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="LIEC_X_EPI" & tabAQ$Annee=="2015"])/
    (som_EPI15),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="LIEC_X_EPI" & tabAQ$Annee=="2016"])/
    (som_EPI16),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="LIEC_X_EPI" & tabAQ$Annee=="2017"])/
    (som_EPI17),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="LIEC_X_EPI" & tabAQ$Annee=="2018"])/
    (som_EPI18))




##### taux pour les _EPA #####


# calcule des sommes

som_EPAtt=sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPI_EPA"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPA_EPA"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="ANS_EPA"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="END_EPA"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="LIEC_X_EPA"])

som_EPA12=sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPI_EPA" & tabAQ$Annee=="2012"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPA_EPA" & tabAQ$Annee=="2012"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="ANS_EPA" & tabAQ$Annee=="2012"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="END_EPA" & tabAQ$Annee=="2012"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="LIEC_X_EPA" & tabAQ$Annee=="2012"])

som_EPA13=sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPI_EPA" & tabAQ$Annee=="2013"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPA_EPA" & tabAQ$Annee=="2013"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="ANS_EPA" & tabAQ$Annee=="2013"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="END_EPA" & tabAQ$Annee=="2013"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="LIEC_X_EPA" & tabAQ$Annee=="2013"])

som_EPA14=sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPI_EPA" & tabAQ$Annee=="2014"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPA_EPA" & tabAQ$Annee=="2014"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="ANS_EPA" & tabAQ$Annee=="2014"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="END_EPA" & tabAQ$Annee=="2014"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="LIEC_X_EPA" & tabAQ$Annee=="2014"])

som_EPA15=sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPI_EPA" & tabAQ$Annee=="2015"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPA_EPA" & tabAQ$Annee=="2015"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="ANS_EPA" & tabAQ$Annee=="2015"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="END_EPA" & tabAQ$Annee=="2015"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="LIEC_X_EPA" & tabAQ$Annee=="2015"])

som_EPA16=sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPI_EPA" & tabAQ$Annee=="2016"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPA_EPA" & tabAQ$Annee=="2016"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="ANS_EPA" & tabAQ$Annee=="2016"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="END_EPA" & tabAQ$Annee=="2016"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="LIEC_X_EPA" & tabAQ$Annee=="2016"])

som_EPA17=sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPI_EPA" & tabAQ$Annee=="2017"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPA_EPA" & tabAQ$Annee=="2017"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="ANS_EPA" & tabAQ$Annee=="2017"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="END_EPA" & tabAQ$Annee=="2017"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="LIEC_X_EPA" & tabAQ$Annee=="2017"])

som_EPA18=sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPI_EPA" & tabAQ$Annee=="2018"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPA_EPA" & tabAQ$Annee=="2018"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="ANS_EPA" & tabAQ$Annee=="2018"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="END_EPA" & tabAQ$Annee=="2018"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="LIEC_X_EPA" & tabAQ$Annee=="2018"])






# vecteurs pour chaque annee

vectEPI_EPA=c(
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPI_EPA" & tabAQ$Annee=="2012"])/
    (som_EPA12),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPI_EPA" & tabAQ$Annee=="2013"])/
    (som_EPA13),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPI_EPA" & tabAQ$Annee=="2014"])/
    (som_EPA14),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPI_EPA" & tabAQ$Annee=="2015"])/
    (som_EPA15),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPI_EPA" & tabAQ$Annee=="2016"])/
    (som_EPA16),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPI_EPA" & tabAQ$Annee=="2017"])/
    (som_EPA17),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPI_EPA" & tabAQ$Annee=="2018"])/
    (som_EPA18))

vectEPA_EPA=c(
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPA_EPA" & tabAQ$Annee=="2012"])/
    (som_EPA12),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPA_EPA" & tabAQ$Annee=="2013"])/
    (som_EPA13),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPA_EPA" & tabAQ$Annee=="2014"])/
    (som_EPA14),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPA_EPA" & tabAQ$Annee=="2015"])/
    (som_EPA15),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPA_EPA" & tabAQ$Annee=="2016"])/
    (som_EPA16),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPA_EPA" & tabAQ$Annee=="2017"])/
    (som_EPA17),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPA_EPA" & tabAQ$Annee=="2018"])/
    (som_EPA18))

vectANS_EPA=c(
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="ANS_EPA" & tabAQ$Annee=="2012"])/
    (som_EPA12),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="ANS_EPA" & tabAQ$Annee=="2013"])/
    (som_EPA13),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="ANS_EPA" & tabAQ$Annee=="2014"])/
    (som_EPA14),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="ANS_EPA" & tabAQ$Annee=="2015"])/
    (som_EPA15),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="ANS_EPA" & tabAQ$Annee=="2016"])/
    (som_EPA16),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="ANS_EPA" & tabAQ$Annee=="2017"])/
    (som_EPA17),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="ANS_EPA" & tabAQ$Annee=="2018"])/
    (som_EPA18))

vectEND_EPA=c(
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="END_EPA" & tabAQ$Annee=="2012"])/
    (som_EPA12),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="END_EPA" & tabAQ$Annee=="2013"])/
    (som_EPA13),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="END_EPA" & tabAQ$Annee=="2014"])/
    (som_EPA14),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="END_EPA" & tabAQ$Annee=="2015"])/
    (som_EPA15),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="END_EPA" & tabAQ$Annee=="2016"])/
    (som_EPA16),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="END_EPA" & tabAQ$Annee=="2017"])/
    (som_EPA17),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="END_EPA" & tabAQ$Annee=="2018"])/
    (som_EPA18))

vectLIEC_X_EPA=c(
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="LIEC_X_EPA" & tabAQ$Annee=="2012"])/
    (som_EPA12),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="LIEC_X_EPA" & tabAQ$Annee=="2013"])/
    (som_EPA13),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="LIEC_X_EPA" & tabAQ$Annee=="2014"])/
    (som_EPA14),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="LIEC_X_EPA" & tabAQ$Annee=="2015"])/
    (som_EPA15),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="LIEC_X_EPA" & tabAQ$Annee=="2016"])/
    (som_EPA16),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="LIEC_X_EPA" & tabAQ$Annee=="2017"])/
    (som_EPA17),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="LIEC_X_EPA" & tabAQ$Annee=="2018"])/
    (som_EPA18))




##### taux pour les _ANS #####


# calcule des sommes

som_ANStt=sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPI_ANS"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPA_ANS"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="ANS_ANS"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="END_ANS"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="LIEC_X_ANS"])

som_ANS12=sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPI_ANS" & tabAQ$Annee=="2012"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPA_ANS" & tabAQ$Annee=="2012"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="ANS_ANS" & tabAQ$Annee=="2012"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="END_ANS" & tabAQ$Annee=="2012"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="LIEC_X_ANS" & tabAQ$Annee=="2012"])

som_ANS13=sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPI_ANS" & tabAQ$Annee=="2013"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPA_ANS" & tabAQ$Annee=="2013"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="ANS_ANS" & tabAQ$Annee=="2013"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="END_ANS" & tabAQ$Annee=="2013"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="LIEC_X_ANS" & tabAQ$Annee=="2013"])

som_ANS14=sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPI_ANS" & tabAQ$Annee=="2014"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPA_ANS" & tabAQ$Annee=="2014"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="ANS_ANS" & tabAQ$Annee=="2014"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="END_ANS" & tabAQ$Annee=="2014"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="LIEC_X_ANS" & tabAQ$Annee=="2014"])

som_ANS15=sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPI_ANS" & tabAQ$Annee=="2015"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPA_ANS" & tabAQ$Annee=="2015"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="ANS_ANS" & tabAQ$Annee=="2015"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="END_ANS" & tabAQ$Annee=="2015"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="LIEC_X_ANS" & tabAQ$Annee=="2015"])

som_ANS16=sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPI_ANS" & tabAQ$Annee=="2016"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPA_ANS" & tabAQ$Annee=="2016"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="ANS_ANS" & tabAQ$Annee=="2016"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="END_ANS" & tabAQ$Annee=="2016"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="LIEC_X_ANS" & tabAQ$Annee=="2016"])

som_ANS17=sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPI_ANS" & tabAQ$Annee=="2017"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPA_ANS" & tabAQ$Annee=="2017"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="ANS_ANS" & tabAQ$Annee=="2017"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="END_ANS" & tabAQ$Annee=="2017"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="LIEC_X_ANS" & tabAQ$Annee=="2017"])

som_ANS18=sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPI_ANS" & tabAQ$Annee=="2018"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPA_ANS" & tabAQ$Annee=="2018"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="ANS_ANS" & tabAQ$Annee=="2018"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="END_ANS" & tabAQ$Annee=="2018"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="LIEC_X_ANS" & tabAQ$Annee=="2018"])






# vecteurs pour chaque annee

vectEPI_ANS=c(
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPI_ANS" & tabAQ$Annee=="2012"])/
    (som_ANS12),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPI_ANS" & tabAQ$Annee=="2013"])/
    (som_ANS13),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPI_ANS" & tabAQ$Annee=="2014"])/
    (som_ANS14),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPI_ANS" & tabAQ$Annee=="2015"])/
    (som_ANS15),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPI_ANS" & tabAQ$Annee=="2016"])/
    (som_ANS16),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPI_ANS" & tabAQ$Annee=="2017"])/
    (som_ANS17),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPI_ANS" & tabAQ$Annee=="2018"])/
    (som_ANS18))

vectEPA_ANS=c(
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPA_ANS" & tabAQ$Annee=="2012"])/
    (som_ANS12),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPA_ANS" & tabAQ$Annee=="2013"])/
    (som_ANS13),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPA_ANS" & tabAQ$Annee=="2014"])/
    (som_ANS14),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPA_ANS" & tabAQ$Annee=="2015"])/
    (som_ANS15),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPA_ANS" & tabAQ$Annee=="2016"])/
    (som_ANS16),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPA_ANS" & tabAQ$Annee=="2017"])/
    (som_ANS17),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPA_ANS" & tabAQ$Annee=="2018"])/
    (som_ANS18))

vectANS_ANS=c(
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="ANS_ANS" & tabAQ$Annee=="2012"])/
    (som_ANS12),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="ANS_ANS" & tabAQ$Annee=="2013"])/
    (som_ANS13),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="ANS_ANS" & tabAQ$Annee=="2014"])/
    (som_ANS14),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="ANS_ANS" & tabAQ$Annee=="2015"])/
    (som_ANS15),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="ANS_ANS" & tabAQ$Annee=="2016"])/
    (som_ANS16),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="ANS_ANS" & tabAQ$Annee=="2017"])/
    (som_ANS17),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="ANS_ANS" & tabAQ$Annee=="2018"])/
    (som_ANS18))

vectEND_ANS=c(
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="END_ANS" & tabAQ$Annee=="2012"])/
    (som_ANS12),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="END_ANS" & tabAQ$Annee=="2013"])/
    (som_ANS13),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="END_ANS" & tabAQ$Annee=="2014"])/
    (som_ANS14),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="END_ANS" & tabAQ$Annee=="2015"])/
    (som_ANS15),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="END_ANS" & tabAQ$Annee=="2016"])/
    (som_ANS16),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="END_ANS" & tabAQ$Annee=="2017"])/
    (som_ANS17),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="END_ANS" & tabAQ$Annee=="2018"])/
    (som_ANS18))

vectLIEC_X_ANS=c(
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="LIEC_X_ANS" & tabAQ$Annee=="2012"])/
    (som_ANS12),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="LIEC_X_ANS" & tabAQ$Annee=="2013"])/
    (som_ANS13),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="LIEC_X_ANS" & tabAQ$Annee=="2014"])/
    (som_ANS14),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="LIEC_X_ANS" & tabAQ$Annee=="2015"])/
    (som_ANS15),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="LIEC_X_ANS" & tabAQ$Annee=="2016"])/
    (som_ANS16),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="LIEC_X_ANS" & tabAQ$Annee=="2017"])/
    (som_ANS17),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="LIEC_X_ANS" & tabAQ$Annee=="2018"])/
    (som_ANS18))




##### taux pour les _END #####


# calcule des sommes

som_ENDtt=sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPI_END"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPA_END"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="ANS_END"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="END_END"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="LIEC_X_END"])

som_END12=sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPI_END" & tabAQ$Annee=="2012"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPA_END" & tabAQ$Annee=="2012"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="ANS_END" & tabAQ$Annee=="2012"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="END_END" & tabAQ$Annee=="2012"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="LIEC_X_END" & tabAQ$Annee=="2012"])

som_END13=sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPI_END" & tabAQ$Annee=="2013"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPA_END" & tabAQ$Annee=="2013"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="ANS_END" & tabAQ$Annee=="2013"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="END_END" & tabAQ$Annee=="2013"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="LIEC_X_END" & tabAQ$Annee=="2013"])

som_END14=sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPI_END" & tabAQ$Annee=="2014"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPA_END" & tabAQ$Annee=="2014"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="ANS_END" & tabAQ$Annee=="2014"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="END_END" & tabAQ$Annee=="2014"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="LIEC_X_END" & tabAQ$Annee=="2014"])

som_END15=sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPI_END" & tabAQ$Annee=="2015"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPA_END" & tabAQ$Annee=="2015"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="ANS_END" & tabAQ$Annee=="2015"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="END_END" & tabAQ$Annee=="2015"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="LIEC_X_END" & tabAQ$Annee=="2015"])

som_END16=sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPI_END" & tabAQ$Annee=="2016"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPA_END" & tabAQ$Annee=="2016"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="ANS_END" & tabAQ$Annee=="2016"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="END_END" & tabAQ$Annee=="2016"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="LIEC_X_END" & tabAQ$Annee=="2016"])

som_END17=sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPI_END" & tabAQ$Annee=="2017"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPA_END" & tabAQ$Annee=="2017"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="ANS_END" & tabAQ$Annee=="2017"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="END_END" & tabAQ$Annee=="2017"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="LIEC_X_END" & tabAQ$Annee=="2017"])

som_END18=sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPI_END" & tabAQ$Annee=="2018"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPA_END" & tabAQ$Annee=="2018"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="ANS_END" & tabAQ$Annee=="2018"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="END_END" & tabAQ$Annee=="2018"])+
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="LIEC_X_END" & tabAQ$Annee=="2018"])






# vecteurs pour chaque annee

vectEPI_END=c(
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPI_END" & tabAQ$Annee=="2012"])/
    (som_END12),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPI_END" & tabAQ$Annee=="2013"])/
    (som_END13),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPI_END" & tabAQ$Annee=="2014"])/
    (som_END14),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPI_END" & tabAQ$Annee=="2015"])/
    (som_END15),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPI_END" & tabAQ$Annee=="2016"])/
    (som_END16),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPI_END" & tabAQ$Annee=="2017"])/
    (som_END17),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPI_END" & tabAQ$Annee=="2018"])/
    (som_END18))

vectEPA_END=c(
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPA_END" & tabAQ$Annee=="2012"])/
    (som_END12),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPA_END" & tabAQ$Annee=="2013"])/
    (som_END13),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPA_END" & tabAQ$Annee=="2014"])/
    (som_END14),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPA_END" & tabAQ$Annee=="2015"])/
    (som_END15),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPA_END" & tabAQ$Annee=="2016"])/
    (som_END16),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPA_END" & tabAQ$Annee=="2017"])/
    (som_END17),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="EPA_END" & tabAQ$Annee=="2018"])/
    (som_END18))

vectANS_END=c(
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="ANS_END" & tabAQ$Annee=="2012"])/
    (som_END12),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="ANS_END" & tabAQ$Annee=="2013"])/
    (som_END13),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="ANS_END" & tabAQ$Annee=="2014"])/
    (som_END14),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="ANS_END" & tabAQ$Annee=="2015"])/
    (som_END15),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="ANS_END" & tabAQ$Annee=="2016"])/
    (som_END16),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="ANS_END" & tabAQ$Annee=="2017"])/
    (som_END17),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="ANS_END" & tabAQ$Annee=="2018"])/
    (som_END18))

vectEND_END=c(
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="END_END" & tabAQ$Annee=="2012"])/
    (som_END12),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="END_END" & tabAQ$Annee=="2013"])/
    (som_END13),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="END_END" & tabAQ$Annee=="2014"])/
    (som_END14),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="END_END" & tabAQ$Annee=="2015"])/
    (som_END15),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="END_END" & tabAQ$Annee=="2016"])/
    (som_END16),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="END_END" & tabAQ$Annee=="2017"])/
    (som_END17),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="END_END" & tabAQ$Annee=="2018"])/
    (som_END18))

vectLIEC_X_END=c(
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="LIEC_X_END" & tabAQ$Annee=="2012"])/
    (som_END12),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="LIEC_X_END" & tabAQ$Annee=="2013"])/
    (som_END13),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="LIEC_X_END" & tabAQ$Annee=="2014"])/
    (som_END14),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="LIEC_X_END" & tabAQ$Annee=="2015"])/
    (som_END15),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="LIEC_X_END" & tabAQ$Annee=="2016"])/
    (som_END16),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="LIEC_X_END" & tabAQ$Annee=="2017"])/
    (som_END17),
  sum(tabAQ$Nbr_EW[tabAQ$FIEC_LIEC=="LIEC_X_END" & tabAQ$Annee=="2018"])/
    (som_END18))



##### Graphiques #####




plot(2013:2018,vectEPI_EPI[2:7]*100,col="#000099",type="l",lwd=4,ylim=c(0,131),
     main="Taux d'identification terrain r?ussi Aquitaine",
     xlab="",
     ylab="",
     yaxt='n')
lines(2013:2018,vectEPA_EPA[2:7]*100,col="#FFFF00",type="l",lwd=4)
lines(2013:2018,vectANS_ANS[2:7]*100,col="#009900",type="l",lwd=4)
lines(2013:2018,vectEND_END[2:7]*100,col="#FF6699",type="l",lwd=4)
lines(2013:2018,vectReussiteAQ*100,col="black",type="l",lwd=5)
axis(side = 2, at=c(0,50,100),las=1,labels=c("0%","50%","100%"))
axis(side = 1, at=2013:2018, pos=-15, tick=F, labels=c(paste("n=",n2013,sep=""),
                                                       paste("n=",n2014,sep=""),paste("n=",n2015,sep=""),
                                                       paste("n=",n2016,sep=""),paste("n=",n2017,sep=""),
                                                       paste("n=",n2018,sep="")))
legend("topleft",
       legend=c("EPI","EPA","ANS","END"),
       fill=c("#000099","#FFFF00","#009900","#FF6699"))


### Appartenances Reelles ###
{
    
  # erreures EPI
  plot(2013:2018,vectEPI_EPI[2:7]*100,col="#000099",type="l",lwd=5,ylim=c(0,131),
       main="Appartenance r?elle des Epig?s identifi?s sur le terrain Aquitaine",
       xlab="",
       ylab="",
       yaxt='n')
  lines(2013:2018,vectEPI_EPA[2:7]*100,col="#FFFF00",type="l",lwd=5)
  lines(2013:2018,vectEPI_ANS[2:7]*100,col="#009900",type="l",lwd=5)
  lines(2013:2018,vectEPI_END[2:7]*100,col="#FF6699",type="l",lwd=5)
  #lines(2013:2018,vectEPI_LIEC_X*100,col="grey",type="l",lwd=5)
  axis(side = 2, at=c(0,50,100),las=1,labels=c("0%","50%","100%"))
  axis(side = 1, at=2013:2018, pos=-15, tick=F, labels=c(paste("n=",n2013,sep=""),
                                                         paste("n=",n2014,sep=""),paste("n=",n2015,sep=""),
                                                         paste("n=",n2016,sep=""),paste("n=",n2017,sep=""),
                                                         paste("n=",n2018,sep="")))
  legend("topleft",
         legend=c("EPI","EPA","ANS","END"),
         fill=c("#000099","#FFFF00","#009900","#FF6699"))
  
  
  
  # erreures EPA
  plot(2013:2018,vectEPA_EPI[2:7]*100,col="#000099",type="l",lwd=5,ylim=c(0,131),
       main="Appartenance r?elle des Epi-An?ciques identifi?s sur le terrain Aquitaine",
       xlab="",
       ylab="",
       las=1,
       yaxt='n')
  lines(2013:2018,vectEPA_EPA[2:7]*100,col="#FFFF00",type="l",lwd=5)
  lines(2013:2018,vectEPA_ANS[2:7]*100,col="#009900",type="l",lwd=5)
  lines(2013:2018,vectEPA_END[2:7]*100,col="#FF6699",type="l",lwd=5)
  #lines(2013:2018,vectEPA_LIEC_X*100,col="grey",type="l",lwd=5)
  axis(side = 2, at=c(0,50,100),las=1,labels=c("0%","50%","100%"))
  axis(side = 1, at=2013:2018, pos=-15, tick=F, labels=c(paste("n=",n2013,sep=""),
                                                         paste("n=",n2014,sep=""),paste("n=",n2015,sep=""),
                                                         paste("n=",n2016,sep=""),paste("n=",n2017,sep=""),
                                                         paste("n=",n2018,sep="")))
  legend("topleft",
         legend=c("EPI","EPA","ANS","END"),
         fill=c("#000099","#FFFF00","#009900","#FF6699"))
  
  
  # erreures ANS
  plot(2013:2018,vectANS_EPI[2:7]*100,col="#000099",type="l",lwd=5,ylim=c(0,131),
       main="Appartenance r?elle des An?ciques-Stricts identifi?s sur le terrain Aquitaine",
       xlab="",
       ylab="",
       las=1,
       yaxt='n')
  lines(2013:2018,vectANS_EPA[2:7]*100,col="#FFFF00",type="l",lwd=5)
  lines(2013:2018,vectANS_ANS[2:7]*100,col="#009900",type="l",lwd=5)
  lines(2013:2018,vectANS_END[2:7]*100,col="#FF6699",type="l",lwd=5)
  #lines(2013:2018,vectANS_LIEC_X*100,col="grey",type="l",lwd=5)
  axis(side = 2, at=c(0,50,100),las=1,labels=c("0%","50%","100%"))
  axis(side = 1, at=2013:2018, pos=-15, tick=F, labels=c(paste("n=",n2013,sep=""),
                                                         paste("n=",n2014,sep=""),paste("n=",n2015,sep=""),
                                                         paste("n=",n2016,sep=""),paste("n=",n2017,sep=""),
                                                         paste("n=",n2018,sep="")))
  legend("topleft",
         legend=c("EPI","EPA","ANS","END"),
         fill=c("#000099","#FFFF00","#009900","#FF6699"))
  
  
  # erreures END
  plot(2013:2018,vectEND_EPI[2:7]*100,col="#000099",type="l",lwd=5,ylim=c(0,131),
       main="Appartenance r?elle des Endog?s identifi?s sur le terrain Aquitaine",
       xlab="",
       ylab="",
       las=1,
       yaxt='n')
  lines(2013:2018,vectEND_EPA[2:7]*100,col="#FFFF00",type="l",lwd=5)
  lines(2013:2018,vectEND_ANS[2:7]*100,col="#009900",type="l",lwd=5)
  lines(2013:2018,vectEND_END[2:7]*100,col="#FF6699",type="l",lwd=5)
  #lines(2013:2018,vectEND_LIEC_X*100,col="grey",type="l",lwd=5)
  axis(side = 2, at=c(0,50,100),las=1,labels=c("0%","50%","100%"))
  axis(side = 1, at=2013:2018, pos=-15, tick=F, labels=c(paste("n=",n2013,sep=""),
                                                         paste("n=",n2014,sep=""),paste("n=",n2015,sep=""),
                                                         paste("n=",n2016,sep=""),paste("n=",n2017,sep=""),
                                                         paste("n=",n2018,sep="")))
  legend("topleft",
         legend=c("EPI","EPA","ANS","END"),
         fill=c("#000099","#FFFF00","#009900","#FF6699"))
  
}


### Erreurs par LIEC ###
{
  
  # erreures EPI
  plot(2013:2018,vectEPI_EPI[2:7]*100,col="#000099",type="l",lwd=5,ylim=c(0,141),
       main="Erreur d'identification terrain des Epig?s Aquitaine",
       xlab="",
       ylab="",
       yaxt='n')
  lines(2013:2018,vectEPA_EPI[2:7]*100,col="#FFFF00",type="l",lwd=5)
  lines(2013:2018,vectANS_EPI[2:7]*100,col="#009900",type="l",lwd=5)
  lines(2013:2018,vectEND_EPI[2:7]*100,col="#FF6699",type="l",lwd=5)
  lines(2013:2018,vectLIEC_X_EPI[2:7]*100,col="grey",type="l",lwd=5)
  axis(side = 2, at=c(0,50,100),las=1,labels=c("0%","50%","100%"))
  axis(side = 1, at=2013:2018, pos=-15, tick=F, labels=c(paste("n=",n2013,sep=""),
                                                         paste("n=",n2014,sep=""),paste("n=",n2015,sep=""),
                                                         paste("n=",n2016,sep=""),paste("n=",n2017,sep=""),
                                                         paste("n=",n2018,sep="")))
  legend("topleft",
         legend=c("EPI","EPA","ANS","END","LIEC_X"),
         fill=c("#000099","#FFFF00","#009900","#FF6699","grey"))
  
  
  
  # erreures EPA
  plot(2013:2018,vectEPI_EPA[2:7]*100,col="#000099",type="l",lwd=5,ylim=c(0,141),
       main="Erreur d'identification terrain des Epi-An?ciques Aquitaine",
       xlab="",
       ylab="",
       las=1,
       yaxt='n')
  lines(2013:2018,vectEPA_EPA[2:7]*100,col="#FFFF00",type="l",lwd=5)
  lines(2013:2018,vectANS_EPA[2:7]*100,col="#009900",type="l",lwd=5)
  lines(2013:2018,vectEND_EPA[2:7]*100,col="#FF6699",type="l",lwd=5)
  lines(2013:2018,vectLIEC_X_EPA[2:7]*100,col="grey",type="l",lwd=5)
  axis(side = 2, at=c(0,50,100),las=1,labels=c("0%","50%","100%"))
  axis(side = 1, at=2013:2018, pos=-15, tick=F, labels=c(paste("n=",n2013,sep=""),
                                                         paste("n=",n2014,sep=""),paste("n=",n2015,sep=""),
                                                         paste("n=",n2016,sep=""),paste("n=",n2017,sep=""),
                                                         paste("n=",n2018,sep="")))
  legend("topleft",
         legend=c("EPI","EPA","ANS","END","LIEC_X"),
         fill=c("#000099","#FFFF00","#009900","#FF6699","grey"))
  
  
  # erreures ANS
  plot(2013:2018,vectEPI_ANS[2:7]*100,col="#000099",type="l",lwd=5,ylim=c(0,141),
       main="Erreur d'identification terrain des An?ciques-Stricts Aquitaine",
       xlab="",
       ylab="",
       las=1,
       yaxt='n')
  lines(2013:2018,vectEPA_ANS[2:7]*100,col="#FFFF00",type="l",lwd=5)
  lines(2013:2018,vectANS_ANS[2:7]*100,col="#009900",type="l",lwd=5)
  lines(2013:2018,vectEND_ANS[2:7]*100,col="#FF6699",type="l",lwd=5)
  lines(2013:2018,vectLIEC_X_ANS[2:7]*100,col="grey",type="l",lwd=5)
  axis(side = 2, at=c(0,50,100),las=1,labels=c("0%","50%","100%"))
  axis(side = 1, at=2013:2018, pos=-15, tick=F, labels=c(paste("n=",n2013,sep=""),
                                                         paste("n=",n2014,sep=""),paste("n=",n2015,sep=""),
                                                         paste("n=",n2016,sep=""),paste("n=",n2017,sep=""),
                                                         paste("n=",n2018,sep="")))
  legend("topleft",
         legend=c("EPI","EPA","ANS","END","LIEC_X"),
         fill=c("#000099","#FFFF00","#009900","#FF6699","grey"))
  
  
  # erreures END
  plot(2013:2018,vectEPI_END[2:7]*100,col="#000099",type="l",lwd=5,ylim=c(0,141),
       main="Erreur d'identification terrain des Endog?s Aquitaine",
       xlab="",
       ylab="",
       las=1,
       yaxt='n')
  lines(2013:2018,vectEPA_END[2:7]*100,col="#FFFF00",type="l",lwd=5)
  lines(2013:2018,vectANS_END[2:7]*100,col="#009900",type="l",lwd=5)
  lines(2013:2018,vectEND_END[2:7]*100,col="#FF6699",type="l",lwd=5)
  lines(2013:2018,vectLIEC_X_END[2:7]*100,col="grey",type="l",lwd=5)
  axis(side = 2, at=c(0,50,100),las=1,labels=c("0%","50%","100%"))
  axis(side = 1, at=2013:2018, pos=-15, tick=F, labels=c(paste("n=",n2013,sep=""),
                                                         paste("n=",n2014,sep=""),paste("n=",n2015,sep=""),
                                                         paste("n=",n2016,sep=""),paste("n=",n2017,sep=""),
                                                         paste("n=",n2018,sep="")))
  legend("topleft",
         legend=c("EPI","EPA","ANS","END","LIEC_X"),
         fill=c("#000099","#FFFF00","#009900","#FF6699","grey"))
  
}








##### Region mal formee : Bourgogne #####

tabBO=Tab2[Tab2$Region%in%"BO",]
nrow(tabBO)==3426
dataID <- unique(sbt_m[,c("cle",colID,"Region")])
nrow(dataID)==1109
dataID=dataID[dataID$Region%in%"BO",]
nrow(dataID)==144

n2013=table(dataID$Annee)[2]
n2014=table(dataID$Annee)[3]
n2015=table(dataID$Annee)[4]
n2016=table(dataID$Annee)[5]
n2017=table(dataID$Annee)[6]
n2018=table(dataID$Annee)[7]

vectReussiteBO=c(
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC%in%c("EPI_EPI","EPA_EPA","ANS_ANS","END_END") & tabBO$Annee%in%2013])/sum(tabBO$Nbr_EW[tabBO$Annee%in%2013]),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC%in%c("EPI_EPI","EPA_EPA","ANS_ANS","END_END") & tabBO$Annee%in%2014])/sum(tabBO$Nbr_EW[tabBO$Annee%in%2014]),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC%in%c("EPI_EPI","EPA_EPA","ANS_ANS","END_END") & tabBO$Annee%in%2015])/sum(tabBO$Nbr_EW[tabBO$Annee%in%2015]),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC%in%c("EPI_EPI","EPA_EPA","ANS_ANS","END_END") & tabBO$Annee%in%2016])/sum(tabBO$Nbr_EW[tabBO$Annee%in%2016]),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC%in%c("EPI_EPI","EPA_EPA","ANS_ANS","END_END") & tabBO$Annee%in%2017])/sum(tabBO$Nbr_EW[tabBO$Annee%in%2017]),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC%in%c("EPI_EPI","EPA_EPA","ANS_ANS","END_END") & tabBO$Annee%in%2018])/sum(tabBO$Nbr_EW[tabBO$Annee%in%2018]))


##### taux pour les _EPI #####


# calcule des sommes

som_EPItt=sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPI_EPI"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPA_EPI"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="ANS_EPI"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="END_EPI"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="LIEC_X_EPI"])

som_EPI12=sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPI_EPI" & tabBO$Annee=="2012"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPA_EPI" & tabBO$Annee=="2012"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="ANS_EPI" & tabBO$Annee=="2012"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="END_EPI" & tabBO$Annee=="2012"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="LIEC_X_EPI" & tabBO$Annee=="2012"])

som_EPI13=sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPI_EPI" & tabBO$Annee=="2013"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPA_EPI" & tabBO$Annee=="2013"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="ANS_EPI" & tabBO$Annee=="2013"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="END_EPI" & tabBO$Annee=="2013"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="LIEC_X_EPI" & tabBO$Annee=="2013"])

som_EPI14=sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPI_EPI" & tabBO$Annee=="2014"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPA_EPI" & tabBO$Annee=="2014"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="ANS_EPI" & tabBO$Annee=="2014"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="END_EPI" & tabBO$Annee=="2014"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="LIEC_X_EPI" & tabBO$Annee=="2014"])

som_EPI15=sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPI_EPI" & tabBO$Annee=="2015"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPA_EPI" & tabBO$Annee=="2015"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="ANS_EPI" & tabBO$Annee=="2015"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="END_EPI" & tabBO$Annee=="2015"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="LIEC_X_EPI" & tabBO$Annee=="2015"])

som_EPI16=sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPI_EPI" & tabBO$Annee=="2016"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPA_EPI" & tabBO$Annee=="2016"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="ANS_EPI" & tabBO$Annee=="2016"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="END_EPI" & tabBO$Annee=="2016"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="LIEC_X_EPI" & tabBO$Annee=="2016"])

som_EPI17=sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPI_EPI" & tabBO$Annee=="2017"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPA_EPI" & tabBO$Annee=="2017"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="ANS_EPI" & tabBO$Annee=="2017"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="END_EPI" & tabBO$Annee=="2017"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="LIEC_X_EPI" & tabBO$Annee=="2017"])

som_EPI18=sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPI_EPI" & tabBO$Annee=="2018"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPA_EPI" & tabBO$Annee=="2018"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="ANS_EPI" & tabBO$Annee=="2018"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="END_EPI" & tabBO$Annee=="2018"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="LIEC_X_EPI" & tabBO$Annee=="2018"])




# vecteurs pour chaque annee

vectEPI_EPI=c(
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPI_EPI" & tabBO$Annee=="2012"])/
    (som_EPI12),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPI_EPI" & tabBO$Annee=="2013"])/
    (som_EPI13),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPI_EPI" & tabBO$Annee=="2014"])/
    (som_EPI14),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPI_EPI" & tabBO$Annee=="2015"])/
    (som_EPI15),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPI_EPI" & tabBO$Annee=="2016"])/
    (som_EPI16),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPI_EPI" & tabBO$Annee=="2017"])/
    (som_EPI17),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPI_EPI" & tabBO$Annee=="2018"])/
    (som_EPI18))

vectEPA_EPI=c(
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPA_EPI" & tabBO$Annee=="2012"])/
    (som_EPI12),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPA_EPI" & tabBO$Annee=="2013"])/
    (som_EPI13),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPA_EPI" & tabBO$Annee=="2014"])/
    (som_EPI14),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPA_EPI" & tabBO$Annee=="2015"])/
    (som_EPI15),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPA_EPI" & tabBO$Annee=="2016"])/
    (som_EPI16),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPA_EPI" & tabBO$Annee=="2017"])/
    (som_EPI17),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPA_EPI" & tabBO$Annee=="2018"])/
    (som_EPI18))

vectANS_EPI=c(
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="ANS_EPI" & tabBO$Annee=="2012"])/
    (som_EPI12),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="ANS_EPI" & tabBO$Annee=="2013"])/
    (som_EPI13),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="ANS_EPI" & tabBO$Annee=="2014"])/
    (som_EPI14),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="ANS_EPI" & tabBO$Annee=="2015"])/
    (som_EPI15),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="ANS_EPI" & tabBO$Annee=="2016"])/
    (som_EPI16),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="ANS_EPI" & tabBO$Annee=="2017"])/
    (som_EPI17),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="ANS_EPI" & tabBO$Annee=="2018"])/
    (som_EPI18))

vectEND_EPI=c(
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="END_EPI" & tabBO$Annee=="2012"])/
    (som_EPI12),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="END_EPI" & tabBO$Annee=="2013"])/
    (som_EPI13),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="END_EPI" & tabBO$Annee=="2014"])/
    (som_EPI14),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="END_EPI" & tabBO$Annee=="2015"])/
    (som_EPI15),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="END_EPI" & tabBO$Annee=="2016"])/
    (som_EPI16),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="END_EPI" & tabBO$Annee=="2017"])/
    (som_EPI17),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="END_EPI" & tabBO$Annee=="2018"])/
    (som_EPI18))

vectLIEC_X_EPI=c(
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="LIEC_X_EPI" & tabBO$Annee=="2012"])/
    (som_EPI12),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="LIEC_X_EPI" & tabBO$Annee=="2013"])/
    (som_EPI13),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="LIEC_X_EPI" & tabBO$Annee=="2014"])/
    (som_EPI14),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="LIEC_X_EPI" & tabBO$Annee=="2015"])/
    (som_EPI15),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="LIEC_X_EPI" & tabBO$Annee=="2016"])/
    (som_EPI16),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="LIEC_X_EPI" & tabBO$Annee=="2017"])/
    (som_EPI17),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="LIEC_X_EPI" & tabBO$Annee=="2018"])/
    (som_EPI18))




##### taux pour les _EPA #####


# calcule des sommes

som_EPAtt=sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPI_EPA"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPA_EPA"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="ANS_EPA"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="END_EPA"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="LIEC_X_EPA"])

som_EPA12=sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPI_EPA" & tabBO$Annee=="2012"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPA_EPA" & tabBO$Annee=="2012"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="ANS_EPA" & tabBO$Annee=="2012"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="END_EPA" & tabBO$Annee=="2012"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="LIEC_X_EPA" & tabBO$Annee=="2012"])

som_EPA13=sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPI_EPA" & tabBO$Annee=="2013"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPA_EPA" & tabBO$Annee=="2013"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="ANS_EPA" & tabBO$Annee=="2013"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="END_EPA" & tabBO$Annee=="2013"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="LIEC_X_EPA" & tabBO$Annee=="2013"])

som_EPA14=sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPI_EPA" & tabBO$Annee=="2014"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPA_EPA" & tabBO$Annee=="2014"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="ANS_EPA" & tabBO$Annee=="2014"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="END_EPA" & tabBO$Annee=="2014"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="LIEC_X_EPA" & tabBO$Annee=="2014"])

som_EPA15=sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPI_EPA" & tabBO$Annee=="2015"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPA_EPA" & tabBO$Annee=="2015"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="ANS_EPA" & tabBO$Annee=="2015"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="END_EPA" & tabBO$Annee=="2015"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="LIEC_X_EPA" & tabBO$Annee=="2015"])

som_EPA16=sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPI_EPA" & tabBO$Annee=="2016"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPA_EPA" & tabBO$Annee=="2016"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="ANS_EPA" & tabBO$Annee=="2016"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="END_EPA" & tabBO$Annee=="2016"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="LIEC_X_EPA" & tabBO$Annee=="2016"])

som_EPA17=sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPI_EPA" & tabBO$Annee=="2017"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPA_EPA" & tabBO$Annee=="2017"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="ANS_EPA" & tabBO$Annee=="2017"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="END_EPA" & tabBO$Annee=="2017"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="LIEC_X_EPA" & tabBO$Annee=="2017"])

som_EPA18=sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPI_EPA" & tabBO$Annee=="2018"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPA_EPA" & tabBO$Annee=="2018"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="ANS_EPA" & tabBO$Annee=="2018"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="END_EPA" & tabBO$Annee=="2018"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="LIEC_X_EPA" & tabBO$Annee=="2018"])






# vecteurs pour chaque annee

vectEPI_EPA=c(
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPI_EPA" & tabBO$Annee=="2012"])/
    (som_EPA12),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPI_EPA" & tabBO$Annee=="2013"])/
    (som_EPA13),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPI_EPA" & tabBO$Annee=="2014"])/
    (som_EPA14),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPI_EPA" & tabBO$Annee=="2015"])/
    (som_EPA15),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPI_EPA" & tabBO$Annee=="2016"])/
    (som_EPA16),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPI_EPA" & tabBO$Annee=="2017"])/
    (som_EPA17),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPI_EPA" & tabBO$Annee=="2018"])/
    (som_EPA18))

vectEPA_EPA=c(
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPA_EPA" & tabBO$Annee=="2012"])/
    (som_EPA12),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPA_EPA" & tabBO$Annee=="2013"])/
    (som_EPA13),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPA_EPA" & tabBO$Annee=="2014"])/
    (som_EPA14),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPA_EPA" & tabBO$Annee=="2015"])/
    (som_EPA15),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPA_EPA" & tabBO$Annee=="2016"])/
    (som_EPA16),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPA_EPA" & tabBO$Annee=="2017"])/
    (som_EPA17),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPA_EPA" & tabBO$Annee=="2018"])/
    (som_EPA18))

vectANS_EPA=c(
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="ANS_EPA" & tabBO$Annee=="2012"])/
    (som_EPA12),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="ANS_EPA" & tabBO$Annee=="2013"])/
    (som_EPA13),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="ANS_EPA" & tabBO$Annee=="2014"])/
    (som_EPA14),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="ANS_EPA" & tabBO$Annee=="2015"])/
    (som_EPA15),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="ANS_EPA" & tabBO$Annee=="2016"])/
    (som_EPA16),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="ANS_EPA" & tabBO$Annee=="2017"])/
    (som_EPA17),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="ANS_EPA" & tabBO$Annee=="2018"])/
    (som_EPA18))

vectEND_EPA=c(
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="END_EPA" & tabBO$Annee=="2012"])/
    (som_EPA12),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="END_EPA" & tabBO$Annee=="2013"])/
    (som_EPA13),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="END_EPA" & tabBO$Annee=="2014"])/
    (som_EPA14),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="END_EPA" & tabBO$Annee=="2015"])/
    (som_EPA15),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="END_EPA" & tabBO$Annee=="2016"])/
    (som_EPA16),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="END_EPA" & tabBO$Annee=="2017"])/
    (som_EPA17),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="END_EPA" & tabBO$Annee=="2018"])/
    (som_EPA18))

vectLIEC_X_EPA=c(
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="LIEC_X_EPA" & tabBO$Annee=="2012"])/
    (som_EPA12),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="LIEC_X_EPA" & tabBO$Annee=="2013"])/
    (som_EPA13),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="LIEC_X_EPA" & tabBO$Annee=="2014"])/
    (som_EPA14),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="LIEC_X_EPA" & tabBO$Annee=="2015"])/
    (som_EPA15),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="LIEC_X_EPA" & tabBO$Annee=="2016"])/
    (som_EPA16),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="LIEC_X_EPA" & tabBO$Annee=="2017"])/
    (som_EPA17),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="LIEC_X_EPA" & tabBO$Annee=="2018"])/
    (som_EPA18))




##### taux pour les _ANS #####


# calcule des sommes

som_ANStt=sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPI_ANS"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPA_ANS"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="ANS_ANS"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="END_ANS"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="LIEC_X_ANS"])

som_ANS12=sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPI_ANS" & tabBO$Annee=="2012"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPA_ANS" & tabBO$Annee=="2012"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="ANS_ANS" & tabBO$Annee=="2012"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="END_ANS" & tabBO$Annee=="2012"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="LIEC_X_ANS" & tabBO$Annee=="2012"])

som_ANS13=sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPI_ANS" & tabBO$Annee=="2013"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPA_ANS" & tabBO$Annee=="2013"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="ANS_ANS" & tabBO$Annee=="2013"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="END_ANS" & tabBO$Annee=="2013"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="LIEC_X_ANS" & tabBO$Annee=="2013"])

som_ANS14=sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPI_ANS" & tabBO$Annee=="2014"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPA_ANS" & tabBO$Annee=="2014"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="ANS_ANS" & tabBO$Annee=="2014"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="END_ANS" & tabBO$Annee=="2014"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="LIEC_X_ANS" & tabBO$Annee=="2014"])

som_ANS15=sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPI_ANS" & tabBO$Annee=="2015"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPA_ANS" & tabBO$Annee=="2015"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="ANS_ANS" & tabBO$Annee=="2015"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="END_ANS" & tabBO$Annee=="2015"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="LIEC_X_ANS" & tabBO$Annee=="2015"])

som_ANS16=sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPI_ANS" & tabBO$Annee=="2016"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPA_ANS" & tabBO$Annee=="2016"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="ANS_ANS" & tabBO$Annee=="2016"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="END_ANS" & tabBO$Annee=="2016"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="LIEC_X_ANS" & tabBO$Annee=="2016"])

som_ANS17=sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPI_ANS" & tabBO$Annee=="2017"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPA_ANS" & tabBO$Annee=="2017"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="ANS_ANS" & tabBO$Annee=="2017"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="END_ANS" & tabBO$Annee=="2017"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="LIEC_X_ANS" & tabBO$Annee=="2017"])

som_ANS18=sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPI_ANS" & tabBO$Annee=="2018"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPA_ANS" & tabBO$Annee=="2018"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="ANS_ANS" & tabBO$Annee=="2018"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="END_ANS" & tabBO$Annee=="2018"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="LIEC_X_ANS" & tabBO$Annee=="2018"])






# vecteurs pour chaque annee

vectEPI_ANS=c(
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPI_ANS" & tabBO$Annee=="2012"])/
    (som_ANS12),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPI_ANS" & tabBO$Annee=="2013"])/
    (som_ANS13),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPI_ANS" & tabBO$Annee=="2014"])/
    (som_ANS14),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPI_ANS" & tabBO$Annee=="2015"])/
    (som_ANS15),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPI_ANS" & tabBO$Annee=="2016"])/
    (som_ANS16),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPI_ANS" & tabBO$Annee=="2017"])/
    (som_ANS17),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPI_ANS" & tabBO$Annee=="2018"])/
    (som_ANS18))

vectEPA_ANS=c(
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPA_ANS" & tabBO$Annee=="2012"])/
    (som_ANS12),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPA_ANS" & tabBO$Annee=="2013"])/
    (som_ANS13),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPA_ANS" & tabBO$Annee=="2014"])/
    (som_ANS14),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPA_ANS" & tabBO$Annee=="2015"])/
    (som_ANS15),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPA_ANS" & tabBO$Annee=="2016"])/
    (som_ANS16),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPA_ANS" & tabBO$Annee=="2017"])/
    (som_ANS17),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPA_ANS" & tabBO$Annee=="2018"])/
    (som_ANS18))

vectANS_ANS=c(
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="ANS_ANS" & tabBO$Annee=="2012"])/
    (som_ANS12),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="ANS_ANS" & tabBO$Annee=="2013"])/
    (som_ANS13),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="ANS_ANS" & tabBO$Annee=="2014"])/
    (som_ANS14),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="ANS_ANS" & tabBO$Annee=="2015"])/
    (som_ANS15),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="ANS_ANS" & tabBO$Annee=="2016"])/
    (som_ANS16),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="ANS_ANS" & tabBO$Annee=="2017"])/
    (som_ANS17),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="ANS_ANS" & tabBO$Annee=="2018"])/
    (som_ANS18))

vectEND_ANS=c(
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="END_ANS" & tabBO$Annee=="2012"])/
    (som_ANS12),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="END_ANS" & tabBO$Annee=="2013"])/
    (som_ANS13),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="END_ANS" & tabBO$Annee=="2014"])/
    (som_ANS14),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="END_ANS" & tabBO$Annee=="2015"])/
    (som_ANS15),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="END_ANS" & tabBO$Annee=="2016"])/
    (som_ANS16),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="END_ANS" & tabBO$Annee=="2017"])/
    (som_ANS17),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="END_ANS" & tabBO$Annee=="2018"])/
    (som_ANS18))

vectLIEC_X_ANS=c(
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="LIEC_X_ANS" & tabBO$Annee=="2012"])/
    (som_ANS12),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="LIEC_X_ANS" & tabBO$Annee=="2013"])/
    (som_ANS13),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="LIEC_X_ANS" & tabBO$Annee=="2014"])/
    (som_ANS14),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="LIEC_X_ANS" & tabBO$Annee=="2015"])/
    (som_ANS15),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="LIEC_X_ANS" & tabBO$Annee=="2016"])/
    (som_ANS16),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="LIEC_X_ANS" & tabBO$Annee=="2017"])/
    (som_ANS17),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="LIEC_X_ANS" & tabBO$Annee=="2018"])/
    (som_ANS18))




##### taux pour les _END #####


# calcule des sommes

som_ENDtt=sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPI_END"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPA_END"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="ANS_END"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="END_END"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="LIEC_X_END"])

som_END12=sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPI_END" & tabBO$Annee=="2012"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPA_END" & tabBO$Annee=="2012"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="ANS_END" & tabBO$Annee=="2012"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="END_END" & tabBO$Annee=="2012"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="LIEC_X_END" & tabBO$Annee=="2012"])

som_END13=sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPI_END" & tabBO$Annee=="2013"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPA_END" & tabBO$Annee=="2013"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="ANS_END" & tabBO$Annee=="2013"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="END_END" & tabBO$Annee=="2013"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="LIEC_X_END" & tabBO$Annee=="2013"])

som_END14=sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPI_END" & tabBO$Annee=="2014"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPA_END" & tabBO$Annee=="2014"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="ANS_END" & tabBO$Annee=="2014"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="END_END" & tabBO$Annee=="2014"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="LIEC_X_END" & tabBO$Annee=="2014"])

som_END15=sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPI_END" & tabBO$Annee=="2015"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPA_END" & tabBO$Annee=="2015"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="ANS_END" & tabBO$Annee=="2015"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="END_END" & tabBO$Annee=="2015"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="LIEC_X_END" & tabBO$Annee=="2015"])

som_END16=sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPI_END" & tabBO$Annee=="2016"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPA_END" & tabBO$Annee=="2016"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="ANS_END" & tabBO$Annee=="2016"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="END_END" & tabBO$Annee=="2016"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="LIEC_X_END" & tabBO$Annee=="2016"])

som_END17=sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPI_END" & tabBO$Annee=="2017"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPA_END" & tabBO$Annee=="2017"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="ANS_END" & tabBO$Annee=="2017"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="END_END" & tabBO$Annee=="2017"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="LIEC_X_END" & tabBO$Annee=="2017"])

som_END18=sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPI_END" & tabBO$Annee=="2018"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPA_END" & tabBO$Annee=="2018"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="ANS_END" & tabBO$Annee=="2018"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="END_END" & tabBO$Annee=="2018"])+
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="LIEC_X_END" & tabBO$Annee=="2018"])






# vecteurs pour chaque annee

vectEPI_END=c(
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPI_END" & tabBO$Annee=="2012"])/
    (som_END12),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPI_END" & tabBO$Annee=="2013"])/
    (som_END13),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPI_END" & tabBO$Annee=="2014"])/
    (som_END14),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPI_END" & tabBO$Annee=="2015"])/
    (som_END15),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPI_END" & tabBO$Annee=="2016"])/
    (som_END16),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPI_END" & tabBO$Annee=="2017"])/
    (som_END17),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPI_END" & tabBO$Annee=="2018"])/
    (som_END18))

vectEPA_END=c(
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPA_END" & tabBO$Annee=="2012"])/
    (som_END12),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPA_END" & tabBO$Annee=="2013"])/
    (som_END13),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPA_END" & tabBO$Annee=="2014"])/
    (som_END14),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPA_END" & tabBO$Annee=="2015"])/
    (som_END15),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPA_END" & tabBO$Annee=="2016"])/
    (som_END16),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPA_END" & tabBO$Annee=="2017"])/
    (som_END17),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="EPA_END" & tabBO$Annee=="2018"])/
    (som_END18))

vectANS_END=c(
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="ANS_END" & tabBO$Annee=="2012"])/
    (som_END12),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="ANS_END" & tabBO$Annee=="2013"])/
    (som_END13),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="ANS_END" & tabBO$Annee=="2014"])/
    (som_END14),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="ANS_END" & tabBO$Annee=="2015"])/
    (som_END15),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="ANS_END" & tabBO$Annee=="2016"])/
    (som_END16),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="ANS_END" & tabBO$Annee=="2017"])/
    (som_END17),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="ANS_END" & tabBO$Annee=="2018"])/
    (som_END18))

vectEND_END=c(
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="END_END" & tabBO$Annee=="2012"])/
    (som_END12),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="END_END" & tabBO$Annee=="2013"])/
    (som_END13),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="END_END" & tabBO$Annee=="2014"])/
    (som_END14),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="END_END" & tabBO$Annee=="2015"])/
    (som_END15),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="END_END" & tabBO$Annee=="2016"])/
    (som_END16),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="END_END" & tabBO$Annee=="2017"])/
    (som_END17),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="END_END" & tabBO$Annee=="2018"])/
    (som_END18))

vectLIEC_X_END=c(
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="LIEC_X_END" & tabBO$Annee=="2012"])/
    (som_END12),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="LIEC_X_END" & tabBO$Annee=="2013"])/
    (som_END13),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="LIEC_X_END" & tabBO$Annee=="2014"])/
    (som_END14),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="LIEC_X_END" & tabBO$Annee=="2015"])/
    (som_END15),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="LIEC_X_END" & tabBO$Annee=="2016"])/
    (som_END16),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="LIEC_X_END" & tabBO$Annee=="2017"])/
    (som_END17),
  sum(tabBO$Nbr_EW[tabBO$FIEC_LIEC=="LIEC_X_END" & tabBO$Annee=="2018"])/
    (som_END18))



##### Graphiques #####




plot(2013:2018,vectEPI_EPI[2:7]*100,col="#000099",type="l",lwd=4,ylim=c(0,131),
     main="Taux d'identification terrain r?ussi Bourgogne",
     xlab="",
     ylab="",
     yaxt='n')
lines(2013:2018,vectEPA_EPA[2:7]*100,col="#FFFF00",type="l",lwd=4)
lines(2013:2018,vectANS_ANS[2:7]*100,col="#009900",type="l",lwd=4)
lines(2013:2018,vectEND_END[2:7]*100,col="#FF6699",type="l",lwd=4)
lines(2013:2018,vectReussiteBO*100,col="black",type="l",lwd=5)
axis(side = 2, at=c(0,50,100),las=1,labels=c("0%","50%","100%"))
axis(side = 1, at=2013:2018, pos=-15, tick=F, labels=c(paste("n=",n2013,sep=""),
                                                       paste("n=",n2014,sep=""),paste("n=",n2015,sep=""),
                                                       paste("n=",n2016,sep=""),paste("n=",n2017,sep=""),
                                                       paste("n=",n2018,sep="")))
legend("topleft",
       legend=c("EPI","EPA","ANS","END"),
       fill=c("#000099","#FFFF00","#009900","#FF6699"))


### Appartenances Reelles ###
{
  
  # erreures EPI
  plot(2013:2018,vectEPI_EPI[2:7]*100,col="#000099",type="l",lwd=5,ylim=c(0,131),
       main="Appartenance r?elle des Epig?s identifi?s sur le terrain Bourgogne",
       xlab="",
       ylab="",
       yaxt='n')
  lines(2013:2018,vectEPI_EPA[2:7]*100,col="#FFFF00",type="l",lwd=5)
  lines(2013:2018,vectEPI_ANS[2:7]*100,col="#009900",type="l",lwd=5)
  lines(2013:2018,vectEPI_END[2:7]*100,col="#FF6699",type="l",lwd=5)
  #lines(2013:2018,vectEPI_LIEC_X*100,col="grey",type="l",lwd=5)
  axis(side = 2, at=c(0,50,100),las=1,labels=c("0%","50%","100%"))
  axis(side = 1, at=2013:2018, pos=-15, tick=F, labels=c(paste("n=",n2013,sep=""),
                                                         paste("n=",n2014,sep=""),paste("n=",n2015,sep=""),
                                                         paste("n=",n2016,sep=""),paste("n=",n2017,sep=""),
                                                         paste("n=",n2018,sep="")))
  legend("topleft",
         legend=c("EPI","EPA","ANS","END"),
         fill=c("#000099","#FFFF00","#009900","#FF6699"))
  
  
  
  # erreures EPA
  plot(2013:2018,vectEPA_EPI[2:7]*100,col="#000099",type="l",lwd=5,ylim=c(0,131),
       main="Appartenance r?elle des Epi-An?ciques identifi?s sur le terrain Bourgogne",
       xlab="",
       ylab="",
       las=1,
       yaxt='n')
  lines(2013:2018,vectEPA_EPA[2:7]*100,col="#FFFF00",type="l",lwd=5)
  lines(2013:2018,vectEPA_ANS[2:7]*100,col="#009900",type="l",lwd=5)
  lines(2013:2018,vectEPA_END[2:7]*100,col="#FF6699",type="l",lwd=5)
  #lines(2013:2018,vectEPA_LIEC_X*100,col="grey",type="l",lwd=5)
  axis(side = 2, at=c(0,50,100),las=1,labels=c("0%","50%","100%"))
  axis(side = 1, at=2013:2018, pos=-15, tick=F, labels=c(paste("n=",n2013,sep=""),
                                                         paste("n=",n2014,sep=""),paste("n=",n2015,sep=""),
                                                         paste("n=",n2016,sep=""),paste("n=",n2017,sep=""),
                                                         paste("n=",n2018,sep="")))
  legend("topleft",
         legend=c("EPI","EPA","ANS","END"),
         fill=c("#000099","#FFFF00","#009900","#FF6699"))
  
  
  # erreures ANS
  plot(2013:2018,vectANS_EPI[2:7]*100,col="#000099",type="l",lwd=5,ylim=c(0,131),
       main="Appartenance r?elle des An?ciques-Stricts identifi?s sur le terrain Bourgogne",
       xlab="",
       ylab="",
       las=1,
       yaxt='n')
  lines(2013:2018,vectANS_EPA[2:7]*100,col="#FFFF00",type="l",lwd=5)
  lines(2013:2018,vectANS_ANS[2:7]*100,col="#009900",type="l",lwd=5)
  lines(2013:2018,vectANS_END[2:7]*100,col="#FF6699",type="l",lwd=5)
  #lines(2013:2018,vectANS_LIEC_X*100,col="grey",type="l",lwd=5)
  axis(side = 2, at=c(0,50,100),las=1,labels=c("0%","50%","100%"))
  axis(side = 1, at=2013:2018, pos=-15, tick=F, labels=c(paste("n=",n2013,sep=""),
                                                         paste("n=",n2014,sep=""),paste("n=",n2015,sep=""),
                                                         paste("n=",n2016,sep=""),paste("n=",n2017,sep=""),
                                                         paste("n=",n2018,sep="")))
  legend("topleft",
         legend=c("EPI","EPA","ANS","END"),
         fill=c("#000099","#FFFF00","#009900","#FF6699"))
  
  
  # erreures END
  plot(2013:2018,vectEND_EPI[2:7]*100,col="#000099",type="l",lwd=5,ylim=c(0,131),
       main="Appartenance r?elle des Endog?s identifi?s sur le terrain Bourgogne",
       xlab="",
       ylab="",
       las=1,
       yaxt='n')
  lines(2013:2018,vectEND_EPA[2:7]*100,col="#FFFF00",type="l",lwd=5)
  lines(2013:2018,vectEND_ANS[2:7]*100,col="#009900",type="l",lwd=5)
  lines(2013:2018,vectEND_END[2:7]*100,col="#FF6699",type="l",lwd=5)
  #lines(2013:2018,vectEND_LIEC_X*100,col="grey",type="l",lwd=5)
  axis(side = 2, at=c(0,50,100),las=1,labels=c("0%","50%","100%"))
  axis(side = 1, at=2013:2018, pos=-15, tick=F, labels=c(paste("n=",n2013,sep=""),
                                                         paste("n=",n2014,sep=""),paste("n=",n2015,sep=""),
                                                         paste("n=",n2016,sep=""),paste("n=",n2017,sep=""),
                                                         paste("n=",n2018,sep="")))
  legend("topleft",
         legend=c("EPI","EPA","ANS","END"),
         fill=c("#000099","#FFFF00","#009900","#FF6699"))
  
}


### Erreurs par LIEC ###
{
  
  # erreures EPI
  plot(2013:2018,vectEPI_EPI[2:7]*100,col="#000099",type="l",lwd=5,ylim=c(0,141),
       main="Erreur d'identification terrain des Epig?s Bourgogne",
       xlab="",
       ylab="",
       yaxt='n')
  lines(2013:2018,vectEPA_EPI[2:7]*100,col="#FFFF00",type="l",lwd=5)
  lines(2013:2018,vectANS_EPI[2:7]*100,col="#009900",type="l",lwd=5)
  lines(2013:2018,vectEND_EPI[2:7]*100,col="#FF6699",type="l",lwd=5)
  lines(2013:2018,vectLIEC_X_EPI[2:7]*100,col="grey",type="l",lwd=5)
  axis(side = 2, at=c(0,50,100),las=1,labels=c("0%","50%","100%"))
  axis(side = 1, at=2013:2018, pos=-15, tick=F, labels=c(paste("n=",n2013,sep=""),
                                                         paste("n=",n2014,sep=""),paste("n=",n2015,sep=""),
                                                         paste("n=",n2016,sep=""),paste("n=",n2017,sep=""),
                                                         paste("n=",n2018,sep="")))
  legend("topleft",
         legend=c("EPI","EPA","ANS","END","LIEC_X"),
         fill=c("#000099","#FFFF00","#009900","#FF6699","grey"))
  
  
  
  # erreures EPA
  plot(2013:2018,vectEPI_EPA[2:7]*100,col="#000099",type="l",lwd=5,ylim=c(0,141),
       main="Erreur d'identification terrain des Epi-An?ciques Bourgogne",
       xlab="",
       ylab="",
       las=1,
       yaxt='n')
  lines(2013:2018,vectEPA_EPA[2:7]*100,col="#FFFF00",type="l",lwd=5)
  lines(2013:2018,vectANS_EPA[2:7]*100,col="#009900",type="l",lwd=5)
  lines(2013:2018,vectEND_EPA[2:7]*100,col="#FF6699",type="l",lwd=5)
  lines(2013:2018,vectLIEC_X_EPA[2:7]*100,col="grey",type="l",lwd=5)
  axis(side = 2, at=c(0,50,100),las=1,labels=c("0%","50%","100%"))
  axis(side = 1, at=2013:2018, pos=-15, tick=F, labels=c(paste("n=",n2013,sep=""),
                                                         paste("n=",n2014,sep=""),paste("n=",n2015,sep=""),
                                                         paste("n=",n2016,sep=""),paste("n=",n2017,sep=""),
                                                         paste("n=",n2018,sep="")))
  legend("topleft",
         legend=c("EPI","EPA","ANS","END","LIEC_X"),
         fill=c("#000099","#FFFF00","#009900","#FF6699","grey"))
  
  
  # erreures ANS
  plot(2013:2018,vectEPI_ANS[2:7]*100,col="#000099",type="l",lwd=5,ylim=c(0,141),
       main="Erreur d'identification terrain des An?ciques-Stricts Bourgogne",
       xlab="",
       ylab="",
       las=1,
       yaxt='n')
  lines(2013:2018,vectEPA_ANS[2:7]*100,col="#FFFF00",type="l",lwd=5)
  lines(2013:2018,vectANS_ANS[2:7]*100,col="#009900",type="l",lwd=5)
  lines(2013:2018,vectEND_ANS[2:7]*100,col="#FF6699",type="l",lwd=5)
  lines(2013:2018,vectLIEC_X_ANS[2:7]*100,col="grey",type="l",lwd=5)
  axis(side = 2, at=c(0,50,100),las=1,labels=c("0%","50%","100%"))
  axis(side = 1, at=2013:2018, pos=-15, tick=F, labels=c(paste("n=",n2013,sep=""),
                                                         paste("n=",n2014,sep=""),paste("n=",n2015,sep=""),
                                                         paste("n=",n2016,sep=""),paste("n=",n2017,sep=""),
                                                         paste("n=",n2018,sep="")))
  legend("topleft",
         legend=c("EPI","EPA","ANS","END","LIEC_X"),
         fill=c("#000099","#FFFF00","#009900","#FF6699","grey"))
  
  
  # erreures END
  plot(2013:2018,vectEPI_END[2:7]*100,col="#000099",type="l",lwd=5,ylim=c(0,141),
       main="Erreur d'identification terrain des Endog?s Bourgogne",
       xlab="",
       ylab="",
       las=1,
       yaxt='n')
  lines(2013:2018,vectEPA_END[2:7]*100,col="#FFFF00",type="l",lwd=5)
  lines(2013:2018,vectANS_END[2:7]*100,col="#009900",type="l",lwd=5)
  lines(2013:2018,vectEND_END[2:7]*100,col="#FF6699",type="l",lwd=5)
  lines(2013:2018,vectLIEC_X_END[2:7]*100,col="grey",type="l",lwd=5)
  axis(side = 2, at=c(0,50,100),las=1,labels=c("0%","50%","100%"))
  axis(side = 1, at=2013:2018, pos=-15, tick=F, labels=c(paste("n=",n2013,sep=""),
                                                         paste("n=",n2014,sep=""),paste("n=",n2015,sep=""),
                                                         paste("n=",n2016,sep=""),paste("n=",n2017,sep=""),
                                                         paste("n=",n2018,sep="")))
  legend("topleft",
         legend=c("EPI","EPA","ANS","END","LIEC_X"),
         fill=c("#000099","#FFFF00","#009900","#FF6699","grey"))
  
}








##### Region mal formee : Franche-Compt? #####

tabFC=Tab2[Tab2$Region%in%"FC",]
nrow(tabFC)==1662
dataID <- unique(sbt_m[,c("cle",colID,"Region")])
nrow(dataID)==1109
dataID=dataID[dataID$Region%in%"FC",]
nrow(dataID)==57

n2013=table(dataID$Annee)[1]
n2014=table(dataID$Annee)[2]
n2015=table(dataID$Annee)[3]
n2016=table(dataID$Annee)[4]
n2017=table(dataID$Annee)[5]
n2018=table(dataID$Annee)[6]


vectReussiteFC=c(
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC%in%c("EPI_EPI","EPA_EPA","ANS_ANS","END_END") & tabFC$Annee%in%2013])/sum(tabFC$Nbr_EW[tabFC$Annee%in%2013]),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC%in%c("EPI_EPI","EPA_EPA","ANS_ANS","END_END") & tabFC$Annee%in%2014])/sum(tabFC$Nbr_EW[tabFC$Annee%in%2014]),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC%in%c("EPI_EPI","EPA_EPA","ANS_ANS","END_END") & tabFC$Annee%in%2015])/sum(tabFC$Nbr_EW[tabFC$Annee%in%2015]),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC%in%c("EPI_EPI","EPA_EPA","ANS_ANS","END_END") & tabFC$Annee%in%2016])/sum(tabFC$Nbr_EW[tabFC$Annee%in%2016]),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC%in%c("EPI_EPI","EPA_EPA","ANS_ANS","END_END") & tabFC$Annee%in%2017])/sum(tabFC$Nbr_EW[tabFC$Annee%in%2017]),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC%in%c("EPI_EPI","EPA_EPA","ANS_ANS","END_END") & tabFC$Annee%in%2018])/sum(tabFC$Nbr_EW[tabFC$Annee%in%2018]))



##### taux pour les _EPI #####


# calcule des sommes

som_EPItt=sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPI_EPI"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPA_EPI"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="ANS_EPI"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="END_EPI"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="LIEC_X_EPI"])

som_EPI12=sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPI_EPI" & tabFC$Annee=="2012"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPA_EPI" & tabFC$Annee=="2012"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="ANS_EPI" & tabFC$Annee=="2012"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="END_EPI" & tabFC$Annee=="2012"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="LIEC_X_EPI" & tabFC$Annee=="2012"])

som_EPI13=sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPI_EPI" & tabFC$Annee=="2013"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPA_EPI" & tabFC$Annee=="2013"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="ANS_EPI" & tabFC$Annee=="2013"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="END_EPI" & tabFC$Annee=="2013"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="LIEC_X_EPI" & tabFC$Annee=="2013"])

som_EPI14=sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPI_EPI" & tabFC$Annee=="2014"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPA_EPI" & tabFC$Annee=="2014"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="ANS_EPI" & tabFC$Annee=="2014"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="END_EPI" & tabFC$Annee=="2014"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="LIEC_X_EPI" & tabFC$Annee=="2014"])

som_EPI15=sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPI_EPI" & tabFC$Annee=="2015"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPA_EPI" & tabFC$Annee=="2015"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="ANS_EPI" & tabFC$Annee=="2015"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="END_EPI" & tabFC$Annee=="2015"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="LIEC_X_EPI" & tabFC$Annee=="2015"])

som_EPI16=sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPI_EPI" & tabFC$Annee=="2016"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPA_EPI" & tabFC$Annee=="2016"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="ANS_EPI" & tabFC$Annee=="2016"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="END_EPI" & tabFC$Annee=="2016"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="LIEC_X_EPI" & tabFC$Annee=="2016"])

som_EPI17=sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPI_EPI" & tabFC$Annee=="2017"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPA_EPI" & tabFC$Annee=="2017"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="ANS_EPI" & tabFC$Annee=="2017"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="END_EPI" & tabFC$Annee=="2017"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="LIEC_X_EPI" & tabFC$Annee=="2017"])

som_EPI18=sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPI_EPI" & tabFC$Annee=="2018"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPA_EPI" & tabFC$Annee=="2018"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="ANS_EPI" & tabFC$Annee=="2018"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="END_EPI" & tabFC$Annee=="2018"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="LIEC_X_EPI" & tabFC$Annee=="2018"])




# vecteurs pour chaque annee

vectEPI_EPI=c(
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPI_EPI" & tabFC$Annee=="2012"])/
    (som_EPI12),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPI_EPI" & tabFC$Annee=="2013"])/
    (som_EPI13),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPI_EPI" & tabFC$Annee=="2014"])/
    (som_EPI14),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPI_EPI" & tabFC$Annee=="2015"])/
    (som_EPI15),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPI_EPI" & tabFC$Annee=="2016"])/
    (som_EPI16),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPI_EPI" & tabFC$Annee=="2017"])/
    (som_EPI17),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPI_EPI" & tabFC$Annee=="2018"])/
    (som_EPI18))

vectEPA_EPI=c(
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPA_EPI" & tabFC$Annee=="2012"])/
    (som_EPI12),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPA_EPI" & tabFC$Annee=="2013"])/
    (som_EPI13),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPA_EPI" & tabFC$Annee=="2014"])/
    (som_EPI14),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPA_EPI" & tabFC$Annee=="2015"])/
    (som_EPI15),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPA_EPI" & tabFC$Annee=="2016"])/
    (som_EPI16),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPA_EPI" & tabFC$Annee=="2017"])/
    (som_EPI17),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPA_EPI" & tabFC$Annee=="2018"])/
    (som_EPI18))

vectANS_EPI=c(
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="ANS_EPI" & tabFC$Annee=="2012"])/
    (som_EPI12),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="ANS_EPI" & tabFC$Annee=="2013"])/
    (som_EPI13),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="ANS_EPI" & tabFC$Annee=="2014"])/
    (som_EPI14),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="ANS_EPI" & tabFC$Annee=="2015"])/
    (som_EPI15),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="ANS_EPI" & tabFC$Annee=="2016"])/
    (som_EPI16),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="ANS_EPI" & tabFC$Annee=="2017"])/
    (som_EPI17),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="ANS_EPI" & tabFC$Annee=="2018"])/
    (som_EPI18))

vectEND_EPI=c(
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="END_EPI" & tabFC$Annee=="2012"])/
    (som_EPI12),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="END_EPI" & tabFC$Annee=="2013"])/
    (som_EPI13),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="END_EPI" & tabFC$Annee=="2014"])/
    (som_EPI14),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="END_EPI" & tabFC$Annee=="2015"])/
    (som_EPI15),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="END_EPI" & tabFC$Annee=="2016"])/
    (som_EPI16),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="END_EPI" & tabFC$Annee=="2017"])/
    (som_EPI17),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="END_EPI" & tabFC$Annee=="2018"])/
    (som_EPI18))

vectLIEC_X_EPI=c(
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="LIEC_X_EPI" & tabFC$Annee=="2012"])/
    (som_EPI12),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="LIEC_X_EPI" & tabFC$Annee=="2013"])/
    (som_EPI13),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="LIEC_X_EPI" & tabFC$Annee=="2014"])/
    (som_EPI14),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="LIEC_X_EPI" & tabFC$Annee=="2015"])/
    (som_EPI15),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="LIEC_X_EPI" & tabFC$Annee=="2016"])/
    (som_EPI16),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="LIEC_X_EPI" & tabFC$Annee=="2017"])/
    (som_EPI17),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="LIEC_X_EPI" & tabFC$Annee=="2018"])/
    (som_EPI18))




##### taux pour les _EPA #####


# calcule des sommes

som_EPAtt=sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPI_EPA"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPA_EPA"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="ANS_EPA"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="END_EPA"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="LIEC_X_EPA"])

som_EPA12=sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPI_EPA" & tabFC$Annee=="2012"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPA_EPA" & tabFC$Annee=="2012"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="ANS_EPA" & tabFC$Annee=="2012"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="END_EPA" & tabFC$Annee=="2012"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="LIEC_X_EPA" & tabFC$Annee=="2012"])

som_EPA13=sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPI_EPA" & tabFC$Annee=="2013"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPA_EPA" & tabFC$Annee=="2013"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="ANS_EPA" & tabFC$Annee=="2013"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="END_EPA" & tabFC$Annee=="2013"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="LIEC_X_EPA" & tabFC$Annee=="2013"])

som_EPA14=sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPI_EPA" & tabFC$Annee=="2014"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPA_EPA" & tabFC$Annee=="2014"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="ANS_EPA" & tabFC$Annee=="2014"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="END_EPA" & tabFC$Annee=="2014"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="LIEC_X_EPA" & tabFC$Annee=="2014"])

som_EPA15=sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPI_EPA" & tabFC$Annee=="2015"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPA_EPA" & tabFC$Annee=="2015"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="ANS_EPA" & tabFC$Annee=="2015"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="END_EPA" & tabFC$Annee=="2015"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="LIEC_X_EPA" & tabFC$Annee=="2015"])

som_EPA16=sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPI_EPA" & tabFC$Annee=="2016"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPA_EPA" & tabFC$Annee=="2016"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="ANS_EPA" & tabFC$Annee=="2016"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="END_EPA" & tabFC$Annee=="2016"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="LIEC_X_EPA" & tabFC$Annee=="2016"])

som_EPA17=sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPI_EPA" & tabFC$Annee=="2017"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPA_EPA" & tabFC$Annee=="2017"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="ANS_EPA" & tabFC$Annee=="2017"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="END_EPA" & tabFC$Annee=="2017"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="LIEC_X_EPA" & tabFC$Annee=="2017"])

som_EPA18=sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPI_EPA" & tabFC$Annee=="2018"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPA_EPA" & tabFC$Annee=="2018"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="ANS_EPA" & tabFC$Annee=="2018"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="END_EPA" & tabFC$Annee=="2018"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="LIEC_X_EPA" & tabFC$Annee=="2018"])






# vecteurs pour chaque annee

vectEPI_EPA=c(
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPI_EPA" & tabFC$Annee=="2012"])/
    (som_EPA12),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPI_EPA" & tabFC$Annee=="2013"])/
    (som_EPA13),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPI_EPA" & tabFC$Annee=="2014"])/
    (som_EPA14),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPI_EPA" & tabFC$Annee=="2015"])/
    (som_EPA15),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPI_EPA" & tabFC$Annee=="2016"])/
    (som_EPA16),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPI_EPA" & tabFC$Annee=="2017"])/
    (som_EPA17),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPI_EPA" & tabFC$Annee=="2018"])/
    (som_EPA18))

vectEPA_EPA=c(
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPA_EPA" & tabFC$Annee=="2012"])/
    (som_EPA12),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPA_EPA" & tabFC$Annee=="2013"])/
    (som_EPA13),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPA_EPA" & tabFC$Annee=="2014"])/
    (som_EPA14),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPA_EPA" & tabFC$Annee=="2015"])/
    (som_EPA15),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPA_EPA" & tabFC$Annee=="2016"])/
    (som_EPA16),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPA_EPA" & tabFC$Annee=="2017"])/
    (som_EPA17),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPA_EPA" & tabFC$Annee=="2018"])/
    (som_EPA18))

vectANS_EPA=c(
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="ANS_EPA" & tabFC$Annee=="2012"])/
    (som_EPA12),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="ANS_EPA" & tabFC$Annee=="2013"])/
    (som_EPA13),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="ANS_EPA" & tabFC$Annee=="2014"])/
    (som_EPA14),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="ANS_EPA" & tabFC$Annee=="2015"])/
    (som_EPA15),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="ANS_EPA" & tabFC$Annee=="2016"])/
    (som_EPA16),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="ANS_EPA" & tabFC$Annee=="2017"])/
    (som_EPA17),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="ANS_EPA" & tabFC$Annee=="2018"])/
    (som_EPA18))

vectEND_EPA=c(
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="END_EPA" & tabFC$Annee=="2012"])/
    (som_EPA12),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="END_EPA" & tabFC$Annee=="2013"])/
    (som_EPA13),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="END_EPA" & tabFC$Annee=="2014"])/
    (som_EPA14),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="END_EPA" & tabFC$Annee=="2015"])/
    (som_EPA15),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="END_EPA" & tabFC$Annee=="2016"])/
    (som_EPA16),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="END_EPA" & tabFC$Annee=="2017"])/
    (som_EPA17),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="END_EPA" & tabFC$Annee=="2018"])/
    (som_EPA18))

vectLIEC_X_EPA=c(
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="LIEC_X_EPA" & tabFC$Annee=="2012"])/
    (som_EPA12),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="LIEC_X_EPA" & tabFC$Annee=="2013"])/
    (som_EPA13),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="LIEC_X_EPA" & tabFC$Annee=="2014"])/
    (som_EPA14),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="LIEC_X_EPA" & tabFC$Annee=="2015"])/
    (som_EPA15),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="LIEC_X_EPA" & tabFC$Annee=="2016"])/
    (som_EPA16),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="LIEC_X_EPA" & tabFC$Annee=="2017"])/
    (som_EPA17),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="LIEC_X_EPA" & tabFC$Annee=="2018"])/
    (som_EPA18))




##### taux pour les _ANS #####


# calcule des sommes

som_ANStt=sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPI_ANS"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPA_ANS"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="ANS_ANS"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="END_ANS"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="LIEC_X_ANS"])

som_ANS12=sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPI_ANS" & tabFC$Annee=="2012"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPA_ANS" & tabFC$Annee=="2012"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="ANS_ANS" & tabFC$Annee=="2012"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="END_ANS" & tabFC$Annee=="2012"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="LIEC_X_ANS" & tabFC$Annee=="2012"])

som_ANS13=sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPI_ANS" & tabFC$Annee=="2013"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPA_ANS" & tabFC$Annee=="2013"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="ANS_ANS" & tabFC$Annee=="2013"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="END_ANS" & tabFC$Annee=="2013"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="LIEC_X_ANS" & tabFC$Annee=="2013"])

som_ANS14=sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPI_ANS" & tabFC$Annee=="2014"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPA_ANS" & tabFC$Annee=="2014"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="ANS_ANS" & tabFC$Annee=="2014"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="END_ANS" & tabFC$Annee=="2014"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="LIEC_X_ANS" & tabFC$Annee=="2014"])

som_ANS15=sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPI_ANS" & tabFC$Annee=="2015"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPA_ANS" & tabFC$Annee=="2015"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="ANS_ANS" & tabFC$Annee=="2015"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="END_ANS" & tabFC$Annee=="2015"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="LIEC_X_ANS" & tabFC$Annee=="2015"])

som_ANS16=sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPI_ANS" & tabFC$Annee=="2016"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPA_ANS" & tabFC$Annee=="2016"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="ANS_ANS" & tabFC$Annee=="2016"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="END_ANS" & tabFC$Annee=="2016"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="LIEC_X_ANS" & tabFC$Annee=="2016"])

som_ANS17=sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPI_ANS" & tabFC$Annee=="2017"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPA_ANS" & tabFC$Annee=="2017"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="ANS_ANS" & tabFC$Annee=="2017"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="END_ANS" & tabFC$Annee=="2017"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="LIEC_X_ANS" & tabFC$Annee=="2017"])

som_ANS18=sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPI_ANS" & tabFC$Annee=="2018"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPA_ANS" & tabFC$Annee=="2018"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="ANS_ANS" & tabFC$Annee=="2018"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="END_ANS" & tabFC$Annee=="2018"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="LIEC_X_ANS" & tabFC$Annee=="2018"])






# vecteurs pour chaque annee

vectEPI_ANS=c(
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPI_ANS" & tabFC$Annee=="2012"])/
    (som_ANS12),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPI_ANS" & tabFC$Annee=="2013"])/
    (som_ANS13),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPI_ANS" & tabFC$Annee=="2014"])/
    (som_ANS14),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPI_ANS" & tabFC$Annee=="2015"])/
    (som_ANS15),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPI_ANS" & tabFC$Annee=="2016"])/
    (som_ANS16),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPI_ANS" & tabFC$Annee=="2017"])/
    (som_ANS17),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPI_ANS" & tabFC$Annee=="2018"])/
    (som_ANS18))

vectEPA_ANS=c(
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPA_ANS" & tabFC$Annee=="2012"])/
    (som_ANS12),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPA_ANS" & tabFC$Annee=="2013"])/
    (som_ANS13),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPA_ANS" & tabFC$Annee=="2014"])/
    (som_ANS14),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPA_ANS" & tabFC$Annee=="2015"])/
    (som_ANS15),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPA_ANS" & tabFC$Annee=="2016"])/
    (som_ANS16),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPA_ANS" & tabFC$Annee=="2017"])/
    (som_ANS17),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPA_ANS" & tabFC$Annee=="2018"])/
    (som_ANS18))

vectANS_ANS=c(
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="ANS_ANS" & tabFC$Annee=="2012"])/
    (som_ANS12),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="ANS_ANS" & tabFC$Annee=="2013"])/
    (som_ANS13),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="ANS_ANS" & tabFC$Annee=="2014"])/
    (som_ANS14),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="ANS_ANS" & tabFC$Annee=="2015"])/
    (som_ANS15),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="ANS_ANS" & tabFC$Annee=="2016"])/
    (som_ANS16),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="ANS_ANS" & tabFC$Annee=="2017"])/
    (som_ANS17),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="ANS_ANS" & tabFC$Annee=="2018"])/
    (som_ANS18))

vectEND_ANS=c(
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="END_ANS" & tabFC$Annee=="2012"])/
    (som_ANS12),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="END_ANS" & tabFC$Annee=="2013"])/
    (som_ANS13),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="END_ANS" & tabFC$Annee=="2014"])/
    (som_ANS14),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="END_ANS" & tabFC$Annee=="2015"])/
    (som_ANS15),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="END_ANS" & tabFC$Annee=="2016"])/
    (som_ANS16),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="END_ANS" & tabFC$Annee=="2017"])/
    (som_ANS17),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="END_ANS" & tabFC$Annee=="2018"])/
    (som_ANS18))

vectLIEC_X_ANS=c(
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="LIEC_X_ANS" & tabFC$Annee=="2012"])/
    (som_ANS12),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="LIEC_X_ANS" & tabFC$Annee=="2013"])/
    (som_ANS13),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="LIEC_X_ANS" & tabFC$Annee=="2014"])/
    (som_ANS14),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="LIEC_X_ANS" & tabFC$Annee=="2015"])/
    (som_ANS15),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="LIEC_X_ANS" & tabFC$Annee=="2016"])/
    (som_ANS16),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="LIEC_X_ANS" & tabFC$Annee=="2017"])/
    (som_ANS17),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="LIEC_X_ANS" & tabFC$Annee=="2018"])/
    (som_ANS18))




##### taux pour les _END #####


# calcule des sommes

som_ENDtt=sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPI_END"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPA_END"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="ANS_END"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="END_END"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="LIEC_X_END"])

som_END12=sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPI_END" & tabFC$Annee=="2012"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPA_END" & tabFC$Annee=="2012"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="ANS_END" & tabFC$Annee=="2012"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="END_END" & tabFC$Annee=="2012"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="LIEC_X_END" & tabFC$Annee=="2012"])

som_END13=sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPI_END" & tabFC$Annee=="2013"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPA_END" & tabFC$Annee=="2013"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="ANS_END" & tabFC$Annee=="2013"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="END_END" & tabFC$Annee=="2013"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="LIEC_X_END" & tabFC$Annee=="2013"])

som_END14=sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPI_END" & tabFC$Annee=="2014"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPA_END" & tabFC$Annee=="2014"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="ANS_END" & tabFC$Annee=="2014"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="END_END" & tabFC$Annee=="2014"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="LIEC_X_END" & tabFC$Annee=="2014"])

som_END15=sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPI_END" & tabFC$Annee=="2015"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPA_END" & tabFC$Annee=="2015"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="ANS_END" & tabFC$Annee=="2015"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="END_END" & tabFC$Annee=="2015"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="LIEC_X_END" & tabFC$Annee=="2015"])

som_END16=sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPI_END" & tabFC$Annee=="2016"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPA_END" & tabFC$Annee=="2016"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="ANS_END" & tabFC$Annee=="2016"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="END_END" & tabFC$Annee=="2016"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="LIEC_X_END" & tabFC$Annee=="2016"])

som_END17=sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPI_END" & tabFC$Annee=="2017"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPA_END" & tabFC$Annee=="2017"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="ANS_END" & tabFC$Annee=="2017"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="END_END" & tabFC$Annee=="2017"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="LIEC_X_END" & tabFC$Annee=="2017"])

som_END18=sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPI_END" & tabFC$Annee=="2018"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPA_END" & tabFC$Annee=="2018"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="ANS_END" & tabFC$Annee=="2018"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="END_END" & tabFC$Annee=="2018"])+
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="LIEC_X_END" & tabFC$Annee=="2018"])






# vecteurs pour chaque annee

vectEPI_END=c(
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPI_END" & tabFC$Annee=="2012"])/
    (som_END12),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPI_END" & tabFC$Annee=="2013"])/
    (som_END13),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPI_END" & tabFC$Annee=="2014"])/
    (som_END14),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPI_END" & tabFC$Annee=="2015"])/
    (som_END15),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPI_END" & tabFC$Annee=="2016"])/
    (som_END16),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPI_END" & tabFC$Annee=="2017"])/
    (som_END17),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPI_END" & tabFC$Annee=="2018"])/
    (som_END18))

vectEPA_END=c(
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPA_END" & tabFC$Annee=="2012"])/
    (som_END12),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPA_END" & tabFC$Annee=="2013"])/
    (som_END13),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPA_END" & tabFC$Annee=="2014"])/
    (som_END14),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPA_END" & tabFC$Annee=="2015"])/
    (som_END15),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPA_END" & tabFC$Annee=="2016"])/
    (som_END16),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPA_END" & tabFC$Annee=="2017"])/
    (som_END17),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPA_END" & tabFC$Annee=="2018"])/
    (som_END18))

vectANS_END=c(
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="ANS_END" & tabFC$Annee=="2012"])/
    (som_END12),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="ANS_END" & tabFC$Annee=="2013"])/
    (som_END13),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="ANS_END" & tabFC$Annee=="2014"])/
    (som_END14),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="ANS_END" & tabFC$Annee=="2015"])/
    (som_END15),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="ANS_END" & tabFC$Annee=="2016"])/
    (som_END16),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="ANS_END" & tabFC$Annee=="2017"])/
    (som_END17),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="ANS_END" & tabFC$Annee=="2018"])/
    (som_END18))

vectEND_END=c(
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="END_END" & tabFC$Annee=="2012"])/
    (som_END12),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="END_END" & tabFC$Annee=="2013"])/
    (som_END13),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="END_END" & tabFC$Annee=="2014"])/
    (som_END14),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="END_END" & tabFC$Annee=="2015"])/
    (som_END15),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="END_END" & tabFC$Annee=="2016"])/
    (som_END16),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="END_END" & tabFC$Annee=="2017"])/
    (som_END17),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="END_END" & tabFC$Annee=="2018"])/
    (som_END18))

vectLIEC_X_END=c(
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="LIEC_X_END" & tabFC$Annee=="2012"])/
    (som_END12),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="LIEC_X_END" & tabFC$Annee=="2013"])/
    (som_END13),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="LIEC_X_END" & tabFC$Annee=="2014"])/
    (som_END14),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="LIEC_X_END" & tabFC$Annee=="2015"])/
    (som_END15),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="LIEC_X_END" & tabFC$Annee=="2016"])/
    (som_END16),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="LIEC_X_END" & tabFC$Annee=="2017"])/
    (som_END17),
  sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="LIEC_X_END" & tabFC$Annee=="2018"])/
    (som_END18))



##### Graphiques #####




plot(2013:2018,vectEPI_EPI[2:7]*100,col="#000099",type="l",lwd=4,ylim=c(0,131),
     main="Taux d'identification terrain r?ussi FC",
     xlab="",
     ylab="",
     yaxt='n')
lines(2013:2018,vectEPA_EPA[2:7]*100,col="#FFFF00",type="l",lwd=4)
lines(2013:2018,vectANS_ANS[2:7]*100,col="#009900",type="l",lwd=4)
lines(2013:2018,vectEND_END[2:7]*100,col="#FF6699",type="l",lwd=4)
lines(2013:2018,vectReussiteFC*100,col="black",type="l",lwd=5)
axis(side = 2, at=c(0,50,100),las=1,labels=c("0%","50%","100%"))
axis(side = 1, at=2013:2018, pos=-15, tick=F, labels=c(paste("n=",n2013,sep=""),
                                                       paste("n=",n2014,sep=""),paste("n=",n2015,sep=""),
                                                       paste("n=",n2016,sep=""),paste("n=",n2017,sep=""),
                                                       paste("n=",n2018,sep="")))
legend("topleft",
       legend=c("EPI","EPA","ANS","END"),
       fill=c("#000099","#FFFF00","#009900","#FF6699"))



yEPI=sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPI_EPI"])/som_EPItt
yEPA=sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="EPA_EPA"])/som_EPAtt
yANS=sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="ANS_ANS"])/som_ANStt
yEND=sum(tabFC$Nbr_EW[tabFC$FIEC_LIEC=="END_END"])/som_ENDtt


barplot(c(yEPI,yEPA,yANS,yEND),ylim=c(0,1))






### Appartenances Reelles ###
{
  
  # erreures EPI
  plot(2013:2018,vectEPI_EPI[2:7]*100,col="#000099",type="l",lwd=5,ylim=c(0,131),
       main="Appartenance r?elle des Epig?s identifi?s sur le terrain FC",
       xlab="",
       ylab="",
       yaxt='n')
  lines(2013:2018,vectEPI_EPA[2:7]*100,col="#FFFF00",type="l",lwd=5)
  lines(2013:2018,vectEPI_ANS[2:7]*100,col="#009900",type="l",lwd=5)
  lines(2013:2018,vectEPI_END[2:7]*100,col="#FF6699",type="l",lwd=5)
  #lines(2013:2018,vectEPI_LIEC_X*100,col="grey",type="l",lwd=5)
  axis(side = 2, at=c(0,50,100),las=1,labels=c("0%","50%","100%"))
  axis(side = 1, at=2013:2018, pos=-15, tick=F, labels=c(paste("n=",n2013,sep=""),
                                                         paste("n=",n2014,sep=""),paste("n=",n2015,sep=""),
                                                         paste("n=",n2016,sep=""),paste("n=",n2017,sep=""),
                                                         paste("n=",n2018,sep="")))
  legend("topleft",
         legend=c("EPI","EPA","ANS","END"),
         fill=c("#000099","#FFFF00","#009900","#FF6699"))
  
  
  
  # erreures EPA
  plot(2013:2018,vectEPA_EPI[2:7]*100,col="#000099",type="l",lwd=5,ylim=c(0,131),
       main="Appartenance r?elle des Epi-An?ciques identifi?s sur le terrain FC",
       xlab="",
       ylab="",
       las=1,
       yaxt='n')
  lines(2013:2018,vectEPA_EPA[2:7]*100,col="#FFFF00",type="l",lwd=5)
  lines(2013:2018,vectEPA_ANS[2:7]*100,col="#009900",type="l",lwd=5)
  lines(2013:2018,vectEPA_END[2:7]*100,col="#FF6699",type="l",lwd=5)
  #lines(2013:2018,vectEPA_LIEC_X*100,col="grey",type="l",lwd=5)
  axis(side = 2, at=c(0,50,100),las=1,labels=c("0%","50%","100%"))
  axis(side = 1, at=2013:2018, pos=-15, tick=F, labels=c(paste("n=",n2013,sep=""),
                                                         paste("n=",n2014,sep=""),paste("n=",n2015,sep=""),
                                                         paste("n=",n2016,sep=""),paste("n=",n2017,sep=""),
                                                         paste("n=",n2018,sep="")))
  legend("topleft",
         legend=c("EPI","EPA","ANS","END"),
         fill=c("#000099","#FFFF00","#009900","#FF6699"))
  
  
  # erreures ANS
  plot(2013:2018,vectANS_EPI[2:7]*100,col="#000099",type="l",lwd=5,ylim=c(0,131),
       main="Appartenance r?elle des An?ciques-Stricts identifi?s sur le terrain FC",
       xlab="",
       ylab="",
       las=1,
       yaxt='n')
  lines(2013:2018,vectANS_EPA[2:7]*100,col="#FFFF00",type="l",lwd=5)
  lines(2013:2018,vectANS_ANS[2:7]*100,col="#009900",type="l",lwd=5)
  lines(2013:2018,vectANS_END[2:7]*100,col="#FF6699",type="l",lwd=5)
  #lines(2013:2018,vectANS_LIEC_X*100,col="grey",type="l",lwd=5)
  axis(side = 2, at=c(0,50,100),las=1,labels=c("0%","50%","100%"))
  axis(side = 1, at=2013:2018, pos=-15, tick=F, labels=c(paste("n=",n2013,sep=""),
                                                         paste("n=",n2014,sep=""),paste("n=",n2015,sep=""),
                                                         paste("n=",n2016,sep=""),paste("n=",n2017,sep=""),
                                                         paste("n=",n2018,sep="")))
  legend("topleft",
         legend=c("EPI","EPA","ANS","END"),
         fill=c("#000099","#FFFF00","#009900","#FF6699"))
  
  
  # erreures END
  plot(2013:2018,vectEND_EPI[2:7]*100,col="#000099",type="l",lwd=5,ylim=c(0,131),
       main="Appartenance r?elle des Endog?s identifi?s sur le terrain FC",
       xlab="",
       ylab="",
       las=1,
       yaxt='n')
  lines(2013:2018,vectEND_EPA[2:7]*100,col="#FFFF00",type="l",lwd=5)
  lines(2013:2018,vectEND_ANS[2:7]*100,col="#009900",type="l",lwd=5)
  lines(2013:2018,vectEND_END[2:7]*100,col="#FF6699",type="l",lwd=5)
  #lines(2013:2018,vectEND_LIEC_X*100,col="grey",type="l",lwd=5)
  axis(side = 2, at=c(0,50,100),las=1,labels=c("0%","50%","100%"))
  axis(side = 1, at=2013:2018, pos=-15, tick=F, labels=c(paste("n=",n2013,sep=""),
                                                         paste("n=",n2014,sep=""),paste("n=",n2015,sep=""),
                                                         paste("n=",n2016,sep=""),paste("n=",n2017,sep=""),
                                                         paste("n=",n2018,sep="")))
  legend("topleft",
         legend=c("EPI","EPA","ANS","END"),
         fill=c("#000099","#FFFF00","#009900","#FF6699"))
  
}


### Erreurs par LIEC ###
{
  
  # erreures EPI
  plot(2013:2018,vectEPI_EPI[2:7]*100,col="#000099",type="l",lwd=5,ylim=c(0,141),
       main="Erreur d'identification terrain des Epig?s FC",
       xlab="",
       ylab="",
       yaxt='n')
  lines(2013:2018,vectEPA_EPI[2:7]*100,col="#FFFF00",type="l",lwd=5)
  lines(2013:2018,vectANS_EPI[2:7]*100,col="#009900",type="l",lwd=5)
  lines(2013:2018,vectEND_EPI[2:7]*100,col="#FF6699",type="l",lwd=5)
  lines(2013:2018,vectLIEC_X_EPI[2:7]*100,col="grey",type="l",lwd=5)
  axis(side = 2, at=c(0,50,100),las=1,labels=c("0%","50%","100%"))
  axis(side = 1, at=2013:2018, pos=-15, tick=F, labels=c(paste("n=",n2013,sep=""),
                                                         paste("n=",n2014,sep=""),paste("n=",n2015,sep=""),
                                                         paste("n=",n2016,sep=""),paste("n=",n2017,sep=""),
                                                         paste("n=",n2018,sep="")))
  legend("topleft",
         legend=c("EPI","EPA","ANS","END","LIEC_X"),
         fill=c("#000099","#FFFF00","#009900","#FF6699","grey"))
  
  
  
  # erreures EPA
  plot(2013:2018,vectEPI_EPA[2:7]*100,col="#000099",type="l",lwd=5,ylim=c(0,141),
       main="Erreur d'identification terrain des Epi-An?ciques FC",
       xlab="",
       ylab="",
       las=1,
       yaxt='n')
  lines(2013:2018,vectEPA_EPA[2:7]*100,col="#FFFF00",type="l",lwd=5)
  lines(2013:2018,vectANS_EPA[2:7]*100,col="#009900",type="l",lwd=5)
  lines(2013:2018,vectEND_EPA[2:7]*100,col="#FF6699",type="l",lwd=5)
  lines(2013:2018,vectLIEC_X_EPA[2:7]*100,col="grey",type="l",lwd=5)
  axis(side = 2, at=c(0,50,100),las=1,labels=c("0%","50%","100%"))
  axis(side = 1, at=2013:2018, pos=-15, tick=F, labels=c(paste("n=",n2013,sep=""),
                                                         paste("n=",n2014,sep=""),paste("n=",n2015,sep=""),
                                                         paste("n=",n2016,sep=""),paste("n=",n2017,sep=""),
                                                         paste("n=",n2018,sep="")))
  legend("topleft",
         legend=c("EPI","EPA","ANS","END","LIEC_X"),
         fill=c("#000099","#FFFF00","#009900","#FF6699","grey"))
  
  
  # erreures ANS
  plot(2013:2018,vectEPI_ANS[2:7]*100,col="#000099",type="l",lwd=5,ylim=c(0,141),
       main="Erreur d'identification terrain des An?ciques-Stricts FC",
       xlab="",
       ylab="",
       las=1,
       yaxt='n')
  lines(2013:2018,vectEPA_ANS[2:7]*100,col="#FFFF00",type="l",lwd=5)
  lines(2013:2018,vectANS_ANS[2:7]*100,col="#009900",type="l",lwd=5)
  lines(2013:2018,vectEND_ANS[2:7]*100,col="#FF6699",type="l",lwd=5)
  lines(2013:2018,vectLIEC_X_ANS[2:7]*100,col="grey",type="l",lwd=5)
  axis(side = 2, at=c(0,50,100),las=1,labels=c("0%","50%","100%"))
  axis(side = 1, at=2013:2018, pos=-15, tick=F, labels=c(paste("n=",n2013,sep=""),
                                                         paste("n=",n2014,sep=""),paste("n=",n2015,sep=""),
                                                         paste("n=",n2016,sep=""),paste("n=",n2017,sep=""),
                                                         paste("n=",n2018,sep="")))
  legend("topleft",
         legend=c("EPI","EPA","ANS","END","LIEC_X"),
         fill=c("#000099","#FFFF00","#009900","#FF6699","grey"))
  
  
  # erreures END
  plot(2013:2018,vectEPI_END[2:7]*100,col="#000099",type="l",lwd=5,ylim=c(0,141),
       main="Erreur d'identification terrain des Endog?s FC",
       xlab="",
       ylab="",
       las=1,
       yaxt='n')
  lines(2013:2018,vectEPA_END[2:7]*100,col="#FFFF00",type="l",lwd=5)
  lines(2013:2018,vectANS_END[2:7]*100,col="#009900",type="l",lwd=5)
  lines(2013:2018,vectEND_END[2:7]*100,col="#FF6699",type="l",lwd=5)
  lines(2013:2018,vectLIEC_X_END[2:7]*100,col="grey",type="l",lwd=5)
  axis(side = 2, at=c(0,50,100),las=1,labels=c("0%","50%","100%"))
  axis(side = 1, at=2013:2018, pos=-15, tick=F, labels=c(paste("n=",n2013,sep=""),
                                                         paste("n=",n2014,sep=""),paste("n=",n2015,sep=""),
                                                         paste("n=",n2016,sep=""),paste("n=",n2017,sep=""),
                                                         paste("n=",n2018,sep="")))
  legend("topleft",
         legend=c("EPI","EPA","ANS","END","LIEC_X"),
         fill=c("#000099","#FFFF00","#009900","#FF6699","grey"))
  
}








##### Region bien formee : PC #####

tabPC=Tab2[Tab2$Region%in%"PC",]
nrow(tabPC)==4179
dataID <- unique(sbt_m[,c("cle",colID,"Region")])
nrow(dataID)==1109
dataID=dataID[dataID$Region%in%"PC",]
nrow(dataID)==189

n2013=table(dataID$Annee)[2]
n2014=table(dataID$Annee)[3]
n2015=table(dataID$Annee)[4]
n2016=table(dataID$Annee)[5]
n2017=table(dataID$Annee)[6]
n2018=table(dataID$Annee)[7]

vectReussitePC=c(
sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC%in%c("EPI_EPI","EPA_EPA","ANS_ANS","END_END") & tabPC$Annee%in%2013])/sum(tabPC$Nbr_EW[tabPC$Annee%in%2013]),
sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC%in%c("EPI_EPI","EPA_EPA","ANS_ANS","END_END") & tabPC$Annee%in%2014])/sum(tabPC$Nbr_EW[tabPC$Annee%in%2014]),
sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC%in%c("EPI_EPI","EPA_EPA","ANS_ANS","END_END") & tabPC$Annee%in%2015])/sum(tabPC$Nbr_EW[tabPC$Annee%in%2015]),
sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC%in%c("EPI_EPI","EPA_EPA","ANS_ANS","END_END") & tabPC$Annee%in%2016])/sum(tabPC$Nbr_EW[tabPC$Annee%in%2016]),
sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC%in%c("EPI_EPI","EPA_EPA","ANS_ANS","END_END") & tabPC$Annee%in%2017])/sum(tabPC$Nbr_EW[tabPC$Annee%in%2017]),
sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC%in%c("EPI_EPI","EPA_EPA","ANS_ANS","END_END") & tabPC$Annee%in%2018])/sum(tabPC$Nbr_EW[tabPC$Annee%in%2018]))


##### taux pour les _EPI #####


# calcule des sommes

som_EPItt=sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPI_EPI"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPA_EPI"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="ANS_EPI"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="END_EPI"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="LIEC_X_EPI"])

som_EPI12=sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPI_EPI" & tabPC$Annee=="2012"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPA_EPI" & tabPC$Annee=="2012"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="ANS_EPI" & tabPC$Annee=="2012"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="END_EPI" & tabPC$Annee=="2012"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="LIEC_X_EPI" & tabPC$Annee=="2012"])

som_EPI13=sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPI_EPI" & tabPC$Annee=="2013"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPA_EPI" & tabPC$Annee=="2013"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="ANS_EPI" & tabPC$Annee=="2013"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="END_EPI" & tabPC$Annee=="2013"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="LIEC_X_EPI" & tabPC$Annee=="2013"])

som_EPI14=sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPI_EPI" & tabPC$Annee=="2014"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPA_EPI" & tabPC$Annee=="2014"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="ANS_EPI" & tabPC$Annee=="2014"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="END_EPI" & tabPC$Annee=="2014"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="LIEC_X_EPI" & tabPC$Annee=="2014"])

som_EPI15=sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPI_EPI" & tabPC$Annee=="2015"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPA_EPI" & tabPC$Annee=="2015"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="ANS_EPI" & tabPC$Annee=="2015"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="END_EPI" & tabPC$Annee=="2015"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="LIEC_X_EPI" & tabPC$Annee=="2015"])

som_EPI16=sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPI_EPI" & tabPC$Annee=="2016"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPA_EPI" & tabPC$Annee=="2016"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="ANS_EPI" & tabPC$Annee=="2016"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="END_EPI" & tabPC$Annee=="2016"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="LIEC_X_EPI" & tabPC$Annee=="2016"])

som_EPI17=sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPI_EPI" & tabPC$Annee=="2017"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPA_EPI" & tabPC$Annee=="2017"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="ANS_EPI" & tabPC$Annee=="2017"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="END_EPI" & tabPC$Annee=="2017"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="LIEC_X_EPI" & tabPC$Annee=="2017"])

som_EPI18=sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPI_EPI" & tabPC$Annee=="2018"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPA_EPI" & tabPC$Annee=="2018"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="ANS_EPI" & tabPC$Annee=="2018"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="END_EPI" & tabPC$Annee=="2018"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="LIEC_X_EPI" & tabPC$Annee=="2018"])




# vecteurs pour chaque annee

vectEPI_EPI=c(
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPI_EPI" & tabPC$Annee=="2012"])/
    (som_EPI12),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPI_EPI" & tabPC$Annee=="2013"])/
    (som_EPI13),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPI_EPI" & tabPC$Annee=="2014"])/
    (som_EPI14),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPI_EPI" & tabPC$Annee=="2015"])/
    (som_EPI15),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPI_EPI" & tabPC$Annee=="2016"])/
    (som_EPI16),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPI_EPI" & tabPC$Annee=="2017"])/
    (som_EPI17),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPI_EPI" & tabPC$Annee=="2018"])/
    (som_EPI18))

vectEPA_EPI=c(
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPA_EPI" & tabPC$Annee=="2012"])/
    (som_EPI12),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPA_EPI" & tabPC$Annee=="2013"])/
    (som_EPI13),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPA_EPI" & tabPC$Annee=="2014"])/
    (som_EPI14),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPA_EPI" & tabPC$Annee=="2015"])/
    (som_EPI15),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPA_EPI" & tabPC$Annee=="2016"])/
    (som_EPI16),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPA_EPI" & tabPC$Annee=="2017"])/
    (som_EPI17),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPA_EPI" & tabPC$Annee=="2018"])/
    (som_EPI18))

vectANS_EPI=c(
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="ANS_EPI" & tabPC$Annee=="2012"])/
    (som_EPI12),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="ANS_EPI" & tabPC$Annee=="2013"])/
    (som_EPI13),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="ANS_EPI" & tabPC$Annee=="2014"])/
    (som_EPI14),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="ANS_EPI" & tabPC$Annee=="2015"])/
    (som_EPI15),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="ANS_EPI" & tabPC$Annee=="2016"])/
    (som_EPI16),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="ANS_EPI" & tabPC$Annee=="2017"])/
    (som_EPI17),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="ANS_EPI" & tabPC$Annee=="2018"])/
    (som_EPI18))

vectEND_EPI=c(
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="END_EPI" & tabPC$Annee=="2012"])/
    (som_EPI12),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="END_EPI" & tabPC$Annee=="2013"])/
    (som_EPI13),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="END_EPI" & tabPC$Annee=="2014"])/
    (som_EPI14),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="END_EPI" & tabPC$Annee=="2015"])/
    (som_EPI15),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="END_EPI" & tabPC$Annee=="2016"])/
    (som_EPI16),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="END_EPI" & tabPC$Annee=="2017"])/
    (som_EPI17),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="END_EPI" & tabPC$Annee=="2018"])/
    (som_EPI18))

vectLIEC_X_EPI=c(
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="LIEC_X_EPI" & tabPC$Annee=="2012"])/
    (som_EPI12),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="LIEC_X_EPI" & tabPC$Annee=="2013"])/
    (som_EPI13),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="LIEC_X_EPI" & tabPC$Annee=="2014"])/
    (som_EPI14),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="LIEC_X_EPI" & tabPC$Annee=="2015"])/
    (som_EPI15),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="LIEC_X_EPI" & tabPC$Annee=="2016"])/
    (som_EPI16),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="LIEC_X_EPI" & tabPC$Annee=="2017"])/
    (som_EPI17),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="LIEC_X_EPI" & tabPC$Annee=="2018"])/
    (som_EPI18))




##### taux pour les _EPA #####


# calcule des sommes

som_EPAtt=sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPI_EPA"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPA_EPA"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="ANS_EPA"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="END_EPA"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="LIEC_X_EPA"])

som_EPA12=sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPI_EPA" & tabPC$Annee=="2012"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPA_EPA" & tabPC$Annee=="2012"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="ANS_EPA" & tabPC$Annee=="2012"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="END_EPA" & tabPC$Annee=="2012"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="LIEC_X_EPA" & tabPC$Annee=="2012"])

som_EPA13=sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPI_EPA" & tabPC$Annee=="2013"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPA_EPA" & tabPC$Annee=="2013"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="ANS_EPA" & tabPC$Annee=="2013"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="END_EPA" & tabPC$Annee=="2013"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="LIEC_X_EPA" & tabPC$Annee=="2013"])

som_EPA14=sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPI_EPA" & tabPC$Annee=="2014"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPA_EPA" & tabPC$Annee=="2014"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="ANS_EPA" & tabPC$Annee=="2014"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="END_EPA" & tabPC$Annee=="2014"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="LIEC_X_EPA" & tabPC$Annee=="2014"])

som_EPA15=sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPI_EPA" & tabPC$Annee=="2015"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPA_EPA" & tabPC$Annee=="2015"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="ANS_EPA" & tabPC$Annee=="2015"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="END_EPA" & tabPC$Annee=="2015"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="LIEC_X_EPA" & tabPC$Annee=="2015"])

som_EPA16=sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPI_EPA" & tabPC$Annee=="2016"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPA_EPA" & tabPC$Annee=="2016"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="ANS_EPA" & tabPC$Annee=="2016"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="END_EPA" & tabPC$Annee=="2016"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="LIEC_X_EPA" & tabPC$Annee=="2016"])

som_EPA17=sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPI_EPA" & tabPC$Annee=="2017"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPA_EPA" & tabPC$Annee=="2017"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="ANS_EPA" & tabPC$Annee=="2017"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="END_EPA" & tabPC$Annee=="2017"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="LIEC_X_EPA" & tabPC$Annee=="2017"])

som_EPA18=sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPI_EPA" & tabPC$Annee=="2018"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPA_EPA" & tabPC$Annee=="2018"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="ANS_EPA" & tabPC$Annee=="2018"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="END_EPA" & tabPC$Annee=="2018"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="LIEC_X_EPA" & tabPC$Annee=="2018"])






# vecteurs pour chaque annee

vectEPI_EPA=c(
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPI_EPA" & tabPC$Annee=="2012"])/
    (som_EPA12),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPI_EPA" & tabPC$Annee=="2013"])/
    (som_EPA13),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPI_EPA" & tabPC$Annee=="2014"])/
    (som_EPA14),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPI_EPA" & tabPC$Annee=="2015"])/
    (som_EPA15),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPI_EPA" & tabPC$Annee=="2016"])/
    (som_EPA16),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPI_EPA" & tabPC$Annee=="2017"])/
    (som_EPA17),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPI_EPA" & tabPC$Annee=="2018"])/
    (som_EPA18))

vectEPA_EPA=c(
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPA_EPA" & tabPC$Annee=="2012"])/
    (som_EPA12),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPA_EPA" & tabPC$Annee=="2013"])/
    (som_EPA13),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPA_EPA" & tabPC$Annee=="2014"])/
    (som_EPA14),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPA_EPA" & tabPC$Annee=="2015"])/
    (som_EPA15),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPA_EPA" & tabPC$Annee=="2016"])/
    (som_EPA16),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPA_EPA" & tabPC$Annee=="2017"])/
    (som_EPA17),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPA_EPA" & tabPC$Annee=="2018"])/
    (som_EPA18))

vectANS_EPA=c(
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="ANS_EPA" & tabPC$Annee=="2012"])/
    (som_EPA12),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="ANS_EPA" & tabPC$Annee=="2013"])/
    (som_EPA13),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="ANS_EPA" & tabPC$Annee=="2014"])/
    (som_EPA14),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="ANS_EPA" & tabPC$Annee=="2015"])/
    (som_EPA15),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="ANS_EPA" & tabPC$Annee=="2016"])/
    (som_EPA16),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="ANS_EPA" & tabPC$Annee=="2017"])/
    (som_EPA17),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="ANS_EPA" & tabPC$Annee=="2018"])/
    (som_EPA18))

vectEND_EPA=c(
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="END_EPA" & tabPC$Annee=="2012"])/
    (som_EPA12),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="END_EPA" & tabPC$Annee=="2013"])/
    (som_EPA13),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="END_EPA" & tabPC$Annee=="2014"])/
    (som_EPA14),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="END_EPA" & tabPC$Annee=="2015"])/
    (som_EPA15),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="END_EPA" & tabPC$Annee=="2016"])/
    (som_EPA16),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="END_EPA" & tabPC$Annee=="2017"])/
    (som_EPA17),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="END_EPA" & tabPC$Annee=="2018"])/
    (som_EPA18))

vectLIEC_X_EPA=c(
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="LIEC_X_EPA" & tabPC$Annee=="2012"])/
    (som_EPA12),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="LIEC_X_EPA" & tabPC$Annee=="2013"])/
    (som_EPA13),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="LIEC_X_EPA" & tabPC$Annee=="2014"])/
    (som_EPA14),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="LIEC_X_EPA" & tabPC$Annee=="2015"])/
    (som_EPA15),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="LIEC_X_EPA" & tabPC$Annee=="2016"])/
    (som_EPA16),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="LIEC_X_EPA" & tabPC$Annee=="2017"])/
    (som_EPA17),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="LIEC_X_EPA" & tabPC$Annee=="2018"])/
    (som_EPA18))




##### taux pour les _ANS #####


# calcule des sommes

som_ANStt=sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPI_ANS"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPA_ANS"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="ANS_ANS"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="END_ANS"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="LIEC_X_ANS"])

som_ANS12=sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPI_ANS" & tabPC$Annee=="2012"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPA_ANS" & tabPC$Annee=="2012"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="ANS_ANS" & tabPC$Annee=="2012"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="END_ANS" & tabPC$Annee=="2012"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="LIEC_X_ANS" & tabPC$Annee=="2012"])

som_ANS13=sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPI_ANS" & tabPC$Annee=="2013"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPA_ANS" & tabPC$Annee=="2013"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="ANS_ANS" & tabPC$Annee=="2013"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="END_ANS" & tabPC$Annee=="2013"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="LIEC_X_ANS" & tabPC$Annee=="2013"])

som_ANS14=sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPI_ANS" & tabPC$Annee=="2014"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPA_ANS" & tabPC$Annee=="2014"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="ANS_ANS" & tabPC$Annee=="2014"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="END_ANS" & tabPC$Annee=="2014"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="LIEC_X_ANS" & tabPC$Annee=="2014"])

som_ANS15=sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPI_ANS" & tabPC$Annee=="2015"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPA_ANS" & tabPC$Annee=="2015"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="ANS_ANS" & tabPC$Annee=="2015"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="END_ANS" & tabPC$Annee=="2015"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="LIEC_X_ANS" & tabPC$Annee=="2015"])

som_ANS16=sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPI_ANS" & tabPC$Annee=="2016"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPA_ANS" & tabPC$Annee=="2016"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="ANS_ANS" & tabPC$Annee=="2016"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="END_ANS" & tabPC$Annee=="2016"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="LIEC_X_ANS" & tabPC$Annee=="2016"])

som_ANS17=sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPI_ANS" & tabPC$Annee=="2017"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPA_ANS" & tabPC$Annee=="2017"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="ANS_ANS" & tabPC$Annee=="2017"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="END_ANS" & tabPC$Annee=="2017"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="LIEC_X_ANS" & tabPC$Annee=="2017"])

som_ANS18=sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPI_ANS" & tabPC$Annee=="2018"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPA_ANS" & tabPC$Annee=="2018"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="ANS_ANS" & tabPC$Annee=="2018"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="END_ANS" & tabPC$Annee=="2018"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="LIEC_X_ANS" & tabPC$Annee=="2018"])






# vecteurs pour chaque annee

vectEPI_ANS=c(
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPI_ANS" & tabPC$Annee=="2012"])/
    (som_ANS12),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPI_ANS" & tabPC$Annee=="2013"])/
    (som_ANS13),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPI_ANS" & tabPC$Annee=="2014"])/
    (som_ANS14),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPI_ANS" & tabPC$Annee=="2015"])/
    (som_ANS15),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPI_ANS" & tabPC$Annee=="2016"])/
    (som_ANS16),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPI_ANS" & tabPC$Annee=="2017"])/
    (som_ANS17),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPI_ANS" & tabPC$Annee=="2018"])/
    (som_ANS18))

vectEPA_ANS=c(
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPA_ANS" & tabPC$Annee=="2012"])/
    (som_ANS12),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPA_ANS" & tabPC$Annee=="2013"])/
    (som_ANS13),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPA_ANS" & tabPC$Annee=="2014"])/
    (som_ANS14),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPA_ANS" & tabPC$Annee=="2015"])/
    (som_ANS15),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPA_ANS" & tabPC$Annee=="2016"])/
    (som_ANS16),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPA_ANS" & tabPC$Annee=="2017"])/
    (som_ANS17),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPA_ANS" & tabPC$Annee=="2018"])/
    (som_ANS18))

vectANS_ANS=c(
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="ANS_ANS" & tabPC$Annee=="2012"])/
    (som_ANS12),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="ANS_ANS" & tabPC$Annee=="2013"])/
    (som_ANS13),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="ANS_ANS" & tabPC$Annee=="2014"])/
    (som_ANS14),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="ANS_ANS" & tabPC$Annee=="2015"])/
    (som_ANS15),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="ANS_ANS" & tabPC$Annee=="2016"])/
    (som_ANS16),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="ANS_ANS" & tabPC$Annee=="2017"])/
    (som_ANS17),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="ANS_ANS" & tabPC$Annee=="2018"])/
    (som_ANS18))

vectEND_ANS=c(
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="END_ANS" & tabPC$Annee=="2012"])/
    (som_ANS12),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="END_ANS" & tabPC$Annee=="2013"])/
    (som_ANS13),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="END_ANS" & tabPC$Annee=="2014"])/
    (som_ANS14),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="END_ANS" & tabPC$Annee=="2015"])/
    (som_ANS15),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="END_ANS" & tabPC$Annee=="2016"])/
    (som_ANS16),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="END_ANS" & tabPC$Annee=="2017"])/
    (som_ANS17),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="END_ANS" & tabPC$Annee=="2018"])/
    (som_ANS18))

vectLIEC_X_ANS=c(
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="LIEC_X_ANS" & tabPC$Annee=="2012"])/
    (som_ANS12),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="LIEC_X_ANS" & tabPC$Annee=="2013"])/
    (som_ANS13),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="LIEC_X_ANS" & tabPC$Annee=="2014"])/
    (som_ANS14),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="LIEC_X_ANS" & tabPC$Annee=="2015"])/
    (som_ANS15),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="LIEC_X_ANS" & tabPC$Annee=="2016"])/
    (som_ANS16),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="LIEC_X_ANS" & tabPC$Annee=="2017"])/
    (som_ANS17),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="LIEC_X_ANS" & tabPC$Annee=="2018"])/
    (som_ANS18))




##### taux pour les _END #####


# calcule des sommes

som_ENDtt=sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPI_END"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPA_END"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="ANS_END"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="END_END"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="LIEC_X_END"])

som_END12=sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPI_END" & tabPC$Annee=="2012"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPA_END" & tabPC$Annee=="2012"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="ANS_END" & tabPC$Annee=="2012"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="END_END" & tabPC$Annee=="2012"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="LIEC_X_END" & tabPC$Annee=="2012"])

som_END13=sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPI_END" & tabPC$Annee=="2013"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPA_END" & tabPC$Annee=="2013"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="ANS_END" & tabPC$Annee=="2013"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="END_END" & tabPC$Annee=="2013"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="LIEC_X_END" & tabPC$Annee=="2013"])

som_END14=sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPI_END" & tabPC$Annee=="2014"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPA_END" & tabPC$Annee=="2014"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="ANS_END" & tabPC$Annee=="2014"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="END_END" & tabPC$Annee=="2014"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="LIEC_X_END" & tabPC$Annee=="2014"])

som_END15=sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPI_END" & tabPC$Annee=="2015"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPA_END" & tabPC$Annee=="2015"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="ANS_END" & tabPC$Annee=="2015"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="END_END" & tabPC$Annee=="2015"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="LIEC_X_END" & tabPC$Annee=="2015"])

som_END16=sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPI_END" & tabPC$Annee=="2016"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPA_END" & tabPC$Annee=="2016"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="ANS_END" & tabPC$Annee=="2016"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="END_END" & tabPC$Annee=="2016"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="LIEC_X_END" & tabPC$Annee=="2016"])

som_END17=sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPI_END" & tabPC$Annee=="2017"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPA_END" & tabPC$Annee=="2017"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="ANS_END" & tabPC$Annee=="2017"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="END_END" & tabPC$Annee=="2017"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="LIEC_X_END" & tabPC$Annee=="2017"])

som_END18=sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPI_END" & tabPC$Annee=="2018"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPA_END" & tabPC$Annee=="2018"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="ANS_END" & tabPC$Annee=="2018"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="END_END" & tabPC$Annee=="2018"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="LIEC_X_END" & tabPC$Annee=="2018"])






# vecteurs pour chaque annee

vectEPI_END=c(
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPI_END" & tabPC$Annee=="2012"])/
    (som_END12),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPI_END" & tabPC$Annee=="2013"])/
    (som_END13),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPI_END" & tabPC$Annee=="2014"])/
    (som_END14),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPI_END" & tabPC$Annee=="2015"])/
    (som_END15),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPI_END" & tabPC$Annee=="2016"])/
    (som_END16),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPI_END" & tabPC$Annee=="2017"])/
    (som_END17),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPI_END" & tabPC$Annee=="2018"])/
    (som_END18))

vectEPA_END=c(
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPA_END" & tabPC$Annee=="2012"])/
    (som_END12),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPA_END" & tabPC$Annee=="2013"])/
    (som_END13),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPA_END" & tabPC$Annee=="2014"])/
    (som_END14),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPA_END" & tabPC$Annee=="2015"])/
    (som_END15),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPA_END" & tabPC$Annee=="2016"])/
    (som_END16),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPA_END" & tabPC$Annee=="2017"])/
    (som_END17),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPA_END" & tabPC$Annee=="2018"])/
    (som_END18))

vectANS_END=c(
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="ANS_END" & tabPC$Annee=="2012"])/
    (som_END12),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="ANS_END" & tabPC$Annee=="2013"])/
    (som_END13),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="ANS_END" & tabPC$Annee=="2014"])/
    (som_END14),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="ANS_END" & tabPC$Annee=="2015"])/
    (som_END15),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="ANS_END" & tabPC$Annee=="2016"])/
    (som_END16),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="ANS_END" & tabPC$Annee=="2017"])/
    (som_END17),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="ANS_END" & tabPC$Annee=="2018"])/
    (som_END18))

vectEND_END=c(
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="END_END" & tabPC$Annee=="2012"])/
    (som_END12),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="END_END" & tabPC$Annee=="2013"])/
    (som_END13),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="END_END" & tabPC$Annee=="2014"])/
    (som_END14),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="END_END" & tabPC$Annee=="2015"])/
    (som_END15),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="END_END" & tabPC$Annee=="2016"])/
    (som_END16),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="END_END" & tabPC$Annee=="2017"])/
    (som_END17),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="END_END" & tabPC$Annee=="2018"])/
    (som_END18))

vectLIEC_X_END=c(
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="LIEC_X_END" & tabPC$Annee=="2012"])/
    (som_END12),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="LIEC_X_END" & tabPC$Annee=="2013"])/
    (som_END13),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="LIEC_X_END" & tabPC$Annee=="2014"])/
    (som_END14),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="LIEC_X_END" & tabPC$Annee=="2015"])/
    (som_END15),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="LIEC_X_END" & tabPC$Annee=="2016"])/
    (som_END16),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="LIEC_X_END" & tabPC$Annee=="2017"])/
    (som_END17),
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="LIEC_X_END" & tabPC$Annee=="2018"])/
    (som_END18))



##### Graphiques #####





plot(2013:2018,vectEPI_EPI[2:7]*100,col="#000099",type="l",lwd=4,ylim=c(0,131),
     main="Taux d'identification terrain r?ussi PC",
     xlab="",
     ylab="",
     yaxt='n')
lines(2013:2018,vectEPA_EPA[2:7]*100,col="#FFFF00",type="l",lwd=4)
lines(2013:2018,vectANS_ANS[2:7]*100,col="#009900",type="l",lwd=4)
lines(2013:2018,vectEND_END[2:7]*100,col="#FF6699",type="l",lwd=4)
lines(2013:2018,vectReussitePC*100,col="black",type="l",lwd=5)
axis(side = 2, at=c(0,50,100),las=1,labels=c("0%","50%","100%"))
axis(side = 1, at=2013:2018, pos=-15, tick=F, labels=c(paste("n=",n2013,sep=""),
                                                       paste("n=",n2014,sep=""),paste("n=",n2015,sep=""),
                                                       paste("n=",n2016,sep=""),paste("n=",n2017,sep=""),
                                                       paste("n=",n2018,sep="")))
legend("topleft",
       legend=c("EPI","EPA","ANS","END"),
       fill=c("#000099","#FFFF00","#009900","#FF6699"))


yEPI=sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPI_EPI"])/som_EPItt
yEPA=sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPA_EPA"])/som_EPAtt
yANS=sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="ANS_ANS"])/som_ANStt
yEND=sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="END_END"])/som_ENDtt


barplot(c(yEPI,yEPA,yANS,yEND),ylim=c(0,1))


tab_bp=data.frame(dt=c(yEPI,yEPA,yANS,yEND),
                  x=rep("FC",4),
                  LIEC=c("EPI","EPA","ANS","END"),
                  pos_text=c(1,1,1,1),
                  Pos=cumsum(c(yEPI,yEPA,yANS,yEND)))

  som_EPItt=sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPI_EPI"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="EPA_EPI"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="ANS_EPI"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="END_EPI"])+
  sum(tabPC$Nbr_EW[tabPC$FIEC_LIEC=="LIEC_X_EPI"])

ggplot(data=tab_bp, aes(x=x,y=dt)) +
       theme_bw(base_size=12) + # Mettre font blanc
       theme(axis.text=element_text(face="bold", size=50),
             axis.title.x = element_text(face="bold", size=50, vjust=0),
             axis.ticks.margin=unit(c(0.5),'cm'),
             axis.ticks=element_line(size=3),
             axis.ticks.length=unit(0.5, "cm"),
             axis.text.x=element_blank(),
             axis.ticks.x=element_blank(),
             panel.grid.major=element_line(size=1.5), # Augmenter taille ligne y
             panel.grid.major.x=element_line(colour=NA), # Suppression lignes axes x
             axis.title=element_text( vjust=5, size=15),
             plot.margin=unit(c(0.5, 0.5, 0.5, 0.5), "cm"), # Augmenter marges
             line=element_line(size=1.5),
             legend.position="none" ) +
       geom_bar(aes(fill=LIEC), stat="identity", width=0.8) + # Barplot ,
       geom_point(aes(y=replace(Pos,which(Nbr_EW==0),NA), x=pos_text ) , 
                  color="grey",  size=63 ) +
       geom_point(aes(y=replace(Pos,which(Nbr_EW==0),NA),  x=pos_text,
                      shape=LIEC, fill=LIEC ), color="grey",  size=60 ) +
       geom_text(aes(y=Pos, colour=LIEC, x=pos_text, label=replace(round(Nbr_EW,1), which(Nbr_EW==0),NA) ) , # Affichage valeur Ab >0
                 size=17, fontface=2) +
       scale_fill_manual(values=rev(c("#7B7B7B", "#FF6699", "#80804D" ,"#009900", "#FFFF00", "#000099")), # Couleur LIEC (Fill)
                         labels=rev(c("Groupe ind?terminable" ,"Endog?s", "NCX" ,"A.T?te Noire", "A.T?te Rouge", "Epig?s"))) + 
       scale_colour_manual(values=rev(c("white", "black", "black", "white", "black", "white")) ) + # couleur texte
       scale_shape_manual(values=rep(21,6)) + # forme rond avec couleur de fond
       scale_y_continuous(expand=c(0.06,0)) +
       scale_x_continuous(NULL) +
       #scale_x_continuous(paste("Parcelle ", data_site$ID_Site, sep=""), limits=c(0.5,1.5) ) +
       ylab(NULL) # SUppression titre axes x et y








### Appartenances Reelles ###
{
  
  # erreures EPI
  plot(2013:2018,vectEPI_EPI[2:7]*100,col="#000099",type="l",lwd=5,ylim=c(0,131),
       main="Appartenance r?elle des Epig?s identifi?s sur le terrain PC",
       xlab="",
       ylab="",
       yaxt='n')
  lines(2013:2018,vectEPI_EPA[2:7]*100,col="#FFFF00",type="l",lwd=5)
  lines(2013:2018,vectEPI_ANS[2:7]*100,col="#009900",type="l",lwd=5)
  lines(2013:2018,vectEPI_END[2:7]*100,col="#FF6699",type="l",lwd=5)
  #lines(2013:2018,vectEPI_LIEC_X*100,col="grey",type="l",lwd=5)
  axis(side = 2, at=c(0,50,100),las=1,labels=c("0%","50%","100%"))
  axis(side = 1, at=2013:2018, pos=-15, tick=F, labels=c(paste("n=",n2013,sep=""),
                                                         paste("n=",n2014,sep=""),paste("n=",n2015,sep=""),
                                                         paste("n=",n2016,sep=""),paste("n=",n2017,sep=""),
                                                         paste("n=",n2018,sep="")))
  legend("topleft",
         legend=c("EPI","EPA","ANS","END"),
         fill=c("#000099","#FFFF00","#009900","#FF6699"))
  
  
  
  # erreures EPA
  plot(2013:2018,vectEPA_EPI[2:7]*100,col="#000099",type="l",lwd=5,ylim=c(0,131),
       main="Appartenance r?elle des Epi-An?ciques identifi?s sur le terrain PC",
       xlab="",
       ylab="",
       las=1,
       yaxt='n')
  lines(2013:2018,vectEPA_EPA[2:7]*100,col="#FFFF00",type="l",lwd=5)
  lines(2013:2018,vectEPA_ANS[2:7]*100,col="#009900",type="l",lwd=5)
  lines(2013:2018,vectEPA_END[2:7]*100,col="#FF6699",type="l",lwd=5)
  #lines(2013:2018,vectEPA_LIEC_X*100,col="grey",type="l",lwd=5)
  axis(side = 2, at=c(0,50,100),las=1,labels=c("0%","50%","100%"))
  axis(side = 1, at=2013:2018, pos=-15, tick=F, labels=c(paste("n=",n2013,sep=""),
                                                         paste("n=",n2014,sep=""),paste("n=",n2015,sep=""),
                                                         paste("n=",n2016,sep=""),paste("n=",n2017,sep=""),
                                                         paste("n=",n2018,sep="")))
  legend("topleft",
         legend=c("EPI","EPA","ANS","END"),
         fill=c("#000099","#FFFF00","#009900","#FF6699"))
  
  
  # erreures ANS
  plot(2013:2018,vectANS_EPI[2:7]*100,col="#000099",type="l",lwd=5,ylim=c(0,131),
       main="Appartenance r?elle des An?ciques-Stricts identifi?s sur le terrain PC",
       xlab="",
       ylab="",
       las=1,
       yaxt='n')
  lines(2013:2018,vectANS_EPA[2:7]*100,col="#FFFF00",type="l",lwd=5)
  lines(2013:2018,vectANS_ANS[2:7]*100,col="#009900",type="l",lwd=5)
  lines(2013:2018,vectANS_END[2:7]*100,col="#FF6699",type="l",lwd=5)
  #lines(2013:2018,vectANS_LIEC_X*100,col="grey",type="l",lwd=5)
  axis(side = 2, at=c(0,50,100),las=1,labels=c("0%","50%","100%"))
  axis(side = 1, at=2013:2018, pos=-15, tick=F, labels=c(paste("n=",n2013,sep=""),
                                                         paste("n=",n2014,sep=""),paste("n=",n2015,sep=""),
                                                         paste("n=",n2016,sep=""),paste("n=",n2017,sep=""),
                                                         paste("n=",n2018,sep="")))
  legend("topleft",
         legend=c("EPI","EPA","ANS","END"),
         fill=c("#000099","#FFFF00","#009900","#FF6699"))
  
  
  # erreures END
  plot(2013:2018,vectEND_EPI[2:7]*100,col="#000099",type="l",lwd=5,ylim=c(0,131),
       main="Appartenance r?elle des Endog?s identifi?s sur le terrain PC",
       xlab="",
       ylab="",
       las=1,
       yaxt='n')
  lines(2013:2018,vectEND_EPA[2:7]*100,col="#FFFF00",type="l",lwd=5)
  lines(2013:2018,vectEND_ANS[2:7]*100,col="#009900",type="l",lwd=5)
  lines(2013:2018,vectEND_END[2:7]*100,col="#FF6699",type="l",lwd=5)
  #lines(2013:2018,vectEND_LIEC_X*100,col="grey",type="l",lwd=5)
  axis(side = 2, at=c(0,50,100),las=1,labels=c("0%","50%","100%"))
  axis(side = 1, at=2013:2018, pos=-15, tick=F, labels=c(paste("n=",n2013,sep=""),
                                                         paste("n=",n2014,sep=""),paste("n=",n2015,sep=""),
                                                         paste("n=",n2016,sep=""),paste("n=",n2017,sep=""),
                                                         paste("n=",n2018,sep="")))
  legend("topleft",
         legend=c("EPI","EPA","ANS","END"),
         fill=c("#000099","#FFFF00","#009900","#FF6699"))
  
}


### Erreurs par LIEC ###
{
  
  # erreures EPI
  plot(2013:2018,vectEPI_EPI[2:7]*100,col="#000099",type="l",lwd=5,ylim=c(0,141),
       main="Erreur d'identification terrain des Epig?s PC",
       xlab="",
       ylab="",
       yaxt='n')
  lines(2013:2018,vectEPA_EPI[2:7]*100,col="#FFFF00",type="l",lwd=5)
  lines(2013:2018,vectANS_EPI[2:7]*100,col="#009900",type="l",lwd=5)
  lines(2013:2018,vectEND_EPI[2:7]*100,col="#FF6699",type="l",lwd=5)
  lines(2013:2018,vectLIEC_X_EPI[2:7]*100,col="grey",type="l",lwd=5)
  axis(side = 2, at=c(0,50,100),las=1,labels=c("0%","50%","100%"))
  axis(side = 1, at=2013:2018, pos=-15, tick=F, labels=c(paste("n=",n2013,sep=""),
                                                         paste("n=",n2014,sep=""),paste("n=",n2015,sep=""),
                                                         paste("n=",n2016,sep=""),paste("n=",n2017,sep=""),
                                                         paste("n=",n2018,sep="")))
  legend("topleft",
         legend=c("EPI","EPA","ANS","END","LIEC_X"),
         fill=c("#000099","#FFFF00","#009900","#FF6699","grey"))
  
  
  
  # erreures EPA
  plot(2013:2018,vectEPI_EPA[2:7]*100,col="#000099",type="l",lwd=5,ylim=c(0,141),
       main="Erreur d'identification terrain des Epi-An?ciques PC",
       xlab="",
       ylab="",
       las=1,
       yaxt='n')
  lines(2013:2018,vectEPA_EPA[2:7]*100,col="#FFFF00",type="l",lwd=5)
  lines(2013:2018,vectANS_EPA[2:7]*100,col="#009900",type="l",lwd=5)
  lines(2013:2018,vectEND_EPA[2:7]*100,col="#FF6699",type="l",lwd=5)
  lines(2013:2018,vectLIEC_X_EPA[2:7]*100,col="grey",type="l",lwd=5)
  axis(side = 2, at=c(0,50,100),las=1,labels=c("0%","50%","100%"))
  axis(side = 1, at=2013:2018, pos=-15, tick=F, labels=c(paste("n=",n2013,sep=""),
                                                         paste("n=",n2014,sep=""),paste("n=",n2015,sep=""),
                                                         paste("n=",n2016,sep=""),paste("n=",n2017,sep=""),
                                                         paste("n=",n2018,sep="")))
  legend("topleft",
         legend=c("EPI","EPA","ANS","END","LIEC_X"),
         fill=c("#000099","#FFFF00","#009900","#FF6699","grey"))
  
  
  # erreures ANS
  plot(2013:2018,vectEPI_ANS[2:7]*100,col="#000099",type="l",lwd=5,ylim=c(0,141),
       main="Erreur d'identification terrain des An?ciques-Stricts PC",
       xlab="",
       ylab="",
       las=1,
       yaxt='n')
  lines(2013:2018,vectEPA_ANS[2:7]*100,col="#FFFF00",type="l",lwd=5)
  lines(2013:2018,vectANS_ANS[2:7]*100,col="#009900",type="l",lwd=5)
  lines(2013:2018,vectEND_ANS[2:7]*100,col="#FF6699",type="l",lwd=5)
  lines(2013:2018,vectLIEC_X_ANS[2:7]*100,col="grey",type="l",lwd=5)
  axis(side = 2, at=c(0,50,100),las=1,labels=c("0%","50%","100%"))
  axis(side = 1, at=2013:2018, pos=-15, tick=F, labels=c(paste("n=",n2013,sep=""),
                                                         paste("n=",n2014,sep=""),paste("n=",n2015,sep=""),
                                                         paste("n=",n2016,sep=""),paste("n=",n2017,sep=""),
                                                         paste("n=",n2018,sep="")))
  legend("topleft",
         legend=c("EPI","EPA","ANS","END","LIEC_X"),
         fill=c("#000099","#FFFF00","#009900","#FF6699","grey"))
  
  
  # erreures END
  plot(2013:2018,vectEPI_END[2:7]*100,col="#000099",type="l",lwd=5,ylim=c(0,141),
       main="Erreur d'identification terrain des Endog?s PC",
       xlab="",
       ylab="",
       las=1,
       yaxt='n')
  lines(2013:2018,vectEPA_END[2:7]*100,col="#FFFF00",type="l",lwd=5)
  lines(2013:2018,vectANS_END[2:7]*100,col="#009900",type="l",lwd=5)
  lines(2013:2018,vectEND_END[2:7]*100,col="#FF6699",type="l",lwd=5)
  lines(2013:2018,vectLIEC_X_END[2:7]*100,col="grey",type="l",lwd=5)
  axis(side = 2, at=c(0,50,100),las=1,labels=c("0%","50%","100%"))
  axis(side = 1, at=2013:2018, pos=-15, tick=F, labels=c(paste("n=",n2013,sep=""),
                                                         paste("n=",n2014,sep=""),paste("n=",n2015,sep=""),
                                                         paste("n=",n2016,sep=""),paste("n=",n2017,sep=""),
                                                         paste("n=",n2018,sep="")))
  legend("topleft",
         legend=c("EPI","EPA","ANS","END","LIEC_X"),
         fill=c("#000099","#FFFF00","#009900","#FF6699","grey"))
  
}








##### Moyennes totale par LIEC #####


moyEPI_EPI=sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPI_EPI"])/(som_EPItt)
moyEPA_EPI=sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPA_EPI"])/(som_EPItt)
moyANS_EPI=sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="ANS_EPI"])/(som_EPItt)
moyEND_EPI=sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="END_EPI"])/(som_EPItt)
moyLIEC_X_EPI=sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="LIEC_X_EPI"])/(som_EPItt)
moyEPI_EPI+moyEPA_EPI+moyANS_EPI+moyEND_EPI+moyLIEC_X_EPI==1
moyEPI_EPI
moyEPA_EPI
moyANS_EPI
moyEND_EPI
moyLIEC_X_EPI


moyEPI_EPA=sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPI_EPA"])/(som_EPAtt)
moyEPA_EPA=sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPA_EPA"])/(som_EPAtt)
moyANS_EPA=sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="ANS_EPA"])/(som_EPAtt)
moyEND_EPA=sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="END_EPA"])/(som_EPAtt)
moyLIEC_X_EPA=sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="LIEC_X_EPA"])/(som_EPAtt)
moyEPI_EPA+moyEPA_EPA+moyANS_EPA+moyEND_EPA+moyLIEC_X_EPA==1
moyEPI_EPA
moyEPA_EPA
moyANS_EPA
moyEND_EPA
moyLIEC_X_EPA


moyEPI_ANS=sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPI_ANS"])/(som_ANStt)
moyEPA_ANS=sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPA_ANS"])/(som_ANStt)
moyANS_ANS=sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="ANS_ANS"])/(som_ANStt)
moyEND_ANS=sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="END_ANS"])/(som_ANStt)
moyLIEC_X_ANS=sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="LIEC_X_ANS"])/(som_ANStt)
moyEPI_ANS+moyEPA_ANS+moyANS_ANS+moyEND_ANS+moyLIEC_X_ANS==1
moyEPI_ANS
moyEPA_ANS
moyANS_ANS
moyEND_ANS
moyLIEC_X_ANS


moyEPI_END=sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPI_END"])/(som_ENDtt)
moyEPA_END=sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPA_END"])/(som_ENDtt)
moyANS_END=sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="ANS_END"])/(som_ENDtt)
moyEND_END=sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="END_END"])/(som_ENDtt)
moyLIEC_X_END=sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="LIEC_X_END"])/(som_ENDtt)
moyEPI_END+moyEPA_END+moyANS_END+moyEND_END+moyLIEC_X_END==1
moyEPI_END
moyEPA_END
moyANS_END
moyEND_END
moyLIEC_X_END





som_LIEC_Xtt=sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPI_LIEC_X"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPA_LIEC_X"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="ANS_LIEC_X"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="END_LIEC_X"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="LIEC_X_LIEC_X"])

som_LIEC_X12=sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPI_LIEC_X" & Tab2$Annee=="2012"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPA_LIEC_X" & Tab2$Annee=="2012"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="ANS_LIEC_X" & Tab2$Annee=="2012"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="END_LIEC_X" & Tab2$Annee=="2012"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="LIEC_X_LIEC_X" & Tab2$Annee=="2012"])

som_LIEC_X13=sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPI_LIEC_X" & Tab2$Annee=="2013"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPA_LIEC_X" & Tab2$Annee=="2013"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="ANS_LIEC_X" & Tab2$Annee=="2013"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="END_LIEC_X" & Tab2$Annee=="2013"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="LIEC_X_LIEC_X" & Tab2$Annee=="2013"])

som_LIEC_X14=sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPI_LIEC_X" & Tab2$Annee=="2014"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPA_LIEC_X" & Tab2$Annee=="2014"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="ANS_LIEC_X" & Tab2$Annee=="2014"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="END_LIEC_X" & Tab2$Annee=="2014"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="LIEC_X_LIEC_X" & Tab2$Annee=="2014"])

som_LIEC_X15=sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPI_LIEC_X" & Tab2$Annee=="2015"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPA_LIEC_X" & Tab2$Annee=="2015"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="ANS_LIEC_X" & Tab2$Annee=="2015"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="END_LIEC_X" & Tab2$Annee=="2015"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="LIEC_X_LIEC_X" & Tab2$Annee=="2015"])

som_LIEC_X16=sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPI_LIEC_X" & Tab2$Annee=="2016"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPA_LIEC_X" & Tab2$Annee=="2016"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="ANS_LIEC_X" & Tab2$Annee=="2016"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="END_LIEC_X" & Tab2$Annee=="2016"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="LIEC_X_LIEC_X" & Tab2$Annee=="2016"])

som_LIEC_X17=sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPI_LIEC_X" & Tab2$Annee=="2017"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPA_LIEC_X" & Tab2$Annee=="2017"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="ANS_LIEC_X" & Tab2$Annee=="2017"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="END_LIEC_X" & Tab2$Annee=="2017"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="LIEC_X_LIEC_X" & Tab2$Annee=="2017"])

som_LIEC_X18=sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPI_LIEC_X" & Tab2$Annee=="2018"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPA_LIEC_X" & Tab2$Annee=="2018"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="ANS_LIEC_X" & Tab2$Annee=="2018"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="END_LIEC_X" & Tab2$Annee=="2018"])+
  sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="LIEC_X_LIEC_X" & Tab2$Annee=="2018"])



moyEPI_LIEC_X=sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPI_LIEC_X"])/(som_LIEC_Xtt)
moyEPA_LIEC_X=sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="EPA_LIEC_X"])/(som_LIEC_Xtt)
moyANS_LIEC_X=sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="ANS_LIEC_X"])/(som_LIEC_Xtt)
moyEND_LIEC_X=sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="END_LIEC_X"])/(som_LIEC_Xtt)
moyLIEC_X_LIEC_X=sum(Tab2$Nbr_EW[Tab2$FIEC_LIEC=="LIEC_X_LIEC_X"])/(som_LIEC_Xtt)
moyEPI_LIEC_X+moyEPA_LIEC_X+moyANS_LIEC_X+moyEND_LIEC_X+moyLIEC_X_LIEC_X==1
moyEPI_LIEC_X
moyEPA_LIEC_X
moyANS_LIEC_X
moyEND_LIEC_X
moyLIEC_X_LIEC_X



save.image("D:/Home/naleveque/Documents/ECOBIO/Equipe/Abstract ISEE/envTauxErr2.Rdata")






##### Graphique 4 regions #####



plot(2013:2018,vectReussiteAQ*100,col=1,type="l",lwd=5,ylim=c(0,131),
     main="Taux d'identification terrain r?ussi",
     xlab="",
     ylab="",
     yaxt='n')
lines(2013:2018,vectReussiteBO*100,col=2,type="l",lwd=5)
lines(2013:2018,vectReussiteFC*100,col=3,type="l",lwd=5)
lines(2013:2018,vectReussitePC*100,col=4,type="l",lwd=5)
axis(side = 2, at=c(0,50,100),las=1,labels=c("0%","50%","100%"))
axis(side = 1, at=2013:2018, pos=-15, tick=F, labels=c(paste("n=",n2013,sep=""),
                                                       paste("n=",n2014,sep=""),paste("n=",n2015,sep=""),
                                                       paste("n=",n2016,sep=""),paste("n=",n2017,sep=""),
                                                       paste("n=",n2018,sep="")))
legend("topleft",
       legend=c("AQ","BO","FC","PC"),
       fill=c(1,2,3,4))


##### PC zoom sur les points critiques - en quoi les gens on classes les LIEC ? #####

#PC 2016 tr?s bons r?sultats : 0.9470069 taux de reussite

tabPC=Tab2[Tab2$Region%in%"PC",]
nrow(tabPC)==4179
dataID <- unique(sbt_m[,c("cle",colID,"Region")])
nrow(dataID)==1109
dataID=dataID[dataID$Region%in%"PC",]
nrow(dataID)==189

tabPC2016=tabPC[tabPC$Annee%in%2016,]
unique(tabPC2016$Stade)#"JV" "SA" "AD" => ok
any(is.na(tabPC2016$Stade))==F

#taux d'id correcte global
sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_EPI","EPA_EPA","ANS_ANS","END_END")])/
sum(tabPC2016$Nbr_EW)#0.9470069 => ok

#les EPI : en quoi ont ils ?t? class?s ?
{
  
  #taux id correcte EPI SA
  t.EPI_EPI.SA=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_EPI") & tabPC2016$Stade%in%"SA"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_EPI","EPA_EPI","ANS_EPI","END_EPI","LIEC_X_EPI")])#0.3811321
  
  #taux id correcte EPI AD
  t.EPI_EPI.AD=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_EPI") & tabPC2016$Stade%in%"AD"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_EPI","EPA_EPI","ANS_EPI","END_EPI","LIEC_X_EPI")])#0.3792453
  
  #taux id correcte EPI JV
  t.EPI_EPI.JV=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_EPI") & tabPC2016$Stade%in%"JV"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_EPI","EPA_EPI","ANS_EPI","END_EPI","LIEC_X_EPI")])#0.1396226
  
  
  
  
  #taux id incorrecte EPA SA
  t.EPA_EPI.SA=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPA_EPI") & tabPC2016$Stade%in%"SA"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_EPI","EPA_EPI","ANS_EPI","END_EPI","LIEC_X_EPI")])#0.04339623
  
  #taux id incorrecte EPA AD
  t.EPA_EPI.AD=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPA_EPI") & tabPC2016$Stade%in%"AD"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_EPI","EPA_EPI","ANS_EPI","END_EPI","LIEC_X_EPI")])#0
  
  #taux id incorrecte EPA JV
  t.EPA_EPI.JV=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPA_EPI") & tabPC2016$Stade%in%"JV"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_EPI","EPA_EPI","ANS_EPI","END_EPI","LIEC_X_EPI")])#0.0490566

  
  
  
  #taux id incorrecte ANS SA
  t.ANS_EPI.SA=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("ANS_EPI") & tabPC2016$Stade%in%"SA"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_EPI","EPA_EPI","ANS_EPI","END_EPI","LIEC_X_EPI")])#0
  
  #taux id incorrecte ANS AD
  t.ANS_EPI.AD=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("ANS_EPI") & tabPC2016$Stade%in%"AD"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_EPI","EPA_EPI","ANS_EPI","END_EPI","LIEC_X_EPI")])#0.001886792
  
  #taux id incorrecte ANS JV
  t.ANS_EPI.JV=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("ANS_EPI") & tabPC2016$Stade%in%"JV"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_EPI","EPA_EPI","ANS_EPI","END_EPI","LIEC_X_EPI")])#0
  
  
  
  
  #taux id incorrecte END SA
  t.END_EPI.SA=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("END_EPI") & tabPC2016$Stade%in%"SA"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_EPI","EPA_EPI","ANS_EPI","END_EPI","LIEC_X_EPI")])#0
  
  #taux id incorrecte END AD
  t.END_EPI.AD=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("END_EPI") & tabPC2016$Stade%in%"AD"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_EPI","EPA_EPI","ANS_EPI","END_EPI","LIEC_X_EPI")])#0.001886792
  
  #taux id incorrecte END JV
  t.END_EPI.JV=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("END_EPI") & tabPC2016$Stade%in%"JV"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_EPI","EPA_EPI","ANS_EPI","END_EPI","LIEC_X_EPI")])#0.003773585
  
  
  
  
  
  #taux id incorrecte LIEC_X SA
  t.LIECX_EPI.SA=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("LIEC_X_EPI") & tabPC2016$Stade%in%"SA"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_EPI","EPA_EPI","ANS_EPI","END_EPI","LIEC_X_EPI")])#0
  
  #taux id incorrecte LIEC_X AD
  t.LIECX_EPI.AD=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("LIEC_X_EPI") & tabPC2016$Stade%in%"AD"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_EPI","EPA_EPI","ANS_EPI","END_EPI","LIEC_X_EPI")])#0
  
  #taux id incorrecte LIEC_X JV
  t.LIECX_EPI.JV=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("LIEC_X_EPI") & tabPC2016$Stade%in%"JV"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_EPI","EPA_EPI","ANS_EPI","END_EPI","LIEC_X_EPI")])#0
  
  
  t.EPI_EPI.AD+t.EPI_EPI.SA+t.EPI_EPI.JV+
  t.EPA_EPI.AD+t.EPA_EPI.SA+t.EPA_EPI.JV+
  t.ANS_EPI.AD+t.ANS_EPI.SA+t.ANS_EPI.JV+
  t.END_EPI.AD+t.END_EPI.SA+t.END_EPI.JV+
  t.LIECX_EPI.AD+t.LIECX_EPI.SA+t.LIECX_EPI.JV==1
  
  
}

#les EPA : en quoi ont ils ?t? class?s ?
{
  
  #taux id EPI SA
  t.EPI_EPA.SA=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_EPA") & tabPC2016$Stade%in%"SA"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_EPA","EPA_EPA","ANS_EPA","END_EPA","LIEC_X_EPA")])
  
  #taux id EPI AD
  t.EPI_EPA.AD=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_EPA") & tabPC2016$Stade%in%"AD"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_EPA","EPA_EPA","ANS_EPA","END_EPA","LIEC_X_EPA")])
  
  #taux id EPI JV
  t.EPI_EPA.JV=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_EPA") & tabPC2016$Stade%in%"JV"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_EPA","EPA_EPA","ANS_EPA","END_EPA","LIEC_X_EPA")])
  
  
  
  
  #taux id EPA SA
  t.EPA_EPA.SA=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPA_EPA") & tabPC2016$Stade%in%"SA"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_EPA","EPA_EPA","ANS_EPA","END_EPA","LIEC_X_EPA")])
  
  #taux id EPA AD
  t.EPA_EPA.AD=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPA_EPA") & tabPC2016$Stade%in%"AD"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_EPA","EPA_EPA","ANS_EPA","END_EPA","LIEC_X_EPA")])
  
  #taux id EPA JV
  t.EPA_EPA.JV=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPA_EPA") & tabPC2016$Stade%in%"JV"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_EPA","EPA_EPA","ANS_EPA","END_EPA","LIEC_X_EPA")])
  
  
  
  
  #taux id ANS SA
  t.ANS_EPA.SA=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("ANS_EPA") & tabPC2016$Stade%in%"SA"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_EPA","EPA_EPA","ANS_EPA","END_EPA","LIEC_X_EPA")])
  
  #taux id ANS AD
  t.ANS_EPA.AD=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("ANS_EPA") & tabPC2016$Stade%in%"AD"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_EPA","EPA_EPA","ANS_EPA","END_EPA","LIEC_X_EPA")])
  
  #taux id ANS JV
  t.ANS_EPA.JV=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("ANS_EPA") & tabPC2016$Stade%in%"JV"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_EPA","EPA_EPA","ANS_EPA","END_EPA","LIEC_X_EPA")])
  
  
  
  
  #taux id END SA
  t.END_EPA.SA=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("END_EPA") & tabPC2016$Stade%in%"SA"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_EPA","EPA_EPA","ANS_EPA","END_EPA","LIEC_X_EPA")])
  
  #taux id END AD
  t.END_EPA.AD=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("END_EPA") & tabPC2016$Stade%in%"AD"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_EPA","EPA_EPA","ANS_EPA","END_EPA","LIEC_X_EPA")])
  
  #taux id END JV
  t.END_EPA.JV=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("END_EPA") & tabPC2016$Stade%in%"JV"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_EPA","EPA_EPA","ANS_EPA","END_EPA","LIEC_X_EPA")])
  
  
  
  
  
  #taux id LIEC_X SA
  t.LIECX_EPA.SA=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("LIEC_X_EPA") & tabPC2016$Stade%in%"SA"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_EPA","EPA_EPA","ANS_EPA","END_EPA","LIEC_X_EPA")])
  
  #taux id LIEC_X AD
  t.LIECX_EPA.AD=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("LIEC_X_EPA") & tabPC2016$Stade%in%"AD"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_EPA","EPA_EPA","ANS_EPA","END_EPA","LIEC_X_EPA")])
  
  #taux id LIEC_X JV
  t.LIECX_EPA.JV=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("LIEC_X_EPA") & tabPC2016$Stade%in%"JV"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_EPA","EPA_EPA","ANS_EPA","END_EPA","LIEC_X_EPA")])
  
  
  t.EPI_EPA.AD+t.EPI_EPA.SA+t.EPI_EPA.JV+
  t.EPA_EPA.AD+t.EPA_EPA.SA+t.EPA_EPA.JV+
  t.ANS_EPA.AD+t.ANS_EPA.SA+t.ANS_EPA.JV+
  t.END_EPA.AD+t.END_EPA.SA+t.END_EPA.JV+
  t.LIECX_EPA.AD+t.LIECX_EPA.SA+t.LIECX_EPA.JV==1
  
  
}

#les ANS : en quoi ont ils ?t? class?s ?
{
  
  #taux id EPI SA
  t.EPI_ANS.SA=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_ANS") & tabPC2016$Stade%in%"SA"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_ANS","EPA_ANS","ANS_ANS","END_ANS","LIEC_X_ANS")])
  
  #taux id EPI AD
  t.EPI_ANS.AD=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_ANS") & tabPC2016$Stade%in%"AD"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_ANS","EPA_ANS","ANS_ANS","END_ANS","LIEC_X_ANS")])
  
  #taux id EPI JV
  t.EPI_ANS.JV=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_ANS") & tabPC2016$Stade%in%"JV"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_ANS","EPA_ANS","ANS_ANS","END_ANS","LIEC_X_ANS")])
  
  
  
  
  #taux id EPA SA
  t.EPA_ANS.SA=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPA_ANS") & tabPC2016$Stade%in%"SA"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_ANS","EPA_ANS","ANS_ANS","END_ANS","LIEC_X_ANS")])
  
  #taux id EPA AD
  t.EPA_ANS.AD=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPA_ANS") & tabPC2016$Stade%in%"AD"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_ANS","EPA_ANS","ANS_ANS","END_ANS","LIEC_X_ANS")])
  
  #taux id EPA JV
  t.EPA_ANS.JV=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPA_ANS") & tabPC2016$Stade%in%"JV"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_ANS","EPA_ANS","ANS_ANS","END_ANS","LIEC_X_ANS")])
  
  
  
  
  #taux id ANS SA
  t.ANS_ANS.SA=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("ANS_ANS") & tabPC2016$Stade%in%"SA"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_ANS","EPA_ANS","ANS_ANS","END_ANS","LIEC_X_ANS")])
  
  #taux id ANS AD
  t.ANS_ANS.AD=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("ANS_ANS") & tabPC2016$Stade%in%"AD"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_ANS","EPA_ANS","ANS_ANS","END_ANS","LIEC_X_ANS")])
  
  #taux id ANS JV
  t.ANS_ANS.JV=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("ANS_ANS") & tabPC2016$Stade%in%"JV"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_ANS","EPA_ANS","ANS_ANS","END_ANS","LIEC_X_ANS")])
  
  
  
  
  #taux id END SA
  t.END_ANS.SA=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("END_ANS") & tabPC2016$Stade%in%"SA"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_ANS","EPA_ANS","ANS_ANS","END_ANS","LIEC_X_ANS")])
  
  #taux id END AD
  t.END_ANS.AD=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("END_ANS") & tabPC2016$Stade%in%"AD"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_ANS","EPA_ANS","ANS_ANS","END_ANS","LIEC_X_ANS")])
  
  #taux id END JV
  t.END_ANS.JV=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("END_ANS") & tabPC2016$Stade%in%"JV"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_ANS","EPA_ANS","ANS_ANS","END_ANS","LIEC_X_ANS")])
  
  
  
  
  
  #taux id LIEC_X SA
  t.LIECX_ANS.SA=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("LIEC_X_ANS") & tabPC2016$Stade%in%"SA"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_ANS","EPA_ANS","ANS_ANS","END_ANS","LIEC_X_ANS")])
  
  #taux id LIEC_X AD
  t.LIECX_ANS.AD=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("LIEC_X_ANS") & tabPC2016$Stade%in%"AD"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_ANS","EPA_ANS","ANS_ANS","END_ANS","LIEC_X_ANS")])
  
  #taux id LIEC_X JV
  t.LIECX_ANS.JV=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("LIEC_X_ANS") & tabPC2016$Stade%in%"JV"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_ANS","EPA_ANS","ANS_ANS","END_ANS","LIEC_X_ANS")])
  
  
  t.EPI_ANS.AD+t.EPI_ANS.SA+t.EPI_ANS.JV+
  t.EPA_ANS.AD+t.EPA_ANS.SA+t.EPA_ANS.JV+
  t.ANS_ANS.AD+t.ANS_ANS.SA+t.ANS_ANS.JV+
  t.END_ANS.AD+t.END_ANS.SA+t.END_ANS.JV+
  t.LIECX_ANS.AD+t.LIECX_ANS.SA+t.LIECX_ANS.JV==1
  
  
}

#les END : en quoi ont ils ?t? class?s ?
{
  
  #taux id EPI SA
  t.EPI_END.SA=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_END") & tabPC2016$Stade%in%"SA"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_END","EPA_END","ANS_END","END_END","LIEC_X_END")])
  
  #taux id EPI AD
  t.EPI_END.AD=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_END") & tabPC2016$Stade%in%"AD"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_END","EPA_END","ANS_END","END_END","LIEC_X_END")])
  
  #taux id EPI JV
  t.EPI_END.JV=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_END") & tabPC2016$Stade%in%"JV"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_END","EPA_END","ANS_END","END_END","LIEC_X_END")])
  
  
  
  
  #taux id EPA SA
  t.EPA_END.SA=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPA_END") & tabPC2016$Stade%in%"SA"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_END","EPA_END","ANS_END","END_END","LIEC_X_END")])
  
  #taux id EPA AD
  t.EPA_END.AD=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPA_END") & tabPC2016$Stade%in%"AD"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_END","EPA_END","ANS_END","END_END","LIEC_X_END")])
  
  #taux id EPA JV
  t.EPA_END.JV=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPA_END") & tabPC2016$Stade%in%"JV"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_END","EPA_END","ANS_END","END_END","LIEC_X_END")])
  
  
  
  
  #taux id ANS SA
  t.ANS_END.SA=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("ANS_END") & tabPC2016$Stade%in%"SA"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_END","EPA_END","ANS_END","END_END","LIEC_X_END")])
  
  #taux id ANS AD
  t.ANS_END.AD=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("ANS_END") & tabPC2016$Stade%in%"AD"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_END","EPA_END","ANS_END","END_END","LIEC_X_END")])
  
  #taux id ANS JV
  t.ANS_END.JV=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("ANS_END") & tabPC2016$Stade%in%"JV"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_END","EPA_END","ANS_END","END_END","LIEC_X_END")])
  
  
  
  
  #taux id END SA
  t.END_END.SA=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("END_END") & tabPC2016$Stade%in%"SA"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_END","EPA_END","ANS_END","END_END","LIEC_X_END")])
  
  #taux id END AD
  t.END_END.AD=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("END_END") & tabPC2016$Stade%in%"AD"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_END","EPA_END","ANS_END","END_END","LIEC_X_END")])
  
  #taux id END JV
  t.END_END.JV=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("END_END") & tabPC2016$Stade%in%"JV"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_END","EPA_END","ANS_END","END_END","LIEC_X_END")])
  
  
  
  
  
  #taux id LIEC_X SA
  t.LIECX_END.SA=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("LIEC_X_END") & tabPC2016$Stade%in%"SA"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_END","EPA_END","ANS_END","END_END","LIEC_X_END")])
  
  #taux id LIEC_X AD
  t.LIECX_END.AD=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("LIEC_X_END") & tabPC2016$Stade%in%"AD"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_END","EPA_END","ANS_END","END_END","LIEC_X_END")])
  
  #taux id LIEC_X JV
  t.LIECX_END.JV=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("LIEC_X_END") & tabPC2016$Stade%in%"JV"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_END","EPA_END","ANS_END","END_END","LIEC_X_END")])
  
  
  t.EPI_END.AD+t.EPI_END.SA+t.EPI_END.JV+
  t.EPA_END.AD+t.EPA_END.SA+t.EPA_END.JV+
  t.ANS_END.AD+t.ANS_END.SA+t.ANS_END.JV+
  t.END_END.AD+t.END_END.SA+t.END_END.JV+
  t.LIECX_END.AD+t.LIECX_END.SA+t.LIECX_END.JV==1
  
  
}



#tabGraphPC=data.frame(
#  EPI=c(t.EPI_EPI.AD,t.EPI_EPI.SA,t.EPI_EPI.JV,
#        t.EPA_EPI.AD,t.EPA_EPI.SA,t.EPA_EPI.JV,
#        t.ANS_EPI.AD,t.ANS_EPI.SA,t.ANS_EPI.JV,
#        t.END_EPI.AD,t.END_EPI.SA,t.END_EPI.JV,
#        t.LIECX_EPI.AD,t.LIECX_EPI.SA,t.LIECX_EPI.JV),
#  EPA=c(t.EPI_EPA.AD,t.EPI_EPA.SA,t.EPI_EPA.JV,
#        t.EPA_EPA.AD,t.EPA_EPA.SA,t.EPA_EPA.JV,
#        t.ANS_EPA.AD,t.ANS_EPA.SA,t.ANS_EPA.JV,
#        t.END_EPA.AD,t.END_EPA.SA,t.END_EPA.JV,
#        t.LIECX_EPA.AD,t.LIECX_EPA.SA,t.LIECX_EPA.JV),
#  ANS=c(t.EPI_ANS.AD,t.EPI_ANS.SA,t.EPI_ANS.JV,
#        t.EPA_ANS.AD,t.EPA_ANS.SA,t.EPA_ANS.JV,
#        t.ANS_ANS.AD,t.ANS_ANS.SA,t.ANS_ANS.JV,
#        t.END_ANS.AD,t.END_ANS.SA,t.END_ANS.JV,
#        t.LIECX_ANS.AD,t.LIECX_ANS.SA,t.LIECX_ANS.JV),
#  END=c(t.EPI_END.AD,t.EPI_END.SA,t.EPI_END.JV,
#        t.EPA_END.AD,t.EPA_END.SA,t.EPA_END.JV,
#        t.ANS_END.AD,t.ANS_END.SA,t.ANS_END.JV,
#        t.END_END.AD,t.END_END.SA,t.END_END.JV,
#        t.LIECX_END.AD,t.LIECX_END.SA,t.LIECX_END.JV))
#tabGraphPC=data.frame(t(tabGraphPC))

tabGraphPC=data.frame(
  Taux=c(
  t.EPI_EPI.AD,t.EPI_EPI.SA,t.EPI_EPI.JV,
  t.EPA_EPI.AD,t.EPA_EPI.SA,t.EPA_EPI.JV,
  t.ANS_EPI.AD,t.ANS_EPI.SA,t.ANS_EPI.JV,
  t.END_EPI.AD,t.END_EPI.SA,t.END_EPI.JV,
  t.LIECX_EPI.AD,t.LIECX_EPI.SA,t.LIECX_EPI.JV,
  t.EPI_EPA.AD,t.EPI_EPA.SA,t.EPI_EPA.JV,
  t.EPA_EPA.AD,t.EPA_EPA.SA,t.EPA_EPA.JV,
  t.ANS_EPA.AD,t.ANS_EPA.SA,t.ANS_EPA.JV,
  t.END_EPA.AD,t.END_EPA.SA,t.END_EPA.JV,
  t.LIECX_EPA.AD,t.LIECX_EPA.SA,t.LIECX_EPA.JV,
  t.EPI_ANS.AD,t.EPI_ANS.SA,t.EPI_ANS.JV,
  t.EPA_ANS.AD,t.EPA_ANS.SA,t.EPA_ANS.JV,
  t.ANS_ANS.AD,t.ANS_ANS.SA,t.ANS_ANS.JV,
  t.END_ANS.AD,t.END_ANS.SA,t.END_ANS.JV,
  t.LIECX_ANS.AD,t.LIECX_ANS.SA,t.LIECX_ANS.JV,
  t.EPI_END.AD,t.EPI_END.SA,t.EPI_END.JV,
  t.EPA_END.AD,t.EPA_END.SA,t.EPA_END.JV,
  t.ANS_END.AD,t.ANS_END.SA,t.ANS_END.JV,
  t.END_END.AD,t.END_END.SA,t.END_END.JV,
  t.LIECX_END.AD,t.LIECX_END.SA,t.LIECX_END.JV),
  Classe=c(
    "t.EPI_EPI.AD","t.EPI_EPI.SA","t.EPI_EPI.JV",
    "t.EPA_EPI.AD","t.EPA_EPI.SA","t.EPA_EPI.JV",
    "t.ANS_EPI.AD","t.ANS_EPI.SA","t.ANS_EPI.JV",
    "t.END_EPI.AD","t.END_EPI.SA","t.END_EPI.JV",
    "t.LIECX_EPI.AD","t.LIECX_EPI.SA","t.LIECX_EPI.JV",
    "t.EPI_EPA.AD","t.EPI_EPA.SA","t.EPI_EPA.JV",
    "t.EPA_EPA.AD","t.EPA_EPA.SA","t.EPA_EPA.JV",
    "t.ANS_EPA.AD","t.ANS_EPA.SA","t.ANS_EPA.JV",
    "t.END_EPA.AD","t.END_EPA.SA","t.END_EPA.JV",
    "t.LIECX_EPA.AD","t.LIECX_EPA.SA","t.LIECX_EPA.JV",
    "t.EPI_ANS.AD","t.EPI_ANS.SA","t.EPI_ANS.JV",
    "t.EPA_ANS.AD","t.EPA_ANS.SA","t.EPA_ANS.JV",
    "t.ANS_ANS.AD","t.ANS_ANS.SA","t.ANS_ANS.JV",
    "t.END_ANS.AD","t.END_ANS.SA","t.END_ANS.JV",
    "t.LIECX_ANS.AD","t.LIECX_ANS.SA","t.LIECX_ANS.JV",
    "t.EPI_END.AD","t.EPI_END.SA","t.EPI_END.JV",
    "t.EPA_END.AD","t.EPA_END.SA","t.EPA_END.JV",
    "t.ANS_END.AD","t.ANS_END.SA","t.ANS_END.JV",
    "t.END_END.AD","t.END_END.SA","t.END_END.JV",
    "t.LIECX_END.AD","t.LIECX_END.SA","t.LIECX_END.JV")
  
  )
tabGraphPC$LIEC=factor(substr(tabGraphPC$Classe,7,9),levels=c("EPI","EPA","ANS","END"))
tabGraphPC$Stade=substr(tabGraphPC$Classe,11,12)
tabGraphPC$LIECid=substr(tabGraphPC$Classe,3,5)
tabGraphPC$LIEC_Stade=paste(tabGraphPC$LIEC,tabGraphPC$Stade,sep="_")

tabGraphPC$LIEC_Stade=factor(tabGraphPC$LIEC_Stade,levels=c(
  "EPI_AD","EPI_SA","EPI_JV","EPA_AD","EPA_SA",
  "EPA_JV","ANS_AD","ANS_SA","ANS_JV","END_AD",
  "END_SA","END_JV","LIECX_AD","LIECX_SA","LIECX_JV"))

ggplot(tabGraphPC, aes(fill=LIEC_Stade, y=Taux, x=LIECid)) + 
  geom_bar(position='stack', stat='identity')+
  theme_classic()+
  scale_fill_manual(values = c(
    rgb(0,0,153,maxColorValue = 255),    rgb(0,0,153,maxColorValue = 255,alpha=150),    rgb(0,0,153,maxColorValue = 255,alpha=100),
    rgb(255,255,0,maxColorValue = 255),  rgb(255,255,0,maxColorValue = 255,alpha=150),  rgb(255,255,0,maxColorValue = 255,alpha=100),
    rgb(0,153,0,maxColorValue = 255),    rgb(0,153,0,maxColorValue = 255,alpha=150),    rgb(0,153,0,maxColorValue = 255,alpha=100),
    rgb(255,102,153,maxColorValue = 255),rgb(255,102,153,maxColorValue = 255,alpha=150),rgb(255,102,153,maxColorValue = 255,alpha=100),
    rgb(123,123,123,maxColorValue = 255),rgb(123,123,123,maxColorValue = 255,alpha=150),rgb(123,123,123,maxColorValue = 255,alpha=100)
  ))



##### PC zoom sur les points critiques - Parmis ce que les gens ont class?s en EPI EPA ANS END : qu'est ce qu'on a r?ellement ? #####


#les EPI identifi?s : que sont ils r?ellement ?
{
  
  #taux id correcte EPI SA
  t2.EPI_EPI.SA=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_EPI") & tabPC2016$Stade%in%"SA"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  #taux id correcte EPI AD
  t2.EPI_EPI.AD=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_EPI") & tabPC2016$Stade%in%"AD"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  #taux id correcte EPI JV
  t2.EPI_EPI.JV=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_EPI") & tabPC2016$Stade%in%"JV"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  
  
  
  #taux id incorrecte EPA SA
  t2.EPI_EPA.SA=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_EPA") & tabPC2016$Stade%in%"SA"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  #taux id incorrecte EPA AD
  t2.EPI_EPA.AD=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_EPA") & tabPC2016$Stade%in%"AD"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  #taux id incorrecte EPA JV
  t2.EPI_EPA.JV=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_EPA") & tabPC2016$Stade%in%"JV"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  
  
  
  #taux id incorrecte ANS SA
  t2.EPI_ANS.SA=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_ANS") & tabPC2016$Stade%in%"SA"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  #taux id incorrecte ANS AD
  t2.EPI_ANS.AD=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_ANS") & tabPC2016$Stade%in%"AD"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  #taux id incorrecte ANS JV
  t2.EPI_ANS.JV=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_ANS") & tabPC2016$Stade%in%"JV"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  
  
  
  #taux id incorrecte END SA
  t2.EPI_END.SA=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_END") & tabPC2016$Stade%in%"SA"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  #taux id incorrecte END AD
  t2.EPI_END.AD=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_END") & tabPC2016$Stade%in%"AD"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  #taux id incorrecte END JV
  t2.EPI_END.JV=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_END") & tabPC2016$Stade%in%"JV"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  
  
  
  
  #taux id incorrecte LIEC_X SA
  t2.EPI_LIECX.SA=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_LIEC_X") & tabPC2016$Stade%in%"SA"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  #taux id incorrecte LIEC_X AD
  t2.EPI_LIECX.AD=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_LIEC_X") & tabPC2016$Stade%in%"AD"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  #taux id incorrecte LIEC_X JV
  t2.EPI_LIECX.JV=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_LIEC_X") & tabPC2016$Stade%in%"JV"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  
  t2.EPI_EPI.AD+t2.EPI_EPI.SA+t2.EPI_EPI.JV+
    t2.EPI_EPA.AD+t2.EPI_EPA.SA+t2.EPI_EPA.JV+
    t2.EPI_ANS.AD+t2.EPI_ANS.SA+t2.EPI_ANS.JV+
    t2.EPI_END.AD+t2.EPI_END.SA+t2.EPI_END.JV+
    t2.EPI_LIECX.AD+t2.EPI_LIECX.SA+t2.EPI_LIECX.JV==1
  
  
}


#les EPA identifi?s : que sont ils r?ellement ?
{
  
  #taux id EPI SA
  t2.EPA_EPI.SA=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPA_EPI") & tabPC2016$Stade%in%"SA"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  #taux id EPI AD
  t2.EPA_EPI.AD=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPA_EPI") & tabPC2016$Stade%in%"AD"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  #taux id EPI JV
  t2.EPA_EPI.JV=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPA_EPI") & tabPC2016$Stade%in%"JV"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  
  
  
  #taux id EPA SA
  t2.EPA_EPA.SA=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPA_EPA") & tabPC2016$Stade%in%"SA"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  #taux id EPA AD
  t2.EPA_EPA.AD=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPA_EPA") & tabPC2016$Stade%in%"AD"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  #taux id EPA JV
  t2.EPA_EPA.JV=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPA_EPA") & tabPC2016$Stade%in%"JV"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  
  
  
  #taux id ANS SA
  t2.EPA_ANS.SA=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPA_ANS") & tabPC2016$Stade%in%"SA"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  #taux id ANS AD
  t2.EPA_ANS.AD=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPA_ANS") & tabPC2016$Stade%in%"AD"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  #taux id ANS JV
  t2.EPA_ANS.JV=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPA_ANS") & tabPC2016$Stade%in%"JV"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  
  
  
  #taux id END SA
  t2.EPA_END.SA=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPA_END") & tabPC2016$Stade%in%"SA"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  #taux id END AD
  t2.EPA_END.AD=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPA_END") & tabPC2016$Stade%in%"AD"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  #taux id END JV
  t2.EPA_END.JV=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPA_END") & tabPC2016$Stade%in%"JV"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  
  
  
  
  #taux id LIEC_X SA
  t2.EPA_LIECX.SA=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPA_LIEC_X") & tabPC2016$Stade%in%"SA"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  #taux id LIEC_X AD
  t2.EPA_LIECX.AD=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPA_LIEC_X") & tabPC2016$Stade%in%"AD"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  #taux id LIEC_X JV
  t2.EPA_LIECX.JV=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPA_LIEC_X") & tabPC2016$Stade%in%"JV"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  
  t2.EPA_EPI.AD+t2.EPA_EPI.SA+t2.EPA_EPI.JV+
    t2.EPA_EPA.AD+t2.EPA_EPA.SA+t2.EPA_EPA.JV+
    t2.EPA_ANS.AD+t2.EPA_ANS.SA+t2.EPA_ANS.JV+
    t2.EPA_END.AD+t2.EPA_END.SA+t2.EPA_END.JV+
    t2.EPA_LIECX.AD+t2.EPA_LIECX.SA+t2.EPA_LIECX.JV==1
  
  
}


#les ANS identifi?s : que sont ils r?ellement ?
{
  
  #taux id EPI SA
  t2.ANS_EPI.SA=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("ANS_EPI") & tabPC2016$Stade%in%"SA"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  #taux id EPI AD
  t2.ANS_EPI.AD=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("ANS_EPI") & tabPC2016$Stade%in%"AD"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  #taux id EPI JV
  t2.ANS_EPI.JV=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("ANS_EPI") & tabPC2016$Stade%in%"JV"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  
  
  
  #taux id EPA SA
  t2.ANS_EPA.SA=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("ANS_EPA") & tabPC2016$Stade%in%"SA"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  #taux id EPA AD
  t2.ANS_EPA.AD=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("ANS_EPA") & tabPC2016$Stade%in%"AD"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  #taux id EPA JV
  t2.ANS_EPA.JV=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("ANS_EPA") & tabPC2016$Stade%in%"JV"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  
  
  
  #taux id ANS SA
  t2.ANS_ANS.SA=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("ANS_ANS") & tabPC2016$Stade%in%"SA"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  #taux id ANS AD
  t2.ANS_ANS.AD=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("ANS_ANS") & tabPC2016$Stade%in%"AD"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  #taux id ANS JV
  t2.ANS_ANS.JV=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("ANS_ANS") & tabPC2016$Stade%in%"JV"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  
  
  
  #taux id END SA
  t2.ANS_END.SA=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("ANS_END") & tabPC2016$Stade%in%"SA"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  #taux id END AD
  t2.ANS_END.AD=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("ANS_END") & tabPC2016$Stade%in%"AD"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  #taux id END JV
  t2.ANS_END.JV=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("ANS_END") & tabPC2016$Stade%in%"JV"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  
  
  
  
  #taux id LIEC_X SA
  t2.ANS_LIECX.SA=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("ANS_LIEC_X") & tabPC2016$Stade%in%"SA"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  #taux id LIEC_X AD
  t2.ANS_LIECX.AD=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("ANS_LIEC_X") & tabPC2016$Stade%in%"AD"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  #taux id LIEC_X JV
  t2.ANS_LIECX.JV=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("ANS_LIEC_X") & tabPC2016$Stade%in%"JV"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  
  t2.ANS_EPI.AD+t2.ANS_EPI.SA+t2.ANS_EPI.JV+
    t2.ANS_EPA.AD+t2.ANS_EPA.SA+t2.ANS_EPA.JV+
    t2.ANS_ANS.AD+t2.ANS_ANS.SA+t2.ANS_ANS.JV+
    t2.ANS_END.AD+t2.ANS_END.SA+t2.ANS_END.JV+
    t2.ANS_LIECX.AD+t2.ANS_LIECX.SA+t2.ANS_LIECX.JV==1
  
  
}


#les END identifi?s : que sont ils r?ellement ?
{
  
  #taux id EPI SA
  t2.END_EPI.SA=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("END_EPI") & tabPC2016$Stade%in%"SA"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  #taux id EPI AD
  t2.END_EPI.AD=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("END_EPI") & tabPC2016$Stade%in%"AD"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  #taux id EPI JV
  t2.END_EPI.JV=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("END_EPI") & tabPC2016$Stade%in%"JV"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  
  
  
  #taux id EPA SA
  t2.END_EPA.SA=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("END_EPA") & tabPC2016$Stade%in%"SA"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  #taux id EPA AD
  t2.END_EPA.AD=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("END_EPA") & tabPC2016$Stade%in%"AD"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  #taux id EPA JV
  t2.END_EPA.JV=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("END_EPA") & tabPC2016$Stade%in%"JV"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  
  
  
  #taux id ANS SA
  t2.END_ANS.SA=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("END_ANS") & tabPC2016$Stade%in%"SA"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  #taux id ANS AD
  t2.END_ANS.AD=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("END_ANS") & tabPC2016$Stade%in%"AD"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  #taux id ANS JV
  t2.END_ANS.JV=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("END_ANS") & tabPC2016$Stade%in%"JV"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  
  
  
  #taux id END SA
  t2.END_END.SA=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("END_END") & tabPC2016$Stade%in%"SA"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  #taux id END AD
  t2.END_END.AD=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("END_END") & tabPC2016$Stade%in%"AD"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  #taux id END JV
  t2.END_END.JV=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("END_END") & tabPC2016$Stade%in%"JV"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  
  
  
  
  #taux id LIEC_X SA
  t2.END_LIECX.SA=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("END_LIEC_X") & tabPC2016$Stade%in%"SA"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  #taux id LIEC_X AD
  t2.END_LIECX.AD=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("END_LIEC_X") & tabPC2016$Stade%in%"AD"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  #taux id LIEC_X JV
  t2.END_LIECX.JV=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("END_LIEC_X") & tabPC2016$Stade%in%"JV"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  
  t2.END_EPI.AD+t2.END_EPI.SA+t2.END_EPI.JV+
    t2.END_EPA.AD+t2.END_EPA.SA+t2.END_EPA.JV+
    t2.END_ANS.AD+t2.END_ANS.SA+t2.END_ANS.JV+
    t2.END_END.AD+t2.END_END.SA+t2.END_END.JV+
    t2.END_LIECX.AD+t2.END_LIECX.SA+t2.END_LIECX.JV==1
  
  
}


#les LIECX identifi?s : que sont ils r?ellement ?
{
  
  #taux id EPI SA
  t2.LIECX_EPI.SA=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("LIEC_X_EPI") & tabPC2016$Stade%in%"SA"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  #taux id EPI AD
  t2.LIECX_EPI.AD=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("LIEC_X_EPI") & tabPC2016$Stade%in%"AD"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  #taux id EPI JV
  t2.LIECX_EPI.JV=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("LIEC_X_EPI") & tabPC2016$Stade%in%"JV"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  
  
  
  #taux id EPA SA
  t2.LIECX_EPA.SA=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("LIEC_X_EPA") & tabPC2016$Stade%in%"SA"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  #taux id EPA AD
  t2.LIECX_EPA.AD=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("LIEC_X_EPA") & tabPC2016$Stade%in%"AD"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  #taux id EPA JV
  t2.LIECX_EPA.JV=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("LIEC_X_EPA") & tabPC2016$Stade%in%"JV"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  
  
  
  #taux id ANS SA
  t2.LIECX_ANS.SA=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("LIEC_X_ANS") & tabPC2016$Stade%in%"SA"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  #taux id ANS AD
  t2.LIECX_ANS.AD=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("LIEC_X_ANS") & tabPC2016$Stade%in%"AD"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  #taux id ANS JV
  t2.LIECX_ANS.JV=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("LIEC_X_ANS") & tabPC2016$Stade%in%"JV"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  
  
  
  #taux id END SA
  t2.LIECX_END.SA=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("LIEC_X_END") & tabPC2016$Stade%in%"SA"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  #taux id END AD
  t2.LIECX_END.AD=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("LIEC_X_END") & tabPC2016$Stade%in%"AD"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  #taux id END JV
  t2.LIECX_END.JV=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("LIEC_X_END") & tabPC2016$Stade%in%"JV"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  
  
  
  
  #taux id LIEC_X SA
  t2.LIECX_LIECX.SA=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("LIEC_X_LIEC_X") & tabPC2016$Stade%in%"SA"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  #taux id LIEC_X AD
  t2.LIECX_LIECX.AD=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("LIEC_X_LIEC_X") & tabPC2016$Stade%in%"AD"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  #taux id LIEC_X JV
  t2.LIECX_LIECX.JV=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("LIEC_X_LIEC_X") & tabPC2016$Stade%in%"JV"])/
    sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  
  t2.LIECX_EPI.AD+t2.LIECX_EPI.SA+t2.LIECX_EPI.JV+
    t2.LIECX_EPA.AD+t2.LIECX_EPA.SA+t2.LIECX_EPA.JV+
    t2.LIECX_ANS.AD+t2.LIECX_ANS.SA+t2.LIECX_ANS.JV+
    t2.LIECX_END.AD+t2.LIECX_END.SA+t2.LIECX_END.JV+
    t2.LIECX_LIECX.AD+t2.LIECX_LIECX.SA+t2.LIECX_LIECX.JV==1
  
  
}


tabGraphPC2=data.frame(
  Taux=c(
    t2.EPI_EPI.AD,t2.EPI_EPI.SA,t2.EPI_EPI.JV,
    t2.EPI_EPA.AD,t2.EPI_EPA.SA,t2.EPI_EPA.JV,
    t2.EPI_ANS.AD,t2.EPI_ANS.SA,t2.EPI_ANS.JV,
    t2.EPI_END.AD,t2.EPI_END.SA,t2.EPI_END.JV,
    t2.EPI_LIECX.AD,t2.EPI_LIECX.SA,t2.EPI_LIECX.JV,
    
    t2.EPA_EPI.AD,t2.EPA_EPI.SA,t2.EPA_EPI.JV,
    t2.EPA_EPA.AD,t2.EPA_EPA.SA,t2.EPA_EPA.JV,
    t2.EPA_ANS.AD,t2.EPA_ANS.SA,t2.EPA_ANS.JV,
    t2.EPA_END.AD,t2.EPA_END.SA,t2.EPA_END.JV,
    t2.EPA_LIECX.AD,t2.EPA_LIECX.SA,t2.EPA_LIECX.JV,
    
    t2.ANS_EPI.AD,t2.ANS_EPI.SA,t2.ANS_EPI.JV,
    t2.ANS_EPA.AD,t2.ANS_EPA.SA,t2.ANS_EPA.JV,
    t2.ANS_ANS.AD,t2.ANS_ANS.SA,t2.ANS_ANS.JV,
    t2.ANS_END.AD,t2.ANS_END.SA,t2.ANS_END.JV,
    t2.ANS_LIECX.AD,t2.ANS_LIECX.SA,t2.ANS_LIECX.JV,
    
    t2.END_EPI.AD,t2.END_EPI.SA,t2.END_EPI.JV,
    t2.END_EPA.AD,t2.END_EPA.SA,t2.END_EPA.JV,
    t2.END_ANS.AD,t2.END_ANS.SA,t2.END_ANS.JV,
    t2.END_END.AD,t2.END_END.SA,t2.END_END.JV,
    t2.END_LIECX.AD,t2.END_LIECX.SA,t2.END_LIECX.JV,
    
    t2.LIECX_EPI.AD,t2.LIECX_EPI.SA,t2.LIECX_EPI.JV,
    t2.LIECX_EPA.AD,t2.LIECX_EPA.SA,t2.LIECX_EPA.JV,
    t2.LIECX_ANS.AD,t2.LIECX_ANS.SA,t2.LIECX_ANS.JV,
    t2.LIECX_END.AD,t2.LIECX_END.SA,t2.LIECX_END.JV,
    t2.LIECX_LIECX.AD,t2.LIECX_LIECX.SA,t2.LIECX_LIECX.JV),
  
  Classe=c(
    "t2.EPI_EPI.AD","t2.EPI_EPI.SA","t2.EPI_EPI.JV",
    "t2.EPI_EPA.AD","t2.EPI_EPA.SA","t2.EPI_EPA.JV",
    "t2.EPI_ANS.AD","t2.EPI_ANS.SA","t2.EPI_ANS.JV",
    "t2.EPI_END.AD","t2.EPI_END.SA","t2.EPI_END.JV",
    "t2.EPI_LIECX.AD","t2.EPI_LIECX.SA","t2.EPI_LIECX.JV",
    
    "t2.EPA_EPI.AD","t2.EPA_EPI.SA","t2.EPA_EPI.JV",
    "t2.EPA_EPA.AD","t2.EPA_EPA.SA","t2.EPA_EPA.JV",
    "t2.EPA_ANS.AD","t2.EPA_ANS.SA","t2.EPA_ANS.JV",
    "t2.EPA_END.AD","t2.EPA_END.SA","t2.EPA_END.JV",
    "t2.EPA_LIECX.AD","t2.EPA_LIECX.SA","t2.EPA_LIECX.JV",
    
    "t2.ANS_EPI.AD","t2.ANS_EPI.SA","t2.ANS_EPI.JV",
    "t2.ANS_EPA.AD","t2.ANS_EPA.SA","t2.ANS_EPA.JV",
    "t2.ANS_ANS.AD","t2.ANS_ANS.SA","t2.ANS_ANS.JV",
    "t2.ANS_END.AD","t2.ANS_END.SA","t2.ANS_END.JV",
    "t2.ANS_LIECX.AD","t2.ANS_LIECX.SA","t2.ANS_LIECX.JV",
    
    "t2.END_EPI.AD","t2.END_EPI.SA","t2.END_EPI.JV",
    "t2.END_EPA.AD","t2.END_EPA.SA","t2.END_EPA.JV",
    "t2.END_ANS.AD","t2.END_ANS.SA","t2.END_ANS.JV",
    "t2.END_END.AD","t2.END_END.SA","t2.END_END.JV",
    "t2.END_LIECX.AD","t2.END_LIECX.SA","t2.END_LIECX.JV",
    
    "t2.LIECX_EPI.AD","t2.LIECX_EPI.SA","t2.LIECX_EPI.JV",
    "t2.LIECX_EPA.AD","t2.LIECX_EPA.SA","t2.LIECX_EPA.JV",
    "t2.LIECX_ANS.AD","t2.LIECX_ANS.SA","t2.LIECX_ANS.JV",
    "t2.LIECX_END.AD","t2.LIECX_END.SA","t2.LIECX_END.JV",
    "t2.LIECX_LIECX.AD","t2.LIECX_LIECX.SA","t2.LIECX_LIECX.JV")
  
)
tabGraphPC2$LIEC=factor(substr(tabGraphPC2$Classe,8,10),levels=c("EPI","EPA","ANS","END","LIECX"))
tabGraphPC2$Stade=substr(tabGraphPC2$Classe,12,13)
tabGraphPC2$LIECid=factor(substr(tabGraphPC2$Classe,4,6),levels=c("EPI","EPA","ANS","END","LIECX"))
tabGraphPC2$LIEC_Stade=paste(tabGraphPC2$LIEC,tabGraphPC2$Stade,sep="_")

tabGraphPC2$LIEC_Stade=factor(tabGraphPC2$LIEC_Stade,levels=c(
  "EPI_AD","EPI_SA","EPI_JV","EPA_AD","EPA_SA",
  "EPA_JV","ANS_AD","ANS_SA","ANS_JV","END_AD",
  "END_SA","END_JV","LIECX_AD","LIECX_SA","LIECX_JV"))

ggplot(tabGraphPC2, aes(fill=LIEC_Stade, y=Taux, x=LIECid)) + 
  geom_bar(position='stack', stat='identity')+
  theme_classic()+
  scale_fill_manual(values = c(
    rgb(0,0,153,maxColorValue = 255),    rgb(0,0,153,maxColorValue = 255,alpha=150),    rgb(0,0,153,maxColorValue = 255,alpha=100),
    rgb(255,255,0,maxColorValue = 255),  rgb(255,255,0,maxColorValue = 255,alpha=100),  rgb(255,255,0,maxColorValue = 255,alpha=25),
    rgb(0,153,0,maxColorValue = 255),    rgb(0,153,0,maxColorValue = 255,alpha=150),    rgb(0,153,0,maxColorValue = 255,alpha=100),
    rgb(255,102,153,maxColorValue = 255),rgb(255,102,153,maxColorValue = 255,alpha=150),rgb(255,102,153,maxColorValue = 255,alpha=100),
    rgb(123,123,123,maxColorValue = 255),rgb(123,123,123,maxColorValue = 255,alpha=150),rgb(123,123,123,maxColorValue = 255,alpha=100)
  ))

#zoom sur les erreurs
tabGraphPC2$Erreur=NA
tabGraphPC2$Erreur[tabGraphPC2$LIEC==tabGraphPC2$LIECid]="Correcte"
tabGraphPC2$Erreur[tabGraphPC2$LIEC!=tabGraphPC2$LIECid]="Incorrecte"

ggplot(tabGraphPC2[tabGraphPC2$Erreur=="Incorrecte",], aes(fill=LIEC_Stade, y=Taux, x=LIECid)) + 
  geom_bar(position='fill', stat='identity')+
  theme_classic()+
  scale_fill_manual(values = c(
    rgb(0,0,153,maxColorValue = 255),    rgb(0,0,153,maxColorValue = 255,alpha=150),    rgb(0,0,153,maxColorValue = 255,alpha=100),
    rgb(255,255,0,maxColorValue = 255),  rgb(255,255,0,maxColorValue = 255,alpha=100),  rgb(255,255,0,maxColorValue = 255,alpha=25),
    rgb(0,153,0,maxColorValue = 255),    rgb(0,153,0,maxColorValue = 255,alpha=150),    rgb(0,153,0,maxColorValue = 255,alpha=100),
    rgb(255,102,153,maxColorValue = 255),rgb(255,102,153,maxColorValue = 255,alpha=150),rgb(255,102,153,maxColorValue = 255,alpha=100),
    rgb(123,123,123,maxColorValue = 255),rgb(123,123,123,maxColorValue = 255,alpha=150),rgb(123,123,123,maxColorValue = 255,alpha=100)
  ))


##### AQ zoom sur les points critiques - Parmis ce que les gens ont class?s en EPI EPA ANS END : qu'est ce qu'on a r?ellement ? #####


tabAQ=Tab2[Tab2$Region%in%"AQ",]
nrow(tabAQ)==1959
dataID <- unique(sbt_m[,c("cle",colID,"Region")])
nrow(dataID)==1109
dataID=dataID[dataID$Region%in%"AQ",]
nrow(dataID)==129



tabAQ2014=tabAQ[tabAQ$Annee%in%2014,]
unique(tabAQ2014$Stade)#"JV" "SA" "AD" => ok
any(is.na(tabAQ2014$Stade))==F

#taux d'id correcte global
sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("EPI_EPI","EPA_EPA","ANS_ANS","END_END")])/
  sum(tabAQ2014$Nbr_EW)#0.4466258 => ok
#vectReussiteAQ


#les EPI identifi?s : que sont ils r?ellement ?
{
  
  #taux id correcte EPI SA
  t2.EPI_EPI.SA=sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("EPI_EPI") & tabAQ2014$Stade%in%"SA"])/
    sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  #taux id correcte EPI AD
  t2.EPI_EPI.AD=sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("EPI_EPI") & tabAQ2014$Stade%in%"AD"])/
    sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  #taux id correcte EPI JV
  t2.EPI_EPI.JV=sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("EPI_EPI") & tabAQ2014$Stade%in%"JV"])/
    sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  
  
  
  #taux id incorrecte EPA SA
  t2.EPI_EPA.SA=sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("EPI_EPA") & tabAQ2014$Stade%in%"SA"])/
    sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  #taux id incorrecte EPA AD
  t2.EPI_EPA.AD=sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("EPI_EPA") & tabAQ2014$Stade%in%"AD"])/
    sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  #taux id incorrecte EPA JV
  t2.EPI_EPA.JV=sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("EPI_EPA") & tabAQ2014$Stade%in%"JV"])/
    sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  
  
  
  #taux id incorrecte ANS SA
  t2.EPI_ANS.SA=sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("EPI_ANS") & tabAQ2014$Stade%in%"SA"])/
    sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  #taux id incorrecte ANS AD
  t2.EPI_ANS.AD=sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("EPI_ANS") & tabAQ2014$Stade%in%"AD"])/
    sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  #taux id incorrecte ANS JV
  t2.EPI_ANS.JV=sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("EPI_ANS") & tabAQ2014$Stade%in%"JV"])/
    sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  
  
  
  #taux id incorrecte END SA
  t2.EPI_END.SA=sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("EPI_END") & tabAQ2014$Stade%in%"SA"])/
    sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  #taux id incorrecte END AD
  t2.EPI_END.AD=sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("EPI_END") & tabAQ2014$Stade%in%"AD"])/
    sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  #taux id incorrecte END JV
  t2.EPI_END.JV=sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("EPI_END") & tabAQ2014$Stade%in%"JV"])/
    sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  
  
  
  
  #taux id incorrecte LIEC_X SA
  t2.EPI_LIECX.SA=sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("EPI_LIEC_X") & tabAQ2014$Stade%in%"SA"])/
    sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  #taux id incorrecte LIEC_X AD
  t2.EPI_LIECX.AD=sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("EPI_LIEC_X") & tabAQ2014$Stade%in%"AD"])/
    sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  #taux id incorrecte LIEC_X JV
  t2.EPI_LIECX.JV=sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("EPI_LIEC_X") & tabAQ2014$Stade%in%"JV"])/
    sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  
  t2.EPI_EPI.AD+t2.EPI_EPI.SA+t2.EPI_EPI.JV+
    t2.EPI_EPA.AD+t2.EPI_EPA.SA+t2.EPI_EPA.JV+
    t2.EPI_ANS.AD+t2.EPI_ANS.SA+t2.EPI_ANS.JV+
    t2.EPI_END.AD+t2.EPI_END.SA+t2.EPI_END.JV+
    t2.EPI_LIECX.AD+t2.EPI_LIECX.SA+t2.EPI_LIECX.JV==1
  
  
  }


#les EPA identifi?s : que sont ils r?ellement ?
{
  
  #taux id EPI SA
  t2.EPA_EPI.SA=sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("EPA_EPI") & tabAQ2014$Stade%in%"SA"])/
    sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  #taux id EPI AD
  t2.EPA_EPI.AD=sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("EPA_EPI") & tabAQ2014$Stade%in%"AD"])/
    sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  #taux id EPI JV
  t2.EPA_EPI.JV=sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("EPA_EPI") & tabAQ2014$Stade%in%"JV"])/
    sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  
  
  
  #taux id EPA SA
  t2.EPA_EPA.SA=sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("EPA_EPA") & tabAQ2014$Stade%in%"SA"])/
    sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  #taux id EPA AD
  t2.EPA_EPA.AD=sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("EPA_EPA") & tabAQ2014$Stade%in%"AD"])/
    sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  #taux id EPA JV
  t2.EPA_EPA.JV=sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("EPA_EPA") & tabAQ2014$Stade%in%"JV"])/
    sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  
  
  
  #taux id ANS SA
  t2.EPA_ANS.SA=sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("EPA_ANS") & tabAQ2014$Stade%in%"SA"])/
    sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  #taux id ANS AD
  t2.EPA_ANS.AD=sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("EPA_ANS") & tabAQ2014$Stade%in%"AD"])/
    sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  #taux id ANS JV
  t2.EPA_ANS.JV=sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("EPA_ANS") & tabAQ2014$Stade%in%"JV"])/
    sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  
  
  
  #taux id END SA
  t2.EPA_END.SA=sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("EPA_END") & tabAQ2014$Stade%in%"SA"])/
    sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  #taux id END AD
  t2.EPA_END.AD=sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("EPA_END") & tabAQ2014$Stade%in%"AD"])/
    sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  #taux id END JV
  t2.EPA_END.JV=sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("EPA_END") & tabAQ2014$Stade%in%"JV"])/
    sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  
  
  
  
  #taux id LIEC_X SA
  t2.EPA_LIECX.SA=sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("EPA_LIEC_X") & tabAQ2014$Stade%in%"SA"])/
    sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  #taux id LIEC_X AD
  t2.EPA_LIECX.AD=sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("EPA_LIEC_X") & tabAQ2014$Stade%in%"AD"])/
    sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  #taux id LIEC_X JV
  t2.EPA_LIECX.JV=sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("EPA_LIEC_X") & tabAQ2014$Stade%in%"JV"])/
    sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  
  round(t2.EPA_EPI.AD+t2.EPA_EPI.SA+t2.EPA_EPI.JV+
    t2.EPA_EPA.AD+t2.EPA_EPA.SA+t2.EPA_EPA.JV+
    t2.EPA_ANS.AD+t2.EPA_ANS.SA+t2.EPA_ANS.JV+
    t2.EPA_END.AD+t2.EPA_END.SA+t2.EPA_END.JV+
    t2.EPA_LIECX.AD+t2.EPA_LIECX.SA+t2.EPA_LIECX.JV,5)==1
  
  
}


#les ANS identifi?s : que sont ils r?ellement ?
{
  
  #taux id EPI SA
  t2.ANS_EPI.SA=sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("ANS_EPI") & tabAQ2014$Stade%in%"SA"])/
    sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  #taux id EPI AD
  t2.ANS_EPI.AD=sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("ANS_EPI") & tabAQ2014$Stade%in%"AD"])/
    sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  #taux id EPI JV
  t2.ANS_EPI.JV=sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("ANS_EPI") & tabAQ2014$Stade%in%"JV"])/
    sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  
  
  
  #taux id EPA SA
  t2.ANS_EPA.SA=sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("ANS_EPA") & tabAQ2014$Stade%in%"SA"])/
    sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  #taux id EPA AD
  t2.ANS_EPA.AD=sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("ANS_EPA") & tabAQ2014$Stade%in%"AD"])/
    sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  #taux id EPA JV
  t2.ANS_EPA.JV=sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("ANS_EPA") & tabAQ2014$Stade%in%"JV"])/
    sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  
  
  
  #taux id ANS SA
  t2.ANS_ANS.SA=sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("ANS_ANS") & tabAQ2014$Stade%in%"SA"])/
    sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  #taux id ANS AD
  t2.ANS_ANS.AD=sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("ANS_ANS") & tabAQ2014$Stade%in%"AD"])/
    sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  #taux id ANS JV
  t2.ANS_ANS.JV=sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("ANS_ANS") & tabAQ2014$Stade%in%"JV"])/
    sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  
  
  
  #taux id END SA
  t2.ANS_END.SA=sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("ANS_END") & tabAQ2014$Stade%in%"SA"])/
    sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  #taux id END AD
  t2.ANS_END.AD=sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("ANS_END") & tabAQ2014$Stade%in%"AD"])/
    sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  #taux id END JV
  t2.ANS_END.JV=sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("ANS_END") & tabAQ2014$Stade%in%"JV"])/
    sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  
  
  
  
  #taux id LIEC_X SA
  t2.ANS_LIECX.SA=sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("ANS_LIEC_X") & tabAQ2014$Stade%in%"SA"])/
    sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  #taux id LIEC_X AD
  t2.ANS_LIECX.AD=sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("ANS_LIEC_X") & tabAQ2014$Stade%in%"AD"])/
    sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  #taux id LIEC_X JV
  t2.ANS_LIECX.JV=sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("ANS_LIEC_X") & tabAQ2014$Stade%in%"JV"])/
    sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  
  round(t2.ANS_EPI.AD+t2.ANS_EPI.SA+t2.ANS_EPI.JV+
    t2.ANS_EPA.AD+t2.ANS_EPA.SA+t2.ANS_EPA.JV+
    t2.ANS_ANS.AD+t2.ANS_ANS.SA+t2.ANS_ANS.JV+
    t2.ANS_END.AD+t2.ANS_END.SA+t2.ANS_END.JV+
    t2.ANS_LIECX.AD+t2.ANS_LIECX.SA+t2.ANS_LIECX.JV,5)==1
  
  
}


#les END identifi?s : que sont ils r?ellement ?
{
  
  #taux id EPI SA
  t2.END_EPI.SA=sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("END_EPI") & tabAQ2014$Stade%in%"SA"])/
    sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  #taux id EPI AD
  t2.END_EPI.AD=sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("END_EPI") & tabAQ2014$Stade%in%"AD"])/
    sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  #taux id EPI JV
  t2.END_EPI.JV=sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("END_EPI") & tabAQ2014$Stade%in%"JV"])/
    sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  
  
  
  #taux id EPA SA
  t2.END_EPA.SA=sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("END_EPA") & tabAQ2014$Stade%in%"SA"])/
    sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  #taux id EPA AD
  t2.END_EPA.AD=sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("END_EPA") & tabAQ2014$Stade%in%"AD"])/
    sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  #taux id EPA JV
  t2.END_EPA.JV=sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("END_EPA") & tabAQ2014$Stade%in%"JV"])/
    sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  
  
  
  #taux id ANS SA
  t2.END_ANS.SA=sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("END_ANS") & tabAQ2014$Stade%in%"SA"])/
    sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  #taux id ANS AD
  t2.END_ANS.AD=sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("END_ANS") & tabAQ2014$Stade%in%"AD"])/
    sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  #taux id ANS JV
  t2.END_ANS.JV=sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("END_ANS") & tabAQ2014$Stade%in%"JV"])/
    sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  
  
  
  #taux id END SA
  t2.END_END.SA=sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("END_END") & tabAQ2014$Stade%in%"SA"])/
    sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  #taux id END AD
  t2.END_END.AD=sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("END_END") & tabAQ2014$Stade%in%"AD"])/
    sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  #taux id END JV
  t2.END_END.JV=sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("END_END") & tabAQ2014$Stade%in%"JV"])/
    sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  
  
  
  
  #taux id LIEC_X SA
  t2.END_LIECX.SA=sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("END_LIEC_X") & tabAQ2014$Stade%in%"SA"])/
    sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  #taux id LIEC_X AD
  t2.END_LIECX.AD=sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("END_LIEC_X") & tabAQ2014$Stade%in%"AD"])/
    sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  #taux id LIEC_X JV
  t2.END_LIECX.JV=sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("END_LIEC_X") & tabAQ2014$Stade%in%"JV"])/
    sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  
  round(t2.END_EPI.AD+t2.END_EPI.SA+t2.END_EPI.JV+
    t2.END_EPA.AD+t2.END_EPA.SA+t2.END_EPA.JV+
    t2.END_ANS.AD+t2.END_ANS.SA+t2.END_ANS.JV+
    t2.END_END.AD+t2.END_END.SA+t2.END_END.JV+
    t2.END_LIECX.AD+t2.END_LIECX.SA+t2.END_LIECX.JV,5)==1
  
  
}


#les LIECX identifi?s : que sont ils r?ellement ?
{
  
  #taux id EPI SA
  t2.LIECX_EPI.SA=sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("LIEC_X_EPI") & tabAQ2014$Stade%in%"SA"])/
    sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  #taux id EPI AD
  t2.LIECX_EPI.AD=sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("LIEC_X_EPI") & tabAQ2014$Stade%in%"AD"])/
    sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  #taux id EPI JV
  t2.LIECX_EPI.JV=sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("LIEC_X_EPI") & tabAQ2014$Stade%in%"JV"])/
    sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  
  
  
  #taux id EPA SA
  t2.LIECX_EPA.SA=sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("LIEC_X_EPA") & tabAQ2014$Stade%in%"SA"])/
    sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  #taux id EPA AD
  t2.LIECX_EPA.AD=sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("LIEC_X_EPA") & tabAQ2014$Stade%in%"AD"])/
    sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  #taux id EPA JV
  t2.LIECX_EPA.JV=sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("LIEC_X_EPA") & tabAQ2014$Stade%in%"JV"])/
    sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  
  
  
  #taux id ANS SA
  t2.LIECX_ANS.SA=sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("LIEC_X_ANS") & tabAQ2014$Stade%in%"SA"])/
    sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  #taux id ANS AD
  t2.LIECX_ANS.AD=sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("LIEC_X_ANS") & tabAQ2014$Stade%in%"AD"])/
    sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  #taux id ANS JV
  t2.LIECX_ANS.JV=sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("LIEC_X_ANS") & tabAQ2014$Stade%in%"JV"])/
    sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  
  
  
  #taux id END SA
  t2.LIECX_END.SA=sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("LIEC_X_END") & tabAQ2014$Stade%in%"SA"])/
    sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  #taux id END AD
  t2.LIECX_END.AD=sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("LIEC_X_END") & tabAQ2014$Stade%in%"AD"])/
    sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  #taux id END JV
  t2.LIECX_END.JV=sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("LIEC_X_END") & tabAQ2014$Stade%in%"JV"])/
    sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  
  
  
  
  #taux id LIEC_X SA
  t2.LIECX_LIECX.SA=sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("LIEC_X_LIEC_X") & tabAQ2014$Stade%in%"SA"])/
    sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  #taux id LIEC_X AD
  t2.LIECX_LIECX.AD=sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("LIEC_X_LIEC_X") & tabAQ2014$Stade%in%"AD"])/
    sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  #taux id LIEC_X JV
  t2.LIECX_LIECX.JV=sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("LIEC_X_LIEC_X") & tabAQ2014$Stade%in%"JV"])/
    sum(tabAQ2014$Nbr_EW[tabAQ2014$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  
  round(t2.LIECX_EPI.AD+t2.LIECX_EPI.SA+t2.LIECX_EPI.JV+
    t2.LIECX_EPA.AD+t2.LIECX_EPA.SA+t2.LIECX_EPA.JV+
    t2.LIECX_ANS.AD+t2.LIECX_ANS.SA+t2.LIECX_ANS.JV+
    t2.LIECX_END.AD+t2.LIECX_END.SA+t2.LIECX_END.JV+
    t2.LIECX_LIECX.AD+t2.LIECX_LIECX.SA+t2.LIECX_LIECX.JV,5)==1
  
  
}


tabGraphAQ2=data.frame(
  Taux=c(
    t2.EPI_EPI.AD,t2.EPI_EPI.SA,t2.EPI_EPI.JV,
    t2.EPI_EPA.AD,t2.EPI_EPA.SA,t2.EPI_EPA.JV,
    t2.EPI_ANS.AD,t2.EPI_ANS.SA,t2.EPI_ANS.JV,
    t2.EPI_END.AD,t2.EPI_END.SA,t2.EPI_END.JV,
    t2.EPI_LIECX.AD,t2.EPI_LIECX.SA,t2.EPI_LIECX.JV,
    
    t2.EPA_EPI.AD,t2.EPA_EPI.SA,t2.EPA_EPI.JV,
    t2.EPA_EPA.AD,t2.EPA_EPA.SA,t2.EPA_EPA.JV,
    t2.EPA_ANS.AD,t2.EPA_ANS.SA,t2.EPA_ANS.JV,
    t2.EPA_END.AD,t2.EPA_END.SA,t2.EPA_END.JV,
    t2.EPA_LIECX.AD,t2.EPA_LIECX.SA,t2.EPA_LIECX.JV,
    
    t2.ANS_EPI.AD,t2.ANS_EPI.SA,t2.ANS_EPI.JV,
    t2.ANS_EPA.AD,t2.ANS_EPA.SA,t2.ANS_EPA.JV,
    t2.ANS_ANS.AD,t2.ANS_ANS.SA,t2.ANS_ANS.JV,
    t2.ANS_END.AD,t2.ANS_END.SA,t2.ANS_END.JV,
    t2.ANS_LIECX.AD,t2.ANS_LIECX.SA,t2.ANS_LIECX.JV,
    
    t2.END_EPI.AD,t2.END_EPI.SA,t2.END_EPI.JV,
    t2.END_EPA.AD,t2.END_EPA.SA,t2.END_EPA.JV,
    t2.END_ANS.AD,t2.END_ANS.SA,t2.END_ANS.JV,
    t2.END_END.AD,t2.END_END.SA,t2.END_END.JV,
    t2.END_LIECX.AD,t2.END_LIECX.SA,t2.END_LIECX.JV,
    
    t2.LIECX_EPI.AD,t2.LIECX_EPI.SA,t2.LIECX_EPI.JV,
    t2.LIECX_EPA.AD,t2.LIECX_EPA.SA,t2.LIECX_EPA.JV,
    t2.LIECX_ANS.AD,t2.LIECX_ANS.SA,t2.LIECX_ANS.JV,
    t2.LIECX_END.AD,t2.LIECX_END.SA,t2.LIECX_END.JV,
    t2.LIECX_LIECX.AD,t2.LIECX_LIECX.SA,t2.LIECX_LIECX.JV),
  
  Classe=c(
    "t2.EPI_EPI.AD","t2.EPI_EPI.SA","t2.EPI_EPI.JV",
    "t2.EPI_EPA.AD","t2.EPI_EPA.SA","t2.EPI_EPA.JV",
    "t2.EPI_ANS.AD","t2.EPI_ANS.SA","t2.EPI_ANS.JV",
    "t2.EPI_END.AD","t2.EPI_END.SA","t2.EPI_END.JV",
    "t2.EPI_LIECX.AD","t2.EPI_LIECX.SA","t2.EPI_LIECX.JV",
    
    "t2.EPA_EPI.AD","t2.EPA_EPI.SA","t2.EPA_EPI.JV",
    "t2.EPA_EPA.AD","t2.EPA_EPA.SA","t2.EPA_EPA.JV",
    "t2.EPA_ANS.AD","t2.EPA_ANS.SA","t2.EPA_ANS.JV",
    "t2.EPA_END.AD","t2.EPA_END.SA","t2.EPA_END.JV",
    "t2.EPA_LIECX.AD","t2.EPA_LIECX.SA","t2.EPA_LIECX.JV",
    
    "t2.ANS_EPI.AD","t2.ANS_EPI.SA","t2.ANS_EPI.JV",
    "t2.ANS_EPA.AD","t2.ANS_EPA.SA","t2.ANS_EPA.JV",
    "t2.ANS_ANS.AD","t2.ANS_ANS.SA","t2.ANS_ANS.JV",
    "t2.ANS_END.AD","t2.ANS_END.SA","t2.ANS_END.JV",
    "t2.ANS_LIECX.AD","t2.ANS_LIECX.SA","t2.ANS_LIECX.JV",
    
    "t2.END_EPI.AD","t2.END_EPI.SA","t2.END_EPI.JV",
    "t2.END_EPA.AD","t2.END_EPA.SA","t2.END_EPA.JV",
    "t2.END_ANS.AD","t2.END_ANS.SA","t2.END_ANS.JV",
    "t2.END_END.AD","t2.END_END.SA","t2.END_END.JV",
    "t2.END_LIECX.AD","t2.END_LIECX.SA","t2.END_LIECX.JV",
    
    "t2.LIECX_EPI.AD","t2.LIECX_EPI.SA","t2.LIECX_EPI.JV",
    "t2.LIECX_EPA.AD","t2.LIECX_EPA.SA","t2.LIECX_EPA.JV",
    "t2.LIECX_ANS.AD","t2.LIECX_ANS.SA","t2.LIECX_ANS.JV",
    "t2.LIECX_END.AD","t2.LIECX_END.SA","t2.LIECX_END.JV",
    "t2.LIECX_LIECX.AD","t2.LIECX_LIECX.SA","t2.LIECX_LIECX.JV")
  
)
tabGraphAQ2$LIEC=factor(substr(tabGraphAQ2$Classe,8,10),levels=c("EPI","EPA","ANS","END","LIECX"))
tabGraphAQ2$Stade=substr(tabGraphAQ2$Classe,12,13)
tabGraphAQ2$LIECid=factor(substr(tabGraphAQ2$Classe,4,6),levels=c("EPI","EPA","ANS","END","LIECX"))
tabGraphAQ2$LIEC_Stade=paste(tabGraphAQ2$LIEC,tabGraphAQ2$Stade,sep="_")

tabGraphAQ2$LIEC_Stade=factor(tabGraphAQ2$LIEC_Stade,levels=c(
  "EPI_AD","EPI_SA","EPI_JV","EPA_AD","EPA_SA",
  "EPA_JV","ANS_AD","ANS_SA","ANS_JV","END_AD",
  "END_SA","END_JV","LIECX_AD","LIECX_SA","LIECX_JV"))

ggplot(tabGraphAQ2, aes(fill=LIEC_Stade, y=Taux, x=LIECid)) + 
  geom_bar(position='stack', stat='identity')+
  theme_classic()+
  scale_fill_manual(values = c(
    rgb(0,0,153,maxColorValue = 255),    rgb(0,0,153,maxColorValue = 255,alpha=150),    rgb(0,0,153,maxColorValue = 255,alpha=100),
    rgb(255,255,0,maxColorValue = 255),  rgb(255,255,0,maxColorValue = 255,alpha=100),  rgb(255,255,0,maxColorValue = 255,alpha=25),
    rgb(0,153,0,maxColorValue = 255),    rgb(0,153,0,maxColorValue = 255,alpha=150),    rgb(0,153,0,maxColorValue = 255,alpha=100),
    rgb(255,102,153,maxColorValue = 255),rgb(255,102,153,maxColorValue = 255,alpha=150),rgb(255,102,153,maxColorValue = 255,alpha=100),
    rgb(123,123,123,maxColorValue = 255),rgb(123,123,123,maxColorValue = 255,alpha=150),rgb(123,123,123,maxColorValue = 255,alpha=100)
  ))

#zoom sur les erreurs
tabGraphAQ2$Erreur=NA
tabGraphAQ2$Erreur[tabGraphAQ2$LIEC==tabGraphAQ2$LIECid]="Correcte"
tabGraphAQ2$Erreur[tabGraphAQ2$LIEC!=tabGraphAQ2$LIECid]="Incorrecte"

ggplot(tabGraphAQ2[tabGraphAQ2$Erreur=="Incorrecte",], aes(fill=LIEC_Stade, y=Taux, x=LIECid)) + 
  geom_bar(position='fill', stat='identity')+
  theme_classic()+
  scale_fill_manual(values = c(
    rgb(0,0,153,maxColorValue = 255),    rgb(0,0,153,maxColorValue = 255,alpha=150),    rgb(0,0,153,maxColorValue = 255,alpha=100),
    rgb(255,255,0,maxColorValue = 255),  rgb(255,255,0,maxColorValue = 255,alpha=100),  rgb(255,255,0,maxColorValue = 255,alpha=25),
    rgb(0,153,0,maxColorValue = 255),    rgb(0,153,0,maxColorValue = 255,alpha=150),    rgb(0,153,0,maxColorValue = 255,alpha=100),
    rgb(255,102,153,maxColorValue = 255),rgb(255,102,153,maxColorValue = 255,alpha=150),rgb(255,102,153,maxColorValue = 255,alpha=100),
    rgb(123,123,123,maxColorValue = 255),rgb(123,123,123,maxColorValue = 255,alpha=150),rgb(123,123,123,maxColorValue = 255,alpha=100)
  ))




##### /!\pas de LIEC_X identifi?s/!\ FC zoom sur les points critiques - Parmis ce que les gens ont class?s en EPI EPA ANS END : qu'est ce qu'on a r?ellement ? #####


tabFC=Tab2[Tab2$Region%in%"FC",]
nrow(tabFC)==1662
dataID <- unique(sbt_m[,c("cle",colID,"Region")])
nrow(dataID)==1109
dataID=dataID[dataID$Region%in%"FC",]
nrow(dataID)==57



tabFC2016=tabFC[tabFC$Annee%in%2016,]
unique(tabFC2016$Stade)#"JV" "SA" "AD" => ok
any(is.na(tabFC2016$Stade))==F

#taux d'id correcte global
sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("EPI_EPI","EPA_EPA","ANS_ANS","END_END")])/
  sum(tabFC2016$Nbr_EW)#0.5334773  => ok
#vectReussiteFC


#les EPI identifi?s : que sont ils r?ellement ?
{
  
  #taux id correcte EPI SA
  t2.EPI_EPI.SA=sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("EPI_EPI") & tabFC2016$Stade%in%"SA"])/
    sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  #taux id correcte EPI AD
  t2.EPI_EPI.AD=sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("EPI_EPI") & tabFC2016$Stade%in%"AD"])/
    sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  #taux id correcte EPI JV
  t2.EPI_EPI.JV=sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("EPI_EPI") & tabFC2016$Stade%in%"JV"])/
    sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  
  
  
  #taux id incorrecte EPA SA
  t2.EPI_EPA.SA=sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("EPI_EPA") & tabFC2016$Stade%in%"SA"])/
    sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  #taux id incorrecte EPA AD
  t2.EPI_EPA.AD=sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("EPI_EPA") & tabFC2016$Stade%in%"AD"])/
    sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  #taux id incorrecte EPA JV
  t2.EPI_EPA.JV=sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("EPI_EPA") & tabFC2016$Stade%in%"JV"])/
    sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  
  
  
  #taux id incorrecte ANS SA
  t2.EPI_ANS.SA=sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("EPI_ANS") & tabFC2016$Stade%in%"SA"])/
    sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  #taux id incorrecte ANS AD
  t2.EPI_ANS.AD=sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("EPI_ANS") & tabFC2016$Stade%in%"AD"])/
    sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  #taux id incorrecte ANS JV
  t2.EPI_ANS.JV=sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("EPI_ANS") & tabFC2016$Stade%in%"JV"])/
    sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  
  
  
  #taux id incorrecte END SA
  t2.EPI_END.SA=sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("EPI_END") & tabFC2016$Stade%in%"SA"])/
    sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  #taux id incorrecte END AD
  t2.EPI_END.AD=sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("EPI_END") & tabFC2016$Stade%in%"AD"])/
    sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  #taux id incorrecte END JV
  t2.EPI_END.JV=sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("EPI_END") & tabFC2016$Stade%in%"JV"])/
    sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  
  
  
  
  #taux id incorrecte LIEC_X SA
  t2.EPI_LIECX.SA=sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("EPI_LIEC_X") & tabFC2016$Stade%in%"SA"])/
    sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  #taux id incorrecte LIEC_X AD
  t2.EPI_LIECX.AD=sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("EPI_LIEC_X") & tabFC2016$Stade%in%"AD"])/
    sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  #taux id incorrecte LIEC_X JV
  t2.EPI_LIECX.JV=sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("EPI_LIEC_X") & tabFC2016$Stade%in%"JV"])/
    sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  
  t2.EPI_EPI.AD+t2.EPI_EPI.SA+t2.EPI_EPI.JV+
    t2.EPI_EPA.AD+t2.EPI_EPA.SA+t2.EPI_EPA.JV+
    t2.EPI_ANS.AD+t2.EPI_ANS.SA+t2.EPI_ANS.JV+
    t2.EPI_END.AD+t2.EPI_END.SA+t2.EPI_END.JV+
    t2.EPI_LIECX.AD+t2.EPI_LIECX.SA+t2.EPI_LIECX.JV==1
  
  
  }


#les EPA identifi?s : que sont ils r?ellement ?
{
  
  #taux id EPI SA
  t2.EPA_EPI.SA=sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("EPA_EPI") & tabFC2016$Stade%in%"SA"])/
    sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  #taux id EPI AD
  t2.EPA_EPI.AD=sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("EPA_EPI") & tabFC2016$Stade%in%"AD"])/
    sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  #taux id EPI JV
  t2.EPA_EPI.JV=sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("EPA_EPI") & tabFC2016$Stade%in%"JV"])/
    sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  
  
  
  #taux id EPA SA
  t2.EPA_EPA.SA=sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("EPA_EPA") & tabFC2016$Stade%in%"SA"])/
    sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  #taux id EPA AD
  t2.EPA_EPA.AD=sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("EPA_EPA") & tabFC2016$Stade%in%"AD"])/
    sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  #taux id EPA JV
  t2.EPA_EPA.JV=sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("EPA_EPA") & tabFC2016$Stade%in%"JV"])/
    sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  
  
  
  #taux id ANS SA
  t2.EPA_ANS.SA=sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("EPA_ANS") & tabFC2016$Stade%in%"SA"])/
    sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  #taux id ANS AD
  t2.EPA_ANS.AD=sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("EPA_ANS") & tabFC2016$Stade%in%"AD"])/
    sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  #taux id ANS JV
  t2.EPA_ANS.JV=sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("EPA_ANS") & tabFC2016$Stade%in%"JV"])/
    sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  
  
  
  #taux id END SA
  t2.EPA_END.SA=sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("EPA_END") & tabFC2016$Stade%in%"SA"])/
    sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  #taux id END AD
  t2.EPA_END.AD=sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("EPA_END") & tabFC2016$Stade%in%"AD"])/
    sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  #taux id END JV
  t2.EPA_END.JV=sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("EPA_END") & tabFC2016$Stade%in%"JV"])/
    sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  
  
  
  
  #taux id LIEC_X SA
  t2.EPA_LIECX.SA=sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("EPA_LIEC_X") & tabFC2016$Stade%in%"SA"])/
    sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  #taux id LIEC_X AD
  t2.EPA_LIECX.AD=sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("EPA_LIEC_X") & tabFC2016$Stade%in%"AD"])/
    sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  #taux id LIEC_X JV
  t2.EPA_LIECX.JV=sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("EPA_LIEC_X") & tabFC2016$Stade%in%"JV"])/
    sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  
  round(t2.EPA_EPI.AD+t2.EPA_EPI.SA+t2.EPA_EPI.JV+
          t2.EPA_EPA.AD+t2.EPA_EPA.SA+t2.EPA_EPA.JV+
          t2.EPA_ANS.AD+t2.EPA_ANS.SA+t2.EPA_ANS.JV+
          t2.EPA_END.AD+t2.EPA_END.SA+t2.EPA_END.JV+
          t2.EPA_LIECX.AD+t2.EPA_LIECX.SA+t2.EPA_LIECX.JV,5)==1
  
  
}


#les ANS identifi?s : que sont ils r?ellement ?
{
  
  #taux id EPI SA
  t2.ANS_EPI.SA=sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("ANS_EPI") & tabFC2016$Stade%in%"SA"])/
    sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  #taux id EPI AD
  t2.ANS_EPI.AD=sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("ANS_EPI") & tabFC2016$Stade%in%"AD"])/
    sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  #taux id EPI JV
  t2.ANS_EPI.JV=sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("ANS_EPI") & tabFC2016$Stade%in%"JV"])/
    sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  
  
  
  #taux id EPA SA
  t2.ANS_EPA.SA=sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("ANS_EPA") & tabFC2016$Stade%in%"SA"])/
    sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  #taux id EPA AD
  t2.ANS_EPA.AD=sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("ANS_EPA") & tabFC2016$Stade%in%"AD"])/
    sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  #taux id EPA JV
  t2.ANS_EPA.JV=sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("ANS_EPA") & tabFC2016$Stade%in%"JV"])/
    sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  
  
  
  #taux id ANS SA
  t2.ANS_ANS.SA=sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("ANS_ANS") & tabFC2016$Stade%in%"SA"])/
    sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  #taux id ANS AD
  t2.ANS_ANS.AD=sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("ANS_ANS") & tabFC2016$Stade%in%"AD"])/
    sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  #taux id ANS JV
  t2.ANS_ANS.JV=sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("ANS_ANS") & tabFC2016$Stade%in%"JV"])/
    sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  
  
  
  #taux id END SA
  t2.ANS_END.SA=sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("ANS_END") & tabFC2016$Stade%in%"SA"])/
    sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  #taux id END AD
  t2.ANS_END.AD=sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("ANS_END") & tabFC2016$Stade%in%"AD"])/
    sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  #taux id END JV
  t2.ANS_END.JV=sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("ANS_END") & tabFC2016$Stade%in%"JV"])/
    sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  
  
  
  
  #taux id LIEC_X SA
  t2.ANS_LIECX.SA=sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("ANS_LIEC_X") & tabFC2016$Stade%in%"SA"])/
    sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  #taux id LIEC_X AD
  t2.ANS_LIECX.AD=sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("ANS_LIEC_X") & tabFC2016$Stade%in%"AD"])/
    sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  #taux id LIEC_X JV
  t2.ANS_LIECX.JV=sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("ANS_LIEC_X") & tabFC2016$Stade%in%"JV"])/
    sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  
  round(t2.ANS_EPI.AD+t2.ANS_EPI.SA+t2.ANS_EPI.JV+
          t2.ANS_EPA.AD+t2.ANS_EPA.SA+t2.ANS_EPA.JV+
          t2.ANS_ANS.AD+t2.ANS_ANS.SA+t2.ANS_ANS.JV+
          t2.ANS_END.AD+t2.ANS_END.SA+t2.ANS_END.JV+
          t2.ANS_LIECX.AD+t2.ANS_LIECX.SA+t2.ANS_LIECX.JV,5)==1
  
  
}


#les END identifi?s : que sont ils r?ellement ?
{
  
  #taux id EPI SA
  t2.END_EPI.SA=sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("END_EPI") & tabFC2016$Stade%in%"SA"])/
    sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  #taux id EPI AD
  t2.END_EPI.AD=sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("END_EPI") & tabFC2016$Stade%in%"AD"])/
    sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  #taux id EPI JV
  t2.END_EPI.JV=sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("END_EPI") & tabFC2016$Stade%in%"JV"])/
    sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  
  
  
  #taux id EPA SA
  t2.END_EPA.SA=sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("END_EPA") & tabFC2016$Stade%in%"SA"])/
    sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  #taux id EPA AD
  t2.END_EPA.AD=sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("END_EPA") & tabFC2016$Stade%in%"AD"])/
    sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  #taux id EPA JV
  t2.END_EPA.JV=sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("END_EPA") & tabFC2016$Stade%in%"JV"])/
    sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  
  
  
  #taux id ANS SA
  t2.END_ANS.SA=sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("END_ANS") & tabFC2016$Stade%in%"SA"])/
    sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  #taux id ANS AD
  t2.END_ANS.AD=sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("END_ANS") & tabFC2016$Stade%in%"AD"])/
    sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  #taux id ANS JV
  t2.END_ANS.JV=sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("END_ANS") & tabFC2016$Stade%in%"JV"])/
    sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  
  
  
  #taux id END SA
  t2.END_END.SA=sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("END_END") & tabFC2016$Stade%in%"SA"])/
    sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  #taux id END AD
  t2.END_END.AD=sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("END_END") & tabFC2016$Stade%in%"AD"])/
    sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  #taux id END JV
  t2.END_END.JV=sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("END_END") & tabFC2016$Stade%in%"JV"])/
    sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  
  
  
  
  #taux id LIEC_X SA
  t2.END_LIECX.SA=sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("END_LIEC_X") & tabFC2016$Stade%in%"SA"])/
    sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  #taux id LIEC_X AD
  t2.END_LIECX.AD=sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("END_LIEC_X") & tabFC2016$Stade%in%"AD"])/
    sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  #taux id LIEC_X JV
  t2.END_LIECX.JV=sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("END_LIEC_X") & tabFC2016$Stade%in%"JV"])/
    sum(tabFC2016$Nbr_EW[tabFC2016$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  
  round(t2.END_EPI.AD+t2.END_EPI.SA+t2.END_EPI.JV+
          t2.END_EPA.AD+t2.END_EPA.SA+t2.END_EPA.JV+
          t2.END_ANS.AD+t2.END_ANS.SA+t2.END_ANS.JV+
          t2.END_END.AD+t2.END_END.SA+t2.END_END.JV+
          t2.END_LIECX.AD+t2.END_LIECX.SA+t2.END_LIECX.JV,5)==1
  
  
}



tabGraphFC2=data.frame(
  Taux=c(
    t2.EPI_EPI.AD,t2.EPI_EPI.SA,t2.EPI_EPI.JV,
    t2.EPI_EPA.AD,t2.EPI_EPA.SA,t2.EPI_EPA.JV,
    t2.EPI_ANS.AD,t2.EPI_ANS.SA,t2.EPI_ANS.JV,
    t2.EPI_END.AD,t2.EPI_END.SA,t2.EPI_END.JV,
    t2.EPI_LIECX.AD,t2.EPI_LIECX.SA,t2.EPI_LIECX.JV,
    
    t2.EPA_EPI.AD,t2.EPA_EPI.SA,t2.EPA_EPI.JV,
    t2.EPA_EPA.AD,t2.EPA_EPA.SA,t2.EPA_EPA.JV,
    t2.EPA_ANS.AD,t2.EPA_ANS.SA,t2.EPA_ANS.JV,
    t2.EPA_END.AD,t2.EPA_END.SA,t2.EPA_END.JV,
    t2.EPA_LIECX.AD,t2.EPA_LIECX.SA,t2.EPA_LIECX.JV,
    
    t2.ANS_EPI.AD,t2.ANS_EPI.SA,t2.ANS_EPI.JV,
    t2.ANS_EPA.AD,t2.ANS_EPA.SA,t2.ANS_EPA.JV,
    t2.ANS_ANS.AD,t2.ANS_ANS.SA,t2.ANS_ANS.JV,
    t2.ANS_END.AD,t2.ANS_END.SA,t2.ANS_END.JV,
    t2.ANS_LIECX.AD,t2.ANS_LIECX.SA,t2.ANS_LIECX.JV,
    
    t2.END_EPI.AD,t2.END_EPI.SA,t2.END_EPI.JV,
    t2.END_EPA.AD,t2.END_EPA.SA,t2.END_EPA.JV,
    t2.END_ANS.AD,t2.END_ANS.SA,t2.END_ANS.JV,
    t2.END_END.AD,t2.END_END.SA,t2.END_END.JV,
    t2.END_LIECX.AD,t2.END_LIECX.SA,t2.END_LIECX.JV),
  
  Classe=c(
    "t2.EPI_EPI.AD","t2.EPI_EPI.SA","t2.EPI_EPI.JV",
    "t2.EPI_EPA.AD","t2.EPI_EPA.SA","t2.EPI_EPA.JV",
    "t2.EPI_ANS.AD","t2.EPI_ANS.SA","t2.EPI_ANS.JV",
    "t2.EPI_END.AD","t2.EPI_END.SA","t2.EPI_END.JV",
    "t2.EPI_LIECX.AD","t2.EPI_LIECX.SA","t2.EPI_LIECX.JV",
    
    "t2.EPA_EPI.AD","t2.EPA_EPI.SA","t2.EPA_EPI.JV",
    "t2.EPA_EPA.AD","t2.EPA_EPA.SA","t2.EPA_EPA.JV",
    "t2.EPA_ANS.AD","t2.EPA_ANS.SA","t2.EPA_ANS.JV",
    "t2.EPA_END.AD","t2.EPA_END.SA","t2.EPA_END.JV",
    "t2.EPA_LIECX.AD","t2.EPA_LIECX.SA","t2.EPA_LIECX.JV",
    
    "t2.ANS_EPI.AD","t2.ANS_EPI.SA","t2.ANS_EPI.JV",
    "t2.ANS_EPA.AD","t2.ANS_EPA.SA","t2.ANS_EPA.JV",
    "t2.ANS_ANS.AD","t2.ANS_ANS.SA","t2.ANS_ANS.JV",
    "t2.ANS_END.AD","t2.ANS_END.SA","t2.ANS_END.JV",
    "t2.ANS_LIECX.AD","t2.ANS_LIECX.SA","t2.ANS_LIECX.JV",
    
    "t2.END_EPI.AD","t2.END_EPI.SA","t2.END_EPI.JV",
    "t2.END_EPA.AD","t2.END_EPA.SA","t2.END_EPA.JV",
    "t2.END_ANS.AD","t2.END_ANS.SA","t2.END_ANS.JV",
    "t2.END_END.AD","t2.END_END.SA","t2.END_END.JV",
    "t2.END_LIECX.AD","t2.END_LIECX.SA","t2.END_LIECX.JV")
  
)
tabGraphFC2$LIEC=factor(substr(tabGraphFC2$Classe,8,10),levels=c("EPI","EPA","ANS","END","LIECX"))
tabGraphFC2$Stade=substr(tabGraphFC2$Classe,12,13)
tabGraphFC2$LIECid=factor(substr(tabGraphFC2$Classe,4,6),levels=c("EPI","EPA","ANS","END","LIECX"))
tabGraphFC2$LIEC_Stade=paste(tabGraphFC2$LIEC,tabGraphFC2$Stade,sep="_")

tabGraphFC2$LIEC_Stade=factor(tabGraphFC2$LIEC_Stade,levels=c(
  "EPI_AD","EPI_SA","EPI_JV","EPA_AD","EPA_SA",
  "EPA_JV","ANS_AD","ANS_SA","ANS_JV","END_AD",
  "END_SA","END_JV","LIECX_AD","LIECX_SA","LIECX_JV"))

ggplot(tabGraphFC2, aes(fill=LIEC_Stade, y=Taux, x=LIECid)) + 
  geom_bar(position='stack', stat='identity')+
  theme_classic()+
  scale_fill_manual(values = c(
    rgb(0,0,153,maxColorValue = 255),    rgb(0,0,153,maxColorValue = 255,alpha=150),    rgb(0,0,153,maxColorValue = 255,alpha=100),
    rgb(255,255,0,maxColorValue = 255),  rgb(255,255,0,maxColorValue = 255,alpha=100),  rgb(255,255,0,maxColorValue = 255,alpha=25),
    rgb(0,153,0,maxColorValue = 255),    rgb(0,153,0,maxColorValue = 255,alpha=150),    rgb(0,153,0,maxColorValue = 255,alpha=100),
    rgb(255,102,153,maxColorValue = 255),rgb(255,102,153,maxColorValue = 255,alpha=150),rgb(255,102,153,maxColorValue = 255,alpha=100),
    rgb(123,123,123,maxColorValue = 255),rgb(123,123,123,maxColorValue = 255,alpha=150),rgb(123,123,123,maxColorValue = 255,alpha=100)
  ))

#zoom sur les erreurs
tabGraphFC2$Erreur=NA
tabGraphFC2$Erreur[tabGraphFC2$LIEC==tabGraphFC2$LIECid]="Correcte"
tabGraphFC2$Erreur[tabGraphFC2$LIEC!=tabGraphFC2$LIECid]="Incorrecte"

ggplot(tabGraphFC2[tabGraphFC2$Erreur=="Incorrecte",], aes(fill=LIEC_Stade, y=Taux, x=LIECid)) + 
  geom_bar(position='fill', stat='identity')+
  theme_classic()+
  scale_fill_manual(values = c(
    rgb(0,0,153,maxColorValue = 255),    rgb(0,0,153,maxColorValue = 255,alpha=150),    rgb(0,0,153,maxColorValue = 255,alpha=100),
    rgb(255,255,0,maxColorValue = 255),  rgb(255,255,0,maxColorValue = 255,alpha=100),  rgb(255,255,0,maxColorValue = 255,alpha=25),
    rgb(0,153,0,maxColorValue = 255),    rgb(0,153,0,maxColorValue = 255,alpha=150),    rgb(0,153,0,maxColorValue = 255,alpha=100),
    rgb(255,102,153,maxColorValue = 255),rgb(255,102,153,maxColorValue = 255,alpha=150),rgb(255,102,153,maxColorValue = 255,alpha=100),
    rgb(123,123,123,maxColorValue = 255),rgb(123,123,123,maxColorValue = 255,alpha=150),rgb(123,123,123,maxColorValue = 255,alpha=100)
  ))




##### BO 2016 zoom sur les points critiques - Parmis ce que les gens ont class?s en EPI EPA ANS END : qu'est ce qu'on a r?ellement ? #####


tabBO=Tab2[Tab2$Region%in%"BO",]
nrow(tabBO)==3426
dataID <- unique(sbt_m[,c("cle",colID,"Region")])
nrow(dataID)==1109
dataID=dataID[dataID$Region%in%"BO",]
nrow(dataID)==144



tabBO2016=tabBO[tabBO$Annee%in%2016,]
unique(tabBO2016$Stade)#"JV" "SA" "AD" => ok
any(is.na(tabBO2016$Stade))==T
nrow(tabBO2016)==488

#taux d'id correcte global
sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPI_EPI","EPA_EPA","ANS_ANS","END_END")])/
  sum(tabBO2016$Nbr_EW)#0.8867060  => ok
#vectReussiteBO

tabBO2016=tabBO2016[!is.na(tabBO2016$Stade),]
any(is.na(tabBO2016$Stade))==F
nrow(tabBO2016)==488-3

#les EPI identifi?s : que sont ils r?ellement ?
{
  
  #taux id correcte EPI SA
  t2.EPI_EPI.SA=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPI_EPI") & tabBO2016$Stade%in%"SA"])/
    sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  #taux id correcte EPI AD
  t2.EPI_EPI.AD=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPI_EPI") & tabBO2016$Stade%in%"AD"])/
    sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  #taux id correcte EPI JV
  t2.EPI_EPI.JV=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPI_EPI") & tabBO2016$Stade%in%"JV"])/
    sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  
  
  
  #taux id incorrecte EPA SA
  t2.EPI_EPA.SA=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPI_EPA") & tabBO2016$Stade%in%"SA"])/
    sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  #taux id incorrecte EPA AD
  t2.EPI_EPA.AD=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPI_EPA") & tabBO2016$Stade%in%"AD"])/
    sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  #taux id incorrecte EPA JV
  t2.EPI_EPA.JV=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPI_EPA") & tabBO2016$Stade%in%"JV"])/
    sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  
  
  
  #taux id incorrecte ANS SA
  t2.EPI_ANS.SA=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPI_ANS") & tabBO2016$Stade%in%"SA"])/
    sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  #taux id incorrecte ANS AD
  t2.EPI_ANS.AD=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPI_ANS") & tabBO2016$Stade%in%"AD"])/
    sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  #taux id incorrecte ANS JV
  t2.EPI_ANS.JV=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPI_ANS") & tabBO2016$Stade%in%"JV"])/
    sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  
  
  
  #taux id incorrecte END SA
  t2.EPI_END.SA=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPI_END") & tabBO2016$Stade%in%"SA"])/
    sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  #taux id incorrecte END AD
  t2.EPI_END.AD=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPI_END") & tabBO2016$Stade%in%"AD"])/
    sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  #taux id incorrecte END JV
  t2.EPI_END.JV=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPI_END") & tabBO2016$Stade%in%"JV"])/
    sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  
  
  
  
  #taux id incorrecte LIEC_X SA
  t2.EPI_LIECX.SA=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPI_LIEC_X") & tabBO2016$Stade%in%"SA"])/
    sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  #taux id incorrecte LIEC_X AD
  t2.EPI_LIECX.AD=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPI_LIEC_X") & tabBO2016$Stade%in%"AD"])/
    sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  #taux id incorrecte LIEC_X JV
  t2.EPI_LIECX.JV=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPI_LIEC_X") & tabBO2016$Stade%in%"JV"])/
    sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  
  t2.EPI_EPI.AD+t2.EPI_EPI.SA+t2.EPI_EPI.JV+
    t2.EPI_EPA.AD+t2.EPI_EPA.SA+t2.EPI_EPA.JV+
    t2.EPI_ANS.AD+t2.EPI_ANS.SA+t2.EPI_ANS.JV+
    t2.EPI_END.AD+t2.EPI_END.SA+t2.EPI_END.JV+
    t2.EPI_LIECX.AD+t2.EPI_LIECX.SA+t2.EPI_LIECX.JV==1
  
  
  }


#les EPA identifi?s : que sont ils r?ellement ?
{
  
  #taux id EPI SA
  t2.EPA_EPI.SA=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPA_EPI") & tabBO2016$Stade%in%"SA"])/
    sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  #taux id EPI AD
  t2.EPA_EPI.AD=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPA_EPI") & tabBO2016$Stade%in%"AD"])/
    sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  #taux id EPI JV
  t2.EPA_EPI.JV=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPA_EPI") & tabBO2016$Stade%in%"JV"])/
    sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  
  
  
  #taux id EPA SA
  t2.EPA_EPA.SA=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPA_EPA") & tabBO2016$Stade%in%"SA"])/
    sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  #taux id EPA AD
  t2.EPA_EPA.AD=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPA_EPA") & tabBO2016$Stade%in%"AD"])/
    sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  #taux id EPA JV
  t2.EPA_EPA.JV=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPA_EPA") & tabBO2016$Stade%in%"JV"])/
    sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  
  
  
  #taux id ANS SA
  t2.EPA_ANS.SA=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPA_ANS") & tabBO2016$Stade%in%"SA"])/
    sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  #taux id ANS AD
  t2.EPA_ANS.AD=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPA_ANS") & tabBO2016$Stade%in%"AD"])/
    sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  #taux id ANS JV
  t2.EPA_ANS.JV=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPA_ANS") & tabBO2016$Stade%in%"JV"])/
    sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  
  
  
  #taux id END SA
  t2.EPA_END.SA=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPA_END") & tabBO2016$Stade%in%"SA"])/
    sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  #taux id END AD
  t2.EPA_END.AD=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPA_END") & tabBO2016$Stade%in%"AD"])/
    sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  #taux id END JV
  t2.EPA_END.JV=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPA_END") & tabBO2016$Stade%in%"JV"])/
    sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  
  
  
  
  #taux id LIEC_X SA
  t2.EPA_LIECX.SA=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPA_LIEC_X") & tabBO2016$Stade%in%"SA"])/
    sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  #taux id LIEC_X AD
  t2.EPA_LIECX.AD=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPA_LIEC_X") & tabBO2016$Stade%in%"AD"])/
    sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  #taux id LIEC_X JV
  t2.EPA_LIECX.JV=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPA_LIEC_X") & tabBO2016$Stade%in%"JV"])/
    sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  
  round(t2.EPA_EPI.AD+t2.EPA_EPI.SA+t2.EPA_EPI.JV+
          t2.EPA_EPA.AD+t2.EPA_EPA.SA+t2.EPA_EPA.JV+
          t2.EPA_ANS.AD+t2.EPA_ANS.SA+t2.EPA_ANS.JV+
          t2.EPA_END.AD+t2.EPA_END.SA+t2.EPA_END.JV+
          t2.EPA_LIECX.AD+t2.EPA_LIECX.SA+t2.EPA_LIECX.JV,5)==1
  
  
}


#les ANS identifi?s : que sont ils r?ellement ?
{
  
  #taux id EPI SA
  t2.ANS_EPI.SA=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("ANS_EPI") & tabBO2016$Stade%in%"SA"])/
    sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  #taux id EPI AD
  t2.ANS_EPI.AD=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("ANS_EPI") & tabBO2016$Stade%in%"AD"])/
    sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  #taux id EPI JV
  t2.ANS_EPI.JV=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("ANS_EPI") & tabBO2016$Stade%in%"JV"])/
    sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  
  
  
  #taux id EPA SA
  t2.ANS_EPA.SA=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("ANS_EPA") & tabBO2016$Stade%in%"SA"])/
    sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  #taux id EPA AD
  t2.ANS_EPA.AD=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("ANS_EPA") & tabBO2016$Stade%in%"AD"])/
    sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  #taux id EPA JV
  t2.ANS_EPA.JV=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("ANS_EPA") & tabBO2016$Stade%in%"JV"])/
    sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  
  
  
  #taux id ANS SA
  t2.ANS_ANS.SA=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("ANS_ANS") & tabBO2016$Stade%in%"SA"])/
    sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  #taux id ANS AD
  t2.ANS_ANS.AD=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("ANS_ANS") & tabBO2016$Stade%in%"AD"])/
    sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  #taux id ANS JV
  t2.ANS_ANS.JV=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("ANS_ANS") & tabBO2016$Stade%in%"JV"])/
    sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  
  
  
  #taux id END SA
  t2.ANS_END.SA=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("ANS_END") & tabBO2016$Stade%in%"SA"])/
    sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  #taux id END AD
  t2.ANS_END.AD=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("ANS_END") & tabBO2016$Stade%in%"AD"])/
    sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  #taux id END JV
  t2.ANS_END.JV=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("ANS_END") & tabBO2016$Stade%in%"JV"])/
    sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  
  
  
  
  #taux id LIEC_X SA
  t2.ANS_LIECX.SA=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("ANS_LIEC_X") & tabBO2016$Stade%in%"SA"])/
    sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  #taux id LIEC_X AD
  t2.ANS_LIECX.AD=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("ANS_LIEC_X") & tabBO2016$Stade%in%"AD"])/
    sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  #taux id LIEC_X JV
  t2.ANS_LIECX.JV=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("ANS_LIEC_X") & tabBO2016$Stade%in%"JV"])/
    sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  
  round(t2.ANS_EPI.AD+t2.ANS_EPI.SA+t2.ANS_EPI.JV+
          t2.ANS_EPA.AD+t2.ANS_EPA.SA+t2.ANS_EPA.JV+
          t2.ANS_ANS.AD+t2.ANS_ANS.SA+t2.ANS_ANS.JV+
          t2.ANS_END.AD+t2.ANS_END.SA+t2.ANS_END.JV+
          t2.ANS_LIECX.AD+t2.ANS_LIECX.SA+t2.ANS_LIECX.JV,5)==1
  
  
}


#les END identifi?s : que sont ils r?ellement ?
{
  
  #taux id EPI SA
  t2.END_EPI.SA=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("END_EPI") & tabBO2016$Stade%in%"SA"])/
    sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  #taux id EPI AD
  t2.END_EPI.AD=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("END_EPI") & tabBO2016$Stade%in%"AD"])/
    sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  #taux id EPI JV
  t2.END_EPI.JV=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("END_EPI") & tabBO2016$Stade%in%"JV"])/
    sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  
  
  
  #taux id EPA SA
  t2.END_EPA.SA=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("END_EPA") & tabBO2016$Stade%in%"SA"])/
    sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  #taux id EPA AD
  t2.END_EPA.AD=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("END_EPA") & tabBO2016$Stade%in%"AD"])/
    sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  #taux id EPA JV
  t2.END_EPA.JV=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("END_EPA") & tabBO2016$Stade%in%"JV"])/
    sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  
  
  
  #taux id ANS SA
  t2.END_ANS.SA=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("END_ANS") & tabBO2016$Stade%in%"SA"])/
    sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  #taux id ANS AD
  t2.END_ANS.AD=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("END_ANS") & tabBO2016$Stade%in%"AD"])/
    sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  #taux id ANS JV
  t2.END_ANS.JV=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("END_ANS") & tabBO2016$Stade%in%"JV"])/
    sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  
  
  
  #taux id END SA
  t2.END_END.SA=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("END_END") & tabBO2016$Stade%in%"SA"])/
    sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  #taux id END AD
  t2.END_END.AD=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("END_END") & tabBO2016$Stade%in%"AD"])/
    sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  #taux id END JV
  t2.END_END.JV=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("END_END") & tabBO2016$Stade%in%"JV"])/
    sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  
  
  
  
  #taux id LIEC_X SA
  t2.END_LIECX.SA=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("END_LIEC_X") & tabBO2016$Stade%in%"SA"])/
    sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  #taux id LIEC_X AD
  t2.END_LIECX.AD=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("END_LIEC_X") & tabBO2016$Stade%in%"AD"])/
    sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  #taux id LIEC_X JV
  t2.END_LIECX.JV=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("END_LIEC_X") & tabBO2016$Stade%in%"JV"])/
    sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  
  round(t2.END_EPI.AD+t2.END_EPI.SA+t2.END_EPI.JV+
          t2.END_EPA.AD+t2.END_EPA.SA+t2.END_EPA.JV+
          t2.END_ANS.AD+t2.END_ANS.SA+t2.END_ANS.JV+
          t2.END_END.AD+t2.END_END.SA+t2.END_END.JV+
          t2.END_LIECX.AD+t2.END_LIECX.SA+t2.END_LIECX.JV,5)==1
  
  
}


#les LIECX identifi?s : que sont ils r?ellement ?
{
  
  #taux id EPI SA
  t2.LIECX_EPI.SA=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("LIEC_X_EPI") & tabBO2016$Stade%in%"SA"])/
    sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  #taux id EPI AD
  t2.LIECX_EPI.AD=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("LIEC_X_EPI") & tabBO2016$Stade%in%"AD"])/
    sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  #taux id EPI JV
  t2.LIECX_EPI.JV=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("LIEC_X_EPI") & tabBO2016$Stade%in%"JV"])/
    sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  
  
  
  #taux id EPA SA
  t2.LIECX_EPA.SA=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("LIEC_X_EPA") & tabBO2016$Stade%in%"SA"])/
    sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  #taux id EPA AD
  t2.LIECX_EPA.AD=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("LIEC_X_EPA") & tabBO2016$Stade%in%"AD"])/
    sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  #taux id EPA JV
  t2.LIECX_EPA.JV=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("LIEC_X_EPA") & tabBO2016$Stade%in%"JV"])/
    sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  
  
  
  #taux id ANS SA
  t2.LIECX_ANS.SA=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("LIEC_X_ANS") & tabBO2016$Stade%in%"SA"])/
    sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  #taux id ANS AD
  t2.LIECX_ANS.AD=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("LIEC_X_ANS") & tabBO2016$Stade%in%"AD"])/
    sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  #taux id ANS JV
  t2.LIECX_ANS.JV=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("LIEC_X_ANS") & tabBO2016$Stade%in%"JV"])/
    sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  
  
  
  #taux id END SA
  t2.LIECX_END.SA=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("LIEC_X_END") & tabBO2016$Stade%in%"SA"])/
    sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  #taux id END AD
  t2.LIECX_END.AD=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("LIEC_X_END") & tabBO2016$Stade%in%"AD"])/
    sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  #taux id END JV
  t2.LIECX_END.JV=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("LIEC_X_END") & tabBO2016$Stade%in%"JV"])/
    sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  
  
  
  
  #taux id LIEC_X SA
  t2.LIECX_LIECX.SA=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("LIEC_X_LIEC_X") & tabBO2016$Stade%in%"SA"])/
    sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  #taux id LIEC_X AD
  t2.LIECX_LIECX.AD=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("LIEC_X_LIEC_X") & tabBO2016$Stade%in%"AD"])/
    sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  #taux id LIEC_X JV
  t2.LIECX_LIECX.JV=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("LIEC_X_LIEC_X") & tabBO2016$Stade%in%"JV"])/
    sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  
  round(t2.LIECX_EPI.AD+t2.LIECX_EPI.SA+t2.LIECX_EPI.JV+
          t2.LIECX_EPA.AD+t2.LIECX_EPA.SA+t2.LIECX_EPA.JV+
          t2.LIECX_ANS.AD+t2.LIECX_ANS.SA+t2.LIECX_ANS.JV+
          t2.LIECX_END.AD+t2.LIECX_END.SA+t2.LIECX_END.JV+
          t2.LIECX_LIECX.AD+t2.LIECX_LIECX.SA+t2.LIECX_LIECX.JV,5)==1
  
  
}


tabGraphBO16_2=data.frame(
  Taux=c(
    t2.EPI_EPI.AD,t2.EPI_EPI.SA,t2.EPI_EPI.JV,
    t2.EPI_EPA.AD,t2.EPI_EPA.SA,t2.EPI_EPA.JV,
    t2.EPI_ANS.AD,t2.EPI_ANS.SA,t2.EPI_ANS.JV,
    t2.EPI_END.AD,t2.EPI_END.SA,t2.EPI_END.JV,
    t2.EPI_LIECX.AD,t2.EPI_LIECX.SA,t2.EPI_LIECX.JV,
    
    t2.EPA_EPI.AD,t2.EPA_EPI.SA,t2.EPA_EPI.JV,
    t2.EPA_EPA.AD,t2.EPA_EPA.SA,t2.EPA_EPA.JV,
    t2.EPA_ANS.AD,t2.EPA_ANS.SA,t2.EPA_ANS.JV,
    t2.EPA_END.AD,t2.EPA_END.SA,t2.EPA_END.JV,
    t2.EPA_LIECX.AD,t2.EPA_LIECX.SA,t2.EPA_LIECX.JV,
    
    t2.ANS_EPI.AD,t2.ANS_EPI.SA,t2.ANS_EPI.JV,
    t2.ANS_EPA.AD,t2.ANS_EPA.SA,t2.ANS_EPA.JV,
    t2.ANS_ANS.AD,t2.ANS_ANS.SA,t2.ANS_ANS.JV,
    t2.ANS_END.AD,t2.ANS_END.SA,t2.ANS_END.JV,
    t2.ANS_LIECX.AD,t2.ANS_LIECX.SA,t2.ANS_LIECX.JV,
    
    t2.END_EPI.AD,t2.END_EPI.SA,t2.END_EPI.JV,
    t2.END_EPA.AD,t2.END_EPA.SA,t2.END_EPA.JV,
    t2.END_ANS.AD,t2.END_ANS.SA,t2.END_ANS.JV,
    t2.END_END.AD,t2.END_END.SA,t2.END_END.JV,
    t2.END_LIECX.AD,t2.END_LIECX.SA,t2.END_LIECX.JV,
    
    t2.LIECX_EPI.AD,t2.LIECX_EPI.SA,t2.LIECX_EPI.JV,
    t2.LIECX_EPA.AD,t2.LIECX_EPA.SA,t2.LIECX_EPA.JV,
    t2.LIECX_ANS.AD,t2.LIECX_ANS.SA,t2.LIECX_ANS.JV,
    t2.LIECX_END.AD,t2.LIECX_END.SA,t2.LIECX_END.JV,
    t2.LIECX_LIECX.AD,t2.LIECX_LIECX.SA,t2.LIECX_LIECX.JV),
  
  Classe=c(
    "t2.EPI_EPI.AD","t2.EPI_EPI.SA","t2.EPI_EPI.JV",
    "t2.EPI_EPA.AD","t2.EPI_EPA.SA","t2.EPI_EPA.JV",
    "t2.EPI_ANS.AD","t2.EPI_ANS.SA","t2.EPI_ANS.JV",
    "t2.EPI_END.AD","t2.EPI_END.SA","t2.EPI_END.JV",
    "t2.EPI_LIECX.AD","t2.EPI_LIECX.SA","t2.EPI_LIECX.JV",
    
    "t2.EPA_EPI.AD","t2.EPA_EPI.SA","t2.EPA_EPI.JV",
    "t2.EPA_EPA.AD","t2.EPA_EPA.SA","t2.EPA_EPA.JV",
    "t2.EPA_ANS.AD","t2.EPA_ANS.SA","t2.EPA_ANS.JV",
    "t2.EPA_END.AD","t2.EPA_END.SA","t2.EPA_END.JV",
    "t2.EPA_LIECX.AD","t2.EPA_LIECX.SA","t2.EPA_LIECX.JV",
    
    "t2.ANS_EPI.AD","t2.ANS_EPI.SA","t2.ANS_EPI.JV",
    "t2.ANS_EPA.AD","t2.ANS_EPA.SA","t2.ANS_EPA.JV",
    "t2.ANS_ANS.AD","t2.ANS_ANS.SA","t2.ANS_ANS.JV",
    "t2.ANS_END.AD","t2.ANS_END.SA","t2.ANS_END.JV",
    "t2.ANS_LIECX.AD","t2.ANS_LIECX.SA","t2.ANS_LIECX.JV",
    
    "t2.END_EPI.AD","t2.END_EPI.SA","t2.END_EPI.JV",
    "t2.END_EPA.AD","t2.END_EPA.SA","t2.END_EPA.JV",
    "t2.END_ANS.AD","t2.END_ANS.SA","t2.END_ANS.JV",
    "t2.END_END.AD","t2.END_END.SA","t2.END_END.JV",
    "t2.END_LIECX.AD","t2.END_LIECX.SA","t2.END_LIECX.JV",
    
    "t2.LIECX_EPI.AD","t2.LIECX_EPI.SA","t2.LIECX_EPI.JV",
    "t2.LIECX_EPA.AD","t2.LIECX_EPA.SA","t2.LIECX_EPA.JV",
    "t2.LIECX_ANS.AD","t2.LIECX_ANS.SA","t2.LIECX_ANS.JV",
    "t2.LIECX_END.AD","t2.LIECX_END.SA","t2.LIECX_END.JV",
    "t2.LIECX_LIECX.AD","t2.LIECX_LIECX.SA","t2.LIECX_LIECX.JV")
  
)
tabGraphBO16_2$LIEC=factor(substr(tabGraphBO16_2$Classe,8,10),levels=c("EPI","EPA","ANS","END","LIECX"))
tabGraphBO16_2$Stade=substr(tabGraphBO16_2$Classe,12,13)
tabGraphBO16_2$LIECid=factor(substr(tabGraphBO16_2$Classe,4,6),levels=c("EPI","EPA","ANS","END","LIECX"))
tabGraphBO16_2$LIEC_Stade=paste(tabGraphBO16_2$LIEC,tabGraphBO16_2$Stade,sep="_")

tabGraphBO16_2$LIEC_Stade=factor(tabGraphBO16_2$LIEC_Stade,levels=c(
  "EPI_AD","EPI_SA","EPI_JV","EPA_AD","EPA_SA",
  "EPA_JV","ANS_AD","ANS_SA","ANS_JV","END_AD",
  "END_SA","END_JV","LIECX_AD","LIECX_SA","LIECX_JV"))

ggplot(tabGraphBO16_2, aes(fill=LIEC_Stade, y=Taux, x=LIECid)) + 
  geom_bar(position='stack', stat='identity')+
  theme_classic()+
  scale_fill_manual(values = c(
    rgb(0,0,153,maxColorValue = 255),    rgb(0,0,153,maxColorValue = 255,alpha=150),    rgb(0,0,153,maxColorValue = 255,alpha=100),
    rgb(255,255,0,maxColorValue = 255),  rgb(255,255,0,maxColorValue = 255,alpha=100),  rgb(255,255,0,maxColorValue = 255,alpha=25),
    rgb(0,153,0,maxColorValue = 255),    rgb(0,153,0,maxColorValue = 255,alpha=150),    rgb(0,153,0,maxColorValue = 255,alpha=100),
    rgb(255,102,153,maxColorValue = 255),rgb(255,102,153,maxColorValue = 255,alpha=150),rgb(255,102,153,maxColorValue = 255,alpha=100),
    rgb(123,123,123,maxColorValue = 255),rgb(123,123,123,maxColorValue = 255,alpha=150),rgb(123,123,123,maxColorValue = 255,alpha=100)
  ))

ggplot(tabGraphBO16_2, aes(fill=LIEC_Stade, y=Taux, x=LIECid)) + 
  geom_bar(position='stack', stat='identity')+
  theme_classic()+
  scale_fill_manual(values = c(
    rgb(0,0,153,maxColorValue = 255),    rgb(0,0,153,maxColorValue = 255),    rgb(0,0,153,maxColorValue = 255,alpha=150),
    rgb(255,255,0,maxColorValue = 255),  rgb(255,255,0,maxColorValue = 255),  rgb(255,255,0,maxColorValue = 255,alpha=100),
    rgb(0,153,0,maxColorValue = 255),    rgb(0,153,0,maxColorValue = 255),    rgb(0,153,0,maxColorValue = 255,alpha=150),
    rgb(255,102,153,maxColorValue = 255),rgb(255,102,153,maxColorValue = 255),rgb(255,102,153,maxColorValue = 255,alpha=150),
    rgb(123,123,123,maxColorValue = 255),rgb(123,123,123,maxColorValue = 255),rgb(123,123,123,maxColorValue = 255,alpha=150)
  ))

#zoom sur les erreurs
tabGraphBO16_2$Erreur=NA
tabGraphBO16_2$Erreur[tabGraphBO16_2$LIEC==tabGraphBO16_2$LIECid]="Correcte"
tabGraphBO16_2$Erreur[tabGraphBO16_2$LIEC!=tabGraphBO16_2$LIECid]="Incorrecte"

ggplot(tabGraphBO16_2[tabGraphBO16_2$Erreur=="Incorrecte",], aes(fill=LIEC_Stade, y=Taux, x=LIECid)) + 
  geom_bar(position='fill', stat='identity')+
  theme_classic()+
  scale_fill_manual(values = c(
    rgb(0,0,153,maxColorValue = 255),    rgb(0,0,153,maxColorValue = 255,alpha=150),    rgb(0,0,153,maxColorValue = 255,alpha=100),
    rgb(255,255,0,maxColorValue = 255),  rgb(255,255,0,maxColorValue = 255,alpha=100),  rgb(255,255,0,maxColorValue = 255,alpha=25),
    rgb(0,153,0,maxColorValue = 255),    rgb(0,153,0,maxColorValue = 255,alpha=150),    rgb(0,153,0,maxColorValue = 255,alpha=100),
    rgb(255,102,153,maxColorValue = 255),rgb(255,102,153,maxColorValue = 255,alpha=150),rgb(255,102,153,maxColorValue = 255,alpha=100),
    rgb(123,123,123,maxColorValue = 255),rgb(123,123,123,maxColorValue = 255,alpha=150),rgb(123,123,123,maxColorValue = 255,alpha=100)
  ))




##### BO 2017 zoom sur les points critiques - Parmis ce que les gens ont class?s en EPI EPA ANS END : qu'est ce qu'on a r?ellement ? #####


tabBO=Tab2[Tab2$Region%in%"BO",]
nrow(tabBO)==3426
dataID <- unique(sbt_m[,c("cle",colID,"Region")])
nrow(dataID)==1109
dataID=dataID[dataID$Region%in%"BO",]
nrow(dataID)==144



tabBO2017=tabBO[tabBO$Annee%in%2017,]
unique(tabBO2017$Stade)#"JV" "SA" "AD" => ok
any(is.na(tabBO2017$Stade))==F
nrow(tabBO2017)==471

#taux d'id correcte global
sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPI_EPI","EPA_EPA","ANS_ANS","END_END")])/
  sum(tabBO2017$Nbr_EW)#0.6357143   => ok
#vectReussiteBO

#les EPI identifi?s : que sont ils r?ellement ?
{
  
  #taux id correcte EPI SA
  t2.EPI_EPI.SA=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPI_EPI") & tabBO2017$Stade%in%"SA"])/
    sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  #taux id correcte EPI AD
  t2.EPI_EPI.AD=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPI_EPI") & tabBO2017$Stade%in%"AD"])/
    sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  #taux id correcte EPI JV
  t2.EPI_EPI.JV=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPI_EPI") & tabBO2017$Stade%in%"JV"])/
    sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  
  
  
  #taux id incorrecte EPA SA
  t2.EPI_EPA.SA=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPI_EPA") & tabBO2017$Stade%in%"SA"])/
    sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  #taux id incorrecte EPA AD
  t2.EPI_EPA.AD=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPI_EPA") & tabBO2017$Stade%in%"AD"])/
    sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  #taux id incorrecte EPA JV
  t2.EPI_EPA.JV=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPI_EPA") & tabBO2017$Stade%in%"JV"])/
    sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  
  
  
  #taux id incorrecte ANS SA
  t2.EPI_ANS.SA=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPI_ANS") & tabBO2017$Stade%in%"SA"])/
    sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  #taux id incorrecte ANS AD
  t2.EPI_ANS.AD=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPI_ANS") & tabBO2017$Stade%in%"AD"])/
    sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  #taux id incorrecte ANS JV
  t2.EPI_ANS.JV=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPI_ANS") & tabBO2017$Stade%in%"JV"])/
    sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  
  
  
  #taux id incorrecte END SA
  t2.EPI_END.SA=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPI_END") & tabBO2017$Stade%in%"SA"])/
    sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  #taux id incorrecte END AD
  t2.EPI_END.AD=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPI_END") & tabBO2017$Stade%in%"AD"])/
    sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  #taux id incorrecte END JV
  t2.EPI_END.JV=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPI_END") & tabBO2017$Stade%in%"JV"])/
    sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  
  
  
  
  #taux id incorrecte LIEC_X SA
  t2.EPI_LIECX.SA=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPI_LIEC_X") & tabBO2017$Stade%in%"SA"])/
    sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  #taux id incorrecte LIEC_X AD
  t2.EPI_LIECX.AD=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPI_LIEC_X") & tabBO2017$Stade%in%"AD"])/
    sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  #taux id incorrecte LIEC_X JV
  t2.EPI_LIECX.JV=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPI_LIEC_X") & tabBO2017$Stade%in%"JV"])/
    sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  
  t2.EPI_EPI.AD+t2.EPI_EPI.SA+t2.EPI_EPI.JV+
    t2.EPI_EPA.AD+t2.EPI_EPA.SA+t2.EPI_EPA.JV+
    t2.EPI_ANS.AD+t2.EPI_ANS.SA+t2.EPI_ANS.JV+
    t2.EPI_END.AD+t2.EPI_END.SA+t2.EPI_END.JV+
    t2.EPI_LIECX.AD+t2.EPI_LIECX.SA+t2.EPI_LIECX.JV==1
  
  
}


#les EPA identifi?s : que sont ils r?ellement ?
{
  
  #taux id EPI SA
  t2.EPA_EPI.SA=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPA_EPI") & tabBO2017$Stade%in%"SA"])/
    sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  #taux id EPI AD
  t2.EPA_EPI.AD=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPA_EPI") & tabBO2017$Stade%in%"AD"])/
    sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  #taux id EPI JV
  t2.EPA_EPI.JV=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPA_EPI") & tabBO2017$Stade%in%"JV"])/
    sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  
  
  
  #taux id EPA SA
  t2.EPA_EPA.SA=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPA_EPA") & tabBO2017$Stade%in%"SA"])/
    sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  #taux id EPA AD
  t2.EPA_EPA.AD=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPA_EPA") & tabBO2017$Stade%in%"AD"])/
    sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  #taux id EPA JV
  t2.EPA_EPA.JV=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPA_EPA") & tabBO2017$Stade%in%"JV"])/
    sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  
  
  
  #taux id ANS SA
  t2.EPA_ANS.SA=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPA_ANS") & tabBO2017$Stade%in%"SA"])/
    sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  #taux id ANS AD
  t2.EPA_ANS.AD=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPA_ANS") & tabBO2017$Stade%in%"AD"])/
    sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  #taux id ANS JV
  t2.EPA_ANS.JV=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPA_ANS") & tabBO2017$Stade%in%"JV"])/
    sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  
  
  
  #taux id END SA
  t2.EPA_END.SA=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPA_END") & tabBO2017$Stade%in%"SA"])/
    sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  #taux id END AD
  t2.EPA_END.AD=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPA_END") & tabBO2017$Stade%in%"AD"])/
    sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  #taux id END JV
  t2.EPA_END.JV=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPA_END") & tabBO2017$Stade%in%"JV"])/
    sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  
  
  
  
  #taux id LIEC_X SA
  t2.EPA_LIECX.SA=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPA_LIEC_X") & tabBO2017$Stade%in%"SA"])/
    sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  #taux id LIEC_X AD
  t2.EPA_LIECX.AD=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPA_LIEC_X") & tabBO2017$Stade%in%"AD"])/
    sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  #taux id LIEC_X JV
  t2.EPA_LIECX.JV=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPA_LIEC_X") & tabBO2017$Stade%in%"JV"])/
    sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  
  round(t2.EPA_EPI.AD+t2.EPA_EPI.SA+t2.EPA_EPI.JV+
          t2.EPA_EPA.AD+t2.EPA_EPA.SA+t2.EPA_EPA.JV+
          t2.EPA_ANS.AD+t2.EPA_ANS.SA+t2.EPA_ANS.JV+
          t2.EPA_END.AD+t2.EPA_END.SA+t2.EPA_END.JV+
          t2.EPA_LIECX.AD+t2.EPA_LIECX.SA+t2.EPA_LIECX.JV,5)==1
  
  
}


#les ANS identifi?s : que sont ils r?ellement ?
{
  
  #taux id EPI SA
  t2.ANS_EPI.SA=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("ANS_EPI") & tabBO2017$Stade%in%"SA"])/
    sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  #taux id EPI AD
  t2.ANS_EPI.AD=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("ANS_EPI") & tabBO2017$Stade%in%"AD"])/
    sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  #taux id EPI JV
  t2.ANS_EPI.JV=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("ANS_EPI") & tabBO2017$Stade%in%"JV"])/
    sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  
  
  
  #taux id EPA SA
  t2.ANS_EPA.SA=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("ANS_EPA") & tabBO2017$Stade%in%"SA"])/
    sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  #taux id EPA AD
  t2.ANS_EPA.AD=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("ANS_EPA") & tabBO2017$Stade%in%"AD"])/
    sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  #taux id EPA JV
  t2.ANS_EPA.JV=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("ANS_EPA") & tabBO2017$Stade%in%"JV"])/
    sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  
  
  
  #taux id ANS SA
  t2.ANS_ANS.SA=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("ANS_ANS") & tabBO2017$Stade%in%"SA"])/
    sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  #taux id ANS AD
  t2.ANS_ANS.AD=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("ANS_ANS") & tabBO2017$Stade%in%"AD"])/
    sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  #taux id ANS JV
  t2.ANS_ANS.JV=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("ANS_ANS") & tabBO2017$Stade%in%"JV"])/
    sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  
  
  
  #taux id END SA
  t2.ANS_END.SA=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("ANS_END") & tabBO2017$Stade%in%"SA"])/
    sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  #taux id END AD
  t2.ANS_END.AD=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("ANS_END") & tabBO2017$Stade%in%"AD"])/
    sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  #taux id END JV
  t2.ANS_END.JV=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("ANS_END") & tabBO2017$Stade%in%"JV"])/
    sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  
  
  
  
  #taux id LIEC_X SA
  t2.ANS_LIECX.SA=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("ANS_LIEC_X") & tabBO2017$Stade%in%"SA"])/
    sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  #taux id LIEC_X AD
  t2.ANS_LIECX.AD=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("ANS_LIEC_X") & tabBO2017$Stade%in%"AD"])/
    sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  #taux id LIEC_X JV
  t2.ANS_LIECX.JV=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("ANS_LIEC_X") & tabBO2017$Stade%in%"JV"])/
    sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  
  round(t2.ANS_EPI.AD+t2.ANS_EPI.SA+t2.ANS_EPI.JV+
          t2.ANS_EPA.AD+t2.ANS_EPA.SA+t2.ANS_EPA.JV+
          t2.ANS_ANS.AD+t2.ANS_ANS.SA+t2.ANS_ANS.JV+
          t2.ANS_END.AD+t2.ANS_END.SA+t2.ANS_END.JV+
          t2.ANS_LIECX.AD+t2.ANS_LIECX.SA+t2.ANS_LIECX.JV,5)==1
  
  
}


#les END identifi?s : que sont ils r?ellement ?
{
  
  #taux id EPI SA
  t2.END_EPI.SA=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("END_EPI") & tabBO2017$Stade%in%"SA"])/
    sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  #taux id EPI AD
  t2.END_EPI.AD=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("END_EPI") & tabBO2017$Stade%in%"AD"])/
    sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  #taux id EPI JV
  t2.END_EPI.JV=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("END_EPI") & tabBO2017$Stade%in%"JV"])/
    sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  
  
  
  #taux id EPA SA
  t2.END_EPA.SA=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("END_EPA") & tabBO2017$Stade%in%"SA"])/
    sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  #taux id EPA AD
  t2.END_EPA.AD=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("END_EPA") & tabBO2017$Stade%in%"AD"])/
    sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  #taux id EPA JV
  t2.END_EPA.JV=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("END_EPA") & tabBO2017$Stade%in%"JV"])/
    sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  
  
  
  #taux id ANS SA
  t2.END_ANS.SA=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("END_ANS") & tabBO2017$Stade%in%"SA"])/
    sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  #taux id ANS AD
  t2.END_ANS.AD=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("END_ANS") & tabBO2017$Stade%in%"AD"])/
    sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  #taux id ANS JV
  t2.END_ANS.JV=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("END_ANS") & tabBO2017$Stade%in%"JV"])/
    sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  
  
  
  #taux id END SA
  t2.END_END.SA=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("END_END") & tabBO2017$Stade%in%"SA"])/
    sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  #taux id END AD
  t2.END_END.AD=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("END_END") & tabBO2017$Stade%in%"AD"])/
    sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  #taux id END JV
  t2.END_END.JV=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("END_END") & tabBO2017$Stade%in%"JV"])/
    sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  
  
  
  
  #taux id LIEC_X SA
  t2.END_LIECX.SA=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("END_LIEC_X") & tabBO2017$Stade%in%"SA"])/
    sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  #taux id LIEC_X AD
  t2.END_LIECX.AD=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("END_LIEC_X") & tabBO2017$Stade%in%"AD"])/
    sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  #taux id LIEC_X JV
  t2.END_LIECX.JV=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("END_LIEC_X") & tabBO2017$Stade%in%"JV"])/
    sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  
  round(t2.END_EPI.AD+t2.END_EPI.SA+t2.END_EPI.JV+
          t2.END_EPA.AD+t2.END_EPA.SA+t2.END_EPA.JV+
          t2.END_ANS.AD+t2.END_ANS.SA+t2.END_ANS.JV+
          t2.END_END.AD+t2.END_END.SA+t2.END_END.JV+
          t2.END_LIECX.AD+t2.END_LIECX.SA+t2.END_LIECX.JV,5)==1
  
  
}


#les LIECX identifi?s : que sont ils r?ellement ?
{
  
  #taux id EPI SA
  t2.LIECX_EPI.SA=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("LIEC_X_EPI") & tabBO2017$Stade%in%"SA"])/
    sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  #taux id EPI AD
  t2.LIECX_EPI.AD=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("LIEC_X_EPI") & tabBO2017$Stade%in%"AD"])/
    sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  #taux id EPI JV
  t2.LIECX_EPI.JV=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("LIEC_X_EPI") & tabBO2017$Stade%in%"JV"])/
    sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  
  
  
  #taux id EPA SA
  t2.LIECX_EPA.SA=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("LIEC_X_EPA") & tabBO2017$Stade%in%"SA"])/
    sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  #taux id EPA AD
  t2.LIECX_EPA.AD=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("LIEC_X_EPA") & tabBO2017$Stade%in%"AD"])/
    sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  #taux id EPA JV
  t2.LIECX_EPA.JV=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("LIEC_X_EPA") & tabBO2017$Stade%in%"JV"])/
    sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  
  
  
  #taux id ANS SA
  t2.LIECX_ANS.SA=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("LIEC_X_ANS") & tabBO2017$Stade%in%"SA"])/
    sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  #taux id ANS AD
  t2.LIECX_ANS.AD=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("LIEC_X_ANS") & tabBO2017$Stade%in%"AD"])/
    sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  #taux id ANS JV
  t2.LIECX_ANS.JV=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("LIEC_X_ANS") & tabBO2017$Stade%in%"JV"])/
    sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  
  
  
  #taux id END SA
  t2.LIECX_END.SA=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("LIEC_X_END") & tabBO2017$Stade%in%"SA"])/
    sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  #taux id END AD
  t2.LIECX_END.AD=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("LIEC_X_END") & tabBO2017$Stade%in%"AD"])/
    sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  #taux id END JV
  t2.LIECX_END.JV=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("LIEC_X_END") & tabBO2017$Stade%in%"JV"])/
    sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  
  
  
  
  #taux id LIEC_X SA
  t2.LIECX_LIECX.SA=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("LIEC_X_LIEC_X") & tabBO2017$Stade%in%"SA"])/
    sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  #taux id LIEC_X AD
  t2.LIECX_LIECX.AD=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("LIEC_X_LIEC_X") & tabBO2017$Stade%in%"AD"])/
    sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  #taux id LIEC_X JV
  t2.LIECX_LIECX.JV=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("LIEC_X_LIEC_X") & tabBO2017$Stade%in%"JV"])/
    sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  
  round(t2.LIECX_EPI.AD+t2.LIECX_EPI.SA+t2.LIECX_EPI.JV+
          t2.LIECX_EPA.AD+t2.LIECX_EPA.SA+t2.LIECX_EPA.JV+
          t2.LIECX_ANS.AD+t2.LIECX_ANS.SA+t2.LIECX_ANS.JV+
          t2.LIECX_END.AD+t2.LIECX_END.SA+t2.LIECX_END.JV+
          t2.LIECX_LIECX.AD+t2.LIECX_LIECX.SA+t2.LIECX_LIECX.JV,5)==1
  
  
}


tabGraphBO17_2=data.frame(
  Taux=c(
    t2.EPI_EPI.AD,t2.EPI_EPI.SA,t2.EPI_EPI.JV,
    t2.EPI_EPA.AD,t2.EPI_EPA.SA,t2.EPI_EPA.JV,
    t2.EPI_ANS.AD,t2.EPI_ANS.SA,t2.EPI_ANS.JV,
    t2.EPI_END.AD,t2.EPI_END.SA,t2.EPI_END.JV,
    t2.EPI_LIECX.AD,t2.EPI_LIECX.SA,t2.EPI_LIECX.JV,
    
    t2.EPA_EPI.AD,t2.EPA_EPI.SA,t2.EPA_EPI.JV,
    t2.EPA_EPA.AD,t2.EPA_EPA.SA,t2.EPA_EPA.JV,
    t2.EPA_ANS.AD,t2.EPA_ANS.SA,t2.EPA_ANS.JV,
    t2.EPA_END.AD,t2.EPA_END.SA,t2.EPA_END.JV,
    t2.EPA_LIECX.AD,t2.EPA_LIECX.SA,t2.EPA_LIECX.JV,
    
    t2.ANS_EPI.AD,t2.ANS_EPI.SA,t2.ANS_EPI.JV,
    t2.ANS_EPA.AD,t2.ANS_EPA.SA,t2.ANS_EPA.JV,
    t2.ANS_ANS.AD,t2.ANS_ANS.SA,t2.ANS_ANS.JV,
    t2.ANS_END.AD,t2.ANS_END.SA,t2.ANS_END.JV,
    t2.ANS_LIECX.AD,t2.ANS_LIECX.SA,t2.ANS_LIECX.JV,
    
    t2.END_EPI.AD,t2.END_EPI.SA,t2.END_EPI.JV,
    t2.END_EPA.AD,t2.END_EPA.SA,t2.END_EPA.JV,
    t2.END_ANS.AD,t2.END_ANS.SA,t2.END_ANS.JV,
    t2.END_END.AD,t2.END_END.SA,t2.END_END.JV,
    t2.END_LIECX.AD,t2.END_LIECX.SA,t2.END_LIECX.JV,
    
    t2.LIECX_EPI.AD,t2.LIECX_EPI.SA,t2.LIECX_EPI.JV,
    t2.LIECX_EPA.AD,t2.LIECX_EPA.SA,t2.LIECX_EPA.JV,
    t2.LIECX_ANS.AD,t2.LIECX_ANS.SA,t2.LIECX_ANS.JV,
    t2.LIECX_END.AD,t2.LIECX_END.SA,t2.LIECX_END.JV,
    t2.LIECX_LIECX.AD,t2.LIECX_LIECX.SA,t2.LIECX_LIECX.JV),
  
  Classe=c(
    "t2.EPI_EPI.AD","t2.EPI_EPI.SA","t2.EPI_EPI.JV",
    "t2.EPI_EPA.AD","t2.EPI_EPA.SA","t2.EPI_EPA.JV",
    "t2.EPI_ANS.AD","t2.EPI_ANS.SA","t2.EPI_ANS.JV",
    "t2.EPI_END.AD","t2.EPI_END.SA","t2.EPI_END.JV",
    "t2.EPI_LIECX.AD","t2.EPI_LIECX.SA","t2.EPI_LIECX.JV",
    
    "t2.EPA_EPI.AD","t2.EPA_EPI.SA","t2.EPA_EPI.JV",
    "t2.EPA_EPA.AD","t2.EPA_EPA.SA","t2.EPA_EPA.JV",
    "t2.EPA_ANS.AD","t2.EPA_ANS.SA","t2.EPA_ANS.JV",
    "t2.EPA_END.AD","t2.EPA_END.SA","t2.EPA_END.JV",
    "t2.EPA_LIECX.AD","t2.EPA_LIECX.SA","t2.EPA_LIECX.JV",
    
    "t2.ANS_EPI.AD","t2.ANS_EPI.SA","t2.ANS_EPI.JV",
    "t2.ANS_EPA.AD","t2.ANS_EPA.SA","t2.ANS_EPA.JV",
    "t2.ANS_ANS.AD","t2.ANS_ANS.SA","t2.ANS_ANS.JV",
    "t2.ANS_END.AD","t2.ANS_END.SA","t2.ANS_END.JV",
    "t2.ANS_LIECX.AD","t2.ANS_LIECX.SA","t2.ANS_LIECX.JV",
    
    "t2.END_EPI.AD","t2.END_EPI.SA","t2.END_EPI.JV",
    "t2.END_EPA.AD","t2.END_EPA.SA","t2.END_EPA.JV",
    "t2.END_ANS.AD","t2.END_ANS.SA","t2.END_ANS.JV",
    "t2.END_END.AD","t2.END_END.SA","t2.END_END.JV",
    "t2.END_LIECX.AD","t2.END_LIECX.SA","t2.END_LIECX.JV",
    
    "t2.LIECX_EPI.AD","t2.LIECX_EPI.SA","t2.LIECX_EPI.JV",
    "t2.LIECX_EPA.AD","t2.LIECX_EPA.SA","t2.LIECX_EPA.JV",
    "t2.LIECX_ANS.AD","t2.LIECX_ANS.SA","t2.LIECX_ANS.JV",
    "t2.LIECX_END.AD","t2.LIECX_END.SA","t2.LIECX_END.JV",
    "t2.LIECX_LIECX.AD","t2.LIECX_LIECX.SA","t2.LIECX_LIECX.JV")
  
)
tabGraphBO17_2$LIEC=factor(substr(tabGraphBO17_2$Classe,8,10),levels=c("EPI","EPA","ANS","END","LIECX"))
tabGraphBO17_2$Stade=substr(tabGraphBO17_2$Classe,12,13)
tabGraphBO17_2$LIECid=factor(substr(tabGraphBO17_2$Classe,4,6),levels=c("EPI","EPA","ANS","END","LIECX"))
tabGraphBO17_2$LIEC_Stade=paste(tabGraphBO17_2$LIEC,tabGraphBO17_2$Stade,sep="_")

tabGraphBO17_2$LIEC_Stade=factor(tabGraphBO17_2$LIEC_Stade,levels=c(
  "EPI_AD","EPI_SA","EPI_JV","EPA_AD","EPA_SA",
  "EPA_JV","ANS_AD","ANS_SA","ANS_JV","END_AD",
  "END_SA","END_JV","LIECX_AD","LIECX_SA","LIECX_JV"))

ggplot(tabGraphBO17_2, aes(fill=LIEC_Stade, y=Taux, x=LIECid)) + 
  geom_bar(position='stack', stat='identity')+
  theme_classic()+
  scale_fill_manual(values = c(
    rgb(0,0,153,maxColorValue = 255),    rgb(0,0,153,maxColorValue = 255,alpha=150),    rgb(0,0,153,maxColorValue = 255,alpha=100),
    rgb(255,255,0,maxColorValue = 255),  rgb(255,255,0,maxColorValue = 255,alpha=100),  rgb(255,255,0,maxColorValue = 255,alpha=25),
    rgb(0,153,0,maxColorValue = 255),    rgb(0,153,0,maxColorValue = 255,alpha=150),    rgb(0,153,0,maxColorValue = 255,alpha=100),
    rgb(255,102,153,maxColorValue = 255),rgb(255,102,153,maxColorValue = 255,alpha=150),rgb(255,102,153,maxColorValue = 255,alpha=100),
    rgb(123,123,123,maxColorValue = 255),rgb(123,123,123,maxColorValue = 255,alpha=150),rgb(123,123,123,maxColorValue = 255,alpha=100)
  ))

ggplot(tabGraphBO17_2, aes(fill=LIEC_Stade, y=Taux, x=LIECid)) + 
  geom_bar(position='stack', stat='identity')+
  theme_classic()+
  scale_fill_manual(values = c(
    rgb(0,0,153,maxColorValue = 255),    rgb(0,0,153,maxColorValue = 255),    rgb(0,0,153,maxColorValue = 255,alpha=150),
    rgb(255,255,0,maxColorValue = 255),  rgb(255,255,0,maxColorValue = 255),  rgb(255,255,0,maxColorValue = 255,alpha=100),
    rgb(0,153,0,maxColorValue = 255),    rgb(0,153,0,maxColorValue = 255),    rgb(0,153,0,maxColorValue = 255,alpha=150),
    rgb(255,102,153,maxColorValue = 255),rgb(255,102,153,maxColorValue = 255),rgb(255,102,153,maxColorValue = 255,alpha=150),
    rgb(123,123,123,maxColorValue = 255),rgb(123,123,123,maxColorValue = 255),rgb(123,123,123,maxColorValue = 255,alpha=150)
  ))

#zoom sur les erreurs
tabGraphBO17_2$Erreur=NA
tabGraphBO17_2$Erreur[tabGraphBO17_2$LIEC==tabGraphBO17_2$LIECid]="Correcte"
tabGraphBO17_2$Erreur[tabGraphBO17_2$LIEC!=tabGraphBO17_2$LIECid]="Incorrecte"

ggplot(tabGraphBO17_2[tabGraphBO17_2$Erreur=="Incorrecte",], aes(fill=LIEC_Stade, y=Taux, x=LIECid)) + 
  geom_bar(position='fill', stat='identity')+
  theme_classic()+
  scale_fill_manual(values = c(
    rgb(0,0,153,maxColorValue = 255),    rgb(0,0,153,maxColorValue = 255,alpha=150),    rgb(0,0,153,maxColorValue = 255,alpha=100),
    rgb(255,255,0,maxColorValue = 255),  rgb(255,255,0,maxColorValue = 255,alpha=100),  rgb(255,255,0,maxColorValue = 255,alpha=25),
    rgb(0,153,0,maxColorValue = 255),    rgb(0,153,0,maxColorValue = 255,alpha=150),    rgb(0,153,0,maxColorValue = 255,alpha=100),
    rgb(255,102,153,maxColorValue = 255),rgb(255,102,153,maxColorValue = 255,alpha=150),rgb(255,102,153,maxColorValue = 255,alpha=100),
    rgb(123,123,123,maxColorValue = 255),rgb(123,123,123,maxColorValue = 255,alpha=150),rgb(123,123,123,maxColorValue = 255,alpha=100)
  ))




##### BO 2016 ABSOLUE - Parmis ce que les gens ont class?s en EPI EPA ANS END : qu'est ce qu'on a r?ellement ? #####




tabBO=Tab2[Tab2$Region%in%"BO",]
nrow(tabBO)==3426
dataID <- unique(sbt_m[,c("cle",colID,"Region")])
nrow(dataID)==1109
dataID=dataID[dataID$Region%in%"BO",]
nrow(dataID)==144



tabBO2016=tabBO[tabBO$Annee%in%2016,]
unique(tabBO2016$Stade)#"JV" "SA" "AD" => ok
any(is.na(tabBO2016$Stade))==T
nrow(tabBO2016)==488

#taux d'id correcte global
sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPI_EPI","EPA_EPA","ANS_ANS","END_END")])/
  sum(tabBO2016$Nbr_EW)#0.8867060  => ok
#vectReussiteBO

tabBO2016=tabBO2016[!is.na(tabBO2016$Stade),]
any(is.na(tabBO2016$Stade))==F
nrow(tabBO2016)==488-3





#les EPI identifi?s : que sont ils r?ellement ?
{
  
  #taux id correcte EPI SA
  ab.EPI_EPI.SA=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPI_EPI") & tabBO2016$Stade%in%"SA"])
  
  #taux id correcte EPI AD
  ab.EPI_EPI.AD=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPI_EPI") & tabBO2016$Stade%in%"AD"])
  
  #taux id correcte EPI JV
  ab.EPI_EPI.JV=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPI_EPI") & tabBO2016$Stade%in%"JV"])
  
  
  
  
  #taux id incorrecte EPA SA
  ab.EPI_EPA.SA=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPI_EPA") & tabBO2016$Stade%in%"SA"])
  
  #taux id incorrecte EPA AD
  ab.EPI_EPA.AD=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPI_EPA") & tabBO2016$Stade%in%"AD"])
  
  #taux id incorrecte EPA JV
  ab.EPI_EPA.JV=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPI_EPA") & tabBO2016$Stade%in%"JV"])
  
  
  
  
  #taux id incorrecte ANS SA
  ab.EPI_ANS.SA=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPI_ANS") & tabBO2016$Stade%in%"SA"])
  
  #taux id incorrecte ANS AD
  ab.EPI_ANS.AD=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPI_ANS") & tabBO2016$Stade%in%"AD"])
  
  #taux id incorrecte ANS JV
  ab.EPI_ANS.JV=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPI_ANS") & tabBO2016$Stade%in%"JV"])
  
  
  
  
  #taux id incorrecte END SA
  ab.EPI_END.SA=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPI_END") & tabBO2016$Stade%in%"SA"])
  
  #taux id incorrecte END AD
  ab.EPI_END.AD=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPI_END") & tabBO2016$Stade%in%"AD"])
  
  #taux id incorrecte END JV
  ab.EPI_END.JV=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPI_END") & tabBO2016$Stade%in%"JV"])
  
  
  
  
  
  #taux id incorrecte LIEC_X SA
  ab.EPI_LIECX.SA=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPI_LIEC_X") & tabBO2016$Stade%in%"SA"])
  
  #taux id incorrecte LIEC_X AD
  ab.EPI_LIECX.AD=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPI_LIEC_X") & tabBO2016$Stade%in%"AD"])
  
  #taux id incorrecte LIEC_X JV
  ab.EPI_LIECX.JV=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPI_LIEC_X") & tabBO2016$Stade%in%"JV"])
  
  
}


#les EPA identifi?s : que sont ils r?ellement ?
{
  
  #taux id EPI SA
  ab.EPA_EPI.SA=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPA_EPI") & tabBO2016$Stade%in%"SA"])
  
  #taux id EPI AD
  ab.EPA_EPI.AD=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPA_EPI") & tabBO2016$Stade%in%"AD"])
  
  #taux id EPI JV
  ab.EPA_EPI.JV=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPA_EPI") & tabBO2016$Stade%in%"JV"])
  
  
  
  
  #taux id EPA SA
  ab.EPA_EPA.SA=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPA_EPA") & tabBO2016$Stade%in%"SA"])
  
  #taux id EPA AD
  ab.EPA_EPA.AD=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPA_EPA") & tabBO2016$Stade%in%"AD"])
  
  #taux id EPA JV
  ab.EPA_EPA.JV=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPA_EPA") & tabBO2016$Stade%in%"JV"])
  
  
  
  
  #taux id ANS SA
  ab.EPA_ANS.SA=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPA_ANS") & tabBO2016$Stade%in%"SA"])
  
  #taux id ANS AD
  ab.EPA_ANS.AD=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPA_ANS") & tabBO2016$Stade%in%"AD"])
  
  #taux id ANS JV
  ab.EPA_ANS.JV=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPA_ANS") & tabBO2016$Stade%in%"JV"])
  
  
  
  
  #taux id END SA
  ab.EPA_END.SA=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPA_END") & tabBO2016$Stade%in%"SA"])
  
  #taux id END AD
  ab.EPA_END.AD=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPA_END") & tabBO2016$Stade%in%"AD"])
  
  #taux id END JV
  ab.EPA_END.JV=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPA_END") & tabBO2016$Stade%in%"JV"])
  
  
  
  
  
  #taux id LIEC_X SA
  ab.EPA_LIECX.SA=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPA_LIEC_X") & tabBO2016$Stade%in%"SA"])
  
  #taux id LIEC_X AD
  ab.EPA_LIECX.AD=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPA_LIEC_X") & tabBO2016$Stade%in%"AD"])
  
  #taux id LIEC_X JV
  ab.EPA_LIECX.JV=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("EPA_LIEC_X") & tabBO2016$Stade%in%"JV"])
  
  
  
  
}


#les ANS identifi?s : que sont ils r?ellement ?
{
  
  #taux id EPI SA
  ab.ANS_EPI.SA=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("ANS_EPI") & tabBO2016$Stade%in%"SA"])
  
  #taux id EPI AD
  ab.ANS_EPI.AD=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("ANS_EPI") & tabBO2016$Stade%in%"AD"])
  
  #taux id EPI JV
  ab.ANS_EPI.JV=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("ANS_EPI") & tabBO2016$Stade%in%"JV"])
  
  
  
  
  #taux id EPA SA
  ab.ANS_EPA.SA=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("ANS_EPA") & tabBO2016$Stade%in%"SA"])
  
  #taux id EPA AD
  ab.ANS_EPA.AD=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("ANS_EPA") & tabBO2016$Stade%in%"AD"])
  
  #taux id EPA JV
  ab.ANS_EPA.JV=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("ANS_EPA") & tabBO2016$Stade%in%"JV"])
  
  
  
  
  #taux id ANS SA
  ab.ANS_ANS.SA=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("ANS_ANS") & tabBO2016$Stade%in%"SA"])
  
  #taux id ANS AD
  ab.ANS_ANS.AD=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("ANS_ANS") & tabBO2016$Stade%in%"AD"])
  
  #taux id ANS JV
  ab.ANS_ANS.JV=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("ANS_ANS") & tabBO2016$Stade%in%"JV"])
  
  
  
  
  #taux id END SA
  ab.ANS_END.SA=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("ANS_END") & tabBO2016$Stade%in%"SA"])
  
  #taux id END AD
  ab.ANS_END.AD=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("ANS_END") & tabBO2016$Stade%in%"AD"])
  
  #taux id END JV
  ab.ANS_END.JV=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("ANS_END") & tabBO2016$Stade%in%"JV"])
  
  
  
  
  
  #taux id LIEC_X SA
  ab.ANS_LIECX.SA=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("ANS_LIEC_X") & tabBO2016$Stade%in%"SA"])
  
  #taux id LIEC_X AD
  ab.ANS_LIECX.AD=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("ANS_LIEC_X") & tabBO2016$Stade%in%"AD"])
  
  #taux id LIEC_X JV
  ab.ANS_LIECX.JV=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("ANS_LIEC_X") & tabBO2016$Stade%in%"JV"])
  
  
  
  
}


#les END identifi?s : que sont ils r?ellement ?
{
  
  #taux id EPI SA
  ab.END_EPI.SA=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("END_EPI") & tabBO2016$Stade%in%"SA"])
  
  #taux id EPI AD
  ab.END_EPI.AD=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("END_EPI") & tabBO2016$Stade%in%"AD"])
  
  #taux id EPI JV
  ab.END_EPI.JV=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("END_EPI") & tabBO2016$Stade%in%"JV"])
  
  
  
  
  #taux id EPA SA
  ab.END_EPA.SA=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("END_EPA") & tabBO2016$Stade%in%"SA"])
  
  #taux id EPA AD
  ab.END_EPA.AD=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("END_EPA") & tabBO2016$Stade%in%"AD"])
  
  #taux id EPA JV
  ab.END_EPA.JV=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("END_EPA") & tabBO2016$Stade%in%"JV"])
  
  
  
  
  #taux id ANS SA
  ab.END_ANS.SA=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("END_ANS") & tabBO2016$Stade%in%"SA"])
  
  #taux id ANS AD
  ab.END_ANS.AD=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("END_ANS") & tabBO2016$Stade%in%"AD"])
  
  #taux id ANS JV
  ab.END_ANS.JV=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("END_ANS") & tabBO2016$Stade%in%"JV"])
  
  
  
  
  #taux id END SA
  ab.END_END.SA=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("END_END") & tabBO2016$Stade%in%"SA"])
  
  #taux id END AD
  ab.END_END.AD=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("END_END") & tabBO2016$Stade%in%"AD"])
  
  #taux id END JV
  ab.END_END.JV=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("END_END") & tabBO2016$Stade%in%"JV"])
  
  
  
  
  #taux id LIEC_X SA
  ab.END_LIECX.SA=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("END_LIEC_X") & tabBO2016$Stade%in%"SA"])
  
  #taux id LIEC_X AD
  ab.END_LIECX.AD=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("END_LIEC_X") & tabBO2016$Stade%in%"AD"])
  
  #taux id LIEC_X JV
  ab.END_LIECX.JV=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("END_LIEC_X") & tabBO2016$Stade%in%"JV"])
  
  
  
  
}


#les LIECX identifi?s : que sont ils r?ellement ?
{
  
  #taux id EPI SA
  ab.LIECX_EPI.SA=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("LIEC_X_EPI") & tabBO2016$Stade%in%"SA"])
  
  #taux id EPI AD
  ab.LIECX_EPI.AD=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("LIEC_X_EPI") & tabBO2016$Stade%in%"AD"])
  
  #taux id EPI JV
  ab.LIECX_EPI.JV=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("LIEC_X_EPI") & tabBO2016$Stade%in%"JV"])
  
  
  
  
  #taux id EPA SA
  ab.LIECX_EPA.SA=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("LIEC_X_EPA") & tabBO2016$Stade%in%"SA"])
  
  #taux id EPA AD
  ab.LIECX_EPA.AD=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("LIEC_X_EPA") & tabBO2016$Stade%in%"AD"])
  
  #taux id EPA JV
  ab.LIECX_EPA.JV=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("LIEC_X_EPA") & tabBO2016$Stade%in%"JV"])
  
  
  
  
  #taux id ANS SA
  ab.LIECX_ANS.SA=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("LIEC_X_ANS") & tabBO2016$Stade%in%"SA"])
  
  #taux id ANS AD
  ab.LIECX_ANS.AD=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("LIEC_X_ANS") & tabBO2016$Stade%in%"AD"])
  
  #taux id ANS JV
  ab.LIECX_ANS.JV=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("LIEC_X_ANS") & tabBO2016$Stade%in%"JV"])
  
  
  
  
  #taux id END SA
  ab.LIECX_END.SA=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("LIEC_X_END") & tabBO2016$Stade%in%"SA"])
  
  #taux id END AD
  ab.LIECX_END.AD=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("LIEC_X_END") & tabBO2016$Stade%in%"AD"])
  
  #taux id END JV
  ab.LIECX_END.JV=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("LIEC_X_END") & tabBO2016$Stade%in%"JV"])
  
  
  
  
  
  #taux id LIEC_X SA
  ab.LIECX_LIECX.SA=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("LIEC_X_LIEC_X") & tabBO2016$Stade%in%"SA"])
  
  #taux id LIEC_X AD
  ab.LIECX_LIECX.AD=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("LIEC_X_LIEC_X") & tabBO2016$Stade%in%"AD"])
  
  #taux id LIEC_X JV
  ab.LIECX_LIECX.JV=sum(tabBO2016$Nbr_EW[tabBO2016$FIEC_LIEC%in%c("LIEC_X_LIEC_X") & tabBO2016$Stade%in%"JV"])
  
  
  
  
}





tabGraphBO16_2=data.frame(
  Taux=c(
    ab.EPI_EPI.AD,ab.EPI_EPI.SA,ab.EPI_EPI.JV,
    ab.EPI_EPA.AD,ab.EPI_EPA.SA,ab.EPI_EPA.JV,
    ab.EPI_ANS.AD,ab.EPI_ANS.SA,ab.EPI_ANS.JV,
    ab.EPI_END.AD,ab.EPI_END.SA,ab.EPI_END.JV,
    ab.EPI_LIECX.AD,ab.EPI_LIECX.SA,ab.EPI_LIECX.JV,
    
    ab.EPA_EPI.AD,ab.EPA_EPI.SA,ab.EPA_EPI.JV,
    ab.EPA_EPA.AD,ab.EPA_EPA.SA,ab.EPA_EPA.JV,
    ab.EPA_ANS.AD,ab.EPA_ANS.SA,ab.EPA_ANS.JV,
    ab.EPA_END.AD,ab.EPA_END.SA,ab.EPA_END.JV,
    ab.EPA_LIECX.AD,ab.EPA_LIECX.SA,ab.EPA_LIECX.JV,
    
    ab.ANS_EPI.AD,ab.ANS_EPI.SA,ab.ANS_EPI.JV,
    ab.ANS_EPA.AD,ab.ANS_EPA.SA,ab.ANS_EPA.JV,
    ab.ANS_ANS.AD,ab.ANS_ANS.SA,ab.ANS_ANS.JV,
    ab.ANS_END.AD,ab.ANS_END.SA,ab.ANS_END.JV,
    ab.ANS_LIECX.AD,ab.ANS_LIECX.SA,ab.ANS_LIECX.JV,
    
    ab.END_EPI.AD,ab.END_EPI.SA,ab.END_EPI.JV,
    ab.END_EPA.AD,ab.END_EPA.SA,ab.END_EPA.JV,
    ab.END_ANS.AD,ab.END_ANS.SA,ab.END_ANS.JV,
    ab.END_END.AD,ab.END_END.SA,ab.END_END.JV,
    ab.END_LIECX.AD,ab.END_LIECX.SA,ab.END_LIECX.JV,
    
    ab.LIECX_EPI.AD,ab.LIECX_EPI.SA,ab.LIECX_EPI.JV,
    ab.LIECX_EPA.AD,ab.LIECX_EPA.SA,ab.LIECX_EPA.JV,
    ab.LIECX_ANS.AD,ab.LIECX_ANS.SA,ab.LIECX_ANS.JV,
    ab.LIECX_END.AD,ab.LIECX_END.SA,ab.LIECX_END.JV,
    ab.LIECX_LIECX.AD,ab.LIECX_LIECX.SA,ab.LIECX_LIECX.JV),
  
  Classe=c(
    "ab.EPI_EPI.AD","ab.EPI_EPI.SA","ab.EPI_EPI.JV",
    "ab.EPI_EPA.AD","ab.EPI_EPA.SA","ab.EPI_EPA.JV",
    "ab.EPI_ANS.AD","ab.EPI_ANS.SA","ab.EPI_ANS.JV",
    "ab.EPI_END.AD","ab.EPI_END.SA","ab.EPI_END.JV",
    "ab.EPI_LIECX.AD","ab.EPI_LIECX.SA","ab.EPI_LIECX.JV",
    
    "ab.EPA_EPI.AD","ab.EPA_EPI.SA","ab.EPA_EPI.JV",
    "ab.EPA_EPA.AD","ab.EPA_EPA.SA","ab.EPA_EPA.JV",
    "ab.EPA_ANS.AD","ab.EPA_ANS.SA","ab.EPA_ANS.JV",
    "ab.EPA_END.AD","ab.EPA_END.SA","ab.EPA_END.JV",
    "ab.EPA_LIECX.AD","ab.EPA_LIECX.SA","ab.EPA_LIECX.JV",
    
    "ab.ANS_EPI.AD","ab.ANS_EPI.SA","ab.ANS_EPI.JV",
    "ab.ANS_EPA.AD","ab.ANS_EPA.SA","ab.ANS_EPA.JV",
    "ab.ANS_ANS.AD","ab.ANS_ANS.SA","ab.ANS_ANS.JV",
    "ab.ANS_END.AD","ab.ANS_END.SA","ab.ANS_END.JV",
    "ab.ANS_LIECX.AD","ab.ANS_LIECX.SA","ab.ANS_LIECX.JV",
    
    "ab.END_EPI.AD","ab.END_EPI.SA","ab.END_EPI.JV",
    "ab.END_EPA.AD","ab.END_EPA.SA","ab.END_EPA.JV",
    "ab.END_ANS.AD","ab.END_ANS.SA","ab.END_ANS.JV",
    "ab.END_END.AD","ab.END_END.SA","ab.END_END.JV",
    "ab.END_LIECX.AD","ab.END_LIECX.SA","ab.END_LIECX.JV",
    
    "ab.LIECX_EPI.AD","ab.LIECX_EPI.SA","ab.LIECX_EPI.JV",
    "ab.LIECX_EPA.AD","ab.LIECX_EPA.SA","ab.LIECX_EPA.JV",
    "ab.LIECX_ANS.AD","ab.LIECX_ANS.SA","ab.LIECX_ANS.JV",
    "ab.LIECX_END.AD","ab.LIECX_END.SA","ab.LIECX_END.JV",
    "ab.LIECX_LIECX.AD","ab.LIECX_LIECX.SA","ab.LIECX_LIECX.JV")
  
)

tabGraphBO16_2$LIEC=factor(substr(tabGraphBO16_2$Classe,8,10),levels=c("EPI","EPA","ANS","END","LIECX"))
tabGraphBO16_2$Stade=substr(tabGraphBO16_2$Classe,12,13)
tabGraphBO16_2$LIECid=factor(substr(tabGraphBO16_2$Classe,4,6),levels=c("EPI","EPA","ANS","END","LIECX"))
tabGraphBO16_2$LIEC_Stade=paste(tabGraphBO16_2$LIEC,tabGraphBO16_2$Stade,sep="_")

tabGraphBO16_2$LIEC_Stade=factor(tabGraphBO16_2$LIEC_Stade,levels=c(
  "EPI_AD","EPI_SA","EPI_JV","EPA_AD","EPA_SA",
  "EPA_JV","ANS_AD","ANS_SA","ANS_JV","END_AD",
  "END_SA","END_JV","LIECX_AD","LIECX_SA","LIECX_JV"))

ggplot(tabGraphBO16_2, aes(fill=LIEC_Stade, y=Taux, x=LIECid)) + 
  geom_bar(position='stack', stat='identity')+
  theme_classic()+ylab("Nombre de VDT")+
  scale_fill_manual(values = c(
    rgb(0,0,153,maxColorValue = 255),    rgb(0,0,153,maxColorValue = 255),    rgb(0,0,153,maxColorValue = 255,alpha=150),
    rgb(255,255,0,maxColorValue = 255),  rgb(255,255,0,maxColorValue = 255),  rgb(255,255,0,maxColorValue = 255,alpha=100),
    rgb(0,153,0,maxColorValue = 255),    rgb(0,153,0,maxColorValue = 255),    rgb(0,153,0,maxColorValue = 255,alpha=150),
    rgb(255,102,153,maxColorValue = 255),rgb(255,102,153,maxColorValue = 255),rgb(255,102,153,maxColorValue = 255,alpha=150),
    rgb(123,123,123,maxColorValue = 255),rgb(123,123,123,maxColorValue = 255),rgb(123,123,123,maxColorValue = 255,alpha=150)
  ))

#zoom sur les erreurs
tabGraphBO16_2$Erreur=NA
tabGraphBO16_2$Erreur[tabGraphBO16_2$LIEC==tabGraphBO16_2$LIECid]="Correcte"
tabGraphBO16_2$Erreur[tabGraphBO16_2$LIEC!=tabGraphBO16_2$LIECid]="Incorrecte"

ggplot(tabGraphBO16_2[tabGraphBO16_2$Erreur=="Incorrecte",], aes(fill=LIEC_Stade, y=Taux, x=LIECid)) + 
  geom_bar(position='fill', stat='identity')+
  theme_classic()+
  scale_fill_manual(values = c(
    rgb(0,0,153,maxColorValue = 255),    rgb(0,0,153,maxColorValue = 255,alpha=150),    rgb(0,0,153,maxColorValue = 255,alpha=100),
    rgb(255,255,0,maxColorValue = 255),  rgb(255,255,0,maxColorValue = 255,alpha=100),  rgb(255,255,0,maxColorValue = 255,alpha=25),
    rgb(0,153,0,maxColorValue = 255),    rgb(0,153,0,maxColorValue = 255,alpha=150),    rgb(0,153,0,maxColorValue = 255,alpha=100),
    rgb(255,102,153,maxColorValue = 255),rgb(255,102,153,maxColorValue = 255,alpha=150),rgb(255,102,153,maxColorValue = 255,alpha=100),
    rgb(123,123,123,maxColorValue = 255),rgb(123,123,123,maxColorValue = 255,alpha=150),rgb(123,123,123,maxColorValue = 255,alpha=100)
  ))




##### BO 2017 ABSOLUE - Parmis ce que les gens ont class?s en EPI EPA ANS END : qu'est ce qu'on a r?ellement ? #####



tabBO=Tab2[Tab2$Region%in%"BO",]
nrow(tabBO)==3426
dataID <- unique(sbt_m[,c("cle",colID,"Region")])
nrow(dataID)==1109
dataID=dataID[dataID$Region%in%"BO",]
nrow(dataID)==144



tabBO2017=tabBO[tabBO$Annee%in%2017,]
unique(tabBO2017$Stade)#"JV" "SA" "AD" => ok
any(is.na(tabBO2017$Stade))==F
nrow(tabBO2017)==471

#taux d'id correcte global
sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPI_EPI","EPA_EPA","ANS_ANS","END_END")])/
  sum(tabBO2017$Nbr_EW)#0.6357143   => ok
#vectReussiteBO

#les EPI identifi?s : que sont ils r?ellement ?
{
  
  #taux id correcte EPI SA
  ab.EPI_EPI.SA=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPI_EPI") & tabBO2017$Stade%in%"SA"])
  
  #taux id correcte EPI AD
  ab.EPI_EPI.AD=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPI_EPI") & tabBO2017$Stade%in%"AD"])
  
  #taux id correcte EPI JV
  ab.EPI_EPI.JV=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPI_EPI") & tabBO2017$Stade%in%"JV"])
  
  
  
  
  #taux id incorrecte EPA SA
  ab.EPI_EPA.SA=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPI_EPA") & tabBO2017$Stade%in%"SA"])
  
  #taux id incorrecte EPA AD
  ab.EPI_EPA.AD=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPI_EPA") & tabBO2017$Stade%in%"AD"])
  
  #taux id incorrecte EPA JV
  ab.EPI_EPA.JV=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPI_EPA") & tabBO2017$Stade%in%"JV"])
  
  
  
  
  #taux id incorrecte ANS SA
  ab.EPI_ANS.SA=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPI_ANS") & tabBO2017$Stade%in%"SA"])
  
  #taux id incorrecte ANS AD
  ab.EPI_ANS.AD=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPI_ANS") & tabBO2017$Stade%in%"AD"])
  
  #taux id incorrecte ANS JV
  ab.EPI_ANS.JV=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPI_ANS") & tabBO2017$Stade%in%"JV"])
  
  
  
  
  #taux id incorrecte END SA
  ab.EPI_END.SA=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPI_END") & tabBO2017$Stade%in%"SA"])
  
  #taux id incorrecte END AD
  ab.EPI_END.AD=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPI_END") & tabBO2017$Stade%in%"AD"])
  
  #taux id incorrecte END JV
  ab.EPI_END.JV=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPI_END") & tabBO2017$Stade%in%"JV"])
  
  
  
  
  
  #taux id incorrecte LIEC_X SA
  ab.EPI_LIECX.SA=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPI_LIEC_X") & tabBO2017$Stade%in%"SA"])
  
  #taux id incorrecte LIEC_X AD
  ab.EPI_LIECX.AD=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPI_LIEC_X") & tabBO2017$Stade%in%"AD"])
  
  #taux id incorrecte LIEC_X JV
  ab.EPI_LIECX.JV=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPI_LIEC_X") & tabBO2017$Stade%in%"JV"])
  
  
  }


#les EPA identifi?s : que sont ils r?ellement ?
{
  
  #taux id EPI SA
  ab.EPA_EPI.SA=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPA_EPI") & tabBO2017$Stade%in%"SA"])
  
  #taux id EPI AD
  ab.EPA_EPI.AD=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPA_EPI") & tabBO2017$Stade%in%"AD"])
  
  #taux id EPI JV
  ab.EPA_EPI.JV=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPA_EPI") & tabBO2017$Stade%in%"JV"])
  
  
  
  
  #taux id EPA SA
  ab.EPA_EPA.SA=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPA_EPA") & tabBO2017$Stade%in%"SA"])
  
  #taux id EPA AD
  ab.EPA_EPA.AD=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPA_EPA") & tabBO2017$Stade%in%"AD"])
  
  #taux id EPA JV
  ab.EPA_EPA.JV=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPA_EPA") & tabBO2017$Stade%in%"JV"])
  
  
  
  
  #taux id ANS SA
  ab.EPA_ANS.SA=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPA_ANS") & tabBO2017$Stade%in%"SA"])
  
  #taux id ANS AD
  ab.EPA_ANS.AD=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPA_ANS") & tabBO2017$Stade%in%"AD"])
  
  #taux id ANS JV
  ab.EPA_ANS.JV=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPA_ANS") & tabBO2017$Stade%in%"JV"])
  
  
  
  
  #taux id END SA
  ab.EPA_END.SA=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPA_END") & tabBO2017$Stade%in%"SA"])
  
  #taux id END AD
  ab.EPA_END.AD=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPA_END") & tabBO2017$Stade%in%"AD"])
  
  #taux id END JV
  ab.EPA_END.JV=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPA_END") & tabBO2017$Stade%in%"JV"])
  
  
  
  
  
  #taux id LIEC_X SA
  ab.EPA_LIECX.SA=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPA_LIEC_X") & tabBO2017$Stade%in%"SA"])
  
  #taux id LIEC_X AD
  ab.EPA_LIECX.AD=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPA_LIEC_X") & tabBO2017$Stade%in%"AD"])
  
  #taux id LIEC_X JV
  ab.EPA_LIECX.JV=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("EPA_LIEC_X") & tabBO2017$Stade%in%"JV"])
  
  
  
  
}


#les ANS identifi?s : que sont ils r?ellement ?
{
  
  #taux id EPI SA
  ab.ANS_EPI.SA=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("ANS_EPI") & tabBO2017$Stade%in%"SA"])
  
  #taux id EPI AD
  ab.ANS_EPI.AD=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("ANS_EPI") & tabBO2017$Stade%in%"AD"])
  
  #taux id EPI JV
  ab.ANS_EPI.JV=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("ANS_EPI") & tabBO2017$Stade%in%"JV"])
  
  
  
  
  #taux id EPA SA
  ab.ANS_EPA.SA=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("ANS_EPA") & tabBO2017$Stade%in%"SA"])
  
  #taux id EPA AD
  ab.ANS_EPA.AD=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("ANS_EPA") & tabBO2017$Stade%in%"AD"])
  
  #taux id EPA JV
  ab.ANS_EPA.JV=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("ANS_EPA") & tabBO2017$Stade%in%"JV"])
  
  
  
  
  #taux id ANS SA
  ab.ANS_ANS.SA=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("ANS_ANS") & tabBO2017$Stade%in%"SA"])
  
  #taux id ANS AD
  ab.ANS_ANS.AD=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("ANS_ANS") & tabBO2017$Stade%in%"AD"])
  
  #taux id ANS JV
  ab.ANS_ANS.JV=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("ANS_ANS") & tabBO2017$Stade%in%"JV"])
  
  
  
  
  #taux id END SA
  ab.ANS_END.SA=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("ANS_END") & tabBO2017$Stade%in%"SA"])
  
  #taux id END AD
  ab.ANS_END.AD=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("ANS_END") & tabBO2017$Stade%in%"AD"])
  
  #taux id END JV
  ab.ANS_END.JV=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("ANS_END") & tabBO2017$Stade%in%"JV"])
  
  
  
  
  
  #taux id LIEC_X SA
  ab.ANS_LIECX.SA=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("ANS_LIEC_X") & tabBO2017$Stade%in%"SA"])
  
  #taux id LIEC_X AD
  ab.ANS_LIECX.AD=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("ANS_LIEC_X") & tabBO2017$Stade%in%"AD"])
  
  #taux id LIEC_X JV
  ab.ANS_LIECX.JV=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("ANS_LIEC_X") & tabBO2017$Stade%in%"JV"])
  
  
  
  
}


#les END identifi?s : que sont ils r?ellement ?
{
  
  #taux id EPI SA
  ab.END_EPI.SA=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("END_EPI") & tabBO2017$Stade%in%"SA"])
  
  #taux id EPI AD
  ab.END_EPI.AD=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("END_EPI") & tabBO2017$Stade%in%"AD"])
  
  #taux id EPI JV
  ab.END_EPI.JV=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("END_EPI") & tabBO2017$Stade%in%"JV"])
  
  
  
  
  #taux id EPA SA
  ab.END_EPA.SA=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("END_EPA") & tabBO2017$Stade%in%"SA"])
  
  #taux id EPA AD
  ab.END_EPA.AD=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("END_EPA") & tabBO2017$Stade%in%"AD"])
  
  #taux id EPA JV
  ab.END_EPA.JV=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("END_EPA") & tabBO2017$Stade%in%"JV"])
  
  
  
  
  #taux id ANS SA
  ab.END_ANS.SA=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("END_ANS") & tabBO2017$Stade%in%"SA"])
  
  #taux id ANS AD
  ab.END_ANS.AD=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("END_ANS") & tabBO2017$Stade%in%"AD"])
  
  #taux id ANS JV
  ab.END_ANS.JV=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("END_ANS") & tabBO2017$Stade%in%"JV"])
  
  
  
  
  #taux id END SA
  ab.END_END.SA=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("END_END") & tabBO2017$Stade%in%"SA"])
  
  #taux id END AD
  ab.END_END.AD=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("END_END") & tabBO2017$Stade%in%"AD"])
  
  #taux id END JV
  ab.END_END.JV=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("END_END") & tabBO2017$Stade%in%"JV"])
  
  
  
  
  #taux id LIEC_X SA
  ab.END_LIECX.SA=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("END_LIEC_X") & tabBO2017$Stade%in%"SA"])
  
  #taux id LIEC_X AD
  ab.END_LIECX.AD=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("END_LIEC_X") & tabBO2017$Stade%in%"AD"])
  
  #taux id LIEC_X JV
  ab.END_LIECX.JV=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("END_LIEC_X") & tabBO2017$Stade%in%"JV"])
  
  
  
  
}


#les LIECX identifi?s : que sont ils r?ellement ?
{
  
  #taux id EPI SA
  ab.LIECX_EPI.SA=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("LIEC_X_EPI") & tabBO2017$Stade%in%"SA"])
  
  #taux id EPI AD
  ab.LIECX_EPI.AD=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("LIEC_X_EPI") & tabBO2017$Stade%in%"AD"])
  
  #taux id EPI JV
  ab.LIECX_EPI.JV=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("LIEC_X_EPI") & tabBO2017$Stade%in%"JV"])
  
  
  
  
  #taux id EPA SA
  ab.LIECX_EPA.SA=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("LIEC_X_EPA") & tabBO2017$Stade%in%"SA"])
  
  #taux id EPA AD
  ab.LIECX_EPA.AD=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("LIEC_X_EPA") & tabBO2017$Stade%in%"AD"])
  
  #taux id EPA JV
  ab.LIECX_EPA.JV=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("LIEC_X_EPA") & tabBO2017$Stade%in%"JV"])
  
  
  
  
  #taux id ANS SA
  ab.LIECX_ANS.SA=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("LIEC_X_ANS") & tabBO2017$Stade%in%"SA"])
  
  #taux id ANS AD
  ab.LIECX_ANS.AD=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("LIEC_X_ANS") & tabBO2017$Stade%in%"AD"])
  
  #taux id ANS JV
  ab.LIECX_ANS.JV=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("LIEC_X_ANS") & tabBO2017$Stade%in%"JV"])
  
  
  
  
  #taux id END SA
  ab.LIECX_END.SA=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("LIEC_X_END") & tabBO2017$Stade%in%"SA"])
  
  #taux id END AD
  ab.LIECX_END.AD=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("LIEC_X_END") & tabBO2017$Stade%in%"AD"])
  
  #taux id END JV
  ab.LIECX_END.JV=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("LIEC_X_END") & tabBO2017$Stade%in%"JV"])
  
  
  
  
  
  #taux id LIEC_X SA
  ab.LIECX_LIECX.SA=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("LIEC_X_LIEC_X") & tabBO2017$Stade%in%"SA"])
  
  #taux id LIEC_X AD
  ab.LIECX_LIECX.AD=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("LIEC_X_LIEC_X") & tabBO2017$Stade%in%"AD"])
  
  #taux id LIEC_X JV
  ab.LIECX_LIECX.JV=sum(tabBO2017$Nbr_EW[tabBO2017$FIEC_LIEC%in%c("LIEC_X_LIEC_X") & tabBO2017$Stade%in%"JV"])
  
  
  
  
}





tabGraphBO17_2=data.frame(
  Taux=c(
    ab.EPI_EPI.AD,ab.EPI_EPI.SA,ab.EPI_EPI.JV,
    ab.EPI_EPA.AD,ab.EPI_EPA.SA,ab.EPI_EPA.JV,
    ab.EPI_ANS.AD,ab.EPI_ANS.SA,ab.EPI_ANS.JV,
    ab.EPI_END.AD,ab.EPI_END.SA,ab.EPI_END.JV,
    ab.EPI_LIECX.AD,ab.EPI_LIECX.SA,ab.EPI_LIECX.JV,
    
    ab.EPA_EPI.AD,ab.EPA_EPI.SA,ab.EPA_EPI.JV,
    ab.EPA_EPA.AD,ab.EPA_EPA.SA,ab.EPA_EPA.JV,
    ab.EPA_ANS.AD,ab.EPA_ANS.SA,ab.EPA_ANS.JV,
    ab.EPA_END.AD,ab.EPA_END.SA,ab.EPA_END.JV,
    ab.EPA_LIECX.AD,ab.EPA_LIECX.SA,ab.EPA_LIECX.JV,
    
    ab.ANS_EPI.AD,ab.ANS_EPI.SA,ab.ANS_EPI.JV,
    ab.ANS_EPA.AD,ab.ANS_EPA.SA,ab.ANS_EPA.JV,
    ab.ANS_ANS.AD,ab.ANS_ANS.SA,ab.ANS_ANS.JV,
    ab.ANS_END.AD,ab.ANS_END.SA,ab.ANS_END.JV,
    ab.ANS_LIECX.AD,ab.ANS_LIECX.SA,ab.ANS_LIECX.JV,
    
    ab.END_EPI.AD,ab.END_EPI.SA,ab.END_EPI.JV,
    ab.END_EPA.AD,ab.END_EPA.SA,ab.END_EPA.JV,
    ab.END_ANS.AD,ab.END_ANS.SA,ab.END_ANS.JV,
    ab.END_END.AD,ab.END_END.SA,ab.END_END.JV,
    ab.END_LIECX.AD,ab.END_LIECX.SA,ab.END_LIECX.JV,
    
    ab.LIECX_EPI.AD,ab.LIECX_EPI.SA,ab.LIECX_EPI.JV,
    ab.LIECX_EPA.AD,ab.LIECX_EPA.SA,ab.LIECX_EPA.JV,
    ab.LIECX_ANS.AD,ab.LIECX_ANS.SA,ab.LIECX_ANS.JV,
    ab.LIECX_END.AD,ab.LIECX_END.SA,ab.LIECX_END.JV,
    ab.LIECX_LIECX.AD,ab.LIECX_LIECX.SA,ab.LIECX_LIECX.JV),
  
  Classe=c(
    "ab.EPI_EPI.AD","ab.EPI_EPI.SA","ab.EPI_EPI.JV",
    "ab.EPI_EPA.AD","ab.EPI_EPA.SA","ab.EPI_EPA.JV",
    "ab.EPI_ANS.AD","ab.EPI_ANS.SA","ab.EPI_ANS.JV",
    "ab.EPI_END.AD","ab.EPI_END.SA","ab.EPI_END.JV",
    "ab.EPI_LIECX.AD","ab.EPI_LIECX.SA","ab.EPI_LIECX.JV",
    
    "ab.EPA_EPI.AD","ab.EPA_EPI.SA","ab.EPA_EPI.JV",
    "ab.EPA_EPA.AD","ab.EPA_EPA.SA","ab.EPA_EPA.JV",
    "ab.EPA_ANS.AD","ab.EPA_ANS.SA","ab.EPA_ANS.JV",
    "ab.EPA_END.AD","ab.EPA_END.SA","ab.EPA_END.JV",
    "ab.EPA_LIECX.AD","ab.EPA_LIECX.SA","ab.EPA_LIECX.JV",
    
    "ab.ANS_EPI.AD","ab.ANS_EPI.SA","ab.ANS_EPI.JV",
    "ab.ANS_EPA.AD","ab.ANS_EPA.SA","ab.ANS_EPA.JV",
    "ab.ANS_ANS.AD","ab.ANS_ANS.SA","ab.ANS_ANS.JV",
    "ab.ANS_END.AD","ab.ANS_END.SA","ab.ANS_END.JV",
    "ab.ANS_LIECX.AD","ab.ANS_LIECX.SA","ab.ANS_LIECX.JV",
    
    "ab.END_EPI.AD","ab.END_EPI.SA","ab.END_EPI.JV",
    "ab.END_EPA.AD","ab.END_EPA.SA","ab.END_EPA.JV",
    "ab.END_ANS.AD","ab.END_ANS.SA","ab.END_ANS.JV",
    "ab.END_END.AD","ab.END_END.SA","ab.END_END.JV",
    "ab.END_LIECX.AD","ab.END_LIECX.SA","ab.END_LIECX.JV",
    
    "ab.LIECX_EPI.AD","ab.LIECX_EPI.SA","ab.LIECX_EPI.JV",
    "ab.LIECX_EPA.AD","ab.LIECX_EPA.SA","ab.LIECX_EPA.JV",
    "ab.LIECX_ANS.AD","ab.LIECX_ANS.SA","ab.LIECX_ANS.JV",
    "ab.LIECX_END.AD","ab.LIECX_END.SA","ab.LIECX_END.JV",
    "ab.LIECX_LIECX.AD","ab.LIECX_LIECX.SA","ab.LIECX_LIECX.JV")
  
)

tabGraphBO17_2$LIEC=factor(substr(tabGraphBO17_2$Classe,8,10),levels=c("EPI","EPA","ANS","END","LIECX"))
tabGraphBO17_2$Stade=substr(tabGraphBO17_2$Classe,12,13)
tabGraphBO17_2$LIECid=factor(substr(tabGraphBO17_2$Classe,4,6),levels=c("EPI","EPA","ANS","END","LIECX"))
tabGraphBO17_2$LIEC_Stade=paste(tabGraphBO17_2$LIEC,tabGraphBO17_2$Stade,sep="_")

tabGraphBO17_2$LIEC_Stade=factor(tabGraphBO17_2$LIEC_Stade,levels=c(
  "EPI_AD","EPI_SA","EPI_JV","EPA_AD","EPA_SA",
  "EPA_JV","ANS_AD","ANS_SA","ANS_JV","END_AD",
  "END_SA","END_JV","LIECX_AD","LIECX_SA","LIECX_JV"))

ggplot(tabGraphBO17_2, aes(fill=LIEC_Stade, y=Taux, x=LIECid)) + 
  geom_bar(position='stack', stat='identity')+
  theme_classic()+ylab("Nombre de VDT")+
  scale_fill_manual(values = c(
    rgb(0,0,153,maxColorValue = 255),    rgb(0,0,153,maxColorValue = 255),    rgb(0,0,153,maxColorValue = 255,alpha=150),
    rgb(255,255,0,maxColorValue = 255),  rgb(255,255,0,maxColorValue = 255),  rgb(255,255,0,maxColorValue = 255,alpha=100),
    rgb(0,153,0,maxColorValue = 255),    rgb(0,153,0,maxColorValue = 255),    rgb(0,153,0,maxColorValue = 255,alpha=150),
    rgb(255,102,153,maxColorValue = 255),rgb(255,102,153,maxColorValue = 255),rgb(255,102,153,maxColorValue = 255,alpha=150),
    rgb(123,123,123,maxColorValue = 255),rgb(123,123,123,maxColorValue = 255),rgb(123,123,123,maxColorValue = 255,alpha=150)
  ))

#zoom sur les erreurs
tabGraphBO17_2$Erreur=NA
tabGraphBO17_2$Erreur[tabGraphBO17_2$LIEC==tabGraphBO17_2$LIECid]="Correcte"
tabGraphBO17_2$Erreur[tabGraphBO17_2$LIEC!=tabGraphBO17_2$LIECid]="Incorrecte"

ggplot(tabGraphBO17_2[tabGraphBO17_2$Erreur=="Incorrecte",], aes(fill=LIEC_Stade, y=Taux, x=LIECid)) + 
  geom_bar(position='fill', stat='identity')+
  theme_classic()+
  scale_fill_manual(values = c(
    rgb(0,0,153,maxColorValue = 255),    rgb(0,0,153,maxColorValue = 255,alpha=150),    rgb(0,0,153,maxColorValue = 255,alpha=100),
    rgb(255,255,0,maxColorValue = 255),  rgb(255,255,0,maxColorValue = 255,alpha=100),  rgb(255,255,0,maxColorValue = 255,alpha=25),
    rgb(0,153,0,maxColorValue = 255),    rgb(0,153,0,maxColorValue = 255,alpha=150),    rgb(0,153,0,maxColorValue = 255,alpha=100),
    rgb(255,102,153,maxColorValue = 255),rgb(255,102,153,maxColorValue = 255,alpha=150),rgb(255,102,153,maxColorValue = 255,alpha=100),
    rgb(123,123,123,maxColorValue = 255),rgb(123,123,123,maxColorValue = 255,alpha=150),rgb(123,123,123,maxColorValue = 255,alpha=100)
  ))




##### NAT 2014 zoom sur les points critiques - Parmis ce que les gens ont class?s en EPI EPA ANS END : qu'est ce qu'on a r?ellement ? #####



nrow(Tab2)==21044
dataID <- unique(sbt_m[,c("cle",colID,"Region")])
nrow(dataID)==1109



tabNAT2014=Tab2[Tab2$Annee%in%2014,]
unique(tabNAT2014$Stade)#"JV" "SA" "AD" => ok
any(is.na(tabNAT2014$Stade))==F
nrow(tabNAT2014)==2938

#taux d'id correcte global
sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPI_EPI","EPA_EPA","ANS_ANS","END_END")])/
  sum(tabNAT2014$Nbr_EW)#0.6466317   => ok
#vectReussiteBO

#les EPI identifi?s : que sont ils r?ellement ?
{
  
  #taux id correcte EPI SA
  t2.EPI_EPI.SA=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPI_EPI") & tabNAT2014$Stade%in%"SA"])/
    sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  #taux id correcte EPI AD
  t2.EPI_EPI.AD=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPI_EPI") & tabNAT2014$Stade%in%"AD"])/
    sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  #taux id correcte EPI JV
  t2.EPI_EPI.JV=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPI_EPI") & tabNAT2014$Stade%in%"JV"])/
    sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  
  
  
  #taux id incorrecte EPA SA
  t2.EPI_EPA.SA=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPI_EPA") & tabNAT2014$Stade%in%"SA"])/
    sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  #taux id incorrecte EPA AD
  t2.EPI_EPA.AD=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPI_EPA") & tabNAT2014$Stade%in%"AD"])/
    sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  #taux id incorrecte EPA JV
  t2.EPI_EPA.JV=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPI_EPA") & tabNAT2014$Stade%in%"JV"])/
    sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  
  
  
  #taux id incorrecte ANS SA
  t2.EPI_ANS.SA=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPI_ANS") & tabNAT2014$Stade%in%"SA"])/
    sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  #taux id incorrecte ANS AD
  t2.EPI_ANS.AD=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPI_ANS") & tabNAT2014$Stade%in%"AD"])/
    sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  #taux id incorrecte ANS JV
  t2.EPI_ANS.JV=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPI_ANS") & tabNAT2014$Stade%in%"JV"])/
    sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  
  
  
  #taux id incorrecte END SA
  t2.EPI_END.SA=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPI_END") & tabNAT2014$Stade%in%"SA"])/
    sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  #taux id incorrecte END AD
  t2.EPI_END.AD=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPI_END") & tabNAT2014$Stade%in%"AD"])/
    sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  #taux id incorrecte END JV
  t2.EPI_END.JV=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPI_END") & tabNAT2014$Stade%in%"JV"])/
    sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  
  
  
  
  #taux id incorrecte LIEC_X SA
  t2.EPI_LIECX.SA=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPI_LIEC_X") & tabNAT2014$Stade%in%"SA"])/
    sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  #taux id incorrecte LIEC_X AD
  t2.EPI_LIECX.AD=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPI_LIEC_X") & tabNAT2014$Stade%in%"AD"])/
    sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  #taux id incorrecte LIEC_X JV
  t2.EPI_LIECX.JV=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPI_LIEC_X") & tabNAT2014$Stade%in%"JV"])/
    sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  
  t2.EPI_EPI.AD+t2.EPI_EPI.SA+t2.EPI_EPI.JV+
    t2.EPI_EPA.AD+t2.EPI_EPA.SA+t2.EPI_EPA.JV+
    t2.EPI_ANS.AD+t2.EPI_ANS.SA+t2.EPI_ANS.JV+
    t2.EPI_END.AD+t2.EPI_END.SA+t2.EPI_END.JV+
    t2.EPI_LIECX.AD+t2.EPI_LIECX.SA+t2.EPI_LIECX.JV==1
  
  
  }


#les EPA identifi?s : que sont ils r?ellement ?
{
  
  #taux id EPI SA
  t2.EPA_EPI.SA=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPA_EPI") & tabNAT2014$Stade%in%"SA"])/
    sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  #taux id EPI AD
  t2.EPA_EPI.AD=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPA_EPI") & tabNAT2014$Stade%in%"AD"])/
    sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  #taux id EPI JV
  t2.EPA_EPI.JV=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPA_EPI") & tabNAT2014$Stade%in%"JV"])/
    sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  
  
  
  #taux id EPA SA
  t2.EPA_EPA.SA=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPA_EPA") & tabNAT2014$Stade%in%"SA"])/
    sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  #taux id EPA AD
  t2.EPA_EPA.AD=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPA_EPA") & tabNAT2014$Stade%in%"AD"])/
    sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  #taux id EPA JV
  t2.EPA_EPA.JV=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPA_EPA") & tabNAT2014$Stade%in%"JV"])/
    sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  
  
  
  #taux id ANS SA
  t2.EPA_ANS.SA=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPA_ANS") & tabNAT2014$Stade%in%"SA"])/
    sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  #taux id ANS AD
  t2.EPA_ANS.AD=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPA_ANS") & tabNAT2014$Stade%in%"AD"])/
    sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  #taux id ANS JV
  t2.EPA_ANS.JV=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPA_ANS") & tabNAT2014$Stade%in%"JV"])/
    sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  
  
  
  #taux id END SA
  t2.EPA_END.SA=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPA_END") & tabNAT2014$Stade%in%"SA"])/
    sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  #taux id END AD
  t2.EPA_END.AD=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPA_END") & tabNAT2014$Stade%in%"AD"])/
    sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  #taux id END JV
  t2.EPA_END.JV=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPA_END") & tabNAT2014$Stade%in%"JV"])/
    sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  
  
  
  
  #taux id LIEC_X SA
  t2.EPA_LIECX.SA=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPA_LIEC_X") & tabNAT2014$Stade%in%"SA"])/
    sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  #taux id LIEC_X AD
  t2.EPA_LIECX.AD=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPA_LIEC_X") & tabNAT2014$Stade%in%"AD"])/
    sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  #taux id LIEC_X JV
  t2.EPA_LIECX.JV=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPA_LIEC_X") & tabNAT2014$Stade%in%"JV"])/
    sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  
  round(t2.EPA_EPI.AD+t2.EPA_EPI.SA+t2.EPA_EPI.JV+
          t2.EPA_EPA.AD+t2.EPA_EPA.SA+t2.EPA_EPA.JV+
          t2.EPA_ANS.AD+t2.EPA_ANS.SA+t2.EPA_ANS.JV+
          t2.EPA_END.AD+t2.EPA_END.SA+t2.EPA_END.JV+
          t2.EPA_LIECX.AD+t2.EPA_LIECX.SA+t2.EPA_LIECX.JV,5)==1
  
  
}


#les ANS identifi?s : que sont ils r?ellement ?
{
  
  #taux id EPI SA
  t2.ANS_EPI.SA=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("ANS_EPI") & tabNAT2014$Stade%in%"SA"])/
    sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  #taux id EPI AD
  t2.ANS_EPI.AD=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("ANS_EPI") & tabNAT2014$Stade%in%"AD"])/
    sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  #taux id EPI JV
  t2.ANS_EPI.JV=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("ANS_EPI") & tabNAT2014$Stade%in%"JV"])/
    sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  
  
  
  #taux id EPA SA
  t2.ANS_EPA.SA=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("ANS_EPA") & tabNAT2014$Stade%in%"SA"])/
    sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  #taux id EPA AD
  t2.ANS_EPA.AD=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("ANS_EPA") & tabNAT2014$Stade%in%"AD"])/
    sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  #taux id EPA JV
  t2.ANS_EPA.JV=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("ANS_EPA") & tabNAT2014$Stade%in%"JV"])/
    sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  
  
  
  #taux id ANS SA
  t2.ANS_ANS.SA=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("ANS_ANS") & tabNAT2014$Stade%in%"SA"])/
    sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  #taux id ANS AD
  t2.ANS_ANS.AD=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("ANS_ANS") & tabNAT2014$Stade%in%"AD"])/
    sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  #taux id ANS JV
  t2.ANS_ANS.JV=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("ANS_ANS") & tabNAT2014$Stade%in%"JV"])/
    sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  
  
  
  #taux id END SA
  t2.ANS_END.SA=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("ANS_END") & tabNAT2014$Stade%in%"SA"])/
    sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  #taux id END AD
  t2.ANS_END.AD=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("ANS_END") & tabNAT2014$Stade%in%"AD"])/
    sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  #taux id END JV
  t2.ANS_END.JV=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("ANS_END") & tabNAT2014$Stade%in%"JV"])/
    sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  
  
  
  
  #taux id LIEC_X SA
  t2.ANS_LIECX.SA=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("ANS_LIEC_X") & tabNAT2014$Stade%in%"SA"])/
    sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  #taux id LIEC_X AD
  t2.ANS_LIECX.AD=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("ANS_LIEC_X") & tabNAT2014$Stade%in%"AD"])/
    sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  #taux id LIEC_X JV
  t2.ANS_LIECX.JV=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("ANS_LIEC_X") & tabNAT2014$Stade%in%"JV"])/
    sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  
  round(t2.ANS_EPI.AD+t2.ANS_EPI.SA+t2.ANS_EPI.JV+
          t2.ANS_EPA.AD+t2.ANS_EPA.SA+t2.ANS_EPA.JV+
          t2.ANS_ANS.AD+t2.ANS_ANS.SA+t2.ANS_ANS.JV+
          t2.ANS_END.AD+t2.ANS_END.SA+t2.ANS_END.JV+
          t2.ANS_LIECX.AD+t2.ANS_LIECX.SA+t2.ANS_LIECX.JV,5)==1
  
  
}


#les END identifi?s : que sont ils r?ellement ?
{
  
  #taux id EPI SA
  t2.END_EPI.SA=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("END_EPI") & tabNAT2014$Stade%in%"SA"])/
    sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  #taux id EPI AD
  t2.END_EPI.AD=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("END_EPI") & tabNAT2014$Stade%in%"AD"])/
    sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  #taux id EPI JV
  t2.END_EPI.JV=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("END_EPI") & tabNAT2014$Stade%in%"JV"])/
    sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  
  
  
  #taux id EPA SA
  t2.END_EPA.SA=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("END_EPA") & tabNAT2014$Stade%in%"SA"])/
    sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  #taux id EPA AD
  t2.END_EPA.AD=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("END_EPA") & tabNAT2014$Stade%in%"AD"])/
    sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  #taux id EPA JV
  t2.END_EPA.JV=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("END_EPA") & tabNAT2014$Stade%in%"JV"])/
    sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  
  
  
  #taux id ANS SA
  t2.END_ANS.SA=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("END_ANS") & tabNAT2014$Stade%in%"SA"])/
    sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  #taux id ANS AD
  t2.END_ANS.AD=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("END_ANS") & tabNAT2014$Stade%in%"AD"])/
    sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  #taux id ANS JV
  t2.END_ANS.JV=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("END_ANS") & tabNAT2014$Stade%in%"JV"])/
    sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  
  
  
  #taux id END SA
  t2.END_END.SA=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("END_END") & tabNAT2014$Stade%in%"SA"])/
    sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  #taux id END AD
  t2.END_END.AD=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("END_END") & tabNAT2014$Stade%in%"AD"])/
    sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  #taux id END JV
  t2.END_END.JV=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("END_END") & tabNAT2014$Stade%in%"JV"])/
    sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  
  
  
  
  #taux id LIEC_X SA
  t2.END_LIECX.SA=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("END_LIEC_X") & tabNAT2014$Stade%in%"SA"])/
    sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  #taux id LIEC_X AD
  t2.END_LIECX.AD=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("END_LIEC_X") & tabNAT2014$Stade%in%"AD"])/
    sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  #taux id LIEC_X JV
  t2.END_LIECX.JV=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("END_LIEC_X") & tabNAT2014$Stade%in%"JV"])/
    sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  
  round(t2.END_EPI.AD+t2.END_EPI.SA+t2.END_EPI.JV+
          t2.END_EPA.AD+t2.END_EPA.SA+t2.END_EPA.JV+
          t2.END_ANS.AD+t2.END_ANS.SA+t2.END_ANS.JV+
          t2.END_END.AD+t2.END_END.SA+t2.END_END.JV+
          t2.END_LIECX.AD+t2.END_LIECX.SA+t2.END_LIECX.JV,5)==1
  
  
}


#les LIECX identifi?s : que sont ils r?ellement ?
{
  
  #taux id EPI SA
  t2.LIECX_EPI.SA=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("LIEC_X_EPI") & tabNAT2014$Stade%in%"SA"])/
    sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  #taux id EPI AD
  t2.LIECX_EPI.AD=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("LIEC_X_EPI") & tabNAT2014$Stade%in%"AD"])/
    sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  #taux id EPI JV
  t2.LIECX_EPI.JV=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("LIEC_X_EPI") & tabNAT2014$Stade%in%"JV"])/
    sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  
  
  
  #taux id EPA SA
  t2.LIECX_EPA.SA=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("LIEC_X_EPA") & tabNAT2014$Stade%in%"SA"])/
    sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  #taux id EPA AD
  t2.LIECX_EPA.AD=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("LIEC_X_EPA") & tabNAT2014$Stade%in%"AD"])/
    sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  #taux id EPA JV
  t2.LIECX_EPA.JV=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("LIEC_X_EPA") & tabNAT2014$Stade%in%"JV"])/
    sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  
  
  
  #taux id ANS SA
  t2.LIECX_ANS.SA=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("LIEC_X_ANS") & tabNAT2014$Stade%in%"SA"])/
    sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  #taux id ANS AD
  t2.LIECX_ANS.AD=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("LIEC_X_ANS") & tabNAT2014$Stade%in%"AD"])/
    sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  #taux id ANS JV
  t2.LIECX_ANS.JV=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("LIEC_X_ANS") & tabNAT2014$Stade%in%"JV"])/
    sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  
  
  
  #taux id END SA
  t2.LIECX_END.SA=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("LIEC_X_END") & tabNAT2014$Stade%in%"SA"])/
    sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  #taux id END AD
  t2.LIECX_END.AD=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("LIEC_X_END") & tabNAT2014$Stade%in%"AD"])/
    sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  #taux id END JV
  t2.LIECX_END.JV=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("LIEC_X_END") & tabNAT2014$Stade%in%"JV"])/
    sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  
  
  
  
  #taux id LIEC_X SA
  t2.LIECX_LIECX.SA=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("LIEC_X_LIEC_X") & tabNAT2014$Stade%in%"SA"])/
    sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  #taux id LIEC_X AD
  t2.LIECX_LIECX.AD=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("LIEC_X_LIEC_X") & tabNAT2014$Stade%in%"AD"])/
    sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  #taux id LIEC_X JV
  t2.LIECX_LIECX.JV=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("LIEC_X_LIEC_X") & tabNAT2014$Stade%in%"JV"])/
    sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  
  round(t2.LIECX_EPI.AD+t2.LIECX_EPI.SA+t2.LIECX_EPI.JV+
          t2.LIECX_EPA.AD+t2.LIECX_EPA.SA+t2.LIECX_EPA.JV+
          t2.LIECX_ANS.AD+t2.LIECX_ANS.SA+t2.LIECX_ANS.JV+
          t2.LIECX_END.AD+t2.LIECX_END.SA+t2.LIECX_END.JV+
          t2.LIECX_LIECX.AD+t2.LIECX_LIECX.SA+t2.LIECX_LIECX.JV,5)==1
  
  
}


tabGraphNAT14_2=data.frame(
  Taux=c(
    t2.EPI_EPI.AD,t2.EPI_EPI.SA,t2.EPI_EPI.JV,
    t2.EPI_EPA.AD,t2.EPI_EPA.SA,t2.EPI_EPA.JV,
    t2.EPI_ANS.AD,t2.EPI_ANS.SA,t2.EPI_ANS.JV,
    t2.EPI_END.AD,t2.EPI_END.SA,t2.EPI_END.JV,
    t2.EPI_LIECX.AD,t2.EPI_LIECX.SA,t2.EPI_LIECX.JV,
    
    t2.EPA_EPI.AD,t2.EPA_EPI.SA,t2.EPA_EPI.JV,
    t2.EPA_EPA.AD,t2.EPA_EPA.SA,t2.EPA_EPA.JV,
    t2.EPA_ANS.AD,t2.EPA_ANS.SA,t2.EPA_ANS.JV,
    t2.EPA_END.AD,t2.EPA_END.SA,t2.EPA_END.JV,
    t2.EPA_LIECX.AD,t2.EPA_LIECX.SA,t2.EPA_LIECX.JV,
    
    t2.ANS_EPI.AD,t2.ANS_EPI.SA,t2.ANS_EPI.JV,
    t2.ANS_EPA.AD,t2.ANS_EPA.SA,t2.ANS_EPA.JV,
    t2.ANS_ANS.AD,t2.ANS_ANS.SA,t2.ANS_ANS.JV,
    t2.ANS_END.AD,t2.ANS_END.SA,t2.ANS_END.JV,
    t2.ANS_LIECX.AD,t2.ANS_LIECX.SA,t2.ANS_LIECX.JV,
    
    t2.END_EPI.AD,t2.END_EPI.SA,t2.END_EPI.JV,
    t2.END_EPA.AD,t2.END_EPA.SA,t2.END_EPA.JV,
    t2.END_ANS.AD,t2.END_ANS.SA,t2.END_ANS.JV,
    t2.END_END.AD,t2.END_END.SA,t2.END_END.JV,
    t2.END_LIECX.AD,t2.END_LIECX.SA,t2.END_LIECX.JV,
    
    t2.LIECX_EPI.AD,t2.LIECX_EPI.SA,t2.LIECX_EPI.JV,
    t2.LIECX_EPA.AD,t2.LIECX_EPA.SA,t2.LIECX_EPA.JV,
    t2.LIECX_ANS.AD,t2.LIECX_ANS.SA,t2.LIECX_ANS.JV,
    t2.LIECX_END.AD,t2.LIECX_END.SA,t2.LIECX_END.JV,
    t2.LIECX_LIECX.AD,t2.LIECX_LIECX.SA,t2.LIECX_LIECX.JV),
  
  Classe=c(
    "t2.EPI_EPI.AD","t2.EPI_EPI.SA","t2.EPI_EPI.JV",
    "t2.EPI_EPA.AD","t2.EPI_EPA.SA","t2.EPI_EPA.JV",
    "t2.EPI_ANS.AD","t2.EPI_ANS.SA","t2.EPI_ANS.JV",
    "t2.EPI_END.AD","t2.EPI_END.SA","t2.EPI_END.JV",
    "t2.EPI_LIECX.AD","t2.EPI_LIECX.SA","t2.EPI_LIECX.JV",
    
    "t2.EPA_EPI.AD","t2.EPA_EPI.SA","t2.EPA_EPI.JV",
    "t2.EPA_EPA.AD","t2.EPA_EPA.SA","t2.EPA_EPA.JV",
    "t2.EPA_ANS.AD","t2.EPA_ANS.SA","t2.EPA_ANS.JV",
    "t2.EPA_END.AD","t2.EPA_END.SA","t2.EPA_END.JV",
    "t2.EPA_LIECX.AD","t2.EPA_LIECX.SA","t2.EPA_LIECX.JV",
    
    "t2.ANS_EPI.AD","t2.ANS_EPI.SA","t2.ANS_EPI.JV",
    "t2.ANS_EPA.AD","t2.ANS_EPA.SA","t2.ANS_EPA.JV",
    "t2.ANS_ANS.AD","t2.ANS_ANS.SA","t2.ANS_ANS.JV",
    "t2.ANS_END.AD","t2.ANS_END.SA","t2.ANS_END.JV",
    "t2.ANS_LIECX.AD","t2.ANS_LIECX.SA","t2.ANS_LIECX.JV",
    
    "t2.END_EPI.AD","t2.END_EPI.SA","t2.END_EPI.JV",
    "t2.END_EPA.AD","t2.END_EPA.SA","t2.END_EPA.JV",
    "t2.END_ANS.AD","t2.END_ANS.SA","t2.END_ANS.JV",
    "t2.END_END.AD","t2.END_END.SA","t2.END_END.JV",
    "t2.END_LIECX.AD","t2.END_LIECX.SA","t2.END_LIECX.JV",
    
    "t2.LIECX_EPI.AD","t2.LIECX_EPI.SA","t2.LIECX_EPI.JV",
    "t2.LIECX_EPA.AD","t2.LIECX_EPA.SA","t2.LIECX_EPA.JV",
    "t2.LIECX_ANS.AD","t2.LIECX_ANS.SA","t2.LIECX_ANS.JV",
    "t2.LIECX_END.AD","t2.LIECX_END.SA","t2.LIECX_END.JV",
    "t2.LIECX_LIECX.AD","t2.LIECX_LIECX.SA","t2.LIECX_LIECX.JV")
  
)

tabGraphNAT14_2$LIEC=factor(substr(tabGraphNAT14_2$Classe,8,10),levels=c("EPI","EPA","ANS","END","LIECX"))
tabGraphNAT14_2$Stade=substr(tabGraphNAT14_2$Classe,12,13)
tabGraphNAT14_2$LIECid=factor(substr(tabGraphNAT14_2$Classe,4,6),levels=c("EPI","EPA","ANS","END","LIECX"))
tabGraphNAT14_2$LIEC_Stade=paste(tabGraphNAT14_2$LIEC,tabGraphNAT14_2$Stade,sep="_")

tabGraphNAT14_2$LIEC_Stade=factor(tabGraphNAT14_2$LIEC_Stade,levels=c(
  "EPI_AD","EPI_SA","EPI_JV","EPA_AD","EPA_SA",
  "EPA_JV","ANS_AD","ANS_SA","ANS_JV","END_AD",
  "END_SA","END_JV","LIECX_AD","LIECX_SA","LIECX_JV"))

ggplot(tabGraphNAT14_2, aes(fill=LIEC_Stade, y=Taux, x=LIECid)) + 
  geom_bar(position='stack', stat='identity')+
  theme_classic()+
  scale_fill_manual(values = c(
    rgb(0,0,153,maxColorValue = 255),    rgb(0,0,153,maxColorValue = 255),    rgb(0,0,153,maxColorValue = 255,alpha=150),
    rgb(255,255,0,maxColorValue = 255),  rgb(255,255,0,maxColorValue = 255),  rgb(255,255,0,maxColorValue = 255,alpha=100),
    rgb(0,153,0,maxColorValue = 255),    rgb(0,153,0,maxColorValue = 255),    rgb(0,153,0,maxColorValue = 255,alpha=150),
    rgb(255,102,153,maxColorValue = 255),rgb(255,102,153,maxColorValue = 255),rgb(255,102,153,maxColorValue = 255,alpha=150),
    rgb(123,123,123,maxColorValue = 255),rgb(123,123,123,maxColorValue = 255),rgb(123,123,123,maxColorValue = 255,alpha=150)
  ))

#zoom sur les erreurs
tabGraphNAT14_2$Erreur=NA
tabGraphNAT14_2$Erreur[tabGraphNAT14_2$LIEC==tabGraphNAT14_2$LIECid]="Correcte"
tabGraphNAT14_2$Erreur[tabGraphNAT14_2$LIEC!=tabGraphNAT14_2$LIECid]="Incorrecte"

ggplot(tabGraphNAT14_2[tabGraphNAT14_2$Erreur=="Incorrecte",], aes(fill=LIEC_Stade, y=Taux, x=LIECid)) + 
  geom_bar(position='fill', stat='identity')+
  theme_classic()+
  scale_fill_manual(values = c(
    rgb(0,0,153,maxColorValue = 255),    rgb(0,0,153,maxColorValue = 255,alpha=150),    rgb(0,0,153,maxColorValue = 255,alpha=100),
    rgb(255,255,0,maxColorValue = 255),  rgb(255,255,0,maxColorValue = 255,alpha=100),  rgb(255,255,0,maxColorValue = 255,alpha=25),
    rgb(0,153,0,maxColorValue = 255),    rgb(0,153,0,maxColorValue = 255,alpha=150),    rgb(0,153,0,maxColorValue = 255,alpha=100),
    rgb(255,102,153,maxColorValue = 255),rgb(255,102,153,maxColorValue = 255,alpha=150),rgb(255,102,153,maxColorValue = 255,alpha=100),
    rgb(123,123,123,maxColorValue = 255),rgb(123,123,123,maxColorValue = 255,alpha=150),rgb(123,123,123,maxColorValue = 255,alpha=100)
  ))




##### NAT 2015 zoom sur les points critiques - Parmis ce que les gens ont class?s en EPI EPA ANS END : qu'est ce qu'on a r?ellement ? #####



nrow(Tab2)==21044
dataID <- unique(sbt_m[,c("cle",colID,"Region")])
nrow(dataID)==1109



tabNAT2015=Tab2[Tab2$Annee%in%2015,]
unique(tabNAT2015$Stade)#"JV" "SA" "AD" "X"
any(is.na(tabNAT2015$Stade))==F
nrow(tabNAT2015)==2751

#taux d'id correcte global
sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPI_EPI","EPA_EPA","ANS_ANS","END_END")])/
  sum(tabNAT2015$Nbr_EW)#0.8171038   => ok
#vectReussiteNAT

tabNAT2015=tabNAT2015[tabNAT2015$Stade!="X",]
nrow(tabNAT2015)==2751-11

#les EPI identifi?s : que sont ils r?ellement ?
{
  
  #taux id correcte EPI SA
  t2.EPI_EPI.SA=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPI_EPI") & tabNAT2015$Stade%in%"SA"])/
    sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  #taux id correcte EPI AD
  t2.EPI_EPI.AD=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPI_EPI") & tabNAT2015$Stade%in%"AD"])/
    sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  #taux id correcte EPI JV
  t2.EPI_EPI.JV=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPI_EPI") & tabNAT2015$Stade%in%"JV"])/
    sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  
  
  
  #taux id incorrecte EPA SA
  t2.EPI_EPA.SA=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPI_EPA") & tabNAT2015$Stade%in%"SA"])/
    sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  #taux id incorrecte EPA AD
  t2.EPI_EPA.AD=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPI_EPA") & tabNAT2015$Stade%in%"AD"])/
    sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  #taux id incorrecte EPA JV
  t2.EPI_EPA.JV=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPI_EPA") & tabNAT2015$Stade%in%"JV"])/
    sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  
  
  
  #taux id incorrecte ANS SA
  t2.EPI_ANS.SA=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPI_ANS") & tabNAT2015$Stade%in%"SA"])/
    sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  #taux id incorrecte ANS AD
  t2.EPI_ANS.AD=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPI_ANS") & tabNAT2015$Stade%in%"AD"])/
    sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  #taux id incorrecte ANS JV
  t2.EPI_ANS.JV=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPI_ANS") & tabNAT2015$Stade%in%"JV"])/
    sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  
  
  
  #taux id incorrecte END SA
  t2.EPI_END.SA=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPI_END") & tabNAT2015$Stade%in%"SA"])/
    sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  #taux id incorrecte END AD
  t2.EPI_END.AD=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPI_END") & tabNAT2015$Stade%in%"AD"])/
    sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  #taux id incorrecte END JV
  t2.EPI_END.JV=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPI_END") & tabNAT2015$Stade%in%"JV"])/
    sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  
  
  
  
  #taux id incorrecte LIEC_X SA
  t2.EPI_LIECX.SA=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPI_LIEC_X") & tabNAT2015$Stade%in%"SA"])/
    sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  #taux id incorrecte LIEC_X AD
  t2.EPI_LIECX.AD=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPI_LIEC_X") & tabNAT2015$Stade%in%"AD"])/
    sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  #taux id incorrecte LIEC_X JV
  t2.EPI_LIECX.JV=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPI_LIEC_X") & tabNAT2015$Stade%in%"JV"])/
    sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  
  t2.EPI_EPI.AD+t2.EPI_EPI.SA+t2.EPI_EPI.JV+
    t2.EPI_EPA.AD+t2.EPI_EPA.SA+t2.EPI_EPA.JV+
    t2.EPI_ANS.AD+t2.EPI_ANS.SA+t2.EPI_ANS.JV+
    t2.EPI_END.AD+t2.EPI_END.SA+t2.EPI_END.JV+
    t2.EPI_LIECX.AD+t2.EPI_LIECX.SA+t2.EPI_LIECX.JV==1
  
  
}


#les EPA identifi?s : que sont ils r?ellement ?
{
  
  #taux id EPI SA
  t2.EPA_EPI.SA=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPA_EPI") & tabNAT2015$Stade%in%"SA"])/
    sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  #taux id EPI AD
  t2.EPA_EPI.AD=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPA_EPI") & tabNAT2015$Stade%in%"AD"])/
    sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  #taux id EPI JV
  t2.EPA_EPI.JV=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPA_EPI") & tabNAT2015$Stade%in%"JV"])/
    sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  
  
  
  #taux id EPA SA
  t2.EPA_EPA.SA=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPA_EPA") & tabNAT2015$Stade%in%"SA"])/
    sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  #taux id EPA AD
  t2.EPA_EPA.AD=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPA_EPA") & tabNAT2015$Stade%in%"AD"])/
    sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  #taux id EPA JV
  t2.EPA_EPA.JV=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPA_EPA") & tabNAT2015$Stade%in%"JV"])/
    sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  
  
  
  #taux id ANS SA
  t2.EPA_ANS.SA=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPA_ANS") & tabNAT2015$Stade%in%"SA"])/
    sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  #taux id ANS AD
  t2.EPA_ANS.AD=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPA_ANS") & tabNAT2015$Stade%in%"AD"])/
    sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  #taux id ANS JV
  t2.EPA_ANS.JV=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPA_ANS") & tabNAT2015$Stade%in%"JV"])/
    sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  
  
  
  #taux id END SA
  t2.EPA_END.SA=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPA_END") & tabNAT2015$Stade%in%"SA"])/
    sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  #taux id END AD
  t2.EPA_END.AD=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPA_END") & tabNAT2015$Stade%in%"AD"])/
    sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  #taux id END JV
  t2.EPA_END.JV=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPA_END") & tabNAT2015$Stade%in%"JV"])/
    sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  
  
  
  
  #taux id LIEC_X SA
  t2.EPA_LIECX.SA=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPA_LIEC_X") & tabNAT2015$Stade%in%"SA"])/
    sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  #taux id LIEC_X AD
  t2.EPA_LIECX.AD=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPA_LIEC_X") & tabNAT2015$Stade%in%"AD"])/
    sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  #taux id LIEC_X JV
  t2.EPA_LIECX.JV=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPA_LIEC_X") & tabNAT2015$Stade%in%"JV"])/
    sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  
  round(t2.EPA_EPI.AD+t2.EPA_EPI.SA+t2.EPA_EPI.JV+
          t2.EPA_EPA.AD+t2.EPA_EPA.SA+t2.EPA_EPA.JV+
          t2.EPA_ANS.AD+t2.EPA_ANS.SA+t2.EPA_ANS.JV+
          t2.EPA_END.AD+t2.EPA_END.SA+t2.EPA_END.JV+
          t2.EPA_LIECX.AD+t2.EPA_LIECX.SA+t2.EPA_LIECX.JV,5)==1
  
  
}


#les ANS identifi?s : que sont ils r?ellement ?
{
  
  #taux id EPI SA
  t2.ANS_EPI.SA=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("ANS_EPI") & tabNAT2015$Stade%in%"SA"])/
    sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  #taux id EPI AD
  t2.ANS_EPI.AD=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("ANS_EPI") & tabNAT2015$Stade%in%"AD"])/
    sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  #taux id EPI JV
  t2.ANS_EPI.JV=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("ANS_EPI") & tabNAT2015$Stade%in%"JV"])/
    sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  
  
  
  #taux id EPA SA
  t2.ANS_EPA.SA=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("ANS_EPA") & tabNAT2015$Stade%in%"SA"])/
    sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  #taux id EPA AD
  t2.ANS_EPA.AD=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("ANS_EPA") & tabNAT2015$Stade%in%"AD"])/
    sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  #taux id EPA JV
  t2.ANS_EPA.JV=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("ANS_EPA") & tabNAT2015$Stade%in%"JV"])/
    sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  
  
  
  #taux id ANS SA
  t2.ANS_ANS.SA=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("ANS_ANS") & tabNAT2015$Stade%in%"SA"])/
    sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  #taux id ANS AD
  t2.ANS_ANS.AD=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("ANS_ANS") & tabNAT2015$Stade%in%"AD"])/
    sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  #taux id ANS JV
  t2.ANS_ANS.JV=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("ANS_ANS") & tabNAT2015$Stade%in%"JV"])/
    sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  
  
  
  #taux id END SA
  t2.ANS_END.SA=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("ANS_END") & tabNAT2015$Stade%in%"SA"])/
    sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  #taux id END AD
  t2.ANS_END.AD=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("ANS_END") & tabNAT2015$Stade%in%"AD"])/
    sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  #taux id END JV
  t2.ANS_END.JV=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("ANS_END") & tabNAT2015$Stade%in%"JV"])/
    sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  
  
  
  
  #taux id LIEC_X SA
  t2.ANS_LIECX.SA=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("ANS_LIEC_X") & tabNAT2015$Stade%in%"SA"])/
    sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  #taux id LIEC_X AD
  t2.ANS_LIECX.AD=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("ANS_LIEC_X") & tabNAT2015$Stade%in%"AD"])/
    sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  #taux id LIEC_X JV
  t2.ANS_LIECX.JV=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("ANS_LIEC_X") & tabNAT2015$Stade%in%"JV"])/
    sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  
  round(t2.ANS_EPI.AD+t2.ANS_EPI.SA+t2.ANS_EPI.JV+
          t2.ANS_EPA.AD+t2.ANS_EPA.SA+t2.ANS_EPA.JV+
          t2.ANS_ANS.AD+t2.ANS_ANS.SA+t2.ANS_ANS.JV+
          t2.ANS_END.AD+t2.ANS_END.SA+t2.ANS_END.JV+
          t2.ANS_LIECX.AD+t2.ANS_LIECX.SA+t2.ANS_LIECX.JV,5)==1
  
  
}


#les END identifi?s : que sont ils r?ellement ?
{
  
  #taux id EPI SA
  t2.END_EPI.SA=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("END_EPI") & tabNAT2015$Stade%in%"SA"])/
    sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  #taux id EPI AD
  t2.END_EPI.AD=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("END_EPI") & tabNAT2015$Stade%in%"AD"])/
    sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  #taux id EPI JV
  t2.END_EPI.JV=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("END_EPI") & tabNAT2015$Stade%in%"JV"])/
    sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  
  
  
  #taux id EPA SA
  t2.END_EPA.SA=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("END_EPA") & tabNAT2015$Stade%in%"SA"])/
    sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  #taux id EPA AD
  t2.END_EPA.AD=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("END_EPA") & tabNAT2015$Stade%in%"AD"])/
    sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  #taux id EPA JV
  t2.END_EPA.JV=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("END_EPA") & tabNAT2015$Stade%in%"JV"])/
    sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  
  
  
  #taux id ANS SA
  t2.END_ANS.SA=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("END_ANS") & tabNAT2015$Stade%in%"SA"])/
    sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  #taux id ANS AD
  t2.END_ANS.AD=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("END_ANS") & tabNAT2015$Stade%in%"AD"])/
    sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  #taux id ANS JV
  t2.END_ANS.JV=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("END_ANS") & tabNAT2015$Stade%in%"JV"])/
    sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  
  
  
  #taux id END SA
  t2.END_END.SA=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("END_END") & tabNAT2015$Stade%in%"SA"])/
    sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  #taux id END AD
  t2.END_END.AD=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("END_END") & tabNAT2015$Stade%in%"AD"])/
    sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  #taux id END JV
  t2.END_END.JV=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("END_END") & tabNAT2015$Stade%in%"JV"])/
    sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  
  
  
  
  #taux id LIEC_X SA
  t2.END_LIECX.SA=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("END_LIEC_X") & tabNAT2015$Stade%in%"SA"])/
    sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  #taux id LIEC_X AD
  t2.END_LIECX.AD=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("END_LIEC_X") & tabNAT2015$Stade%in%"AD"])/
    sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  #taux id LIEC_X JV
  t2.END_LIECX.JV=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("END_LIEC_X") & tabNAT2015$Stade%in%"JV"])/
    sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  
  round(t2.END_EPI.AD+t2.END_EPI.SA+t2.END_EPI.JV+
          t2.END_EPA.AD+t2.END_EPA.SA+t2.END_EPA.JV+
          t2.END_ANS.AD+t2.END_ANS.SA+t2.END_ANS.JV+
          t2.END_END.AD+t2.END_END.SA+t2.END_END.JV+
          t2.END_LIECX.AD+t2.END_LIECX.SA+t2.END_LIECX.JV,5)==1
  
  
}


#les LIECX identifi?s : que sont ils r?ellement ?
{
  
  #taux id EPI SA
  t2.LIECX_EPI.SA=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("LIEC_X_EPI") & tabNAT2015$Stade%in%"SA"])/
    sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  #taux id EPI AD
  t2.LIECX_EPI.AD=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("LIEC_X_EPI") & tabNAT2015$Stade%in%"AD"])/
    sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  #taux id EPI JV
  t2.LIECX_EPI.JV=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("LIEC_X_EPI") & tabNAT2015$Stade%in%"JV"])/
    sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  
  
  
  #taux id EPA SA
  t2.LIECX_EPA.SA=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("LIEC_X_EPA") & tabNAT2015$Stade%in%"SA"])/
    sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  #taux id EPA AD
  t2.LIECX_EPA.AD=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("LIEC_X_EPA") & tabNAT2015$Stade%in%"AD"])/
    sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  #taux id EPA JV
  t2.LIECX_EPA.JV=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("LIEC_X_EPA") & tabNAT2015$Stade%in%"JV"])/
    sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  
  
  
  #taux id ANS SA
  t2.LIECX_ANS.SA=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("LIEC_X_ANS") & tabNAT2015$Stade%in%"SA"])/
    sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  #taux id ANS AD
  t2.LIECX_ANS.AD=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("LIEC_X_ANS") & tabNAT2015$Stade%in%"AD"])/
    sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  #taux id ANS JV
  t2.LIECX_ANS.JV=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("LIEC_X_ANS") & tabNAT2015$Stade%in%"JV"])/
    sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  
  
  
  #taux id END SA
  t2.LIECX_END.SA=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("LIEC_X_END") & tabNAT2015$Stade%in%"SA"])/
    sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  #taux id END AD
  t2.LIECX_END.AD=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("LIEC_X_END") & tabNAT2015$Stade%in%"AD"])/
    sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  #taux id END JV
  t2.LIECX_END.JV=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("LIEC_X_END") & tabNAT2015$Stade%in%"JV"])/
    sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  
  
  
  
  #taux id LIEC_X SA
  t2.LIECX_LIECX.SA=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("LIEC_X_LIEC_X") & tabNAT2015$Stade%in%"SA"])/
    sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  #taux id LIEC_X AD
  t2.LIECX_LIECX.AD=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("LIEC_X_LIEC_X") & tabNAT2015$Stade%in%"AD"])/
    sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  #taux id LIEC_X JV
  t2.LIECX_LIECX.JV=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("LIEC_X_LIEC_X") & tabNAT2015$Stade%in%"JV"])/
    sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  
  round(t2.LIECX_EPI.AD+t2.LIECX_EPI.SA+t2.LIECX_EPI.JV+
          t2.LIECX_EPA.AD+t2.LIECX_EPA.SA+t2.LIECX_EPA.JV+
          t2.LIECX_ANS.AD+t2.LIECX_ANS.SA+t2.LIECX_ANS.JV+
          t2.LIECX_END.AD+t2.LIECX_END.SA+t2.LIECX_END.JV+
          t2.LIECX_LIECX.AD+t2.LIECX_LIECX.SA+t2.LIECX_LIECX.JV,5)==1
  
  
}


tabGraphNAT15_2=data.frame(
  Taux=c(
    t2.EPI_EPI.AD,t2.EPI_EPI.SA,t2.EPI_EPI.JV,
    t2.EPI_EPA.AD,t2.EPI_EPA.SA,t2.EPI_EPA.JV,
    t2.EPI_ANS.AD,t2.EPI_ANS.SA,t2.EPI_ANS.JV,
    t2.EPI_END.AD,t2.EPI_END.SA,t2.EPI_END.JV,
    t2.EPI_LIECX.AD,t2.EPI_LIECX.SA,t2.EPI_LIECX.JV,
    
    t2.EPA_EPI.AD,t2.EPA_EPI.SA,t2.EPA_EPI.JV,
    t2.EPA_EPA.AD,t2.EPA_EPA.SA,t2.EPA_EPA.JV,
    t2.EPA_ANS.AD,t2.EPA_ANS.SA,t2.EPA_ANS.JV,
    t2.EPA_END.AD,t2.EPA_END.SA,t2.EPA_END.JV,
    t2.EPA_LIECX.AD,t2.EPA_LIECX.SA,t2.EPA_LIECX.JV,
    
    t2.ANS_EPI.AD,t2.ANS_EPI.SA,t2.ANS_EPI.JV,
    t2.ANS_EPA.AD,t2.ANS_EPA.SA,t2.ANS_EPA.JV,
    t2.ANS_ANS.AD,t2.ANS_ANS.SA,t2.ANS_ANS.JV,
    t2.ANS_END.AD,t2.ANS_END.SA,t2.ANS_END.JV,
    t2.ANS_LIECX.AD,t2.ANS_LIECX.SA,t2.ANS_LIECX.JV,
    
    t2.END_EPI.AD,t2.END_EPI.SA,t2.END_EPI.JV,
    t2.END_EPA.AD,t2.END_EPA.SA,t2.END_EPA.JV,
    t2.END_ANS.AD,t2.END_ANS.SA,t2.END_ANS.JV,
    t2.END_END.AD,t2.END_END.SA,t2.END_END.JV,
    t2.END_LIECX.AD,t2.END_LIECX.SA,t2.END_LIECX.JV,
    
    t2.LIECX_EPI.AD,t2.LIECX_EPI.SA,t2.LIECX_EPI.JV,
    t2.LIECX_EPA.AD,t2.LIECX_EPA.SA,t2.LIECX_EPA.JV,
    t2.LIECX_ANS.AD,t2.LIECX_ANS.SA,t2.LIECX_ANS.JV,
    t2.LIECX_END.AD,t2.LIECX_END.SA,t2.LIECX_END.JV,
    t2.LIECX_LIECX.AD,t2.LIECX_LIECX.SA,t2.LIECX_LIECX.JV),
  
  Classe=c(
    "t2.EPI_EPI.AD","t2.EPI_EPI.SA","t2.EPI_EPI.JV",
    "t2.EPI_EPA.AD","t2.EPI_EPA.SA","t2.EPI_EPA.JV",
    "t2.EPI_ANS.AD","t2.EPI_ANS.SA","t2.EPI_ANS.JV",
    "t2.EPI_END.AD","t2.EPI_END.SA","t2.EPI_END.JV",
    "t2.EPI_LIECX.AD","t2.EPI_LIECX.SA","t2.EPI_LIECX.JV",
    
    "t2.EPA_EPI.AD","t2.EPA_EPI.SA","t2.EPA_EPI.JV",
    "t2.EPA_EPA.AD","t2.EPA_EPA.SA","t2.EPA_EPA.JV",
    "t2.EPA_ANS.AD","t2.EPA_ANS.SA","t2.EPA_ANS.JV",
    "t2.EPA_END.AD","t2.EPA_END.SA","t2.EPA_END.JV",
    "t2.EPA_LIECX.AD","t2.EPA_LIECX.SA","t2.EPA_LIECX.JV",
    
    "t2.ANS_EPI.AD","t2.ANS_EPI.SA","t2.ANS_EPI.JV",
    "t2.ANS_EPA.AD","t2.ANS_EPA.SA","t2.ANS_EPA.JV",
    "t2.ANS_ANS.AD","t2.ANS_ANS.SA","t2.ANS_ANS.JV",
    "t2.ANS_END.AD","t2.ANS_END.SA","t2.ANS_END.JV",
    "t2.ANS_LIECX.AD","t2.ANS_LIECX.SA","t2.ANS_LIECX.JV",
    
    "t2.END_EPI.AD","t2.END_EPI.SA","t2.END_EPI.JV",
    "t2.END_EPA.AD","t2.END_EPA.SA","t2.END_EPA.JV",
    "t2.END_ANS.AD","t2.END_ANS.SA","t2.END_ANS.JV",
    "t2.END_END.AD","t2.END_END.SA","t2.END_END.JV",
    "t2.END_LIECX.AD","t2.END_LIECX.SA","t2.END_LIECX.JV",
    
    "t2.LIECX_EPI.AD","t2.LIECX_EPI.SA","t2.LIECX_EPI.JV",
    "t2.LIECX_EPA.AD","t2.LIECX_EPA.SA","t2.LIECX_EPA.JV",
    "t2.LIECX_ANS.AD","t2.LIECX_ANS.SA","t2.LIECX_ANS.JV",
    "t2.LIECX_END.AD","t2.LIECX_END.SA","t2.LIECX_END.JV",
    "t2.LIECX_LIECX.AD","t2.LIECX_LIECX.SA","t2.LIECX_LIECX.JV")
  
)

tabGraphNAT15_2$LIEC=factor(substr(tabGraphNAT15_2$Classe,8,10),levels=c("EPI","EPA","ANS","END","LIECX"))
tabGraphNAT15_2$Stade=substr(tabGraphNAT15_2$Classe,12,13)
tabGraphNAT15_2$LIECid=factor(substr(tabGraphNAT15_2$Classe,4,6),levels=c("EPI","EPA","ANS","END","LIECX"))
tabGraphNAT15_2$LIEC_Stade=paste(tabGraphNAT15_2$LIEC,tabGraphNAT15_2$Stade,sep="_")

tabGraphNAT15_2$LIEC_Stade=factor(tabGraphNAT15_2$LIEC_Stade,levels=c(
  "EPI_AD","EPI_SA","EPI_JV","EPA_AD","EPA_SA",
  "EPA_JV","ANS_AD","ANS_SA","ANS_JV","END_AD",
  "END_SA","END_JV","LIECX_AD","LIECX_SA","LIECX_JV"))

ggplot(tabGraphNAT15_2, aes(fill=LIEC_Stade, y=Taux, x=LIECid)) + 
  geom_bar(position='stack', stat='identity')+
  theme_classic()+
  scale_fill_manual(values = c(
    rgb(0,0,153,maxColorValue = 255),    rgb(0,0,153,maxColorValue = 255),    rgb(0,0,153,maxColorValue = 255,alpha=150),
    rgb(255,255,0,maxColorValue = 255),  rgb(255,255,0,maxColorValue = 255),  rgb(255,255,0,maxColorValue = 255,alpha=100),
    rgb(0,153,0,maxColorValue = 255),    rgb(0,153,0,maxColorValue = 255),    rgb(0,153,0,maxColorValue = 255,alpha=150),
    rgb(255,102,153,maxColorValue = 255),rgb(255,102,153,maxColorValue = 255),rgb(255,102,153,maxColorValue = 255,alpha=150),
    rgb(123,123,123,maxColorValue = 255),rgb(123,123,123,maxColorValue = 255),rgb(123,123,123,maxColorValue = 255,alpha=150)
  ))

#zoom sur les erreurs
tabGraphNAT15_2$Erreur=NA
tabGraphNAT15_2$Erreur[tabGraphNAT15_2$LIEC==tabGraphNAT15_2$LIECid]="Correcte"
tabGraphNAT15_2$Erreur[tabGraphNAT15_2$LIEC!=tabGraphNAT15_2$LIECid]="Incorrecte"

ggplot(tabGraphNAT15_2[tabGraphNAT15_2$Erreur=="Incorrecte",], aes(fill=LIEC_Stade, y=Taux, x=LIECid)) + 
  geom_bar(position='fill', stat='identity')+
  theme_classic()+
  scale_fill_manual(values = c(
    rgb(0,0,153,maxColorValue = 255),    rgb(0,0,153,maxColorValue = 255,alpha=150),    rgb(0,0,153,maxColorValue = 255,alpha=100),
    rgb(255,255,0,maxColorValue = 255),  rgb(255,255,0,maxColorValue = 255,alpha=100),  rgb(255,255,0,maxColorValue = 255,alpha=25),
    rgb(0,153,0,maxColorValue = 255),    rgb(0,153,0,maxColorValue = 255,alpha=150),    rgb(0,153,0,maxColorValue = 255,alpha=100),
    rgb(255,102,153,maxColorValue = 255),rgb(255,102,153,maxColorValue = 255,alpha=150),rgb(255,102,153,maxColorValue = 255,alpha=100),
    rgb(123,123,123,maxColorValue = 255),rgb(123,123,123,maxColorValue = 255,alpha=150),rgb(123,123,123,maxColorValue = 255,alpha=100)
  ))




##### NAT 2014 ABSOLUE - Parmis ce que les gens ont class?s en EPI EPA ANS END : qu'est ce qu'on a r?ellement ? #####



nrow(Tab2)==21044
dataID <- unique(sbt_m[,c("cle",colID,"Region")])
nrow(dataID)==1109



tabNAT2014=Tab2[Tab2$Annee%in%2014,]
unique(tabNAT2014$Stade)#"JV" "SA" "AD" "X"
any(is.na(tabNAT2014$Stade))==F
nrow(tabNAT2014)==2938

#taux d'id correcte global
sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPI_EPI","EPA_EPA","ANS_ANS","END_END")])/
  sum(tabNAT2014$Nbr_EW)#0.6466317   => ok
#vectReussiteNAT

#les EPI identifi?s : que sont ils r?ellement ?
{
  
  #taux id correcte EPI SA
  ab.EPI_EPI.SA=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPI_EPI") & tabNAT2014$Stade%in%"SA"])
  
  #taux id correcte EPI AD
  ab.EPI_EPI.AD=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPI_EPI") & tabNAT2014$Stade%in%"AD"])
  
  #taux id correcte EPI JV
  ab.EPI_EPI.JV=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPI_EPI") & tabNAT2014$Stade%in%"JV"])
  
  
  
  
  #taux id incorrecte EPA SA
  ab.EPI_EPA.SA=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPI_EPA") & tabNAT2014$Stade%in%"SA"])
  
  #taux id incorrecte EPA AD
  ab.EPI_EPA.AD=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPI_EPA") & tabNAT2014$Stade%in%"AD"])
  
  #taux id incorrecte EPA JV
  ab.EPI_EPA.JV=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPI_EPA") & tabNAT2014$Stade%in%"JV"])
  
  
  
  
  #taux id incorrecte ANS SA
  ab.EPI_ANS.SA=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPI_ANS") & tabNAT2014$Stade%in%"SA"])
  
  #taux id incorrecte ANS AD
  ab.EPI_ANS.AD=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPI_ANS") & tabNAT2014$Stade%in%"AD"])
  
  #taux id incorrecte ANS JV
  ab.EPI_ANS.JV=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPI_ANS") & tabNAT2014$Stade%in%"JV"])
  
  
  
  
  #taux id incorrecte END SA
  ab.EPI_END.SA=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPI_END") & tabNAT2014$Stade%in%"SA"])
  
  #taux id incorrecte END AD
  ab.EPI_END.AD=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPI_END") & tabNAT2014$Stade%in%"AD"])
  
  #taux id incorrecte END JV
  ab.EPI_END.JV=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPI_END") & tabNAT2014$Stade%in%"JV"])
  
  
  
  
  
  #taux id incorrecte LIEC_X SA
  ab.EPI_LIECX.SA=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPI_LIEC_X") & tabNAT2014$Stade%in%"SA"])
  
  #taux id incorrecte LIEC_X AD
  ab.EPI_LIECX.AD=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPI_LIEC_X") & tabNAT2014$Stade%in%"AD"])
  
  #taux id incorrecte LIEC_X JV
  ab.EPI_LIECX.JV=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPI_LIEC_X") & tabNAT2014$Stade%in%"JV"])
  
  
  }


#les EPA identifi?s : que sont ils r?ellement ?
{
  
  #taux id EPI SA
  ab.EPA_EPI.SA=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPA_EPI") & tabNAT2014$Stade%in%"SA"])
  
  #taux id EPI AD
  ab.EPA_EPI.AD=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPA_EPI") & tabNAT2014$Stade%in%"AD"])
  
  #taux id EPI JV
  ab.EPA_EPI.JV=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPA_EPI") & tabNAT2014$Stade%in%"JV"])
  
  
  
  
  #taux id EPA SA
  ab.EPA_EPA.SA=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPA_EPA") & tabNAT2014$Stade%in%"SA"])
  
  #taux id EPA AD
  ab.EPA_EPA.AD=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPA_EPA") & tabNAT2014$Stade%in%"AD"])
  
  #taux id EPA JV
  ab.EPA_EPA.JV=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPA_EPA") & tabNAT2014$Stade%in%"JV"])
  
  
  
  
  #taux id ANS SA
  ab.EPA_ANS.SA=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPA_ANS") & tabNAT2014$Stade%in%"SA"])
  
  #taux id ANS AD
  ab.EPA_ANS.AD=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPA_ANS") & tabNAT2014$Stade%in%"AD"])
  
  #taux id ANS JV
  ab.EPA_ANS.JV=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPA_ANS") & tabNAT2014$Stade%in%"JV"])
  
  
  
  
  #taux id END SA
  ab.EPA_END.SA=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPA_END") & tabNAT2014$Stade%in%"SA"])
  
  #taux id END AD
  ab.EPA_END.AD=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPA_END") & tabNAT2014$Stade%in%"AD"])
  
  #taux id END JV
  ab.EPA_END.JV=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPA_END") & tabNAT2014$Stade%in%"JV"])
  
  
  
  
  
  #taux id LIEC_X SA
  ab.EPA_LIECX.SA=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPA_LIEC_X") & tabNAT2014$Stade%in%"SA"])
  
  #taux id LIEC_X AD
  ab.EPA_LIECX.AD=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPA_LIEC_X") & tabNAT2014$Stade%in%"AD"])
  
  #taux id LIEC_X JV
  ab.EPA_LIECX.JV=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("EPA_LIEC_X") & tabNAT2014$Stade%in%"JV"])
  
  
  
  
}


#les ANS identifi?s : que sont ils r?ellement ?
{
  
  #taux id EPI SA
  ab.ANS_EPI.SA=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("ANS_EPI") & tabNAT2014$Stade%in%"SA"])
  
  #taux id EPI AD
  ab.ANS_EPI.AD=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("ANS_EPI") & tabNAT2014$Stade%in%"AD"])
  
  #taux id EPI JV
  ab.ANS_EPI.JV=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("ANS_EPI") & tabNAT2014$Stade%in%"JV"])
  
  
  
  
  #taux id EPA SA
  ab.ANS_EPA.SA=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("ANS_EPA") & tabNAT2014$Stade%in%"SA"])
  
  #taux id EPA AD
  ab.ANS_EPA.AD=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("ANS_EPA") & tabNAT2014$Stade%in%"AD"])
  
  #taux id EPA JV
  ab.ANS_EPA.JV=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("ANS_EPA") & tabNAT2014$Stade%in%"JV"])
  
  
  
  
  #taux id ANS SA
  ab.ANS_ANS.SA=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("ANS_ANS") & tabNAT2014$Stade%in%"SA"])
  
  #taux id ANS AD
  ab.ANS_ANS.AD=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("ANS_ANS") & tabNAT2014$Stade%in%"AD"])
  
  #taux id ANS JV
  ab.ANS_ANS.JV=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("ANS_ANS") & tabNAT2014$Stade%in%"JV"])
  
  
  
  
  #taux id END SA
  ab.ANS_END.SA=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("ANS_END") & tabNAT2014$Stade%in%"SA"])
  
  #taux id END AD
  ab.ANS_END.AD=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("ANS_END") & tabNAT2014$Stade%in%"AD"])
  
  #taux id END JV
  ab.ANS_END.JV=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("ANS_END") & tabNAT2014$Stade%in%"JV"])
  
  
  
  
  
  #taux id LIEC_X SA
  ab.ANS_LIECX.SA=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("ANS_LIEC_X") & tabNAT2014$Stade%in%"SA"])
  
  #taux id LIEC_X AD
  ab.ANS_LIECX.AD=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("ANS_LIEC_X") & tabNAT2014$Stade%in%"AD"])
  
  #taux id LIEC_X JV
  ab.ANS_LIECX.JV=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("ANS_LIEC_X") & tabNAT2014$Stade%in%"JV"])
  
  
  
  
}


#les END identifi?s : que sont ils r?ellement ?
{
  
  #taux id EPI SA
  ab.END_EPI.SA=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("END_EPI") & tabNAT2014$Stade%in%"SA"])
  
  #taux id EPI AD
  ab.END_EPI.AD=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("END_EPI") & tabNAT2014$Stade%in%"AD"])
  
  #taux id EPI JV
  ab.END_EPI.JV=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("END_EPI") & tabNAT2014$Stade%in%"JV"])
  
  
  
  
  #taux id EPA SA
  ab.END_EPA.SA=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("END_EPA") & tabNAT2014$Stade%in%"SA"])
  
  #taux id EPA AD
  ab.END_EPA.AD=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("END_EPA") & tabNAT2014$Stade%in%"AD"])
  
  #taux id EPA JV
  ab.END_EPA.JV=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("END_EPA") & tabNAT2014$Stade%in%"JV"])
  
  
  
  
  #taux id ANS SA
  ab.END_ANS.SA=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("END_ANS") & tabNAT2014$Stade%in%"SA"])
  
  #taux id ANS AD
  ab.END_ANS.AD=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("END_ANS") & tabNAT2014$Stade%in%"AD"])
  
  #taux id ANS JV
  ab.END_ANS.JV=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("END_ANS") & tabNAT2014$Stade%in%"JV"])
  
  
  
  
  #taux id END SA
  ab.END_END.SA=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("END_END") & tabNAT2014$Stade%in%"SA"])
  
  #taux id END AD
  ab.END_END.AD=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("END_END") & tabNAT2014$Stade%in%"AD"])
  
  #taux id END JV
  ab.END_END.JV=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("END_END") & tabNAT2014$Stade%in%"JV"])
  
  
  
  
  #taux id LIEC_X SA
  ab.END_LIECX.SA=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("END_LIEC_X") & tabNAT2014$Stade%in%"SA"])
  
  #taux id LIEC_X AD
  ab.END_LIECX.AD=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("END_LIEC_X") & tabNAT2014$Stade%in%"AD"])
  
  #taux id LIEC_X JV
  ab.END_LIECX.JV=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("END_LIEC_X") & tabNAT2014$Stade%in%"JV"])
  
  
  
  
}


#les LIECX identifi?s : que sont ils r?ellement ?
{
  
  #taux id EPI SA
  ab.LIECX_EPI.SA=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("LIEC_X_EPI") & tabNAT2014$Stade%in%"SA"])
  
  #taux id EPI AD
  ab.LIECX_EPI.AD=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("LIEC_X_EPI") & tabNAT2014$Stade%in%"AD"])
  
  #taux id EPI JV
  ab.LIECX_EPI.JV=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("LIEC_X_EPI") & tabNAT2014$Stade%in%"JV"])
  
  
  
  
  #taux id EPA SA
  ab.LIECX_EPA.SA=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("LIEC_X_EPA") & tabNAT2014$Stade%in%"SA"])
  
  #taux id EPA AD
  ab.LIECX_EPA.AD=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("LIEC_X_EPA") & tabNAT2014$Stade%in%"AD"])
  
  #taux id EPA JV
  ab.LIECX_EPA.JV=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("LIEC_X_EPA") & tabNAT2014$Stade%in%"JV"])
  
  
  
  
  #taux id ANS SA
  ab.LIECX_ANS.SA=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("LIEC_X_ANS") & tabNAT2014$Stade%in%"SA"])
  
  #taux id ANS AD
  ab.LIECX_ANS.AD=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("LIEC_X_ANS") & tabNAT2014$Stade%in%"AD"])
  
  #taux id ANS JV
  ab.LIECX_ANS.JV=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("LIEC_X_ANS") & tabNAT2014$Stade%in%"JV"])
  
  
  
  
  #taux id END SA
  ab.LIECX_END.SA=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("LIEC_X_END") & tabNAT2014$Stade%in%"SA"])
  
  #taux id END AD
  ab.LIECX_END.AD=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("LIEC_X_END") & tabNAT2014$Stade%in%"AD"])
  
  #taux id END JV
  ab.LIECX_END.JV=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("LIEC_X_END") & tabNAT2014$Stade%in%"JV"])
  
  
  
  
  
  #taux id LIEC_X SA
  ab.LIECX_LIECX.SA=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("LIEC_X_LIEC_X") & tabNAT2014$Stade%in%"SA"])
  
  #taux id LIEC_X AD
  ab.LIECX_LIECX.AD=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("LIEC_X_LIEC_X") & tabNAT2014$Stade%in%"AD"])
  
  #taux id LIEC_X JV
  ab.LIECX_LIECX.JV=sum(tabNAT2014$Nbr_EW[tabNAT2014$FIEC_LIEC%in%c("LIEC_X_LIEC_X") & tabNAT2014$Stade%in%"JV"])
  
  
  
  
}


tabGraphNAT14_2=data.frame(
  Taux=c(
    ab.EPI_EPI.AD,ab.EPI_EPI.SA,ab.EPI_EPI.JV,
    ab.EPI_EPA.AD,ab.EPI_EPA.SA,ab.EPI_EPA.JV,
    ab.EPI_ANS.AD,ab.EPI_ANS.SA,ab.EPI_ANS.JV,
    ab.EPI_END.AD,ab.EPI_END.SA,ab.EPI_END.JV,
    ab.EPI_LIECX.AD,ab.EPI_LIECX.SA,ab.EPI_LIECX.JV,
    
    ab.EPA_EPI.AD,ab.EPA_EPI.SA,ab.EPA_EPI.JV,
    ab.EPA_EPA.AD,ab.EPA_EPA.SA,ab.EPA_EPA.JV,
    ab.EPA_ANS.AD,ab.EPA_ANS.SA,ab.EPA_ANS.JV,
    ab.EPA_END.AD,ab.EPA_END.SA,ab.EPA_END.JV,
    ab.EPA_LIECX.AD,ab.EPA_LIECX.SA,ab.EPA_LIECX.JV,
    
    ab.ANS_EPI.AD,ab.ANS_EPI.SA,ab.ANS_EPI.JV,
    ab.ANS_EPA.AD,ab.ANS_EPA.SA,ab.ANS_EPA.JV,
    ab.ANS_ANS.AD,ab.ANS_ANS.SA,ab.ANS_ANS.JV,
    ab.ANS_END.AD,ab.ANS_END.SA,ab.ANS_END.JV,
    ab.ANS_LIECX.AD,ab.ANS_LIECX.SA,ab.ANS_LIECX.JV,
    
    ab.END_EPI.AD,ab.END_EPI.SA,ab.END_EPI.JV,
    ab.END_EPA.AD,ab.END_EPA.SA,ab.END_EPA.JV,
    ab.END_ANS.AD,ab.END_ANS.SA,ab.END_ANS.JV,
    ab.END_END.AD,ab.END_END.SA,ab.END_END.JV,
    ab.END_LIECX.AD,ab.END_LIECX.SA,ab.END_LIECX.JV,
    
    ab.LIECX_EPI.AD,ab.LIECX_EPI.SA,ab.LIECX_EPI.JV,
    ab.LIECX_EPA.AD,ab.LIECX_EPA.SA,ab.LIECX_EPA.JV,
    ab.LIECX_ANS.AD,ab.LIECX_ANS.SA,ab.LIECX_ANS.JV,
    ab.LIECX_END.AD,ab.LIECX_END.SA,ab.LIECX_END.JV,
    ab.LIECX_LIECX.AD,ab.LIECX_LIECX.SA,ab.LIECX_LIECX.JV),
  
  Classe=c(
    "ab.EPI_EPI.AD","ab.EPI_EPI.SA","ab.EPI_EPI.JV",
    "ab.EPI_EPA.AD","ab.EPI_EPA.SA","ab.EPI_EPA.JV",
    "ab.EPI_ANS.AD","ab.EPI_ANS.SA","ab.EPI_ANS.JV",
    "ab.EPI_END.AD","ab.EPI_END.SA","ab.EPI_END.JV",
    "ab.EPI_LIECX.AD","ab.EPI_LIECX.SA","ab.EPI_LIECX.JV",
    
    "ab.EPA_EPI.AD","ab.EPA_EPI.SA","ab.EPA_EPI.JV",
    "ab.EPA_EPA.AD","ab.EPA_EPA.SA","ab.EPA_EPA.JV",
    "ab.EPA_ANS.AD","ab.EPA_ANS.SA","ab.EPA_ANS.JV",
    "ab.EPA_END.AD","ab.EPA_END.SA","ab.EPA_END.JV",
    "ab.EPA_LIECX.AD","ab.EPA_LIECX.SA","ab.EPA_LIECX.JV",
    
    "ab.ANS_EPI.AD","ab.ANS_EPI.SA","ab.ANS_EPI.JV",
    "ab.ANS_EPA.AD","ab.ANS_EPA.SA","ab.ANS_EPA.JV",
    "ab.ANS_ANS.AD","ab.ANS_ANS.SA","ab.ANS_ANS.JV",
    "ab.ANS_END.AD","ab.ANS_END.SA","ab.ANS_END.JV",
    "ab.ANS_LIECX.AD","ab.ANS_LIECX.SA","ab.ANS_LIECX.JV",
    
    "ab.END_EPI.AD","ab.END_EPI.SA","ab.END_EPI.JV",
    "ab.END_EPA.AD","ab.END_EPA.SA","ab.END_EPA.JV",
    "ab.END_ANS.AD","ab.END_ANS.SA","ab.END_ANS.JV",
    "ab.END_END.AD","ab.END_END.SA","ab.END_END.JV",
    "ab.END_LIECX.AD","ab.END_LIECX.SA","ab.END_LIECX.JV",
    
    "ab.LIECX_EPI.AD","ab.LIECX_EPI.SA","ab.LIECX_EPI.JV",
    "ab.LIECX_EPA.AD","ab.LIECX_EPA.SA","ab.LIECX_EPA.JV",
    "ab.LIECX_ANS.AD","ab.LIECX_ANS.SA","ab.LIECX_ANS.JV",
    "ab.LIECX_END.AD","ab.LIECX_END.SA","ab.LIECX_END.JV",
    "ab.LIECX_LIECX.AD","ab.LIECX_LIECX.SA","ab.LIECX_LIECX.JV")
  
)

tabGraphNAT14_2$LIEC=factor(substr(tabGraphNAT14_2$Classe,8,10),levels=c("EPI","EPA","ANS","END","LIECX"))
tabGraphNAT14_2$Stade=substr(tabGraphNAT14_2$Classe,12,13)
tabGraphNAT14_2$LIECid=factor(substr(tabGraphNAT14_2$Classe,4,6),levels=c("EPI","EPA","ANS","END","LIECX"))
tabGraphNAT14_2$LIEC_Stade=paste(tabGraphNAT14_2$LIEC,tabGraphNAT14_2$Stade,sep="_")

tabGraphNAT14_2$LIEC_Stade=factor(tabGraphNAT14_2$LIEC_Stade,levels=c(
  "EPI_AD","EPI_SA","EPI_JV","EPA_AD","EPA_SA",
  "EPA_JV","ANS_AD","ANS_SA","ANS_JV","END_AD",
  "END_SA","END_JV","LIECX_AD","LIECX_SA","LIECX_JV"))

ggplot(tabGraphNAT14_2, aes(fill=LIEC_Stade, y=Taux, x=LIECid)) + 
  geom_bar(position='fill', stat='identity')+
  theme_classic()+ylab("Nombre de VDT")+
  scale_fill_manual(values = c(
    rgb(0,0,153,maxColorValue = 255),    rgb(0,0,153,maxColorValue = 255),    rgb(0,0,153,maxColorValue = 255,alpha=150),
    rgb(255,255,0,maxColorValue = 255),  rgb(255,255,0,maxColorValue = 255),  rgb(255,255,0,maxColorValue = 255,alpha=100),
    rgb(0,153,0,maxColorValue = 255),    rgb(0,153,0,maxColorValue = 255),    rgb(0,153,0,maxColorValue = 255,alpha=150),
    rgb(255,102,153,maxColorValue = 255),rgb(255,102,153,maxColorValue = 255),rgb(255,102,153,maxColorValue = 255,alpha=150),
    rgb(123,123,123,maxColorValue = 255),rgb(123,123,123,maxColorValue = 255),rgb(123,123,123,maxColorValue = 255,alpha=150)
  ))

#zoom sur les erreurs
tabGraphNAT14_2$Erreur=NA
tabGraphNAT14_2$Erreur[tabGraphNAT14_2$LIEC==tabGraphNAT14_2$LIECid]="Correcte"
tabGraphNAT14_2$Erreur[tabGraphNAT14_2$LIEC!=tabGraphNAT14_2$LIECid]="Incorrecte"

ggplot(tabGraphNAT14_2[tabGraphNAT14_2$Erreur=="Incorrecte",], aes(fill=LIEC_Stade, y=Taux, x=LIECid)) + 
  geom_bar(position='fill', stat='identity')+
  theme_classic()+
  scale_fill_manual(values = c(
    rgb(0,0,153,maxColorValue = 255),    rgb(0,0,153,maxColorValue = 255,alpha=150),    rgb(0,0,153,maxColorValue = 255,alpha=100),
    rgb(255,255,0,maxColorValue = 255),  rgb(255,255,0,maxColorValue = 255,alpha=100),  rgb(255,255,0,maxColorValue = 255,alpha=25),
    rgb(0,153,0,maxColorValue = 255),    rgb(0,153,0,maxColorValue = 255,alpha=150),    rgb(0,153,0,maxColorValue = 255,alpha=100),
    rgb(255,102,153,maxColorValue = 255),rgb(255,102,153,maxColorValue = 255,alpha=150),rgb(255,102,153,maxColorValue = 255,alpha=100),
    rgb(123,123,123,maxColorValue = 255),rgb(123,123,123,maxColorValue = 255,alpha=150),rgb(123,123,123,maxColorValue = 255,alpha=100)
  ))




##### NAT 2015 ABSOLUE - Parmis ce que les gens ont class?s en EPI EPA ANS END : qu'est ce qu'on a r?ellement ? #####



nrow(Tab2)==21044
dataID <- unique(sbt_m[,c("cle",colID,"Region")])
nrow(dataID)==1109



tabNAT2015=Tab2[Tab2$Annee%in%2015,]
unique(tabNAT2015$Stade)#"JV" "SA" "AD" "X"
any(is.na(tabNAT2015$Stade))==F
nrow(tabNAT2015)==2751

#taux d'id correcte global
sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPI_EPI","EPA_EPA","ANS_ANS","END_END")])/
  sum(tabNAT2015$Nbr_EW)#0.8171038   => ok
#vectReussiteNAT

tabNAT2015=tabNAT2015[tabNAT2015$Stade!="X",]
nrow(tabNAT2015)==2751-11

#les EPI identifi?s : que sont ils r?ellement ?
{
  
  #taux id correcte EPI SA
  ab.EPI_EPI.SA=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPI_EPI") & tabNAT2015$Stade%in%"SA"])
  
  #taux id correcte EPI AD
  ab.EPI_EPI.AD=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPI_EPI") & tabNAT2015$Stade%in%"AD"])
  
  #taux id correcte EPI JV
  ab.EPI_EPI.JV=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPI_EPI") & tabNAT2015$Stade%in%"JV"])
  
  
  
  
  #taux id incorrecte EPA SA
  ab.EPI_EPA.SA=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPI_EPA") & tabNAT2015$Stade%in%"SA"])
  
  #taux id incorrecte EPA AD
  ab.EPI_EPA.AD=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPI_EPA") & tabNAT2015$Stade%in%"AD"])
  
  #taux id incorrecte EPA JV
  ab.EPI_EPA.JV=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPI_EPA") & tabNAT2015$Stade%in%"JV"])
  
  
  
  
  #taux id incorrecte ANS SA
  ab.EPI_ANS.SA=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPI_ANS") & tabNAT2015$Stade%in%"SA"])
  
  #taux id incorrecte ANS AD
  ab.EPI_ANS.AD=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPI_ANS") & tabNAT2015$Stade%in%"AD"])
  
  #taux id incorrecte ANS JV
  ab.EPI_ANS.JV=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPI_ANS") & tabNAT2015$Stade%in%"JV"])
  
  
  
  
  #taux id incorrecte END SA
  ab.EPI_END.SA=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPI_END") & tabNAT2015$Stade%in%"SA"])
  
  #taux id incorrecte END AD
  ab.EPI_END.AD=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPI_END") & tabNAT2015$Stade%in%"AD"])
  
  #taux id incorrecte END JV
  ab.EPI_END.JV=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPI_END") & tabNAT2015$Stade%in%"JV"])
  
  
  
  
  
  #taux id incorrecte LIEC_X SA
  ab.EPI_LIECX.SA=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPI_LIEC_X") & tabNAT2015$Stade%in%"SA"])
  
  #taux id incorrecte LIEC_X AD
  ab.EPI_LIECX.AD=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPI_LIEC_X") & tabNAT2015$Stade%in%"AD"])
  
  #taux id incorrecte LIEC_X JV
  ab.EPI_LIECX.JV=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPI_LIEC_X") & tabNAT2015$Stade%in%"JV"])
  
  
}


#les EPA identifi?s : que sont ils r?ellement ?
{
  
  #taux id EPI SA
  ab.EPA_EPI.SA=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPA_EPI") & tabNAT2015$Stade%in%"SA"])
  
  #taux id EPI AD
  ab.EPA_EPI.AD=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPA_EPI") & tabNAT2015$Stade%in%"AD"])
  
  #taux id EPI JV
  ab.EPA_EPI.JV=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPA_EPI") & tabNAT2015$Stade%in%"JV"])
  
  
  
  
  #taux id EPA SA
  ab.EPA_EPA.SA=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPA_EPA") & tabNAT2015$Stade%in%"SA"])
  
  #taux id EPA AD
  ab.EPA_EPA.AD=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPA_EPA") & tabNAT2015$Stade%in%"AD"])
  
  #taux id EPA JV
  ab.EPA_EPA.JV=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPA_EPA") & tabNAT2015$Stade%in%"JV"])
  
  
  
  
  #taux id ANS SA
  ab.EPA_ANS.SA=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPA_ANS") & tabNAT2015$Stade%in%"SA"])
  
  #taux id ANS AD
  ab.EPA_ANS.AD=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPA_ANS") & tabNAT2015$Stade%in%"AD"])
  
  #taux id ANS JV
  ab.EPA_ANS.JV=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPA_ANS") & tabNAT2015$Stade%in%"JV"])
  
  
  
  
  #taux id END SA
  ab.EPA_END.SA=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPA_END") & tabNAT2015$Stade%in%"SA"])
  
  #taux id END AD
  ab.EPA_END.AD=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPA_END") & tabNAT2015$Stade%in%"AD"])
  
  #taux id END JV
  ab.EPA_END.JV=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPA_END") & tabNAT2015$Stade%in%"JV"])
  
  
  
  
  
  #taux id LIEC_X SA
  ab.EPA_LIECX.SA=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPA_LIEC_X") & tabNAT2015$Stade%in%"SA"])
  
  #taux id LIEC_X AD
  ab.EPA_LIECX.AD=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPA_LIEC_X") & tabNAT2015$Stade%in%"AD"])
  
  #taux id LIEC_X JV
  ab.EPA_LIECX.JV=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("EPA_LIEC_X") & tabNAT2015$Stade%in%"JV"])
  
  
  
  
}


#les ANS identifi?s : que sont ils r?ellement ?
{
  
  #taux id EPI SA
  ab.ANS_EPI.SA=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("ANS_EPI") & tabNAT2015$Stade%in%"SA"])
  
  #taux id EPI AD
  ab.ANS_EPI.AD=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("ANS_EPI") & tabNAT2015$Stade%in%"AD"])
  
  #taux id EPI JV
  ab.ANS_EPI.JV=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("ANS_EPI") & tabNAT2015$Stade%in%"JV"])
  
  
  
  
  #taux id EPA SA
  ab.ANS_EPA.SA=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("ANS_EPA") & tabNAT2015$Stade%in%"SA"])
  
  #taux id EPA AD
  ab.ANS_EPA.AD=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("ANS_EPA") & tabNAT2015$Stade%in%"AD"])
  
  #taux id EPA JV
  ab.ANS_EPA.JV=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("ANS_EPA") & tabNAT2015$Stade%in%"JV"])
  
  
  
  
  #taux id ANS SA
  ab.ANS_ANS.SA=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("ANS_ANS") & tabNAT2015$Stade%in%"SA"])
  
  #taux id ANS AD
  ab.ANS_ANS.AD=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("ANS_ANS") & tabNAT2015$Stade%in%"AD"])
  
  #taux id ANS JV
  ab.ANS_ANS.JV=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("ANS_ANS") & tabNAT2015$Stade%in%"JV"])
  
  
  
  
  #taux id END SA
  ab.ANS_END.SA=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("ANS_END") & tabNAT2015$Stade%in%"SA"])
  
  #taux id END AD
  ab.ANS_END.AD=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("ANS_END") & tabNAT2015$Stade%in%"AD"])
  
  #taux id END JV
  ab.ANS_END.JV=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("ANS_END") & tabNAT2015$Stade%in%"JV"])
  
  
  
  
  
  #taux id LIEC_X SA
  ab.ANS_LIECX.SA=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("ANS_LIEC_X") & tabNAT2015$Stade%in%"SA"])
  
  #taux id LIEC_X AD
  ab.ANS_LIECX.AD=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("ANS_LIEC_X") & tabNAT2015$Stade%in%"AD"])
  
  #taux id LIEC_X JV
  ab.ANS_LIECX.JV=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("ANS_LIEC_X") & tabNAT2015$Stade%in%"JV"])
  
  
  
  
}


#les END identifi?s : que sont ils r?ellement ?
{
  
  #taux id EPI SA
  ab.END_EPI.SA=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("END_EPI") & tabNAT2015$Stade%in%"SA"])
  
  #taux id EPI AD
  ab.END_EPI.AD=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("END_EPI") & tabNAT2015$Stade%in%"AD"])
  
  #taux id EPI JV
  ab.END_EPI.JV=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("END_EPI") & tabNAT2015$Stade%in%"JV"])
  
  
  
  
  #taux id EPA SA
  ab.END_EPA.SA=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("END_EPA") & tabNAT2015$Stade%in%"SA"])
  
  #taux id EPA AD
  ab.END_EPA.AD=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("END_EPA") & tabNAT2015$Stade%in%"AD"])
  
  #taux id EPA JV
  ab.END_EPA.JV=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("END_EPA") & tabNAT2015$Stade%in%"JV"])
  
  
  
  
  #taux id ANS SA
  ab.END_ANS.SA=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("END_ANS") & tabNAT2015$Stade%in%"SA"])
  
  #taux id ANS AD
  ab.END_ANS.AD=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("END_ANS") & tabNAT2015$Stade%in%"AD"])
  
  #taux id ANS JV
  ab.END_ANS.JV=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("END_ANS") & tabNAT2015$Stade%in%"JV"])
  
  
  
  
  #taux id END SA
  ab.END_END.SA=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("END_END") & tabNAT2015$Stade%in%"SA"])
  
  #taux id END AD
  ab.END_END.AD=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("END_END") & tabNAT2015$Stade%in%"AD"])
  
  #taux id END JV
  ab.END_END.JV=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("END_END") & tabNAT2015$Stade%in%"JV"])
  
  
  
  
  #taux id LIEC_X SA
  ab.END_LIECX.SA=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("END_LIEC_X") & tabNAT2015$Stade%in%"SA"])
  
  #taux id LIEC_X AD
  ab.END_LIECX.AD=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("END_LIEC_X") & tabNAT2015$Stade%in%"AD"])
  
  #taux id LIEC_X JV
  ab.END_LIECX.JV=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("END_LIEC_X") & tabNAT2015$Stade%in%"JV"])
  
  
  
  
}


#les LIECX identifi?s : que sont ils r?ellement ?
{
  
  #taux id EPI SA
  ab.LIECX_EPI.SA=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("LIEC_X_EPI") & tabNAT2015$Stade%in%"SA"])
  
  #taux id EPI AD
  ab.LIECX_EPI.AD=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("LIEC_X_EPI") & tabNAT2015$Stade%in%"AD"])
  
  #taux id EPI JV
  ab.LIECX_EPI.JV=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("LIEC_X_EPI") & tabNAT2015$Stade%in%"JV"])
  
  
  
  
  #taux id EPA SA
  ab.LIECX_EPA.SA=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("LIEC_X_EPA") & tabNAT2015$Stade%in%"SA"])
  
  #taux id EPA AD
  ab.LIECX_EPA.AD=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("LIEC_X_EPA") & tabNAT2015$Stade%in%"AD"])
  
  #taux id EPA JV
  ab.LIECX_EPA.JV=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("LIEC_X_EPA") & tabNAT2015$Stade%in%"JV"])
  
  
  
  
  #taux id ANS SA
  ab.LIECX_ANS.SA=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("LIEC_X_ANS") & tabNAT2015$Stade%in%"SA"])
  
  #taux id ANS AD
  ab.LIECX_ANS.AD=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("LIEC_X_ANS") & tabNAT2015$Stade%in%"AD"])
  
  #taux id ANS JV
  ab.LIECX_ANS.JV=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("LIEC_X_ANS") & tabNAT2015$Stade%in%"JV"])
  
  
  
  
  #taux id END SA
  ab.LIECX_END.SA=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("LIEC_X_END") & tabNAT2015$Stade%in%"SA"])
  
  #taux id END AD
  ab.LIECX_END.AD=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("LIEC_X_END") & tabNAT2015$Stade%in%"AD"])
  
  #taux id END JV
  ab.LIECX_END.JV=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("LIEC_X_END") & tabNAT2015$Stade%in%"JV"])
  
  
  
  
  
  #taux id LIEC_X SA
  ab.LIECX_LIECX.SA=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("LIEC_X_LIEC_X") & tabNAT2015$Stade%in%"SA"])
  
  #taux id LIEC_X AD
  ab.LIECX_LIECX.AD=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("LIEC_X_LIEC_X") & tabNAT2015$Stade%in%"AD"])
  
  #taux id LIEC_X JV
  ab.LIECX_LIECX.JV=sum(tabNAT2015$Nbr_EW[tabNAT2015$FIEC_LIEC%in%c("LIEC_X_LIEC_X") & tabNAT2015$Stade%in%"JV"])
  
  
  
  
}


tabGraphNAT15_2=data.frame(
  Taux=c(
    ab.EPI_EPI.AD,ab.EPI_EPI.SA,ab.EPI_EPI.JV,
    ab.EPI_EPA.AD,ab.EPI_EPA.SA,ab.EPI_EPA.JV,
    ab.EPI_ANS.AD,ab.EPI_ANS.SA,ab.EPI_ANS.JV,
    ab.EPI_END.AD,ab.EPI_END.SA,ab.EPI_END.JV,
    ab.EPI_LIECX.AD,ab.EPI_LIECX.SA,ab.EPI_LIECX.JV,
    
    ab.EPA_EPI.AD,ab.EPA_EPI.SA,ab.EPA_EPI.JV,
    ab.EPA_EPA.AD,ab.EPA_EPA.SA,ab.EPA_EPA.JV,
    ab.EPA_ANS.AD,ab.EPA_ANS.SA,ab.EPA_ANS.JV,
    ab.EPA_END.AD,ab.EPA_END.SA,ab.EPA_END.JV,
    ab.EPA_LIECX.AD,ab.EPA_LIECX.SA,ab.EPA_LIECX.JV,
    
    ab.ANS_EPI.AD,ab.ANS_EPI.SA,ab.ANS_EPI.JV,
    ab.ANS_EPA.AD,ab.ANS_EPA.SA,ab.ANS_EPA.JV,
    ab.ANS_ANS.AD,ab.ANS_ANS.SA,ab.ANS_ANS.JV,
    ab.ANS_END.AD,ab.ANS_END.SA,ab.ANS_END.JV,
    ab.ANS_LIECX.AD,ab.ANS_LIECX.SA,ab.ANS_LIECX.JV,
    
    ab.END_EPI.AD,ab.END_EPI.SA,ab.END_EPI.JV,
    ab.END_EPA.AD,ab.END_EPA.SA,ab.END_EPA.JV,
    ab.END_ANS.AD,ab.END_ANS.SA,ab.END_ANS.JV,
    ab.END_END.AD,ab.END_END.SA,ab.END_END.JV,
    ab.END_LIECX.AD,ab.END_LIECX.SA,ab.END_LIECX.JV,
    
    ab.LIECX_EPI.AD,ab.LIECX_EPI.SA,ab.LIECX_EPI.JV,
    ab.LIECX_EPA.AD,ab.LIECX_EPA.SA,ab.LIECX_EPA.JV,
    ab.LIECX_ANS.AD,ab.LIECX_ANS.SA,ab.LIECX_ANS.JV,
    ab.LIECX_END.AD,ab.LIECX_END.SA,ab.LIECX_END.JV,
    ab.LIECX_LIECX.AD,ab.LIECX_LIECX.SA,ab.LIECX_LIECX.JV),
  
  Classe=c(
    "ab.EPI_EPI.AD","ab.EPI_EPI.SA","ab.EPI_EPI.JV",
    "ab.EPI_EPA.AD","ab.EPI_EPA.SA","ab.EPI_EPA.JV",
    "ab.EPI_ANS.AD","ab.EPI_ANS.SA","ab.EPI_ANS.JV",
    "ab.EPI_END.AD","ab.EPI_END.SA","ab.EPI_END.JV",
    "ab.EPI_LIECX.AD","ab.EPI_LIECX.SA","ab.EPI_LIECX.JV",
    
    "ab.EPA_EPI.AD","ab.EPA_EPI.SA","ab.EPA_EPI.JV",
    "ab.EPA_EPA.AD","ab.EPA_EPA.SA","ab.EPA_EPA.JV",
    "ab.EPA_ANS.AD","ab.EPA_ANS.SA","ab.EPA_ANS.JV",
    "ab.EPA_END.AD","ab.EPA_END.SA","ab.EPA_END.JV",
    "ab.EPA_LIECX.AD","ab.EPA_LIECX.SA","ab.EPA_LIECX.JV",
    
    "ab.ANS_EPI.AD","ab.ANS_EPI.SA","ab.ANS_EPI.JV",
    "ab.ANS_EPA.AD","ab.ANS_EPA.SA","ab.ANS_EPA.JV",
    "ab.ANS_ANS.AD","ab.ANS_ANS.SA","ab.ANS_ANS.JV",
    "ab.ANS_END.AD","ab.ANS_END.SA","ab.ANS_END.JV",
    "ab.ANS_LIECX.AD","ab.ANS_LIECX.SA","ab.ANS_LIECX.JV",
    
    "ab.END_EPI.AD","ab.END_EPI.SA","ab.END_EPI.JV",
    "ab.END_EPA.AD","ab.END_EPA.SA","ab.END_EPA.JV",
    "ab.END_ANS.AD","ab.END_ANS.SA","ab.END_ANS.JV",
    "ab.END_END.AD","ab.END_END.SA","ab.END_END.JV",
    "ab.END_LIECX.AD","ab.END_LIECX.SA","ab.END_LIECX.JV",
    
    "ab.LIECX_EPI.AD","ab.LIECX_EPI.SA","ab.LIECX_EPI.JV",
    "ab.LIECX_EPA.AD","ab.LIECX_EPA.SA","ab.LIECX_EPA.JV",
    "ab.LIECX_ANS.AD","ab.LIECX_ANS.SA","ab.LIECX_ANS.JV",
    "ab.LIECX_END.AD","ab.LIECX_END.SA","ab.LIECX_END.JV",
    "ab.LIECX_LIECX.AD","ab.LIECX_LIECX.SA","ab.LIECX_LIECX.JV")
  
)

tabGraphNAT15_2$LIEC=factor(substr(tabGraphNAT15_2$Classe,8,10),levels=c("EPI","EPA","ANS","END","LIECX"))
tabGraphNAT15_2$Stade=substr(tabGraphNAT15_2$Classe,12,13)
tabGraphNAT15_2$LIECid=factor(substr(tabGraphNAT15_2$Classe,4,6),levels=c("EPI","EPA","ANS","END","LIECX"))
tabGraphNAT15_2$LIEC_Stade=paste(tabGraphNAT15_2$LIEC,tabGraphNAT15_2$Stade,sep="_")

tabGraphNAT15_2$LIEC_Stade=factor(tabGraphNAT15_2$LIEC_Stade,levels=c(
  "EPI_AD","EPI_SA","EPI_JV","EPA_AD","EPA_SA",
  "EPA_JV","ANS_AD","ANS_SA","ANS_JV","END_AD",
  "END_SA","END_JV","LIECX_AD","LIECX_SA","LIECX_JV"))

ggplot(tabGraphNAT15_2, aes(fill=LIEC_Stade, y=Taux, x=LIECid)) + 
  geom_bar(position='stack', stat='identity')+
  theme_classic()+ylab("Nombre de VDT")+
  scale_fill_manual(values = c(
    rgb(0,0,153,maxColorValue = 255),    rgb(0,0,153,maxColorValue = 255),    rgb(0,0,153,maxColorValue = 255,alpha=150),
    rgb(255,255,0,maxColorValue = 255),  rgb(255,255,0,maxColorValue = 255),  rgb(255,255,0,maxColorValue = 255,alpha=100),
    rgb(0,153,0,maxColorValue = 255),    rgb(0,153,0,maxColorValue = 255),    rgb(0,153,0,maxColorValue = 255,alpha=150),
    rgb(255,102,153,maxColorValue = 255),rgb(255,102,153,maxColorValue = 255),rgb(255,102,153,maxColorValue = 255,alpha=150),
    rgb(123,123,123,maxColorValue = 255),rgb(123,123,123,maxColorValue = 255),rgb(123,123,123,maxColorValue = 255,alpha=150)
  ))

#zoom sur les erreurs
tabGraphNAT15_2$Erreur=NA
tabGraphNAT15_2$Erreur[tabGraphNAT15_2$LIEC==tabGraphNAT15_2$LIECid]="Correcte"
tabGraphNAT15_2$Erreur[tabGraphNAT15_2$LIEC!=tabGraphNAT15_2$LIECid]="Incorrecte"

ggplot(tabGraphNAT15_2[tabGraphNAT15_2$Erreur=="Incorrecte",], aes(fill=LIEC_Stade, y=Taux, x=LIECid)) + 
  geom_bar(position='fill', stat='identity')+
  theme_classic()+
  scale_fill_manual(values = c(
    rgb(0,0,153,maxColorValue = 255),    rgb(0,0,153,maxColorValue = 255,alpha=150),    rgb(0,0,153,maxColorValue = 255,alpha=100),
    rgb(255,255,0,maxColorValue = 255),  rgb(255,255,0,maxColorValue = 255,alpha=100),  rgb(255,255,0,maxColorValue = 255,alpha=25),
    rgb(0,153,0,maxColorValue = 255),    rgb(0,153,0,maxColorValue = 255,alpha=150),    rgb(0,153,0,maxColorValue = 255,alpha=100),
    rgb(255,102,153,maxColorValue = 255),rgb(255,102,153,maxColorValue = 255,alpha=150),rgb(255,102,153,maxColorValue = 255,alpha=100),
    rgb(123,123,123,maxColorValue = 255),rgb(123,123,123,maxColorValue = 255,alpha=150),rgb(123,123,123,maxColorValue = 255,alpha=100)
  ))




##### Graphique BO vs NAT #####



plot(2013:2018,vectReussiteNAT*100,col="red",type="l",lwd=4,ylim=c(50,100),
     main="Taux d'identification terrain r?ussi Aquitaine",
     xlab="",
     ylab="",
     yaxt='n')
lines(2013:2018,vectReussiteBO*100,col="black",type="l",lwd=4)
axis(side = 2, at=c(0,50,100),las=1,labels=c("0%","50%","100%"))
axis(side = 1, at=2013:2018, pos=-15, tick=F, labels=c(paste("n=",n2013,sep=""),
                                                       paste("n=",n2014,sep=""),paste("n=",n2015,sep=""),
                                                       paste("n=",n2016,sep=""),paste("n=",n2017,sep=""),
                                                       paste("n=",n2018,sep="")))



plot(2013:2018,rep(mean(vectReussiteNAT)*100,6),col="red",type="l",lwd=4,ylim=c(50,100),
     main="Taux d'identification terrain r?ussi Aquitaine",
     xlab="",
     ylab="",
     yaxt='n')
lines(2013:2018,vectReussiteBO*100,col="black",type="l",lwd=4)
axis(side = 2, at=c(0,50,100),las=1,labels=c("0%","50%","100%"))
axis(side = 1, at=2013:2018, pos=-15, tick=F, labels=c(paste("n=",n2013,sep=""),
                                                       paste("n=",n2014,sep=""),paste("n=",n2015,sep=""),
                                                       paste("n=",n2016,sep=""),paste("n=",n2017,sep=""),
                                                       paste("n=",n2018,sep="")))


##### Boxplot intro #####

vAnnee=c()
vRegion=c()
vCombi=c()
vAbs=c()
vtauxErr=c()


for(i in unique(Tab2$Region)){
  dat_i=Tab2[Tab2$Region%in%i,]
  for(j in unique(dat_i$Annee)){
    dat_ij=dat_i[dat_i$Annee==j,]
    vAnnee=c(vAnnee,rep(j,25))
    vRegion=c(vRegion,rep(i,25))
    vtauxErr=c(vtauxErr,
               sum(dat_ij$Nbr_EW[dat_ij$FIEC_LIEC=="EPI_EPI"])/sum(dat_ij$Nbr_EW),
               sum(dat_ij$Nbr_EW[dat_ij$FIEC_LIEC=="EPA_EPA"])/sum(dat_ij$Nbr_EW),
               sum(dat_ij$Nbr_EW[dat_ij$FIEC_LIEC=="ANS_ANS"])/sum(dat_ij$Nbr_EW),
               sum(dat_ij$Nbr_EW[dat_ij$FIEC_LIEC=="END_END"])/sum(dat_ij$Nbr_EW),
               sum(dat_ij$Nbr_EW[dat_ij$FIEC_LIEC=="LIEC_X_LIEC_X"])/sum(dat_ij$Nbr_EW),
               sum(dat_ij$Nbr_EW[dat_ij$FIEC_LIEC=="EPI_EPA"])/sum(dat_ij$Nbr_EW),
               sum(dat_ij$Nbr_EW[dat_ij$FIEC_LIEC=="EPI_ANS"])/sum(dat_ij$Nbr_EW),
               sum(dat_ij$Nbr_EW[dat_ij$FIEC_LIEC=="EPI_END"])/sum(dat_ij$Nbr_EW),
               sum(dat_ij$Nbr_EW[dat_ij$FIEC_LIEC=="EPI_LIEC_X"])/sum(dat_ij$Nbr_EW),
               sum(dat_ij$Nbr_EW[dat_ij$FIEC_LIEC=="EPA_EPI"])/sum(dat_ij$Nbr_EW),
               sum(dat_ij$Nbr_EW[dat_ij$FIEC_LIEC=="EPA_ANS"])/sum(dat_ij$Nbr_EW),
               sum(dat_ij$Nbr_EW[dat_ij$FIEC_LIEC=="EPA_END"])/sum(dat_ij$Nbr_EW),
               sum(dat_ij$Nbr_EW[dat_ij$FIEC_LIEC=="EPA_LIEC_X"])/sum(dat_ij$Nbr_EW),
               sum(dat_ij$Nbr_EW[dat_ij$FIEC_LIEC=="ANS_EPI"])/sum(dat_ij$Nbr_EW),
               sum(dat_ij$Nbr_EW[dat_ij$FIEC_LIEC=="ANS_EPA"])/sum(dat_ij$Nbr_EW),
               sum(dat_ij$Nbr_EW[dat_ij$FIEC_LIEC=="ANS_END"])/sum(dat_ij$Nbr_EW),
               sum(dat_ij$Nbr_EW[dat_ij$FIEC_LIEC=="ANS_LIEC_X"])/sum(dat_ij$Nbr_EW),
               sum(dat_ij$Nbr_EW[dat_ij$FIEC_LIEC=="END_EPI"])/sum(dat_ij$Nbr_EW),
               sum(dat_ij$Nbr_EW[dat_ij$FIEC_LIEC=="END_EPA"])/sum(dat_ij$Nbr_EW),
               sum(dat_ij$Nbr_EW[dat_ij$FIEC_LIEC=="END_ANS"])/sum(dat_ij$Nbr_EW),
               sum(dat_ij$Nbr_EW[dat_ij$FIEC_LIEC=="END_LIEC_X"])/sum(dat_ij$Nbr_EW),
               sum(dat_ij$Nbr_EW[dat_ij$FIEC_LIEC=="LIEC_X_EPI"])/sum(dat_ij$Nbr_EW),
               sum(dat_ij$Nbr_EW[dat_ij$FIEC_LIEC=="LIEC_X_EPA"])/sum(dat_ij$Nbr_EW),
               sum(dat_ij$Nbr_EW[dat_ij$FIEC_LIEC=="LIEC_X_ANS"])/sum(dat_ij$Nbr_EW),
               sum(dat_ij$Nbr_EW[dat_ij$FIEC_LIEC=="LIEC_X_END"])/sum(dat_ij$Nbr_EW))
    vAbs=c(vAbs,
               sum(dat_ij$Nbr_EW[dat_ij$FIEC_LIEC=="EPI_EPI"]),
               sum(dat_ij$Nbr_EW[dat_ij$FIEC_LIEC=="EPA_EPA"]),
               sum(dat_ij$Nbr_EW[dat_ij$FIEC_LIEC=="ANS_ANS"]),
               sum(dat_ij$Nbr_EW[dat_ij$FIEC_LIEC=="END_END"]),
               sum(dat_ij$Nbr_EW[dat_ij$FIEC_LIEC=="LIEC_X_LIEC_X"]),
               sum(dat_ij$Nbr_EW[dat_ij$FIEC_LIEC=="EPI_EPA"]),
               sum(dat_ij$Nbr_EW[dat_ij$FIEC_LIEC=="EPI_ANS"]),
               sum(dat_ij$Nbr_EW[dat_ij$FIEC_LIEC=="EPI_END"]),
               sum(dat_ij$Nbr_EW[dat_ij$FIEC_LIEC=="EPI_LIEC_X"]),
               sum(dat_ij$Nbr_EW[dat_ij$FIEC_LIEC=="EPA_EPI"]),
               sum(dat_ij$Nbr_EW[dat_ij$FIEC_LIEC=="EPA_ANS"]),
               sum(dat_ij$Nbr_EW[dat_ij$FIEC_LIEC=="EPA_END"]),
               sum(dat_ij$Nbr_EW[dat_ij$FIEC_LIEC=="EPA_LIEC_X"]),
               sum(dat_ij$Nbr_EW[dat_ij$FIEC_LIEC=="ANS_EPI"]),
               sum(dat_ij$Nbr_EW[dat_ij$FIEC_LIEC=="ANS_EPA"]),
               sum(dat_ij$Nbr_EW[dat_ij$FIEC_LIEC=="ANS_END"]),
               sum(dat_ij$Nbr_EW[dat_ij$FIEC_LIEC=="ANS_LIEC_X"]),
               sum(dat_ij$Nbr_EW[dat_ij$FIEC_LIEC=="END_EPI"]),
               sum(dat_ij$Nbr_EW[dat_ij$FIEC_LIEC=="END_EPA"]),
               sum(dat_ij$Nbr_EW[dat_ij$FIEC_LIEC=="END_ANS"]),
               sum(dat_ij$Nbr_EW[dat_ij$FIEC_LIEC=="END_LIEC_X"]),
               sum(dat_ij$Nbr_EW[dat_ij$FIEC_LIEC=="LIEC_X_EPI"]),
               sum(dat_ij$Nbr_EW[dat_ij$FIEC_LIEC=="LIEC_X_EPA"]),
               sum(dat_ij$Nbr_EW[dat_ij$FIEC_LIEC=="LIEC_X_ANS"]),
               sum(dat_ij$Nbr_EW[dat_ij$FIEC_LIEC=="LIEC_X_END"]))
    vCombi=c(vCombi,
             "EPI_EPI",
             "EPA_EPA",
             "ANS_ANS",
             "END_END",
             "LIEC_X_LIEC_X",
             "EPI_EPA",
             "EPI_ANS",
             "EPI_END",
             "EPI_LIEC_X",
             "EPA_EPI",
             "EPA_ANS",
             "EPA_END",
             "EPA_LIEC_X",
             "ANS_EPI",
             "ANS_EPA",
             "ANS_END",
             "ANS_LIEC_X",
             "END_EPI",
             "END_EPA",
             "END_ANS",
             "END_LIEC_X",
             "LIEC_X_EPI",
             "LIEC_X_EPA",
             "LIEC_X_ANS",
             "LIEC_X_END")

  }
}



Tab_err=data.frame(
  Annee=vAnnee,
  Region=vRegion,
  Combi=vCombi,
  NbrVDT=vAbs,
  tauxErr=vtauxErr
)

#test

str(Tab_err)
str(Tab2[,c("Annee","Region","FIEC_LIEC")])

sum(Tab_err$tauxErr)==(nrow(Tab_err)/25)
nrow(anti_join(Tab2[,c("Annee","Region","FIEC_LIEC")],Tab_err[,c("Annee","Region","Combi")]))==0

#pour chaque Annee/Region, on garde que les r?ussites (sauf les LIEC_X)

Tab_reu=Tab_err[Tab_err$Combi%in%c("EPI_EPI","EPA_EPA","ANS_ANS","END_END"),]
nrow(unique(Tab_reu[,c("Annee","Region")]))*4==
  nrow(unique(Tab_reu[,c("Annee","Region","Combi")]))
nrow(unique(Tab_reu[,c("Annee","Region")]))==67

str(Tab_reu)
Tab_reufin=aggregate(Tab_reu$tauxErr,
                     by=Tab_reu[,c("Annee","Region")],
                     sum)
colnames(Tab_reufin)[colnames(Tab_reufin)=="x"]="tauxErr"

nrow(unique(Tab_reu[,c("Annee","Region")]))==nrow(Tab_reufin)
sum(Tab_reufin$tauxErr)==sum(Tab_reu$tauxErr)
sum(Tab_reufin$tauxErr)==sum(Tab_err$tauxErr[Tab_err$Combi%in%c("EPI_EPI","EPA_EPA","ANS_ANS","END_END")])


#j'ajoute le nombre de VDT, car on enlevera les nombre trop petits

Tab_NbrVDT=aggregate(Tab_err$NbrVDT,
                     by=Tab_err[,c("Annee","Region")],
                     sum)
colnames(Tab_NbrVDT)[colnames(Tab_NbrVDT)=="x"]="NbrVDT"




nrow(anti_join(unique(Tab_NbrVDT[,c("Annee","Region")]),
               unique(Tab_reufin[,c("Annee","Region")])))==0
nrow(anti_join(unique(Tab_reufin[,c("Annee","Region")]),
               unique(Tab_NbrVDT[,c("Annee","Region")])))==0

Tab_reu_Final=merge(Tab_reufin,Tab_NbrVDT,by=c("Annee","Region"))

nrow(anti_join(unique(Tab_reu_Final[,c("Annee","Region","tauxErr")]),
               unique(Tab_reufin[,c("Annee","Region","tauxErr")])))==0
nrow(anti_join(unique(Tab_reufin[,c("Annee","Region","tauxErr")]),
                unique(Tab_reu_Final[,c("Annee","Region","tauxErr")])))==0

nrow(anti_join(unique(Tab_reu_Final[,c("Annee","Region","NbrVDT")]),
               unique(Tab_NbrVDT[,c("Annee","Region","NbrVDT")])))==0
nrow(anti_join(unique(Tab_NbrVDT[,c("Annee","Region","NbrVDT")]),
              unique(Tab_reu_Final[,c("Annee","Region","NbrVDT")])))==0



#j'ajoute le nombre de LIEC_X id sur terrain, car on enlevera les nombre trop elev?s//NbrVDT (les gens s'en foutaient)

Tab_LIECX=Tab_err[Tab_err$Combi%in%c("LIEC_X_LIEC_X","LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END"),]
nrow(Tab_LIECX)==nrow(Tab_err)/5
Tab_LIECX=aggregate(Tab_LIECX$NbrVDT,
                  by=Tab_LIECX[,c("Annee","Region")],
                  sum)
colnames(Tab_LIECX)[colnames(Tab_LIECX)=="x"]="NbrLIECX"

nrow(Tab_LIECX)==nrow(Tab_reu_Final)

nrow(anti_join(unique(Tab_LIECX[,c("Annee","Region")]),
               unique(Tab_reu_Final[,c("Annee","Region")])))==0
nrow(anti_join(unique(Tab_reu_Final[,c("Annee","Region")]),
               unique(Tab_LIECX[,c("Annee","Region")])))==0



Tab_reu_Final2=merge(Tab_reu_Final,Tab_LIECX,by=c("Annee","Region"))

nrow(anti_join(unique(Tab_reu_Final2[,c("Annee","Region","tauxErr","NbrVDT")]),
               unique(Tab_reu_Final[,c("Annee","Region","tauxErr","NbrVDT")])))==0
nrow(anti_join(unique(Tab_reu_Final[,c("Annee","Region","tauxErr","NbrVDT")]),
               unique(Tab_reu_Final2[,c("Annee","Region","tauxErr","NbrVDT")])))==0

nrow(anti_join(unique(Tab_reu_Final2[,c("Annee","Region","NbrLIECX")]),
               unique(Tab_LIECX[,c("Annee","Region","NbrLIECX")])))==0
nrow(anti_join(unique(Tab_LIECX[,c("Annee","Region","NbrLIECX")]),
               unique(Tab_reu_Final2[,c("Annee","Region","NbrLIECX")])))==0


Tab_reu_Final2$titre="A"
Tab_reu_Final2$tauxLIECX=Tab_reu_Final2$NbrLIECX/Tab_reu_Final2$NbrVDT



Tab_reu_Final2Save=Tab_reu_Final2

#j'enl?ve les lignes avec un taux de LIECX trop ?lev?
nrow(Tab_reu_Final2)==67
Tab_reu_Final2=Tab_reu_Final2[Tab_reu_Final2$tauxLIECX<0.5,]
nrow(Tab_reu_Final2)==62

#j'enl?ve les lignes avec trop peu de VDT
nrow(Tab_reu_Final2)==62
Tab_reu_Final2=Tab_reu_Final2[Tab_reu_Final2$NbrVDT>62,]
nrow(Tab_reu_Final2)==57


ggplot(Tab_reu_Final2, aes(x=titre,y=tauxErr)) + 
  geom_jitter(color=rgb(76,159,104,max=255),width=0.2,size=5,height=0) +
  geom_boxplot(size=1.5,alpha = 0) + scale_x_discrete(drop = FALSE)+
  # scale_fill_manual(values=c("#BB4AFF", "#619CFF", "#9E9E9E"))+ 
  theme_bw(base_size=12) + # Mettre font blanc
  theme(title=element_text(face="bold", size=15),
        axis.text=element_text(face="bold", size=15),
        axis.title.x = element_text(face="bold", vjust=0),#size=25
        # axis.ticks.margin=unit(c(0.5),'cm'),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_line(size=2),
        axis.ticks.length=unit(0.5, "cm"),
        axis.text.x=element_blank(),
        panel.grid=element_blank(), # Augmenter taille ligne y
        panel.grid.major.x=element_blank(), # Suppression lignes axes x
        axis.title=element_text( vjust=5, size=15),
        line=element_line(size=1.5),
        legend.position="right",
        legend.title = element_text(face="bold", size=15),
        legend.text = element_text(size=15)) +
  stat_summary(fun.y=mean, geom="point", shape=18, size=10, col="black", fill = "black")+ #???shape=18
  scale_y_continuous(labels = scales::percent, limits=c(0,1)) + 
  ylab(NULL) + 
  xlab(NULL) 



ggplot(Tab_reu_Final2, aes(x=titre,y=tauxErr,fill=titre)) + 
  geom_boxplot(size=1.5) + scale_x_discrete(drop = FALSE)+
  # scale_fill_manual(values=c("#BB4AFF", "#619CFF", "#9E9E9E"))+ 
  theme_bw(base_size=12) + # Mettre font blanc
  theme(title=element_text(face="bold", size=15),
        axis.text=element_text(face="bold", size=15),
        axis.title.x = element_text(face="bold", vjust=0),#size=25
        # axis.ticks.margin=unit(c(0.5),'cm'),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_line(size=2),
        axis.ticks.length=unit(0.5, "cm"),
        axis.text.x=element_blank(),
        panel.grid=element_blank(), # Augmenter taille ligne y
        panel.grid.major.x=element_blank(), # Suppression lignes axes x
        axis.title=element_text( vjust=5, size=15),
        line=element_line(size=1.5),
        legend.position="none",
        legend.title = element_text(face="bold", size=15),
        legend.text = element_text(size=15)) +
  stat_summary(fun.y=mean, geom="point", shape=18, size=10, col="black", fill = "black")+ #???shape=18
  scale_y_continuous(labels = scales::percent, limits=c(0,1)) + 
  scale_fill_manual(values=c("#4C9F68")) + 
  ylab(NULL) + 
  xlab(NULL) 


ggplot(Tab_reu_Final2, aes(x=titre,y=tauxErr,fill=titre)) + 
  geom_boxplot(size=1.5) + scale_x_discrete(drop = FALSE)+
  # scale_fill_manual(values=c("#BB4AFF", "#619CFF", "#9E9E9E"))+ 
  theme_bw(base_size=12) + # Mettre font blanc
  theme(title=element_text(face="bold", size=15),
        axis.text=element_text(face="bold", size=15),
        axis.title.x = element_text(face="bold", vjust=0),#size=25
        # axis.ticks.margin=unit(c(0.5),'cm'),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_line(size=2),
        axis.ticks.length=unit(0.5, "cm"),
        axis.text.x=element_blank(),
        panel.grid=element_blank(), # Augmenter taille ligne y
        panel.grid.major.x=element_blank(), # Suppression lignes axes x
        axis.title=element_text( vjust=5, size=15),
        line=element_line(size=1.5),
        legend.position="none",
        legend.title = element_text(face="bold", size=15),
        legend.text = element_text(size=15)) +
  stat_summary(fun.y=mean, geom="point", shape=18, size=10, col="black", fill = "black")+ #???shape=18
  scale_y_continuous(labels = scales::percent, limits=c(0,1)) + 
  scale_fill_manual(values=c("grey")) + 
  ylab(NULL) + 
  xlab(NULL) 



##### MAUVAIS GROUPE AL 2012 #####



tabAL=Tab2[Tab2$Region%in%"AL",]
nrow(tabAL)==1309
dataID <- unique(sbt_m[,c("cle",colID,"Region")])
nrow(dataID)==1109
dataID=dataID[dataID$Region%in%"AL",]
nrow(dataID)==99




tabAL2012=tabAL[tabAL$Annee%in%2012,]
nrow(tabAL2012)==244

#taux d'id correcte global
sum(tabAL2012$Nbr_EW[tabAL2012$FIEC_LIEC%in%c("EPI_EPI","EPA_EPA","ANS_ANS","END_END")])/
  sum(tabAL2012$Nbr_EW)#0.4008499  => ok
#Tab_reu_Final2$tauxErr[Tab_reu_Final2$Region=="AL" & Tab_reu_Final2$Annee==2012]


#les EPI identifi?s : que sont ils r?ellement ?
{
  #nombre id correcte EPI
  n2.EPI_EPI=sum(tabAL2012$Nbr_EW[tabAL2012$FIEC_LIEC%in%c("EPI_EPI")])
  
  #nombre id incorrecte EPA 
  n2.EPI_EPA=sum(tabAL2012$Nbr_EW[tabAL2012$FIEC_LIEC%in%c("EPI_EPA")])
  
  #nombre id incorrecte ANS 
  n2.EPI_ANS=sum(tabAL2012$Nbr_EW[tabAL2012$FIEC_LIEC%in%c("EPI_ANS")])
  
  #nombre id incorrecte END 
  n2.EPI_END=sum(tabAL2012$Nbr_EW[tabAL2012$FIEC_LIEC%in%c("EPI_END")])
  
  #nombre id incorrecte LIEC_X 
  n2.EPI_LIECX=sum(tabAL2012$Nbr_EW[tabAL2012$FIEC_LIEC%in%c("EPI_LIEC_X")])
  
  
  
  
  #taux id correcte EPI
  t2.EPI_EPI=sum(tabAL2012$Nbr_EW[tabAL2012$FIEC_LIEC%in%c("EPI_EPI")])/
    sum(tabAL2012$Nbr_EW[tabAL2012$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])

  #taux id incorrecte EPA 
  t2.EPI_EPA=sum(tabAL2012$Nbr_EW[tabAL2012$FIEC_LIEC%in%c("EPI_EPA")])/
    sum(tabAL2012$Nbr_EW[tabAL2012$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])

  #taux id incorrecte ANS 
  t2.EPI_ANS=sum(tabAL2012$Nbr_EW[tabAL2012$FIEC_LIEC%in%c("EPI_ANS")])/
    sum(tabAL2012$Nbr_EW[tabAL2012$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])

  #taux id incorrecte END 
  t2.EPI_END=sum(tabAL2012$Nbr_EW[tabAL2012$FIEC_LIEC%in%c("EPI_END")])/
    sum(tabAL2012$Nbr_EW[tabAL2012$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])

  #taux id incorrecte LIEC_X 
  t2.EPI_LIECX=sum(tabAL2012$Nbr_EW[tabAL2012$FIEC_LIEC%in%c("EPI_LIEC_X")])/
    sum(tabAL2012$Nbr_EW[tabAL2012$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
  
  t2.EPI_EPI+
    t2.EPI_EPA+
    t2.EPI_ANS+
    t2.EPI_END+
    t2.EPI_LIECX==1
  
  
}

#les EPA identifi?s : que sont ils r?ellement ?
{
  
  #nombre id correcte EPI
  n2.EPA_EPI=sum(tabAL2012$Nbr_EW[tabAL2012$FIEC_LIEC%in%c("EPA_EPI")])
  
  #nombre id incorrecte EPA 
  n2.EPA_EPA=sum(tabAL2012$Nbr_EW[tabAL2012$FIEC_LIEC%in%c("EPA_EPA")])
  
  #nombre id incorrecte ANS 
  n2.EPA_ANS=sum(tabAL2012$Nbr_EW[tabAL2012$FIEC_LIEC%in%c("EPA_ANS")])
  
  #nombre id incorrecte END 
  n2.EPA_END=sum(tabAL2012$Nbr_EW[tabAL2012$FIEC_LIEC%in%c("EPA_END")])
  
  #nombre id incorrecte LIEC_X 
  n2.EPA_LIECX=sum(tabAL2012$Nbr_EW[tabAL2012$FIEC_LIEC%in%c("EPA_LIEC_X")])
  
  
  
  
  
  
  #taux id correcte EPI
  t2.EPA_EPI=sum(tabAL2012$Nbr_EW[tabAL2012$FIEC_LIEC%in%c("EPA_EPI")])/
    sum(tabAL2012$Nbr_EW[tabAL2012$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  
  
  
  #taux id incorrecte EPA 
  t2.EPA_EPA=sum(tabAL2012$Nbr_EW[tabAL2012$FIEC_LIEC%in%c("EPA_EPA")])/
    sum(tabAL2012$Nbr_EW[tabAL2012$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  
  
  
  #taux id incorrecte ANS 
  t2.EPA_ANS=sum(tabAL2012$Nbr_EW[tabAL2012$FIEC_LIEC%in%c("EPA_ANS")])/
    sum(tabAL2012$Nbr_EW[tabAL2012$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  
  
  
  
  #taux id incorrecte END 
  t2.EPA_END=sum(tabAL2012$Nbr_EW[tabAL2012$FIEC_LIEC%in%c("EPA_END")])/
    sum(tabAL2012$Nbr_EW[tabAL2012$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  
  
  
  #taux id incorrecte LIEC_X 
  t2.EPA_LIECX=sum(tabAL2012$Nbr_EW[tabAL2012$FIEC_LIEC%in%c("EPA_LIEC_X")])/
    sum(tabAL2012$Nbr_EW[tabAL2012$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
  
  
  
  t2.EPA_EPI+
    t2.EPA_EPA+
    t2.EPA_ANS+
    t2.EPA_END+
    t2.EPA_LIECX==1
  
  
}

#les ANS identifi?s : que sont ils r?ellement ?
{
  
  #nombre id correcte EPI
  n2.ANS_EPI=sum(tabAL2012$Nbr_EW[tabAL2012$FIEC_LIEC%in%c("ANS_EPI")])
  
  #nombre id incorrecte EPA 
  n2.ANS_EPA=sum(tabAL2012$Nbr_EW[tabAL2012$FIEC_LIEC%in%c("ANS_EPA")])
  
  #nombre id incorrecte ANS 
  n2.ANS_ANS=sum(tabAL2012$Nbr_EW[tabAL2012$FIEC_LIEC%in%c("ANS_ANS")])
  
  #nombre id incorrecte END 
  n2.ANS_END=sum(tabAL2012$Nbr_EW[tabAL2012$FIEC_LIEC%in%c("ANS_END")])
  
  #nombre id incorrecte LIEC_X 
  n2.ANS_LIECX=sum(tabAL2012$Nbr_EW[tabAL2012$FIEC_LIEC%in%c("ANS_LIEC_X")])
  
  
  
  
  
  
  #taux id correcte EPI
  t2.ANS_EPI=sum(tabAL2012$Nbr_EW[tabAL2012$FIEC_LIEC%in%c("ANS_EPI")])/
    sum(tabAL2012$Nbr_EW[tabAL2012$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  
  
  
  #taux id incorrecte EPA 
  t2.ANS_EPA=sum(tabAL2012$Nbr_EW[tabAL2012$FIEC_LIEC%in%c("ANS_EPA")])/
    sum(tabAL2012$Nbr_EW[tabAL2012$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  
  
  
  #taux id incorrecte ANS 
  t2.ANS_ANS=sum(tabAL2012$Nbr_EW[tabAL2012$FIEC_LIEC%in%c("ANS_ANS")])/
    sum(tabAL2012$Nbr_EW[tabAL2012$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  
  
  
  
  #taux id incorrecte END 
  t2.ANS_END=sum(tabAL2012$Nbr_EW[tabAL2012$FIEC_LIEC%in%c("ANS_END")])/
    sum(tabAL2012$Nbr_EW[tabAL2012$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  
  
  
  #taux id incorrecte LIEC_X 
  t2.ANS_LIECX=sum(tabAL2012$Nbr_EW[tabAL2012$FIEC_LIEC%in%c("ANS_LIEC_X")])/
    sum(tabAL2012$Nbr_EW[tabAL2012$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
  
  
  
  t2.ANS_EPI+
    t2.ANS_EPA+
    t2.ANS_ANS+
    t2.ANS_END+
    t2.ANS_LIECX==1
  
  
}

#les END identifi?s : que sont ils r?ellement ?
{
  
  #nombre id correcte EPI
  n2.END_EPI=sum(tabAL2012$Nbr_EW[tabAL2012$FIEC_LIEC%in%c("END_EPI")])
  
  #nombre id incorrecte EPA 
  n2.END_EPA=sum(tabAL2012$Nbr_EW[tabAL2012$FIEC_LIEC%in%c("END_EPA")])
  
  #nombre id incorrecte ANS 
  n2.END_ANS=sum(tabAL2012$Nbr_EW[tabAL2012$FIEC_LIEC%in%c("END_ANS")])
  
  #nombre id incorrecte END 
  n2.END_END=sum(tabAL2012$Nbr_EW[tabAL2012$FIEC_LIEC%in%c("END_END")])
  
  #nombre id incorrecte LIEC_X 
  n2.END_LIECX=sum(tabAL2012$Nbr_EW[tabAL2012$FIEC_LIEC%in%c("END_LIEC_X")])
  
  
  
  
  
  
  #taux id correcte EPI
  t2.END_EPI=sum(tabAL2012$Nbr_EW[tabAL2012$FIEC_LIEC%in%c("END_EPI")])/
    sum(tabAL2012$Nbr_EW[tabAL2012$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  
  
  
  #taux id incorrecte EPA 
  t2.END_EPA=sum(tabAL2012$Nbr_EW[tabAL2012$FIEC_LIEC%in%c("END_EPA")])/
    sum(tabAL2012$Nbr_EW[tabAL2012$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  
  
  
  #taux id incorrecte ANS 
  t2.END_ANS=sum(tabAL2012$Nbr_EW[tabAL2012$FIEC_LIEC%in%c("END_ANS")])/
    sum(tabAL2012$Nbr_EW[tabAL2012$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  
  
  
  
  #taux id incorrecte END 
  t2.END_END=sum(tabAL2012$Nbr_EW[tabAL2012$FIEC_LIEC%in%c("END_END")])/
    sum(tabAL2012$Nbr_EW[tabAL2012$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  
  
  
  #taux id incorrecte LIEC_X 
  t2.END_LIECX=sum(tabAL2012$Nbr_EW[tabAL2012$FIEC_LIEC%in%c("END_LIEC_X")])/
    sum(tabAL2012$Nbr_EW[tabAL2012$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
  
  
  
  t2.END_EPI+
    t2.END_EPA+
    t2.END_ANS+
    t2.END_END+
    t2.END_LIECX==1
  
  
}

#les LIEC_X identifi?s : que sont ils r?ellement ?
{
  
  #nombre id correcte EPI
  n2.LIECX_EPI=sum(tabAL2012$Nbr_EW[tabAL2012$FIEC_LIEC%in%c("LIEC_X_EPI")])
  
  #nombre id incorrecte EPA 
  n2.LIECX_EPA=sum(tabAL2012$Nbr_EW[tabAL2012$FIEC_LIEC%in%c("LIEC_X_EPA")])
  
  #nombre id incorrecte ANS 
  n2.LIECX_ANS=sum(tabAL2012$Nbr_EW[tabAL2012$FIEC_LIEC%in%c("LIEC_X_ANS")])
  
  #nombre id incorrecte END 
  n2.LIECX_END=sum(tabAL2012$Nbr_EW[tabAL2012$FIEC_LIEC%in%c("LIEC_X_END")])
  
  #nombre id incorrecte LIEC_X 
  n2.LIECX_LIECX=sum(tabAL2012$Nbr_EW[tabAL2012$FIEC_LIEC%in%c("LIEC_X_LIEC_X")])
  
  
  
  
  
  #taux id correcte EPI
  t2.LIECX_EPI=sum(tabAL2012$Nbr_EW[tabAL2012$FIEC_LIEC%in%c("LIEC_X_EPI")])/
    sum(tabAL2012$Nbr_EW[tabAL2012$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  
  
  
  #taux id incorrecte EPA 
  t2.LIECX_EPA=sum(tabAL2012$Nbr_EW[tabAL2012$FIEC_LIEC%in%c("LIEC_X_EPA")])/
    sum(tabAL2012$Nbr_EW[tabAL2012$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  
  
  
  #taux id incorrecte ANS 
  t2.LIECX_ANS=sum(tabAL2012$Nbr_EW[tabAL2012$FIEC_LIEC%in%c("LIEC_X_ANS")])/
    sum(tabAL2012$Nbr_EW[tabAL2012$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  
  
  
  
  #taux id incorrecte END 
  t2.LIECX_END=sum(tabAL2012$Nbr_EW[tabAL2012$FIEC_LIEC%in%c("LIEC_X_END")])/
    sum(tabAL2012$Nbr_EW[tabAL2012$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  
  
  
  #taux id incorrecte LIEC_X 
  t2.LIECX_LIECX=sum(tabAL2012$Nbr_EW[tabAL2012$FIEC_LIEC%in%c("LIEC_X_LIEC_X")])/
    sum(tabAL2012$Nbr_EW[tabAL2012$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
  
  
  
  t2.LIECX_EPI+
    t2.LIECX_EPA+
    t2.LIECX_ANS+
    t2.LIECX_END+
    t2.LIECX_LIECX==1
  
  
}



tabGraphAL12_2=data.frame(
  Nbr=c(
    n2.EPI_EPI,
    n2.EPI_EPA,
    n2.EPI_ANS,
    n2.EPI_END,
    n2.EPI_LIECX,
    n2.EPA_EPI,
    n2.EPA_EPA,
    n2.EPA_ANS,
    n2.EPA_END,
    n2.EPA_LIECX,
    n2.ANS_EPI,
    n2.ANS_EPA,
    n2.ANS_ANS,
    n2.ANS_END,
    n2.ANS_LIECX,
    n2.END_EPI,
    n2.END_EPA,
    n2.END_ANS,
    n2.END_END,
    n2.END_LIECX,
    n2.LIECX_EPI,
    n2.LIECX_EPA,
    n2.LIECX_ANS,
    n2.LIECX_END,
    n2.LIECX_LIECX),
  
  Taux=c(
    t2.EPI_EPI,
    t2.EPI_EPA,
    t2.EPI_ANS,
    t2.EPI_END,
    t2.EPI_LIECX,
    t2.EPA_EPI,
    t2.EPA_EPA,
    t2.EPA_ANS,
    t2.EPA_END,
    t2.EPA_LIECX,
    t2.ANS_EPI,
    t2.ANS_EPA,
    t2.ANS_ANS,
    t2.ANS_END,
    t2.ANS_LIECX,
    t2.END_EPI,
    t2.END_EPA,
    t2.END_ANS,
    t2.END_END,
    t2.END_LIECX,
    t2.LIECX_EPI,
    t2.LIECX_EPA,
    t2.LIECX_ANS,
    t2.LIECX_END,
    t2.LIECX_LIECX),
  
  Classe=c(
    "t2.EPI_EPI",
    "t2.EPI_EPA",
    "t2.EPI_ANS",
    "t2.EPI_END",
    "t2.EPI_LIECX",
    "t2.EPA_EPI",
    "t2.EPA_EPA",
    "t2.EPA_ANS",
    "t2.EPA_END",
    "t2.EPA_LIECX",
    "t2.ANS_EPI",
    "t2.ANS_EPA",
    "t2.ANS_ANS",
    "t2.ANS_END",
    "t2.ANS_LIECX",
    "t2.END_EPI",
    "t2.END_EPA",
    "t2.END_ANS",
    "t2.END_END",
    "t2.END_LIECX",
    "t2.LIECX_EPI",
    "t2.LIECX_EPA",
    "t2.LIECX_ANS",
    "t2.LIECX_END",
    "t2.LIECX_LIECX")
  
)



tabGraphAL12_2=tabGraphAL12_2[c(5,4,3,2,1,10,9,8,7,6,15,14,13,12,11,20,19,18,17,16,25,24,23,22,21),]
tabGraphAL12_2$LIEC=factor(substr(tabGraphAL12_2$Classe,8,10),levels=c("EPI","EPA","ANS","END","LIECX"))
tabGraphAL12_2$LIECid=factor(substr(tabGraphAL12_2$Classe,4,6),levels=c("EPI","EPA","ANS","END","LIECX"))
anti_join(tabGraphAL12_2,ddply(tabGraphAL12_2, .(Classe), mutate, Pos_Y=cumsum(Nbr)-Nbr/2 ))
anti_join(ddply(tabGraphAL12_2, .(Classe), mutate, Pos_Y=cumsum(Nbr)-Nbr/2 ),tabGraphAL12_2)
tabGraphAL12_2 <- ddply(tabGraphAL12_2, .(LIECid), mutate, Pos_Y=cumsum(Nbr)-Nbr/2 )

tabGraphAL12_2=tabGraphAL12_2[tabGraphAL12_2$Nbr>0,]
tabGraphAL12_2=tabGraphAL12_2[tabGraphAL12_2$LIEC!="LIECX",]

vpos_X=c()
for(iLIEC in levels(tabGraphAL12_2$LIECid)){
  intmd=seq(0.4,-0.4,length.out = nrow(tabGraphAL12_2[tabGraphAL12_2$LIECid%in%iLIEC,]))
  vpos_X=c(vpos_X,intmd+which(levels(tabGraphAL12_2$LIECid)%in%iLIEC))
}
tabGraphAL12_2$pos_X <- vpos_X

tabGraphAL12_2$couleur=rep(NA,nrow(tabGraphAL12_2))
tabGraphAL12_2$couleur[tabGraphAL12_2$LIEC%in%"EPI"]="#000099"
tabGraphAL12_2$couleur[tabGraphAL12_2$LIEC%in%"EPA"]="#FFFF00"
tabGraphAL12_2$couleur[tabGraphAL12_2$LIEC%in%"ANS"]="#009900"
tabGraphAL12_2$couleur[tabGraphAL12_2$LIEC%in%"END"]="#FF6699"
tabGraphAL12_2$couleur[tabGraphAL12_2$LIEC%in%"LIECX"]="#7B7B7B"



tabGraphAL12_2$colecr=rep(NA,nrow(tabGraphAL12_2))
tabGraphAL12_2$colecr[tabGraphAL12_2$LIEC%in%"EPI"]="#FFFFFF"
tabGraphAL12_2$colecr[tabGraphAL12_2$LIEC%in%"EPA"]="#000000"
tabGraphAL12_2$colecr[tabGraphAL12_2$LIEC%in%"ANS"]="#FFFFFF"
tabGraphAL12_2$colecr[tabGraphAL12_2$LIEC%in%"END"]="#000000"
tabGraphAL12_2$colecr[tabGraphAL12_2$LIEC%in%"LIECX"]="#FFFFFF"






ggplot(tabGraphAL12_2, aes(fill=LIEC, y=Taux, x=LIECid)) + 
  geom_bar(position='stack', stat='identity')+
  theme_classic()+
  scale_fill_manual(values = c("#000099","#FFFF00","#009900","#FF6699","#7B7B7B"))

ggplot(tabGraphAL12_2, aes(fill=LIEC, y=Taux, x=LIECid)) + 
  geom_bar(position='fill', stat='identity')+
  theme_classic()+
  scale_fill_manual(values = c("#000099","#FFFF00","#009900","#FF6699","#7B7B7B"))

ggplot(tabGraphAL12_2, aes(fill=LIEC, y=Nbr, x=LIECid)) + 
  geom_bar(position='fill', stat='identity')+
  theme_classic()+
  scale_fill_manual(values = c("#000099","#FFFF00","#009900","#FF6699","#7B7B7B"))

ggplot(tabGraphAL12_2, aes(fill=LIEC, y=Nbr, x=LIECid)) + 
  geom_bar(position='stack', stat='identity')+
  theme_classic()+
  scale_fill_manual(values = c("#000099","#FFFF00","#009900","#FF6699","#7B7B7B"))







#zoom sur les erreurs
tabGraphAL12_2$Erreur=NA
tabGraphAL12_2$Erreur[tabGraphAL12_2$LIEC==tabGraphAL12_2$LIECid]="Correcte"
tabGraphAL12_2$Erreur[tabGraphAL12_2$LIEC!=tabGraphAL12_2$LIECid]="Incorrecte"



ggplot(tabGraphAL12_2[tabGraphAL12_2$Erreur=="Incorrecte" & tabGraphAL12_2$LIECid!="LIECX",], 
       aes(fill=LIEC, y=Taux, x=LIECid)) + 
  geom_bar(position='stack', stat='identity')+
  theme_classic()+
  theme(title=element_text(face="bold", size=15),
        axis.text=element_text(face="bold", size=15),
        axis.title.x = element_text(face="bold", vjust=0),#size=25
        # axis.ticks.margin=unit(c(0.5),'cm'),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_line(size=2),
        axis.ticks.length=unit(0.5, "cm"),
        axis.text.x=element_blank(),
        panel.grid=element_blank(), # Augmenter taille ligne y
        panel.grid.major.x=element_blank(), # Suppression lignes axes x
        axis.title=element_text( vjust=5, size=15),
        line=element_line(size=1.5),
        legend.position="none",
        legend.title = element_text(face="bold", size=15),
        legend.text = element_text(size=15)) +
  scale_y_continuous(labels = scales::percent, limits=c(0,1)) + 
  scale_fill_manual(values = c("#000099","#FFFF00","#009900","#FF6699","#7B7B7B"))+
  xlab(NULL)+ylab(NULL)



ggplot(tabGraphAL12_2[tabGraphAL12_2$Erreur=="Incorrecte" & tabGraphAL12_2$LIECid!="LIECX",], 
       aes(fill=LIEC, y=Nbr, x=LIECid)) + 
  geom_bar(position='stack', stat='identity')+
  theme_classic()+
  theme(title=element_text(face="bold", size=15),
        axis.text=element_text(face="bold", size=15),
        axis.title.x = element_text(face="bold", vjust=0),#size=25
        # axis.ticks.margin=unit(c(0.5),'cm'),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_line(size=2),
        axis.ticks.length=unit(0.5, "cm"),
        axis.text.x=element_blank(),
        panel.grid=element_blank(), # Augmenter taille ligne y
        panel.grid.major.x=element_blank(), # Suppression lignes axes x
        axis.title=element_text( vjust=5, size=15),
        line=element_line(size=1.5),
        legend.position="none",
        legend.title = element_text(face="bold", size=15),
        legend.text = element_text(size=15)) +
  scale_fill_manual(values = c("#000099","#FFFF00","#009900","#FF6699","#7B7B7B"))+
  xlab(NULL)+ylab(NULL)


ggplot(tabGraphAL12_2[tabGraphAL12_2$LIECid!="LIECX",], 
       aes(fill=LIEC, y=Nbr, x=LIECid)) + 
  geom_bar(position='stack', stat='identity')+
  theme_classic()+
  theme(title=element_text(face="bold", size=15),
        axis.text=element_text(face="bold", size=15),
        axis.title.x = element_text(face="bold", vjust=0),#size=25
        # axis.ticks.margin=unit(c(0.5),'cm'),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_line(size=2),
        axis.ticks.length=unit(0.5, "cm"),
        axis.text.x=element_blank(),
        panel.grid=element_blank(), # Augmenter taille ligne y
        panel.grid.major.x=element_blank(), # Suppression lignes axes x
        axis.title=element_text( vjust=5, size=15),
        line=element_line(size=1.5),
        legend.position="none",
        legend.title = element_text(face="bold", size=15),
        legend.text = element_text(size=15)) +
  scale_fill_manual(values = c("#000099","#FFFF00","#009900","#FF6699","#7B7B7B"))+
  xlab(NULL)+ylab(NULL)+ylim(0,300)+


  geom_point(aes(y=Pos_Y, x=pos_X ),color="grey", size=14, show.legend=F)+  # Point pour entourer le nb de VDT
  geom_point(aes(y=Pos_Y, x=pos_X ),color=tabGraphAL12_2$couleur[tabGraphAL12_2$LIECid!="LIECX"], size=12, show.legend=F)+  # Point pour entourer le nb de VDT
  geom_text(aes(y=Pos_Y, x=pos_X ), label=tabGraphAL12_2$Nbr[tabGraphAL12_2$LIECid!="LIECX"], # Affichage valeur Ab >0
            size=5, fontface=2, color=tabGraphAL12_2$colecr[tabGraphAL12_2$LIECid!="LIECX"], show.legend=F)    # Augmenter taille et font (gras)


  
##### BON GROUPE NP 2015 #####
  
  tabPC=Tab2[Tab2$Region%in%"NP",]
  nrow(tabPC)==1208
  dataID <- unique(sbt_m[,c("cle",colID,"Region")])
  nrow(dataID)==1109
  dataID=dataID[dataID$Region%in%"NP",]
  nrow(dataID)==78
  
  
  
  
  tabNP2015=tabPC[tabPC$Annee%in%2015,]
  nrow(tabNP2015)==150
  
  #taux d'id correcte global
  sum(tabNP2015$Nbr_EW[tabNP2015$FIEC_LIEC%in%c("EPI_EPI","EPA_EPA","ANS_ANS","END_END")])/
    sum(tabNP2015$Nbr_EW)#0.9653179  => ok
  #Tab_reu_Final2$tauxErr[Tab_reu_Final2$Region=="NP" & Tab_reu_Final2$Annee==2015]
  
  
  #les EPI identifi?s : que sont ils r?ellement ?
  {
    #nombre id correcte EPI
    n2.EPI_EPI=sum(tabNP2015$Nbr_EW[tabNP2015$FIEC_LIEC%in%c("EPI_EPI")])
    
    #nombre id incorrecte EPA 
    n2.EPI_EPA=sum(tabNP2015$Nbr_EW[tabNP2015$FIEC_LIEC%in%c("EPI_EPA")])
    
    #nombre id incorrecte ANS 
    n2.EPI_ANS=sum(tabNP2015$Nbr_EW[tabNP2015$FIEC_LIEC%in%c("EPI_ANS")])
    
    #nombre id incorrecte END 
    n2.EPI_END=sum(tabNP2015$Nbr_EW[tabNP2015$FIEC_LIEC%in%c("EPI_END")])
    
    #nombre id incorrecte LIEC_X 
    n2.EPI_LIECX=sum(tabNP2015$Nbr_EW[tabNP2015$FIEC_LIEC%in%c("EPI_LIEC_X")])
    
    
    
    
    #taux id correcte EPI
    t2.EPI_EPI=sum(tabNP2015$Nbr_EW[tabNP2015$FIEC_LIEC%in%c("EPI_EPI")])/
      sum(tabNP2015$Nbr_EW[tabNP2015$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
    
    #taux id incorrecte EPA 
    t2.EPI_EPA=sum(tabNP2015$Nbr_EW[tabNP2015$FIEC_LIEC%in%c("EPI_EPA")])/
      sum(tabNP2015$Nbr_EW[tabNP2015$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
    
    #taux id incorrecte ANS 
    t2.EPI_ANS=sum(tabNP2015$Nbr_EW[tabNP2015$FIEC_LIEC%in%c("EPI_ANS")])/
      sum(tabNP2015$Nbr_EW[tabNP2015$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
    
    #taux id incorrecte END 
    t2.EPI_END=sum(tabNP2015$Nbr_EW[tabNP2015$FIEC_LIEC%in%c("EPI_END")])/
      sum(tabNP2015$Nbr_EW[tabNP2015$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
    
    #taux id incorrecte LIEC_X 
    t2.EPI_LIECX=sum(tabNP2015$Nbr_EW[tabNP2015$FIEC_LIEC%in%c("EPI_LIEC_X")])/
      sum(tabNP2015$Nbr_EW[tabNP2015$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
    
    t2.EPI_EPI+
      t2.EPI_EPA+
      t2.EPI_ANS+
      t2.EPI_END+
      t2.EPI_LIECX==1
    
    
    }
  
  #les EPA identifi?s : que sont ils r?ellement ?
  {
    
    #nombre id correcte EPI
    n2.EPA_EPI=sum(tabNP2015$Nbr_EW[tabNP2015$FIEC_LIEC%in%c("EPA_EPI")])
    
    #nombre id incorrecte EPA 
    n2.EPA_EPA=sum(tabNP2015$Nbr_EW[tabNP2015$FIEC_LIEC%in%c("EPA_EPA")])
    
    #nombre id incorrecte ANS 
    n2.EPA_ANS=sum(tabNP2015$Nbr_EW[tabNP2015$FIEC_LIEC%in%c("EPA_ANS")])
    
    #nombre id incorrecte END 
    n2.EPA_END=sum(tabNP2015$Nbr_EW[tabNP2015$FIEC_LIEC%in%c("EPA_END")])
    
    #nombre id incorrecte LIEC_X 
    n2.EPA_LIECX=sum(tabNP2015$Nbr_EW[tabNP2015$FIEC_LIEC%in%c("EPA_LIEC_X")])
    
    
    
    
    
    
    #taux id correcte EPI
    t2.EPA_EPI=sum(tabNP2015$Nbr_EW[tabNP2015$FIEC_LIEC%in%c("EPA_EPI")])/
      sum(tabNP2015$Nbr_EW[tabNP2015$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
    
    
    
    
    #taux id incorrecte EPA 
    t2.EPA_EPA=sum(tabNP2015$Nbr_EW[tabNP2015$FIEC_LIEC%in%c("EPA_EPA")])/
      sum(tabNP2015$Nbr_EW[tabNP2015$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
    
    
    
    
    #taux id incorrecte ANS 
    t2.EPA_ANS=sum(tabNP2015$Nbr_EW[tabNP2015$FIEC_LIEC%in%c("EPA_ANS")])/
      sum(tabNP2015$Nbr_EW[tabNP2015$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
    
    
    
    
    
    #taux id incorrecte END 
    t2.EPA_END=sum(tabNP2015$Nbr_EW[tabNP2015$FIEC_LIEC%in%c("EPA_END")])/
      sum(tabNP2015$Nbr_EW[tabNP2015$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
    
    
    
    
    #taux id incorrecte LIEC_X 
    t2.EPA_LIECX=sum(tabNP2015$Nbr_EW[tabNP2015$FIEC_LIEC%in%c("EPA_LIEC_X")])/
      sum(tabNP2015$Nbr_EW[tabNP2015$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
    
    
    
    t2.EPA_EPI+
      t2.EPA_EPA+
      t2.EPA_ANS+
      t2.EPA_END+
      t2.EPA_LIECX==1
    
    
  }
  
  #les ANS identifi?s : que sont ils r?ellement ?
  {
    
    #nombre id correcte EPI
    n2.ANS_EPI=sum(tabNP2015$Nbr_EW[tabNP2015$FIEC_LIEC%in%c("ANS_EPI")])
    
    #nombre id incorrecte EPA 
    n2.ANS_EPA=sum(tabNP2015$Nbr_EW[tabNP2015$FIEC_LIEC%in%c("ANS_EPA")])
    
    #nombre id incorrecte ANS 
    n2.ANS_ANS=sum(tabNP2015$Nbr_EW[tabNP2015$FIEC_LIEC%in%c("ANS_ANS")])
    
    #nombre id incorrecte END 
    n2.ANS_END=sum(tabNP2015$Nbr_EW[tabNP2015$FIEC_LIEC%in%c("ANS_END")])
    
    #nombre id incorrecte LIEC_X 
    n2.ANS_LIECX=sum(tabNP2015$Nbr_EW[tabNP2015$FIEC_LIEC%in%c("ANS_LIEC_X")])
    
    
    
    
    
    
    #taux id correcte EPI
    t2.ANS_EPI=sum(tabNP2015$Nbr_EW[tabNP2015$FIEC_LIEC%in%c("ANS_EPI")])/
      sum(tabNP2015$Nbr_EW[tabNP2015$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
    
    
    
    
    #taux id incorrecte EPA 
    t2.ANS_EPA=sum(tabNP2015$Nbr_EW[tabNP2015$FIEC_LIEC%in%c("ANS_EPA")])/
      sum(tabNP2015$Nbr_EW[tabNP2015$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
    
    
    
    
    #taux id incorrecte ANS 
    t2.ANS_ANS=sum(tabNP2015$Nbr_EW[tabNP2015$FIEC_LIEC%in%c("ANS_ANS")])/
      sum(tabNP2015$Nbr_EW[tabNP2015$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
    
    
    
    
    
    #taux id incorrecte END 
    t2.ANS_END=sum(tabNP2015$Nbr_EW[tabNP2015$FIEC_LIEC%in%c("ANS_END")])/
      sum(tabNP2015$Nbr_EW[tabNP2015$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
    
    
    
    
    #taux id incorrecte LIEC_X 
    t2.ANS_LIECX=sum(tabNP2015$Nbr_EW[tabNP2015$FIEC_LIEC%in%c("ANS_LIEC_X")])/
      sum(tabNP2015$Nbr_EW[tabNP2015$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
    
    
    
    t2.ANS_EPI+
      t2.ANS_EPA+
      t2.ANS_ANS+
      t2.ANS_END+
      t2.ANS_LIECX==1
    
    
  }
  
  #les END identifi?s : que sont ils r?ellement ?
  {
    
    #nombre id correcte EPI
    n2.END_EPI=sum(tabNP2015$Nbr_EW[tabNP2015$FIEC_LIEC%in%c("END_EPI")])
    
    #nombre id incorrecte EPA 
    n2.END_EPA=sum(tabNP2015$Nbr_EW[tabNP2015$FIEC_LIEC%in%c("END_EPA")])
    
    #nombre id incorrecte ANS 
    n2.END_ANS=sum(tabNP2015$Nbr_EW[tabNP2015$FIEC_LIEC%in%c("END_ANS")])
    
    #nombre id incorrecte END 
    n2.END_END=sum(tabNP2015$Nbr_EW[tabNP2015$FIEC_LIEC%in%c("END_END")])
    
    #nombre id incorrecte LIEC_X 
    n2.END_LIECX=sum(tabNP2015$Nbr_EW[tabNP2015$FIEC_LIEC%in%c("END_LIEC_X")])
    
    
    
    
    
    
    #taux id correcte EPI
    t2.END_EPI=sum(tabNP2015$Nbr_EW[tabNP2015$FIEC_LIEC%in%c("END_EPI")])/
      sum(tabNP2015$Nbr_EW[tabNP2015$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
    
    
    
    
    #taux id incorrecte EPA 
    t2.END_EPA=sum(tabNP2015$Nbr_EW[tabNP2015$FIEC_LIEC%in%c("END_EPA")])/
      sum(tabNP2015$Nbr_EW[tabNP2015$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
    
    
    
    
    #taux id incorrecte ANS 
    t2.END_ANS=sum(tabNP2015$Nbr_EW[tabNP2015$FIEC_LIEC%in%c("END_ANS")])/
      sum(tabNP2015$Nbr_EW[tabNP2015$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
    
    
    
    
    
    #taux id incorrecte END 
    t2.END_END=sum(tabNP2015$Nbr_EW[tabNP2015$FIEC_LIEC%in%c("END_END")])/
      sum(tabNP2015$Nbr_EW[tabNP2015$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
    
    
    
    
    #taux id incorrecte LIEC_X 
    t2.END_LIECX=sum(tabNP2015$Nbr_EW[tabNP2015$FIEC_LIEC%in%c("END_LIEC_X")])/
      sum(tabNP2015$Nbr_EW[tabNP2015$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
    
    
    
    t2.END_EPI+
      t2.END_EPA+
      t2.END_ANS+
      t2.END_END+
      t2.END_LIECX==1
    
    
  }
  
  #les LIEC_X identifi?s : que sont ils r?ellement ?
  {
    
    #nombre id correcte EPI
    n2.LIECX_EPI=sum(tabNP2015$Nbr_EW[tabNP2015$FIEC_LIEC%in%c("LIEC_X_EPI")])
    
    #nombre id incorrecte EPA 
    n2.LIECX_EPA=sum(tabNP2015$Nbr_EW[tabNP2015$FIEC_LIEC%in%c("LIEC_X_EPA")])
    
    #nombre id incorrecte ANS 
    n2.LIECX_ANS=sum(tabNP2015$Nbr_EW[tabNP2015$FIEC_LIEC%in%c("LIEC_X_ANS")])
    
    #nombre id incorrecte END 
    n2.LIECX_END=sum(tabNP2015$Nbr_EW[tabNP2015$FIEC_LIEC%in%c("LIEC_X_END")])
    
    #nombre id incorrecte LIEC_X 
    n2.LIECX_LIECX=sum(tabNP2015$Nbr_EW[tabNP2015$FIEC_LIEC%in%c("LIEC_X_LIEC_X")])
    
    
    
    
    
    #taux id correcte EPI
    t2.LIECX_EPI=sum(tabNP2015$Nbr_EW[tabNP2015$FIEC_LIEC%in%c("LIEC_X_EPI")])/
      sum(tabNP2015$Nbr_EW[tabNP2015$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
    
    
    
    
    #taux id incorrecte EPA 
    t2.LIECX_EPA=sum(tabNP2015$Nbr_EW[tabNP2015$FIEC_LIEC%in%c("LIEC_X_EPA")])/
      sum(tabNP2015$Nbr_EW[tabNP2015$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
    
    
    
    
    #taux id incorrecte ANS 
    t2.LIECX_ANS=sum(tabNP2015$Nbr_EW[tabNP2015$FIEC_LIEC%in%c("LIEC_X_ANS")])/
      sum(tabNP2015$Nbr_EW[tabNP2015$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
    
    
    
    
    
    #taux id incorrecte END 
    t2.LIECX_END=sum(tabNP2015$Nbr_EW[tabNP2015$FIEC_LIEC%in%c("LIEC_X_END")])/
      sum(tabNP2015$Nbr_EW[tabNP2015$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
    
    
    
    
    #taux id incorrecte LIEC_X 
    t2.LIECX_LIECX=sum(tabNP2015$Nbr_EW[tabNP2015$FIEC_LIEC%in%c("LIEC_X_LIEC_X")])/
      sum(tabNP2015$Nbr_EW[tabNP2015$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
    
    
    
    t2.LIECX_EPI+
      t2.LIECX_EPA+
      t2.LIECX_ANS+
      t2.LIECX_END+
      t2.LIECX_LIECX==1
    
    
  }
  
  
  
  tabGraphNP15_2=data.frame(
    Nbr=c(
      n2.EPI_EPI,
      n2.EPI_EPA,
      n2.EPI_ANS,
      n2.EPI_END,
      n2.EPI_LIECX,
      n2.EPA_EPI,
      n2.EPA_EPA,
      n2.EPA_ANS,
      n2.EPA_END,
      n2.EPA_LIECX,
      n2.ANS_EPI,
      n2.ANS_EPA,
      n2.ANS_ANS,
      n2.ANS_END,
      n2.ANS_LIECX,
      n2.END_EPI,
      n2.END_EPA,
      n2.END_ANS,
      n2.END_END,
      n2.END_LIECX,
      n2.LIECX_EPI,
      n2.LIECX_EPA,
      n2.LIECX_ANS,
      n2.LIECX_END,
      n2.LIECX_LIECX),
    
    Taux=c(
      t2.EPI_EPI,
      t2.EPI_EPA,
      t2.EPI_ANS,
      t2.EPI_END,
      t2.EPI_LIECX,
      t2.EPA_EPI,
      t2.EPA_EPA,
      t2.EPA_ANS,
      t2.EPA_END,
      t2.EPA_LIECX,
      t2.ANS_EPI,
      t2.ANS_EPA,
      t2.ANS_ANS,
      t2.ANS_END,
      t2.ANS_LIECX,
      t2.END_EPI,
      t2.END_EPA,
      t2.END_ANS,
      t2.END_END,
      t2.END_LIECX,
      t2.LIECX_EPI,
      t2.LIECX_EPA,
      t2.LIECX_ANS,
      t2.LIECX_END,
      t2.LIECX_LIECX),
    
    Classe=c(
      "t2.EPI_EPI",
      "t2.EPI_EPA",
      "t2.EPI_ANS",
      "t2.EPI_END",
      "t2.EPI_LIECX",
      "t2.EPA_EPI",
      "t2.EPA_EPA",
      "t2.EPA_ANS",
      "t2.EPA_END",
      "t2.EPA_LIECX",
      "t2.ANS_EPI",
      "t2.ANS_EPA",
      "t2.ANS_ANS",
      "t2.ANS_END",
      "t2.ANS_LIECX",
      "t2.END_EPI",
      "t2.END_EPA",
      "t2.END_ANS",
      "t2.END_END",
      "t2.END_LIECX",
      "t2.LIECX_EPI",
      "t2.LIECX_EPA",
      "t2.LIECX_ANS",
      "t2.LIECX_END",
      "t2.LIECX_LIECX")
    
  )
  
  
  tabGraphNP15_2$LIEC=factor(substr(tabGraphNP15_2$Classe,8,10),levels=c("EPI","EPA","ANS","END","LIECX"))
  tabGraphNP15_2$LIECid=factor(substr(tabGraphNP15_2$Classe,4,6),levels=c("EPI","EPA","ANS","END","LIECX"))
  anti_join(tabGraphNP15_2,ddply(tabGraphNP15_2, .(Classe), mutate, Pos_Y=cumsum(Nbr)-Nbr/2 ))
  anti_join(ddply(tabGraphNP15_2, .(Classe), mutate, Pos_Y=cumsum(Nbr)-Nbr/2 ),tabGraphNP15_2)
  tabGraphNP15_2 <- ddply(tabGraphNP15_2, .(Classe), mutate, Pos_Y=cumsum(Nbr)-Nbr/2 )
  pos_X <- NULL
  for(j in 1:nlevels(tabGraphNP15_2$Classe) ){
    # Pour centrer les ?tiqu?tes en fonction du nombre de Classe par LIEC
    ss_tab_mod <- subset(tabGraphNP15_2, Classe==levels(tabGraphNP15_2$Classe)[j] )
    
    # Si la Classe a des vdt
    if( sum(ss_tab_mod$Nbr,na.rm=T)>0 ){
      ss_tab_mod_bon <- ss_tab_mod[ss_tab_mod$Nbr!=0,] #supprime ligne vide pour Nbr
      ss_tab_mod_bon$LIEC <- factor(ss_tab_mod_bon$LIEC ) # ?quivaut ? un droplevels()
      pos_X2 <- NULL#position etiquette
      
      if( nlevels(ss_tab_mod_bon$LIEC)==1 ){ pos_X2 <- j #un LIEC : au centre 
      }else{
        if( Disper_LIEC==0 ){#pour d?caler les etiquettes
          pos_X2 <- rev(j+seq(-0.4,0.4,0.8/(nlevels(ss_tab_mod_bon$LIEC)-1) ) )
        }else{#pour position de l'ecart
          pos_X2 <- rev(j+seq(-0.3,0.3,0.6/(nlevels(ss_tab_mod_bon$LIEC)-1) ) )
        }
      }
      
      if( "AB_NCX" %in% colnames(data_VDT_input) ){
        pos_X2 <- c(pos_X2,  rep(NA, 6-length(pos_X2) ) )
      }else{
        pos_X2 <- c(pos_X2,  rep(NA, 5-length(pos_X2) ) )
      }
      
      # Si la Classe n'a pas de vdt 
    }else{
      ss_tab_mod_bon <- NULL
      if( "AB_NCX" %in% colnames(data_VDT_input) ){
        pos_X2 <- rep(NA, 6) 
      }else{
        pos_X2 <- rep(NA, 5) 
      }
    }
    
    if( "AB_NCX" %in% colnames(data_VDT_input) ){
      names(pos_X2) <- c( levels(ss_tab_mod_bon$LIEC), setdiff( c("LIEC_X", "END", "AnS_END", "ANS", "EPA", "EPI"), 
                                                              levels(ss_tab_mod_bon$LIEC)) )
    }else{
      names(pos_X2) <- c( levels(ss_tab_mod_bon$LIEC), setdiff( c("LIEC_X","END", "ANS", "EPA", "EPI"), 
                                                              levels(ss_tab_mod_bon$LIEC)) )
    }
    #remet dans l'odre des levels
    pos_X2 <- pos_X2[levels(tabGraphNP15_2$LIEC)]
    pos_X <- c( pos_X,  pos_X2 )
  }
  tabGraphNP15_2$pos_X <- pos_X
  
  
  
  
  
  
  
  ggplot(tabGraphNP15_2, aes(fill=LIEC, y=Taux, x=LIECid)) + 
    geom_bar(position='stack', stat='identity')+
    theme_classic()+
    scale_fill_manual(values = c("#000099","#FFFF00","#009900","#FF6699","#7B7B7B"))
  
  ggplot(tabGraphNP15_2, aes(fill=LIEC, y=Taux, x=LIECid)) + 
    geom_bar(position='fill', stat='identity')+
    theme_classic()+
    scale_fill_manual(values = c("#000099","#FFFF00","#009900","#FF6699","#7B7B7B"))
  
  ggplot(tabGraphNP15_2, aes(fill=LIEC, y=Nbr, x=LIECid)) + 
    geom_bar(position='fill', stat='identity')+
    theme_classic()+
    scale_fill_manual(values = c("#000099","#FFFF00","#009900","#FF6699","#7B7B7B"))
  
  ggplot(tabGraphNP15_2, aes(fill=LIEC, y=Nbr, x=LIECid)) + 
    geom_bar(position='stack', stat='identity')+
    theme_classic()+
    scale_fill_manual(values = c("#000099","#FFFF00","#009900","#FF6699","#7B7B7B"))
  
  
  
  
  
  
  
  #zoom sur les erreurs
  tabGraphNP15_2$Erreur=NA
  tabGraphNP15_2$Erreur[tabGraphNP15_2$LIEC==tabGraphNP15_2$LIECid]="Correcte"
  tabGraphNP15_2$Erreur[tabGraphNP15_2$LIEC!=tabGraphNP15_2$LIECid]="Incorrecte"
  
  
  
  ggplot(tabGraphNP15_2[tabGraphNP15_2$Erreur=="Incorrecte" & tabGraphNP15_2$LIECid!="LIECX",], 
         aes(fill=LIEC, y=Taux, x=LIECid)) + 
    geom_bar(position='stack', stat='identity')+
    theme_classic()+
    theme(title=element_text(face="bold", size=15),
          axis.text=element_text(face="bold", size=15),
          axis.title.x = element_text(face="bold", vjust=0),#size=25
          # axis.ticks.margin=unit(c(0.5),'cm'),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_line(size=2),
          axis.ticks.length=unit(0.5, "cm"),
          axis.text.x=element_blank(),
          panel.grid=element_blank(), # Augmenter taille ligne y
          panel.grid.major.x=element_blank(), # Suppression lignes axes x
          axis.title=element_text( vjust=5, size=15),
          line=element_line(size=1.5),
          legend.position="none",
          legend.title = element_text(face="bold", size=15),
          legend.text = element_text(size=15)) +
    scale_y_continuous(labels = scales::percent, limits=c(0,1)) + 
    scale_fill_manual(values = c("#000099","#FFFF00","#009900","#FF6699","#7B7B7B"))+
    xlab(NULL)+ylab(NULL)
  
  
  
  ggplot(tabGraphNP15_2[tabGraphNP15_2$Erreur=="Incorrecte" & tabGraphNP15_2$LIECid!="LIECX",], 
         aes(fill=LIEC, y=Taux, x=LIECid)) + 
    geom_bar(position='stack', stat='identity')+
    theme_classic()+
    theme(title=element_text(face="bold", size=15),
          axis.text=element_text(face="bold", size=15),
          axis.title.x = element_text(face="bold", vjust=0),#size=25
          # axis.ticks.margin=unit(c(0.5),'cm'),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_line(size=2),
          axis.ticks.length=unit(0.5, "cm"),
          axis.text.x=element_blank(),
          panel.grid=element_blank(), # Augmenter taille ligne y
          panel.grid.major.x=element_blank(), # Suppression lignes axes x
          axis.title=element_text( vjust=5, size=15),
          line=element_line(size=1.5),
          legend.position="none",
          legend.title = element_text(face="bold", size=15),
          legend.text = element_text(size=15)) +
    scale_y_continuous(labels = scales::percent, limits=c(0,0.15)) + 
    scale_fill_manual(values = c("#000099","#FFFF00","#009900","#FF6699","#7B7B7B"))+
    xlab(NULL)+ylab(NULL)
  
  
  
  ggplot(tabGraphNP15_2[tabGraphNP15_2$Erreur=="Incorrecte" & tabGraphNP15_2$LIECid!="LIECX",], 
         aes(fill=LIEC, y=Nbr, x=LIECid)) + 
    geom_bar(position='stack', stat='identity')+
    theme_classic()+
    theme(title=element_text(face="bold", size=15),
          axis.text=element_text(face="bold", size=15),
          axis.title.x = element_text(face="bold", vjust=0),#size=25
          # axis.ticks.margin=unit(c(0.5),'cm'),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_line(size=2),
          axis.ticks.length=unit(0.5, "cm"),
          axis.text.x=element_blank(),
          panel.grid=element_blank(), # Augmenter taille ligne y
          panel.grid.major.x=element_blank(), # Suppression lignes axes x
          axis.title=element_text( vjust=5, size=15),
          line=element_line(size=1.5),
          legend.position="none",
          legend.title = element_text(face="bold", size=15),
          legend.text = element_text(size=15)) +
    scale_fill_manual(values = c("#000099","#FFFF00","#009900","#FF6699","#7B7B7B"))+
    xlab(NULL)+ylab(NULL)
  
  
  ggplot(tabGraphNP15_2[tabGraphNP15_2$LIECid!="LIECX",], 
         aes(fill=LIEC, y=Nbr, x=LIECid)) + 
    geom_bar(position='stack', stat='identity')+
    theme_classic()+
    theme(title=element_text(face="bold", size=15),
          axis.text=element_text(face="bold", size=15),
          axis.title.x = element_text(face="bold", vjust=0),#size=25
          # axis.ticks.margin=unit(c(0.5),'cm'),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_line(size=2),
          axis.ticks.length=unit(0.5, "cm"),
          axis.text.x=element_blank(),
          panel.grid=element_blank(), # Augmenter taille ligne y
          panel.grid.major.x=element_blank(), # Suppression lignes axes x
          axis.title=element_text( vjust=5, size=15),
          line=element_line(size=1.5),
          legend.position="none",
          legend.title = element_text(face="bold", size=15),
          legend.text = element_text(size=15)) +
    scale_fill_manual(values = c("#000099","#FFFF00","#009900","#FF6699","#7B7B7B"))+
    xlab(NULL)+ylab(NULL)
  
  
  geom_point(aes(y=replace(Pos_Y,which(AB==0),NA), x=pos_X ),
             color="grey", size=Size_Rond, show.legend=F) + # Point pour entourer le nb de VDT
    geom_point(aes(y=replace(Pos_Y, which(AB==0),NA),  x=pos_X, fill=LIEC, shape=LIEC), 
               color="grey", size=Size_Rond+1, show.legend=F) + # Point pour entourer le nb de VDT
    geom_text(aes(y=replace(Pos_Y, which(AB==0),NA), colour=LIEC, x=pos_X), label=label_AB$AB, # Affichage valeur Ab >0
              size=Size_Eti, fontface=2, show.legend=F)    # Augmenter taille et font (gras)
  
  
  
##### BON GROUPE NP 2014 #####
  
  tabNP=Tab2[Tab2$Region%in%"NP",]
  nrow(tabNP)==1208
  dataID <- unique(sbt_m[,c("cle",colID,"Region")])
  nrow(dataID)==1109
  dataID=dataID[dataID$Region%in%"NP",]
  nrow(dataID)==78
  
  
  
  
  tabNP2014=tabNP[tabNP$Annee%in%2014,]
  nrow(tabNP2014)==183
  
  #taux d'id correcte global
  sum(tabNP2014$Nbr_EW[tabNP2014$FIEC_LIEC%in%c("EPI_EPI","EPA_EPA","ANS_ANS","END_END")])/
    sum(tabNP2014$Nbr_EW)#0.9615385  => ok
  #Tab_reu_Final2$tauxErr[Tab_reu_Final2$Region=="NP" & Tab_reu_Final2$Annee==2014]
  
  
  #les EPI identifi?s : que sont ils r?ellement ?
  {
    #nombre id correcte EPI
    n2.EPI_EPI=sum(tabNP2014$Nbr_EW[tabNP2014$FIEC_LIEC%in%c("EPI_EPI")])
    
    #nombre id incorrecte EPA 
    n2.EPI_EPA=sum(tabNP2014$Nbr_EW[tabNP2014$FIEC_LIEC%in%c("EPI_EPA")])
    
    #nombre id incorrecte ANS 
    n2.EPI_ANS=sum(tabNP2014$Nbr_EW[tabNP2014$FIEC_LIEC%in%c("EPI_ANS")])
    
    #nombre id incorrecte END 
    n2.EPI_END=sum(tabNP2014$Nbr_EW[tabNP2014$FIEC_LIEC%in%c("EPI_END")])
    
    #nombre id incorrecte LIEC_X 
    n2.EPI_LIECX=sum(tabNP2014$Nbr_EW[tabNP2014$FIEC_LIEC%in%c("EPI_LIEC_X")])
    
    
    
    
    #taux id correcte EPI
    t2.EPI_EPI=sum(tabNP2014$Nbr_EW[tabNP2014$FIEC_LIEC%in%c("EPI_EPI")])/
      sum(tabNP2014$Nbr_EW[tabNP2014$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
    
    #taux id incorrecte EPA 
    t2.EPI_EPA=sum(tabNP2014$Nbr_EW[tabNP2014$FIEC_LIEC%in%c("EPI_EPA")])/
      sum(tabNP2014$Nbr_EW[tabNP2014$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
    
    #taux id incorrecte ANS 
    t2.EPI_ANS=sum(tabNP2014$Nbr_EW[tabNP2014$FIEC_LIEC%in%c("EPI_ANS")])/
      sum(tabNP2014$Nbr_EW[tabNP2014$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
    
    #taux id incorrecte END 
    t2.EPI_END=sum(tabNP2014$Nbr_EW[tabNP2014$FIEC_LIEC%in%c("EPI_END")])/
      sum(tabNP2014$Nbr_EW[tabNP2014$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
    
    #taux id incorrecte LIEC_X 
    t2.EPI_LIECX=sum(tabNP2014$Nbr_EW[tabNP2014$FIEC_LIEC%in%c("EPI_LIEC_X")])/
      sum(tabNP2014$Nbr_EW[tabNP2014$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
    
    t2.EPI_EPI+
      t2.EPI_EPA+
      t2.EPI_ANS+
      t2.EPI_END+
      t2.EPI_LIECX==1
    
    
    }
  
  #les EPA identifi?s : que sont ils r?ellement ?
  {
    
    #nombre id correcte EPI
    n2.EPA_EPI=sum(tabNP2014$Nbr_EW[tabNP2014$FIEC_LIEC%in%c("EPA_EPI")])
    
    #nombre id incorrecte EPA 
    n2.EPA_EPA=sum(tabNP2014$Nbr_EW[tabNP2014$FIEC_LIEC%in%c("EPA_EPA")])
    
    #nombre id incorrecte ANS 
    n2.EPA_ANS=sum(tabNP2014$Nbr_EW[tabNP2014$FIEC_LIEC%in%c("EPA_ANS")])
    
    #nombre id incorrecte END 
    n2.EPA_END=sum(tabNP2014$Nbr_EW[tabNP2014$FIEC_LIEC%in%c("EPA_END")])
    
    #nombre id incorrecte LIEC_X 
    n2.EPA_LIECX=sum(tabNP2014$Nbr_EW[tabNP2014$FIEC_LIEC%in%c("EPA_LIEC_X")])
    
    
    
    
    
    
    #taux id correcte EPI
    t2.EPA_EPI=sum(tabNP2014$Nbr_EW[tabNP2014$FIEC_LIEC%in%c("EPA_EPI")])/
      sum(tabNP2014$Nbr_EW[tabNP2014$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
    
    
    
    
    #taux id incorrecte EPA 
    t2.EPA_EPA=sum(tabNP2014$Nbr_EW[tabNP2014$FIEC_LIEC%in%c("EPA_EPA")])/
      sum(tabNP2014$Nbr_EW[tabNP2014$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
    
    
    
    
    #taux id incorrecte ANS 
    t2.EPA_ANS=sum(tabNP2014$Nbr_EW[tabNP2014$FIEC_LIEC%in%c("EPA_ANS")])/
      sum(tabNP2014$Nbr_EW[tabNP2014$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
    
    
    
    
    
    #taux id incorrecte END 
    t2.EPA_END=sum(tabNP2014$Nbr_EW[tabNP2014$FIEC_LIEC%in%c("EPA_END")])/
      sum(tabNP2014$Nbr_EW[tabNP2014$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
    
    
    
    
    #taux id incorrecte LIEC_X 
    t2.EPA_LIECX=sum(tabNP2014$Nbr_EW[tabNP2014$FIEC_LIEC%in%c("EPA_LIEC_X")])/
      sum(tabNP2014$Nbr_EW[tabNP2014$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
    
    
    
    t2.EPA_EPI+
      t2.EPA_EPA+
      t2.EPA_ANS+
      t2.EPA_END+
      t2.EPA_LIECX==1
    
    
  }
  
  #les ANS identifi?s : que sont ils r?ellement ?
  {
    
    #nombre id correcte EPI
    n2.ANS_EPI=sum(tabNP2014$Nbr_EW[tabNP2014$FIEC_LIEC%in%c("ANS_EPI")])
    
    #nombre id incorrecte EPA 
    n2.ANS_EPA=sum(tabNP2014$Nbr_EW[tabNP2014$FIEC_LIEC%in%c("ANS_EPA")])
    
    #nombre id incorrecte ANS 
    n2.ANS_ANS=sum(tabNP2014$Nbr_EW[tabNP2014$FIEC_LIEC%in%c("ANS_ANS")])
    
    #nombre id incorrecte END 
    n2.ANS_END=sum(tabNP2014$Nbr_EW[tabNP2014$FIEC_LIEC%in%c("ANS_END")])
    
    #nombre id incorrecte LIEC_X 
    n2.ANS_LIECX=sum(tabNP2014$Nbr_EW[tabNP2014$FIEC_LIEC%in%c("ANS_LIEC_X")])
    
    
    
    
    
    
    #taux id correcte EPI
    t2.ANS_EPI=sum(tabNP2014$Nbr_EW[tabNP2014$FIEC_LIEC%in%c("ANS_EPI")])/
      sum(tabNP2014$Nbr_EW[tabNP2014$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
    
    
    
    
    #taux id incorrecte EPA 
    t2.ANS_EPA=sum(tabNP2014$Nbr_EW[tabNP2014$FIEC_LIEC%in%c("ANS_EPA")])/
      sum(tabNP2014$Nbr_EW[tabNP2014$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
    
    
    
    
    #taux id incorrecte ANS 
    t2.ANS_ANS=sum(tabNP2014$Nbr_EW[tabNP2014$FIEC_LIEC%in%c("ANS_ANS")])/
      sum(tabNP2014$Nbr_EW[tabNP2014$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
    
    
    
    
    
    #taux id incorrecte END 
    t2.ANS_END=sum(tabNP2014$Nbr_EW[tabNP2014$FIEC_LIEC%in%c("ANS_END")])/
      sum(tabNP2014$Nbr_EW[tabNP2014$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
    
    
    
    
    #taux id incorrecte LIEC_X 
    t2.ANS_LIECX=sum(tabNP2014$Nbr_EW[tabNP2014$FIEC_LIEC%in%c("ANS_LIEC_X")])/
      sum(tabNP2014$Nbr_EW[tabNP2014$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
    
    
    
    t2.ANS_EPI+
      t2.ANS_EPA+
      t2.ANS_ANS+
      t2.ANS_END+
      t2.ANS_LIECX==1
    
    
  }
  
  #les END identifi?s : que sont ils r?ellement ?
  {
    
    #nombre id correcte EPI
    n2.END_EPI=sum(tabNP2014$Nbr_EW[tabNP2014$FIEC_LIEC%in%c("END_EPI")])
    
    #nombre id incorrecte EPA 
    n2.END_EPA=sum(tabNP2014$Nbr_EW[tabNP2014$FIEC_LIEC%in%c("END_EPA")])
    
    #nombre id incorrecte ANS 
    n2.END_ANS=sum(tabNP2014$Nbr_EW[tabNP2014$FIEC_LIEC%in%c("END_ANS")])
    
    #nombre id incorrecte END 
    n2.END_END=sum(tabNP2014$Nbr_EW[tabNP2014$FIEC_LIEC%in%c("END_END")])
    
    #nombre id incorrecte LIEC_X 
    n2.END_LIECX=sum(tabNP2014$Nbr_EW[tabNP2014$FIEC_LIEC%in%c("END_LIEC_X")])
    
    
    
    
    
    
    #taux id correcte EPI
    t2.END_EPI=sum(tabNP2014$Nbr_EW[tabNP2014$FIEC_LIEC%in%c("END_EPI")])/
      sum(tabNP2014$Nbr_EW[tabNP2014$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
    
    
    
    
    #taux id incorrecte EPA 
    t2.END_EPA=sum(tabNP2014$Nbr_EW[tabNP2014$FIEC_LIEC%in%c("END_EPA")])/
      sum(tabNP2014$Nbr_EW[tabNP2014$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
    
    
    
    
    #taux id incorrecte ANS 
    t2.END_ANS=sum(tabNP2014$Nbr_EW[tabNP2014$FIEC_LIEC%in%c("END_ANS")])/
      sum(tabNP2014$Nbr_EW[tabNP2014$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
    
    
    
    
    
    #taux id incorrecte END 
    t2.END_END=sum(tabNP2014$Nbr_EW[tabNP2014$FIEC_LIEC%in%c("END_END")])/
      sum(tabNP2014$Nbr_EW[tabNP2014$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
    
    
    
    
    #taux id incorrecte LIEC_X 
    t2.END_LIECX=sum(tabNP2014$Nbr_EW[tabNP2014$FIEC_LIEC%in%c("END_LIEC_X")])/
      sum(tabNP2014$Nbr_EW[tabNP2014$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
    
    
    
    round(t2.END_EPI+
      t2.END_EPA+
      t2.END_ANS+
      t2.END_END+
      t2.END_LIECX,6)==1
    
    
  }
  
  #les LIEC_X identifi?s : que sont ils r?ellement ?
  {
    
    #nombre id correcte EPI
    n2.LIECX_EPI=sum(tabNP2014$Nbr_EW[tabNP2014$FIEC_LIEC%in%c("LIEC_X_EPI")])
    
    #nombre id incorrecte EPA 
    n2.LIECX_EPA=sum(tabNP2014$Nbr_EW[tabNP2014$FIEC_LIEC%in%c("LIEC_X_EPA")])
    
    #nombre id incorrecte ANS 
    n2.LIECX_ANS=sum(tabNP2014$Nbr_EW[tabNP2014$FIEC_LIEC%in%c("LIEC_X_ANS")])
    
    #nombre id incorrecte END 
    n2.LIECX_END=sum(tabNP2014$Nbr_EW[tabNP2014$FIEC_LIEC%in%c("LIEC_X_END")])
    
    #nombre id incorrecte LIEC_X 
    n2.LIECX_LIECX=sum(tabNP2014$Nbr_EW[tabNP2014$FIEC_LIEC%in%c("LIEC_X_LIEC_X")])
    
    
    
    
    
    #taux id correcte EPI
    t2.LIECX_EPI=sum(tabNP2014$Nbr_EW[tabNP2014$FIEC_LIEC%in%c("LIEC_X_EPI")])/
      sum(tabNP2014$Nbr_EW[tabNP2014$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
    
    
    
    
    #taux id incorrecte EPA 
    t2.LIECX_EPA=sum(tabNP2014$Nbr_EW[tabNP2014$FIEC_LIEC%in%c("LIEC_X_EPA")])/
      sum(tabNP2014$Nbr_EW[tabNP2014$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
    
    
    
    
    #taux id incorrecte ANS 
    t2.LIECX_ANS=sum(tabNP2014$Nbr_EW[tabNP2014$FIEC_LIEC%in%c("LIEC_X_ANS")])/
      sum(tabNP2014$Nbr_EW[tabNP2014$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
    
    
    
    
    
    #taux id incorrecte END 
    t2.LIECX_END=sum(tabNP2014$Nbr_EW[tabNP2014$FIEC_LIEC%in%c("LIEC_X_END")])/
      sum(tabNP2014$Nbr_EW[tabNP2014$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
    
    
    
    
    #taux id incorrecte LIEC_X 
    t2.LIECX_LIECX=sum(tabNP2014$Nbr_EW[tabNP2014$FIEC_LIEC%in%c("LIEC_X_LIEC_X")])/
      sum(tabNP2014$Nbr_EW[tabNP2014$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
    
    
    
    t2.LIECX_EPI+
      t2.LIECX_EPA+
      t2.LIECX_ANS+
      t2.LIECX_END+
      t2.LIECX_LIECX==1
    
    
  }
  
  
  
  tabGraphNP14_2=data.frame(
    Nbr=c(
      n2.EPI_EPI,
      n2.EPI_EPA,
      n2.EPI_ANS,
      n2.EPI_END,
      n2.EPI_LIECX,
      n2.EPA_EPI,
      n2.EPA_EPA,
      n2.EPA_ANS,
      n2.EPA_END,
      n2.EPA_LIECX,
      n2.ANS_EPI,
      n2.ANS_EPA,
      n2.ANS_ANS,
      n2.ANS_END,
      n2.ANS_LIECX,
      n2.END_EPI,
      n2.END_EPA,
      n2.END_ANS,
      n2.END_END,
      n2.END_LIECX,
      n2.LIECX_EPI,
      n2.LIECX_EPA,
      n2.LIECX_ANS,
      n2.LIECX_END,
      n2.LIECX_LIECX),
    
    Taux=c(
      t2.EPI_EPI,
      t2.EPI_EPA,
      t2.EPI_ANS,
      t2.EPI_END,
      t2.EPI_LIECX,
      t2.EPA_EPI,
      t2.EPA_EPA,
      t2.EPA_ANS,
      t2.EPA_END,
      t2.EPA_LIECX,
      t2.ANS_EPI,
      t2.ANS_EPA,
      t2.ANS_ANS,
      t2.ANS_END,
      t2.ANS_LIECX,
      t2.END_EPI,
      t2.END_EPA,
      t2.END_ANS,
      t2.END_END,
      t2.END_LIECX,
      t2.LIECX_EPI,
      t2.LIECX_EPA,
      t2.LIECX_ANS,
      t2.LIECX_END,
      t2.LIECX_LIECX),
    
    Classe=c(
      "t2.EPI_EPI",
      "t2.EPI_EPA",
      "t2.EPI_ANS",
      "t2.EPI_END",
      "t2.EPI_LIECX",
      "t2.EPA_EPI",
      "t2.EPA_EPA",
      "t2.EPA_ANS",
      "t2.EPA_END",
      "t2.EPA_LIECX",
      "t2.ANS_EPI",
      "t2.ANS_EPA",
      "t2.ANS_ANS",
      "t2.ANS_END",
      "t2.ANS_LIECX",
      "t2.END_EPI",
      "t2.END_EPA",
      "t2.END_ANS",
      "t2.END_END",
      "t2.END_LIECX",
      "t2.LIECX_EPI",
      "t2.LIECX_EPA",
      "t2.LIECX_ANS",
      "t2.LIECX_END",
      "t2.LIECX_LIECX")
    
  )
  
  
  tabGraphNP14_2$LIEC=factor(substr(tabGraphNP14_2$Classe,8,10),levels=c("EPI","EPA","ANS","END","LIECX"))
  tabGraphNP14_2$LIECid=factor(substr(tabGraphNP14_2$Classe,4,6),levels=c("EPI","EPA","ANS","END","LIECX"))
  anti_join(tabGraphNP14_2,ddply(tabGraphNP14_2, .(Classe), mutate, Pos_Y=cumsum(Nbr)-Nbr/2 ))
  anti_join(ddply(tabGraphNP14_2, .(Classe), mutate, Pos_Y=cumsum(Nbr)-Nbr/2 ),tabGraphNP14_2)
  tabGraphNP14_2 <- ddply(tabGraphNP14_2, .(Classe), mutate, Pos_Y=cumsum(Nbr)-Nbr/2 )
  pos_X <- NULL
  for(j in 1:nlevels(tabGraphNP14_2$Classe) ){
    # Pour centrer les ?tiqu?tes en fonction du nombre de Classe par LIEC
    ss_tab_mod <- subset(tabGraphNP14_2, Classe==levels(tabGraphNP14_2$Classe)[j] )
    
    # Si la Classe a des vdt
    if( sum(ss_tab_mod$Nbr,na.rm=T)>0 ){
      ss_tab_mod_bon <- ss_tab_mod[ss_tab_mod$Nbr!=0,] #supprime ligne vide pour Nbr
      ss_tab_mod_bon$LIEC <- factor(ss_tab_mod_bon$LIEC ) # ?quivaut ? un droplevels()
      pos_X2 <- NULL#position etiquette
      
      if( nlevels(ss_tab_mod_bon$LIEC)==1 ){ pos_X2 <- j #un LIEC : au centre 
      }else{
        if( Disper_LIEC==0 ){#pour d?caler les etiquettes
          pos_X2 <- rev(j+seq(-0.4,0.4,0.8/(nlevels(ss_tab_mod_bon$LIEC)-1) ) )
        }else{#pour position de l'ecart
          pos_X2 <- rev(j+seq(-0.3,0.3,0.6/(nlevels(ss_tab_mod_bon$LIEC)-1) ) )
        }
      }
      
      if( "AB_NCX" %in% colnames(data_VDT_input) ){
        pos_X2 <- c(pos_X2,  rep(NA, 6-length(pos_X2) ) )
      }else{
        pos_X2 <- c(pos_X2,  rep(NA, 5-length(pos_X2) ) )
      }
      
      # Si la Classe n'a pas de vdt 
    }else{
      ss_tab_mod_bon <- NULL
      if( "AB_NCX" %in% colnames(data_VDT_input) ){
        pos_X2 <- rep(NA, 6) 
      }else{
        pos_X2 <- rep(NA, 5) 
      }
    }
    
    if( "AB_NCX" %in% colnames(data_VDT_input) ){
      names(pos_X2) <- c( levels(ss_tab_mod_bon$LIEC), setdiff( c("LIEC_X", "END", "AnS_END", "ANS", "EPA", "EPI"), 
                                                              levels(ss_tab_mod_bon$LIEC)) )
    }else{
      names(pos_X2) <- c( levels(ss_tab_mod_bon$LIEC), setdiff( c("LIEC_X","END", "ANS", "EPA", "EPI"), 
                                                              levels(ss_tab_mod_bon$LIEC)) )
    }
    #remet dans l'odre des levels
    pos_X2 <- pos_X2[levels(tabGraphNP14_2$LIEC)]
    pos_X <- c( pos_X,  pos_X2 )
  }
  tabGraphNP14_2$pos_X <- pos_X
  
  
  
  
  
  
  
  ggplot(tabGraphNP14_2, aes(fill=LIEC, y=Taux, x=LIECid)) + 
    geom_bar(position='stack', stat='identity')+
    theme_classic()+
    scale_fill_manual(values = c("#000099","#FFFF00","#009900","#FF6699","#7B7B7B"))
  
  ggplot(tabGraphNP14_2, aes(fill=LIEC, y=Taux, x=LIECid)) + 
    geom_bar(position='fill', stat='identity')+
    theme_classic()+
    scale_fill_manual(values = c("#000099","#FFFF00","#009900","#FF6699","#7B7B7B"))
  
  ggplot(tabGraphNP14_2, aes(fill=LIEC, y=Nbr, x=LIECid)) + 
    geom_bar(position='fill', stat='identity')+
    theme_classic()+
    scale_fill_manual(values = c("#000099","#FFFF00","#009900","#FF6699","#7B7B7B"))
  
  ggplot(tabGraphNP14_2, aes(fill=LIEC, y=Nbr, x=LIECid)) + 
    geom_bar(position='stack', stat='identity')+
    theme_classic()+
    scale_fill_manual(values = c("#000099","#FFFF00","#009900","#FF6699","#7B7B7B"))
  
  
  
  
  
  
  
  #zoom sur les erreurs
  tabGraphNP14_2$Erreur=NA
  tabGraphNP14_2$Erreur[tabGraphNP14_2$LIEC==tabGraphNP14_2$LIECid]="Correcte"
  tabGraphNP14_2$Erreur[tabGraphNP14_2$LIEC!=tabGraphNP14_2$LIECid]="Incorrecte"
  
  
  
  ggplot(tabGraphNP14_2[tabGraphNP14_2$Erreur=="Incorrecte" & tabGraphNP14_2$LIECid!="LIECX",], 
         aes(fill=LIEC, y=Taux, x=LIECid)) + 
    geom_bar(position='stack', stat='identity')+
    theme_classic()+
    theme(title=element_text(face="bold", size=15),
          axis.text=element_text(face="bold", size=15),
          axis.title.x = element_text(face="bold", vjust=0),#size=25
          # axis.ticks.margin=unit(c(0.5),'cm'),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_line(size=2),
          axis.ticks.length=unit(0.5, "cm"),
          axis.text.x=element_blank(),
          panel.grid=element_blank(), # Augmenter taille ligne y
          panel.grid.major.x=element_blank(), # Suppression lignes axes x
          axis.title=element_text( vjust=5, size=15),
          line=element_line(size=1.5),
          legend.position="none",
          legend.title = element_text(face="bold", size=15),
          legend.text = element_text(size=15)) +
    scale_y_continuous(labels = scales::percent, limits=c(0,1)) + 
    scale_fill_manual(values = c("#000099","#FFFF00","#009900","#FF6699","#7B7B7B"))+
    xlab(NULL)+ylab(NULL)
  
  
  
  ggplot(tabGraphNP14_2[tabGraphNP14_2$Erreur=="Incorrecte" & tabGraphNP14_2$LIECid!="LIECX",], 
         aes(fill=LIEC, y=Taux, x=LIECid)) + 
    geom_bar(position='stack', stat='identity')+
    theme_classic()+
    theme(title=element_text(face="bold", size=15),
          axis.text=element_text(face="bold", size=15),
          axis.title.x = element_text(face="bold", vjust=0),#size=25
          # axis.ticks.margin=unit(c(0.5),'cm'),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_line(size=2),
          axis.ticks.length=unit(0.5, "cm"),
          axis.text.x=element_blank(),
          panel.grid=element_blank(), # Augmenter taille ligne y
          panel.grid.major.x=element_blank(), # Suppression lignes axes x
          axis.title=element_text( vjust=5, size=15),
          line=element_line(size=1.5),
          legend.position="none",
          legend.title = element_text(face="bold", size=15),
          legend.text = element_text(size=15)) +
    scale_y_continuous(labels = scales::percent, limits=c(0,0.25)) + 
    scale_fill_manual(values = c("#000099","#FFFF00","#009900","#FF6699","#7B7B7B"))+
    xlab(NULL)+ylab(NULL)
  
  
  
  ggplot(tabGraphNP14_2[tabGraphNP14_2$Erreur=="Incorrecte" & tabGraphNP14_2$LIECid!="LIECX",], 
         aes(fill=LIEC, y=Nbr, x=LIECid)) + 
    geom_bar(position='stack', stat='identity')+
    theme_classic()+
    theme(title=element_text(face="bold", size=15),
          axis.text=element_text(face="bold", size=15),
          axis.title.x = element_text(face="bold", vjust=0),#size=25
          # axis.ticks.margin=unit(c(0.5),'cm'),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_line(size=2),
          axis.ticks.length=unit(0.5, "cm"),
          axis.text.x=element_blank(),
          panel.grid=element_blank(), # Augmenter taille ligne y
          panel.grid.major.x=element_blank(), # Suppression lignes axes x
          axis.title=element_text( vjust=5, size=15),
          line=element_line(size=1.5),
          legend.position="none",
          legend.title = element_text(face="bold", size=15),
          legend.text = element_text(size=15)) +
    scale_fill_manual(values = c("#000099","#FFFF00","#009900","#FF6699","#7B7B7B"))+
    xlab(NULL)+ylab(NULL)
  
  
  ggplot(tabGraphNP14_2[tabGraphNP14_2$LIECid!="LIECX",], 
         aes(fill=LIEC, y=Nbr, x=LIECid)) + 
    geom_bar(position='stack', stat='identity')+
    theme_classic()+
    theme(title=element_text(face="bold", size=15),
          axis.text=element_text(face="bold", size=15),
          axis.title.x = element_text(face="bold", vjust=0),#size=25
          # axis.ticks.margin=unit(c(0.5),'cm'),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_line(size=2),
          axis.ticks.length=unit(0.5, "cm"),
          axis.text.x=element_blank(),
          panel.grid=element_blank(), # Augmenter taille ligne y
          panel.grid.major.x=element_blank(), # Suppression lignes axes x
          axis.title=element_text( vjust=5, size=15),
          line=element_line(size=1.5),
          legend.position="none",
          legend.title = element_text(face="bold", size=15),
          legend.text = element_text(size=15)) +
    scale_fill_manual(values = c("#000099","#FFFF00","#009900","#FF6699","#7B7B7B"))+
    xlab(NULL)+ylab(NULL)
  
  
  geom_point(aes(y=replace(Pos_Y,which(AB==0),NA), x=pos_X ),
             color="grey", size=Size_Rond, show.legend=F) + # Point pour entourer le nb de VDT
    geom_point(aes(y=replace(Pos_Y, which(AB==0),NA),  x=pos_X, fill=LIEC, shape=LIEC), 
               color="grey", size=Size_Rond+1, show.legend=F) + # Point pour entourer le nb de VDT
    geom_text(aes(y=replace(Pos_Y, which(AB==0),NA), colour=LIEC, x=pos_X), label=label_AB$AB, # Affichage valeur Ab >0
              size=Size_Eti, fontface=2, show.legend=F)    # Augmenter taille et font (gras)
  
  
  
##### BON GROUPE PC 2016 #####
  
  tabPC=Tab2[Tab2$Region%in%"PC",]
  nrow(tabPC)==4179
  dataID <- unique(sbt_m[,c("cle",colID,"Region")])
  nrow(dataID)==1109
  dataID=dataID[dataID$Region%in%"PC",]
  nrow(dataID)==189
  
  
  
  
  tabPC2016=tabPC[tabPC$Annee%in%2016,]
  nrow(tabPC2016)==511
  
  #taux d'id correcte global
  sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_EPI","EPA_EPA","ANS_ANS","END_END")])/
    sum(tabPC2016$Nbr_EW)#0.9470069  => ok
  #Tab_reu_Final2$tauxErr[Tab_reu_Final2$Region=="PC" & Tab_reu_Final2$Annee==2016]
  
  
  #les EPI identifi?s : que sont ils r?ellement ?
  {
    #nombre id correcte EPI
    n2.EPI_EPI=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_EPI")])
    
    #nombre id incorrecte EPA 
    n2.EPI_EPA=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_EPA")])
    
    #nombre id incorrecte ANS 
    n2.EPI_ANS=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_ANS")])
    
    #nombre id incorrecte END 
    n2.EPI_END=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_END")])
    
    #nombre id incorrecte LIEC_X 
    n2.EPI_LIECX=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_LIEC_X")])
    
    
    
    
    #taux id correcte EPI
    t2.EPI_EPI=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_EPI")])/
      sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
    
    #taux id incorrecte EPA 
    t2.EPI_EPA=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_EPA")])/
      sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
    
    #taux id incorrecte ANS 
    t2.EPI_ANS=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_ANS")])/
      sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
    
    #taux id incorrecte END 
    t2.EPI_END=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_END")])/
      sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
    
    #taux id incorrecte LIEC_X 
    t2.EPI_LIECX=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_LIEC_X")])/
      sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
    
    t2.EPI_EPI+
      t2.EPI_EPA+
      t2.EPI_ANS+
      t2.EPI_END+
      t2.EPI_LIECX==1
    
    
    }
  
  #les EPA identifi?s : que sont ils r?ellement ?
  {
    
    #nombre id correcte EPI
    n2.EPA_EPI=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPA_EPI")])
    
    #nombre id incorrecte EPA 
    n2.EPA_EPA=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPA_EPA")])
    
    #nombre id incorrecte ANS 
    n2.EPA_ANS=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPA_ANS")])
    
    #nombre id incorrecte END 
    n2.EPA_END=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPA_END")])
    
    #nombre id incorrecte LIEC_X 
    n2.EPA_LIECX=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPA_LIEC_X")])
    
    
    
    
    
    
    #taux id correcte EPI
    t2.EPA_EPI=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPA_EPI")])/
      sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
    
    
    
    
    #taux id incorrecte EPA 
    t2.EPA_EPA=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPA_EPA")])/
      sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
    
    
    
    
    #taux id incorrecte ANS 
    t2.EPA_ANS=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPA_ANS")])/
      sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
    
    
    
    
    
    #taux id incorrecte END 
    t2.EPA_END=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPA_END")])/
      sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
    
    
    
    
    #taux id incorrecte LIEC_X 
    t2.EPA_LIECX=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPA_LIEC_X")])/
      sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
    
    
    
    t2.EPA_EPI+
      t2.EPA_EPA+
      t2.EPA_ANS+
      t2.EPA_END+
      t2.EPA_LIECX==1
    
    
  }
  
  #les ANS identifi?s : que sont ils r?ellement ?
  {
    
    #nombre id correcte EPI
    n2.ANS_EPI=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("ANS_EPI")])
    
    #nombre id incorrecte EPA 
    n2.ANS_EPA=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("ANS_EPA")])
    
    #nombre id incorrecte ANS 
    n2.ANS_ANS=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("ANS_ANS")])
    
    #nombre id incorrecte END 
    n2.ANS_END=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("ANS_END")])
    
    #nombre id incorrecte LIEC_X 
    n2.ANS_LIECX=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("ANS_LIEC_X")])
    
    
    
    
    
    
    #taux id correcte EPI
    t2.ANS_EPI=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("ANS_EPI")])/
      sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
    
    
    
    
    #taux id incorrecte EPA 
    t2.ANS_EPA=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("ANS_EPA")])/
      sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
    
    
    
    
    #taux id incorrecte ANS 
    t2.ANS_ANS=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("ANS_ANS")])/
      sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
    
    
    
    
    
    #taux id incorrecte END 
    t2.ANS_END=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("ANS_END")])/
      sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
    
    
    
    
    #taux id incorrecte LIEC_X 
    t2.ANS_LIECX=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("ANS_LIEC_X")])/
      sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
    
    
    
    t2.ANS_EPI+
      t2.ANS_EPA+
      t2.ANS_ANS+
      t2.ANS_END+
      t2.ANS_LIECX==1
    
    
  }
  
  #les END identifi?s : que sont ils r?ellement ?
  {
    
    #nombre id correcte EPI
    n2.END_EPI=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("END_EPI")])
    
    #nombre id incorrecte EPA 
    n2.END_EPA=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("END_EPA")])
    
    #nombre id incorrecte ANS 
    n2.END_ANS=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("END_ANS")])
    
    #nombre id incorrecte END 
    n2.END_END=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("END_END")])
    
    #nombre id incorrecte LIEC_X 
    n2.END_LIECX=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("END_LIEC_X")])
    
    
    
    
    
    
    #taux id correcte EPI
    t2.END_EPI=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("END_EPI")])/
      sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
    
    
    
    
    #taux id incorrecte EPA 
    t2.END_EPA=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("END_EPA")])/
      sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
    
    
    
    
    #taux id incorrecte ANS 
    t2.END_ANS=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("END_ANS")])/
      sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
    
    
    
    
    
    #taux id incorrecte END 
    t2.END_END=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("END_END")])/
      sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
    
    
    
    
    #taux id incorrecte LIEC_X 
    t2.END_LIECX=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("END_LIEC_X")])/
      sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
    
    
    
    round(t2.END_EPI+
            t2.END_EPA+
            t2.END_ANS+
            t2.END_END+
            t2.END_LIECX,6)==1
    
    
  }
  
  #les LIEC_X identifi?s : que sont ils r?ellement ?
  {
    
    #nombre id correcte EPI
    n2.LIECX_EPI=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("LIEC_X_EPI")])
    
    #nombre id incorrecte EPA 
    n2.LIECX_EPA=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("LIEC_X_EPA")])
    
    #nombre id incorrecte ANS 
    n2.LIECX_ANS=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("LIEC_X_ANS")])
    
    #nombre id incorrecte END 
    n2.LIECX_END=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("LIEC_X_END")])
    
    #nombre id incorrecte LIEC_X 
    n2.LIECX_LIECX=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("LIEC_X_LIEC_X")])
    
    
    
    
    
    #taux id correcte EPI
    t2.LIECX_EPI=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("LIEC_X_EPI")])/
      sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
    
    
    
    
    #taux id incorrecte EPA 
    t2.LIECX_EPA=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("LIEC_X_EPA")])/
      sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
    
    
    
    
    #taux id incorrecte ANS 
    t2.LIECX_ANS=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("LIEC_X_ANS")])/
      sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
    
    
    
    
    
    #taux id incorrecte END 
    t2.LIECX_END=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("LIEC_X_END")])/
      sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
    
    
    
    
    #taux id incorrecte LIEC_X 
    t2.LIECX_LIECX=sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("LIEC_X_LIEC_X")])/
      sum(tabPC2016$Nbr_EW[tabPC2016$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
    
    
    
    t2.LIECX_EPI+
      t2.LIECX_EPA+
      t2.LIECX_ANS+
      t2.LIECX_END+
      t2.LIECX_LIECX==1
    
    
  }
  
  
  
  tabGraphPC16_2=data.frame(
    Nbr=c(
      n2.EPI_EPI,
      n2.EPI_EPA,
      n2.EPI_ANS,
      n2.EPI_END,
      n2.EPI_LIECX,
      n2.EPA_EPI,
      n2.EPA_EPA,
      n2.EPA_ANS,
      n2.EPA_END,
      n2.EPA_LIECX,
      n2.ANS_EPI,
      n2.ANS_EPA,
      n2.ANS_ANS,
      n2.ANS_END,
      n2.ANS_LIECX,
      n2.END_EPI,
      n2.END_EPA,
      n2.END_ANS,
      n2.END_END,
      n2.END_LIECX,
      n2.LIECX_EPI,
      n2.LIECX_EPA,
      n2.LIECX_ANS,
      n2.LIECX_END,
      n2.LIECX_LIECX),
    
    Taux=c(
      t2.EPI_EPI,
      t2.EPI_EPA,
      t2.EPI_ANS,
      t2.EPI_END,
      t2.EPI_LIECX,
      t2.EPA_EPI,
      t2.EPA_EPA,
      t2.EPA_ANS,
      t2.EPA_END,
      t2.EPA_LIECX,
      t2.ANS_EPI,
      t2.ANS_EPA,
      t2.ANS_ANS,
      t2.ANS_END,
      t2.ANS_LIECX,
      t2.END_EPI,
      t2.END_EPA,
      t2.END_ANS,
      t2.END_END,
      t2.END_LIECX,
      t2.LIECX_EPI,
      t2.LIECX_EPA,
      t2.LIECX_ANS,
      t2.LIECX_END,
      t2.LIECX_LIECX),
    
    Classe=c(
      "t2.EPI_EPI",
      "t2.EPI_EPA",
      "t2.EPI_ANS",
      "t2.EPI_END",
      "t2.EPI_LIECX",
      "t2.EPA_EPI",
      "t2.EPA_EPA",
      "t2.EPA_ANS",
      "t2.EPA_END",
      "t2.EPA_LIECX",
      "t2.ANS_EPI",
      "t2.ANS_EPA",
      "t2.ANS_ANS",
      "t2.ANS_END",
      "t2.ANS_LIECX",
      "t2.END_EPI",
      "t2.END_EPA",
      "t2.END_ANS",
      "t2.END_END",
      "t2.END_LIECX",
      "t2.LIECX_EPI",
      "t2.LIECX_EPA",
      "t2.LIECX_ANS",
      "t2.LIECX_END",
      "t2.LIECX_LIECX")
    
  )
  
  
  
  tabGraphPC16_2=tabGraphPC16_2[c(5,4,3,2,1,10,9,8,7,6,15,14,13,12,11,20,19,18,17,16,25,24,23,22,21),]
  tabGraphPC16_2$LIEC=factor(substr(tabGraphPC16_2$Classe,8,10),levels=c("EPI","EPA","ANS","END","LIECX"))
  tabGraphPC16_2$LIECid=factor(substr(tabGraphPC16_2$Classe,4,6),levels=c("EPI","EPA","ANS","END","LIECX"))
  #tabGraphPC16_2$Nbr=round(tabGraphPC16_2$Nbr/4)#On cheat le nombre pour que ce soit du mm ordre que AL 12
  
  anti_join(tabGraphPC16_2,ddply(tabGraphPC16_2, .(Classe), mutate, Pos_Y=cumsum(Nbr)-Nbr/2 ))
  anti_join(ddply(tabGraphPC16_2, .(Classe), mutate, Pos_Y=cumsum(Nbr)-Nbr/2 ),tabGraphPC16_2)
  tabGraphPC16_2 <- ddply(tabGraphPC16_2, .(LIECid), mutate, Pos_Y=cumsum(Nbr)-Nbr/2 )
  
  tabGraphPC16_2=tabGraphPC16_2[tabGraphPC16_2$Nbr>0,]
  tabGraphPC16_2=tabGraphPC16_2[tabGraphPC16_2$LIEC!="LIECX",]
  
  vpos_X=c()
  for(iLIEC in levels(tabGraphPC16_2$LIECid)){
    intmd=seq(0.4,-0.4,length.out = nrow(tabGraphPC16_2[tabGraphPC16_2$LIECid%in%iLIEC,]))
    vpos_X=c(vpos_X,intmd+which(levels(tabGraphPC16_2$LIECid)%in%iLIEC))
  }
  tabGraphPC16_2$pos_X <- vpos_X
  
  tabGraphPC16_2$couleur=rep(NA,nrow(tabGraphPC16_2))
  tabGraphPC16_2$couleur[tabGraphPC16_2$LIEC%in%"EPI"]="#000099"
  tabGraphPC16_2$couleur[tabGraphPC16_2$LIEC%in%"EPA"]="#FFFF00"
  tabGraphPC16_2$couleur[tabGraphPC16_2$LIEC%in%"ANS"]="#009900"
  tabGraphPC16_2$couleur[tabGraphPC16_2$LIEC%in%"END"]="#FF6699"
  tabGraphPC16_2$couleur[tabGraphPC16_2$LIEC%in%"LIECX"]="#7B7B7B"
  
  
  
  tabGraphPC16_2$colecr=rep(NA,nrow(tabGraphPC16_2))
  tabGraphPC16_2$colecr[tabGraphPC16_2$LIEC%in%"EPI"]="#FFFFFF"
  tabGraphPC16_2$colecr[tabGraphPC16_2$LIEC%in%"EPA"]="#000000"
  tabGraphPC16_2$colecr[tabGraphPC16_2$LIEC%in%"ANS"]="#FFFFFF"
  tabGraphPC16_2$colecr[tabGraphPC16_2$LIEC%in%"END"]="#000000"
  tabGraphPC16_2$colecr[tabGraphPC16_2$LIEC%in%"LIECX"]="#FFFFFF"
  
  

  
  
  
  ggplot(tabGraphPC16_2, aes(fill=LIEC, y=Taux, x=LIECid)) + 
    geom_bar(position='stack', stat='identity')+
    theme_classic()+
    scale_fill_manual(values = c("#000099","#FFFF00","#009900","#FF6699","#7B7B7B"))
  
  ggplot(tabGraphPC16_2, aes(fill=LIEC, y=Taux, x=LIECid)) + 
    geom_bar(position='fill', stat='identity')+
    theme_classic()+
    scale_fill_manual(values = c("#000099","#FFFF00","#009900","#FF6699","#7B7B7B"))
  
  ggplot(tabGraphPC16_2, aes(fill=LIEC, y=Nbr, x=LIECid)) + 
    geom_bar(position='fill', stat='identity')+
    theme_classic()+
    scale_fill_manual(values = c("#000099","#FFFF00","#009900","#FF6699","#7B7B7B"))
  
  ggplot(tabGraphPC16_2, aes(fill=LIEC, y=Nbr, x=LIECid)) + 
    geom_bar(position='stack', stat='identity')+
    theme_classic()+
    scale_fill_manual(values = c("#000099","#FFFF00","#009900","#FF6699","#7B7B7B"))
  
  
  
  
  
  
  
  #zoom sur les erreurs
  tabGraphPC16_2$Erreur=NA
  tabGraphPC16_2$Erreur[tabGraphPC16_2$LIEC==tabGraphPC16_2$LIECid]="Correcte"
  tabGraphPC16_2$Erreur[tabGraphPC16_2$LIEC!=tabGraphPC16_2$LIECid]="Incorrecte"
  
  
  
  ggplot(tabGraphPC16_2[tabGraphPC16_2$Erreur=="Incorrecte" & tabGraphPC16_2$LIECid!="LIECX",], 
         aes(fill=LIEC, y=Taux, x=LIECid)) + 
    geom_bar(position='stack', stat='identity')+
    theme_classic()+
    theme(title=element_text(face="bold", size=15),
          axis.text=element_text(face="bold", size=15),
          axis.title.x = element_text(face="bold", vjust=0),#size=25
          # axis.ticks.margin=unit(c(0.5),'cm'),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_line(size=2),
          axis.ticks.length=unit(0.5, "cm"),
          axis.text.x=element_blank(),
          panel.grid=element_blank(), # Augmenter taille ligne y
          panel.grid.major.x=element_blank(), # Suppression lignes axes x
          axis.title=element_text( vjust=5, size=15),
          line=element_line(size=1.5),
          legend.position="none",
          legend.title = element_text(face="bold", size=15),
          legend.text = element_text(size=15)) +
    scale_y_continuous(labels = scales::percent, limits=c(0,1)) + 
    scale_fill_manual(values = c("#000099","#FFFF00","#009900","#FF6699","#7B7B7B"))+
    xlab(NULL)+ylab(NULL)
  
  
  
  ggplot(tabGraphPC16_2[tabGraphPC16_2$Erreur=="Incorrecte" & tabGraphPC16_2$LIECid!="LIECX",], 
         aes(fill=LIEC, y=Taux, x=LIECid)) + 
    geom_bar(position='stack', stat='identity')+
    theme_classic()+
    theme(title=element_text(face="bold", size=15),
          axis.text=element_text(face="bold", size=15),
          axis.title.x = element_text(face="bold", vjust=0),#size=25
          # axis.ticks.margin=unit(c(0.5),'cm'),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_line(size=2),
          axis.ticks.length=unit(0.5, "cm"),
          axis.text.x=element_blank(),
          panel.grid=element_blank(), # Augmenter taille ligne y
          panel.grid.major.x=element_blank(), # Suppression lignes axes x
          axis.title=element_text( vjust=5, size=15),
          line=element_line(size=1.5),
          legend.position="none",
          legend.title = element_text(face="bold", size=15),
          legend.text = element_text(size=15)) +
    scale_y_continuous(labels = scales::percent, limits=c(0,0.25)) + 
    scale_fill_manual(values = c("#000099","#FFFF00","#009900","#FF6699","#7B7B7B"))+
    xlab(NULL)+ylab(NULL)
  
  
  
  ggplot(tabGraphPC16_2[tabGraphPC16_2$Erreur=="Incorrecte" & tabGraphPC16_2$LIECid!="LIECX",], 
         aes(fill=LIEC, y=Nbr, x=LIECid)) + 
    geom_bar(position='stack', stat='identity')+
    theme_classic()+
    theme(title=element_text(face="bold", size=15),
          axis.text=element_text(face="bold", size=15),
          axis.title.x = element_text(face="bold", vjust=0),#size=25
          # axis.ticks.margin=unit(c(0.5),'cm'),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_line(size=2),
          axis.ticks.length=unit(0.5, "cm"),
          axis.text.x=element_blank(),
          panel.grid=element_blank(), # Augmenter taille ligne y
          panel.grid.major.x=element_blank(), # Suppression lignes axes x
          axis.title=element_text( vjust=5, size=15),
          line=element_line(size=1.5),
          legend.position="none",
          legend.title = element_text(face="bold", size=15),
          legend.text = element_text(size=15)) +
    scale_fill_manual(values = c("#000099","#FFFF00","#009900","#FF6699","#7B7B7B"))+
    xlab(NULL)+ylab(NULL)
  
  
  ggplot(tabGraphPC16_2[tabGraphPC16_2$LIECid!="LIECX",], 
         aes(fill=LIEC, y=Nbr, x=LIECid)) + 
    geom_bar(position='stack', stat='identity')+
    theme_classic()+
    theme(title=element_text(face="bold", size=15),
          axis.text=element_text(face="bold", size=15),
          axis.title.x = element_text(face="bold", vjust=0),#size=25
          # axis.ticks.margin=unit(c(0.5),'cm'),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_line(size=2),
          axis.ticks.length=unit(0.5, "cm"),
          axis.text.x=element_blank(),
          panel.grid=element_blank(), # Augmenter taille ligne y
          panel.grid.major.x=element_blank(), # Suppression lignes axes x
          axis.title=element_text( vjust=5, size=15),
          line=element_line(size=1.5),
          legend.position="none",
          legend.title = element_text(face="bold", size=15),
          legend.text = element_text(size=15)) +
    scale_fill_manual(values = c("#000099","#FFFF00","#009900","#FF6699","#7B7B7B"))+
    xlab(NULL)+ylab(NULL)+ylim(0,1000)+
  
  
  geom_point(aes(y=Pos_Y, x=pos_X ),color="grey", size=14, show.legend=F)+  # Point pour entourer le nb de VDT
  geom_point(aes(y=Pos_Y, x=pos_X ),color=tabGraphPC16_2$couleur[tabGraphPC16_2$LIECid!="LIECX"], size=12, show.legend=F)+  # Point pour entourer le nb de VDT
  geom_text(aes(y=Pos_Y, x=pos_X ), label=tabGraphPC16_2$Nbr[tabGraphPC16_2$LIECid!="LIECX"], # Affichage valeur Ab >0
            size=5, fontface=2, color=tabGraphPC16_2$colecr[tabGraphPC16_2$LIECid!="LIECX"], show.legend=F)    # Augmenter taille et font (gras)

  
  
##### GROUPE nat #####
  
  nrow(Tab2)==21044
  dataID <- unique(sbt_m[,c("cle",colID,"Region")])
  nrow(dataID)==1109
 
  
  Tabnat=Tab2[!(Tab2$Region%in%"CO" & Tab2$Annee%in%2014),]
  Tabnat=Tab2[!(Tab2$Region%in%"CO" & Tab2$Annee%in%2016),]
  Tabnat=Tab2[!(Tab2$Region%in%"AU" & Tab2$Annee%in%2014),]
  Tabnat=Tab2[!(Tab2$Region%in%"AU" & Tab2$Annee%in%2017),]
  Tabnat=Tab2[!(Tab2$Region%in%"AU" & Tab2$Annee%in%2018),]

  Tabnat=Tab2[!(Tab2$Region%in%"PA" & Tab2$Annee%in%2015),]
  Tabnat=Tab2[!(Tab2$Region%in%"PA" & Tab2$Annee%in%2016),]
  Tabnat=Tab2[!(Tab2$Region%in%"PA" & Tab2$Annee%in%2017),]
  Tabnat=Tab2[!(Tab2$Region%in%"PA" & Tab2$Annee%in%2018),]
  Tabnat=Tab2[!(Tab2$Region%in%"AU" & Tab2$Annee%in%2015),]
  
  
  
  nrow(Tabnat)==21042
  
  #taux d'id correcte global
  sum(Tabnat$Nbr_EW[Tabnat$FIEC_LIEC%in%c("EPI_EPI","EPA_EPA","ANS_ANS","END_END")])/
    sum(Tabnat$Nbr_EW)#0.7204936  

  
  #les EPI identifi?s : que sont ils r?ellement ?
  {
    #nombre id correcte EPI
    n2.EPI_EPI=sum(Tabnat$Nbr_EW[Tabnat$FIEC_LIEC%in%c("EPI_EPI")])
    
    #nombre id incorrecte EPA 
    n2.EPI_EPA=sum(Tabnat$Nbr_EW[Tabnat$FIEC_LIEC%in%c("EPI_EPA")])
    
    #nombre id incorrecte ANS 
    n2.EPI_ANS=sum(Tabnat$Nbr_EW[Tabnat$FIEC_LIEC%in%c("EPI_ANS")])
    
    #nombre id incorrecte END 
    n2.EPI_END=sum(Tabnat$Nbr_EW[Tabnat$FIEC_LIEC%in%c("EPI_END")])
    
    #nombre id incorrecte LIEC_X 
    n2.EPI_LIECX=sum(Tabnat$Nbr_EW[Tabnat$FIEC_LIEC%in%c("EPI_LIEC_X")])
    
    
    
    
    #taux id correcte EPI
    t2.EPI_EPI=sum(Tabnat$Nbr_EW[Tabnat$FIEC_LIEC%in%c("EPI_EPI")])/
      sum(Tabnat$Nbr_EW[Tabnat$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
    
    #taux id incorrecte EPA 
    t2.EPI_EPA=sum(Tabnat$Nbr_EW[Tabnat$FIEC_LIEC%in%c("EPI_EPA")])/
      sum(Tabnat$Nbr_EW[Tabnat$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
    
    #taux id incorrecte ANS 
    t2.EPI_ANS=sum(Tabnat$Nbr_EW[Tabnat$FIEC_LIEC%in%c("EPI_ANS")])/
      sum(Tabnat$Nbr_EW[Tabnat$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
    
    #taux id incorrecte END 
    t2.EPI_END=sum(Tabnat$Nbr_EW[Tabnat$FIEC_LIEC%in%c("EPI_END")])/
      sum(Tabnat$Nbr_EW[Tabnat$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
    
    #taux id incorrecte LIEC_X 
    t2.EPI_LIECX=sum(Tabnat$Nbr_EW[Tabnat$FIEC_LIEC%in%c("EPI_LIEC_X")])/
      sum(Tabnat$Nbr_EW[Tabnat$FIEC_LIEC%in%c("EPI_EPI","EPI_EPA","EPI_ANS","EPI_END","EPI_LIEC_X")])
    
    t2.EPI_EPI+
      t2.EPI_EPA+
      t2.EPI_ANS+
      t2.EPI_END+
      t2.EPI_LIECX==1
    
    
    }
  
  #les EPA identifi?s : que sont ils r?ellement ?
  {
    
    #nombre id correcte EPI
    n2.EPA_EPI=sum(Tabnat$Nbr_EW[Tabnat$FIEC_LIEC%in%c("EPA_EPI")])
    
    #nombre id incorrecte EPA 
    n2.EPA_EPA=sum(Tabnat$Nbr_EW[Tabnat$FIEC_LIEC%in%c("EPA_EPA")])
    
    #nombre id incorrecte ANS 
    n2.EPA_ANS=sum(Tabnat$Nbr_EW[Tabnat$FIEC_LIEC%in%c("EPA_ANS")])
    
    #nombre id incorrecte END 
    n2.EPA_END=sum(Tabnat$Nbr_EW[Tabnat$FIEC_LIEC%in%c("EPA_END")])
    
    #nombre id incorrecte LIEC_X 
    n2.EPA_LIECX=sum(Tabnat$Nbr_EW[Tabnat$FIEC_LIEC%in%c("EPA_LIEC_X")])
    
    
    
    
    
    
    #taux id correcte EPI
    t2.EPA_EPI=sum(Tabnat$Nbr_EW[Tabnat$FIEC_LIEC%in%c("EPA_EPI")])/
      sum(Tabnat$Nbr_EW[Tabnat$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
    
    
    
    
    #taux id incorrecte EPA 
    t2.EPA_EPA=sum(Tabnat$Nbr_EW[Tabnat$FIEC_LIEC%in%c("EPA_EPA")])/
      sum(Tabnat$Nbr_EW[Tabnat$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
    
    
    
    
    #taux id incorrecte ANS 
    t2.EPA_ANS=sum(Tabnat$Nbr_EW[Tabnat$FIEC_LIEC%in%c("EPA_ANS")])/
      sum(Tabnat$Nbr_EW[Tabnat$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
    
    
    
    
    
    #taux id incorrecte END 
    t2.EPA_END=sum(Tabnat$Nbr_EW[Tabnat$FIEC_LIEC%in%c("EPA_END")])/
      sum(Tabnat$Nbr_EW[Tabnat$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
    
    
    
    
    #taux id incorrecte LIEC_X 
    t2.EPA_LIECX=sum(Tabnat$Nbr_EW[Tabnat$FIEC_LIEC%in%c("EPA_LIEC_X")])/
      sum(Tabnat$Nbr_EW[Tabnat$FIEC_LIEC%in%c("EPA_EPI","EPA_EPA","EPA_ANS","EPA_END","EPA_LIEC_X")])
    
    
    
    round(t2.EPA_EPI+
      t2.EPA_EPA+
      t2.EPA_ANS+
      t2.EPA_END+
      t2.EPA_LIECX,6)==1
    
    
  }
  
  #les ANS identifi?s : que sont ils r?ellement ?
  {
    
    #nombre id correcte EPI
    n2.ANS_EPI=sum(Tabnat$Nbr_EW[Tabnat$FIEC_LIEC%in%c("ANS_EPI")])
    
    #nombre id incorrecte EPA 
    n2.ANS_EPA=sum(Tabnat$Nbr_EW[Tabnat$FIEC_LIEC%in%c("ANS_EPA")])
    
    #nombre id incorrecte ANS 
    n2.ANS_ANS=sum(Tabnat$Nbr_EW[Tabnat$FIEC_LIEC%in%c("ANS_ANS")])
    
    #nombre id incorrecte END 
    n2.ANS_END=sum(Tabnat$Nbr_EW[Tabnat$FIEC_LIEC%in%c("ANS_END")])
    
    #nombre id incorrecte LIEC_X 
    n2.ANS_LIECX=sum(Tabnat$Nbr_EW[Tabnat$FIEC_LIEC%in%c("ANS_LIEC_X")])
    
    
    
    
    
    
    #taux id correcte EPI
    t2.ANS_EPI=sum(Tabnat$Nbr_EW[Tabnat$FIEC_LIEC%in%c("ANS_EPI")])/
      sum(Tabnat$Nbr_EW[Tabnat$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
    
    
    
    
    #taux id incorrecte EPA 
    t2.ANS_EPA=sum(Tabnat$Nbr_EW[Tabnat$FIEC_LIEC%in%c("ANS_EPA")])/
      sum(Tabnat$Nbr_EW[Tabnat$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
    
    
    
    
    #taux id incorrecte ANS 
    t2.ANS_ANS=sum(Tabnat$Nbr_EW[Tabnat$FIEC_LIEC%in%c("ANS_ANS")])/
      sum(Tabnat$Nbr_EW[Tabnat$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
    
    
    
    
    
    #taux id incorrecte END 
    t2.ANS_END=sum(Tabnat$Nbr_EW[Tabnat$FIEC_LIEC%in%c("ANS_END")])/
      sum(Tabnat$Nbr_EW[Tabnat$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
    
    
    
    
    #taux id incorrecte LIEC_X 
    t2.ANS_LIECX=sum(Tabnat$Nbr_EW[Tabnat$FIEC_LIEC%in%c("ANS_LIEC_X")])/
      sum(Tabnat$Nbr_EW[Tabnat$FIEC_LIEC%in%c("ANS_EPI","ANS_EPA","ANS_ANS","ANS_END","ANS_LIEC_X")])
    
    
    
    t2.ANS_EPI+
      t2.ANS_EPA+
      t2.ANS_ANS+
      t2.ANS_END+
      t2.ANS_LIECX==1
    
    
  }
  
  #les END identifi?s : que sont ils r?ellement ?
  {
    
    #nombre id correcte EPI
    n2.END_EPI=sum(Tabnat$Nbr_EW[Tabnat$FIEC_LIEC%in%c("END_EPI")])
    
    #nombre id incorrecte EPA 
    n2.END_EPA=sum(Tabnat$Nbr_EW[Tabnat$FIEC_LIEC%in%c("END_EPA")])
    
    #nombre id incorrecte ANS 
    n2.END_ANS=sum(Tabnat$Nbr_EW[Tabnat$FIEC_LIEC%in%c("END_ANS")])
    
    #nombre id incorrecte END 
    n2.END_END=sum(Tabnat$Nbr_EW[Tabnat$FIEC_LIEC%in%c("END_END")])
    
    #nombre id incorrecte LIEC_X 
    n2.END_LIECX=sum(Tabnat$Nbr_EW[Tabnat$FIEC_LIEC%in%c("END_LIEC_X")])
    
    
    
    
    
    
    #taux id correcte EPI
    t2.END_EPI=sum(Tabnat$Nbr_EW[Tabnat$FIEC_LIEC%in%c("END_EPI")])/
      sum(Tabnat$Nbr_EW[Tabnat$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
    
    
    
    
    #taux id incorrecte EPA 
    t2.END_EPA=sum(Tabnat$Nbr_EW[Tabnat$FIEC_LIEC%in%c("END_EPA")])/
      sum(Tabnat$Nbr_EW[Tabnat$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
    
    
    
    
    #taux id incorrecte ANS 
    t2.END_ANS=sum(Tabnat$Nbr_EW[Tabnat$FIEC_LIEC%in%c("END_ANS")])/
      sum(Tabnat$Nbr_EW[Tabnat$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
    
    
    
    
    
    #taux id incorrecte END 
    t2.END_END=sum(Tabnat$Nbr_EW[Tabnat$FIEC_LIEC%in%c("END_END")])/
      sum(Tabnat$Nbr_EW[Tabnat$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
    
    
    
    
    #taux id incorrecte LIEC_X 
    t2.END_LIECX=sum(Tabnat$Nbr_EW[Tabnat$FIEC_LIEC%in%c("END_LIEC_X")])/
      sum(Tabnat$Nbr_EW[Tabnat$FIEC_LIEC%in%c("END_EPI","END_EPA","END_ANS","END_END","END_LIEC_X")])
    
    
    
    round(t2.END_EPI+
            t2.END_EPA+
            t2.END_ANS+
            t2.END_END+
            t2.END_LIECX,6)==1
    
    
  }
  
  #les LIEC_X identifi?s : que sont ils r?ellement ?
  {
    
    #nombre id correcte EPI
    n2.LIECX_EPI=sum(Tabnat$Nbr_EW[Tabnat$FIEC_LIEC%in%c("LIEC_X_EPI")])
    
    #nombre id incorrecte EPA 
    n2.LIECX_EPA=sum(Tabnat$Nbr_EW[Tabnat$FIEC_LIEC%in%c("LIEC_X_EPA")])
    
    #nombre id incorrecte ANS 
    n2.LIECX_ANS=sum(Tabnat$Nbr_EW[Tabnat$FIEC_LIEC%in%c("LIEC_X_ANS")])
    
    #nombre id incorrecte END 
    n2.LIECX_END=sum(Tabnat$Nbr_EW[Tabnat$FIEC_LIEC%in%c("LIEC_X_END")])
    
    #nombre id incorrecte LIEC_X 
    n2.LIECX_LIECX=sum(Tabnat$Nbr_EW[Tabnat$FIEC_LIEC%in%c("LIEC_X_LIEC_X")])
    
    
    
    
    
    #taux id correcte EPI
    t2.LIECX_EPI=sum(Tabnat$Nbr_EW[Tabnat$FIEC_LIEC%in%c("LIEC_X_EPI")])/
      sum(Tabnat$Nbr_EW[Tabnat$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
    
    
    
    
    #taux id incorrecte EPA 
    t2.LIECX_EPA=sum(Tabnat$Nbr_EW[Tabnat$FIEC_LIEC%in%c("LIEC_X_EPA")])/
      sum(Tabnat$Nbr_EW[Tabnat$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
    
    
    
    
    #taux id incorrecte ANS 
    t2.LIECX_ANS=sum(Tabnat$Nbr_EW[Tabnat$FIEC_LIEC%in%c("LIEC_X_ANS")])/
      sum(Tabnat$Nbr_EW[Tabnat$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
    
    
    
    
    
    #taux id incorrecte END 
    t2.LIECX_END=sum(Tabnat$Nbr_EW[Tabnat$FIEC_LIEC%in%c("LIEC_X_END")])/
      sum(Tabnat$Nbr_EW[Tabnat$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
    
    
    
    
    #taux id incorrecte LIEC_X 
    t2.LIECX_LIECX=sum(Tabnat$Nbr_EW[Tabnat$FIEC_LIEC%in%c("LIEC_X_LIEC_X")])/
      sum(Tabnat$Nbr_EW[Tabnat$FIEC_LIEC%in%c("LIEC_X_EPI","LIEC_X_EPA","LIEC_X_ANS","LIEC_X_END","LIEC_X_LIEC_X")])
    
    
    
    t2.LIECX_EPI+
      t2.LIECX_EPA+
      t2.LIECX_ANS+
      t2.LIECX_END+
      t2.LIECX_LIECX==1
    
    
  }
  
  
  
  tabGraphNAT_2=data.frame(
    Nbr=c(
      n2.EPI_EPI,
      n2.EPI_EPA,
      n2.EPI_ANS,
      n2.EPI_END,
      n2.EPI_LIECX,
      n2.EPA_EPI,
      n2.EPA_EPA,
      n2.EPA_ANS,
      n2.EPA_END,
      n2.EPA_LIECX,
      n2.ANS_EPI,
      n2.ANS_EPA,
      n2.ANS_ANS,
      n2.ANS_END,
      n2.ANS_LIECX,
      n2.END_EPI,
      n2.END_EPA,
      n2.END_ANS,
      n2.END_END,
      n2.END_LIECX,
      n2.LIECX_EPI,
      n2.LIECX_EPA,
      n2.LIECX_ANS,
      n2.LIECX_END,
      n2.LIECX_LIECX),
    
    Taux=c(
      t2.EPI_EPI,
      t2.EPI_EPA,
      t2.EPI_ANS,
      t2.EPI_END,
      t2.EPI_LIECX,
      t2.EPA_EPI,
      t2.EPA_EPA,
      t2.EPA_ANS,
      t2.EPA_END,
      t2.EPA_LIECX,
      t2.ANS_EPI,
      t2.ANS_EPA,
      t2.ANS_ANS,
      t2.ANS_END,
      t2.ANS_LIECX,
      t2.END_EPI,
      t2.END_EPA,
      t2.END_ANS,
      t2.END_END,
      t2.END_LIECX,
      t2.LIECX_EPI,
      t2.LIECX_EPA,
      t2.LIECX_ANS,
      t2.LIECX_END,
      t2.LIECX_LIECX),
    
    Classe=c(
      "t2.EPI_EPI",
      "t2.EPI_EPA",
      "t2.EPI_ANS",
      "t2.EPI_END",
      "t2.EPI_LIECX",
      "t2.EPA_EPI",
      "t2.EPA_EPA",
      "t2.EPA_ANS",
      "t2.EPA_END",
      "t2.EPA_LIECX",
      "t2.ANS_EPI",
      "t2.ANS_EPA",
      "t2.ANS_ANS",
      "t2.ANS_END",
      "t2.ANS_LIECX",
      "t2.END_EPI",
      "t2.END_EPA",
      "t2.END_ANS",
      "t2.END_END",
      "t2.END_LIECX",
      "t2.LIECX_EPI",
      "t2.LIECX_EPA",
      "t2.LIECX_ANS",
      "t2.LIECX_END",
      "t2.LIECX_LIECX")
    
  )
  
  
  tabGraphNAT_2=tabGraphNAT_2[c(5,4,3,2,1,10,9,8,7,6,15,14,13,12,11,20,19,18,17,16,25,24,23,22,21),]
  tabGraphNAT_2$LIEC=factor(substr(tabGraphNAT_2$Classe,8,10),levels=c("EPI","EPA","ANS","END","LIECX"))
  tabGraphNAT_2$LIECid=factor(substr(tabGraphNAT_2$Classe,4,6),levels=c("EPI","EPA","ANS","END","LIECX"))
  tabGraphNAT_2=tabGraphNAT_2[tabGraphNAT_2$Nbr>0,]
  tabGraphNAT_2=tabGraphNAT_2[tabGraphNAT_2$LIEC!="LIECX",]

  anti_join(tabGraphNAT_2,ddply(tabGraphNAT_2, .(Classe), mutate, Pos_Y=cumsum(Nbr)-Nbr/2 ))
  anti_join(ddply(tabGraphNAT_2, .(Classe), mutate, Pos_Y=cumsum(Nbr)-Nbr/2 ),tabGraphNAT_2)
  tabGraphNAT_2 <- ddply(tabGraphNAT_2, .(LIECid), mutate, Pos_Y=cumsum(Nbr)-Nbr/2 )
  
  
  vpos_X=c()
  for(iLIEC in levels(tabGraphNAT_2$LIECid)){
    intmd=seq(0.4,-0.4,length.out = nrow(tabGraphNAT_2[tabGraphNAT_2$LIECid%in%iLIEC,]))
    vpos_X=c(vpos_X,intmd+which(levels(tabGraphNAT_2$LIECid)%in%iLIEC))
  }
  tabGraphNAT_2$pos_X <- vpos_X
  
  tabGraphNAT_2$couleur=rep(NA,nrow(tabGraphNAT_2))
  tabGraphNAT_2$couleur[tabGraphNAT_2$LIEC%in%"EPI"]="#000099"
  tabGraphNAT_2$couleur[tabGraphNAT_2$LIEC%in%"EPA"]="#FFFF00"
  tabGraphNAT_2$couleur[tabGraphNAT_2$LIEC%in%"ANS"]="#009900"
  tabGraphNAT_2$couleur[tabGraphNAT_2$LIEC%in%"END"]="#FF6699"
  tabGraphNAT_2$couleur[tabGraphNAT_2$LIEC%in%"LIECX"]="#7B7B7B"
  
  
  
  tabGraphNAT_2$colecr=rep(NA,nrow(tabGraphNAT_2))
  tabGraphNAT_2$colecr[tabGraphNAT_2$LIEC%in%"EPI"]="#FFFFFF"
  tabGraphNAT_2$colecr[tabGraphNAT_2$LIEC%in%"EPA"]="#000000"
  tabGraphNAT_2$colecr[tabGraphNAT_2$LIEC%in%"ANS"]="#FFFFFF"
  tabGraphNAT_2$colecr[tabGraphNAT_2$LIEC%in%"END"]="#000000"
  tabGraphNAT_2$colecr[tabGraphNAT_2$LIEC%in%"LIECX"]="#FFFFFF"
  
  
  
  
  
  
  ggplot(tabGraphNAT_2, aes(fill=LIEC, y=Taux, x=LIECid)) + 
    geom_bar(position='stack', stat='identity')+
    theme_classic()+
    scale_fill_manual(values = c("#000099","#FFFF00","#009900","#FF6699","#7B7B7B"))
  
  ggplot(tabGraphNAT_2, aes(fill=LIEC, y=Taux, x=LIECid)) + 
    geom_bar(position='fill', stat='identity')+
    theme_classic()+
    scale_fill_manual(values = c("#000099","#FFFF00","#009900","#FF6699","#7B7B7B"))
  
  ggplot(tabGraphNAT_2, aes(fill=LIEC, y=Nbr, x=LIECid)) + 
    geom_bar(position='fill', stat='identity')+
    theme_classic()+
    scale_fill_manual(values = c("#000099","#FFFF00","#009900","#FF6699","#7B7B7B"))
  
  ggplot(tabGraphNAT_2, aes(fill=LIEC, y=Nbr, x=LIECid)) + 
    geom_bar(position='stack', stat='identity')+
    theme_classic()+
    scale_fill_manual(values = c("#000099","#FFFF00","#009900","#FF6699","#7B7B7B"))
  
  
  
  
  
  
  
  #zoom sur les erreurs
  tabGraphNAT_2$Erreur=NA
  tabGraphNAT_2$Erreur[tabGraphNAT_2$LIEC==tabGraphNAT_2$LIECid]="Correcte"
  tabGraphNAT_2$Erreur[tabGraphNAT_2$LIEC!=tabGraphNAT_2$LIECid]="Incorrecte"
  
  
  
  ggplot(tabGraphNAT_2[tabGraphNAT_2$Erreur=="Incorrecte" & tabGraphNAT_2$LIECid!="LIECX",], 
         aes(fill=LIEC, y=Taux, x=LIECid)) + 
    geom_bar(position='stack', stat='identity')+
    theme_classic()+
    theme(title=element_text(face="bold", size=15),
          axis.text=element_text(face="bold", size=15),
          axis.title.x = element_text(face="bold", vjust=0),#size=25
          # axis.ticks.margin=unit(c(0.5),'cm'),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_line(size=2),
          axis.ticks.length=unit(0.5, "cm"),
          axis.text.x=element_blank(),
          panel.grid=element_blank(), # Augmenter taille ligne y
          panel.grid.major.x=element_blank(), # Suppression lignes axes x
          axis.title=element_text( vjust=5, size=15),
          line=element_line(size=1.5),
          legend.position="none",
          legend.title = element_text(face="bold", size=15),
          legend.text = element_text(size=15)) +
    scale_y_continuous(labels = scales::percent, limits=c(0,1)) + 
    scale_fill_manual(values = c("#000099","#FFFF00","#009900","#FF6699","#7B7B7B"))+
    xlab(NULL)+ylab(NULL)
  
  
  
  ggplot(tabGraphNAT_2[tabGraphNAT_2$Erreur=="Incorrecte" & tabGraphNAT_2$LIECid!="LIECX",], 
         aes(fill=LIEC, y=Taux, x=LIECid)) + 
    geom_bar(position='stack', stat='identity')+
    theme_classic()+
    theme(title=element_text(face="bold", size=15),
          axis.text=element_text(face="bold", size=15),
          axis.title.x = element_text(face="bold", vjust=0),#size=25
          # axis.ticks.margin=unit(c(0.5),'cm'),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_line(size=2),
          axis.ticks.length=unit(0.5, "cm"),
          axis.text.x=element_blank(),
          panel.grid=element_blank(), # Augmenter taille ligne y
          panel.grid.major.x=element_blank(), # Suppression lignes axes x
          axis.title=element_text( vjust=5, size=15),
          line=element_line(size=1.5),
          legend.position="none",
          legend.title = element_text(face="bold", size=15),
          legend.text = element_text(size=15)) +
    scale_y_continuous(labels = scales::percent, limits=c(0,0.25)) + 
    scale_fill_manual(values = c("#000099","#FFFF00","#009900","#FF6699","#7B7B7B"))+
    xlab(NULL)+ylab(NULL)
  
  
  
  ggplot(tabGraphNAT_2[tabGraphNAT_2$Erreur=="Incorrecte" & tabGraphNAT_2$LIECid!="LIECX",], 
         aes(fill=LIEC, y=Nbr, x=LIECid)) + 
    geom_bar(position='stack', stat='identity')+
    theme_classic()+
    theme(title=element_text(face="bold", size=15),
          axis.text=element_text(face="bold", size=15),
          axis.title.x = element_text(face="bold", vjust=0),#size=25
          # axis.ticks.margin=unit(c(0.5),'cm'),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_line(size=2),
          axis.ticks.length=unit(0.5, "cm"),
          axis.text.x=element_blank(),
          panel.grid=element_blank(), # Augmenter taille ligne y
          panel.grid.major.x=element_blank(), # Suppression lignes axes x
          axis.title=element_text( vjust=5, size=15),
          line=element_line(size=1.5),
          legend.position="none",
          legend.title = element_text(face="bold", size=15),
          legend.text = element_text(size=15)) +
    scale_fill_manual(values = c("#000099","#FFFF00","#009900","#FF6699","#7B7B7B"))+
    xlab(NULL)+ylab(NULL)
  
  
  ggplot(tabGraphNAT_2[tabGraphNAT_2$LIECid!="LIECX",], 
         aes(fill=LIEC, y=Nbr, x=LIECid)) + 
    geom_bar(position='stack', stat='identity')+
    theme_classic()+
    theme(title=element_text(face="bold", size=15),
          axis.text=element_text(face="bold", size=15),
          axis.title.x = element_text(face="bold", vjust=0),#size=25
          # axis.ticks.margin=unit(c(0.5),'cm'),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_line(size=2),
          axis.ticks.length=unit(0.5, "cm"),
          axis.text.x=element_blank(),
          panel.grid=element_blank(), # Augmenter taille ligne y
          panel.grid.major.x=element_blank(), # Suppression lignes axes x
          axis.title=element_text( vjust=5, size=15),
          line=element_line(size=1.5),
          legend.position="none",
          legend.title = element_text(face="bold", size=15),
          legend.text = element_text(size=15)) +
    scale_fill_manual(values = c("#000099","#FFFF00","#009900","#FF6699","#7B7B7B"))+
    xlab(NULL)+ylab(NULL)+
    geom_point(aes(y=Pos_Y, x=pos_X ),color="grey", size=17, show.legend=F)+  # Point pour entourer le nb de VDT
    geom_point(aes(y=Pos_Y, x=pos_X ),color=tabGraphNAT_2$couleur[tabGraphNAT_2$LIECid!="LIECX"], size=15, show.legend=F)+  # Point pour entourer le nb de VDT
    geom_text(aes(y=Pos_Y, x=pos_X ), label=tabGraphNAT_2$Nbr[tabGraphNAT_2$LIECid!="LIECX"], # Affichage valeur Ab >0
              size=5, fontface=2, color=tabGraphNAT_2$colecr[tabGraphNAT_2$LIECid!="LIECX"], show.legend=F)    # Augmenter taille et font (gras)
  
  
    
  
  
#####
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  

#load("C:/Users/Nathan/Documents/ECOBIO/Equipe/Abstract ISEE/envTauxErr_v2022.06.13.Rdata")




































