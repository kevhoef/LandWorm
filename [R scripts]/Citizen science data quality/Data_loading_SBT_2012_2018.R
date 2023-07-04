#### PACKAGE LOADING ####
library(readxl)
library(dplyr)

#### DATASET LOADING ####
sbt_m <- read_excel("./[Database]/Raw datasets/National/SBT_ENI/ENI _ Data Compilation Moutarde 2012-2018 V2023.03.17.xlsx", 
      sheet = "Compilation data_Brut", col_types = c("text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text",  "text", "text", "text", "text", "text", "text", "text", "text", "numeric",  "numeric", "text", "text", "text", "text", "text"))


#### DATASET CURATION AND FORMATING  ####
##### Change column name, taxon, to match common template ####
#names(sbt_m)
#str(sbt_m)
sbt_m = rename(sbt_m,'CE.Deter.Terrain'='CE Déter Terrain','Incomplet'='Imcomplet','Interet'='Intérêt') %>%
  mutate(Code_Taxon =replace(Code_Taxon, Code_Taxon== 'Ali', 'ALI'),
         Annee = factor(Annee),
         ID_Site = factor(ID_Site),
         CE.Deter.Terrain <- factor(CE.Deter.Terrain),
         Code_Taxon <- factor(Code_Taxon),
         GF <- factor(GF),
         Stade <- factor(Stade))
##### Changing the order of the columns, simply for aesthetic reasons ####
enteteVDT=c(
  "Programme",	"Protocole",	"Code_Methode", "Annee", "Date_Prelevement",	"ID_Site", "Code_Parcelle", "Site", "Parcelle", 
  "Modalite", "Cadre", "Sous.cadre", "Bloc", "Repetition",	
  "CE.Deter.Terrain",	"Deter.Code.Taxon",	"Code_Taxon",	"GF", "Stade", "Incomplet",	
  "Nbr_VDT",	"Pds",	"Pds_GF",	
  "Stockage", "Phrase_Fiche_P1",	"Phrase_Fiche_P2", "Commentaires" 
)

#setdiff(enteteVDT,colnames(sbt_m))
#sbt_m[,setdiff(enteteVDT,colnames(sbt_m))]=NA
#sbt_m=sbt_m[,c(enteteVDT)]

# Création AB_cor et BM_cor
#ActiveDataSet("sbt_m")
#Create_T_Brut_Final(data_VDT_input = sbt_m)


