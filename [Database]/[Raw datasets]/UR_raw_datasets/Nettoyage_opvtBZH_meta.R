opvtBZH_rd = read_excel("./[Database]/Raw datasets/opvtBZH_rd.xlsx")


opvtBZH_rd <- opvtBZH_rd %>%
  mutate(Temp = ID_Site,
         ID_Site = Code_Parcelle,
         Code_Parcelle = Temp) %>%
  select(-Temp) 

opvtBZH_rd <- opvtBZH_rd %>%
  mutate(ID_Site = str_replace(ID_Site, "STAC", "STAC-"))

write.csv(opvtBZH_rd, "[Database]/Raw datasets/opvtBZH_rd.csv", row.names = FALSE)




opvtBZH_meta = read_excel("[Database]/[Raw datasets]/Plot_description_datasets/opvtBZH_meta.xlsx", sheet = "meta_19-21")

opvtBZH_meta$GPS_X <- as.numeric(gsub(",", ".", opvtBZH_meta$GPS_X))
opvtBZH_meta$GPS_Y <- as.numeric(gsub(",", ".", opvtBZH_meta$GPS_Y))

opvtBZH_meta <- opvtBZH_meta %>% 
  dplyr::rename(GPS_Y = GPS_X, GPS_X = GPS_Y)

opvtBZH_meta <- opvtBZH_meta %>%
  mutate(Temp = ID_Site,
         ID_Site = Code_Parcelle,
         Code_Parcelle = Temp) %>%
  select(-Temp) 

write.csv(opvtBZH_meta, "[Database]/[Raw datasets]/Plot_description_datasets/opvtBZH_meta.csv", row.names = FALSE)


