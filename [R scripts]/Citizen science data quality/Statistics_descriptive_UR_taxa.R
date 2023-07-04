########################################################################
##### Calcul UR des ESP par STADE de Developement at national scale ####


sbt_colerror_esp_national = CS_error_simpl %>%  
  group_by(FIEC_LIEC_Stade,Code_Taxon,LIEC)  %>% 
  summarize(Nbr_EW=sum(Nbr_EW) ) 

sbt_colerror_esp = CS_error_simpl %>%  
  group_by(Code_Taxon)  %>% 
  summarize(Nbr_EW=sum(Nbr_EW) ) 
view(sbt_colerror_esp)

sbt_colerror_esp_national_trans = pivot_wider(sbt_colerror_esp_national,
                                        id_cols = c( "Code_Taxon","LIEC"),
                                        names_from = FIEC_LIEC_Stade, 
                                        values_from = Nbr_EW,
                                        values_fill = list(Nbr_EW = 0))


summary(sbt_colerror_esp_national_trans)
sbt_colerror_esp_national_trans = mutate(sbt_colerror_esp_national_trans, 
                                   UR = ifelse (LIEC =="EPI",((EPI_EPI_JV+EPI_EPI_AD+ANE_RH_EPI_JV+ANE_RH_EPI_AD+ANE_BH_EPI_JV+ANE_BH_EPI_AD+END_EPI_JV+END_EPI_AD+GF_X_EPI_JV+GF_X_EPI_AD)-(EPI_EPI_JV+EPI_EPI_AD))/(EPI_EPI_JV+EPI_EPI_AD+ANE_RH_EPI_JV+ANE_RH_EPI_AD+ANE_BH_EPI_JV+ANE_BH_EPI_AD+END_EPI_JV+END_EPI_AD+GF_X_EPI_JV+GF_X_EPI_AD),
                                                ifelse (LIEC=="ANE_RH", ((EPI_ANE_RH_JV+EPI_ANE_RH_AD+ANE_RH_ANE_RH_JV+ANE_RH_ANE_RH_AD+ANE_BH_ANE_RH_JV+ANE_BH_ANE_RH_AD+END_ANE_RH_JV+END_ANE_RH_AD+GF_X_ANE_RH_JV+GF_X_ANE_RH_AD)-(ANE_RH_ANE_RH_JV+ANE_RH_ANE_RH_AD))/(EPI_ANE_RH_JV+EPI_ANE_RH_AD+ANE_RH_ANE_RH_JV+ANE_RH_ANE_RH_AD+ANE_BH_ANE_RH_JV+ANE_BH_ANE_RH_AD+END_ANE_RH_JV+END_ANE_RH_AD+GF_X_ANE_RH_JV+GF_X_ANE_RH_AD),
                                                        ifelse (LIEC == "ANE_BH", ((EPI_ANE_BH_JV+EPI_ANE_BH_AD+ANE_RH_ANE_BH_JV+ANE_RH_ANE_BH_AD+ANE_BH_ANE_BH_JV+ANE_BH_ANE_BH_AD+END_ANE_BH_JV+END_ANE_BH_AD+GF_X_ANE_BH_JV+GF_X_ANE_BH_AD)-(ANE_BH_ANE_BH_JV+ANE_BH_ANE_BH_AD))/(EPI_ANE_BH_JV+EPI_ANE_BH_AD+ANE_RH_ANE_BH_JV+ANE_RH_ANE_BH_AD+ANE_BH_ANE_BH_JV+ANE_BH_ANE_BH_AD+END_ANE_BH_JV+END_ANE_BH_AD+GF_X_ANE_BH_JV+GF_X_ANE_BH_AD),
                                                                ifelse (LIEC == "END", ((EPI_END_JV+EPI_END_AD+ANE_RH_END_JV+ANE_RH_END_AD+ANE_BH_END_JV+ANE_BH_END_AD+END_END_JV+END_END_AD+GF_X_END_JV+GF_X_END_AD)-(END_END_JV+END_END_AD))/(EPI_END_JV+EPI_END_AD+ANE_RH_END_JV+ANE_RH_END_AD+ANE_BH_END_JV+ANE_BH_END_AD+END_END_JV+END_END_AD+GF_X_END_JV+GF_X_END_AD),
                                                                        ifelse(LIEC =="GF_X", NA,NA))))))

sbt_colerror_esp_national_trans = mutate(sbt_colerror_esp_national_trans, 
                                   UR_JV = ifelse (LIEC =="EPI",((EPI_EPI_JV+ANE_RH_EPI_JV+ANE_BH_EPI_JV+END_EPI_JV+GF_X_EPI_JV)-(EPI_EPI_JV))/(EPI_EPI_JV+ANE_RH_EPI_JV+ANE_BH_EPI_JV+END_EPI_JV+GF_X_EPI_JV),
                                                   ifelse (LIEC=="ANE_RH", ((EPI_ANE_RH_JV+ANE_RH_ANE_RH_JV+ANE_BH_ANE_RH_JV+END_ANE_RH_JV+GF_X_ANE_RH_JV)-(ANE_RH_ANE_RH_JV))/(EPI_ANE_RH_JV+ANE_RH_ANE_RH_JV+ANE_BH_ANE_RH_JV+END_ANE_RH_JV+GF_X_ANE_RH_JV),
                                                           ifelse (LIEC == "ANE_BH", ((EPI_ANE_BH_JV+ANE_RH_ANE_BH_JV+ANE_BH_ANE_BH_JV+END_ANE_BH_JV+GF_X_ANE_BH_JV)-(ANE_BH_ANE_BH_JV))/(EPI_ANE_BH_JV+ANE_RH_ANE_BH_JV+ANE_BH_ANE_BH_JV+END_ANE_BH_JV+GF_X_ANE_BH_JV),
                                                                   ifelse (LIEC == "END", ((EPI_END_JV+ANE_RH_END_JV+ANE_BH_END_JV+END_END_JV+GF_X_END_JV)-(END_END_JV))/(EPI_END_JV+ANE_RH_END_JV+ANE_BH_END_JV+END_END_JV+GF_X_END_JV),
                                                                           ifelse(LIEC =="GF_X", NA,NA))))))
sbt_colerror_esp_national_trans = mutate(sbt_colerror_esp_national_trans, 
                                   UR_SA = ifelse (LIEC =="EPI",((EPI_EPI_SA+ANE_RH_EPI_SA+ANE_BH_EPI_SA+END_EPI_SA+GF_X_EPI_SA)-(EPI_EPI_SA))/(EPI_EPI_SA+ANE_RH_EPI_SA+ANE_BH_EPI_SA+END_EPI_SA+GF_X_EPI_SA),
                                                   ifelse (LIEC=="ANE_RH", ((EPI_ANE_RH_SA+ANE_RH_ANE_RH_SA+ANE_BH_ANE_RH_SA+END_ANE_RH_SA+GF_X_ANE_RH_SA)-(ANE_RH_ANE_RH_SA))/(EPI_ANE_RH_SA+ANE_RH_ANE_RH_SA+ANE_BH_ANE_RH_SA+END_ANE_RH_SA+GF_X_ANE_RH_SA),
                                                           ifelse (LIEC == "ANE_BH", ((EPI_ANE_BH_SA+ANE_RH_ANE_BH_SA+ANE_BH_ANE_BH_SA+END_ANE_BH_SA+GF_X_ANE_BH_SA)-(ANE_BH_ANE_BH_SA))/(EPI_ANE_BH_SA+ANE_RH_ANE_BH_SA+ANE_BH_ANE_BH_SA+END_ANE_BH_SA+GF_X_ANE_BH_SA),
                                                                   ifelse (LIEC == "END", ((EPI_END_SA+ANE_RH_END_SA+ANE_BH_END_SA+END_END_SA+GF_X_END_SA)-(END_END_SA))/(EPI_END_SA+ANE_RH_END_SA+ANE_BH_END_SA+END_END_SA+GF_X_END_SA),
                                                                           ifelse(LIEC =="GF_X", NA,NA))))))

sbt_colerror_esp_national_trans = mutate(sbt_colerror_esp_national_trans, 
                                   UR_AD = ifelse (LIEC =="EPI",((EPI_EPI_AD+ANE_RH_EPI_AD+ANE_BH_EPI_AD+END_EPI_AD+GF_X_EPI_AD)-(EPI_EPI_AD))/(EPI_EPI_AD+ANE_RH_EPI_AD+ANE_BH_EPI_AD+END_EPI_AD+GF_X_EPI_AD),
                                                   ifelse (LIEC=="ANE_RH", ((EPI_ANE_RH_AD+ANE_RH_ANE_RH_AD+ANE_BH_ANE_RH_AD+END_ANE_RH_AD+GF_X_ANE_RH_AD)-(ANE_RH_ANE_RH_AD))/(EPI_ANE_RH_AD+ANE_RH_ANE_RH_AD+ANE_BH_ANE_RH_AD+END_ANE_RH_AD+GF_X_ANE_RH_AD),
                                                           ifelse (LIEC == "ANE_BH", ((EPI_ANE_BH_AD+ANE_RH_ANE_BH_AD+ANE_BH_ANE_BH_AD+END_ANE_BH_AD+GF_X_ANE_BH_AD)-(ANE_BH_ANE_BH_AD))/(EPI_ANE_BH_AD+ANE_RH_ANE_BH_AD+ANE_BH_ANE_BH_AD+END_ANE_BH_AD+GF_X_ANE_BH_AD),
                                                                   ifelse (LIEC == "END", ((EPI_END_AD+ANE_RH_END_AD+ANE_BH_END_AD+END_END_AD+GF_X_END_AD)-(END_END_AD))/(EPI_END_AD+ANE_RH_END_AD+ANE_BH_END_AD+END_END_AD+GF_X_END_AD),
                                                                           ifelse(LIEC =="GF_X", NA,NA))))))

sbt_colerror_esp_national_UR <- right_join(sbt_colerror_esp_national_trans, sbt_colerror_esp, by = "Code_Taxon")
sbt_colerror_esp_national_UR =select(sbt_colerror_esp_national_UR, Code_Taxon,UR, UR_JV,UR_SA, UR_AD,Nbr_EW)


write.table(sbt_colerror_esp_national_UR, "D:/Home/khoeffner/Downloads/sbt_colerror_esp_national_UR.csv", sep=",", col.names=TRUE, row.names=TRUE, quote=TRUE, na="NA")
#view(sbt_colerror_esp_national_UR)





sbt_colerror_esp_national_UR