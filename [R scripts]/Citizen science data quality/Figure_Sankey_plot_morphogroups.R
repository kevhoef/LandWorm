#### PREPARATION DATA ####

data_all <- CS_error %>% 
  mutate(Stade = fct_recode(Stade, "AD" = "SA")) %>%
  select(`FIEC`,LIEC, Stade,Nbr_EW) %>% 
  group_by(`FIEC`,LIEC,Stade) %>% 
  mutate(sum=sum(Nbr_EW)) %>% 
  select(-Nbr_EW) %>% 
  distinct() %>% 
  filter(!(LIEC %in% "GF_X"))  %>%
  #filter(!(FIEC %in% "GF_X"))%>%
  arrange(match(LIEC,c("EPI", "ANE_RH", "ANE_BH","END")))%>%
  arrange(match(FIEC,c("EPI", "ANE_RH", "ANE_BH","END")))


tableau_freq <- table(CS_error$FIEC)

# Niveau de facteur à compter
niveau <- "GF_X"

# Vérification si le niveau de facteur existe dans le tableau de fréquences
if (niveau %in% names(tableau_freq)) {
  # Nombre de lignes avec le niveau de facteur spécifié
  nombre_lignes <- tableau_freq[niveau]
  
  # Affichage du résultat
  print(nombre_lignes)
} else {
  print("Le niveau de facteur spécifié n'existe pas.")
}

data_tot <- CS_error %>% 
  select(FIEC,LIEC,Nbr_EW) %>% 
  group_by(FIEC,LIEC) %>% 
  mutate(sum=sum(Nbr_EW)) %>% 
  select(-Nbr_EW) %>% 
  distinct() %>% 
  filter(!(LIEC %in% "GF_X"))  %>%
  #filter(!(FIEC %in% "GF_X"))%>%
  arrange(match(LIEC,c("EPI", "ANE_RH", "ANE_BH","END")))%>%
  arrange(match(FIEC,c("EPI", "ANE_RH", "ANE_BH","END")))%>%
  ungroup() %>%
  pivot_wider(names_from = LIEC,values_from = sum )


data_tot <- as.data.frame(data_tot)
rownames(data_tot) <- data_tot[,1]
data_tot <- data_tot[,-1]
#data<- t(data) ??? dans quel sens ?

#graphes
data_tot_long <- data_tot %>%
  rownames_to_column %>%
  gather(key = 'key', value = 'value', -rowname) %>%
  filter(value > 0)
colnames(data_tot_long) <- c("target", "source", "value")
data_tot_long$target <- paste(data_tot_long$target, " ", sep="")


nodes_ad <- data.frame(name=c(as.character(data_tot_long$source), as.character(data_tot_long$target)) %>% unique())


data_tot_long$IDsource=match(data_tot_long$source, nodes_ad$name)-1 
data_tot_long$IDtarget=match(data_tot_long$target, nodes_ad$name)-1

#ColourScal ='d3.scaleOrdinal() .range(["#FDE725FF","#B4DE2CFF","#6DCD59FF","#35B779FF","grey"])'

# Add a 'group' column to each connection:
data_tot_long$group <- as.factor(c("type_a","type_b","type_c","type_d","type_e","type_a","type_b","type_c","type_d","type_e","type_a","type_b","type_c","type_d","type_e","type_a","type_b","type_c","type_d","type_e"))

# Give a color for each group:
my_color_ad <- 'd3.scaleOrdinal() .domain(["type_a", "type_b", "type_c","type_d","type_e"]) .range(["darkblue", "gold", "seagreen","deeppink", "grey"])'




# Make the Network
sankeyNetwork(Links = data_tot_long, Nodes = nodes_ad,
              Source = "IDtarget", Target = "IDsource",
              Value = "value", NodeID = "name", 
              sinksRight=FALSE,colourScale=my_color_ad, LinkGroup="group",nodeWidth=40, fontSize=13, nodePadding=20)






#### SANKEY DATA AD ####

data_ad <- data_all %>% 
  filter(Stade=="AD") %>% 
  ungroup() %>% 
  select(-Stade) %>% 
  pivot_wider(names_from = LIEC,values_from = sum )


data_ad <- as.data.frame(data_ad)
rownames(data_ad) <- data_ad[,1]
data_ad <- data_ad[,-1]
#data<- t(data) ??? dans quel sens ?

#graphes
data_ad_long <- data_ad %>%
  rownames_to_column %>%
  gather(key = 'key', value = 'value', -rowname) %>%
  filter(value > 0)
colnames(data_ad_long) <- c("source","target", "value")
data_ad_long$target <- paste(data_ad_long$target, " ", sep="")


nodes_ad <- data.frame(name=c(as.character(data_ad_long$source), as.character(data_ad_long$target)) %>% unique())


data_ad_long$IDsource=match(data_ad_long$source, nodes_ad$name)-1 
data_ad_long$IDtarget=match(data_ad_long$target, nodes_ad$name)-1

#ColourScal ='d3.scaleOrdinal() .range(["#FDE725FF","#B4DE2CFF","#6DCD59FF","#35B779FF","grey"])'

# Add a 'group' column to each connection:
data_ad_long$group <- as.factor(c("type_a","type_b","type_c","type_d","type_e","type_a","type_b","type_c","type_d","type_e","type_a","type_b","type_c","type_d","type_e","type_a","type_b","type_c","type_d","type_e"))

# Give a color for each group:
my_color_ad <- 'd3.scaleOrdinal() .domain(["type_a", "type_b", "type_c","type_d","type_e"]) .range(["darkblue", "gold", "seagreen","deeppink", "grey"])'




# Make the Network
sankeyNetwork(Links = data_ad_long, Nodes = nodes_ad,
              Source = "IDsource", Target = "IDtarget",
              Value = "value", NodeID = "name", 
              sinksRight=FALSE,colourScale=my_color_ad, LinkGroup="group",nodeWidth=40, fontSize=13, nodePadding=20)



#### SANKEY DATA JV ####

data_JV <- data_all %>% 
  filter(Stade=="JV") %>% 
  ungroup() %>% 
  select(-Stade) %>% 
  pivot_wider(names_from = LIEC,values_from = sum )

data_JV <- as.data.frame(data_JV)
rownames(data_JV) <- data_JV[,1]
data_JV <- data_JV[,-1]
#data<- t(data) ??? dans quel sens ?

#graphes
data_JV_long <- data_JV %>%
  rownames_to_column %>%
  gather(key = 'key', value = 'value', -rowname) %>%
  filter(value > 0)
colnames(data_JV_long) <- c("source", "target", "value")
data_JV_long$target <- paste(data_JV_long$target, " ", sep="")


nodes_JV <- data.frame(name=c(as.character(data_JV_long$source), as.character(data_JV_long$target)) %>% unique())


data_JV_long$IDsource=match(data_JV_long$source, nodes_JV$name)-1 
data_JV_long$IDtarget=match(data_JV_long$target, nodes_JV$name)-1

#ColourScal ='d3.scaleOrdinal() .range(["#FDE725FF","#B4DE2CFF","#6DCD59FF","#35B779FF","grey"])'

# Add a JV 'group' column to each connection:
data_JV_long$group <- as.factor(c("type_a","type_b","type_c","type_d","type_e","type_a","type_b","type_c","type_d","type_e","type_a","type_b","type_c","type_d","type_e","type_a","type_b","type_c","type_d","type_e"))

# Give a color for each group:
my_color_JV <- 'd3.scaleOrdinal() .domain(["type_a", "type_b", "type_c","type_d","type_e"]) .range(["darkblue", "gold", "seagreen","deeppink", "grey"])'





# Make the Network
sankeyNetwork(Links = data_JV_long, Nodes = nodes_JV,
              Source = "IDsource", Target = "IDtarget",
              Value = "value", NodeID = "name", 
              sinksRight=FALSE,colourScale=my_color_JV, LinkGroup="group",nodeWidth=40, fontSize=13, nodePadding=20, iterations = 0)

