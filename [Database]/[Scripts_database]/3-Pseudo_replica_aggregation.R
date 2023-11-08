#### PACKAGE LOADING ####
library(plyr)

#### DATA AGGREGATION BY PSEUDO REPLICA ####
Create_T_synth_rep <- function( data_VDT_input, nom_tab){ 

  #Picking up the table name
  table_name <- paste0(deparse(substitute(data_VDT_input)), "_rep")
  
  #Identification variable
  Variable_Bandeau <- c("Programme", "Protocole" ,"ID_Site", "Code_Parcelle" ,"Annee", "Modalite", "Bloc", "Repetition")#, "Cadre" j'enleve car cadre est situ? ? l'int?rieur d'une r?p?tition
  
  #Transformation of variables into factors
  data_VDT_input$Programme <- factor(data_VDT_input$Programme) 
  data_VDT_input$Protocole <- factor(data_VDT_input$Protocole) 
  data_VDT_input$Annee <- factor(data_VDT_input$Annee) 
  data_VDT_input$ID_Site <- factor(data_VDT_input$ID_Site) 
  data_VDT_input$Repetition <- factor(data_VDT_input$Repetition)  
  data_VDT_input$Code_Taxon <- factor(data_VDT_input$Code_Taxon)  

  # Create an ID by pasting together the values of columns specified in Variable_Bandeau
  data_VDT_input <- data_VDT_input %>%
    mutate(ID = apply(select(., !!Variable_Bandeau), 1, paste, collapse = "/-/")) %>%
    mutate(ID = as.factor(ID))
  
  #####   Total_AB  ##### 
  # Cross-tabulation
  T_AB_tot_CD <- xtabs(AB_cor ~ ID, data = data_VDT_input, na.action = na.pass)
  # Convert to data frame
  T_AB_tot <- data.frame(matrix(T_AB_tot_CD, ncol = 1, 
                                dimnames = list(rownames(T_AB_tot_CD), "Total_AB")),
                         ID = rownames(T_AB_tot_CD))
  
  #####  Total_AB_by_Stage  ##### 
  # If 'Stage' column exists
  if (!is.null(data_VDT_input$Stade)) {
    # Perform cross-tabulation
    T_AB_STAD_CD <- xtabs(AB_cor ~ ID + Stade, data = data_VDT_input, na.action = na.pass)
    
    # Convert to data frame
    T_AB_STAD <- data.frame(matrix(T_AB_STAD_CD, ncol = ncol(T_AB_STAD_CD),
                                   dimnames = list(rownames(T_AB_STAD_CD),
                                                   paste("AB", colnames(T_AB_STAD_CD), sep = "_"))),
                            ID = rownames(T_AB_STAD_CD))
    
    # Find columns named "AB_X" and prefix them with "STAD_" to distinguish from species "X"
    colnames(T_AB_STAD)[which(colnames(T_AB_STAD) == "AB_X")] <- "AB_STAD_X"
  }
  

  ####  Total_AB_by_Species  ##### 
# If 'Code_Taxon' column exists
if (!is.null(data_VDT_input$Code_Taxon)) {
  # Perform cross-tabulation
  T_AB_ESP_CD <- xtabs(AB_cor ~ ID + Code_Taxon, data = data_VDT_input, na.action = na.pass)
  
  # Convert to data frame
  T_AB_ESP <- data.frame(matrix(T_AB_ESP_CD, ncol = ncol(T_AB_ESP_CD),
                                dimnames = list(rownames(T_AB_ESP_CD),
                                                paste("AB", colnames(T_AB_ESP_CD), sep = "_"))),
                         ID = rownames(T_AB_ESP_CD))
  
  # Find columns named "AB_X" and change them to "AB_SPX" to distinguish from stage "X"
  colnames(T_AB_ESP)[which(colnames(T_AB_ESP) == "AB_X")] <- "AB_SPX"
}

  ### ############################################################  ###
  ### If no Biomass exists (Merge AB tables by column)               ###
  ### ############################################################  ###
  
  # When Biomass (BM_cor) column is not present in the data
  if (is.null(data_VDT_input$BM_cor)) {
    
    # Case: Both Stade and Code_Taxon columns are present
    if (!is.null(data_VDT_input$Stade) & !is.null(data_VDT_input$Code_Taxon)) {
      T_synth_rep_output <- join_all(list(T_AB_tot, T_AB_STAD, T_AB_ESP),
                                     by = "ID")
    } 
    
    # Case: Only Code_Taxon is present, and Stade is not
    else if (!is.null(data_VDT_input$Code_Taxon)) {
      T_synth_rep_output <- join_all(list(T_AB_tot, T_AB_ESP),
                                     by = "ID")
    }
    
    # Case: Only Stade is present, and Code_Taxon is not
    else if (!is.null(data_VDT_input$Stade)) {
      T_synth_rep_output <- join_all(list(T_AB_tot, T_AB_STAD),
                                     by = "ID")
    }
    
    # Case: Neither Stade nor Code_Taxon is present
    else {
      T_synth_rep_output <- join_all(list(T_AB_tot),
                                     by = "ID")
    }

    
    ### *********************************** ###   
    ### If Biomass exists: Cross-tabulate BM ###
    ### *********************************** ###  
    # If Biomass (BM_cor) column is present in the data
  } else {
    
    #####  ------  ##### 
    #####  BM_tot  ##### 
    #####  ------  #####  
    # Cross-tabulation
    T_BM_tot_CD <- xtabs(BM_cor ~ ID, data = data_VDT_input, na.action = na.pass) # na.rm = TRUE
    # Convert to data frame
    T_BM_tot <- data.frame(matrix(T_BM_tot_CD, ncol = 1, dimnames = list(rownames(T_BM_tot_CD), "BM_tot")),
                           ID = rownames(T_BM_tot_CD))

  
    
    #####  ------  ##### 
    #####  BM_ESP  ##### 
    #####  ------  #####  
    # If the Code_Taxon column exists
    if(is.null(data_VDT_input$Code_Taxon) == FALSE) {
      # Perform cross-tabulation
      T_BM_ESP_CD <- xtabs(BM_cor ~ ID + Code_Taxon, data = data_VDT_input, na.action = na.pass)
      
      # Convert to data frame
      T_BM_ESP <- data.frame(matrix(T_BM_ESP_CD, ncol = ncol(T_BM_ESP_CD), 
                                    dimnames = list(rownames(T_BM_ESP_CD),
                                                    paste("BM", colnames(T_BM_ESP_CD), sep = "_"))),
                             ID = rownames(T_BM_ESP_CD))
      
      # Rename columns with 'X' to distinguish them from species columns
      colnames(T_BM_ESP)[which(colnames(T_BM_ESP) == "BM_X")] <- "BM_SPX"
    }
    
    
    ### -----------------------------------------------------------------------------------------  ###
    # Combine the tables to create a synthetic table of results by repetition
    ### -----------------------------------------------------------------------------------------  ###
    
    # If both Stade and Code_Taxon columns exist
    if (is.null(data_VDT_input$Stade) == FALSE & is.null(data_VDT_input$Code_Taxon) == FALSE) {
      T_synth_rep_output <- join_all(list(T_AB_tot, T_BM_tot, T_AB_STAD, T_AB_ESP, T_BM_ESP), 
                                     by = "ID")
      
      # If only Code_Taxon column exists and not Stade
    } else if (is.null(data_VDT_input$Code_Taxon) == FALSE) {
      T_synth_rep_output <- join_all(list(T_AB_tot, T_BM_tot, T_AB_ESP, T_BM_ESP), 
                                     by = "ID")
      
      # If only Stade column exists and not Code_Taxon
    } else if (is.null(data_VDT_input$Stade) == FALSE) {
      T_synth_rep_output <- join_all(list(T_AB_tot, T_BM_tot, T_AB_STAD), 
                                     by = "ID")
      
      # If neither Stade nor Code_Taxon columns exist
    } else { 
      T_synth_rep_output <- join_all(list(T_AB_tot, T_BM_tot), by = "ID") 
    } 
    # END of the BM cross-tabulation
  }
  
  
  ### ########################################################### ###
  # Create the Tab_Bandeau_rep table and merge it with T_synth_rep_output
  ### ########################################################### ###
  
  # Create the table based on the identifier
  Tab_Bandeau_rep <- data.frame(matrix(unlist(strsplit(as.character(T_synth_rep_output$ID), '/-/')), 
                                       dimnames = list(rep("", length(T_synth_rep_output$ID)), Variable_Bandeau),
                                       nrow = length(T_synth_rep_output$ID), byrow = TRUE),
                                stringsAsFactors = TRUE)
  
  # Merge Tab_Bandeau_rep with T_synth_rep_output
  T_synth_rep_output <- cbind(Tab_Bandeau_rep, T_synth_rep_output)
  
  # Delete the ID column
  T_synth_rep_output <- T_synth_rep_output %>% 
    select(-ID)
  
  # Assign the final data frame in the environment
  assign(nom_tab, T_synth_rep_output, envir=parent.frame()) 
  
  # Return the name of the variable where the final data frame is stored
  return(nom_tab)
}

#### LOOP #####
# Get the names of data frames in the global environment that end with "_rd"
data_frames_names <- ls(pattern = "_rd$", envir = .GlobalEnv)

# Create an empty list to store the output tables
output_tables_list <- list()

# Loop through each data frame and apply the function
for (table_name in data_frames_names) {
  nom_tab <- paste0(table_name, "_rep")  # Set the name for the new data frame

  # Use tryCatch to catch any errors that might occur during the function call
  tryCatch({
    # Print the name of the current table being processed
    cat("Processing table:", table_name, "\n")
    
    # Call the Create_T_synth_rep function
    output_table <- Create_T_synth_rep(get(table_name), nom_tab)
    
    # Check if the function returns a valid output
    if (!is.null(nom_tab)) {
      cat("Table", nom_tab, "has been created.\n")
      # Store the output table in the list with a unique name based on the original table name
      output_tables_list[[nom_tab]] <- get(output_table)
    } else {
      cat("Error when creating the table", nom_tab, "\n")
    }
  }, error = function(e) {
    # If an error occurs, print the error message
    cat("Error when processing the table", table_name, ":", conditionMessage(e), "\n")
  })
  
  # Loop through the names of the tables in the list
  for (table_name in names(output_tables_list)) {
    # Use 'assign' to create a variable with the name 'table_name' in the global environment
    # and assign the corresponding table to it
    assign(table_name, output_tables_list[[table_name]], envir = .GlobalEnv)
  }
}
