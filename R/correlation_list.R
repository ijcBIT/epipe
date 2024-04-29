#' Correlation analysis
#' 
#' @param clean_object Rgset object
#' @param path Path to save the correlation plot
#' 


# path<-custom_paths$qc_folder
# 
#x<-correlation_analysis(clean_EPICv2,
 #                        path=custom_paths$qc_folder)


# Cal passar param path, perque el volem a correlation_plot. path--qc

correlation_analysis<-function(clean_object,path){
  data_prepared<-preparation(clean_object)
  listed_variables<-correlated_list(data_prepared)
  correlation_plot(listed_variables,data_prepared,path)
  return(listed_variables)
}

## EXAMPLES: 

#x<-correlation_analysis(clean_EPICv2)
# 
# data_prepared<-preparation(clean_EPICv2)
# listed_variables<-correlated_list(data_prepared)
# correlation_plot(listed_variables,data_prepared)






#' Sample sheet check and preparation 
#' 
#' @description Check and prepare the sample sheet for correlation analysis
#' 
#' @param object Clean object (Rgset)
#' @import S4Vectors
#' @import Biobase

# Punt de partida-->clean object (includes sex variable)

preparation<-function(object){
  library(S4Vectors)
  library(Biobase)
  # Just select the variables that may be correlated (not all the variables that we have)
  ###### Ex: Sample_name is irrelevant 
  
  #metadata<-colData(object)
  metadata = object@colData
  
  df<-data.frame(metadata)
  ss_new<-get_category(df)
  
  variables<-names(ss_new)
  category<-attributes(ss_new)$category
  
  # Select variables that are not 'ids'
  non_ids_variables <- variables[category != "ids"]
  
  # Subset the dataframe with selected variables
  ss_new_subset <- ss_new[, non_ids_variables, with = FALSE]
  
  ##################################################################
  
  # Transform character variables into a factor: 
  for (col in names(ss_new_subset)) {
    if (is.character(ss_new_subset[[col]])) {
      ss_new_subset[[col]] <- as.factor(ss_new_subset[[col]])
    }
  }
  
  ###############################################################
  
  # Can happen that some variables have more levels than what they truly have (ex: Type variable: covs, case,control)-->covs is not a level
  
  library(data.table)
  
  # Iterate over each column in the data.table
  for (col in names(ss_new_subset)) {
    # Check if the column is a factor and if "covs" is one of its levels
    if (is.factor(ss_new_subset[[col]]) && "covs" %in% levels(ss_new_subset[[col]])) {
      # Remove the level "covs" from the factor variable
      ss_new_subset[[col]] <- droplevels(ss_new_subset[[col]], exclude = "covs")
    }
    if (is.factor(ss_new_subset[[col]]) && "batch" %in% levels(ss_new_subset[[col]])) {
      # Remove the level "covs" from the factor variable
      ss_new_subset[[col]] <- droplevels(ss_new_subset[[col]], exclude = "batch")
    }
  }
  
  return(ss_new_subset)
  
}







#' List of correlated variables
#'
#' @description Obtain a list of variables that are correlated or its design is not well balanced, depending on each type of variables(numerical/categorical)
#' 
#' @param dataset Sample sheet 
#'
#' @return A list with correlated variables


correlated_list<-function(dataset){


  cor_data<-list()
  p.values<-NULL
  # Loop over each column in the data.table
  for (i in 1:length(names(dataset))) {
    variable1 <- names(dataset)[i]  # Get the name of the first variable
    
    # Start the inner loop from the next column to enforce the order
    for (j in (i + 1):length(names(dataset))) {
      variable2 <- names(dataset)[j]  # Get the name of the second variable
      
      
      ##################################################################
      # # If both variables are numeric--> Pearson correlation test
      if (is.numeric(dataset[[variable1]]) & is.numeric(dataset[[variable2]]) & variable1!=variable2){
        
        cor_score<-cor.test(x=dataset[[variable1]],y=dataset[[variable2]])$estimate
        #print(cor_score)
        if (abs(cor_score)>0.5){
          cor_data <- append(cor_data, list(c(variable1, variable2,cor_score)))
          #cor_data<-c(cor_data,paste0(variable1,' and ', variable2,' are correlated'))
          #p.values<-c(p.values,cor_score)
        }
        
      }
      
      ###########################################################
      # #If variable 1 is numerical and variable 2 is caregorical
      if (is.numeric(dataset[[variable1]]) & is.factor(dataset[[variable2]])){
        
        
        # If categorical variable has two groups --> t.test
        # If categorical variable has more than 2 groups --> anova test
        if (length(levels(dataset[[variable2]]))==2){
          p.value<-t.test(dataset[[variable1]]~dataset[[variable2]])$p.value
          
        }else{
          anova<-aov(dataset[[variable1]]~dataset[[variable2]])
          p.value<-summary(anova)[[1]]$`Pr(>F)`[1]
          
        }
        
        if (p.value<0.05){
          cor_data <- append(cor_data, list(c(variable1, variable2,p.value)))
          #cor_data<-c(cor_data,paste0(variable1,'and', variable2, 'design is not balanced'))
          #p.values<-c(p.values,p.value)
        }
  
      }
      
      #############################################################
      # #If variable 1 is categorical and variable 2 is numerical
      
      if (is.numeric(dataset[[variable2]]) & is.factor(dataset[[variable1]])){
        
        
        # If categorical variable has two groups --> t.test
        # If categorical variable has more than 2 groups --> anova test
        if (length(levels(dataset[[variable1]]))==2){
          p.value<-t.test(dataset[[variable2]]~dataset[[variable1]])$p.value
          
        }else{
          anova<-aov(dataset[[variable2]]~dataset[[variable1]])
          p.value<-summary(anova)[[1]]$`Pr(>F)`[1]
          
        }
        
        if (p.value<0.05){
          cor_data <- append(cor_data, list(c(variable1, variable2,p.value)))
          #cor_data<-c(cor_data,paste0(variable1,'and', variable2, 'design is not balanced'))
          #p.values<-c(p.values,p.value)
        }
        
      }
      
      
      ##########################################################################
      #If both variables are categorical---> chisq.test (multinomial) or test.equal proportions (binomial)
      
      
       if (is.factor(dataset[[variable1]]) & is.factor(dataset[[variable2]]) & variable1!=variable2){
         
         
         # if both variables have 2 levels --> test equal proportions 
         
         if (length(levels(dataset[[variable1]]))==2 & length(levels(dataset[[variable2]]))==2){
           contingency_table<-table(dataset[[variable1]],dataset[[variable2]])
           p.value<-prop.test(contingency_table)$p.value
           if (p.value<0.05){
             cor_data <- append(cor_data, list(c(variable1, variable2,p.value)))
                 #p.values<-c(p.values,p.value)
               }
           
         }
         
         
         # if at least one variable has more than 2 levels --> chisq.test 
         
         if (length(levels(dataset[[variable1]])) > 2 | length(levels(dataset[[variable2]])) > 2) {
           
           # Perform chi-square test of independence
           contingency_table <- table(dataset[[variable1]], dataset[[variable2]])
           chisq_test_result <- chisq.test(contingency_table)
           p.value <- chisq_test_result$p.value
           if (p.value < 0.05) {
             cor_data <- append(cor_data, list(c(variable1, variable2,p.value)))
             #p.values <- c(p.values, p.value)
           }
           
         }
         
         
        
         
       }
      
      
    }
    
  }
  
  return(cor_data)

}



# How to use: 
#coorrelated_l<-correlated_list(tips)
