library(data.table)
### metadata
author <- "Izar de Villasante and Marina Vilardell"
description <- "Test project"

### Raw data paths (Genomics Unit):
path_to_Samples_on_array_excel = "W://GENOMICS_UNIT/SampleDB/Samples_on_array_last.xls"
path_to_SDNA.xls = "W://GENOMICS_UNIT/DBB/SDNA.xls"
path_to_storage = "/ijc/LABS/GENOMICS/RAW/Arrays/"
path_to_LTS = "/ijc/LABS/GENOMICS/LTS"


### Folder paths:
results_folder <- "results/"
analysis_folder <- "analysis/"

### Values:
# Path to input sample_sheets:
data_paths <- c(
  ex_EPIC = "inst/extdata/EPIC/sample_sheet_EPIC.rds",
  ex_EPICv2 = "inst/extdata/EPICv2/sample_sheet_EPICv2.rds"
  )

data_names <- c(
  names(data_paths)
)

arraytype <- c(
  "EPIC",
  "EPICv2"
)

project_name<-'Epipe test'

values <- data.table(cbind(data_names,data_paths,arraytype,project_name))
#values <- data.table(data_names="",data_paths="inst/extdata/EPICv2/sample_sheet_EPICv2.rds", arraytype="EPICv2")



# Sample sheet:
# Sentrix_ID must be a variable name in the sample sheet or must be included in the barcode (ex barcode: 207505140099_R08C01)



# Normalization:
norm_function <-c("noob") #

# Filters:
remove_sex <- TRUE
crossreactive=T
snps=T

# Model:
group_var = "Sample_Group" # contains the values of the variable of interest that we want to model.

# Parameters:
## PCA top variable:
topN <- c(100,1000,5000,10000)


values <- values[,.(norm=rlang::syms(norm_function)),by=data_names][values,on=.(data_names)]


### Advanced functions:

### Metadata:
idcol= "Sample_Name" # column with unique sample ids.

### Color palettes:

#"#1B9E77" "#D95F02" "#7570B3" "#E7298A" "#66A61E" "#E6AB02" "#A6761D" "#666666"

pal_discrete =  c("#1B9E77","#7570B3","#D95F02", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666")


#pal_discrete =  c(
 # "#191919", "#0075DC", "#F0A0FF", "#993F00", "#005C31", "#FFE100", "#FF0010",
 # "#FFCC99", "#808080", "#94FFB5", "#8F7C00", "#9DCC00", "#C20088", "#003380",
 #  "#FFA405", "#FFA8BB", "#426600", "#2BCE48", "#4C005C", "#00998F", "#E0FF66",
 # "#740AFF", "#990000", "#FFFF80", "#5EF1F2", "#FF5005"
#)


# models (In case a more specific model is desired). Overrides group_var & covs:

# A) Specific model for each sample_sheet:
# models<- c(
#   m1 = "~ 0 + Condition",
#   m2 = "~ 0 + Type"
# )
# values <- values$models<-models

# B) Multiple models for each sample_sheet:
# models<- c(
#   m1 = "~ 0 + Type * Condition",
#   m2 = "~ 0 + Condition",
#   m3 = "~ 0 + Type"
# )
# values <- values[,.(models=models),by=data_names][values,on=.(data_names)]

# Singular is one against all the rest of variables in Sample_Group
# Pairwise is all against all
# gr: groups together multiple values making the mean.
# Subset: Type == Control, subset the sample_sheet so only a subgroup is left

# Contrasts:
Contrasts <- NULL

# Samp_Group:

# covs:
# Report options:
# Add as many combinations as you would like

report_colab <- tibble::tibble(
  report_name = "report_colab",
  qc_include = "false",
  corrplot_include = "false",
  parameter_tunning_plots = "false",
  code_fold = TRUE,
  show_warning = FALSE,
  shows_message = FALSE

)

report_analyst <- tibble::tibble(
  report_name = "report_analyst",
  qc_include = "true",
  corrplot_include = "true",
  parameter_tunning_plots = "true",
  code_fold = FALSE,
  show_warning = TRUE,
  shows_message = TRUE

)

report_opts <- rbind(report_colab,report_analyst)

report_parameters <- tibble::tibble(
  report_opts,
  path_to_Samples_on_array_excel = path_to_Samples_on_array_excel,
  path_to_SDNA.xls = path_to_SDNA.xls,
  path_to_storage = path_to_storage,
  path_to_LTS = path_to_LTS
)

