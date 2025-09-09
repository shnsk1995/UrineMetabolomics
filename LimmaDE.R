# Load limma
library(limma)

# Extract metadata and metabolite data
data <- read.csv(data, file = "data/FormattedData.csv",row.names = FALSE)
metab_data <- data
metadata <- metab_data[, 1:7]  # Columns: SpecimenID, SampleGroup, etc.
metabolites <- metab_data[, 8:ncol(metab_data)]
rownames(metabolites) <- metab_data$SpecimenID

# Ensure SampleGroup is a factor with "Unexposed" as reference
metadata$SampleGroup <- factor(metadata$SampleGroup, levels = c("Unexposed", "Exposed"))

# Define target SG (usually 1.020)
target_SG <- 1.020

# Convert SG to numeric (in case it's read as character)
SG_sample <- as.numeric(metadata$SpecificGravity)

# Apply SG normalization
SG_normalized <- sweep(metabolites, 1, (target_SG - 1) / (SG_sample - 1), FUN = "*")


# Log transformation (add small constant to avoid log(0))
log_data <- log10(SG_normalized + 1e-6)

# Transpose metabolite data so that features are rows and samples are columns
metabolite_matrix <- t(as.matrix(log_data))

# Confirm row and column correspondence
# rows = metabolites, columns = samples (should match metadata rows)
stopifnot(ncol(metabolite_matrix) == nrow(metadata))

# Build design matrix with SampleGroup and DaysAfterExposure as covariate
design <- model.matrix(~ SampleGroup + DaysAfterExposure, data = metadata)

# Fit model
fit <- lmFit(metabolite_matrix, design)

# Apply empirical Bayes moderation
fit <- eBayes(fit)

# Get results for the group comparison (Exposed vs Unexposed)
results <- topTable(fit, coef = "SampleGroupExposed", number = Inf, adjust.method = "BH")

# View top differentially abundant metabolites
head(results)


#Save results
write.csv(results, file = "data/LimmaDEResults.csv", row.names = TRUE)






#############################################LIMMA-After input from Metaboanalyst###############################################
datasets <- c("RAW","SPA","VPA","SGN")
covariates <- c("NoCov","DOE","WOE","Age")

for(i in 1:length(datasets)){
  for (j in 1:length(covariates)) {
    DoLimma(dataset = datasets[i],covariate = covariates[j])
  }
}

DoLimma <- function(dataset,covariate){
  
  
  # Extract metadata and metabolite data
  data <- read.csv(file = paste0("data/Metaboanalyst/",dataset,"/data_normalized.csv"), row.names = 1)
  metabolites <- data[-1,]
  metabolites <- metabolites %>% mutate(across(everything(), as.numeric))
  metadata <- read.csv("data/SampleData.csv")
  colnames(metabolites) <- gsub("\\.","-",colnames(metabolites))
  rownames(metadata) <- metadata$Sample
  
  # Ensure SampleGroup is a factor with "Unexposed" as reference
  metadata$SampleGroup <- factor(metadata$SampleGroup, levels = c("Unexposed", "Exposed"))
  
  #Order samples in metabolites and meta data
  metadata <- metadata[match(colnames(metabolites),rownames(metadata)),]
  
  # Build design matrix with SampleGroup and DaysAfterExposure as covariate
  if(covariate == "NoCov"){
    design <- model.matrix(~ SampleGroup , data = metadata)
  }else if(covariate == "DOE"){
    design <- model.matrix(~ SampleGroup + Days.after.exposure , data = metadata)
  }else if(covariate == "WOE"){
    design <- model.matrix(~ SampleGroup + Weeks.after.exposure, data = metadata)
  }else if(covariate == "Age"){
    design <- model.matrix(~ SampleGroup + Age, data = metadata)
  }
  
  # Fit model
  fit <- lmFit(metabolites, design)
  
  # Apply empirical Bayes moderation
  fit <- eBayes(fit)
  
  # Get results for the group comparison (Exposed vs Unexposed)
  results <- topTable(fit, coef = "SampleGroupExposed", number = Inf, adjust.method = "BH")
  
  # View top differentially abundant metabolites
  head(results)
  
  if(!dir.exists("data/LimmaResults")){
    dir.create("data/LimmaResults")
  }
  
  if(!dir.exists(paste0("data/LimmaResults/",covariate))){
    dir.create(paste0("data/LimmaResults/",covariate))
  }
  
  if(!dir.exists(paste0("data/LimmaResults/",covariate,"/",dataset))){
    dir.create(paste0("data/LimmaResults/",covariate,"/",dataset))
  }
  
  #Save results
  write.csv(results, file = paste0("data/LimmaResults/",covariate,"/",dataset,"/",dataset,"_",covariate,"_LimmaResults.csv"), row.names = TRUE)
}
