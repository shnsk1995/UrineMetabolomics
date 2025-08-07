# Load limma
library(limma)

# Extract metadata and metabolite data
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

