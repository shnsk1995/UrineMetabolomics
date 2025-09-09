

###Fetching and cleaning Raw data###########
plotTitle <- "Raw data"
rawData <- read.xlsx("data/Metabolomics\ result\ for\ oxidative\ stress\ study.xlsx", sheet = "Raw peak intensity")
rawDataOSM <- read.xlsx("data/OXIDATIVE\ STRESS\ MARKERS\ IN\ URINE\ SAMPLES.xlsx", sheet = "Raw peak intensities")
colnames(rawData) <- rawData[1,]
rawData <- rawData[-c(1, 245, 246),-c(126,127) ]
rawDataOSM <- rawDataOSM[-c(17, 18, 14, 9), ]

identical(colnames(rawData[,6:ncol(rawData)]), colnames(rawDataOSM[,6:ncol(rawDataOSM)]))

cnames <- c("Compound","Compound ID","Formula","Mass/Charge","Rt.(min)",colnames(rawData[,6:ncol(rawData)]))

colnames(rawData) <- cnames
colnames(rawDataOSM) <- cnames

identical(colnames(rawData), colnames(rawDataOSM))

rawData <- rbind(rawData,rawDataOSM)

rawData <- rawData[!duplicated(rawData$Compound),]

rownames(rawData) <- rawData$Compound



###Fetching and cleaning SumPeakArea data###########
plotTitle <- "Sum Peak Area"
sumPeakArea <- read.xlsx("data/Metabolomics\ result\ for\ oxidative\ stress\ study.xlsx", sheet = "normalized by sum peak area")
sumPeakAreaOSM <- read.xlsx("data/OXIDATIVE\ STRESS\ MARKERS\ IN\ URINE\ SAMPLES.xlsx", sheet = "Normalized by sum peak area")

identical(colnames(sumPeakArea[,6:ncol(sumPeakArea)]), colnames(sumPeakAreaOSM[,6:ncol(sumPeakAreaOSM)]))

cnames <- c("Compound","Compound ID","Formula","Mass/Charge","Rt.(min)",colnames(sumPeakArea[,6:ncol(sumPeakArea)]))

colnames(sumPeakArea) <- cnames
colnames(sumPeakAreaOSM) <- cnames

identical(colnames(sumPeakArea), colnames(sumPeakAreaOSM))

sumPeakArea <- rbind(sumPeakArea,sumPeakAreaOSM)

sumPeakArea <- sumPeakArea[!duplicated(sumPeakArea$Compound),]

rownames(sumPeakArea) <- sumPeakArea$Compound



###Fetching and cleaning volume data###########
plotTitle <- "Volume"
volume <- read.xlsx("data/Metabolomics\ result\ for\ oxidative\ stress\ study.xlsx", sheet = "normalized by volume")
volumeOSM <- read.xlsx("data/OXIDATIVE\ STRESS\ MARKERS\ IN\ URINE\ SAMPLES.xlsx", sheet = "Normalized by volume peak area")

identical(colnames(volume[,6:ncol(volume)]), colnames(volumeOSM[,6:ncol(volumeOSM)]))

cnames <- c("Compound","Compound ID","Formula","Mass/Charge","Rt.(min)",colnames(volume[,6:ncol(volume)]))

colnames(volume) <- cnames
colnames(volumeOSM) <- cnames

identical(colnames(volume), colnames(volumeOSM))

volume <- rbind(volume,volumeOSM)

volume <- volume[!duplicated(volume$Compound),]

rownames(volume) <- volume$Compound




###################Fetch Sample Data########################
sampleData <- read.xlsx("data/FFCCS Oxidative Stress Sample Manifest 2-17-2025.xlsx")
colnames(sampleData) <- sampleData[7,]
sampleData <- sampleData[-(1:7),]
sampleData <- sampleData %>%
  dplyr::mutate(SampleGroup = if_else(is.na(`Date of Fire Response`),"Unexposed","Exposed"))
colnames(sampleData)
write.csv(sampleData,"data/SampleData.csv", row.names = FALSE)



#######SG Normalization######################
sgData <- rawData

# -----------------------------
# 1. Extract metabolite data (sample columns)
# -----------------------------
meta_cols <- sgData[, 1:5]          # first 5 columns: Compound info
sample_data <- sgData[, 6:ncol(sgData)]  # sample columns
cnames <- colnames(sample_data)

# Make sure sample_data is numeric
sample_data <- as.data.frame(lapply(sample_data, as.numeric))
colnames(sample_data) <- cnames

# -----------------------------
# 2. Extract Specific Gravity values
# -----------------------------
sg_values <- sampleData$`Specific Gravity`[match(colnames(sample_data), sampleData$`Specimen ID`)]

# Ensure numeric
sg_values <- as.numeric(sg_values)

# Check alignment
stopifnot(all(colnames(sample_data) == sampleData$`Specimen ID`[match(colnames(sample_data), sampleData$`Specimen ID`)]))

# -----------------------------
# 3. Standard SG normalization
# -----------------------------
SG_ref <- 1.020

sg_normalized <- sweep(sample_data, 2, sg_values, FUN = function(x, sg) {
  x * (SG_ref - 1) / (sg - 1)
})

# -----------------------------
# 4. Combine with metabolite metadata
# -----------------------------
sgData <- cbind(meta_cols, sg_normalized)


#####Save all data into a excel workbook
wb <- createWorkbook()

# Add sheets and write each data frame
addWorksheet(wb, "Raw Data")
writeData(wb, "Raw Data", rawData)

addWorksheet(wb, "SPA Data")
writeData(wb, "SPA Data", sumPeakArea)

addWorksheet(wb, "VPA Data")
writeData(wb, "VPA Data", volume)

addWorksheet(wb, "SGN Data")
writeData(wb, "SGN Data", sgData)

# Save workbook
saveWorkbook(wb, "data/CleanedData.xlsx", overwrite = TRUE)

