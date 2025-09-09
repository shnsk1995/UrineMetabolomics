

#Fetching Sample Data##############
sampleData <- read.xlsx("data/FFCCS Oxidative Stress Sample Manifest 2-17-2025.xlsx")
colnames(sampleData) <- sampleData[7,]
sampleData <- sampleData[-(1:7),]


###Fetching Raw data###########
nomralizationMethod <- "Raw peak intensity"
plotTitle <- "Raw data"
rawData <- read.xlsx("data/Metabolomics\ result\ for\ oxidative\ stress\ study.xlsx", sheet = nomralizationMethod)
colnames(rawData) <- rawData[1,]
rawData <- rawData[-c(1, 245, 246), ]

sampleData <- sampleData %>%
  dplyr::mutate(SampleGroup = if_else(is.na(`Date of Fire Response`),"Unexposed","Exposed"))

columnNames <- c("SpecimenID","SampleGroup","Sex","Age","SpecificGravity",rawData$compound)


################Formulate DF###############################################
qDataRaw <- rawData


data <- as.data.frame(matrix(ncol = 248,nrow = 0))
colnames(data) <- columnNames


for(i in 1:nrow(sampleData)){
  
  
  datarow <- as.data.frame(matrix(ncol = 248,nrow = 1))
  colnames(datarow) <- columnNames
  
  datarow$SpecimenID[1] <- sampleData$`Specimen ID`[i]
  datarow$SampleGroup[1] <- sampleData$SampleGroup[i]
  datarow$Sex[1] <- sampleData$Sex[i]
  datarow$Age[1] <- sampleData$Age[i]
  datarow$SpecificGravity[1] <- sampleData$`Specific Gravity`[i]
  for (j in 1:nrow(qDataRaw)) {
    
    datarow[1,qDataRaw$compound[j]] <- qDataRaw[qDataRaw$compound == qDataRaw$compound[j],sampleData$`Specimen ID`[i]]
    
  }
  
  data <- rbind(data,datarow)
  
}

# Convert numeric columns (from column 5 onward)
# You can adapt this range as needed
for (i in 4:ncol(data)) {
  # Remove any non-numeric characters (like "*", ">", etc.)
  data[[i]] <- gsub("[^0-9\\.]", "", data[[i]])
  data[[i]] <- as.numeric(data[[i]])
}


##############Normalize raw data for SG#################
##Skip for other datatypes

plotTitle <- "by specific gravity"

# Define target SG (usually 1.020)
target_SG <- 1.020

# Convert SG to numeric (in case it's read as character)
SG_sample <- as.numeric(data$SpecificGravity)

# Extract metabolite data (assuming from column 6 onward)
metabolite_data <- data[, 6:ncol(data)]

# Apply SG normalization
SG_normalized <- sweep(metabolite_data, 1, (target_SG - 1) / (SG_sample - 1), FUN = "*")
rownames(SG_normalized) <- data$SpecimenID

SG_normalized <- log10(SG_normalized + 1e-6)


#####Median_Sample_BoxPlot####################

SG_normalized$SpecimenID <- data$SpecimenID
SG_normalized$SampleGroup <- data$SampleGroup
SG_normalized$Sex <- data$Sex
SG_normalized$Age <- data$Age
SG_normalized$SpecificGravity <- data$SpecificGravity
SG_normalized <- SG_normalized[, c((ncol(SG_normalized)-4):ncol(SG_normalized), 1:(ncol(SG_normalized)-5))]
data <- SG_normalized



samData <- data[,1:5]
data <- log10(data[,6:ncol(data)] + 1e-6)
data <- cbind(samData,data)
rownames(data) <- data$SpecimenID

df <- data


# Assuming first 5 columns are sample info and the rest are metabolites
sample_info_cols <- 1:5
metabolite_cols <- 6:ncol(df)

# Pivot data to long format
df_long <- df %>%
  pivot_longer(
    cols = 6:ncol(df),
    names_to = "Metabolite",
    values_to = "Value"
  )

# Compute median per SpecimenID and make sure it's numeric
medians <- df_long %>%
  group_by(SpecimenID) %>%
  summarise(Median = as.numeric(median(Value, na.rm = TRUE)))

# Plot: bars with a line on top
p <- ggplot(medians, aes(x = SpecimenID, y = Median, group = 1)) +
  geom_col(fill = "lightblue", width = 0.7) +        # ensure visible bars
  geom_line(color = "red", size = 1) +               # line tracing top
  geom_point(color = "red", size = 2) +              # points on top
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = paste0("Median ",plotTitle," Normalized Metabolite per Specimen"),
       x = "Specimen ID",
       y = paste0("Median ",plotTitle," Normalized Metabolite Value"),)+
  theme(
    plot.title = element_text(hjust = 0.5)
  )

# Show plot
print(p)

# Save high-resolution plots
ggsave(paste0("data/",plotTitle,"_median_bar_line_plot.pdf"), p, width = 18, height = 6)
ggsave(paste0("data/",plotTitle,"_median_bar_line_plot.jpeg"), p, width = 18, height = 6, dpi = 300)
