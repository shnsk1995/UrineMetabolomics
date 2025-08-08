#Script to perform PCA and plot the PC1 and PC2



#Required Libraries
library(dplyr)
library(limma)
library(edgeR)
library(data.table)
library(ggplot2)
library(openxlsx)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)
library(factoextra)
library(pheatmap)



#Read sample and metabolomics dataset
qDataRaw <- read.xlsx("data/Metabolomics\ result\ for\ oxidative\ stress\ study.xlsx", sheet = "Raw peak intensity")
colnames(qDataRaw) <- qDataRaw[1,]
qDataRaw <- qDataRaw[-c(1, 245, 246), ]

sampleData <- read.xlsx("data/FFCCS Oxidative Stress Sample Manifest 2-17-2025.xlsx")
colnames(sampleData) <- sampleData[7,]
sampleData <- sampleData[-(1:7),]
sampleData <- sampleData %>%
  dplyr::mutate(SampleGroup = if_else(is.na(`Date of Fire Response`),"Unexposed","Exposed"))



#Create a dataframe with both sample and metabolomics dataset
columnNames <- c("SpecimenID","SampleGroup","Sex","Age","SpecificGravity","WeeksAfterExposure","DaysAfterExposure",qDataRaw$compound)


data <- as.data.frame(matrix(ncol = 250,nrow = 0))
colnames(data) <- columnNames


for(i in 1:nrow(sampleData)){
  
  
  datarow <- as.data.frame(matrix(ncol = 250,nrow = 1))
  colnames(datarow) <- columnNames
  
  datarow$SpecimenID[1] <- sampleData$`Specimen ID`[i]
  datarow$SampleGroup[1] <- sampleData$SampleGroup[i]
  datarow$Sex[1] <- sampleData$Sex[i]
  datarow$Age[1] <- sampleData$Age[i]
  datarow$SpecificGravity[1] <- sampleData$`Specific Gravity`[i]
  datarow$WeeksAfterExposure[1] <- sampleData$`Weeks after exposure`[i]
  datarow$DaysAfterExposure[1] <- sampleData$`Days after exposure`[i]
  for (j in 1:nrow(qDataRaw)) {
    
    datarow[1,qDataRaw$compound[j]] <- qDataRaw[qDataRaw$compound == qDataRaw$compound[j],sampleData$`Specimen ID`[i]]
    
  }
  
  data <- rbind(data,datarow)
  
}


# Convert numeric columns (from column 5 onward)
# You can adapt this range as needed
for (i in 5:ncol(data)) {
  # Remove any non-numeric characters (like "*", ">", etc.)
  data[[i]] <- gsub("[^0-9\\.]", "", data[[i]])
  data[[i]] <- as.numeric(data[[i]])
}

write.csv(data, file = "data/FormattedData.csv",row.names = FALSE)

# Define target SG (usually 1.020)
target_SG <- 1.020

# Convert SG to numeric (in case it's read as character)
SG_sample <- as.numeric(data$SpecificGravity)

# Extract metabolite data (assuming from column 8 onward)
metabolite_data <- data[, 8:ncol(data)]

# Apply SG normalization
SG_normalized <- sweep(metabolite_data, 1, (target_SG - 1) / (SG_sample - 1), FUN = "*")


# Log transformation (add small constant to avoid log(0))
log_data <- log10(SG_normalized + 1e-6)

# Autoscaling (mean-center and scale)
scaled_data <- scale(log_data)

# Run PCA
pca_result <- prcomp(scaled_data)

# Combine PCA scores with metadata
pca_df <- as.data.frame(pca_result$x)
pca_df$SpecimenID <- data$SpecimenID
pca_df$SampleGroup <- data$SampleGroup
pca_df$Sex <- data$Sex
pca_df$WeeksAfterExposure <- data$WeeksAfterExposure
pca_df$DaysAfterExposure <- data$DaysAfterExposure



# Define custom colors for PCA plot
custom_colors <- c("1" = "blue", "2" = "red", "3" = "violet", "4" = "orange")

#PCA plot colored by groups based on duration between exposure and collection
#Even for unexposed samples the duration is calculated with the reference to the exposure date
weeks <- ggplot(pca_df, aes(x = PC1, y = PC2, 
                   color = factor(WeeksAfterExposure), 
                   shape = SampleGroup, 
                   label = SpecimenID)) +
  geom_point(size = 3) +
  geom_text_repel(size = 3, max.overlaps = 100) +
  theme_minimal() +
  labs(title = "PCA of SG-Normalized Urine Metabolomics Data",
       x = paste0("PC1 (", round(summary(pca_result)$importance[2, 1] * 100, 1), "% variance)"),
       y = paste0("PC2 (", round(summary(pca_result)$importance[2, 2] * 100, 1), "% variance)"),
       color = "Weeks After Exposure",
       shape = "Sample Group") +
  scale_color_manual(values = custom_colors) +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5, size = 25),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 12))


#PCA plot colored by number of days passed after exposure date.
#Even for unexposed samples the duration is calculated with the reference to the exposure date
Days <- ggplot(pca_df, aes(x = PC1, y = PC2, 
                   color = DaysAfterExposure, 
                   shape = SampleGroup, 
                   label = SpecimenID)) +
  geom_point(size = 3) +
  geom_text_repel(size = 3, max.overlaps = 100) +
  theme_minimal() +
  labs(title = "PCA of SG-Normalized Urine Metabolomics Data",
       x = paste0("PC1 (", round(summary(pca_result)$importance[2, 1] * 100, 1), "% variance)"),
       y = paste0("PC2 (", round(summary(pca_result)$importance[2, 2] * 100, 1), "% variance)"),
       color = "Days After Exposure",
       shape = "Sample Group") +
  scale_color_gradient(low = "blue", high = "red")+
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5, size = 25),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 12))


#Save the plots in plots folder
if(!dir.exists("plots")){
  dir.create("plots")
}
pdf(file = "plots/PCA_SampleGroups_Weeks.pdf", width = 11.69, height = 8.27)

print(weeks)

dev.off()


pdf(file = "plots/PCA_SampleGroups_Days.pdf", width = 11.69, height = 8.27)

print(Days)

dev.off()


jpeg(file = "plots/PCA_SampleGroups_Weeks.jpeg",units = "in", width = 11.69, height = 8.27, quality = 100, res = 600)

print(weeks)

dev.off()


jpeg(file = "plots/PCA_SampleGroups_Days.jpeg",units = "in", width = 11.69, height = 8.27, quality = 100, res = 600)

print(Days)

dev.off()

