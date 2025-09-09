

# Function to create median plot
PlotMedian <- function(data, sample_info, plot_title, file_prefix) {
  
  # Extract sample columns (columns 6 onwards)
  sample_data <- data[, 6:ncol(data)]
  
  # Ensure numeric
  cnames <- colnames(sample_data)
  sample_data <- as.data.frame(lapply(sample_data, as.numeric))
  sample_data <- log10(sample_data + 1e-6)
  colnames(sample_data) <- cnames
  
  # Compute median per sample
  medians <- apply(sample_data, 2, median, na.rm = TRUE)
  
  # Create dataframe for plotting
  plot_df <- data.frame(
    SpecimenID = colnames(sample_data),
    Median = medians
  )
  
  # Optional: order by SpecimenID as in sampleData
  plot_df$SpecimenID <- factor(plot_df$SpecimenID, levels = sample_info$`Specimen ID`)
  
  # Create ggplot
  p <- ggplot(plot_df, aes(x = SpecimenID, y = Median, group = 1)) +
    geom_col(fill = "skyblue") +
    geom_line(color = "red", size = 1) +
    geom_point(color = "red", size = 2) +
    labs(
      title = plot_title,
      x = "Specimen ID",
      y = "Median Intensity"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 20, face = "bold"), # title in middle
      axis.title = element_text(size = 16, face = "bold"),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10)
    )
  
  # Save as JPEG
  ggsave(paste0(file_prefix, ".jpeg"), p, width = 18, height = 6, units = "in", dpi = 300)
  
  # Save as PDF
  ggsave(paste0(file_prefix, ".pdf"), p, width = 18, height = 6)
  
  return(p)
}

# -----------------------------
# Generate plots for each dataset
# -----------------------------
p_raw <- PlotMedian(rawData, sampleData, "Median per Sample - Raw Data", "data/median_rawData")
p_sum <- PlotMedian(sumPeakArea, sampleData, "Median per Sample - Sum Peak Area Normalized", "data/median_sumPeakArea")
p_vol <- PlotMedian(volume, sampleData, "Median per Sample - Volume Normalized", "data/median_volume")
p_sg <- PlotMedian(sgData, sampleData, "Median per Sample - SG Normalized", "data/median_sgData")

