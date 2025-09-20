
# ========================
# Read input files
# ========================
limma_results <- read.csv("data/LimmaResults/NoCov/RAW/RAW_NoCov_LimmaResults.csv", header = TRUE, stringsAsFactors = FALSE, row.names = 1)
limma_results <- data.frame(Metabolite = rownames(limma_results), limma_results)
norm_data <- read.csv("data/MetaboAnalyst/Raw/data_normalized.csv", header = TRUE, stringsAsFactors = FALSE, row.names = 1)

# ========================
# Prepare data
# ========================
# Get top 5 metabolites based on adjusted p-value
top_metabolites <- limma_results %>%
  arrange(P.Value) %>%
  slice(1:5) %>%
  pull(Metabolite)

# Extract labels (first row is group labels)
sample_labels <- as.character(norm_data[1, ])
norm_data <- norm_data[-1, ]  # Remove label row

# ========================
# Create output folders
# ========================
dir.create("data/Box_Plots", recursive = TRUE, showWarnings = FALSE)
dir.create("data/Violin_Plots", recursive = TRUE, showWarnings = FALSE)

# ========================
# Plotting function
# ========================
plot_and_save <- function(metabolite, data_matrix, labels) {
  # Prepare tidy data
  values <- as.numeric(data_matrix[metabolite, ])
  df <- data.frame(
    Sample = colnames(data_matrix),
    Group = labels,
    Value = values
  )
  
  # Titles
  title_text <- paste("Expression of", metabolite)
  
  # ------------------------
  # Box plot
  # ------------------------
  p1 <- ggplot(df, aes(x = Group, y = Value, fill = Group)) +
    geom_boxplot(outlier.shape = 21) +
    geom_jitter(width = 0.2, alpha = 0.7, size = 3) +
    labs(title = title_text, x = "Group", y = "Normalized Expression") +
    theme_bw(base_size = 16) +
    theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
  
  ggsave(filename = paste0("data/Box_Plots/", metabolite, ".jpeg"), plot = p1, width = 8, height = 6, dpi = 300)
  ggsave(filename = paste0("data/Box_Plots/", metabolite, ".pdf"), plot = p1, width = 8, height = 6)
  
  # ------------------------
  # Violin plot
  # ------------------------
  p2 <- ggplot(df, aes(x = Group, y = Value, fill = Group)) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_boxplot(width = 0.1, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.7, size = 3) +
    labs(title = title_text, x = "Group", y = "Normalized Expression") +
    theme_bw(base_size = 16) +
    theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
  
  ggsave(filename = paste0("data/Violin_Plots/", metabolite, ".jpeg"), plot = p2, width = 8, height = 6, dpi = 300)
  ggsave(filename = paste0("data/Violin_Plots/", metabolite, ".pdf"), plot = p2, width = 8, height = 6)
}

# ========================
# Generate plots for top 5 metabolites
# ========================
for (met in top_metabolites) {
  plot_and_save(met, norm_data, sample_labels)
}

