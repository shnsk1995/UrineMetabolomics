

# Function: Perform venn comparisons and save plots + Excel
performVennComparisons <- function(datasets, covariates, 
                                   results_dir = "data/LimmaResults", 
                                   out_dir = "data/VennComparison",
                                   p_cutoff = 0.05, 
                                   logFC_cutoff = 0) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Helper to load significant metabolites
  load_sig <- function(dataset, covariate) {
    file <- file.path(results_dir, covariate, dataset,
                      paste0(dataset, "_", covariate, "_LimmaResults.csv"))
    df <- read.csv(file, row.names = 1)
    sig <- rownames(df[df$P.Value < p_cutoff & abs(df$logFC) > logFC_cutoff, ])
    return(sig)
  }
  
  
  
  # Storage for all sig sets
  all_sets <- list()
  for (cov in covariates) {
    for (ds in datasets) {
      all_sets[[paste(ds, cov, sep = "_")]] <- load_sig(ds, cov)
    }
  }
  
  # Function to plot & save venn + Excel
  save_venn <- function(sets, name) {
    # Plot
    p <- ggvenn(sets, show_percentage = FALSE, 
                set_name_size = 5, text_size = 5) +
      ggtitle(name) +
      theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"))
    
    jpeg(file.path(out_dir, paste0("Venn_", name, ".jpeg")), width = 3000, height = 2000, res = 300)
    print(p)
    dev.off()
    ggsave(file.path(out_dir, paste0("Venn_", name, ".pdf")), plot = p, width = 15, height = 10)
    
    # Save common elements to Excel
    wb <- createWorkbook()
    for (i in 1:length(sets)) {
      addWorksheet(wb, names(sets)[i])
      writeData(wb, sheet = names(sets)[i], sets[[i]])
    }
    # Intersection sheet
    addWorksheet(wb, "Common")
    common <- Reduce(intersect, sets)
    writeData(wb, sheet = "Common", common)
    
    saveWorkbook(wb, file.path(out_dir, paste0("Venn_", name, ".xlsx")), overwrite = TRUE)
  }
  
  # 1. Across datasets per covariate
  for (cov in covariates) {
    sets <- all_sets[grep(paste0("_", cov, "$"), names(all_sets))]
    save_venn(sets, paste0("Datasets_", cov))
  }
  
  # 2. Across covariates per dataset
  for (ds in datasets) {
    sets <- all_sets[grep(paste0("^", ds, "_"), names(all_sets))]
    save_venn(sets, paste0("Covariates_", ds))
  }
  
  message("✅ Venn comparisons completed. Results saved in ", out_dir)

  
  return(all_sets)
  
}


datasets <- c("RAW","SPA","VPA","SGN")
covariates <- c("NoCov","DOE","WOE","Age")

allSets <- performVennComparisons(datasets, covariates)


out_dir <- "data/VennComparison"


# 1. Compute overlap matrix
set_names <- names(allSets)
n <- length(allSets)
overlap_mat <- matrix(0, n, n, dimnames = list(set_names, set_names))

for (i in seq_len(n)) {
  for (j in seq_len(n)) {
    overlap_mat[i, j] <- length(intersect(allSets[[i]], allSets[[j]]))
  }
}

# 2. Convert to data.frame for ggplot
df <- as.data.frame(as.table(overlap_mat))
colnames(df) <- c("Set1", "Set2", "Overlap")

# 3. Plot heatmap
g <- ggplot(df, aes(x = Set1, y = Set2, fill = Overlap)) +
  geom_tile(color = "black", size = 0.6) +  # black border
  geom_text(aes(label = Overlap), 
            color = "black", 
            fontface = "bold", 
            size = 5) +  # increase text size
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme_minimal(base_size = 14) +
  ggtitle("Heatmap showing common number of significant genes (pval < 0.05) across all comparisons")+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    axis.title = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

jpeg(file.path(out_dir, paste0("HeatMap_AllResults.jpeg")), width = 4000, height = 2000, res = 300)
print(g)
dev.off()
ggsave(file.path(out_dir, paste0("HeatMap_AllResults.pdf")), plot = g, width = 15, height = 10)



library(ggplot2)
library(stringr)

# Example: allSets is your list of 16 named vectors
# names(allSets) like "RNA_Up", "RNA_Down", "ATAC_Up", "ATAC_Down"

# 1. Create a dataframe with set sizes
df <- data.frame(
  Set = names(allSets),
  Size = sapply(allSets, length)
)

# 2. Extract prefix and suffix from names
df$Dataset <- str_extract(df$Set, "^[^_]+")   # part before "_"
df$Covariate <- str_extract(df$Set, "[^_]+$")   # part after "_"

# Reorder Set factor by Covariate, then by Size within each Covariate (optional)
df$Set <- factor(df$Set, levels = df$Set[order(df$Covariate, df$Size, decreasing = TRUE)])

# ---------- PLOT 1: Group by suffix ----------
p_covariate <- ggplot(df, aes(x = Set, y = Size, fill = Covariate)) +
  geom_col(color = "black") +
  theme_minimal(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(title = "Barplot grouped by Covariate", x = " Dataset", y = "Number of significant metabolites (pval < 0.05)")

jpeg(file.path(out_dir, paste0("BarPlot_AllResults_Covariate.jpeg")), width = 3000, height = 2000, res = 300)
print(p_covariate)
dev.off()
ggsave(file.path(out_dir, paste0("BarPlot_AllResults_Covariate.pdf")), plot = p_covariate, width = 14, height = 10)

# Reorder Set factor by Covariate, then by Size within each Covariate (optional)
df$Set <- factor(df$Set, levels = df$Set[order(df$Dataset, df$Size, decreasing = TRUE)])


# ---------- PLOT 2: Group by prefix ----------
p_dataset <- ggplot(df, aes(x = Set, y = Size, fill = Dataset)) +
  geom_col(color = "black") +
  theme_minimal(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(title = "Barplot grouped by dataset", x = "Dataset", y = "Number of significant metabolites (pval < 0.05)")

jpeg(file.path(out_dir, paste0("BarPlot_AllResults_Dataset.jpeg")), width = 3000, height = 2000, res = 300)
print(p_dataset)
dev.off()
ggsave(file.path(out_dir, paste0("BarPlot_AllResults_Dataset.pdf")), plot = p_dataset, width = 10, height = 10)

