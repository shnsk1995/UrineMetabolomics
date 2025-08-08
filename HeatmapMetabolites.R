library(pheatmap)
library(dplyr)

# Step 1: Filter top 30 significant metabolites (nominal p < 0.05)
results <- read.csv(file = "data/LimmaDEResults.csv", row.names = 1)
sig_metabs <- results %>%
  filter(P.Value < 0.05) %>%
  arrange(P.Value) %>%
  head(30)

top30_metabs <- rownames(sig_metabs)
writeLines(top30_metabs,"data/sig_meta.txt")

# Step 2: Subset expression matrix
expr_top30 <- metabolite_matrix[top30_metabs, , drop = FALSE]

# Step 3: Match and order metadata with expression matrix columns
# Ensure sample order consistency
metadata <- metadata %>%
  filter(SpecimenID %in% colnames(expr_top30)) %>%
  arrange(SampleGroup)

# Reorder expression matrix columns to match metadata order
expr_top30 <- expr_top30[, metadata$SpecimenID]

# Step 4: Optional Z-score normalization (row-wise)
expr_scaled <- t(scale(t(expr_top30)))
#expr_scaled <- expr_top30

# Step 5: Create top annotation for sample groups
annotation_col <- data.frame(SampleGroup = metadata$SampleGroup)
rownames(annotation_col) <- metadata$SpecimenID

pdf("plots/Top_Heatmap.pdf", width = 20, height = 5)

pheatmap(expr_scaled,
         show_rownames = TRUE,
         show_colnames = TRUE,
         border_color = "black",   # Enables borders
         cellheight = 10,          # Adjust for visibility
         cellwidth = 8,           # Adjust for your number of samples
         fontsize_row = 8,
         fontsize_col = 8,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         annotation_col = annotation_col,
         main = "Top Differential Metabolites (p < 0.05)",
         color = colorRampPalette(c("blue", "white", "red"))(100)
         )
dev.off()


png("plots/Top_Heatmap.png", width = 10500, height = 3000, res = 600)  # Adjust size & resolution

pheatmap(expr_scaled,
         show_rownames = TRUE,
         show_colnames = TRUE,
         border_color = "black",   # Enables borders
         cellheight = 10,          # Adjust for visibility
         cellwidth = 8,           # Adjust for your number of samples
         fontsize_row = 8,
         fontsize_col = 8,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         annotation_col = annotation_col,
         main = "Top Differential Metabolites (p < 0.05)",
         color = colorRampPalette(c("blue", "white", "red"))(100)
)
dev.off()



