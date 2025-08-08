library(ggplot2)
library(ggrepel)  # for nice non-overlapping labels

# Example: assuming your limma results dataframe is called "res"
results <- read.csv(file = "data/LimmaDEResults.csv", row.names = 1)
df <- data.frame(Metabolite = rownames(results), results)

# Thresholds
logFC_cutoff <- 0.1
p_cutoff <- 0.05

# Add a significance flag
df$Significant <- with(df,
                       ifelse(abs(logFC) >= logFC_cutoff & P.Value <= p_cutoff,
                              "Significant", "Not Significant"))

# Subset for labeling
df$Label <- ifelse(df$Significant == "Significant", df$Metabolite, NA)

# Create main title and subtitle with thresholds
main_title <- "Volcano Plot"
sub_title <- paste0("Thresholds: |log2FC| ≥ ", logFC_cutoff,
                    ", p ≤ ", p_cutoff)

# Plot
p <- ggplot(df, aes(x = logFC, y = -log10(P.Value), color = Significant)) +
  geom_point(alpha = 0.8, size = 2) +
  geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff),
             linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(p_cutoff),
             linetype = "dashed", color = "black") +
  geom_text_repel(aes(label = Label), size = 4, max.overlaps = Inf) +
  scale_color_manual(values = c("Not Significant" = "grey60",
                                "Significant" = "red")) +
  labs(title = main_title,
       subtitle = sub_title,
       x = "Log2 Fold Change",
       y = "-Log10 P-value") +
  theme_classic(base_size = 18) +
  theme(
    plot.title = element_text(size = 22, face = "bold", hjust = 0.5),  # centered
    plot.subtitle = element_text(size = 18, hjust = 0.5),
    axis.title = element_text(size = 20, face = "bold"),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16)
  )

# Save plot
ggsave("plots/volcano_plot.png", p, width = 8, height = 6, dpi = 300)
ggsave("plots/volcano_plot.pdf", p, width = 8, height = 6)
