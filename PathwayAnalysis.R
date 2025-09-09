

# Set paths
input_file <- "data/Raw_NoCov_PathwayResults/pathway_results.csv"  # replace with your CSV file path
output_dir <- "data/PathwayResults"      # replace with your desired output folder
dir.create(output_dir, showWarnings = FALSE)

# Read the data
df <- read.csv(input_file, row.names = 1)
df$Pathway <- rownames(df)


# ------------------ 2. Dot Plot ------------------
dot_plot <- ggplot(df, aes(x = Impact, y = reorder(Pathway, Impact), size = Hits, color = -log10(Raw.p))) +
  geom_point() +
  theme_minimal(base_size = 18) +
  labs(
    title = "Pathway Enrichment - Dot Plot",
    x = "Impact",
    y = "Pathway",
    size = "Hits",
    color = "-log10(pvalue)"
  ) +
  scale_color_viridis_c() +
  theme(
    axis.text.y = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

# Save dot plot
ggsave(filename = file.path(output_dir, "dot_plot.pdf"), plot = dot_plot, width = 10, height = 8)
ggsave(filename = file.path(output_dir, "dot_plot.jpeg"), plot = dot_plot, width = 10, height = 8, dpi = 300)

