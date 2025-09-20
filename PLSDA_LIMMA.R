
# ========================
# Read input files
# ========================
limma_results <- read.csv("data/LimmaResults/NoCov/RAW/RAW_NoCov_LimmaResults.csv", header = TRUE, stringsAsFactors = FALSE, row.names = 1)
limma_results <- data.frame(Metabolite = rownames(limma_results), limma_results)
plsda_results <- read.csv("data/PLSDA_VIP_RAW_NoCov.csv", header = TRUE, stringsAsFactors = FALSE, row.names = 1)

# ========================
# Step 1: Filter limma significant metabolites
# ========================
sig_metabolites <- limma_results %>%
  filter(P.Value < 0.05) %>%
  pull(Metabolite)

num_sig <- length(sig_metabolites)
cat("Number of significant metabolites from limma:", num_sig, "\n")

# ========================
# Step 2: Get top VIP metabolites based on Comp.1
# ========================
top_vip_metabolites <- plsda_results %>%
  arrange(desc(Comp..1)) %>%
  head(num_sig) %>%
  rownames()

cat("Number of top VIP metabolites selected:", length(top_vip_metabolites), "\n")

# ========================
# Step 3: Compare overlap
# ========================
overlap <- intersect(sig_metabolites, top_vip_metabolites)
only_limma <- setdiff(sig_metabolites, top_vip_metabolites)
only_vip <- setdiff(top_vip_metabolites, sig_metabolites)

# ========================
# Step 4: Save Venn Diagram (PNG + PDF)
# ========================
dir.create("data/Venn_Plots", recursive = TRUE, showWarnings = FALSE)

venn_plot <- venn.diagram(
  x = list(Limma = sig_metabolites, VIP = top_vip_metabolites),
  category.names = c("Limma (P<0.05)", "PLSDA (Top VIP)"),
  filename = NULL,
  output = TRUE,
  imagetype = "png",
  col = "black",
  fill = c("#1f77b4", "#ff7f0e"),
  alpha = 0.5,
  cex = 4,
  cat.cex = 1.8,
  cat.fontface = "bold",
  main = "Venn Diagram: Limma vs VIP",
  main.cex = 2
)

png("data/Venn_Plots/Limma_vs_VIP.png", width = 7500, height = 4500, res = 300)
grid.draw(venn_plot)
dev.off()

pdf("data/Venn_Plots/Limma_vs_VIP.pdf", width = 18, height = 10)
grid.draw(venn_plot)
dev.off()

# ========================
# Step 5: Save Excel file with three sheets
# ========================
wb <- createWorkbook()

addWorksheet(wb, "Overlap")
writeData(wb, "Overlap", data.frame(Overlap = overlap))

addWorksheet(wb, "Only_Limma")
writeData(wb, "Only_Limma", data.frame(Only_Limma = only_limma))

addWorksheet(wb, "Only_VIP")
writeData(wb, "Only_VIP", data.frame(Only_VIP = only_vip))

saveWorkbook(wb, "data/Limma_VIP_Comparison.xlsx", overwrite = TRUE)

cat("✅ Venn diagram saved in data/Venn_Plots/ (PNG + PDF)\n")
cat("✅ Excel file saved: Limma_VIP_Comparison.xlsx\n")
