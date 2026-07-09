library(Seurat)
library(dplyr)
library(ggplot2)

source("C:/Users/SAMSUNG/Desktop/project4/scripts/config.R")

liver_combined <- readRDS(file.path(CACHE_DIR, "03_liver_annotated.rds"))

genes_present <- intersect(CANDIDATE_GENES, rownames(liver_combined))
genes_missing <- setdiff(CANDIDATE_GENES, rownames(liver_combined))
if (length(genes_missing) > 0) cat("WARNING: not detected in data:", paste(genes_missing, collapse = ", "), "\n")

liver_combined$cell_type_condition <- paste(liver_combined$cell_type, liver_combined$condition, sep = "_")

# DotPlot across annotated cell types, split by condition
p_dot <- DotPlot(liver_combined, features = genes_present, group.by = "cell_type", split.by = "condition",
                  cols = c("lightblue", "red")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(FIG_DIR, "04_candidate_genes_dotplot.png"), p_dot, width = 10, height = 8, dpi = 150)

p_vln <- VlnPlot(liver_combined, features = genes_present, group.by = "cell_type", split.by = "condition",
                  pt.size = 0, ncol = 2)
ggsave(file.path(FIG_DIR, "04_candidate_genes_violin.png"), p_vln, width = 14, height = 16, dpi = 150)

p_feat <- FeaturePlot(liver_combined, features = genes_present, reduction = "umap", ncol = 3)
ggsave(file.path(FIG_DIR, "04_candidate_genes_featureplot.png"), p_feat, width = 15, height = 10, dpi = 150)

# Summary table: % expressing + average expression per cell type x condition x gene
expr_mat <- FetchData(liver_combined, vars = c(genes_present, "cell_type", "condition"))
summary_df <- expr_mat %>%
  tidyr::pivot_longer(cols = all_of(genes_present), names_to = "gene", values_to = "expr") %>%
  group_by(gene, cell_type, condition) %>%
  summarise(pct_expressing = 100 * mean(expr > 0), avg_expr = mean(expr), n_cells = n(), .groups = "drop") %>%
  arrange(gene, cell_type, condition)

write.csv(summary_df, file.path(TABLE_DIR, "04_candidate_gene_localization_summary.csv"), row.names = FALSE)

cat("\nTop localization per gene (cell type with highest pct_expressing, cirrhotic):\n")
top_call <- summary_df %>%
  filter(condition == "cirrhotic") %>%
  group_by(gene) %>%
  slice_max(pct_expressing, n = 1) %>%
  select(gene, cell_type, pct_expressing, avg_expr)
print(as.data.frame(top_call))

saveRDS(liver_combined, file.path(CACHE_DIR, "04_liver_final.rds"))
