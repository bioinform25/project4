library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)

source("C:/Users/SAMSUNG/Desktop/project4/scripts/config.R")

liver_combined <- readRDS(file.path(CACHE_DIR, "04_liver_final.rds"))

# LPS recognition complex (gut-derived LPS sensing): LBP shuttles LPS to CD14,
# which hands off to the TLR4/LY96(MD-2) receptor complex -> NF-kB activation.
# Motivation: an independent 2026 mouse study (LXN-THBS2 axis, GSE174748) reports
# THBS2 activates hepatic stellate cells via the TLR4-TGF-beta/FAK pathway; a 2025
# human study on this same GSE136103 cohort examined gut-derived LPS effects on
# hepatocytes (via AOAH) but did not check LPS-receptor expression on stellate
# cells. This script tests whether the Mesenchyme population driving the
# THY1/LUM/THBS2 fibrogenic signal (04_candidate_gene_localization.R) also
# expresses the machinery to directly sense gut-derived LPS.
LPS_PATHWAY_GENES <- c("LBP", "CD14", "TLR4", "LY96")

genes_present <- intersect(LPS_PATHWAY_GENES, rownames(liver_combined))
genes_missing <- setdiff(LPS_PATHWAY_GENES, rownames(liver_combined))
if (length(genes_missing) > 0) cat("WARNING: not detected in data:", paste(genes_missing, collapse = ", "), "\n")

p_dot <- DotPlot(liver_combined, features = genes_present, group.by = "cell_type", split.by = "condition",
                  cols = c("lightblue", "red")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(FIG_DIR, "06_lps_pathway_dotplot.png"), p_dot, width = 8, height = 8, dpi = 150)

p_vln <- VlnPlot(liver_combined, features = genes_present, group.by = "cell_type", split.by = "condition",
                  pt.size = 0, ncol = 2)
ggsave(file.path(FIG_DIR, "06_lps_pathway_violin.png"), p_vln, width = 14, height = 10, dpi = 150)

expr_mat <- FetchData(liver_combined, vars = c(genes_present, "cell_type", "condition"))
summary_df <- expr_mat %>%
  pivot_longer(cols = all_of(genes_present), names_to = "gene", values_to = "expr") %>%
  group_by(gene, cell_type, condition) %>%
  summarise(pct_expressing = 100 * mean(expr > 0), avg_expr = mean(expr), n_cells = n(), .groups = "drop") %>%
  arrange(gene, cell_type, condition)
write.csv(summary_df, file.path(TABLE_DIR, "06_lps_pathway_summary.csv"), row.names = FALSE)

cat("LPS pathway gene expression, Mesenchyme vs MP (classic LPS responder, positive control):\n")
print(summary_df %>% filter(cell_type %in% c("Mesenchyme", "MP")) %>% as.data.frame())

# Within-Mesenchyme co-expression: do cells with higher THBS2/LUM/THY1 also
# show higher TLR4-pathway expression? (cell-level correlation, cirrhotic only,
# where the fibrogenic signature is active)
mesenchyme_cirr <- subset(liver_combined, cell_type == "Mesenchyme" & condition == "cirrhotic")
fibro_genes <- intersect(c("THBS2", "LUM", "THY1"), rownames(mesenchyme_cirr))
coexpr_mat <- FetchData(mesenchyme_cirr, vars = c(genes_present, fibro_genes))
cor_mat <- cor(coexpr_mat, method = "spearman")
write.csv(round(cor_mat, 3), file.path(TABLE_DIR, "06_lps_fibrogenic_correlation_mesenchyme_cirrhotic.csv"))
cat("\nSpearman correlation, LPS-pathway genes vs fibrogenic genes, within cirrhotic Mesenchyme cells (n=",
    ncol(mesenchyme_cirr), "):\n", sep = "")
print(round(cor_mat[genes_present, fibro_genes, drop = FALSE], 3))
