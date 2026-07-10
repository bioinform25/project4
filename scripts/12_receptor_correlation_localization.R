# Deep-dive on the 4 druggable receptors that emerged from the CellChat
# "signaling gained in cirrhotic Mesenchyme" screen (09/11): CD44, ITGA1,
# ITGB1, EDNRB. Two questions:
#
# 1. Do these receptors track with the fibrogenic signature (LUM/THY1/THBS2)
#    at the single-cell level within the same (cirrhotic Mesenchyme) cells,
#    the same way 06_gutliver_lps_pathway.R tested for the LPS receptor genes?
# 2. The fibrosis literature's best-established stellate-cell-specific
#    integrin is alpha11beta1 (ITGA11), not alpha1beta1 (ITGA1) -- both are
#    present in this dataset, so compare their cell-type localization
#    directly rather than assuming the literature's target is the same as
#    ours.

source("scripts/config.R")
suppressMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(pheatmap)
})

liver_combined <- readRDS(file.path(CACHE_DIR, "04_liver_final.rds"))

RECEPTOR_GENES <- c("CD44", "ITGA1", "ITGB1", "EDNRB", "EDNRA")
FIBRO_GENES <- c("LUM", "THY1", "THBS2")
INTEGRIN_COMPARE <- c("ITGA1", "ITGA8", "ITGA11", "ITGB1")

# ---------------------------------------------------------------------------
# 1) cell-level correlation within cirrhotic Mesenchyme
# ---------------------------------------------------------------------------
mesenchyme_cirr <- subset(liver_combined, cell_type == "Mesenchyme" & condition == "cirrhotic")
cat("cirrhotic Mesenchyme cells:", ncol(mesenchyme_cirr), "\n")

coexpr_mat <- FetchData(mesenchyme_cirr, vars = c(RECEPTOR_GENES, FIBRO_GENES))
cor_mat <- cor(coexpr_mat, method = "spearman")
write.csv(round(cor_mat, 3),
          file.path(TABLE_DIR, "12_receptor_fibrogenic_correlation_mesenchyme_cirrhotic.csv"))

cat("\nSpearman correlation, candidate receptors vs LUM/THY1/THBS2, within cirrhotic Mesenchyme:\n")
print(round(cor_mat[RECEPTOR_GENES, FIBRO_GENES], 3))

p_cor <- pheatmap(cor_mat, display_numbers = TRUE, number_format = "%.2f",
                   main = "Receptor-fibrogenic correlation (cirrhotic Mesenchyme)",
                   filename = file.path(FIG_DIR, "12_receptor_correlation_heatmap.png"),
                   width = 7.0, height = 5.5)

# ---------------------------------------------------------------------------
# 2) ITGA1 vs ITGA8 vs ITGA11 -- which integrin alpha subunit is actually
# most Mesenchyme/cirrhosis-specific in THIS cohort? (ITGA11 is the subunit
# with the strongest stellate-cell-specific fibrosis literature, not ITGA1)
# ---------------------------------------------------------------------------
p_dot <- DotPlot(liver_combined, features = INTEGRIN_COMPARE, group.by = "cell_type",
                  split.by = "condition", cols = c("lightblue", "red")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("ITGA1 vs ITGA8 vs ITGA11 (+ITGB1) localization by cell type")
ggsave(file.path(FIG_DIR, "12_integrin_subunit_comparison_dotplot.png"), p_dot, width = 9, height = 8, dpi = 150)

expr_mat <- FetchData(liver_combined, vars = c(INTEGRIN_COMPARE, "cell_type", "condition"))
integrin_summary <- expr_mat %>%
  pivot_longer(cols = all_of(INTEGRIN_COMPARE), names_to = "gene", values_to = "expr") %>%
  group_by(gene, cell_type, condition) %>%
  summarise(pct_expressing = 100 * mean(expr > 0), avg_expr = mean(expr), n_cells = n(), .groups = "drop") %>%
  arrange(gene, cell_type, condition)
write.csv(integrin_summary, file.path(TABLE_DIR, "12_integrin_subunit_localization_summary.csv"), row.names = FALSE)

cat("\nTop localization per integrin gene (cirrhotic):\n")
top_call <- integrin_summary %>% filter(condition == "cirrhotic") %>%
  group_by(gene) %>% slice_max(pct_expressing, n = 1) %>%
  select(gene, cell_type, pct_expressing, avg_expr)
print(as.data.frame(top_call))

cat("\nMesenchyme-specific pct_expressing, healthy vs cirrhotic, per gene:\n")
print(integrin_summary %>% filter(cell_type == "Mesenchyme") %>%
        select(gene, condition, pct_expressing) %>% as.data.frame())

cat("\nDone. results/tables/12_*.csv, results/figures/12_*.png\n")
