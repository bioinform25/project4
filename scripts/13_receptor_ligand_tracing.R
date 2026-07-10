# Which cell type actually SENDS the signal that CD44 / ITGA1+ITGB1 / EDNRB
# receive on Mesenchyme cells, and how does that change healthy -> cirrhotic?
# Reuses the already-computed 11_mesenchyme_incoming_signaling_all.csv
# (built in 11_cellchat_compare.R from the cached CellChat objects) rather
# than reloading the CellChat objects themselves.

source("scripts/config.R")
suppressMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
})

all_sig <- read_csv(file.path(TABLE_DIR, "11_mesenchyme_incoming_signaling_all.csv"), show_col_types = FALSE)

TARGET_RECEPTORS <- c("CD44", "ITGA1_ITGB1", "EDNRB", "EDNRA")
sub <- all_sig %>% filter(receptor %in% TARGET_RECEPTORS)
cat("edges into Mesenchyme via CD44 / ITGA1+ITGB1 / EDNRB / EDNRA:", nrow(sub), "\n")

# ---------------------------------------------------------------------------
# per receptor: which source cell type(s), which ligand(s), which condition
# ---------------------------------------------------------------------------
trace_summary <- sub %>%
  group_by(receptor, source, condition) %>%
  summarise(n_ligands = n_distinct(ligand), ligands = paste(sort(unique(ligand)), collapse = "; "),
            .groups = "drop") %>%
  arrange(receptor, desc(n_ligands))
write_csv(trace_summary, file.path(TABLE_DIR, "13_receptor_ligand_source_tracing.csv"))

cat("\n--- CD44: source cell types feeding Mesenchyme ---\n")
print(trace_summary %>% filter(receptor == "CD44") %>% as.data.frame())

cat("\n--- ITGA1+ITGB1 (collagen/laminin receptor): source cell types ---\n")
print(trace_summary %>% filter(receptor == "ITGA1_ITGB1") %>% as.data.frame())

cat("\n--- EDN1-EDNRA vs EDN1-EDNRB: does the SOURCE change, or just the receptor? ---\n")
edn <- sub %>% filter(receptor %in% c("EDNRA", "EDNRB")) %>%
  select(source, ligand, receptor, condition, prob, pval)
print(as.data.frame(edn))
cat("\n-> both draw on the same ligand (EDN1) from the same source (Endothelia); only the\n",
    "   RECEPTOR side switches from EDNRA (already significant in healthy) to EDNRB\n",
    "   (newly significant in cirrhotic) -- i.e. Mesenchyme's endothelin-sensing repertoire\n",
    "   itself changes with disease, not just how much EDN1 is around.\n", sep = "")

# ---------------------------------------------------------------------------
# figure: number of distinct incoming ligand signals per source, per receptor,
# healthy vs cirrhotic
# ---------------------------------------------------------------------------
plot_df <- trace_summary %>%
  mutate(receptor = factor(receptor, levels = TARGET_RECEPTORS))
p <- ggplot(plot_df, aes(x = source, y = n_ligands, fill = condition)) +
  geom_col(position = "dodge") +
  facet_wrap(~receptor, scales = "free_x") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = "# distinct ligand signals into Mesenchyme",
       title = "Ligand-source tracing for the 4 candidate receptors")
ggsave(file.path(FIG_DIR, "13_receptor_ligand_source_barplot.png"), p, width = 9, height = 7, dpi = 300)

cat("\nDone. results/tables/13_receptor_ligand_source_tracing.csv, results/figures/13_receptor_ligand_source_barplot.png\n")
