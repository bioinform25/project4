# Upstream transcription-factor regulator analysis for the Mesenchyme
# fibrogenic signature (LUM/THY1/THBS2) found in 04_candidate_gene_localization.R.
#
# Rationale: those 3 genes localize to, and are strongly upregulated within,
# activated Mesenchyme cells (hepatic stellate cells / portal fibroblasts) in
# cirrhosis. This script asks which transcription factors are differentially
# active in that same cirrhotic-vs-healthy Mesenchyme contrast, and which of
# those TFs directly regulate LUM/THY1/THBS2 per the CollecTRI network -- i.e.
# which TFs are plausible upstream drivers of the fibrogenic program, rather
# than downstream effector genes themselves.

source("scripts/config.R")
suppressMessages({
  library(Seurat)
  library(decoupleR)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(httr)
  library(readr)
})

# decoupleR::get_collectri() / OmnipathR::transcriptional() currently error
# ("Join columns in `y` must be present in the data... ncbi_tax_id") because
# the installed OmnipathR release's post-processing code expects a response
# schema the live server no longer returns. Fetching the same CollecTRI data
# directly from the OmniPath REST API and building the source/target/mor
# triplet ourselves (mor = consensus_stimulation - consensus_inhibition, the
# same convention get_collectri() uses internally) avoids the bug entirely.
fetch_collectri <- function() {
  url <- "https://omnipathdb.org/interactions?resources=CollecTRI&genesymbols=1&fields=sources&format=tsv"
  r <- GET(url, timeout(60))
  stopifnot(status_code(r) == 200)
  raw <- read_tsv(I(content(r, "text", encoding = "UTF-8")), show_col_types = FALSE)
  raw %>%
    transmute(source = source_genesymbol, target = target_genesymbol,
              mor = as.integer(consensus_stimulation) - as.integer(consensus_inhibition)) %>%
    filter(mor != 0) %>%
    distinct(source, target, .keep_all = TRUE)
}

obj <- readRDS(file.path(CACHE_DIR, "04_liver_final.rds"))
mes <- subset(obj, cell_type == "Mesenchyme")
cat("Mesenchyme cells:", ncol(mes), "\n")
print(table(mes$patient, mes$condition))

# ---------------------------------------------------------------------------
# Patient-level pseudobulk (avoids treating cells as independent replicates --
# same rationale as project1's power/MDE correction: n=5 vs n=5 patients is
# the real sample size here, not n=cells)
# ---------------------------------------------------------------------------
counts <- GetAssayData(mes, assay = "RNA", layer = "counts")
patient_id <- mes$patient
patients <- unique(patient_id)

pb <- sapply(patients, function(p) Matrix::rowSums(counts[, patient_id == p, drop = FALSE]))
colnames(pb) <- patients
cat("\npseudobulk matrix:", nrow(pb), "genes x", ncol(pb), "patients\n")

# keep genes with reasonable detection across the pseudobulk samples
keep <- rowSums(pb >= 5) >= 5
pb <- pb[keep, ]
cat("genes retained after filtering:", nrow(pb), "\n")

# CPM + log, then per-gene z-score across patients (standard input for
# decoupleR's linear-model (ulm) TF activity method)
cpm <- t(t(pb) / colSums(pb)) * 1e6
logcpm <- log1p(cpm)
z <- t(scale(t(logcpm)))
z[is.na(z)] <- 0

condition <- ifelse(grepl("^cirrhotic", colnames(pb)), "cirrhotic", "healthy")
cat("\npseudobulk samples:", paste(colnames(pb), "=", condition, collapse = ", "), "\n")

# ---------------------------------------------------------------------------
# CollecTRI TF regulon activity (decoupleR)
# ---------------------------------------------------------------------------
net <- fetch_collectri()
cat("\nCollecTRI network: ", nrow(net), "TF-target edges,", length(unique(net$source)), "TFs\n")

acts <- run_ulm(mat = z, net = net, .source = "source", .target = "target",
                 .mor = "mor", minsize = 5)

act_wide <- acts %>%
  filter(statistic == "ulm") %>%
  select(source, condition_sample = condition, score) %>%
  pivot_wider(names_from = condition_sample, values_from = score) %>%
  column_to_rownames("source") %>%
  as.matrix()
# decoupleR's `condition` column here holds the pseudobulk sample name; realign
colnames(act_wide) <- colnames(pb)[match(colnames(act_wide), colnames(pb))]

cirr_cols <- colnames(act_wide)[condition[match(colnames(act_wide), colnames(pb))] == "cirrhotic"]
heal_cols <- colnames(act_wide)[condition[match(colnames(act_wide), colnames(pb))] == "healthy"]

wilcox_p <- apply(act_wide, 1, function(row) {
  tryCatch(wilcox.test(row[cirr_cols], row[heal_cols])$p.value, error = function(e) NA)
})
mean_diff <- rowMeans(act_wide[, cirr_cols, drop = FALSE]) - rowMeans(act_wide[, heal_cols, drop = FALSE])

tf_results <- tibble(
  TF = rownames(act_wide),
  mean_activity_cirrhotic = rowMeans(act_wide[, cirr_cols, drop = FALSE]),
  mean_activity_healthy = rowMeans(act_wide[, heal_cols, drop = FALSE]),
  activity_diff = mean_diff,
  pvalue = wilcox_p
) %>%
  mutate(padj = p.adjust(pvalue, method = "BH")) %>%
  arrange(pvalue)

# which TFs directly target our 3 fibrogenic genes in CollecTRI?
fibro_genes <- c("LUM", "THY1", "THBS2")
tf_to_fibro <- net %>% filter(target %in% fibro_genes)
cat("\nCollecTRI edges into LUM/THY1/THBS2:", nrow(tf_to_fibro), "\n")
print(tf_to_fibro)

tf_results <- tf_results %>%
  left_join(
    tf_to_fibro %>% group_by(source) %>%
      summarise(regulates = paste(target, mor, sep = "(", collapse = "; ")) %>%
      mutate(regulates = paste0(regulates, ")")) %>%
      rename(TF = source),
    by = "TF"
  ) %>%
  mutate(regulates_fibrogenic_gene = !is.na(regulates))

write.csv(tf_results, file.path(TABLE_DIR, "07_tf_activity_mesenchyme.csv"), row.names = FALSE)
write.csv(tf_to_fibro, file.path(TABLE_DIR, "07_collectri_edges_to_candidate_genes.csv"), row.names = FALSE)

cat("\nsignificant TFs (padj<0.05):", sum(tf_results$padj < 0.05, na.rm = TRUE), "\n")
cat("TFs among these that also directly regulate LUM/THY1/THBS2:\n")
print(tf_results %>% filter(padj < 0.05, regulates_fibrogenic_gene) %>%
        select(TF, activity_diff, padj, regulates))

# ---------------------------------------------------------------------------
# Figure: top differentially active TFs, flagging direct regulators of the
# fibrogenic genes
# ---------------------------------------------------------------------------
top_tf <- tf_results %>% filter(!is.na(padj)) %>% arrange(padj) %>% head(20)
p <- ggplot(top_tf, aes(x = reorder(TF, activity_diff), y = activity_diff,
                         fill = regulates_fibrogenic_gene)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = "#d62728", "FALSE" = "grey60"),
                     labels = c("TRUE" = "Regulates LUM/THY1/THBS2 (CollecTRI)",
                                "FALSE" = "Other differentially active TF"),
                     name = NULL) +
  labs(x = NULL, y = "TF activity difference (cirrhotic - healthy)",
       title = "Top differentially active TFs, Mesenchyme (cirrhotic vs healthy)",
       subtitle = "CollecTRI regulon activity (decoupleR::run_ulm), patient-level pseudobulk (n=5 vs 5)") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom")
ggsave(file.path(FIG_DIR, "07_tf_activity_barplot.png"), p, width = 8, height = 7, dpi = 300)

cat("\nDone. Results in results/tables/07_*.csv, figure in results/figures/07_tf_activity_barplot.png\n")
