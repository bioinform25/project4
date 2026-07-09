# Integrate the TF (07) and miRNA (08) upstream-regulator results with
# druggability, and build a combined regulator -> target-gene network figure.
#
# TF druggability is queried from DGIdb (same GraphQL approach as
# project3/scripts/04_target_discovery.R). DGIdb does not meaningfully cover
# miRNA-targeting agents (it is a gene/protein-drug database), so miRNA
# "druggability" here is a small hand-curated table of clinical-stage
# miRNA-modulating compounds relevant to the specific miRNAs that came out of
# 08_mirna_target_analysis.R, verified via literature search rather than a
# programmatic query (see comments for sources).

source("scripts/config.R")
suppressMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(httr)
  library(jsonlite)
  library(igraph)
})

# ---------------------------------------------------------------------------
# TF side: candidates = TFs with a direct CollecTRI edge into LUM/THY1/THBS2
# ---------------------------------------------------------------------------
tf_results <- read_csv(file.path(TABLE_DIR, "07_tf_activity_mesenchyme.csv"), show_col_types = FALSE)
tf_candidates <- tf_results %>% filter(regulates_fibrogenic_gene) %>% arrange(pvalue)
cat("TF candidates (direct CollecTRI regulators of LUM/THY1/THBS2):", nrow(tf_candidates), "\n")
print(tf_candidates %>% select(TF, activity_diff, pvalue, padj, regulates))

gql_query <- function(gene_names) {
  q <- list(
    query = "query genes($names: [String!]) { genes(names: $names) { nodes { name interactions { drug { name } } } } }",
    variables = list(names = as.list(gene_names))
  )
  r <- POST("https://dgidb.org/api/graphql",
            body = toJSON(q, auto_unbox = TRUE),
            content_type("application/json"))
  stopifnot(status_code(r) == 200)
  fromJSON(content(r, "text", encoding = "UTF-8"), simplifyVector = FALSE)
}

dgidb_res <- gql_query(unique(tf_candidates$TF))
gene_nodes <- dgidb_res$data$genes$nodes
tf_drug_summary <- purrr::map_dfr(gene_nodes, function(n) {
  drugs <- unique(purrr::map_chr(n$interactions, ~ .x$drug$name))
  tibble(TF = n$name, n_drugs = length(drugs),
         example_drugs = paste(head(drugs, 5), collapse = "; "))
})

tf_candidates <- tf_candidates %>%
  left_join(tf_drug_summary, by = "TF") %>%
  mutate(n_drugs = ifelse(is.na(n_drugs), 0, n_drugs),
         example_drugs = ifelse(is.na(example_drugs), "", example_drugs),
         nominal_sig = pvalue < 0.05)

write_csv(tf_candidates, file.path(TABLE_DIR, "09_tf_candidates_druggability.csv"))
cat("\nTF druggability:\n")
print(tf_candidates %>% select(TF, activity_diff, pvalue, regulates, n_drugs, example_drugs))

# ---------------------------------------------------------------------------
# miRNA side: hand-curated clinical-stage compound check for the
# fibrosis-family miRNAs flagged in 08_mirna_target_analysis.R. Verified via
# web search 2026-07-10 (ClinicalTrials.gov / primary literature), not a
# database query -- DGIdb does not cover miRNA-targeting drugs.
# ---------------------------------------------------------------------------
mirna_drug_notes <- tribble(
  ~mirna_family, ~compound, ~modality, ~indication, ~status,
  "miR-29", "Remlarsen (MRG-201)", "miR-29 mimic", "Cutaneous fibrosis / keloid (Phase 2, NCT03601052)",
    "Discontinued -- Phase 2 showed severe immune-related adverse events despite on-target ECM/collagen repression in Phase 1",
  "miR-21", "Lademirsen (RG-012)", "anti-miR-21 (LNA)", "Alport syndrome kidney fibrosis (Phase 2, HERA study)",
    "Discontinued 2022 -- stopped for futility (did not meet eGFR-decline efficacy endpoint; safety was acceptable)",
  "miR-34a", "MRX34", "liposomal miR-34a mimic", "Advanced solid tumors incl. liver cancer (Phase 1)",
    "Halted -- serious immune-mediated adverse events including patient deaths; not developed for fibrosis, listed for context since miR-34a targets THBS2/THY1 here",
  "miR-214", "none identified", "-", "-", "No clinical-stage compound found in this search"
)
write_csv(mirna_drug_notes, file.path(TABLE_DIR, "09_mirna_clinical_stage_compounds.csv"))

mirna_hits <- read_csv(file.path(TABLE_DIR, "08_mirna_validated_targets.csv"), show_col_types = FALSE) %>%
  filter(!is.na(known_fibrosis_family))
cat("\nmiRNA hits in known fibrosis families, with clinical-stage compound context:\n")
print(mirna_hits %>% select(mature_mirna_id, target_symbol, known_fibrosis_family))
print(mirna_drug_notes)

# ---------------------------------------------------------------------------
# combined regulator -> target network figure
# ---------------------------------------------------------------------------
tf_edges <- tf_candidates %>%
  separate_rows(regulates, sep = "; ") %>%
  mutate(target = gsub("\\(.*\\)", "", regulates),
         regulator = TF, reg_type = "TF",
         highlight = nominal_sig) %>%
  select(regulator, target, reg_type, highlight)

mirna_edges <- mirna_hits %>%
  distinct(mature_mirna_id, target_symbol, known_fibrosis_family) %>%
  transmute(regulator = mature_mirna_id, target = target_symbol, reg_type = "miRNA",
            highlight = TRUE)

all_edges <- bind_rows(tf_edges, mirna_edges)
write_csv(all_edges, file.path(TABLE_DIR, "09_combined_regulator_network_edges.csv"))

g <- graph_from_data_frame(all_edges[, c("regulator", "target")], directed = TRUE)
node_type <- ifelse(V(g)$name %in% c("LUM", "THY1", "THBS2"), "target",
                     ifelse(V(g)$name %in% tf_edges$regulator, "TF", "miRNA"))
node_color <- c("target" = "#2ca02c", "TF" = "#1f77b4", "miRNA" = "#d62728")[node_type]
edge_highlight <- all_edges$highlight[match(
  paste(get.edgelist(g)[, 1], get.edgelist(g)[, 2]),
  paste(all_edges$regulator, all_edges$target)
)]

png(file.path(FIG_DIR, "09_regulator_network.png"), width = 2600, height = 2200, res = 300)
par(mar = c(0, 0, 2, 0))
plot(g,
     vertex.size = ifelse(node_type == "target", 14, 8),
     vertex.color = node_color,
     vertex.label.cex = ifelse(node_type == "target", 0.9, 0.55),
     vertex.label.color = "black",
     vertex.frame.color = "white",
     edge.color = ifelse(edge_highlight, "black", "grey80"),
     edge.width = ifelse(edge_highlight, 1.6, 0.6),
     edge.arrow.size = 0.25,
     layout = layout_with_fr(g),
     main = "Upstream TF / miRNA regulators of the Mesenchyme fibrogenic signature")
legend("bottomleft",
       legend = c("Candidate gene (LUM/THY1/THBS2)", "Transcription factor", "miRNA"),
       col = c("#2ca02c", "#1f77b4", "#d62728"), pch = 19, bty = "n", cex = 0.8)
dev.off()

cat("\nDone. results/tables/09_*.csv, results/figures/09_regulator_network.png\n")
