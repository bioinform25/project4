# GSE136103 — Cell-type localization of 7 fibrosis candidate genes in cirrhotic liver scRNA-seq

Re-analysis of Ramachandran et al. 2019, *Nature* ("Resolving the fibrotic
niche of human liver cirrhosis using single-cell transcriptomics", GSE136103)
to answer a question raised by [project3](https://github.com/bioinform25/project3)
(GSE135251 bulk RNA-seq): the druggability screen there surfaced 7 candidate
genes outside the canonical fibrosis gene list — `CCL21`, `CXCL8`, `CCL20`,
`EPCAM`, `LUM`, `THY1`, `THBS2` — but bulk RNA-seq cannot say *which liver
cell type* actually expresses them. This project uses single-cell data from
the same disease (human liver cirrhosis) to localize each gene.

## Design

- **Data**: 20 human liver 10x samples (5 healthy + 5 cirrhotic patients,
  CD45+/CD45- sorted fractions), 60,925 cells post-QC
  (`nFeature_RNA > 300`, `percent.mt < 30`, per paper Methods). Blood and
  mouse samples in the GEO series were excluded (not relevant to this
  question).
- **Integration**: samples were normalized/scaled/PCA'd together and
  batch-corrected with **Harmony** (`RunHarmony`, batch = individual 10x
  sample), converging in 3 iterations. This is a deliberate upgrade over a
  practice-level prior script that used a naive `merge()` across 10 patients
  with no batch correction — post-Harmony UMAP shows samples well-mixed
  within each cluster (`results/figures/02_umap_sample.png`), confirming
  cluster identity is driven by cell type rather than patient-of-origin.
- **Clustering**: `FindClusters` (Harmony embedding, dims 1:15, resolution
  0.6) → 20 clusters, matching the paper's reported cluster count at the
  same resolution.
- **Annotation**: cluster identity was assigned by scoring each cluster
  against the paper's own supplementary lineage signature gene sets
  (`cell_lineage_signature_genes.xlsx`, parsed programmatically rather than
  hand-retyped) via `AddModuleScore`, taking the highest-scoring lineage per
  cluster. Every cluster had one lineage score clearly separated from the
  rest, and assignments were cross-checked against a marker-gene dotplot
  (`results/figures/03_marker_dotplot_by_cluster.png`) before being finalized
  — the prior practice script stopped at an unlabeled dotplot and never
  actually assigned cell types.
- **Candidate gene localization**: for the 7 genes, percent-expressing and
  average expression were computed per annotated cell type, split by
  condition (healthy vs cirrhotic).
- **Composition shift** (stretch analysis): cell-type proportions compared
  healthy vs cirrhotic (chi-square + per-type proportion tests, BH-adjusted),
  echoing the paper's own "fibrotic niche" composition-shift finding.

## Results

**Cell types identified** (12 lineages, 60,925 cells): T cell, ILC, MP
(mononuclear phagocyte), B cell, Plasma cell, pDC, Mast cell, Endothelia,
Mesenchyme, Hepatocyte, Cholangiocyte, Cycling
(`results/tables/03_cluster_annotation.csv`,
`results/figures/03_umap_annotated.png`).

**Gene localization** (`results/tables/04_candidate_gene_localization_summary.csv`,
`results/figures/04_candidate_genes_dotplot.png`):

| Gene | Localizes to | Cirrhotic vs healthy (% expressing) |
|---|---|---|
| `LUM` | Mesenchyme | 5.1% → 34.3% |
| `THY1` | Mesenchyme | 5.4% → 41.9% |
| `THBS2` | Mesenchyme | 2.4% → 23.9% |
| `EPCAM` | Cholangiocyte | 52.5% → 87.5% |
| `CXCL8` | Cholangiocyte / MP (no clear single winner) | ~22% → ~15-20% (flat/slight decrease) |
| `CCL21` | Mesenchyme / Endothelia (Hepatocyte nominally highest but n=49-262 cells, too small to trust) | 26.6% → (Mesenchyme, healthy not in top 6) |
| `CCL20` | No clear localization — weak signal (≤8%) in every cell type | not condition-differential |

Three of the seven genes (`LUM`, `THY1`, `THBS2`) give a clean, coherent
story: all three localize specifically to the **Mesenchyme** cluster
(hepatic stellate cells / portal fibroblasts, identified by
`PDGFRB`/`ACTA2`/`COL1A1-3`/`DCN`) and are each dramatically more expressed
in cirrhotic than healthy liver within that same cell type. This is
consistent with known stellate-cell activation biology (`THY1`/CD90 is a
established activated-stellate-cell marker) and mechanistically explains
why project3 found these genes elevated in advanced-fibrosis bulk tissue —
they mark or are produced by activated stellate cells, which expand and
activate during fibrosis. `EPCAM` localizes just as clearly to
**Cholangiocytes**, consistent with the cirrhosis-associated "ductular
reaction" (biliary epithelial expansion). `CXCL8` and `CCL20` (chemokines)
do not resolve to one clear producer cell type in this dataset — expression
is diffuse across epithelial and myeloid populations, and `CCL21`'s top hit
sits in a very small hepatocyte subpopulation that should not be
over-interpreted.

**Composition shift** (`results/tables/05_composition_proportion_test.csv`,
`results/figures/05_composition_barplot.png`): Cholangiocyte proportion
nearly triples in cirrhotic liver (3.2% → 9.1% of cells, BH-adj p < 1e-206),
consistent with ductular reaction; Endothelia more than doubles (8.9% →
19.3%), consistent with the paper's reported vascular remodeling.
Notably, **Mesenchyme's own proportion slightly *decreases*** (4.6% → 3.0%)
even as `LUM`/`THY1`/`THBS2` expression *within* that population rises
sharply — i.e., the fibrogenic signal in cirrhosis reflects activation of
existing stellate cells rather than net stellate-cell expansion, a
distinction bulk RNA-seq alone cannot make.

## Pipeline

```
scripts/config.R                       # shared paths + candidate gene list
scripts/00_load_qc.R                   # Read10X, QC filter, per-sample metadata
scripts/01_integration.R               # merge, normalize, PCA, Harmony
scripts/02_clustering.R                # FindClusters + UMAP
scripts/03_annotation.R                # lineage scoring -> cell-type labels
scripts/04_candidate_gene_localization.R  # the core deliverable
scripts/05_composition_shift.R         # healthy vs cirrhotic composition test
```

Run in order from the repo root with `Rscript scripts/0N_....R`. Each script
reads the previous step's cached `.rds` from `data_cache/` (gitignored, not
committed — regenerate by re-running from `00_load_qc.R`).

## Data

Raw 10x files (20 human liver samples, ~270MB) and the paper's supplementary
`cell_lineage_signature_genes.xlsx` are not committed; scripts read them from
a local path (`RAW_DATA_DIR` in `scripts/config.R`) pointing at the extracted
GEO supplementary `GSE136103_RAW.tar`. Blood (4 samples) and mouse (2
samples) folders in the same GEO series were excluded from this analysis.

## Limitations

- Harmony batch correction was run at the level of individual 10x
  sample/chip, not patient — this is the correct granularity for a technical
  batch effect, but means CD45+/CD45- fractions from the same patient are
  still treated as separate batches (as they are separate 10x runs).
- Cell-type assignment is score-based (highest average `AddModuleScore`
  across the paper's own marker sets) rather than manual expert curation;
  it agreed with the marker dotplot for all 20 clusters, but a genuine
  ambiguous population (cluster 3: elevated in both ILC and T-cell markers,
  called ILC) exists and could plausibly be split further (e.g. NKT cells).
- `CCL21`, `CXCL8`, and `CCL20` did not give as clean a localization signal
  as the other 4 genes — this is reported honestly rather than forced into
  a tidy narrative. `CCL21`'s Hepatocyte signal in particular rests on only
  49 (healthy) / 262 (cirrhotic) cells and should be treated as
  hypothesis-generating, not conclusive (per prior lesson on not
  over-trusting small-n subpopulations).
- This dataset and project3's GSE135251 are different patient cohorts (no
  paired samples) — the mechanistic link drawn here (bulk DEG signal ↔
  single-cell source population) is a plausible biological inference, not a
  directly matched within-patient validation.
