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
- **Gut-liver axis extension**: a literature check (below) found that GSE136103
  has already been reused for a gut-liver-axis question (LPS-detoxifying
  enzyme AOAH in hepatocytes/macrophages), and separately that `THBS2`
  activates stellate cells via the TLR4-TGF-β/FAK pathway in an independent
  mouse study — but nobody had checked LPS-receptor expression on the
  Mesenchyme/stellate population in this human cohort. `scripts/06_gutliver_lps_pathway.R`
  checks LPS recognition-complex genes (`LBP`, `CD14`, `TLR4`, `LY96`/MD-2)
  across cell types and tests their cell-level correlation with
  `THBS2`/`LUM`/`THY1` within cirrhotic Mesenchyme cells.

## Literature context (checked before extending the analysis)

Before deciding where to take this project next, recent (2023-2026) literature
using GSE136103 and each of the 7 candidate genes was reviewed:

- GSE136103 has been reused for NAFLD-fibrosis subsetting, an HSC
  heterogeneity study, a fibrogenic-macrophage atlas, deconvolution
  references, a drug-repurposing patent (cathepsin B/H inhibitors — unrelated
  genes), and the gut-liver-axis/AOAH paper above. **None localize this
  specific 7-gene panel** or connect it back to a bulk MASLD fibrosis
  screen — this project's core angle is not duplicated elsewhere.
- `THY1` (CD90) and `EPCAM` are well-established canonical markers (activated
  stellate cell / ductular reaction, respectively) — consistent with, not
  novel relative to, existing literature.
- `THBS2`: an independent 2026 mouse study (LXN-THBS2 axis, GSE174748) shows
  it activates HSCs via TLR4-TGF-β/FAK — supports this project's Mesenchyme
  localization finding.
- `LUM`: known since 2012 as a fibrosis-associated ECM gene, confirmed as a
  stage-3 (myofibroblast) activated-HSC marker in a mouse scRNA-seq study —
  this project's human Mesenchyme localization is a cross-species
  confirmation, not a new discovery.
- `CCL20`: a 2018 in vitro study (LX-2 stellate cell line + palmitic acid
  loading) reported HSCs as the primary CCL20 source in NAFLD. **This
  project's in vivo data disagrees** — Mesenchyme cells show near-zero CCL20
  detection (0.3-0.7% expressing, see Limitations) in this cirrhosis cohort.
  This discrepancy is reported rather than smoothed over; plausible
  explanations include immortalized-cell-line artifact, scRNA-seq dropout for
  a lowly-expressed secreted cytokine, or a genuine difference between the
  NAFLD/palmitic-acid stimulus and mixed-etiology cirrhosis.

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

**Gut-liver axis / LPS receptor expression**
(`results/tables/06_lps_pathway_summary.csv`,
`results/figures/06_lps_pathway_dotplot.png`): as expected, MP
(macrophages/monocytes) is the dominant LPS-sensing population (`CD14`
64-79%, `LY96` 57-69% expressing) and Hepatocyte is the dominant `LBP`
source — both textbook biology, serving as a positive-control sanity check.
In the **Mesenchyme** population, `CD14` and `LY96`(MD-2) expression rise
markedly from healthy to cirrhotic (`CD14` 8.2%→16.9%, `LY96` 4.1%→19.8%);
`TLR4` itself stays low in absolute terms (0.4%→1.3%) but also increases
~3.5-fold. Within cirrhotic Mesenchyme cells, `LY96` shows the clearest
(modest) positive cell-level correlation with the fibrogenic genes (Spearman
ρ = 0.18-0.25 vs `THBS2`/`LUM`/`THY1`;
`results/tables/06_lps_fibrogenic_correlation_mesenchyme_cirrhotic.csv`),
`TLR4` more weakly so (ρ = 0.05-0.15). This is consistent with, and extends
to a human cirrhosis cohort, the published mouse finding that `THBS2` acts
through TLR4 signaling in stellate cells — i.e., activated stellate cells in
cirrhotic liver upregulate the accessory machinery to sense gut-derived LPS
directly, alongside their fibrogenic gene program. The correlations are
modest and `TLR4` transcript levels are low (a known scRNA-seq limitation —
see Limitations), so this is reported as suggestive support for a
gut-liver-axis mechanism, not a confirmed causal link.

## Discussion — what is (and isn't) new here

Most of this project's results **confirm existing literature** rather than
discover something unpublished — worth stating plainly rather than overselling
a re-analysis project as novel science:

- `THY1` → activated stellate cell and `EPCAM` → cholangiocyte are established
  canonical markers.
- `THBS2` and `LUM`'s localization to Mesenchyme/HSCs is already reported in
  independent mouse studies (2026 LXN-THBS2-TLR4 axis paper, 2012+ LUM
  fibrosis literature, mouse HSC-activation scRNA atlases) — this project's
  contribution there is a **human cross-species confirmation**, not a new
  finding.
- The general concept that stellate cells "activate" without necessarily
  expanding in number during fibrosis is also not new to the field, even
  though this project's specific quantification of it for `LUM`/`THY1`/`THBS2`
  on this cohort had not been published before.

Two things found during this project genuinely have not been published
elsewhere, as far as the literature check above could determine:

1. **The CCL20 discrepancy.** A 2018 in vitro study (LX-2 stellate-cell line +
   palmitic acid loading) reported hepatic stellate cells as the primary
   source of CCL20 in NAFLD. This project's in vivo human cirrhosis data
   directly disagrees: Mesenchyme cells show near-zero `CCL20` detection
   (0.3-0.7% expressing) in GSE136103. Nobody appears to have checked this
   specific claim against in vivo single-cell data before. This is reported
   as an **honest contradiction worth flagging**, not a resolved answer —
   plausible explanations (immortalized-cell-line artifact vs. scRNA dropout
   vs. NAFLD/palmitic-acid-specific stimulus not present in mixed-etiology
   cirrhosis) are listed above but not adjudicated here.
2. **HSC-expressed LPS-receptor machinery, in this human cohort.** A 2025
   eLife paper already reused this exact GSE136103 cohort for a
   gut-liver-axis question but only examined AOAH in macrophages/hepatocytes;
   a separate 2026 mouse study showed `THBS2` activates HSCs via TLR4
   signaling, but never checked LPS-receptor expression on stellate cells
   directly. This project is, as far as could be determined, the first to
   check `CD14`/`LY96`/`TLR4` expression specifically in the Mesenchyme
   population of a human cirrhosis single-cell dataset, finding it rises
   sharply in cirrhotic tissue and correlates modestly (Spearman ρ = 0.18-0.25)
   with the `THBS2`/`LUM`/`THY1` fibrogenic signature.

**Both of these should be framed as hypothesis-generating pilot observations,
not established discoveries.** The LPS-receptor correlation is modest in
strength and rests on sparse `TLR4` transcript detection (a known scRNA-seq
limitation, see Limitations) — it suggests the published mouse TLR4-HSC
mechanism *could* extend to human cirrhotic tissue, not that it has been
proven to. The CCL20 contradiction is a real disagreement between two data
sources, not evidence that either one is "wrong." Presenting these as "the
published mouse mechanism might extend to humans" and "a cell-line finding
may not hold in situ" is the honest, defensible framing for a fair
presentation — not "this project discovered a new fibrosis pathway."

## Pipeline

```
scripts/config.R                       # shared paths + candidate gene list
scripts/00_load_qc.R                   # Read10X, QC filter, per-sample metadata
scripts/01_integration.R               # merge, normalize, PCA, Harmony
scripts/02_clustering.R                # FindClusters + UMAP
scripts/03_annotation.R                # lineage scoring -> cell-type labels
scripts/04_candidate_gene_localization.R  # the core deliverable
scripts/05_composition_shift.R         # healthy vs cirrhotic composition test
scripts/06_gutliver_lps_pathway.R      # LPS receptor pathway, gut-liver axis extension
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
- `TLR4` transcript detection is sparse even in MP (the expected
  high-expressing population, ~19-20%), a known scRNA-seq limitation for
  this gene (low mRNA copy number, receptor recycling) — absolute `TLR4` pct
  in Mesenchyme (≤1.3%) likely understates true receptor presence, and the
  LPS-pathway/fibrogenic-gene correlations (ρ ≤ 0.25) should be read as
  suggestive, not confirmatory.
