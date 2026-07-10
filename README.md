# GSE136103 — Cell-type localization + upstream regulator analysis of fibrosis candidate genes in cirrhotic liver scRNA-seq

Re-analysis of Ramachandran et al. 2019, *Nature* ("Resolving the fibrotic
niche of human liver cirrhosis using single-cell transcriptomics", GSE136103)
to answer a question raised by [project3](https://github.com/bioinform25/project3)
(GSE135251 bulk RNA-seq): the druggability screen there surfaced 7 candidate
genes outside the canonical fibrosis gene list — `CCL21`, `CXCL8`, `CCL20`,
`EPCAM`, `LUM`, `THY1`, `THBS2` — but bulk RNA-seq cannot say *which liver
cell type* actually expresses them. This project uses single-cell data from
the same disease (human liver cirrhosis) to localize each gene, then (in a
2026-07-10 extension, see "Upstream TF / miRNA regulator analysis" below)
asks which transcription factors and miRNAs regulate the 3 genes that
localize most cleanly, as a more pharmacologically-oriented (protein/TF/
miRNA target) follow-up, and finally (same-day second extension, see
"CellChat ligand-receptor cell-cell communication analysis" below) asks
which *other* liver cell types signal into the Mesenchyme population via
ligand-receptor pairs, and how that changes in cirrhosis. A third extension
(see "Receptor deep-dive" below) then digs specifically into the 4 druggable
receptors that extension surfaced (`CD44`, `ITGA1`, `ITGB1`, `EDNRB`).

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

## Upstream TF / miRNA regulator analysis (2026-07-10 extension)

A second follow-up question, motivated by wanting a pharmacologically-oriented
(protein/TF/miRNA) mechanism angle rather than just downstream effector genes:
of the 3 genes that cleanly localize to activated Mesenchyme with a clean
condition effect (`LUM`, `THY1`, `THBS2`), what regulates them, and are any of
those regulators more tractable drug targets than the effector genes
themselves?

- **TF activity** (`scripts/07_tf_regulon_analysis.R`): Mesenchyme cells were
  aggregated to **patient-level pseudobulk** (n=5 healthy, n=5 cirrhotic --
  the real biological replicate count, not the cell count) to avoid
  pseudoreplication. TF regulon activity was inferred with
  `decoupleR::run_ulm()` against the **CollecTRI** transcription-factor
  network (1,188 TFs, fetched directly from the OmniPath REST API -- the
  `decoupleR::get_collectri()` / `OmnipathR::transcriptional()` R wrappers
  currently error against the live server's response schema, a package/server
  version mismatch, not a data problem). TF activity was compared cirrhotic
  vs healthy per TF (Wilcoxon, BH-adjusted), and cross-referenced against
  which TFs have a direct CollecTRI edge into `LUM`/`THY1`/`THBS2`.
- **miRNA targets** (`scripts/08_mirna_target_analysis.R`): `multiMiR` was
  queried for validated (miRTarBase/miRecords/TarBase) and predicted
  (TargetScan) interactions targeting the same 3 genes, then cross-checked
  against a literature-curated list of established fibrosis-associated miRNA
  families (miR-29, -21, -192, -200, -214, -34a, -122, -155).
- **Druggability** (`scripts/09_regulator_druggability.R`): TF candidates
  were queried against DGIdb (same approach as project3). DGIdb does not
  meaningfully cover miRNA-targeting agents, so miRNA "druggability" instead
  used a small hand-curated table of clinical-stage miRNA-modulating
  compounds relevant to the specific miRNA families found, verified via live
  web search (source links in that script's comments), not a database query.

### Results

**TF activity**: with only 5 patients per group, **no TF survives BH
correction** (min padj ≈ 0.52) -- this is an honest limitation of the sample
size, not a null result to hide. Reporting nominal (uncorrected) p-values
instead, among the 22 TFs with a direct CollecTRI edge into `LUM`/`THY1`/
`THBS2` (`results/tables/07_tf_activity_mesenchyme.csv`,
`results/figures/07_tf_activity_barplot.png`):

| TF | Target | Direction (mor) | Activity diff (cirr-healthy) | nominal p | DGIdb drugs |
|---|---|---|---|---|---|
| EOMES | THY1 | activates | +1.43 | 0.0079 | 1 (glatiramer) |
| ZFP42 | THY1 | activates | +3.48 | 0.0079 | 0 |
| SMAD4 | THBS2 | activates | +2.42 | 0.0159 | 5 (incl. cetuximab) |

`SMAD4` is the most mechanistically satisfying hit: it is the canonical
co-SMAD of the TGF-β pathway (already in project3's *canonical* fibrosis gene
list), directly activates `THBS2` per CollecTRI, and is nominally more active
in cirrhotic Mesenchyme -- i.e., this re-derives, rather than discovers, the
textbook TGF-β/SMAD → fibrogenic-ECM axis, which is a useful positive-control
sanity check on the pseudobulk/decoupleR pipeline itself. `EOMES` and `ZFP42`
are the statistically strongest hits but biologically less expected in this
context (EOMES: T-box factor, mainly T/NK-cell differentiation; ZFP42/REX1:
pluripotency factor) -- flagged as an intriguing, not confirmed, lead per the
Limitations below. `SMAD3`/`SMAD2` (the other TGF-β SMADs) also target
`THBS2` in CollecTRI but were not nominally significant here (p=0.22, 0.69).

**miRNA targets**: `LUM`/`THY1`/`THBS2` had 424 validated interactions
combined; the multi-gene ones landing in an established fibrosis-miRNA family
are the most interpretable (`results/tables/08_mirna_validated_targets.csv`,
`09_combined_regulator_network_edges.csv`, `figures/09_regulator_network.png`):

| miRNA | Targets (of LUM/THY1/THBS2) | Family | Clinical-stage compound |
|---|---|---|---|
| hsa-miR-29a/b/c-3p | miR-29b: all 3; miR-29a/c: 2 of 3 | miR-29 (classic anti-fibrotic/anti-collagen family) | Remlarsen (MRG-201), miR-29 mimic -- Phase 2 in cutaneous fibrosis, **discontinued** for immune-related adverse events |
| hsa-miR-21-5p | LUM, THBS2 | miR-21 (classic pro-fibrotic family) | Lademirsen (RG-012), anti-miR-21 -- Phase 2 in Alport kidney fibrosis, **discontinued 2022** for futility (safety was acceptable, efficacy endpoint not met) |
| hsa-miR-34a-5p | THY1, THBS2 | miR-34a | MRX34, liposomal mimic -- Phase 1 in solid tumors (not fibrosis), **halted** after serious immune-mediated adverse events incl. deaths |
| hsa-miR-214-3p | THBS2 | miR-214 | none identified |

Note also that `let-7` family members and `miR-1-3p` technically hit all 3
genes with validated evidence too (`results/tables/08_mirna_multi_target_
candidates.csv`), but these are broadly-expressed, extensively-profiled
miRNAs with large target repertoires in general -- landing atop a
multi-target list is expected for them regardless of fibrosis relevance, so
they are not treated as fibrosis-specific leads the way the miR-29/-21/-34a/
-214 family hits are.

**Taken together**, this analysis's most defensible, novel-ish observation is
that **`miR-29b-3p` has validated evidence for repressing all three
Mesenchyme fibrogenic genes at once** (`LUM`, `THY1`, `THBS2`) -- consistent
with, and mechanistically explaining why, a miR-29 mimic (remlarsen) was
pursued clinically for fibrosis in the first place, even though that specific
compound did not reach approval. The `SMAD4`→`THBS2` TF finding is a
positive-control confirmation of known TGF-β biology rather than a new
finding. The `EOMES`/`ZFP42`→`THY1` TF findings and the `miR-214`→`THBS2`
finding are the least-explored leads but rest on the weakest statistical
footing (nominal p only, n=5 vs 5, no independent replication) and no
clinical-stage compound was identified for `miR-214` in this search.

### Limitations (regulator analysis)

- n=5 vs 5 patients for the pseudobulk TF-activity test is small; the exact
  Wilcoxon test at this sample size has a minimum possible p-value of
  ~0.0079, and nothing survives BH correction across ~1,000 TFs tested. All
  TF findings here are nominal/exploratory, not confirmed.
- CollecTRI and multiMiR (TargetScan/miRTarBase) edges are literature- and
  prediction-derived, not validated in this specific cirrhotic Mesenchyme
  population -- they say a regulatory relationship is plausible/reported
  elsewhere, not that it is active in this exact cell state.
- The miRNA-drug clinical-stage compound table is a small, manually verified
  set focused on the specific miRNA families that came out of this analysis,
  not an exhaustive or database-driven survey -- absence of a listed compound
  (e.g., for miR-214) means none was found in this search, not that none
  exists.
- All three clinical-stage miRNA-modulating compounds identified
  (remlarsen, lademirsen, MRX34) failed to reach approval, for different
  reasons (immune tolerability x2, efficacy futility x1) -- this is reported
  as an important caveat on the drug-class level: "druggable mechanism" here
  should not be read as "clinically de-risked modality."

## CellChat ligand-receptor cell-cell communication analysis (2026-07-10, second extension)

The TF/miRNA analysis above asks what drives `LUM`/`THY1`/`THBS2` expression
from *inside* the Mesenchyme cell. This extension asks the complementary
question with **CellChat**: which *other* liver cell types signal into
Mesenchyme via ligand-receptor pairs, and does that signaling repertoire
change between healthy and cirrhotic liver? This was the CellChat direction
discussed earlier in this project's history and deliberately deferred until
the TF/miRNA analysis was done first.

### Design

- `scripts/10_cellchat_prep.R`: the Seurat object was split by `condition`
  (healthy: 35,074 cells; cirrhotic: 25,851 cells) and a separate CellChat
  object was built for each using **all three** CellChatDB.human categories
  (Secreted Signaling, ECM-Receptor, Cell-Cell Contact — ECM-Receptor
  specifically because `LUM`/`THY1`/`THBS2` are ECM-associated genes).
  `population.size = TRUE` was used because this cohort mixes CD45+/CD45-
  sorted fractions, so raw cell-type proportions in the data do not reflect
  true tissue composition. Standard pipeline: `subsetData` →
  `identifyOverExpressedGenes`/`Interactions` → `computeCommunProb` →
  `filterCommunication(min.cells=10)` → `computeCommunProbPathway` →
  `aggregateNet`.
- `scripts/11_cellchat_compare.R`: merged the two conditions for comparison,
  checked whether `LUM`/`THY1`/`THBS2` themselves appear as ligands/receptors
  in CellChatDB, then specifically extracted all significant signaling
  **into Mesenchyme** (`targets.use = "Mesenchyme"`) per condition, computed
  which specific ligand-receptor pairs are gained/lost in cirrhosis, and
  queried DGIdb for the receptor genes driving the gained signaling.

### Results

**`THBS2` and `THY1` are themselves signaling ligands in CellChatDB** (`LUM`
is not present in CellChatDB as either a ligand or receptor) —
`results/tables/10_candidate_genes_in_cellchatdb.csv`:

| Ligand | Category | Receptors |
|---|---|---|
| `THBS2` | ECM-Receptor | ITGA3+ITGB1, ITGAV+ITGB3, SDC1, SDC4, CD36, **CD47** |
| `THY1` | Cell-Cell Contact | ITGAM+ITGB2, ITGAX+ITGB2, ITGAV+ITGB3 |

This means the Mesenchyme cells expressing more `THBS2`/`THY1` in cirrhosis
are not just passively marked by these genes — they may actively be
*signaling outward* through them (e.g. the `THBS2`-`CD47` "don't-eat-me"
axis is a recognized immune-evasion signal in other CD47-expressing
disease contexts, though this project did not test that specific
hypothesis here — it is flagged as a lead, not demonstrated).

**Global network summary** (`results/figures/10_interaction_count_strength.png`):
total inferred interactions were similar between conditions (healthy 1,512
vs cirrhotic 1,551), but overall interaction *strength* was lower in
cirrhotic (0.518 vs 0.354) — i.e., cirrhotic liver does not have a globally
"louder" communication network, it has a **reorganized** one (61 vs 68
significant pathways detected per condition; `results/figures/
11_pathway_information_flow.png` shows which pathways are healthy- vs
cirrhotic-dominant across the whole network).

**Signaling into Mesenchyme specifically** grows substantially: 46
significant ligand-receptor pairs in healthy vs **76 in cirrhotic**
(`results/tables/11_mesenchyme_incoming_signaling_all.csv`) — 56 pairs are
newly significant in cirrhotic ("gained"), 26 drop out ("lost")
(`results/tables/11_mesenchyme_incoming_{gained,lost}_in_cirrhotic.csv`).
The gained signaling is dominated by **`COLLAGEN` (20 pairs) and `LAMININ`
(14 pairs)**, and critically, 18 of those specific collagen/laminin pairs
have Mesenchyme itself as the *source* — i.e., activated stellate cells
increasingly signal to each other (autocrine/paracrine) via the very ECM
components they produce, a self-reinforcing fibrotic amplification loop.
The remaining gained signal comes from Endothelia (11 pairs, consistent
with this project's earlier finding of near-tripled Endothelia proportion
in cirrhosis), Cholangiocyte (10, consistent with the ductular-reaction
composition shift), and Hepatocyte (9), via `MK`(midkine)/`FN1`/`SPP1`
signaling (`results/figures/11_mesenchyme_incoming_bubble.png`).

**Druggability of the gained-signal receptors**
(`results/tables/11_gained_receptor_druggability.csv`), ranked by how many
gained edges each receptor gene participates in:

| Receptor | Gained edges | DGIdb drugs | Note |
|---|---|---|---|
| `CD44` | 24 | 12 | hyaluronan/`SPP1`/collagen receptor |
| `ITGA1` | 18 | 1 (SAN-300) | collagen receptor subunit |
| `ITGB1` | 18 | 15 | collagen/laminin receptor subunit |
| `LRP1` | 5 | 3 | |
| `EDNRB` | 1 | **29** (incl. **bosentan, ambrisentan** — approved endothelin-receptor antagonists) | only 1 gained edge, but an already-approved drug class exists |
| `PDGFRA` | 1 | 77 (incl. imatinib) | extensively drugged in oncology already |
| `TGFBR2` | 1 | 10 | canonical TGF-β axis, third independent time it surfaces this session (also hit in the TF-regulon and canonical-gene-list analyses) |

`ITGA1`/`ITGB1` (the collagen-receptor integrin pair) and `CD44` are the
clearest mechanistic story here — they are the direct receptors for the
gained collagen/laminin autocrine signal identified above, and both are at
least partially druggable already (an anti-ITGA1 candidate, SAN-300,
exists; CD44 and ITGB1 have broader existing pharmacology). `EDNRB`
appearing at all is notable less for its edge count (only 1) and more
because **approved endothelin-receptor-antagonist drugs already exist**
(bosentan, ambrisentan, for pulmonary arterial hypertension) — endothelin
signaling in liver fibrosis has prior literature support, so this is a
plausible repurposing lead rather than a new discovery.

### Limitations (CellChat analysis)

- CellChat infers communication from average expression per cell-type group
  per condition (not per patient), so — unlike the TF-regulon analysis above
  — this is **not** patient-level pseudoreplication-safe; a result could in
  principle be driven by one or two patients with more cells in a given
  type/condition. This is a known property of the standard CellChat workflow
  (patient-stratified CellChat runs were not attempted here) and should be
  kept in mind alongside the patient-level cell-count table in
  `04_candidate_gene_localization.R`'s composition analysis.
  `Hepatocyte`/`Mast cell`/`pDC` have the fewest cells overall and their
  CellChat edges should be read cautiously for the same reason flagged
  earlier for `CCL21`'s Hepatocyte localization.
- "Gained"/"lost" here means *newly crossing CellChat's significance
  threshold*, not necessarily "absent" vs "present" biologically — a pair
  just below threshold in healthy and just above in cirrhotic is technically
  "gained" but the underlying biology is more continuous than that binary
  framing suggests.
- CellChatDB ligand-receptor annotations (e.g., `THBS2`-`CD47`) are curated
  from the general literature, not validated in this cohort or cell state —
  their appearance here means the interaction is plausible per the
  database, not that it was directly demonstrated to occur.
- As with the TF/miRNA section, these are hypothesis-generating
  observations from a single re-analyzed public cohort, not new
  experimentally validated biology.

## Receptor deep-dive: CD44 / ITGA1+ITGB1 / EDNRB (2026-07-10, third extension)

The CellChat analysis above surfaced 4 druggable receptors receiving newly
gained signaling in cirrhotic Mesenchyme (`CD44`, `ITGA1`, `ITGB1`, `EDNRB`).
The user asked to dig deeper into these specifically. Three questions, using
only data already computed in this repo (no new datasets):

### Design

- `scripts/12_receptor_correlation_localization.R`: (a) within cirrhotic
  Mesenchyme cells, Spearman-correlated `CD44`/`ITGA1`/`ITGB1`/`EDNRB`/`EDNRA`
  against the fibrogenic signature (`LUM`/`THY1`/`THBS2`), same method as the
  LPS-pathway correlation in `06_gutliver_lps_pathway.R`; (b) compared
  `ITGA1` against `ITGA8` and `ITGA11` (both present in this dataset) by
  cell-type localization -- `ITGA11` is the integrin alpha subunit with the
  strongest stellate-cell-specific fibrosis literature, not `ITGA1`, so this
  checks whether the literature's preferred target is even detectable here.
- `scripts/13_receptor_ligand_tracing.R`: reuses the already-computed
  `11_mesenchyme_incoming_signaling_all.csv` (no need to reload the large
  CellChat objects) to trace which cell type sends each receptor's ligand(s),
  and whether that source changes between healthy and cirrhotic.
- `scripts/14_receptor_clinical_context.R`: hand-curated, web-search-verified
  (2026-07-10) clinical/druggability context per receptor -- same rigor as
  the miRNA compound check in `09_regulator_druggability.R`, i.e. checking
  whether a real compound targets the *exact* gene/subtype that changed here,
  not just "a drug exists for this gene family."

### Results

**Correlation** (`results/tables/12_receptor_fibrogenic_correlation_mesenchyme_cirrhotic.csv`,
`figures/12_receptor_correlation_heatmap.png`): `EDNRB` is the only receptor
with a clean positive correlation with all 3 fibrogenic genes (LUM 0.22, THY1
0.16, THBS2 0.07), while `EDNRA` is *negatively* correlated with all 3 (-0.13
to -0.08) -- an independent confirmation, via a completely different method
(cell-level correlation vs. CellChat's differential-significance test), that
`EDNRB` specifically (not `EDNRA`) tracks with fibrogenic activation. `CD44`,
`ITGA1`, `ITGB1` show only weak correlations (0.01-0.29) with the fibrogenic
genes and with each other.

**Integrin subunit comparison** (`results/tables/12_integrin_subunit_localization_summary.csv`,
`figures/12_integrin_subunit_comparison_dotplot.png`) -- a genuinely
unexpected result: `ITGA11`, the subunit with by far the strongest
stellate-cell-specific fibrosis literature (hedgehog-ITGA11 axis,
miR-12135/ITGA11 axis), is barely detected in this dataset (5.2% of cirrhotic
Mesenchyme cells, up from 1.5% healthy) -- far below `ITGA1` (41.4%, up from
23.9%). `ITGA8` actually *decreases* in Mesenchyme (15.0% -> 4.6%) and its
top-localizing cell type is Plasma cell, not Mesenchyme, directly
contradicting the assumption it would behave like a stellate marker here.
This is reported as an honest discrepancy: either (a) `ITGA11`'s established
role doesn't transfer well from the mouse/culture-activated-HSC literature to
this human in-vivo cirrhosis cohort, or (b) `ITGA11` suffers the same
low-copy-number scRNA-seq dropout already flagged for `TLR4` in
`06_gutliver_lps_pathway.R`. Either way, this project's own hit (`ITGA1`) is
*not* the field's preferred integrin target, and that gap is stated plainly
rather than glossed over.

**Ligand-source tracing** (`results/tables/13_receptor_ligand_source_tracing.csv`,
`figures/13_receptor_ligand_source_barplot.png`): all `CD44` and
`ITGA1+ITGB1` signaling into Mesenchyme is exclusively cirrhotic (zero
healthy-condition edges for either receptor) -- Mesenchyme itself is the
largest single source (autocrine, 9-10 distinct ligands: collagens,
laminins, `FN1`), with Endothelia/Cholangiocyte/Hepatocyte contributing the
rest paracrinely. For the endothelin axis, both `EDN1-EDNRA` (healthy) and
`EDN1-EDNRB` (cirrhotic) draw on the *same* ligand from the *same* source
(Endothelia) -- only the receptor side switches, meaning Mesenchyme's
endothelin-sensing repertoire itself changes with disease, not simply how
much `EDN1` is around.

**Clinical context** (`results/tables/14_receptor_clinical_context.csv`) --
the honest bottom line: **none of the 4 receptors has a compound that is
both actively in clinical development and precisely matched to the specific
gene/subtype that actually changed in this analysis**:

| Gene (our hit) | Compound checked | Gap |
|---|---|---|
| `ITGA1`/`ITGB1` | SAN-300 (anti-VLA-1) | Phase 1 in rheumatoid arthritis only, never tested in fibrosis |
| `ITGA1`/`ITGB1` | (n/a) | The field's real target is `ITGA11`, a different alpha subunit, barely detected in our data |
| `CD44` | RG7356 (anti-CD44) | Terminated early in oncology (no dose-response); not a fibrosis compound |
| `EDNRB` | Bosentan / Ambrisentan | Bosentan carries hepatotoxicity risk; ambrisentan is ETA-selective, doesn't engage EDNRB (ETB) |
| `EDNRA` | Zibotentan | Actively in an ongoing cirrhosis trial (ZEAL-UNLOCK) -- but targets the receptor already active in *healthy* tissue, not the EDNRB signal that is new in cirrhosis |

### Limitations (receptor deep-dive)

- The correlation and ligand-tracing analyses are internally consistent with
  each other and with the earlier CellChat/TF results, but all derive from
  the same single GSE136103 cohort -- this is convergent evidence within one
  dataset, not independent replication.
- `ITGA11`/`TLR4`-type low-count scRNA-seq dropout is a recurring caveat in
  this project (also flagged in `06_gutliver_lps_pathway.R`) and cannot be
  fully distinguished from genuine low expression without deeper-coverage
  data (e.g., targeted qPCR or a snRNA-seq platform with better capture for
  low-copy transcripts).
- The clinical-context table is a snapshot of a small, manually verified set
  of compounds (2026-07-10), not an exhaustive drug-repurposing screen.

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
scripts/07_tf_regulon_analysis.R       # CollecTRI TF activity, Mesenchyme pseudobulk
scripts/08_mirna_target_analysis.R     # multiMiR miRNA targets of LUM/THY1/THBS2
scripts/09_regulator_druggability.R    # DGIdb (TF) + literature (miRNA) druggability, network figure
scripts/10_cellchat_prep.R             # per-condition CellChat objects (healthy, cirrhotic)
scripts/11_cellchat_compare.R          # merged comparison, signaling into Mesenchyme, receptor druggability
scripts/12_receptor_correlation_localization.R  # CD44/ITGA1/ITGB1/EDNRB vs fibrogenic genes; ITGA1 vs ITGA8 vs ITGA11
scripts/13_receptor_ligand_tracing.R   # which cell type sends each receptor's ligand, healthy vs cirrhotic
scripts/14_receptor_clinical_context.R # literature-verified compound/trial status per receptor
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
