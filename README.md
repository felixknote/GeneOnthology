# GeneOntology — Multi-Ontology GO Enrichment for *E. coli* K-12

R pipeline for Gene Ontology enrichment across all three namespaces
(Biological Process, Molecular Function, Cellular Component) for a
user-supplied gene list. Produces publication-ready hierarchical
treeplots, gene–concept networks, single-class gene dendrograms, and
gene × gene distance matrices. No biological information is hardcoded —
only gene symbols.

## Scripts

| Script | Purpose |
|---|---|
| `go_analysis.R` | **Main script** — runs BP, MF, CC separately + combined analysis |
| `go_tree.R` | Legacy single-ontology script (BP only) |

---

## `go_analysis.R` — Multi-Ontology Pipeline

### What it does

1. **Symbol → Entrez** via `org.EcK12.eg.db` (`clusterProfiler::bitr`) — once, shared
2. **Per-ontology loop** over BP, MF, CC — each written to its own subfolder:
   - GO enrichment with BH-adjusted p < 0.05, q < 0.2 (`enrichGO`)
   - Redundancy reduction by semantic similarity ≥ 0.7 (`simplify`)
   - Pairwise term similarity via Jaccard (`pairwise_termsim`)
   - Data-driven cluster number k via the gap statistic (Tibshirani et al. 2001)
   - Hierarchical treeplot of non-redundant GO terms (`treeplot`)
   - Gene–concept network (`cnetplot`)
   - Single-class gene dendrogram (block distance: Jaccard within cluster,
     1 + semantic distance between clusters)
   - Gene × gene distance matrix saved as CSV
3. **Combined analysis** in `combined/` — merges all three enrichment sets:
   - Gene × gene Jaccard distance matrix across BP + MF + CC terms
   - Gene dendrogram colored by dominant ontology
   - Per-ontology cnetplots in the combined context

### Output folder structure

```
Gene Ontology/
├── BP/
│   ├── go_treeplot.{pdf,png}
│   ├── go_cnetplot.{pdf,png}
│   ├── go_single_class.{pdf,png}
│   └── gene_distance_matrix.csv
├── MF/
│   └── … (same files)
├── CC/
│   └── … (same files)
└── combined/
    ├── go_cnetplot_BP.{pdf,png}
    ├── go_cnetplot_MF.{pdf,png}
    ├── go_cnetplot_CC.{pdf,png}
    ├── go_single_class.{pdf,png}
    └── gene_distance_matrix.csv   ← Jaccard across all GO terms
```

All plots are 16:9; treeplot height scales with term count.

---

## Setup

R ≥ 4.3 with Bioconductor. Install once:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "org.EcK12.eg.db",
  "clusterProfiler",
  "enrichplot",
  "GOSemSim"
))

install.packages(c(
  "ggplot2", "dplyr", "tidyr",
  "ggdendro", "RColorBrewer", "cluster"
))
```

## Usage

```bash
Rscript go_analysis.R
```

or source `go_analysis.R` from RStudio / VSCode. Subfolders are
created automatically; plots render to the active graphics device and
are also saved to file.

## Customising the gene list

Edit the `gene_names` vector at the top of `go_analysis.R`:

```r
gene_names <- c(
  "ftsZ", "gyrA", "rpoB",   # add or remove gene symbols
  ...
)
```

Everything downstream — enrichment, clustering, pathway labels,
colors, plot dimensions — is resolved automatically at runtime.

## Distance matrix

Each subfolder contains `gene_distance_matrix.csv` — a symmetric
gene × gene matrix where 0 = identical GO term profile and 1 = no
shared terms. Within the per-ontology runs the distance uses a block
structure (Jaccard within cluster; 1 + inter-cluster semantic distance
between clusters). The combined matrix uses pure Jaccard across the
union of all GO terms from BP, MF, and CC.

## Methods

| Step | Reference |
|---|---|
| `clusterProfiler` 4.0 | Wu et al. (2021) *The Innovation* 2:100141 |
| GO semantic similarity | Yu et al. (2015) *Bioinformatics* 26:976 |
| DOSE framework | Yu et al. (2015) *Bioinformatics* 31:608 |
| Gap statistic (k) | Tibshirani et al. (2001) *JRSS-B* 63:411 |

Statistical thresholds: BH-adjusted p < 0.05, semantic similarity ≥ 0.7,
q-value < 0.2.

## Organism

*Escherichia coli* K-12 MG1655 — NCBI taxonomy ID 83333,
Bioconductor annotation package `org.EcK12.eg.db`.
