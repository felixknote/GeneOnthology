# GeneOntology — GO & KEGG Enrichment for *E. coli* K-12

R pipeline for Gene Ontology (Biological Process) and KEGG pathway
enrichment of a user-supplied gene list, with publication-ready
hierarchical treeplots, gene–concept networks, and a single-class
gene dendrogram. No biological information is hardcoded — only gene
symbols.

## What it does

Given a list of *E. coli* K-12 MG1655 gene symbols in `go_tree.R`, the
script:

1. **Symbol → Entrez** via `org.EcK12.eg.db` (`clusterProfiler::bitr`)
2. **GO:BP enrichment** with BH-adjusted p < 0.05, q < 0.2 (`enrichGO`)
3. **Redundancy reduction** by semantic similarity ≥ 0.7 (`simplify`)
4. **Pairwise term similarity** (Jaccard; `pairwise_termsim`)
5. **Data-driven cluster number k** via the gap statistic
   (Tibshirani et al. 2001) — no manual tuning
6. **Hierarchical treeplot** of GO:BP terms (`treeplot`)
7. **Gene–concept network** linking genes to their GO terms (`cnetplot`)
8. **Single-class gene dendrogram** — each gene assigned to its
   dominant GO cluster; block distance matrix guarantees the tree
   topology mirrors GO pathway structure
9. **KEGG pathway enrichment** (`enrichKEGG`, organism `eco`)

### Outputs

| File | Contents |
|---|---|
| `go_treeplot.{pdf,png}` | Hierarchical dendrogram of non-redundant GO:BP terms |
| `go_cnetplot.{pdf,png}` | Gene–concept network (gene ↔ GO term) |
| `go_single_class.{pdf,png}` | Gene dendrogram, one pathway label per gene |
| *(console)* | Significant KEGG pathways with GeneRatio and adjusted p |

All plots are 16:9; sizes scale with the number of terms so dense
figures don't overlap.

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
Rscript go_tree.R
```

or source `go_tree.R` from RStudio / VSCode. Plots render to the
active graphics device and are saved alongside the script.

## Customising the gene list

Edit the `gene_names` vector at the top of `go_tree.R`:

```r
gene_names <- c(
  "ftsZ", "gyrA", "rpoB",   # add or remove gene symbols
  ...
)
```

Everything downstream — enrichment, clustering, pathway labels,
colors, plot dimensions — is resolved automatically at runtime.

## Methods

| Step | Reference |
|---|---|
| `clusterProfiler` 4.0 | Wu et al. (2021) *The Innovation* 2:100141 |
| GO semantic similarity | Yu et al. (2015) *Bioinformatics* 26:976 |
| DOSE framework | Yu et al. (2015) *Bioinformatics* 31:608 |
| Gap statistic (k) | Tibshirani et al. (2001) *JRSS-B* 63:411 |
| KEGG pathways | Kanehisa & Goto (2000) *Nucleic Acids Res* 28:27 |

Statistical thresholds follow Bioinformatics community standards:
BH-adjusted p < 0.05, semantic similarity ≥ 0.7, q-value < 0.2.

## Organism

*Escherichia coli* K-12 MG1655 — NCBI taxonomy ID 83333, KEGG prefix
`eco`, Bioconductor annotation package `org.EcK12.eg.db`.
