# GeneOntology — E. coli Gene Classifier

Fetches functional annotations for a list of *E. coli* K-12 MG1655 genes from
three public databases and produces a hierarchy tree and an annotated Excel
workbook. No biological information is hardcoded — only gene names.

## What it does

Given a flat list of gene names in `classify_genes.py`, the script:

1. **UniProt** — resolves accession, protein name, GO terms, and KEGG cross-reference
2. **EcoCyc/BioCyc** — resolves protein/gene frame IDs, fetches curated GO annotations,
   pathway memberships, pathway hierarchy, and molecular weight
3. **KEGG** — fetches specific pathway associations

Genes are grouped by the most informative level of the EcoCyc super-pathway
hierarchy. When no pathway data is available the most-annotated GO Biological
Process term (majority vote across EcoCyc + UniProt sources) is used as fallback.

### Outputs

| File | Contents |
|---|---|
| `gene_hierarchy.txt` | ASCII tree grouped by biological function |
| `gene_classification.xlsx` | One row per gene — accession, protein name, group, EcoCyc pathways, KEGG pathways, GO terms (BP / MF / CC), molecular weight |

## Setup

```bash
pip install requests openpyxl
```

Create a `.env` file next to the script with your
[EcoCyc](https://ecocyc.org) account credentials
(free registration required for authenticated API access):

```
ECOCYC_EMAIL=you@example.com
ECOCYC_PASSWORD=yourpassword
```

> `.env` is listed in `.gitignore` and will never be committed.

## Usage

```bash
python classify_genes.py
```

The script prints per-gene progress as it runs (~5–10 min for 28 genes due to
API rate limits). Outputs are written to the same directory as the script.

## Customising the gene list

Edit `GENE_NAMES` at the top of `classify_genes.py`. Add or remove gene names;
everything else — accessions, groups, colors, pathway hierarchy — is resolved
automatically at runtime.

```python
GENE_NAMES: list[str] = [
    "ftsZ", "gyrA", "rpoB",   # add or remove genes here
    ...
]
```

## Data sources

| Source | URL | Auth |
|---|---|---|
| EcoCyc / BioCyc | https://websvc.biocyc.org | Free account (Tier 1) |
| UniProt REST | https://rest.uniprot.org | None |
| KEGG REST | https://rest.kegg.jp | None |

## Organism

*Escherichia coli* K-12 MG1655 — NCBI taxonomy ID 83333, KEGG prefix `eco`,
BioCyc database ID `ECOLI`.
