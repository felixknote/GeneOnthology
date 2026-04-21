# =============================================================================
# GO + KEGG Enrichment — E. coli K12
#
# Methods follow:
#   [1] Wu et al. (2021) clusterProfiler 4.0: A universal enrichment tool for
#       interpreting omics data. The Innovation 2(3):100141.
#       doi:10.1016/j.xinn.2021.100141
#   [2] Yu et al. (2015) GOSemSim: an R package for measuring semantic
#       similarity among GO terms and gene products. Bioinformatics 26(7):976.
#       doi:10.1093/bioinformatics/btq064
#   [3] Yu et al. (2015) DOSE: an R/Bioconductor package for disease ontology
#       semantic and enrichment analysis. Bioinformatics 31(4):608.
#       doi:10.1093/bioinformatics/btu684
#   [4] Tibshirani et al. (2001) Estimating the number of clusters in a data
#       set via the gap statistic. JRSS-B 63(2):411.
#       doi:10.1111/1467-9868.00293
#   [5] Kanehisa & Goto (2000) KEGG: Kyoto Encyclopedia of Genes and Genomes.
#       Nucleic Acids Res 28(1):27. doi:10.1093/nar/28.1.27
#
# Statistical thresholds: BH-adjusted p < 0.05, semantic similarity >= 0.7
# Cluster number k: data-driven via gap statistic [4] — no manual tuning.
#
# Outputs (all 16:9):
#   go_treeplot.pdf / .png    — hierarchical dendrogram of GO:BP terms
#   go_cnetplot.pdf / .png    — gene-concept network (gene → GO term)
#   go_single_class.pdf / .png — gene dendrogram, leaves = genes
#   kegg_dotplot.pdf / .png   — KEGG pathway enrichment dot plot
#   kegg_cnetplot.pdf / .png  — KEGG gene-concept network
#
# Install once: run install_packages.R
# =============================================================================

# ---- make user library visible in VSCode R sessions -------------------------
local({
  user_lib <- file.path(Sys.getenv("APPDATA"), "R", "win-library",
                        paste0(R.version$major, ".",
                               substr(R.version$minor, 1, 1)))
  if (dir.exists(user_lib) && !(user_lib %in% .libPaths()))
    .libPaths(c(user_lib, .libPaths()))
})

# -----------------------------------------------------------------------------
# CONFIGURATION — edit only this block
# -----------------------------------------------------------------------------
gene_names <- c(
  "ftsZ", "ftsI", "mrdA", "mrcB", "mrcA", "murA", "murC",
  "rpsA", "rpsL", "rplA", "rplC",
  "folA", "folP",
  "dnaE", "dnaB", "gyrA", "gyrB", "parC", "parE",
  "rpoA", "rpoB",
  "secA", "secY",
  "lpxA", "lpxC", "lptA", "lptC", "msbA"
)

# Output directory — NULL saves next to this script
output_dir <- NULL
# -----------------------------------------------------------------------------


# ---- Setup ------------------------------------------------------------------
suppressPackageStartupMessages({
  library(org.EcK12.eg.db)
  library(clusterProfiler)
  library(enrichplot)
  library(GOSemSim)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

if (is.null(output_dir)) {
  # Works when run via Rscript directly
  args        <- commandArgs(trailingOnly = FALSE)
  script_arg  <- grep("^--file=", args, value = TRUE)
  output_dir  <- if (length(script_arg) == 1) {
    dirname(normalizePath(sub("^--file=", "", script_arg)))
  } else {
    # Fallback: sourced interactively
    tryCatch(dirname(normalizePath(sys.frame(1)$ofile)),
             error = function(e) getwd())
  }
}
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)


# ---- 1. Symbol → Entrez ID --------------------------------------------------
message("Converting ", length(gene_names), " gene symbols to Entrez IDs ...")

gene_df <- tryCatch(
  clusterProfiler::bitr(gene_names,
                        fromType = "SYMBOL", toType = "ENTREZID",
                        OrgDb    = org.EcK12.eg.db),
  error = function(e) stop("bitr() failed: ", conditionMessage(e))
)

missed <- setdiff(gene_names, gene_df$SYMBOL)
if (length(missed) > 0)
  warning(length(missed), " gene(s) not found in org.EcK12.eg.db and skipped: ",
          paste(missed, collapse = ", "))
if (nrow(gene_df) == 0) stop("No genes mapped. Aborting.")

message("  Mapped: ", nrow(gene_df), " / ", length(gene_names))


# ---- 2. GO:BP enrichment (publication-standard thresholds) ------------------
# Thresholds follow Wu et al. 2021 (clusterProfiler 4.0, The Innovation)
# and standard Bioinformatics practice:
#   p.adj (BH) < 0.05,  q-value < 0.2
message("Running GO:BP enrichment (BH-adjusted p < 0.05) ...")

ego <- enrichGO(
  gene          = gene_df$ENTREZID,
  universe      = keys(org.EcK12.eg.db, keytype = "ENTREZID"),
  OrgDb         = org.EcK12.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE          # gene symbols instead of Entrez IDs in plots
)

n_raw <- sum(ego@result$p.adjust < 0.05)
message("  Significant terms before simplification: ", n_raw)

if (n_raw == 0) stop("No significant GO:BP terms found. Try a larger gene set.")


# ---- 3. Remove redundant GO terms (semantic similarity ≥ 0.7) ---------------
# simplify() keeps the most significant representative from each cluster of
# semantically similar terms. Cutoff 0.7 is the published default.
# Reference: Yu et al. (2015) OMICS doi:10.1089/omi.2014.0083
message("Removing redundant GO terms (similarity cutoff 0.7) ...")

ego_s <- clusterProfiler::simplify(ego,
                                   cutoff   = 0.7,
                                   by       = "p.adjust",
                                   select_fun = min)

n_simplified <- nrow(ego_s@result)
message("  Terms after simplification: ", n_simplified)


# ---- 4. Pairwise semantic similarity (for treeplot dendrogram) --------------
message("Computing pairwise Jaccard similarity ...")
ego_sim <- pairwise_termsim(ego_s, method = "JC")


# ---- 5. Data-driven cluster number via gap statistic on GO dendrogram -------
# Build the same hclust object that treeplot() uses internally, then use the
# gap statistic to find the natural number of clusters (no user input needed).
.auto_k <- function(sim_mat, k_max = 10) {
  # sim_mat: pairwise similarity matrix from pairwise_termsim()
  mat <- as.matrix(sim_mat@termsim)
  mat[is.na(mat)] <- 0
  dist_mat <- as.dist(1 - mat)
  if (nrow(mat) < 3) return(1L)

  hc <- hclust(dist_mat, method = "ward.D")
  k_max <- min(k_max, nrow(mat) - 1L)

  # Gap statistic (Tibshirani et al. 2001)
  set.seed(42)
  gap <- cluster::clusGap(as.matrix(dist_mat),
                          FUN       = function(x, k) list(cluster = cutree(hc, k)),
                          K.max     = k_max,
                          B         = 50,
                          d.power   = 2)
  k_opt <- cluster::maxSE(gap$Tab[, "gap"],
                           gap$Tab[, "SE.sim"],
                           method = "Tibs2001SEmax")
  as.integer(k_opt)
}

# cluster package ships with base R — no extra install needed
if (!requireNamespace("cluster", quietly = TRUE))
  stop("'cluster' package is missing — it normally ships with R base.")

k <- .auto_k(ego_sim, k_max = min(10L, n_simplified - 1L))
message("  Optimal cluster count (gap statistic): ", k)


# ---- 6. Treeplot ------------------------------------------------------------
message("Building hierarchical treeplot ...")

p_tree <- treeplot(
  ego_sim,
  nCluster      = k,
  showCategory  = n_simplified,   # show all non-redundant terms
  color         = "p.adjust",
  hexpand       = 0.85,           # right-side space for labels + bar + group names
  tiplab_offset = 0.3,            # gap between tree tip and label text
  extend        = 3.5,            # how far the highlight extends past tip labels → vertical bar position
  cladelab_offset = 9.0           # group name offset from vertical bar
) +
  labs(
    title    = "GO Biological Process — hierarchical term tree",
    subtitle = paste0(
      nrow(gene_df), " E. coli K12 genes  |  ",
      n_simplified, " non-redundant terms (from ", n_raw, " significant)  |  ",
      "BH p.adj < 0.05  |  k = ", k, " clusters (gap statistic)"
    )
  ) +
  theme(
    plot.title    = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 11, color = "grey40"),
    legend.position = "none"
  )

# Post-hoc: shrink GO term label text so dense clusters don't overlap.
# Also nudge labels right so they don't sit on top of the dot nodes.
for (.i in seq_along(p_tree$layers)) {
  lyr <- p_tree$layers[[.i]]
  if (inherits(lyr$geom, "GeomText")) {
    lyr$aes_params$size    <- 4.2      # larger, legible labels
    lyr$aes_params$hjust   <- 0        # left-align so text reads away from dots
  }
}
rm(.i, lyr)


# ---- 7. Save + display ------------------------------------------------------
# print() sends the plot to the active R graphics device (RStudio / VSCode viewer)
# ggsave() writes to file — it does NOT display; both are needed.
save_plot <- function(p, base, w = 16, h = 11) {
  print(p)                                          # show in R device
  for (ext in c("pdf", "png")) {
    path <- file.path(output_dir, paste0(base, ".", ext))
    ggplot2::ggsave(path, plot = p, width = w, height = h,
                    dpi = if (ext == "png") 300 else 72)
    message("  Saved: ", path)
  }
}

# Height driven by term count so content fills the canvas (no empty bottom).
# 0.38" per term works well at font size 4.2; min 10" for titles/legend.
tree_h <- max(10, n_simplified * 0.38)
tree_w <- round(tree_h * 16 / 9, 1)
save_plot(p_tree, "go_treeplot", w = tree_w, h = tree_h)             # 16:9


# ---- 8. Gene-concept network (shows which gene drives which GO term) ---------
# cnetplot() is the standard companion to treeplot: it draws GO term nodes
# connected to their member gene nodes so you can read off which genes are
# responsible for each enriched term.  Genes shared across multiple terms
# appear as hubs — common in multi-functional proteins.
message("Building gene-concept network (cnetplot) ...")

set.seed(42)  # reproducible force-directed layout
p_cnet <- cnetplot(
  ego_s,
  showCategory = n_simplified    # all non-redundant terms
) +
  labs(
    title    = "Gene-concept network — GO Biological Process",
    subtitle = paste0(
      nrow(gene_df), " E. coli K12 genes  |  ",
      n_simplified, " non-redundant terms  |  BH p.adj < 0.05"
    )
  ) +
  theme(
    plot.title    = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 11, color = "grey40")
  )

# Increase node-label text size post-hoc.
# cnetplot uses GeomTextRepel (ggrepel) — check by class name, not inheritance.
for (.i in seq_along(p_cnet$layers)) {
  .cls <- class(p_cnet$layers[[.i]]$geom)[1]
  if (.cls %in% c("GeomTextRepel", "GeomText", "GeomLabel", "GeomLabelRepel"))
    p_cnet$layers[[.i]]$aes_params$size <- 6
}
rm(.i, .cls)

save_plot(p_cnet, "go_cnetplot", w = 32, h = 18)   # 16:9


# ---- 9. Gene dendrogram — leaves = genes, branches = GO pathways -----------
# Design follows Eisen et al. (1998) hierarchical clustering dendrogram
# (PNAS 95:14863) adapted for GO-term space.
#
# Distance matrix (Prathipati et al. 2016, BMC Bioinformatics approach):
#   Within-cluster gene pair  → Jaccard distance on shared GO terms ∈ [0, 1]
#   Between-cluster gene pair → 1 + inter-cluster semantic distance ∈ [1, 2]
# The +1 offset guarantees every within-cluster distance is strictly less
# than every between-cluster distance, so the dendrogram topology mirrors
# the GO pathway structure exactly — genes sharing a pathway always form a
# clade together before merging with genes from other pathways.
#
# Single-class assignment (dominant-cluster, Alexa & Rahnenfuhrer 2009):
#   For each gene, count enriched GO terms per treeplot cluster.
#   Assign to the cluster with the highest count; tiebreak = lowest p.adj.
message("Building gene dendrogram (single-class assignment) ...")

# Install ggdendro if missing (CRAN, ships as ggplot2 companion)
if (!requireNamespace("ggdendro", quietly = TRUE))
  install.packages("ggdendro", repos = "https://cloud.r-project.org",
                   lib = .libPaths()[1], quiet = TRUE)
library(ggdendro)

# -- 9a. Reproduce treeplot's internal clustering (term → cluster ID) ---------
# enrichplot >= 1.19 stores termsim with Description as rownames, not GO IDs.
# Handle both formats by checking whether rownames look like GO IDs.
mat_sim <- as.matrix(ego_sim@termsim)
mat_sim[is.na(mat_sim)] <- 0
term_cut <- cutree(hclust(as.dist(1 - mat_sim), method = "ward.D"), k = k)

rn <- rownames(mat_sim)
if (all(grepl("^GO:", rn))) {
  term_to_cluster <- setNames(term_cut, rn)           # keyed by GO ID
} else {
  # keyed by Description — remap to GO ID via ego_s@result
  desc2id <- setNames(ego_s@result$ID, ego_s@result$Description)
  term_to_cluster <- setNames(term_cut, desc2id[rn])  # keyed by GO ID
}

# Representative cluster label = lowest p.adj GO term in each cluster
cluster_labels <- ego_s@result %>%
  mutate(cluster = term_to_cluster[ID]) %>%
  group_by(cluster) %>%
  slice_min(p.adjust, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
             transmute(cluster,
                    cluster_label = Description)

# -- 9b. Dominant-cluster assignment per gene ---------------------------------
gene_long <- ego_s@result %>%
  mutate(cluster = term_to_cluster[ID]) %>%
  select(cluster, geneID, p.adjust) %>%
  separate_rows(geneID, sep = "/")

gene_class <- gene_long %>%
  filter(!is.na(cluster)) %>%
  group_by(geneID, cluster) %>%
  summarise(n_terms = n(), best_padj = min(p.adjust, na.rm = TRUE), .groups = "drop") %>%
  group_by(geneID) %>%
  arrange(desc(n_terms), best_padj) %>%
  slice(1) %>%
  ungroup() %>%
  left_join(cluster_labels, by = "cluster") %>%
  filter(!is.na(cluster), !is.na(geneID))

# -- 9c. Gene × GO-term binary membership matrix ------------------------------
gene_go_pairs <- ego_s@result %>%
  select(ID, geneID) %>%
  separate_rows(geneID, sep = "/")

g_vec <- gene_class$geneID          # genes present in enrichment
t_vec <- unique(gene_go_pairs$ID)

gene_term_mat <- matrix(0L, nrow = length(g_vec), ncol = length(t_vec),
                        dimnames = list(g_vec, t_vec))
for (i in seq_len(nrow(gene_go_pairs)))
  if (gene_go_pairs$geneID[i] %in% g_vec)
    gene_term_mat[gene_go_pairs$geneID[i], gene_go_pairs$ID[i]] <- 1L

# -- 9d. Inter-cluster semantic distance (centroid average) -------------------
# termsim may be indexed by Description rather than GO ID — create a
# GO-ID → matrix-key lookup so we can always index the matrix correctly.
id_to_matkey <- if (all(grepl("^GO:", rn))) {
  setNames(ego_s@result$ID, ego_s@result$ID)          # identity
} else {
  setNames(ego_s@result$Description, ego_s@result$ID) # ID → Description
}

cl_ids  <- sort(unique(gene_class$cluster))
cl_dist <- matrix(0, length(cl_ids), length(cl_ids),
                  dimnames = list(cl_ids, cl_ids))
for (i in seq_along(cl_ids))
  for (j in seq_along(cl_ids)) {
    if (i == j) next
    ti <- id_to_matkey[names(term_to_cluster)[term_to_cluster == cl_ids[i]]]
    tj <- id_to_matkey[names(term_to_cluster)[term_to_cluster == cl_ids[j]]]
    ti <- ti[!is.na(ti)]; tj <- tj[!is.na(tj)]
    cl_dist[i, j] <- 1 - mean(ego_sim@termsim[ti, tj, drop = FALSE], na.rm = TRUE)
  }

# -- 9e. Block gene-gene distance matrix --------------------------------------
gene_cv  <- setNames(gene_class$cluster, gene_class$geneID)
n_g      <- length(g_vec)
D        <- matrix(0, n_g, n_g, dimnames = list(g_vec, g_vec))

for (i in seq_len(n_g))
  for (j in seq_len(i - 1L)) {
    gi <- g_vec[i]; gj <- g_vec[j]
    ci <- gene_cv[gi]; cj <- gene_cv[gj]
    vi <- gene_term_mat[gi, ]; vj <- gene_term_mat[gj, ]
    u  <- sum(vi | vj)
    jac <- if (u == 0) 0 else 1 - sum(vi & vj) / u
    d   <- if (!is.na(ci) && !is.na(cj) && ci == cj) jac
           else if (!is.na(ci) && !is.na(cj)) 1 + cl_dist[as.character(ci), as.character(cj)]
           else 1
    D[i, j] <- D[j, i] <- d
  }

gene_hc   <- hclust(as.dist(D), method = "ward.D")
dend_data <- ggdendro::dendro_data(as.dendrogram(gene_hc), type = "rectangle")

# -- 9f. Annotate leaves with cluster info ------------------------------------
leaf_df <- dend_data$labels %>%
  left_join(gene_class %>% select(geneID, cluster, cluster_label),
            by = c("label" = "geneID"))

# In the HORIZONTAL layout (left-to-right), the ggdendro coordinate roles swap:
#   plot-x  ← segment y  (Ward distance; reversed so root is at LEFT)
#   plot-y  ← segment x  (leaf positions along vertical axis)
# Cluster span is now along the VERTICAL axis (original x)
bracket_df <- leaf_df %>%
  group_by(cluster, cluster_label) %>%
  summarise(y_min = min(x), y_max = max(x), .goals = "drop", .groups = "drop") %>%
  mutate(y_mid = (y_min + y_max) / 2)

max_h <- max(dend_data$segments$y, dend_data$segments$yend)

# Cluster bar: vertical line to the RIGHT of the leaves (x = negative in reversed scale)
bar_x   <- -0.10 * max_h   # colored bar just right of leaf labels
label_x <- -0.22 * max_h   # pathway name text further right

# Color palette — one color per cluster
n_cl   <- length(cl_ids)
cl_pal <- setNames(
  if (n_cl <= 8) RColorBrewer::brewer.pal(max(3, n_cl), "Dark2")[seq_len(n_cl)]
  else           grDevices::hcl.colors(n_cl, palette = "Dark 3"),
  cl_ids)

# Right-side x-limit: enough room for the longest cluster label text.
# At size 4.5 (ggplot mm units ≈ 1.6 mm/pt), a 50-char label ≈ 55 mm.
# 1 data unit ≈ canvas_width_mm / x_range; at 32" × 25.4 = 812 mm and
# x_range ≈ max_h*1.05 + |label_x|, we budget label_x - max_label_chars*0.09.
max_label_chars <- max(nchar(bracket_df$cluster_label), na.rm = TRUE)
x_right_lim     <- label_x - max_label_chars * 0.09 * (max_h / 10)

# -- 9g. Build horizontal dendrogram plot (left = root, right = leaves) -------
p_gene_tree <- ggplot() +
  # Tree branches (swap x↔y so height is horizontal, leaf order is vertical)
  geom_segment(data = dend_data$segments,
               aes(x = y, y = x, xend = yend, yend = xend),
               color = "grey35", linewidth = 0.5) +
  # Gene names — just right of each leaf (anchor at x=0, text extends right)
  geom_text(data = leaf_df,
            aes(x = -0.02 * max_h, y = x, label = label,
                color = factor(cluster)),
            angle = 0, hjust = 0, vjust = 0.5, size = 5) +
  # Colored vertical bar spanning each cluster's genes
  geom_segment(data = bracket_df,
               aes(x = bar_x, xend = bar_x,
                   y = y_min - 0.35, yend = y_max + 0.35,
                   color = factor(cluster)),
               linewidth = 5, lineend = "round") +
  # Cluster pathway label — horizontal text to the right of the bar
  geom_text(data = bracket_df,
            aes(x = label_x, y = y_mid, label = cluster_label,
                color = factor(cluster)),
            angle = 0, hjust = 0, vjust = 0.5, size = 4.5, fontface = "italic") +
  scale_color_manual(values = cl_pal, guide = "none") +
  # Reversed x: root (max_h) at LEFT, leaves (0) at RIGHT; extend for labels
  scale_x_reverse(
    limits = c(max_h * 1.06, x_right_lim),
    name   = "Ward distance"
  ) +
  scale_y_continuous(expand = expansion(add = 0.8)) +
  labs(
    title    = "Gene dendrogram — single pathway assignment per gene",
    subtitle = paste0(
      "Leaves = genes  |  Branch length = shared GO term distance (Jaccard)  |  ",
      "Cluster separation enforced by inter-cluster semantic distance  |  k = ",
      k, " GO clusters"
    ),
    y = NULL
  ) +
  theme_classic(base_size = 16) +
  theme(
    plot.title    = element_text(size = 20, face = "bold"),
    plot.subtitle = element_text(size = 13, color = "grey40"),
    axis.title.x  = element_text(size = 14),
    axis.text.x   = element_text(size = 13),
    axis.text.y   = element_blank(),
    axis.ticks.y  = element_blank(),
    axis.line.y   = element_blank()
  )

save_plot(p_gene_tree, "go_single_class", w = 32, h = 18)   # 16:9


# ---- 10. KEGG pathway enrichment --------------------------------------------
# enrichKEGG() uses KEGG REST API (requires internet). Organism code "eco"
# = E. coli K12 MG1655 (KEGG organism T00007).
# KEGG expects its own locus IDs (b-numbers), not NCBI Entrez IDs.
# bitr_kegg() converts ncbi-geneid → kegg via the KEGG REST API.
# Reference: Kanehisa & Goto 2000 [5]
message("Running KEGG pathway enrichment ...")

kegg_map <- tryCatch(
  clusterProfiler::bitr_kegg(gene_df$ENTREZID,
                              fromType = "ncbi-geneid",
                              toType   = "kegg",
                              organism = "eco"),
  error = function(e) { message("  bitr_kegg() failed: ", conditionMessage(e)); NULL }
)

if (is.null(kegg_map) || nrow(kegg_map) == 0) {
  message("  Could not convert Entrez IDs to KEGG locus IDs — skipping KEGG.")
  ekegg <- NULL
} else {
  message("  Mapped ", nrow(kegg_map), " / ", nrow(gene_df), " genes to KEGG IDs")
  ekegg <- tryCatch(
    clusterProfiler::enrichKEGG(
      gene          = kegg_map$kegg,
      organism      = "eco",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.2
    ),
    error = function(e) {
      message("  KEGG enrichment failed: ", conditionMessage(e))
      NULL
    }
  )
}

if (is.null(ekegg) || nrow(ekegg@result[ekegg@result$p.adjust < 0.05, ]) == 0) {
  message("  No significant KEGG pathways found.")
} else {
  n_kegg  <- nrow(ekegg@result[ekegg@result$p.adjust < 0.05, ])
  message("  Significant KEGG pathways (", n_kegg, "):")
  sig <- ekegg@result[ekegg@result$p.adjust < 0.05,
                       c("Description", "GeneRatio", "p.adjust", "geneID")]
  for (i in seq_len(nrow(sig)))
    message("    ", sig$Description[i],
            "  GeneRatio=", sig$GeneRatio[i],
            "  p.adj=", signif(sig$p.adjust[i], 3))
}

message("Done.")
