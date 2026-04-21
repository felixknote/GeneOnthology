# =============================================================================
# GO Enrichment — E. coli K12 — Multi-Ontology
#
# Runs the full enrichment pipeline for each GO namespace separately
# (BP, MF, CC) and then a combined analysis that merges all three.
#
# Methods:
#   Wu et al. (2021) clusterProfiler 4.0. The Innovation 2:100141.
#   Yu et al. (2015) GOSemSim. Bioinformatics 26:976.
#   Tibshirani et al. (2001) Gap statistic. JRSS-B 63:411.
#
# Outputs per ontology subfolder (BP/, MF/, CC/):
#   go_treeplot.{pdf,png}      hierarchical GO term dendrogram
#   go_cnetplot.{pdf,png}      gene-concept network
#   go_single_class.{pdf,png}  gene dendrogram (one pathway label per gene)
#   gene_distance_matrix.csv   gene x gene Jaccard distance matrix
#
# Outputs in combined/:
#   go_cnetplot_BP/MF/CC.{pdf,png}  per-ontology cnetplots (shared gene set)
#   go_single_class.{pdf,png}       gene dendrogram colored by dominant ontology
#   gene_distance_matrix.csv        gene x gene Jaccard matrix across BP+MF+CC
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

# Root output directory — NULL saves next to this script
output_base <- NULL
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
  library(RColorBrewer)
})

if (!requireNamespace("ggdendro", quietly = TRUE))
  install.packages("ggdendro", repos = "https://cloud.r-project.org",
                   lib = .libPaths()[1], quiet = TRUE)
library(ggdendro)

if (!requireNamespace("cluster", quietly = TRUE))
  stop("'cluster' package missing — it normally ships with R base.")

if (is.null(output_base)) {
  args       <- commandArgs(trailingOnly = FALSE)
  script_arg <- grep("^--file=", args, value = TRUE)
  output_base <- if (length(script_arg) == 1) {
    dirname(normalizePath(sub("^--file=", "", script_arg)))
  } else {
    tryCatch(dirname(normalizePath(sys.frame(1)$ofile)),
             error = function(e) getwd())
  }
}


# ---- Shared helpers ---------------------------------------------------------

ont_label <- function(ont) {
  lut <- c(BP = "GO Biological Process",
           MF = "GO Molecular Function",
           CC = "GO Cellular Component")
  ifelse(ont %in% names(lut), lut[ont], ont)
}

save_plot <- function(p, base, output_dir, w = 16, h = 11) {
  print(p)
  for (ext in c("pdf", "png")) {
    path <- file.path(output_dir, paste0(base, ".", ext))
    ggplot2::ggsave(path, plot = p, width = w, height = h,
                    dpi = if (ext == "png") 300 else 72)
    message("  Saved: ", path)
  }
}

.auto_k <- function(sim_mat, k_max = 10) {
  mat <- as.matrix(sim_mat@termsim)
  mat[is.na(mat)] <- 0
  dist_mat <- as.dist(1 - mat)
  if (nrow(mat) < 3) return(1L)
  hc    <- hclust(dist_mat, method = "ward.D")
  k_max <- min(k_max, nrow(mat) - 1L)
  set.seed(42)
  gap <- cluster::clusGap(as.matrix(dist_mat),
                          FUN   = function(x, k) list(cluster = cutree(hc, k)),
                          K.max = k_max, B = 50, d.power = 2)
  as.integer(cluster::maxSE(gap$Tab[, "gap"], gap$Tab[, "SE.sim"],
                             method = "Tibs2001SEmax"))
}

# Builds the block gene-gene distance matrix (Prathipati et al. 2016 approach).
# Same within a cluster: Jaccard on GO term membership [0,1].
# Different clusters:   1 + inter-cluster semantic distance [1,2].
# Returns the numeric matrix D with gene names as dimnames.
.build_gene_dist <- function(ego_s, ego_sim, k) {

  # -- term → cluster mapping --------------------------------------------------
  mat_sim  <- as.matrix(ego_sim@termsim)
  mat_sim[is.na(mat_sim)] <- 0
  term_cut <- cutree(hclust(as.dist(1 - mat_sim), method = "ward.D"), k = k)
  rn       <- rownames(mat_sim)

  if (all(grepl("^GO:", rn))) {
    term_to_cluster <- setNames(term_cut, rn)
  } else {
    desc2id         <- setNames(ego_s@result$ID, ego_s@result$Description)
    term_to_cluster <- setNames(term_cut, desc2id[rn])
  }

  # -- representative cluster label (lowest p.adj term) -----------------------
  cluster_labels <- ego_s@result %>%
    mutate(cluster = term_to_cluster[ID]) %>%
    group_by(cluster) %>%
    slice_min(p.adjust, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    transmute(cluster, cluster_label = Description)

  # -- gene → dominant cluster ------------------------------------------------
  gene_long <- ego_s@result %>%
    mutate(cluster = term_to_cluster[ID]) %>%
    select(cluster, geneID, p.adjust) %>%
    separate_rows(geneID, sep = "/")

  gene_class <- gene_long %>%
    filter(!is.na(cluster)) %>%
    group_by(geneID, cluster) %>%
    summarise(n_terms = n(), best_padj = min(p.adjust, na.rm = TRUE),
              .groups = "drop") %>%
    group_by(geneID) %>%
    arrange(desc(n_terms), best_padj) %>%
    slice(1) %>%
    ungroup() %>%
    left_join(cluster_labels, by = "cluster") %>%
    filter(!is.na(cluster), !is.na(geneID))

  # -- gene x term binary matrix ----------------------------------------------
  gene_go_pairs <- ego_s@result %>%
    select(ID, geneID) %>%
    separate_rows(geneID, sep = "/")

  g_vec <- gene_class$geneID
  t_vec <- unique(gene_go_pairs$ID)

  gene_term_mat <- matrix(0L, nrow = length(g_vec), ncol = length(t_vec),
                          dimnames = list(g_vec, t_vec))
  for (i in seq_len(nrow(gene_go_pairs)))
    if (gene_go_pairs$geneID[i] %in% g_vec)
      gene_term_mat[gene_go_pairs$geneID[i], gene_go_pairs$ID[i]] <- 1L

  # -- inter-cluster semantic distance ----------------------------------------
  id_to_matkey <- if (all(grepl("^GO:", rn))) {
    setNames(ego_s@result$ID, ego_s@result$ID)
  } else {
    setNames(ego_s@result$Description, ego_s@result$ID)
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

  # -- block gene-gene distance matrix ----------------------------------------
  gene_cv <- setNames(gene_class$cluster, gene_class$geneID)
  n_g     <- length(g_vec)
  D       <- matrix(0, n_g, n_g, dimnames = list(g_vec, g_vec))

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

  list(D = D, gene_class = gene_class, gene_term_mat = gene_term_mat,
       term_to_cluster = term_to_cluster, cluster_labels = cluster_labels,
       cl_ids = cl_ids, cl_dist = cl_dist)
}

# Builds the gene dendrogram ggplot from a precomputed distance matrix.
# gene_class must have columns: geneID, cluster, cluster_label.
.plot_gene_dendrogram <- function(D, gene_class, title_str, subtitle_str,
                                  palette = NULL) {
  gene_hc   <- hclust(as.dist(D), method = "ward.D")
  dend_data <- ggdendro::dendro_data(as.dendrogram(gene_hc), type = "rectangle")

  leaf_df <- dend_data$labels %>%
    left_join(gene_class %>% select(geneID, cluster, cluster_label),
              by = c("label" = "geneID"))

  bracket_df <- leaf_df %>%
    group_by(cluster, cluster_label) %>%
    summarise(y_min = min(x), y_max = max(x), .goals = "drop",
              .groups = "drop") %>%
    mutate(y_mid = (y_min + y_max) / 2)

  max_h   <- max(dend_data$segments$y, dend_data$segments$yend)
  bar_x   <- -0.10 * max_h
  label_x <- -0.22 * max_h

  cl_ids <- sort(unique(gene_class$cluster))
  n_cl   <- length(cl_ids)
  if (is.null(palette))
    palette <- setNames(
      if (n_cl <= 8) RColorBrewer::brewer.pal(max(3, n_cl), "Dark2")[seq_len(n_cl)]
      else           grDevices::hcl.colors(n_cl, palette = "Dark 3"),
      cl_ids)

  max_label_chars <- max(nchar(bracket_df$cluster_label), na.rm = TRUE)
  x_right_lim     <- label_x - max_label_chars * 0.09 * (max_h / 10)

  p <- ggplot() +
    geom_segment(data = dend_data$segments,
                 aes(x = y, y = x, xend = yend, yend = xend),
                 color = "grey35", linewidth = 0.5) +
    geom_text(data = leaf_df,
              aes(x = -0.02 * max_h, y = x, label = label,
                  color = factor(cluster)),
              angle = 0, hjust = 0, vjust = 0.5, size = 5) +
    geom_segment(data = bracket_df,
                 aes(x = bar_x, xend = bar_x,
                     y = y_min - 0.35, yend = y_max + 0.35,
                     color = factor(cluster)),
                 linewidth = 5, lineend = "round") +
    geom_text(data = bracket_df,
              aes(x = label_x, y = y_mid, label = cluster_label,
                  color = factor(cluster)),
              angle = 0, hjust = 0, vjust = 0.5, size = 4.5, fontface = "italic") +
    scale_color_manual(values = palette, guide = "none") +
    scale_x_reverse(limits = c(max_h * 1.06, x_right_lim),
                    name   = "Ward distance") +
    scale_y_continuous(expand = expansion(add = 0.8)) +
    labs(title = title_str, subtitle = subtitle_str, y = NULL) +
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
  p
}


# =============================================================================
# run_single_ont()  — full enrichment pipeline for one GO namespace
# =============================================================================
run_single_ont <- function(gene_df, ont, output_dir) {
  label <- ont_label(ont)
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # -- 1. GO enrichment -------------------------------------------------------
  message("\n=== ", label, " ===")
  message("Running enrichment (BH-adjusted p < 0.05) ...")

  ego <- enrichGO(
    gene          = gene_df$ENTREZID,
    universe      = keys(org.EcK12.eg.db, keytype = "ENTREZID"),
    OrgDb         = org.EcK12.eg.db,
    ont           = ont,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2,
    readable      = TRUE
  )

  n_raw <- sum(ego@result$p.adjust < 0.05)
  message("  Significant terms before simplification: ", n_raw)
  if (n_raw == 0) {
    message("  No significant terms — skipping ", ont, ".")
    return(invisible(NULL))
  }

  # -- 2. Simplify redundant terms --------------------------------------------
  message("Removing redundant terms (similarity cutoff 0.7) ...")
  ego_s <- clusterProfiler::simplify(ego, cutoff = 0.7, by = "p.adjust",
                                     select_fun = min)
  n_simplified <- nrow(ego_s@result)
  message("  Terms after simplification: ", n_simplified)

  # -- 3. Pairwise similarity -------------------------------------------------
  message("Computing pairwise Jaccard similarity ...")
  ego_sim <- pairwise_termsim(ego_s, method = "JC")

  # -- 4. Gap statistic for k -------------------------------------------------
  k <- .auto_k(ego_sim, k_max = min(10L, n_simplified - 1L))
  message("  Optimal cluster count (gap statistic): ", k)

  # -- 5. Treeplot ------------------------------------------------------------
  message("Building treeplot ...")
  p_tree <- treeplot(
    ego_sim,
    nCluster      = k,
    showCategory  = n_simplified,
    color         = "p.adjust",
    hexpand       = 0.85,
    tiplab_offset = 0.3,
    extend        = 3.5,
    cladelab_offset = 9.0
  ) +
    labs(
      title    = paste0(label, " — hierarchical term tree"),
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

  for (.i in seq_along(p_tree$layers)) {
    lyr <- p_tree$layers[[.i]]
    if (inherits(lyr$geom, "GeomText")) {
      lyr$aes_params$size  <- 4.2
      lyr$aes_params$hjust <- 0
    }
  }
  rm(.i, lyr)

  tree_h <- max(10, n_simplified * 0.38)
  tree_w <- round(tree_h * 16 / 9, 1)
  save_plot(p_tree, "go_treeplot", output_dir, w = tree_w, h = tree_h)

  # -- 6. Gene-concept network ------------------------------------------------
  message("Building cnetplot ...")
  set.seed(42)
  p_cnet <- cnetplot(ego_s, showCategory = n_simplified) +
    labs(
      title    = paste0("Gene-concept network — ", label),
      subtitle = paste0(
        nrow(gene_df), " E. coli K12 genes  |  ",
        n_simplified, " non-redundant terms  |  BH p.adj < 0.05"
      )
    ) +
    theme(
      plot.title    = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 11, color = "grey40")
    )

  for (.i in seq_along(p_cnet$layers)) {
    .cls <- class(p_cnet$layers[[.i]]$geom)[1]
    if (.cls %in% c("GeomTextRepel", "GeomText", "GeomLabel", "GeomLabelRepel"))
      p_cnet$layers[[.i]]$aes_params$size <- 6
  }
  rm(.i, .cls)
  save_plot(p_cnet, "go_cnetplot", output_dir, w = 32, h = 18)

  # -- 7. Gene dendrogram + distance matrix -----------------------------------
  message("Building gene dendrogram ...")
  dist_res <- .build_gene_dist(ego_s, ego_sim, k)
  D        <- dist_res$D

  p_gene <- .plot_gene_dendrogram(
    D          = D,
    gene_class = dist_res$gene_class,
    title_str  = paste0("Gene dendrogram — ", label),
    subtitle_str = paste0(
      "Leaves = genes  |  Jaccard distance on shared GO terms  |  k = ",
      k, " GO clusters"
    )
  )
  save_plot(p_gene, "go_single_class", output_dir, w = 32, h = 18)

  # -- 8. Save distance matrix ------------------------------------------------
  csv_path <- file.path(output_dir, "gene_distance_matrix.csv")
  write.csv(as.data.frame(D), csv_path)
  message("  Saved: ", csv_path)

  invisible(list(simplified = ego_s, full = ego))
}


# =============================================================================
# run_combined()  — merge BP + MF + CC, cluster genes by full GO profile
# =============================================================================
run_combined <- function(ego_list, gene_df, output_dir) {
  message("\n=== Combined (BP + MF + CC) ===")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  ego_list <- Filter(Negate(is.null), ego_list)
  if (length(ego_list) == 0) {
    message("  No enrichment results — skipping combined.")
    return(invisible(NULL))
  }

  ont_names <- names(ego_list)

  # -- 1. Merge result rows -----------------------------------------------------
  # Simplified: used for distance matrix and cnetplot (non-redundant terms).
  # Full (unsimplified): used for cluster labeling only — simplify() removes
  # specific terms like "peptidoglycan biosynthetic process" in favour of their
  # parents, so we go back to the full enrichment to find precise labels.
  all_results <- bind_rows(lapply(ont_names, function(o)
    ego_list[[o]]$simplified@result %>% mutate(ont = o)
  ))
  all_results_full <- bind_rows(lapply(ont_names, function(o)
    ego_list[[o]]$full@result %>% mutate(ont = o)
  ))

  # -- 2. Gene × term binary matrix ---------------------------------------------
  message("Building combined gene x term binary matrix ...")
  all_pairs <- all_results %>%
    select(ID, geneID) %>%
    separate_rows(geneID, sep = "/")

  g_all <- sort(unique(all_pairs$geneID))
  t_all <- unique(all_results$ID)
  n_g   <- length(g_all)

  gene_term_mat <- matrix(0L, nrow = n_g, ncol = length(t_all),
                          dimnames = list(g_all, t_all))
  for (i in seq_len(nrow(all_pairs)))
    gene_term_mat[all_pairs$geneID[i], all_pairs$ID[i]] <- 1L

  # -- 3. Jaccard gene-gene distance matrix ------------------------------------
  message("Computing Jaccard distance matrix (", n_g, " genes × ",
          length(t_all), " terms) ...")
  D_comb <- matrix(0, n_g, n_g, dimnames = list(g_all, g_all))
  for (i in seq_len(n_g))
    for (j in seq_len(i - 1L)) {
      vi <- gene_term_mat[i, ]; vj <- gene_term_mat[j, ]
      u  <- sum(vi | vj)
      D_comb[i, j] <- D_comb[j, i] <- if (u == 0) 1 else 1 - sum(vi & vj) / u
    }

  # -- 4. Gap statistic → k gene clusters --------------------------------------
  message("Finding optimal gene cluster count (gap statistic) ...")
  gene_hc <- hclust(as.dist(D_comb), method = "ward.D")
  k_max   <- min(10L, n_g - 1L)
  set.seed(42)
  gap <- cluster::clusGap(
    as.matrix(D_comb),
    FUN   = function(x, k) list(cluster = cutree(gene_hc, k)),
    K.max = k_max, B = 50, d.power = 2
  )
  k_genes <- as.integer(cluster::maxSE(gap$Tab[, "gap"], gap$Tab[, "SE.sim"],
                                        method = "Tibs2001SEmax"))
  message("  Optimal gene cluster count: ", k_genes)
  gene_clusters <- cutree(gene_hc, k = k_genes)

  # -- 5. Label each gene cluster with its most representative GO term ---------
  # Sort priority:
  #   1. BP before MF/CC — Biological Process terms are the most interpretable
  #      labels; MF/CC are used only when no BP term covers the cluster.
  #   2. Specificity = n_cl / Count descending — fraction of the term's
  #      annotated genes that fall in this cluster; favours specific terms over
  #      broad parents that happen to cover every gene.
  #   3. p.adjust ascending — tiebreak by enrichment significance.
  cluster_labels_df <- bind_rows(lapply(sort(unique(gene_clusters)), function(cl) {
    cl_genes <- names(gene_clusters)[gene_clusters == cl]
    scored <- all_results_full %>%
      rowwise() %>%
      mutate(n_cl        = length(intersect(strsplit(geneID, "/")[[1]], cl_genes)),
             specificity = n_cl / Count) %>%
      ungroup() %>%
      filter(n_cl > 0) %>%
      arrange(ont != "BP", desc(specificity), p.adjust) %>%
      slice(1)
    data.frame(
      cluster       = cl,
      cluster_label = if (nrow(scored) > 0) scored$Description[1]
                      else paste("Cluster", cl)
    )
  }))

  gene_class <- data.frame(geneID  = names(gene_clusters),
                            cluster = as.integer(gene_clusters)) %>%
    left_join(cluster_labels_df, by = "cluster")

  # -- 6. Gene dendrogram -------------------------------------------------------
  message("Building combined gene dendrogram ...")
  p_dend <- .plot_gene_dendrogram(
    D            = D_comb,
    gene_class   = gene_class,
    title_str    = "Gene dendrogram — combined GO (BP + MF + CC)",
    subtitle_str = paste0(
      "Leaves = genes  |  Jaccard distance across ", length(t_all),
      " GO terms (BP + MF + CC)  |  k = ", k_genes,
      " gene clusters (gap statistic)"
    )
  )
  save_plot(p_dend, "go_single_class", output_dir, w = 32, h = 18)

  # -- 7. Single unified cnetplot (all three ontologies merged) ----------------
  # Inject the merged result table into a copy of the BP enrichResult so that
  # cnetplot() can render it — it only reads @result, OrgDb, and readable.
  message("Building unified combined cnetplot ...")
  combined_ego         <- ego_list[[1]]$simplified
  combined_ego@result  <- all_results %>% select(-ont)
  n_show               <- nrow(combined_ego@result)
  set.seed(42)
  p_cnet <- cnetplot(combined_ego, showCategory = n_show) +
    labs(
      title    = "Gene-concept network — combined GO (BP + MF + CC)",
      subtitle = paste0(
        n_g, " E. coli K12 genes  |  ", n_show,
        " non-redundant terms across BP, MF, CC  |  BH p.adj < 0.05"
      )
    ) +
    theme(
      plot.title    = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 11, color = "grey40")
    )
  for (.i in seq_along(p_cnet$layers)) {
    .cls <- class(p_cnet$layers[[.i]]$geom)[1]
    if (.cls %in% c("GeomTextRepel", "GeomText", "GeomLabel", "GeomLabelRepel"))
      p_cnet$layers[[.i]]$aes_params$size <- 5
  }
  rm(.i, .cls)
  save_plot(p_cnet, "go_cnetplot_combined", output_dir, w = 36, h = 20)

  # -- 8. Save distance matrix --------------------------------------------------
  csv_path <- file.path(output_dir, "gene_distance_matrix.csv")
  write.csv(as.data.frame(D_comb), csv_path)
  message("  Saved: ", csv_path)

  invisible(D_comb)
}


# =============================================================================
# Main
# =============================================================================

# Convert gene symbols to Entrez IDs once
message("Converting ", length(gene_names), " gene symbols to Entrez IDs ...")
gene_df <- tryCatch(
  clusterProfiler::bitr(gene_names,
                        fromType = "SYMBOL", toType = "ENTREZID",
                        OrgDb    = org.EcK12.eg.db),
  error = function(e) stop("bitr() failed: ", conditionMessage(e))
)
missed <- setdiff(gene_names, gene_df$SYMBOL)
if (length(missed) > 0)
  warning(length(missed), " gene(s) not mapped and skipped: ",
          paste(missed, collapse = ", "))
if (nrow(gene_df) == 0) stop("No genes mapped. Aborting.")
message("  Mapped: ", nrow(gene_df), " / ", length(gene_names))

# Run per-ontology analyses
ego_results <- list()
for (ont in c("BP", "MF", "CC")) {
  ego_results[[ont]] <- run_single_ont(
    gene_df    = gene_df,
    ont        = ont,
    output_dir = file.path(output_base, ont)
  )
}

# Run combined analysis
run_combined(
  ego_list   = ego_results,
  gene_df    = gene_df,
  output_dir = file.path(output_base, "combined")
)

message("\nAll analyses complete.")
