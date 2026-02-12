# ============================================================
# PIPoL: Tree-aware decontamination via label propagation
# ============================================================


# ------------------------------------------------------------
# Internal: Estimate smoothing scale tau via BTS
# ------------------------------------------------------------
.BTS_tau <- function(
    physeq,
    level,
    alpha = 0.05,
    quantile_cutoff = 0.5
) {
  tree <- phyloseq::phy_tree(physeq)
  taxa_df <- as.data.frame(phyloseq::tax_table(physeq))
  D <- stats::cophenetic(tree)
  
  col_lc <- tolower(colnames(taxa_df))
  idx <- which(col_lc == tolower(level))
  if (length(idx) == 0) {
    stop("Taxonomic level not found: ", level)
  }
  
  tax_label <- taxa_df[[idx]]
  groups <- split(rownames(taxa_df), tax_label)
  
  d_max <- vapply(groups, function(g) {
    if (length(g) < 2) return(NA_real_)
    subD <- D[g, g, drop = FALSE]
    max(subD[upper.tri(subD)], na.rm = TRUE)
  }, numeric(1))
  
  d_cut <- stats::quantile(d_max, probs = quantile_cutoff, na.rm = TRUE)
  tau <- d_cut / (-log(alpha))
  
  as.numeric(tau)
}


# ------------------------------------------------------------
# Internal: PIPoL tree-aware label propagation
# ------------------------------------------------------------
.pipol_tree_core <- function(
    physeq,
    outcome_var,
    tau,
    confounder_vars = NULL,
    weighted = TRUE
) {
  otu <- as(phyloseq::otu_table(physeq), "matrix")
  if (!phyloseq::taxa_are_rows(physeq)) otu <- t(otu)
  
  meta <- as(phyloseq::sample_data(physeq), "data.frame")
  tree <- phyloseq::phy_tree(physeq)
  
  outcome <- factor(meta[[outcome_var]])
  y_mat <- stats::model.matrix(~ outcome - 1)
  colnames(y_mat) <- levels(outcome)
  
  # IPW or centering
  if (!is.null(confounder_vars)) {
    X <- stats::model.matrix(
      as.formula(paste("~", paste(confounder_vars, collapse = "+"))),
      data = meta
    )
    ps_mat <- sapply(seq_len(ncol(y_mat)), function(k) {
      fit <- stats::glm(y_mat[, k] ~ X - 1, family = stats::binomial())
      pmax(stats::predict(fit, type = "response"), 1e-6)
    })
    y_eff <- y_mat / ps_mat
  } else {
    y_eff <- scale(y_mat, center = TRUE, scale = FALSE)
  }
  
  # Host–taxon matrix
  if (weighted) {
    W <- sweep(otu, 2, colSums(otu), FUN = "/")
    W[is.na(W)] <- 0
  } else {
    W <- (otu > 0) * 1
  }
  
  R0 <- W %*% y_eff
  
  D <- stats::cophenetic(tree)
  K <- exp(-D / tau)
  K <- sweep(K, 1, rowSums(K), FUN = "/")
  
  R <- K %*% R0
  rownames(R) <- rownames(otu)
  colnames(R) <- colnames(y_mat)
  
  R
}


# ------------------------------------------------------------
# Public API
# ------------------------------------------------------------


#' Enrichment–Prevalence Analysis (Tree-aware)
#'
#' Perform group-wise co-enrichment and co-existence analysis
#' using tree-aware PIPoL propagation.
#'
#' @param physeq A phyloseq object.
#' @param outcome_var Sample-level outcome variable.
#' @param group Group level to analyze.
#' @param tau_level Taxonomic level used for smoothing.
#' @param alpha Significance level controlling tree decay scale.
#' @param confounders Optional confounders for IPW.
#' @param p_upper Upper-tail p-value cutoff for enrichment.
#' @param p_lower Lower-tail p-value cutoff for depletion.
#' @param min_prevalence Minimum prevalence required.
#' @export
enrichment_prevalence <- function(
    physeq,
    outcome_var,
    group,
    tau_level = c("Genus", "Family", "Species"),
    alpha = 0.05,
    confounders = NULL,
    p_upper = 0.05,
    p_lower = 0.05,
    min_prevalence = 0,
    vline = NULL,
    hline = NULL
) {
  
  requireNamespace("phyloseq", quietly = TRUE)
  requireNamespace("ggplot2", quietly = TRUE)
  
  tau_level <- match.arg(tau_level)
  
  # ------------------------------------------------------------
  #  Compute tau
  # ------------------------------------------------------------
  tau_val <- .BTS_tau(
    physeq,
    level = tau_level,
    alpha = alpha
  )
  
  # ------------------------------------------------------------
  #  Generate Rtrees
  # ------------------------------------------------------------
  weighted_Rtree <- .pipol_tree_core(
    physeq,
    outcome_var = outcome_var,
    tau = tau_val,
    confounder_vars = confounders,
    weighted = TRUE
  )
  
  unweighted_Rtree <- .pipol_tree_core(
    physeq,
    outcome_var = outcome_var,
    tau = tau_val,
    confounder_vars = confounders,
    weighted = FALSE
  )
  
  if (!(group %in% colnames(weighted_Rtree))) {
    stop("Group not found in outcome levels.")
  }
  
  # ------------------------------------------------------------
  #  OTU alignment
  # ------------------------------------------------------------
  otu <- as(phyloseq::otu_table(physeq), "matrix")
  if (!phyloseq::taxa_are_rows(physeq)) otu <- t(otu)
  
  common_features <- intersect(rownames(weighted_Rtree), rownames(otu))
  if (length(common_features) == 0) {
    stop("No overlapping features.")
  }
  
  otu <- otu[common_features, , drop = FALSE]
  weighted_Rtree <- weighted_Rtree[common_features, , drop = FALSE]
  unweighted_Rtree <- unweighted_Rtree[common_features, , drop = FALSE]
  
  # ------------------------------------------------------------
  #  Extract group samples
  # ------------------------------------------------------------
  meta <- as(phyloseq::sample_data(physeq), "data.frame")
  samples <- rownames(meta[meta[[outcome_var]] == group, , drop = FALSE])
  
  if (length(samples) == 0) {
    stop("No samples found for group.")
  }
  
  otu_sub <- otu[, samples, drop = FALSE]
  
  prevalence <- rowSums(otu_sub > 0) / ncol(otu_sub)
  
  denom <- colSums(otu_sub)
  denom[denom == 0] <- NA
  otu_rel <- sweep(otu_sub, 2, denom, "/")
  mean_abundance <- rowMeans(otu_rel, na.rm = TRUE)
  
  keep <- prevalence > min_prevalence &
    is.finite(prevalence) &
    is.finite(mean_abundance)
  
  if (!any(keep)) {
    stop("No taxa remain after filtering.")
  }
  
  prevalence <- prevalence[keep]
  mean_abundance <- mean_abundance[keep]
  features <- names(prevalence)
  
  weighted_Rtree <- weighted_Rtree[features, , drop = FALSE]
  unweighted_Rtree <- unweighted_Rtree[features, , drop = FALSE]
  
  # ------------------------------------------------------------
  #  Z-scores
  # ------------------------------------------------------------
  wz <- scale(weighted_Rtree)
  uz <- scale(unweighted_Rtree)
  
  pw_u <- pnorm(wz, lower.tail = FALSE)
  pu_u <- pnorm(uz, lower.tail = FALSE)
  pw_l <- pnorm(wz, lower.tail = TRUE)
  pu_l <- pnorm(uz, lower.tail = TRUE)
  
  category <- ifelse(
    pw_u[, group] < p_upper & pu_u[, group] < p_upper,
    "Core",
    ifelse(
      pw_l[, group] < p_lower & pu_l[, group] < p_lower,
      "Transient",
      "Intermediate"
    )
  )
  
  df <- data.frame(
    Feature = features,
    Prevalence = prevalence,
    Mean_abundance = mean_abundance,
    Weighted_z = wz[, group],
    Unweighted_z = uz[, group],
    Category = factor(category,
                      levels = c("Core", "Intermediate", "Transient")),
    row.names = features
  )
  
  # ------------------------------------------------------------
  #  Plots
  # ------------------------------------------------------------
  color_map <- c(
    Core = "#D7263D",
    Intermediate = "gray70",
    Transient = "#1E90FF"
  )
  
  shape_map <- c(
    Core = 17,
    Intermediate = 16,
    Transient = 15
  )
  
  # 1️⃣ Z vs Z
  p1 <- ggplot2::ggplot(
    df,
    ggplot2::aes(Unweighted_z, Weighted_z,
                 color = Category,
                 shape = Category)
  ) +
    ggplot2::geom_point(alpha = 0.8, size = 2.2) +
    ggplot2::scale_color_manual(values = color_map, drop = FALSE) +
    ggplot2::scale_shape_manual(values = shape_map, drop = FALSE) +
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::labs(
      x = "Coexist Z Score",
      y = "Coenrichment Z Score"
    )
  
  # 2️⃣ Prevalence vs Coenrichment
  p2 <- ggplot2::ggplot(
    df,
    ggplot2::aes(Prevalence, Weighted_z)
  ) +
    ggplot2::geom_point(color = "gray30",
                        alpha = 0.7,
                        size = 2) +
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::labs(
      x = "Prevalence",
      y = "Coenrichment Z Score"
    )
  
  # 3️⃣ Mean abundance vs Coenrichment
  p3 <- ggplot2::ggplot(
    df,
    ggplot2::aes(Mean_abundance, Weighted_z)
  ) +
    ggplot2::geom_point(color = "gray30",
                        alpha = 0.7,
                        size = 2) +
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::labs(
      x = "Mean Relative Abundance",
      y = "Coenrichment Z Score"
    )
  
  # 4️⃣ Prevalence vs Mean abundance
  p4 <- ggplot2::ggplot(
    df,
    ggplot2::aes(Prevalence, Mean_abundance,
                 color = Category,
                 shape = Category)
  ) +
    ggplot2::geom_point(alpha = 0.8, size = 2.2) +
    ggplot2::scale_color_manual(values = color_map, drop = FALSE) +
    ggplot2::scale_shape_manual(values = shape_map, drop = FALSE) +
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::labs(
      x = "Prevalence",
      y = "Mean Relative Abundance"
    )
  
  if (!is.null(vline)) {
    p4 <- p4 +
      ggplot2::geom_vline(
        xintercept = vline,
        linetype = "dashed",
        color = "gray50"
      )
  }
  
  if (!is.null(hline)) {
    p4 <- p4 +
      ggplot2::geom_hline(
        yintercept = hline,
        linetype = "dashed",
        color = "gray50"
      )
  }
  
  return(list(
    table = df,
    plot_enrichment = p1,
    plot_prev_enrich = p2,
    plot_abund_enrich = p3,
    plot_prev_abund = p4,
    tau = tau_val
  ))
  
}
