#' Enrichment–Prevalence Analysis
#'
#' Perform group-wise co-enrichment and co-existence analysis
#' with automatic handling of sparse and dense compositional data.
#'
#' @param physeq A phyloseq object.
#' @param weighted_Rtree Weighted association matrix (features x groups).
#' @param unweighted_Rtree Unweighted association matrix (features x groups).
#' @param group_var Sample metadata variable defining groups.
#' @param group Group level to analyze.
#' @param p_upper Upper-tail p-value cutoff for enrichment.
#' @param p_lower Lower-tail p-value cutoff for depletion.
#' @param min_prevalence Minimum prevalence required to retain a feature.
#' @param vline Optional vertical cutoff for prevalence-abundance plot.
#' @param hline Optional horizontal cutoff for prevalence-abundance plot.
#'
#' @return A list containing:
#' \itemize{
#'   \item table: data.frame with enrichment statistics
#'   \item plot_enrichment: co-existence vs co-enrichment Z-scores
#'   \item plot_prev_enrich: prevalence vs co-enrichment
#'   \item plot_abund_enrich: abundance vs co-enrichment
#'   \item plot_prev_abund: prevalence vs abundance
#' }
#'
#' @export
enrichment_prevalence <- function(
    physeq,
    weighted_Rtree,
    unweighted_Rtree,
    group_var,
    group,
    p_upper = 0.05,
    p_lower = 0.05,
    min_prevalence = 0,
    vline = NULL,
    hline = NULL
) {
  
  requireNamespace("phyloseq", quietly = TRUE)
  requireNamespace("ggplot2", quietly = TRUE)
  
  ## ----------------------------
  ## Input checks
  ## ----------------------------
  stopifnot(
    identical(dim(weighted_Rtree), dim(unweighted_Rtree)),
    identical(rownames(weighted_Rtree), rownames(unweighted_Rtree)),
    identical(colnames(weighted_Rtree), colnames(unweighted_Rtree))
  )
  
  if (!(group %in% colnames(weighted_Rtree))) {
    stop("`group` not found in Rtree columns.")
  }
  
  meta <- as(phyloseq::sample_data(physeq), "data.frame")
  if (!(group_var %in% colnames(meta))) {
    stop("`group_var` not found in sample_data.")
  }
  
  samples <- rownames(meta[meta[[group_var]] == group, , drop = FALSE])
  if (length(samples) == 0) {
    stop("No samples found for group = ", group)
  }
  
  ## ----------------------------
  ## OTU alignment
  ## ----------------------------
  otu <- as(phyloseq::otu_table(physeq), "matrix")
  if (!phyloseq::taxa_are_rows(physeq)) otu <- t(otu)
  
  common_features <- intersect(rownames(weighted_Rtree), rownames(otu))
  if (length(common_features) == 0) {
    stop("No overlapping features between physeq and Rtree.")
  }
  
  otu <- otu[common_features, , drop = FALSE]
  weighted_Rtree <- weighted_Rtree[common_features, , drop = FALSE]
  unweighted_Rtree <- unweighted_Rtree[common_features, , drop = FALSE]
  
  ## ----------------------------
  ## Prevalence & abundance
  ## ----------------------------
  otu_sub <- otu[, samples, drop = FALSE]
  
  prevalence <- rowSums(otu_sub > 0) / ncol(otu_sub)
  
  denom <- colSums(otu_sub)
  denom[denom == 0] <- NA
  otu_rel <- sweep(otu_sub, 2, denom, "/")
  mean_abundance <- rowMeans(otu_rel, na.rm = TRUE)
  
  keep <- is.finite(prevalence) &
    is.finite(mean_abundance) &
    prevalence > min_prevalence
  
  if (!any(keep)) {
    stop("No features remain after prevalence filtering.")
  }
  
  prevalence <- prevalence[keep]
  mean_abundance <- mean_abundance[keep]
  features <- names(prevalence)
  
  weighted_Rtree <- weighted_Rtree[features, , drop = FALSE]
  unweighted_Rtree <- unweighted_Rtree[features, , drop = FALSE]
  
  ## ----------------------------
  ## Z-scores & p-values
  ## ----------------------------
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
    Category = factor(category, levels = c("Core", "Intermediate", "Transient")),
    row.names = features
  )
  
  ## ----------------------------
  ## Plots
  ## ----------------------------
  color_map <- c(
    Core = "#D7263D",
    Intermediate = "gray70",
    Transient = "#1E90FF"
  )
  shape_map <- c(Core = 17, Intermediate = 16, Transient = 15)
  
  p_enrich <- ggplot2::ggplot(
    df,
    ggplot2::aes(Unweighted_z, Weighted_z, color = Category, shape = Category)
  ) +
    ggplot2::geom_point(alpha = 0.8, size = 2.2) +
    ggplot2::scale_color_manual(values = color_map, drop = FALSE) +
    ggplot2::scale_shape_manual(values = shape_map, drop = FALSE) +
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::labs(
      x = "Co-existence Z-score",
      y = "Co-enrichment Z-score"
    )
  
  p_prev_enrich <- ggplot2::ggplot(df, ggplot2::aes(Prevalence, Weighted_z)) +
    ggplot2::geom_point(color = "gray30", alpha = 0.7, size = 2) +
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::labs(x = "Prevalence", y = "Co-enrichment Z-score")
  
  p_abund_enrich <- ggplot2::ggplot(df, ggplot2::aes(Mean_abundance, Weighted_z)) +
    ggplot2::geom_point(color = "gray30", alpha = 0.7, size = 2) +
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::labs(x = "Mean Relative Abundance", y = "Co-enrichment Z-score")
  
  p_prev_abund <- ggplot2::ggplot(
    df,
    ggplot2::aes(Prevalence, Mean_abundance, color = Category, shape = Category)
  ) +
    ggplot2::geom_point(alpha = 0.8, size = 2.2) +
    ggplot2::scale_color_manual(values = color_map, drop = FALSE) +
    ggplot2::scale_shape_manual(values = shape_map, drop = FALSE) +
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::labs(x = "Prevalence", y = "Mean Relative Abundance")
  
  if (!is.null(vline)) {
    p_prev_abund <- p_prev_abund +
      ggplot2::geom_vline(xintercept = vline, linetype = "dashed", color = "gray50")
  }
  if (!is.null(hline)) {
    p_prev_abund <- p_prev_abund +
      ggplot2::geom_hline(yintercept = hline, linetype = "dashed", color = "gray50")
  }
  
  list(
    table = df,
    plot_enrichment = p_enrich,
    plot_prev_enrich = p_prev_enrich,
    plot_abund_enrich = p_abund_enrich,
    plot_prev_abund = p_prev_abund
  )
}
