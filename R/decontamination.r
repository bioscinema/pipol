# ============================================================
# RIPoL: Tree-aware decontamination via label propagation
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
# Internal: RIPoL tree-aware label propagation
# ------------------------------------------------------------
.ripol_tree_core <- function(
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
# Internal: Core / transient classification from Rtrees
# ------------------------------------------------------------
.decontam_from_rtrees <- function(
    Rtree_weighted,
    Rtree_unweighted,
    label,
    gamma_weighted,
    gamma_unweighted,
    physeq = NULL
) {
  stopifnot(
    all(dim(Rtree_weighted) == dim(Rtree_unweighted)),
    all(rownames(Rtree_weighted) == rownames(Rtree_unweighted))
  )
  
  z_w <- scale(Rtree_weighted[, label])
  z_u <- scale(Rtree_unweighted[, label])
  
  p_w <- stats::pnorm(z_w, lower.tail = FALSE)
  p_u <- stats::pnorm(z_u, lower.tail = FALSE)
  
  thresh_w <- stats::qnorm(1 - gamma_weighted)
  thresh_u <- stats::qnorm(1 - gamma_unweighted)
  
  category <- ifelse(
    z_w > thresh_w & z_u > thresh_u, "High risk",
    ifelse(z_w < -thresh_w & z_u < -thresh_u, "Low risk", "Medium risk")
  )
  
  df <- data.frame(
    OTU = rownames(Rtree_weighted),
    CoexistZ = z_u,
    CoenrichZ = z_w,
    Category = factor(category,
                      levels = c("High risk", "Medium risk", "Low risk"))
  )
  
  p_scatter <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = CoexistZ, y = CoenrichZ,
                 color = Category, shape = Category)
  ) +
    ggplot2::geom_point(alpha = 0.8, size = 2.2) +
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::scale_color_manual(
      values = c("High risk" = "#D7263D",
                 "Medium risk" = "gray70",
                 "Low risk" = "#1E90FF")
    ) +
    ggplot2::labs(
      x = "Coexistence Z score",
      y = "Co-enrichment Z score"
    )
  
  list(
    table = df,
    plot = p_scatter
  )
}


# ------------------------------------------------------------
# Public API
# ------------------------------------------------------------

#' Tree-aware microbiome decontamination via RIPoL
#'
#' Identifies core, intermediate, and transient taxa by jointly modeling
#' prevalence-based co-existence and abundance-based co-enrichment using
#' phylogeny-aware label propagation.
#'
#' @param physeq A phyloseq object
#' @param outcome_var Sample-level outcome variable
#' @param label Target group label
#' @param tau Taxonomic level used for smoothing ("Genus", "Family", "Species")
#' @param alpha Significance level controlling tree decay scale
#' @param gamma_weighted Upper-tail cutoff for weighted propagation
#' @param gamma_unweighted Upper-tail cutoff for unweighted propagation
#' @param confounders Optional confounder variables for IPW
#'
#' @return A list containing a results table and ggplot object
#' @export
decontamination <- function(
    physeq,
    outcome_var,
    label,
    tau = c("Genus", "Family", "Species"),
    alpha = 0.05,
    gamma_weighted = 0.05,
    gamma_unweighted = 0.05,
    confounders = NULL
) {
  tau <- match.arg(tau)
  
  tau_val <- .BTS_tau(
    physeq,
    level = tau,
    alpha = alpha
  )
  
  R_w <- .ripol_tree_core(
    physeq,
    outcome_var = outcome_var,
    tau = tau_val,
    confounder_vars = confounders,
    weighted = TRUE
  )
  
  R_u <- .ripol_tree_core(
    physeq,
    outcome_var = outcome_var,
    tau = tau_val,
    confounder_vars = confounders,
    weighted = FALSE
  )
  
  .decontam_from_rtrees(
    R_w,
    R_u,
    label = label,
    gamma_weighted = gamma_weighted,
    gamma_unweighted = gamma_unweighted,
    physeq = physeq
  )
}
