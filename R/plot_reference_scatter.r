#' Reference-based scatter plot for RIPoL results
#'
#' Visualize coexistence vs. co-enrichment scores produced by RIPoL,
#' colored by reference classification (e.g. Oral / Contaminant),
#' optionally faceted by prevalence strata defined by user-specified
#' bins or quantile-based grouping.
#'
#' @param results_df Data frame containing RIPoL results.
#' @param score_x Column name for x-axis score (default: "coexists.p").
#' @param score_y Column name for y-axis score (default: "coenrich.p").
#' @param group_col Column name for reference classification.
#' @param abundance_col Column name for abundance (used as point size).
#' @param prevalence_col Column name for prevalence (numeric or factor).
#' @param facet Logical; whether to facet by prevalence.
#' @param prevalence_bins Numeric vector of cut points for prevalence.
#'   If provided, overrides \code{n_prevalence_bins}.
#' @param n_prevalence_bins Integer; number of quantile-based prevalence bins.
#' @param palette Named vector of colors for reference groups.
#' @param alpha Point transparency.
#'
#' @return A ggplot object.
#' @export
plot_reference_scatter <- function(
    results_df,
    score_x = "coexists.p",
    score_y = "coenrich.p",
    group_col = "Group",
    abundance_col = "Abundance",
    prevalence_col = "Prevalence",
    facet = TRUE,
    prevalence_bins = NULL,
    n_prevalence_bins = NULL,
    palette = c(
      "Oral" = "blue",
      "Contaminant" = "red",
      "Ambiguous" = "gray80"
    ),
    alpha = 0.3
) {
  
  stopifnot(
    score_x %in% colnames(results_df),
    score_y %in% colnames(results_df),
    group_col %in% colnames(results_df)
  )
  
  df <- results_df
  
  # --------------------------------------------------
  # Construct prevalence facet variable if requested
  # --------------------------------------------------
  if (facet && prevalence_col %in% colnames(df)) {
    
    prev <- df[[prevalence_col]]
    
    if (!is.null(prevalence_bins)) {
      df$.facet_prevalence <- cut(
        prev,
        breaks = prevalence_bins,
        include.lowest = TRUE,
        right = FALSE
      )
      
    } else if (!is.null(n_prevalence_bins)) {
      qs <- stats::quantile(
        prev,
        probs = seq(0, 1, length.out = n_prevalence_bins + 1),
        na.rm = TRUE
      )
      qs <- unique(qs)
      
      df$.facet_prevalence <- cut(
        prev,
        breaks = qs,
        include.lowest = TRUE
      )
      
    } else {
      df$.facet_prevalence <- factor(prev)
    }
  }
  
  # --------------------------------------------------
  # Build plot
  # --------------------------------------------------
  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = .data[[score_x]],
      y = .data[[score_y]],
      color = .data[[group_col]],
      size = if (abundance_col %in% colnames(df)) .data[[abundance_col]] else NULL
    )
  ) +
    ggplot2::geom_point(alpha = alpha) +
    ggplot2::scale_color_manual(
      values = palette,
      name = "Reference\nClassification"
    ) +
    ggplot2::labs(
      x = paste0("Score (", score_x, ")"),
      y = paste0("Score (", score_y, ")")
    ) +
    ggplot2::guides(size = "none") +
    ggplot2::theme_classic(base_size = 14)
  
  if (facet && ".facet_prevalence" %in% colnames(df)) {
    p <- p + ggplot2::facet_grid(
      stats::as.formula("~ .facet_prevalence")
    )
  }
  
  return(p)
}
