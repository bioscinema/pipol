BTS <- function(physeq, level = c("Species", "Genus", "Family"), alpha = 0.05, quantile_cutoff = 0.5) {
  require(phyloseq)
  require(ape)
  require(dplyr)
  
  tree <- phy_tree(physeq)
  taxa_df <- as.data.frame(tax_table(physeq))
  D <- cophenetic(tree)  # pairwise patristic distance matrix
  
  result <- data.frame(Taxonomic_Level = character(), Recommended_Tau = numeric())
  
  # lowercase column names for matching
  colnames_lc <- tolower(colnames(taxa_df))
  
  for (L in level) {
    L_match <- which(colnames_lc == tolower(L))
    if (length(L_match) == 0) next  # skip if not found
    L_real <- colnames(taxa_df)[L_match]  # get original column name
    
    taxa_df$tax_label <- taxa_df[[L_real]]
    taxon_groups <- split(rownames(taxa_df), taxa_df$tax_label)
    
    d_max_vec <- sapply(taxon_groups, function(grp) {
      if (length(grp) < 2) return(NA)
      subD <- D[grp, grp, drop = FALSE]
      max(subD[upper.tri(subD)], na.rm = TRUE)
    })
    
    d_cut <- quantile(d_max_vec, probs = quantile_cutoff, na.rm = TRUE)
    tau_L <- d_cut / (-log(alpha))
    
    result <- rbind(result, data.frame(Taxonomic_Level = L_real, Recommended_Tau = tau_L))
  }
  
  return(result)
}


RIPoL_tree <- function(physeq, outcome_var, confounder_vars = NULL, tau, weighted = TRUE) {
  library(phyloseq)
  library(ape)
  library(dplyr)
  library(stats)
  
  # 1. read data
  otu <- as(otu_table(physeq),"matrix")
  if (!taxa_are_rows(physeq)) otu <- t(otu)
  meta <- as(sample_data(physeq),"data.frame")
  tree <- phy_tree(physeq)
  
  # 2. one hot encoding
  outcome <- factor(meta[[outcome_var]])
  y_mat <- model.matrix(~ outcome - 1)  
  colnames(y_mat) <- gsub("^outcome", "", colnames(y_mat))
  K <- ncol(y_mat)
  
  # 3. Estimate propensity score and compute IPW-weighted label
  if (!is.null(confounder_vars)) {
    X <- model.matrix(as.formula(paste("~", paste(confounder_vars, collapse = "+"))), data = meta)
    ps_mat <- sapply(1:K, function(k) {
      fit <- glm(y_mat[, k] ~ X - 1, family = binomial)
      pmax(predict(fit, type = "response"), 1e-6)
    })
    colnames(ps_mat) <- colnames(y_mat)
    y_weighted <- y_mat / ps_mat
  } else {
    # No confounder: center outcome
    y_weighted <- scale(y_mat, center = TRUE, scale = FALSE)
  }
  
  # 4. Host–taxon bipartite matrix: weighted or binary
  if (weighted) {
    W <- sweep(otu, 2, colSums(otu), FUN = "/")  # TSS normalization
    W[is.na(W)] <- 0
  } else {
    W <- (otu > 0) * 1  # Binary presence/absence
  }
  
  # 5. Host to Taxon 
  R0 <- W %*% y_weighted 
  
  # 6. Tree kernel smoothing
  D <- cophenetic(tree)
  K <- exp(-D / tau)
  K <- sweep(K, 1, rowSums(K), FUN = "/")  
  

  R_final <- K %*% R0  # m × K matrix
  
  rownames(R_final) <- rownames(otu)
  colnames(R_final) <- colnames(y_mat)
  
  return(R_final)
}

decontamination <- function(Rtree_weighted, Rtree_unweighted, 
                               gamma_weighted = 0.05, gamma_unweighted = 0.05,
                               physeq = NULL,label) {########################################

  stopifnot(all(dim(Rtree_weighted) == dim(Rtree_unweighted)))
  stopifnot(all(rownames(Rtree_weighted) == rownames(Rtree_unweighted)))
  stopifnot(all(colnames(Rtree_weighted) == colnames(Rtree_unweighted)))
  
  taxa_ids <- rownames(Rtree_weighted)
  class_labels <- colnames(Rtree_weighted)
  core_taxa_list <- list()
  enrichment_list <- list()
  existence_list <- list()
#
  score_w <- Rtree_weighted[, label]
  score_u <- Rtree_unweighted[, label]
  score_w<-scale(score_w)
  score_u<-scale(score_u)
  p_w<-1-pnorm(score_w)
  p_u<-1-pnorm(score_u)

#
  # core_taxa <- intersect(top_w, top_u)
  
  # core_taxa_list[[label]] <- core_taxa
  enrichment_list<-row.names(Rtree_weighted[p_w<gamma_weighted,])
  existence_list <- row.names(Rtree_unweighted[p_u<gamma_unweighted,])

  
  # all_core_taxa <- unique(unlist(core_taxa_list))
  enrichment_core_taxa <- unique(enrichment_list)
  existence_core_taxa <- unique(existence_list)
  if (!is.null(physeq)) {
    tax_tab <- as(tax_table(physeq), "matrix")
    # taxonomy_subset <- tax_tab[all_core_taxa, , drop = FALSE]
    enrichment_subset <- tax_tab[enrichment_core_taxa, , drop = FALSE]
    existence_subset <- tax_tab[existence_core_taxa, , drop = FALSE]
  } else {
    enrichment_subset <- NULL
    existence_subset <- NULL
  }
  
  df_plot <- data.frame(
    OTU = taxa_ids,
    coexist_z = as.numeric(score_u),
    coenrich_z = as.numeric(score_w)
  )
  
  # 计算 95% 分位数阈值
  thresh_x <- qnorm(1-gamma_unweighted)
  thresh_y <- qnorm(1-gamma_weighted)
  
  df_plot$risk <- "Low risk"
  df_plot$risk[df_plot$coexist_z > thresh_x & df_plot$coenrich_z > thresh_y] <- "High risk"
  df_plot$risk[
    (df_plot$coexist_z > thresh_x & df_plot$coenrich_z <= thresh_y) |
      (df_plot$coexist_z <= thresh_x & df_plot$coenrich_z > thresh_y)
  ] <- "Medium risk"
  
  # 映射颜色和形状
  color_map <- c(
    "High risk" = "#D7263D",
    "Medium risk" = "gray70",
    "Low risk" = "#1E90FF"
  )
  
  shape_map <- c(
    "High risk" = 17,
    "Medium risk" = 16,
    "Low risk" = 15
  )
  library(ggplot2)
  plot_scatter <- ggplot(df_plot, aes(x = coexist_z, y = coenrich_z)) +
    geom_point(aes(color = risk, shape = risk), size = 2) +
    scale_color_manual(values = color_map) +
    scale_shape_manual(values = shape_map) +
    theme_classic(base_size = 14) +
    labs(
      x = "Coexist Z Score",
      y = "Coenrichment Z Score",
      color = "Category",
      shape = "Category"
    )
  return(list(enrich_list = enrichment_list,
              exist_list = existence_list,
              enrichment_table = enrichment_subset,
              existence_table = existence_subset,
              scatter_plot = plot_scatter
              ))
}

generate_core_taxa <- function(Rtree_weighted, Rtree_unweighted, 
                               gamma_weighted = 0.05, gamma_unweighted = 0.05,
                               physeq = NULL) {

  stopifnot(all(dim(Rtree_weighted) == dim(Rtree_unweighted)))
  stopifnot(all(rownames(Rtree_weighted) == rownames(Rtree_unweighted)))
  stopifnot(all(colnames(Rtree_weighted) == colnames(Rtree_unweighted)))
  
  taxa_ids <- rownames(Rtree_weighted)
  class_labels <- colnames(Rtree_weighted)
  core_taxa_list <- list()
  
  for (label in class_labels) {
    score_w <- Rtree_weighted[, label]
    score_u <- Rtree_unweighted[, label]

    n_w <- max(1, floor(length(score_w) * gamma_weighted))
    n_u <- max(1, floor(length(score_u) * gamma_unweighted))
    
    top_w <- names(sort(score_w, decreasing = TRUE))[1:n_w]
    top_u <- names(sort(score_u, decreasing = TRUE))[1:n_u]
    
    core_taxa <- intersect(top_w, top_u)
    core_taxa_list[[label]] <- core_taxa
  }
  
  all_core_taxa <- unique(unlist(core_taxa_list))
  
  if (!is.null(physeq)) {
    tax_tab <- as(tax_table(physeq), "matrix")
    taxonomy_subset <- tax_tab[all_core_taxa, , drop = FALSE]
  } else {
    taxonomy_subset <- NULL
  }
  
  return(list(core_taxa_by_class = core_taxa_list,
              taxonomy_table = taxonomy_subset))
}

train_HINT_model <- function(W, Rtree, covariates, labels, 
                             n_iter = 1000, lr = 0.01, seed = 123, verbose = TRUE) {
  set.seed(seed)
  
  n <- nrow(W)
  m <- ncol(W)
  p <- ncol(covariates)
  K <- length(unique(labels))
  d <- ncol(Rtree)
  
  # One-hot encode labels
  y_mat <- model.matrix(~ factor(labels) - 1)
  colnames(y_mat) <- levels(factor(labels))
  
  # Initialize parameters
  B <- matrix(rnorm(p * K, mean = 0, sd = 0.01), nrow = p)
  V <- matrix(rnorm(d * K, mean = 0, sd = 0.01), nrow = d)
  
  # Softmax function
  softmax <- function(x) {
    exp_x <- exp(x - apply(x, 1, max))
    exp_x / rowSums(exp_x)
  }
  
  loss_history <- numeric(n_iter)
  
  for (iter in 1:n_iter) {
    # Forward pass
    E <- W %*% Rtree           # ecological embedding: n × d
    eta <- covariates %*% B + E %*% V  # total logit: n × K
    p_hat <- softmax(eta)
    
    # Loss (cross-entropy)
    loss <- -sum(y_mat * log(p_hat + 1e-10)) / n
    loss_history[iter] <- loss
    
    if (verbose && iter %% 100 == 0) {
      cat("Iter:", iter, "Loss:", round(loss, 4), "\n")
    }
    
    # Gradients
    grad_eta <- (p_hat - y_mat) / n  # n × K
    grad_B <- t(covariates) %*% grad_eta  # p × K
    grad_V <- t(E) %*% grad_eta           # d × K
    
    # Gradient descent update
    B <- B - lr * grad_B
    V <- V - lr * grad_V
  }
  
  return(list(B = B, V = V, loss = loss_history))
}




ContamScore <- function(physeq, Rtree) {
  library(phyloseq)
  library(matrixStats)
  
  # 1. OTU table
  otu <- as.matrix(otu_table(physeq))
  if (!taxa_are_rows(physeq)) otu <- t(otu)
  prevalence <- rowSums(otu > 0) / ncol(otu)
  
  # 2. Normalize to TSS (sample-wise)
  W <- sweep(otu, 2, colSums(otu), FUN = "/")
  W[is.na(W)] <- 0  # handle 0/0
  
  # 3. For each taxon, compute ||r_j||^2
  strength_vec <- rowSums(Rtree^2)  # length m
  
  # 4. Construct influence matrix: M_ij = W_ij * ||r_j||^2
  influence_matrix <- t(t(W) * strength_vec)  # W: m × n; -> t(W): n × m; strength_vec: m
  
  # 5. For each OTU j: strength = sum over i, dispersion = CV over i
  influence_strength <- rowSums(influence_matrix)  # m-vector
  influence_mean <- rowMeans(influence_matrix)
  influence_sd <- apply(influence_matrix, 1, sd)
  influence_cv <- influence_sd / (influence_mean + 1e-8)
  
  # 6. Return data.frame
  df <- data.frame(
    OTU = rownames(Rtree),
    InfluenceStrength = influence_strength,
    InfluenceCV = influence_cv,
    Prevalence = prevalence[rownames(Rtree)]
  )
  
  return(df)
}

plot_cv_vs_strength <- function(contam_df) {
  library(ggplot2)
  
  p <- ggplot(contam_df, aes(
    x = -log10(InfluenceStrength),
    y = InfluenceCV,
    size = Prevalence + 0.25
  )) +
    geom_point(alpha = 0.6, color = "#003366") +
    theme_minimal(base_size = 14) +
    labs(
      title = "Contamination Signal: Strength vs. Dispersion",
      x = expression(-log[10]~"(Total Influence Strength)"),
      y = "Influence Coefficient of Variation (CV)",
      size = "Prevalence"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  return(p)
}

OTU_to_Taxonomy <- function(physeq, otu_ids) {
  tax_tab <- as(tax_table(physeq),"matrix")
  tax_tab <- as.data.frame(tax_tab)
  tax_subset <- tax_tab[rownames(tax_tab) %in% otu_ids, , drop = FALSE]
  tax_subset <- tibble::tibble(OTU = rownames(tax_subset), tax_subset)
  return(tax_subset)
}

compute_enrichment_vs_prevalence <- function(physeq, weighted_Rtree, unweighted_Rtree,group1, group2,weighted_gamma,unweighted_gamma) {
  library(phyloseq)
  library(ggplot2)
  
  # OTU matrix
  otu <- as.matrix(otu_table(physeq))
  if (!taxa_are_rows(physeq)) otu <- t(otu)
  #prevalence <- rowSums(otu > 0) / ncol(otu)
  
  # Enrichment Score: group1 - group2
  if (!(group1 %in% colnames(weighted_Rtree)) || !(group2 %in% colnames(weighted_Rtree))) {
    stop("Group names not found in Rtree columns.")
  }
  
  enrichment_score <- weighted_Rtree[, group1] - weighted_Rtree[, group2]
  occupation_score <- unweighted_Rtree[, group1] - unweighted_Rtree[, group2]
  df <- data.frame(
    OTU = rownames(weighted_Rtree),
    OccupationScore = log10(abs(occupation_score))*ifelse(occupation_score>=0,1,-1),
    EnrichmentScore = log10(abs(enrichment_score))*ifelse(enrichment_score>=0,1,-1)
  )
  
  # Plot
  p <- ggplot(df, aes(x = df$OccupationScore, y = df$EnrichmentScore)) +
    geom_point(alpha = 0.6, color = "#003366") +
    annotate("rect",
             xmin = quantile(df$OccupationScore,unweighted_gamma),
             xmax = quantile(df$OccupationScore,1-unweighted_gamma),
             ymin = quantile(df$EnrichmentScore,weighted_gamma),
             ymax = quantile(df$EnrichmentScore,1-weighted_gamma),
             color = "red",
             fill = NA,
             linetype = "dashed",
             size = 1) +
    theme_minimal(base_size = 14) +
    labs(
      title = "OccupationScore vs. Enrichment Score",
      x = "OccupationScore",
      y = paste("Enrichment Score:", group1, "–", group2)
    ) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  return(list(
    table = df,
    plot = p
  ))
}


compute_enrichment_vs_prevalence2 <- function(
    physeq, 
    weighted_Rtree, 
    unweighted_Rtree,
    groupvariable,
    group1, 
    lower_thresh = 0.05, 
    upper_thresh = 0.05
) {
  library(phyloseq)
  library(ggplot2)
  
  # OTU matrix
  otu <- as(otu_table(physeq),"matrix")
  #otu <- as.matrix(otu)
  if (!taxa_are_rows(physeq)) otu <- t(otu)
  
  meta <- as(sample_data(physeq), "data.frame")
  
  # 找到满足条件的样本名
  group_samples <- rownames(meta[meta[[groupvariable]] == group1, ])
  
  # 过滤 OTU 表
  otu_sub <- otu[, group_samples, drop = FALSE]
  
  # 计算 prevalence（每个 OTU 在组内出现的样本数 / 样本总数）
  prevalence <- rowSums(otu_sub > 0) / ncol(otu_sub)
  
  otu_tss <- sweep(otu_sub, 2, colSums(otu_sub), FUN = "/")
  
  # 计算每个 OTU 在该组内的平均相对丰度
  mean_abundance <- rowMeans(otu_tss, na.rm = TRUE)
  
  # Normalize Rtrees
  if (!(group1 %in% colnames(weighted_Rtree)) ) {
    stop("Group name not found in Rtree columns.")
  }
  weighted_Rtree <- scale(weighted_Rtree, center = TRUE, scale = TRUE)
  unweighted_Rtree <- scale(unweighted_Rtree, center = TRUE, scale = TRUE)
  
  p_weighted <- pnorm(weighted_Rtree, lower.tail = FALSE)
  p_unweighted <- pnorm(unweighted_Rtree, lower.tail = FALSE)
  
  #p_adj_weighted <- apply(p_weighted, 2, p.adjust, method = "BH")
  #p_adj_unweighted <- apply(p_unweighted, 2, p.adjust, method = "BH")
  
  p_weighted_lower <- pnorm(weighted_Rtree, lower.tail = TRUE)
  p_unweighted_lower <- pnorm(unweighted_Rtree, lower.tail = TRUE)
  
  #p_adj_weighted_lower <- apply(p_weighted_lower, 2, p.adjust, method = "BH")
  #p_adj_unweighted_lower <- apply(p_unweighted_lower, 2, p.adjust, method = "BH")


  # 分类标签
  category <- ifelse(
    p_weighted[,group1] < upper_thresh & p_unweighted[,group1] < upper_thresh, "Core",
    ifelse(p_weighted_lower[,group1] < lower_thresh & p_unweighted_lower[,group1] < lower_thresh, "Transient", "Intermediate")
  )
  
  df <- data.frame(
    OTU = rownames(weighted_Rtree),
    Unweighted_upper = p_unweighted[,group1],
    Unweighted_lower =p_unweighted_lower[,group1],
    Weighted_upper = p_weighted[,group1],
    Weighted_lower = p_weighted_lower[,group1],
    Prevalence = prevalence,
    Category = factor(category, levels = c("Core", "Intermediate", "Transient")),
    Unweighted_z = unweighted_Rtree[,group1],
    Weighted_z = weighted_Rtree[,group1],
    Mean_abundance = mean_abundance
  )
  color_map <- c("Core" = "#D7263D", "Intermediate" = "gray70", "Transient" = "#1E90FF")
  shape_map <- c("Core" = 17, "Intermediate" = 16, "Transient" = 15)
  p1 <- ggplot(df, aes(x = Unweighted_z, y = Weighted_z, color = Category, shape = Category)) +
    geom_point(alpha = 0.8, size = 2.2) +
    scale_color_manual(values = color_map) +
    scale_shape_manual(values = shape_map) +
    theme_classic(base_size = 14) +
    labs(
      x = "Coexist Z Score",
      y = "Coenrichment Z Score"
    ) +
    theme(plot.title = element_text(hjust = 0.5))
  
  # 图2: Prevalence vs Coenrichment
  p2 <- ggplot(df, aes(x = Prevalence, y = Weighted_z)) +
    geom_point(color = "gray30", alpha = 0.7, size = 2)+
    scale_color_manual(values = color_map) +
    theme_classic(base_size = 14) +
    labs(
      x = "Prevalence",
      y = "Coenrichment Z Score"
    )
  
  # 图3: Mean abundance vs Coenrichment
  p3 <- ggplot(df, aes(x = Mean_abundance, y = Weighted_z)) +
    geom_point(color = "gray30", alpha = 0.7, size = 2)+
    scale_color_manual(values = color_map) +
    theme_classic(base_size = 14) +
    labs(
      x = "Mean Relative Abundance",
      y = "Coenrichment Z Score"
    )
  p4 <- ggplot(df, aes(x = Prevalence, y = Mean_abundance, color = Category, shape = Category)) +
    geom_point(alpha = 0.8, size = 2.2) +
    scale_color_manual(values = color_map) +
    scale_shape_manual(values = shape_map) +
    theme_classic(base_size = 14) +
    geom_vline(xintercept = 0.6, color = "gray50", linetype = "dashed", linewidth = 0.5) +
    geom_hline(yintercept = 0.01, color = "gray50", linetype = "dashed", linewidth = 0.5)+
    labs(
      x = "Prevalence",
      y = "Mean Relative Abundance"
    )
  return(list(
    table = df,
    plot1 = p1,
    plot2 = p2,
    plot3 = p3,
    plot4 = p4
  ))
}

oa_plot<-function(rel_abund,occupancy_high,occupancy_low,threshold){
  
  otu_rel <- as(otu_table(rel_abund), "matrix")
  if (!taxa_are_rows(rel_abund)) {
    otu_rel <- t(otu_rel)
  }
  occupancy <- rowMeans(otu_rel > 0)
  mean_abundance <- rowMeans(otu_rel)
  
  occ_abund_df <- data.frame(
    Taxon = rownames(otu_rel),
    Occupancy = occupancy,
    MeanAbundance = mean_abundance
  )
  
  occ_abund_df$Category <- case_when(
    occupancy >= occupancy_high & mean_abundance >= threshold ~ "Core",
    occupancy <= occupancy_low & mean_abundance < threshold ~ "Transient",
    TRUE ~ "Intermediate"
  )
  
  category_colors <- c("Core" = "#0072B2", "Transient" = "#D55E00", "Intermediate" = "gray60")
  
  ggplot(occ_abund_df, aes(x = Occupancy, y = MeanAbundance, color = Category)) +
    geom_point(alpha = 0.8, size = 2.2) +
    scale_y_log10(labels = scales::scientific, breaks = scales::log_breaks()) +
    scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.2)) +
    scale_color_manual(values = category_colors) +
    labs(
      title = "Occupancy-Abundance Plot",
      subtitle = "Core taxa: Occ ≥ 0.6 & Abund ≥ 0.1%",
      x = "Occupancy (Proportion of Samples Present)",
      y = "Mean Relative Abundance (log scale)",
      color = "Taxon Category"
    ) +
    geom_hline(yintercept = threshold, linetype = "dashed", color = "black") +
    geom_vline(xintercept = occupancy_high, linetype = "dashed", color = "black") +
    geom_vline(xintercept = occupancy_low, linetype = "dashed", color = "black") +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.subtitle = element_markdown(size = 12, hjust = 0.5),
      legend.position = "right"
    )
}

get_tax_label <- function(tax_table_subset, levels = c("Genus", "Family", "Order", "Class", "Phylum", "Kingdom")) {
  label_vec <- character(nrow(tax_table_subset))
  
  for (i in seq_len(nrow(tax_table_subset))) {
    for (level in levels) {
      val <- as.character(tax_table_subset[i, level])
      if (!is.na(val) && val != "") {
        label_vec[i] <- val
        break
      }
    }
  }
  
  return(label_vec)
}



check_taxa_presence_matrix <- function(weighted, unweighted, level = "Genus") {
  
  HOMD <- read.csv("./decontam reference/HOMD_taxon.xls.csv")
  Saltercontam <- read.csv("./decontam reference/contamination_database.csv")# 清洗输入向量
  taxa_w <- tolower(trimws(weighted[!is.na(weighted)]))
  taxa_u <- tolower(trimws(unweighted[!is.na(unweighted)]))
  
  # 合并唯一分类名
  all_taxa <- sort(unique(c(taxa_w, taxa_u)))
  
  # 处理 HOMD：flatten + 标准化
  ref_homd <- tolower(trimws(unlist(HOMD)))
  ref_homd <- ref_homd[!is.na(ref_homd)]
  
  # 处理 Saltercontam：确保是字符向量
  if (is.data.frame(Saltercontam)) {
    ref_salter <- tolower(trimws(unlist(Saltercontam)))
  } else {
    ref_salter <- tolower(trimws(Saltercontam))
  }
  ref_salter <- ref_salter[!is.na(ref_salter)]
  
  # 构建结果表格
  df <- data.frame(
    Taxon = all_taxa,
    Weighted = as.integer(all_taxa %in% taxa_w),
    Unweighted = as.integer(all_taxa %in% taxa_u),
    HOMD = as.integer(all_taxa %in% ref_homd),
    Saltercontam = as.integer(all_taxa %in% ref_salter),
    stringsAsFactors = FALSE
  )
  
  return(df)
}


compute_enrichment_vs_prevalence_ibd <- function(
    physeq,
    weighted_Rtree,
    unweighted_Rtree,
    groupvariable,
    group1,
    lower_thresh = 0.05,
    upper_thresh = 0.05,
    vline = 0.2,
    hline = 0.01
) {
  library(phyloseq)
  library(ggplot2)
  
  # ----------------------------
  # 0) 基础检查 + OTU 对齐
  # ----------------------------
  stopifnot(all(dim(weighted_Rtree) == dim(unweighted_Rtree)))
  stopifnot(all(rownames(weighted_Rtree) == rownames(unweighted_Rtree)))
  stopifnot(all(colnames(weighted_Rtree) == colnames(unweighted_Rtree)))
  
  if (!(group1 %in% colnames(weighted_Rtree))) {
    stop("Group name not found in Rtree columns: ", group1)
  }
  
  # physeq OTU matrix
  otu <- as(otu_table(physeq), "matrix")
  if (!taxa_are_rows(physeq)) otu <- t(otu)
  
  meta <- as(sample_data(physeq), "data.frame")
  if (!(groupvariable %in% colnames(meta))) {
    stop("groupvariable not found in sample_data: ", groupvariable)
  }
  
  group_samples <- rownames(meta[meta[[groupvariable]] == group1, , drop = FALSE])
  if (length(group_samples) == 0) {
    stop("No samples found for group1 = ", group1, " in ", groupvariable)
  }
  
  # 强制以共同 OTU 为准（避免任何命名/顺序问题）
  common_otus <- intersect(rownames(weighted_Rtree), rownames(otu))
  if (length(common_otus) == 0) stop("No overlapping OTUs between physeq and Rtree.")
  
  weighted_Rtree <- weighted_Rtree[common_otus, , drop = FALSE]
  unweighted_Rtree <- unweighted_Rtree[common_otus, , drop = FALSE]
  otu <- otu[common_otus, , drop = FALSE]
  
  # group 内子表
  otu_sub <- otu[, group_samples, drop = FALSE]
  
  # ----------------------------
  # 1) prevalence + mean abundance（组内）
  # ----------------------------
  prevalence <- rowSums(otu_sub > 0) / ncol(otu_sub)
  
  denom <- colSums(otu_sub)
  denom[denom == 0] <- NA  # 防止 0/0
  otu_tss <- sweep(otu_sub, 2, denom, FUN = "/")
  mean_abundance <- rowMeans(otu_tss, na.rm = TRUE)
  
  # 只保留在该组内出现过的 OTU（IBD blank 很重要）
  keep <- is.finite(prevalence) & is.finite(mean_abundance) & (prevalence > 0)
  # 如果你希望更宽松（只要 prevalence>0 不管 mean_abundance），可把 mean_abundance 条件去掉
  prevalence <- prevalence[keep]
  mean_abundance <- mean_abundance[keep]
  
  # 同步裁剪 Rtree
  keep_otus <- names(prevalence)
  weighted_Rtree <- weighted_Rtree[keep_otus, , drop = FALSE]
  unweighted_Rtree <- unweighted_Rtree[keep_otus, , drop = FALSE]
  
  if (length(keep_otus) == 0) {
    stop("After filtering prevalence>0, no OTUs remain for plotting. (Group may be too sparse.)")
  }
  
  # ----------------------------
  # 2) Rtree 标准化 + p 值 + 分类
  # ----------------------------
  weighted_z_all <- scale(weighted_Rtree, center = TRUE, scale = TRUE)
  unweighted_z_all <- scale(unweighted_Rtree, center = TRUE, scale = TRUE)
  
  p_weighted_upper <- pnorm(weighted_z_all, lower.tail = FALSE)
  p_unweighted_upper <- pnorm(unweighted_z_all, lower.tail = FALSE)
  p_weighted_lower <- pnorm(weighted_z_all, lower.tail = TRUE)
  p_unweighted_lower <- pnorm(unweighted_z_all, lower.tail = TRUE)
  
  category <- ifelse(
    p_weighted_upper[, group1] < upper_thresh & p_unweighted_upper[, group1] < upper_thresh, "High risk",
    ifelse(p_weighted_lower[, group1] < lower_thresh & p_unweighted_lower[, group1] < lower_thresh, "Low risk", "Medium risk")
  )
  names(category) <- rownames(weighted_z_all)
  
  df <- data.frame(
    OTU = rownames(weighted_z_all),
    Unweighted_upper = p_unweighted_upper[, group1],
    Unweighted_lower = p_unweighted_lower[, group1],
    Weighted_upper = p_weighted_upper[, group1],
    Weighted_lower = p_weighted_lower[, group1],
    Prevalence = prevalence[rownames(weighted_z_all)],
    Mean_abundance = mean_abundance[rownames(weighted_z_all)],
    Unweighted_z = unweighted_z_all[, group1],
    Weighted_z = weighted_z_all[, group1],
    Category = factor(category[rownames(weighted_z_all)], levels = c("High risk", "Low risk", "Medium risk")),
    row.names = rownames(weighted_z_all)
  )
  
  # 最后再保险一次：去掉任何 NA
  df <- df[is.finite(df$Prevalence) & is.finite(df$Mean_abundance) & !is.na(df$Category), , drop = FALSE]
  if (nrow(df) == 0) stop("df became empty after final NA filtering; please inspect prevalence/mean_abundance.")
  
  # ----------------------------
  # 3) 作图（与 tumor 的风格一致）
  # ----------------------------
  color_map <- c("High risk" = "#D7263D", "Low risk" = "#1E90FF", "Medium risk" = "gray70")
  shape_map <- c("High risk" = 17, "Low risk" = 15, "Medium risk" = 16)
  
  p1 <- ggplot(df, aes(x = Unweighted_z, y = Weighted_z, color = Category, shape = Category)) +
    geom_point(alpha = 0.8, size = 2.2) +
    scale_color_manual(values = color_map, drop = FALSE) +
    scale_shape_manual(values = shape_map, drop = FALSE) +
    theme_classic(base_size = 14) +
    labs(
      title = paste("Group:", group1),
      x = "Coexist Z Score",
      y = "Coenrichment Z Score"
    ) +
    theme(plot.title = element_text(hjust = 0.5))
  
  p2 <- ggplot(df, aes(x = Prevalence, y = Weighted_z)) +
    geom_point(color = "gray30", alpha = 0.7, size = 2) +
    theme_classic(base_size = 14) +
    labs(x = "Prevalence", y = "Coenrichment Z Score")
  
  p3 <- ggplot(df, aes(x = Mean_abundance, y = Weighted_z)) +
    geom_point(color = "gray30", alpha = 0.7, size = 2) +
    theme_classic(base_size = 14) +
    labs(x = "Mean Relative Abundance", y = "Coenrichment Z Score")
  
  p4 <- ggplot(df, aes(x = Prevalence, y = Mean_abundance, color = Category, shape = Category)) +
    geom_point(alpha = 0.8, size = 2.2) +
    scale_color_manual(values = color_map, drop = FALSE) +
    scale_shape_manual(values = shape_map, drop = FALSE) +
    theme_classic(base_size = 14) +
    geom_vline(xintercept = vline, color = "gray50", linetype = "dashed", linewidth = 0.5) +
    geom_hline(yintercept = hline, color = "gray50", linetype = "dashed", linewidth = 0.5) +
    labs(x = "Prevalence", y = "Mean Relative Abundance")
  
  return(list(
    table = df,
    plot1 = p1,
    plot2 = p2,
    plot3 = p3,
    plot4 = p4
  ))
}
