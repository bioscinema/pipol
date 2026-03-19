# PIPoL **Phylogeny-Informed Propagation of Labels for Microbiome Decontamination**

PIPoL is an R package for identifying potential contaminant taxa in microbiome sequencing data using **tree-aware label propagation**.  
The method integrates **co-existence (presence/absence)** and **co-enrichment (abundance-weighted)** signals on a **phylogenetic tree** to distinguish likely contaminants from true biological taxa.

---

## Key Features

- **Tree-aware propagation** using phylogenetic distances  
- **Dual modeling of contamination signals**
  - Co-existence (unweighted, presence/absence)
  - Co-enrichment (weighted, abundance-based)
- **Explicit support for negative controls** (e.g. blanks)
- **Integrated visualization** for enrichment, prevalence, and risk stratification
- Designed for **real microbiome datasets** (IBD, oral, tumor, etc.)

---

## Installation

Install the development version from GitHub:

```r
remotes::install_github("bioscinema/PIPoL")
```

Load the package:

```r
library(PIPoL)
```

---

## Conceptual Overview

PIPoL operates in three main steps:

1. Estimate a smoothing scale (tau) from the phylogenetic tree

2. Propagate sample-level labels across taxa on the tree

3. Integrate co-existence and co-enrichment signals to assess contamination risk

Pipeline overview:

```scss
phyloseq object
      ↓
   BTS (tau)
      ↓
PIPoL_tree (weighted / unweighted)
      ↓
decontamination
      ↓
risk scores + visualization
```

---

## Basic Workflow

### 1. Prepare a phyloseq object

Your phyloseq object should contain:

- OTU / ASV table

- Taxonomy table

- Sample metadata with a group variable (e.g. "blank", "case", "control")

```r
physeq
```

---

### 2. Estimate tree smoothing parameter

```r
tau_lib <- BTS(physeq)
```

---

### 3. Run tree-aware label propagation (is embedded in each functions)

Abundance-weighted (co-enrichment):

```r
Rtree_weighted <- PIPoL_tree(
  physeq,
  outcome_var = "GroupLabel",
  tau = tau_lib,
  weighted = TRUE
)
```

Presence/absence (co-existence):

```r
Rtree_unweighted <- PIPoL_tree(
  physeq,
  outcome_var = "GroupLabel",
  tau = tau_lib,
  weighted = FALSE
)
```

---

### 4. Identify contaminant taxa

```r
res <- decontamination(
    physeq,
    outcome_var,
    label,
    tau = c("Genus", "Family", "Species"),
    alpha = 0.05,
    gamma_weighted = 0.05,
    gamma_unweighted = 0.05,
    confounders = NULL
)
```

---

## Output

The decontamination() function returns a list containing:

```r
names(res)
```

- enrich_list
Taxa enriched in the target label (abundance-based)

- exist_list
Taxa enriched in the target label (presence-based)

- enrichment_table
Taxonomy table for enrichment hits

- existence_table
Taxonomy table for co-existence hits

- scatter_plot
Joint co-existence / co-enrichment visualization

Display the main plot:

```r
res$scatter_plot
```

---

## Enrichment–Prevalence Analysis

### Core function (method-level)

```r
enrichment_prevalence(
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
)
```

This function provides a generic, method-level characterization of how PIPoL enrichment signals relate to:

- prevalence

- mean relative abundance

- co-existence versus co-enrichment structure

---


## Notes on Interpretation

- High-risk taxa show strong co-existence and strong co-enrichment with negative controls

- Low-risk taxa show depletion in both signals

- Intermediate taxa require further inspection

Thresholds (gamma, prevalence cutoffs) should be adjusted based on study design, sequencing depth, and availability of negative controls.

