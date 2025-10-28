# ICB Response Explorer

A Shiny web application for exploring **immunotherapy response** at both **patient-level (pseudobulk)** and **cell-level (single-cell DE analysis)** across multiple datasets, including **PAN-T** and **CAR-T** studies.

---

## Overview

**ICB Response Explorer** allows users to interactively visualize and compare **Responder vs Non-responder** differences in gene expression across studies and cell types.

The app supports:
- **Pan-T Patient-level**  
- **Pan-T Cell-level**
- **CAR-T patient-level**
- **CAR-T cell-level** 
- (Planned) **CRISPR** 

---

## Key Features

| Feature | Description |
|----------|-------------|
| **Cancer Type Coverage** | Includes samples from **Skin (melanoma)**, **Liver (HCC)**, **Breast**, and **Kidney (RCC)** cancers. |
| **Dynamic Gene Index** | Builds searchable gene lists on-demand for each dataset (PAN-T and CAR-T). |
| **Interactive Visualizations** | Generates heatmaps, boxplots, CPM bar plots, dot plots, and lollipop plots. |
| **Download Support** | Allows exporting filtered DE results as `.tsv` tables. |

---

## Dataset Source

Pan-T Datasets are derived from the Nature Medicine (2023) publication:  
> **"A single-cell immune atlas of solid tumors treated with immune checkpoint blockade"**  
> [https://www.nature.com/articles/s41591-023-02371-y#data-availability](https://www.nature.com/articles/s41591-023-02371-y#data-availability)

Data were accessed from the GEO/SCP repositories linked in the publication and reformatted for Shiny-based visualization.

---

## Run Locally

### Prerequisites

- R (≥ 4.3)
- RStudio (recommended)
- Required packages:
  ```r
  install.packages(c(
    "shiny", "shinyjqui", "tidyverse", "data.table",
    "DT", "scales", "ggpubr", "memoise", "cachem"
  ))
  library(shiny)
  runApp("path/to/ICB_Response_Explorer")
