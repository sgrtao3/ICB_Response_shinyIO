# ICB Response Explorer

A Shiny web application for exploring **immunotherapy response** at both **patient-level (pseudobulk)** and **cell-level (single-cell DE analysis)** across multiple datasets, including **PAN-T** and **CAR-T** studies.

---

## Overview

**ICB Response Explorer** allows users to interactively visualize and compare **Responder vs Non-responder** differences in gene expression across studies and cell types.

The app supports:
- **Patient-level (pseudobulk)** analysis with boxplots and heatmaps  
- **Cell-level** analysis from precomputed DE tables with filtering, visualization, and export  
- **CAR-T patient-level** and **cell-level** modules  
- (Planned) **CRISPR** analysis module  

---

## Key Features

| Feature | Description |
|----------|-------------|
| **Lazy Loading** | Large datasets (e.g. expression matrices) are only loaded when needed, preventing crashes on limited-memory servers. |
| **Caching (`cachem`)** | Uses `cachem::cache_mem` to store computed results in memory for faster re-renders. |
| **Dynamic Gene Index** | Builds searchable gene lists on-demand for each dataset (PAN-T and CAR-T). |
| **Interactive Visualizations** | Generates heatmaps, boxplots, CPM bar plots, dot plots, and lollipop plots. |
| **Download Support** | Allows exporting filtered DE results as `.tsv` tables. |

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
