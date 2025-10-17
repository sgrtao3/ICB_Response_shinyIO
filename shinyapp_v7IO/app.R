## ===================== app.R — Combined Explorer (low-RAM) =====================
## Patient-level (pseudobulk) + Cell-level (DE) with on-demand file reading
## ============================================================================

## ----- Packages -----
library(shiny)
library(shinyjqui)
library(tidyverse)   # dplyr, ggplot2, forcats, readr, tibble, tidyr
library(data.table)
library(DT)
library(scales)
library(ggpubr)

options(shiny.sanitize.errors = FALSE)
`%||%` <- function(a,b) if (is.null(a)) b else a

## =========================
## = Patient-level helpers =
## =========================

# Expression loader (RDS or TSV) with duplicate column collapse
load_expr <- function(path) {
  if (!file.exists(path)) stop("Expression file not found: ", path)
  if (grepl("\\.rds$", path, ignore.case = TRUE)) {
    message("Reading RDS file: ", path)
    return(readRDS(path))
  }
  message("Reading TSV file: ", path)
  expr_dt <- data.table::fread(path, sep = "\t", header = TRUE, showProgress = FALSE)
  data.table::setnames(expr_dt, old = names(expr_dt)[1], new = "gene")
  
  fix_ids <- function(x) {
    x %>% stringr::str_trim() %>%
      stringr::str_replace_all('^[\'"]|[\'"]$', '') %>%
      stringr::str_replace_all("`", "") %>%
      stringr::str_replace("(\\.x|\\.y)$", "") %>%
      stringr::str_replace("\\.\\d+$", "")
  }
  colnames(expr_dt) <- c("gene", fix_ids(colnames(expr_dt)[-1]))
  
  long_dt <- data.table::melt(
    data.table::as.data.table(expr_dt),
    id.vars = "gene", variable.name = "ID", value.name = "value"
  )
  long_dt <- long_dt[, .(value = mean(value, na.rm = TRUE)), by = .(gene, ID)]
  wide_dt <- data.table::dcast(long_dt, gene ~ ID, value.var = "value")
  
  as.data.frame(wide_dt) %>%
    mutate(
      gene = as.character(gene),
      gene = stringr::str_replace_all(gene, '^[\'"]|[\'"]$', ''),
      gene = stringr::str_trim(gene),
      gene = dplyr::na_if(gene, "")
    ) %>%
    filter(!is.na(gene)) %>%
    group_by(gene) %>%
    summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
    tibble::column_to_rownames("gene")
}

normalize_id <- function(x) {
  x %>% stringr::str_trim() %>%
    stringr::str_replace_all('^[\'"]|[\'"]$', '') %>%
    stringr::str_replace_all("`", "") %>%
    stringr::str_replace("(\\.x|\\.y)$", "") %>%
    stringr::str_replace("\\.\\d+$", "")
}

filter_meta <- function(meta_df, prepost, celltypes) {
  meta_f <- meta_df %>%
    dplyr::filter(is.na(response_group) | response_group %in% c("Responder","Non-responder"))
  
  if (!is.null(prepost) && prepost %in% c("Pre","Post") && "pre_post" %in% names(meta_f)) {
    meta_f <- meta_f %>% dplyr::filter(pre_post %in% prepost)
  }
  if (!is.null(celltypes) && length(celltypes) > 0 && "cell_type" %in% names(meta_f)) {
    meta_f <- meta_f %>% dplyr::filter(cell_type %in% celltypes)
  }
  meta_f %>% dplyr::select(dplyr::any_of(c("ID","donor_id","cell_type","pre_post","response_group","study")))
}

safe_wilcox <- function(df) {
  df <- df %>% dplyr::filter(response_group %in% c("Responder","Non-responder"))
  n_R  <- sum(df$response_group == "Responder", na.rm = TRUE)
  n_NR <- sum(df$response_group == "Non-responder", na.rm = TRUE)
  if (isTRUE(n_R > 0 & n_NR > 0)) {
    tryCatch(wilcox.test(expr ~ response_group, data = df, exact = FALSE)$p.value,
             error = function(e) NA_real_)
  } else NA_real_
}

compute_one_study <- function(study_label, gene, prepost, celltypes, meta_all_df, expr_all_mat) {
  mf <- meta_all_df %>% dplyr::filter(study == study_label)
  if (nrow(mf) == 0) return(NULL)
  
  mf <- mf %>% dplyr::filter(response_group %in% c("Responder","Non-responder"))
  if (!is.null(prepost) && prepost %in% c("Pre","Post") && "pre_post" %in% names(mf)) {
    mf <- mf %>% dplyr::filter(pre_post %in% prepost)
  }
  if (!is.null(celltypes) && length(celltypes) > 0 && "cell_type" %in% names(mf)) {
    mf <- mf %>% dplyr::filter(cell_type %in% celltypes)
  }
  
  common_ids <- intersect(colnames(expr_all_mat), mf$ID)
  if (length(common_ids) == 0) return(NULL)
  if (!(gene %in% rownames(expr_all_mat))) return(NULL)
  
  sample_df <- tibble::tibble(
    ID   = common_ids,
    expr = as.numeric(expr_all_mat[gene, common_ids, drop = TRUE])
  ) %>% dplyr::left_join(mf, by = "ID") %>%
    dplyr::filter(!is.na(expr)) %>%
    dplyr::mutate(response_group = factor(response_group, levels = c("Non-responder","Responder")))
  
  if (nrow(sample_df) == 0) return(NULL)
  
  patient_df <- sample_df %>%
    dplyr::group_by(cell_type, donor_id, response_group) %>%
    dplyr::summarise(expr = mean(expr, na.rm = TRUE), .groups = "drop")
  
  stats_tbl <- patient_df %>%
    dplyr::group_split(cell_type) %>%
    purrr::map_dfr(function(df) {
      df <- df %>% dplyr::filter(response_group %in% c("Responder","Non-responder"))
      n_R  <- sum(df$response_group == "Responder", na.rm = TRUE)
      n_NR <- sum(df$response_group == "Non-responder", na.rm = TRUE)
      medR  <- suppressWarnings(stats::median(df$expr[df$response_group == "Responder"], na.rm = TRUE))
      medNR <- suppressWarnings(stats::median(df$expr[df$response_group == "Non-responder"], na.rm = TRUE))
      dir <- dplyr::case_when(
        is.na(medR) | is.na(medNR) ~ 0L,
        medR == medNR              ~ 0L,
        medR >  medNR              ~ 1L,
        TRUE                       ~ -1L
      )
      p <- safe_wilcox(df)
      tibble::tibble(
        Study        = study_label,
        cell_type    = df$cell_type[1],
        n_R          = n_R,
        n_NR         = n_NR,
        median_R     = medR,
        median_NR    = medNR,
        direction    = dir,
        p_value      = p,
        signed_score = dplyr::if_else(is.na(p), NA_real_, -log10(p) * dir)
      )
    })
  
  list(patient = patient_df %>% dplyr::mutate(Study = study_label),
       stats   = stats_tbl)
}

## =======================
## = Cell-level helpers  =
## =======================

# Resolve path (tries both raw and data/<dir>)
resolve_dir <- function(dir_path) {
  candidates <- c(dir_path, file.path("data", dir_path))
  existing <- candidates[dir.exists(candidates)]
  if (length(existing) == 0) return(NULL)
  existing[1]
}

list_de_files <- function(dir_path) {
  if (is.null(dir_path)) return(character(0))
  list.files(dir_path,
             pattern = "\\.(tsv|txt|tsv\\.gz|txt\\.gz)$",
             full.names = TRUE, ignore.case = TRUE, recursive = TRUE)
}

# Read *only rows for one gene* from a single file (keeps RAM tiny)
read_gene_from_file <- function(fp, gene_sel, use_symbol, KEEP_COLS) {
  cols_to_read <- intersect(KEEP_COLS, c(
    "Study","cluster","gene","hgnc_symbol","gene_biotype",
    "p_val","avg_log2FC","pct.1","pct.2","p_val_adj",
    "Responder_CPM","Non_responder_CPM"
  ))
  dt <- tryCatch(
    fread(fp, sep = "\t", header = TRUE, select = cols_to_read,
          data.table = TRUE, showProgress = FALSE),
    error = function(e) NULL
  )
  if (is.null(dt) || !nrow(dt)) return(NULL)
  
  # Harmonize pct column names if needed
  if (!("pct.1" %in% names(dt)) && ("pct_1" %in% names(dt))) data.table::setnames(dt, "pct_1", "pct.1")
  if (!("pct.2" %in% names(dt)) && ("pct_2" %in% names(dt))) data.table::setnames(dt, "pct_2", "pct.2")
  
  # Filter by gene/symbol
  if (use_symbol && "hgnc_symbol" %in% names(dt)) {
    dt <- dt[tolower(hgnc_symbol) == tolower(gene_sel) | tolower(gene) == tolower(gene_sel)]
  } else {
    dt <- dt[tolower(gene) == tolower(gene_sel) | ("hgnc_symbol" %in% names(dt) &&
                                                     tolower(hgnc_symbol) == tolower(gene_sel))]
  }
  if (!nrow(dt)) return(NULL)
  
  # Coerce numerics (prevents character-vs-numeric comparisons)
  num_cols <- intersect(
    c("p_val","avg_log2FC","p_val_adj","Responder_CPM","Non_responder_CPM","pct.1","pct.2"),
    names(dt)
  )
  for (cc in num_cols) dt[[cc]] <- suppressWarnings(as.numeric(dt[[cc]]))
  
  tibble::as_tibble(dt)
}


# Build a small index of unique gene IDs for the dropdown (reads just 2 cols)
build_gene_index <- function(files) {
  if (!length(files)) return(character(0))
  uniq <- new.env(parent = emptyenv())
  add_vec <- function(v) { for (z in v[!is.na(v) & nzchar(v)]) uniq[[tolower(z)]] <- TRUE }
  for (fp in files) {
    dt <- tryCatch(
      fread(fp, sep = "\t", header = TRUE, select = c("gene","hgnc_symbol"),
            data.table = TRUE, showProgress = FALSE, nThread = 1L),
      error = function(e) NULL
    )
    if (is.null(dt)) next
    if (!("hgnc_symbol" %in% names(dt))) dt[, hgnc_symbol := NA_character_]
    add_vec(dt$gene); add_vec(dt$hgnc_symbol)
    # Free as we go
    rm(dt); gc(FALSE)
  }
  sort(names(as.list(uniq)))
}

## =========================
## = CONFIG (relative)    =
## =========================
p_meta_path <- file.path("data", "combined_meta_using_Combined_outcome.tsv")
p_expr_path <- file.path("data", "combined_expr_using_Combined_outcome.rds")

dir_pre  <- resolve_dir("DE_per_study_pre")
dir_post <- resolve_dir("DE_per_study_post")

KEEP_COLS <- c("Study","cluster","hgnc_symbol","gene","gene_biotype",
               "p_val","avg_log2FC","pct.1","pct.2","p_val_adj",
               "Responder_CPM","Non_responder_CPM")

validate <- shiny::validate; need <- shiny::need

## =========== Load patient-level data ===========
if (!file.exists(p_meta_path))
  stop("Missing file: ", p_meta_path, "\nBundle it under the app's data/ folder.")
if (!file.exists(p_expr_path))
  stop("Missing file: ", p_expr_path, "\nBundle it under the app's data/ folder.")

meta_all <- readr::read_tsv(p_meta_path, show_col_types = FALSE) %>%
  mutate(
    across(c(ID, donor_id, pre_post, cell_type, study), ~ as.character(.x)),
    ID        = normalize_id(ID),
    donor_id  = ifelse(is.na(donor_id) | donor_id == "", ID, donor_id),
    Combined_outcome = trimws(Combined_outcome),
    response_group = case_when(
      Combined_outcome %in% c("Favourable")   ~ "Responder",
      Combined_outcome %in% c("Unfavourable") ~ "Non-responder",
      Combined_outcome %in% c("UT","n/a", NA) ~ NA_character_,
      TRUE                                     ~ NA_character_
    )
  )

study_keep <- meta_all %>%
  dplyr::group_by(study) %>%
  dplyr::summarise(any_labeled = any(!is.na(response_group)), .groups = "drop") %>%
  dplyr::filter(any_labeled) %>% dplyr::pull(study)
meta_all <- meta_all %>% dplyr::filter(study %in% study_keep)

expr_mat_all     <- load_expr(p_expr_path)
p_all_genes      <- sort(rownames(expr_mat_all))
p_all_celltypes  <- sort(unique(meta_all$cell_type))

## =========== Prepare cell-level file lists + small gene index ===========
files_pre  <- list_de_files(dir_pre)
files_post <- list_de_files(dir_post)
all_files  <- unique(c(files_pre, files_post))

# Small dropdown index (only names, not whole tables)
gene_index <- build_gene_index(all_files)

## ======================
## =         UI         =
## ======================

ui <- navbarPage(
  title = "ICB Response Explorer",
  id = "top_tabs",
  
  ## ---- Patient-level tab ----
  tabPanel(
    "Patient-level",
    sidebarLayout(
      sidebarPanel(
        width = 3,
        h4("Responder vs Non-responder (Patient-level)"),
        selectizeInput(
          "p_gene", "Gene:", choices = p_all_genes,
          selected = "A1BG",
          options = list(placeholder = "Type a gene symbol…", create = FALSE)
        ),
        radioButtons(
          "p_prepost", "Timepoint:",
          choices = c("Pre","Post","Both"), selected = "Post", inline = TRUE
        ),
        selectizeInput(
          "p_celltypes", "Cell types (optional):",
          choices = p_all_celltypes, multiple = TRUE
        ),
        uiOutput("p_plot_study_picker"),
        actionButton("p_run", "Update", class = "btn-primary")
      ),
      mainPanel(
        tabsetPanel(
          id = "p_subtabs",
          tabPanel(
            "Heatmaps",
            br(),
            fluidRow(
              column(
                12, h4("Score"),
                jqui_resizable(
                  div(id = "p_heat_signed_box", style = "height:600px;",
                      plotOutput("p_heat_signed", height = "100%")),
                  options = list(handles = "se")
                )
              ),
              column(
                12, h4("P-value"),
                jqui_resizable(
                  div(id = "p_heat_pval_box", style = "height:600px;",
                      plotOutput("p_heat_pval", height = "100%")),
                  options = list(handles = "se")
                )
              ),
              column(
                12, h4("Direction"),
                jqui_resizable(
                  div(id = "p_heat_dir_box", style = "height:600px;",
                      plotOutput("p_heat_dir", height = "100%")),
                  options = list(handles = "se")
                )
              )
            ),
            br(),
            h5("Notes:"),
            tags$ul(
              tags$li("Signed score = -log10(p) × direction, where direction = +1 if Responder median > Non-responder median, -1 if lower, 0 if tie."),
              tags$li("P-value heatmap shows raw p-values comparing Responder vs Non-responder."),
              tags$li("Direction heatmap shows +1 (Responder higher), 0 (tie), -1 (Non-responder higher).")
            )
          ),
          tabPanel(
            "Per-cell-type plot",
            br(),
            jqui_resizable(
              div(id = "p_facet_plot_box", style = "height:700px;",
                  plotOutput("p_facet_plot", height = "100%")),
              options = list(handles = "se")
            ),
            br(),
            h4("Per-cell-type statistics"),
            DTOutput("p_stats_table")
          )
        )
      )
    )
  ),
  
  ## ---- Cell-level tab ----
  tabPanel(
    "Cell-level",
    sidebarLayout(
      sidebarPanel(
        width = 3,
        h4("Cross-study per-gene viewer (Pre / Post folders)"),
        selectizeInput(
          "c_gene_sel",
          "Gene (HGNC symbol or Ensembl ID):",
          choices = gene_index, multiple = FALSE,
          options = list(placeholder = "Type to search...", create = FALSE)
        ),
        checkboxInput(
          "c_use_symbol",
          "Search/display by HGNC symbol (fallback to Ensembl if missing)",
          TRUE
        ),
        radioButtons(
          "c_prepost", "Timepoint:",
          choices = c("Pre","Post"), selected = "Post", inline = TRUE
        ),
        tags$hr(),
        h4("Filters"),
        sliderInput(
          "c_fdr_thr", "FDR (p_val_adj) threshold:",
          min = 1e-10, max = 0.25, value = 0.05, step = 0.005, ticks = FALSE
        ),
        sliderInput(
          "c_pct_thr", "Min pct expressed in either group (pct.1 or pct.2):",
          min = 0, max = 1, value = 0, step = 0.05
        ),
        checkboxInput("c_show_only_sig", "Show only significant rows (p_val_adj < FDR)", FALSE),
        tags$hr(),
        downloadButton("c_dl_table", "Download table (.tsv)")
      ),
      mainPanel(
        tabsetPanel(
          id = "c_subtabs",
          tabPanel("Table", br(), DTOutput("c_res_table")),
          tabPanel(
            "Dot heatmap", br(),
            jqui_resizable(
              div(id = "c_dot_heatmap_box", style = "height:520px;",
                  plotOutput("c_dot_heatmap", height = "100%")),
              options = list(handles = "se")
            )
          ),
          tabPanel(
            "Paired CPM bars", br(),
            jqui_resizable(
              div(id = "c_cpm_bars_box", style = "height:620px;",
                  plotOutput("c_cpm_bars", height = "100%")),
              options = list(handles = "se")
            )
          ),
          tabPanel(
            "Lollipop (log2FC)", br(),
            jqui_resizable(
              div(id = "c_lollipop_box", style = "height:620px;",
                  plotOutput("c_lollipop", height = "100%")),
              options = list(handles = "se")
            )
          )
        )
      )
    )
  )
)

## ======================
## =       SERVER       =
## ======================

server <- function(input, output, session) {
  h <- function(x, default) if (is.null(x)) default else x
  
  ## -----------------------------
  ## Patient-level state & logic
  ## -----------------------------
  p_params <- reactiveVal(list(gene = "A1BG", prepost = "Post", celltypes = NULL))
  
  observeEvent(input$p_run, {
    p_params(list(
      gene      = input$p_gene,
      prepost   = if (input$p_prepost == "Both") NULL else input$p_prepost,
      celltypes = input$p_celltypes
    ))
  }, ignoreInit = FALSE)
  
  p_meta_filt <- reactive({
    par <- p_params(); req(par)
    filter_meta(meta_all, par$prepost %||% "Both", par$celltypes)
  })
  
  p_studies_available <- reactive({
    mf <- p_meta_filt(); req(nrow(mf) > 0)
    sort(unique(mf$study))
  })
  
  output$p_plot_study_picker <- renderUI({
    ch <- p_studies_available()
    selectInput("p_plot_study", "Dataset for Per-cell-type plot:", choices = ch, selected = ch[1])
  })
  
  p_combined_stats <- reactive({
    par <- p_params(); req(par)
    mf <- p_meta_filt(); req(nrow(mf) > 0)
    purrr::map_dfr(p_studies_available(), function(st) {
      res <- compute_one_study(st, par$gene, par$prepost, par$celltypes, mf, expr_mat_all)
      if (is.null(res)) return(NULL)
      res$stats
    })
  })
  
  p_chosen_patient_df <- reactive({
    par <- p_params(); req(par, input$p_plot_study)
    mf <- p_meta_filt()
    res <- compute_one_study(input$p_plot_study, par$gene, par$prepost, par$celltypes, mf, expr_mat_all)
    validate(need(!is.null(res), "No data for this study with current filters."))
    res$patient
  })
  
  p_chosen_stats_tbl <- reactive({
    par <- p_params(); req(par, input$p_plot_study)
    mf <- p_meta_filt()
    res <- compute_one_study(input$p_plot_study, par$gene, par$prepost, par$celltypes, mf, expr_mat_all)
    validate(need(!is.null(res), "No data for this study with current filters."))
    res$stats
  })
  
  ## Patient plots
  output$p_heat_signed <- renderPlot({
    df <- p_combined_stats() %>% dplyr::select(CellType = cell_type, Study, ThisStudy = signed_score)
    req(nrow(df) > 0)
    ggplot(df, aes(x = Study, y = forcats::fct_rev(factor(CellType)), fill = ThisStudy)) +
      geom_tile(color = "white") +
      geom_text(aes(label = ifelse(is.na(ThisStudy), "NA", sprintf("%.2f", ThisStudy))), size = 3) +
      scale_fill_gradient2(low = scales::muted("blue"), mid = "white", high = scales::muted("red"),
                           midpoint = 0, na.value = "grey90") +
      labs(x = NULL, y = NULL, fill = "Signed\nscore") +
      theme_minimal(base_size = 13) +
      theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  }, height = function() h(input$p_heat_signed_box_height, 600))
  
  output$p_heat_pval <- renderPlot({
    df <- p_combined_stats() %>%
      dplyr::select(CellType = cell_type, Study, p_raw = p_value)
    req(nrow(df) > 0)
    
    df <- df %>%
      dplyr::mutate(
        p_raw    = suppressWarnings(as.numeric(p_raw)),
        sig_star = ifelse(!is.na(p_raw) & p_raw < 0.05, "*", ""),
        lab      = ifelse(
          is.na(p_raw), "NA",
          paste0(format(p_raw, digits = 2, scientific = TRUE), sig_star)
        )
      )
    
    ggplot(df, aes(x = Study, y = forcats::fct_rev(factor(CellType)), fill = p_raw)) +
      geom_tile(color = "white") +
      # numeric labels
      geom_text(aes(label = ifelse(is.na(p_raw),
                                   "NA",
                                   format(p_raw, digits = 2, scientific = TRUE))),
                size = 3, color = "black") +
      # bold white stars, slightly shifted
      geom_text(aes(label = ifelse(!is.na(p_raw) & p_raw < 0.05, "*", "")),
                color = "white", fontface = "bold", size = 6,
                nudge_x = 0.15,  # move right
                nudge_y = -0.15) +  # move down
      scale_fill_gradient(
        low = scales::muted("blue"),
        high = scales::muted("red"),
        limits = c(0, 1),
        oob = scales::squish,
        na.value = "grey90",
        guide = guide_colorbar(reverse = TRUE)
      ) +
      labs(x = NULL, y = NULL, fill = "p-value") +
      theme_minimal(base_size = 13) +
      theme(panel.grid = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1))
  }, height = function() h(input$p_heat_pval_box_height, 600))
  
  output$p_heat_dir <- renderPlot({
    df <- p_combined_stats() %>%
      dplyr::select(CellType = cell_type, Study, ThisStudy = direction) %>%
      dplyr::mutate(Direction = factor(ThisStudy, levels = c(-1,0,1), labels = c("-1","0","+1")))
    req(nrow(df) > 0)
    ggplot(df, aes(x = Study, y = forcats::fct_rev(factor(CellType)), fill = Direction)) +
      geom_tile(color = "white") +
      geom_text(aes(label = as.character(Direction)), size = 3) +
      scale_fill_manual(values = c("-1" = scales::muted("blue"), "0" = "grey85", "+1" = scales::muted("red")),
                        na.value = "grey90") +
      labs(x = NULL, y = NULL, fill = "Direction") +
      theme_minimal(base_size = 13) +
      theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  }, height = function() h(input$p_heat_dir_box_height, 600))
  
  output$p_facet_plot <- renderPlot({
    pec <- p_chosen_patient_df()
    par <- p_params(); req(nrow(pec) > 0, par$gene)
    ggplot(pec, aes(x = response_group, y = expr, fill = response_group)) +
      geom_boxplot(width = 0.35, outlier.shape = NA, alpha = 0.6) +
      geom_jitter(width = 0.08, size = 1.6, alpha = 0.7) +
      stat_summary(fun = median, geom = "point", size = 5, colour = "black", shape = 95) +
      facet_wrap(~ cell_type, scales = "free_y") +
      labs(
        title = paste(par$gene, "Expression (patient-level)",
                      ifelse(is.null(par$prepost), "Pre+Post", par$prepost),
                      input$p_plot_study),
        x = NULL, y = "Expression level"
      ) +
      theme_minimal(base_size = 14) +
      theme(legend.position = "none") +
      ggpubr::stat_compare_means(comparisons = list(c("Responder","Non-responder")),
                                 method = "wilcox.test", label = "p.signif")
  }, height = function() h(input$p_facet_plot_box_height, 700))
  
  output$p_stats_table <- DT::renderDT({
    st <- p_chosen_stats_tbl(); req(nrow(st) > 0)
    st %>%
      dplyr::mutate(study = Study, p_value = signif(p_value, 3), signed_score = round(signed_score, 3)) %>%
      dplyr::select(study, cell_type, n_R, n_NR, median_R, median_NR, direction, p_value, signed_score) %>%
      DT::datatable(options = list(pageLength = 10), rownames = FALSE)
  })
  
  ## ---------------------------
  ## Cell-level (on-demand)    -
  ## ---------------------------
  
  # Choose file set based on timepoint
  current_files <- reactive({
    if (req(input$c_prepost) == "Pre") files_pre else files_post
  })
  
  # Read just the selected gene from all files for the current timepoint
  c_filtered <- reactive({
    # Explicit dependencies so any change re-runs the reactive:
    tp        <- req(input$c_prepost)
    gene_sel  <- req(input$c_gene_sel)
    fdr_thr   <- req(input$c_fdr_thr)
    only_sig  <- isTRUE(input$c_show_only_sig)
    use_sym   <- isTRUE(input$c_use_symbol)
    
    files <- if (tp == "Pre") current_files() else current_files()  # current_files() already reacts to c_prepost
    validate(need(length(files) > 0, "No cell-level files found for this timepoint."))
    
    # Stream per file and keep only matching rows
    pieces <- lapply(files, function(fp) {
      read_gene_from_file(fp, gene_sel, use_sym, KEEP_COLS)
    })
    dat <- dplyr::bind_rows(pieces)
    validate(need(nrow(dat) > 0, "Selected gene not found in these files."))
    
    # Safety: coerce numerics after bind (some files might slip odd types)
    num_cols <- intersect(
      c("p_val","avg_log2FC","p_val_adj","Responder_CPM","Non_responder_CPM","pct.1","pct.2"),
      names(dat)
    )
    for (cc in num_cols) dat[[cc]] <- suppressWarnings(as.numeric(dat[[cc]]))
    
    # Derive helper columns
    dat <- dat %>%
      mutate(
        pre_post = tp,
        `pct.1`  = `pct.1` %||% NA_real_,
        `pct.2`  = `pct.2` %||% NA_real_,
        pct_min  = pmin(`pct.1`, `pct.2`),
        pct_diff = `pct.1` - `pct.2`,
        cpm_diff = (Responder_CPM %||% NA_real_) - (Non_responder_CPM %||% NA_real_),
        direction = case_when(
          avg_log2FC >  0 ~ "Responder up",
          avg_log2FC <  0 ~ "Non-responder up",
          TRUE            ~ "No change"
        ),
        signif = !is.na(p_val_adj) & p_val_adj < fdr_thr
      ) %>%
      filter(is.na(pct_min) | pct_min >= input$c_pct_thr)
    
    # >>> Apply the checkbox filter here <<<
    if (only_sig) {
      dat <- dat %>% filter(signif)
    }
    
    dat %>% arrange(Study, cluster)
  })
  
  output$c_res_table <- DT::renderDT({
    dat <- c_filtered(); req(nrow(dat) > 0)
    show_cols <- c("Study","pre_post","cluster","gene","hgnc_symbol","gene_biotype",
                   "p_val","avg_log2FC","pct.1","pct.2","p_val_adj",
                   "Responder_CPM","Non_responder_CPM","pct_diff","cpm_diff","direction")
    show_cols <- intersect(show_cols, names(dat))
    DT::datatable(
      dat[, show_cols, drop = FALSE],
      rownames = FALSE, filter = "top",
      extensions = "Buttons",
      options = list(pageLength = 25, dom = "Bfrtip",
                     buttons = c("copy", "csv", "excel"),
                     scrollX = TRUE)
    ) %>%
      DT::formatRound(c("p_val","p_val_adj"), 3) %>%
      DT::formatRound(c("avg_log2FC","pct.1","pct.2","pct_diff"), 3) %>%
      DT::formatRound(c("Responder_CPM","Non_responder_CPM","cpm_diff"), 3)
  })
  
  output$c_dl_table <- downloadHandler(
    filename = function() paste0(
      "per_gene_table_",
      gsub("[^A-Za-z0-9]+","_", input$c_gene_sel),
      "_", tolower(input$c_prepost), ".tsv"
    ),
    content = function(file) data.table::fwrite(c_filtered(), file, sep = "\t", quote = FALSE, na = "NA")
  )
  
  output$c_dot_heatmap <- renderPlot({
    dat <- c_filtered(); req(nrow(dat) > 0)
    ttl_suffix <- paste0(" — ", input$c_prepost)
    ggplot(dat, aes(x = cluster, y = Study)) +
      geom_point(aes(size = -log10(p_val_adj + 1e-300), fill = avg_log2FC),
                 shape = 21,
                 color = ifelse(dat$signif, "black", "grey70"),
                 alpha = 0.9) +
      scale_fill_gradient2(name = "log2FC\n(Resp/Nonresp)",
                           low = scales::muted("blue"), mid = "white", high = scales::muted("red"),
                           midpoint = 0) +
      scale_size_continuous(name = expression(-log[10]("FDR")), range = c(2, 10)) +
      labs(x = "Cell type (cluster)", y = "Study",
           title = paste0("Dot heatmap — ", input$c_gene_sel, ttl_suffix)) +
      theme_minimal(base_size = 12) +
      theme(axis.text.x = element_text(angle = 40, hjust = 1),
            panel.grid.minor = element_blank())
  }, height = function() h(input$c_dot_heatmap_box_height, 520))
  
  output$c_cpm_bars <- renderPlot({
    dat <- c_filtered(); req(nrow(dat) > 0)
    validate(need(all(c("Responder_CPM","Non_responder_CPM") %in% names(dat)),
                  "CPM columns not available in the loaded files."))
    dat_long <- dat %>%
      select(Study, pre_post, cluster, Responder_CPM, Non_responder_CPM) %>%
      tidyr::pivot_longer(c(Responder_CPM, Non_responder_CPM), names_to = "Group", values_to = "CPM")
    ttl_suffix <- paste0(" — ", input$c_prepost)
    ggplot(dat_long, aes(x = cluster, y = CPM, fill = Group)) +
      geom_col(position = position_dodge(width = 0.75), width = 0.7) +
      facet_wrap(~ Study, scales = "free_y") +
      labs(x = "Cell type (cluster)", y = "CPM",
           title = paste0("Responder vs Non-responder CPM — ", input$c_gene_sel, ttl_suffix)) +
      theme_minimal(base_size = 12) +
      theme(axis.text.x = element_text(angle = 40, hjust = 1))
  }, height = function() h(input$c_cpm_bars_box_height, 620))
  
  output$c_lollipop <- renderPlot({
    dat <- c_filtered(); req(nrow(dat) > 0)
    dat <- dat %>% dplyr::mutate(cluster_f = forcats::fct_reorder(cluster, avg_log2FC, .fun = median))
    ttl_suffix <- paste0(" — ", input$c_prepost)
    ggplot(dat, aes(y = cluster_f, x = avg_log2FC, color = Study)) +
      geom_segment(aes(yend = cluster_f, x = 0, xend = avg_log2FC),
                   linewidth = 0.7, alpha = 0.6) +
      geom_point(size = 3) +
      geom_vline(xintercept = 0, linetype = "dashed") +
      labs(x = "log2FC (Responder / Non-responder)", y = "Cell type (cluster)",
           title = paste0("Effect sizes by cell type — ", input$c_gene_sel, ttl_suffix)) +
      theme_minimal(base_size = 12)
  }, height = function() h(input$c_lollipop_box_height, 620))
}

shinyApp(ui, server)
