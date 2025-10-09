# app.R  — Multi-study heatmaps from combined_meta / combined_expr
## ----- Packages -----
library(shiny)
library(tidyverse)
library(data.table)
library(ggpubr)
library(stringr)
library(purrr)
library(scales)
library(DT)

message("== start ==")
message("reading data ...");  df <- readRDS("data/combined_expr.rds")
message("data read OK")
message("reading META: data/combined_meta.tsv")
meta_all <- readr::read_tsv("data/combined_meta.tsv", show_col_types = FALSE)
# Show columns the server actually sees
message("META columns: ", paste(names(meta_all), collapse = ", "))
message("building UI ...")

## ----- Config: file paths -----
meta_path  <- "data/combined_meta.tsv"
#expr_path  <- "data/combined_expr.tsv"
expr_path <- "data/combined_expr.rds"
app_title  <- "Responder vs Non-responder"

## ----- Helpers -----
# Expression loader that:
#  - normalizes sample IDs (strip .x/.y/.1 suffixes, quotes/backticks, trim)
#  - collapses duplicate sample columns by row-mean
#  - returns genes x samples matrix (rownames = gene, colnames = sample IDs)
load_expr <- function(path) {
  #expr_dt <- data.table::fread(path, sep = "\t", header = TRUE)
  #data.table::setnames(expr_dt, old = names(expr_dt)[1], new = "gene")
  # detect file type
  if (grepl("\\.rds$", path, ignore.case = TRUE)) {
    message("Reading RDS file: ", path)
    return(readRDS(path))
  }
  
  message("Reading TSV file: ", path)
  expr_dt <- data.table::fread(path, sep = "\t", header = TRUE)
  data.table::setnames(expr_dt, old = names(expr_dt)[1], new = "gene")
  
  
  fix_ids <- function(x) {
    x %>%
      stringr::str_trim() %>%
      stringr::str_replace_all('^[\'"]|[\'"]$', '') %>%
      stringr::str_replace_all("`", "") %>%
      stringr::str_replace("(\\.x|\\.y)$", "") %>%
      stringr::str_replace("\\.\\d+$", "")
  }
  colnames(expr_dt) <- c("gene", fix_ids(colnames(expr_dt)[-1]))
  
  
  # melt -> collapse duplicates -> cast back wide
  long_dt <- data.table::melt(
    data.table::as.data.table(expr_dt),
    id.vars = "gene",
    variable.name = "ID",
    value.name = "value"
  )
  long_dt <- long_dt[, .(value = mean(value, na.rm = TRUE)), by = .(gene, ID)]
  wide_dt <- data.table::dcast(long_dt, gene ~ ID, value.var = "value")
  
  wide_df <- as.data.frame(wide_dt) %>%
    mutate(
      gene = as.character(gene),
      gene = stringr::str_replace_all(gene, '^[\'"]|[\'"]$', ''),
      gene = stringr::str_trim(gene),
      gene = dplyr::na_if(gene, "")
    ) %>%
    filter(!is.na(gene)) %>%
    group_by(gene) %>%
    summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
  
  wide_df %>% tibble::column_to_rownames("gene")
}

normalize_id <- function(x) {
  x %>%
    stringr::str_trim() %>%
    stringr::str_replace_all('^[\'"]|[\'"]$', '') %>%  # strip outer quotes
    stringr::str_replace_all("`", "") %>%              # drop backticks
    stringr::str_replace("(\\.x|\\.y)$", "") %>%       # <<< important: drop .x/.y
    stringr::str_replace("\\.\\d+$", "")               # drop trailing .1/.2/...
}

# Filter/meta normalizer; preserves `study`
filter_meta <- function(meta_df, prepost, celltypes) {
  meta_f <- meta_df
  # response_group is created in meta_base; keep rows with known response if you intended that
  meta_f <- meta_f %>% dplyr::filter(is.na(response_group) | response_group %in% c("Responder","Non-responder"))
  if (!is.null(prepost) && prepost %in% c("Pre","Post") && "pre_post" %in% names(meta_f)) {
    meta_f <- meta_f %>% dplyr::filter(pre_post %in% prepost)
  }
  if (!is.null(celltypes) && length(celltypes) > 0 && "cell_type" %in% names(meta_f)) {
    meta_f <- meta_f %>% dplyr::filter(cell_type %in% celltypes)
  }
  meta_f %>% dplyr::select(dplyr::any_of(c("ID","donor_id","cell_type","pre_post","response_group","study")))
}

# Compute patient-level df and stats for ONE study label
safe_wilcox <- function(df) {
  # df has columns: expr, response_group
  # keep only R/NR rows
  df <- df %>% dplyr::filter(response_group %in% c("Responder","Non-responder"))
  n_R  <- sum(df$response_group == "Responder",    na.rm = TRUE)
  n_NR <- sum(df$response_group == "Non-responder", na.rm = TRUE)
  if (isTRUE(n_R > 0 & n_NR > 0)) {
    # guard against NA expr and zero-length numeric vectors
    tryCatch(
      {
        wilcox.test(expr ~ response_group, data = df, exact = FALSE)$p.value
      },
      error = function(e) NA_real_
    )
  } else {
    NA_real_
  }
}

compute_one_study <- function(study_label, gene, prepost, celltypes, meta_all_df, expr_all_mat) {
  mf_all <- meta_all_df
  # keep only this study
  mf <- mf_all %>% dplyr::filter(study == study_label)
  # ensure canonical columns exist
  # (ID, donor_id, cell_type, pre_post, response_group, study)
  if (nrow(mf) == 0) return(NULL)
  # drop unlabeled outcomes right here (matches your local behavior)
  mf <- mf %>% dplyr::filter(response_group %in% c("Responder","Non-responder"))
  # optional filters
  if (!is.null(prepost) && prepost %in% c("Pre","Post") && "pre_post" %in% names(mf)) {
    mf <- mf %>% dplyr::filter(pre_post %in% prepost)
  }
  if (!is.null(celltypes) && length(celltypes) > 0 && "cell_type" %in% names(mf)) {
    mf <- mf %>% dplyr::filter(cell_type %in% celltypes)
  }
  # intersect with expression columns
  common_ids <- intersect(colnames(expr_all_mat), mf$ID)
  if (length(common_ids) == 0) return(NULL)
  if (!(gene %in% rownames(expr_all_mat))) return(NULL)
  # build sample-level df (drop NA expr)
  sample_df <- tibble::tibble(
    ID   = common_ids,
    expr = as.numeric(expr_all_mat[gene, common_ids, drop = TRUE])
  ) %>%
    dplyr::left_join(mf, by = "ID") %>%
    dplyr::filter(!is.na(expr)) %>%
    dplyr::mutate(
      response_group = factor(response_group, levels = c("Non-responder","Responder"))
    )
  if (nrow(sample_df) == 0) return(NULL)
  # patient-level (mean per donor × cell_type × response)
  patient_df <- sample_df %>%
    dplyr::group_by(cell_type, donor_id, response_group) %>%
    dplyr::summarise(expr = mean(expr, na.rm = TRUE), .groups = "drop")
  if (nrow(patient_df) == 0) {
    return(list(
      patient = tibble::tibble(cell_type=character(), donor_id=character(),
                               response_group=factor(levels=c("Non-responder","Responder")),
                               expr=numeric(), Study=character()),
      stats   = tibble::tibble(Study=character(), cell_type=character(),
                               n_R=integer(), n_NR=integer(),
                               median_R=numeric(), median_NR=numeric(),
                               direction=integer(), p_value=numeric(), signed_score=numeric())
    ))
  }
  # stats per cell_type with all guards
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
        Study = study_label,
        cell_type = df$cell_type[1],
        n_R = n_R, n_NR = n_NR,
        median_R = medR, median_NR = medNR,
        direction = dir,
        p_value = p,
        signed_score = dplyr::if_else(is.na(p), NA_real_, -log10(p) * dir)
      )
    })
  list(
    patient = patient_df %>% dplyr::mutate(Study = study_label),
    stats   = stats_tbl
  )
}


## ----- Load metadata & expression -----
meta_all <- readr::read_tsv(meta_path, show_col_types = FALSE) %>%
  mutate(
    across(c(ID, donor_id, pre_post, cell_type, study), ~ as.character(.x)),
    ID = normalize_id(ID) %>%
      str_trim() %>%
      str_replace_all('^[\'"]|[\'"]$', '') %>%
      str_replace_all("`", "") %>%
      str_replace("\\.\\d+$", ""),
    donor_id = ifelse(is.na(donor_id) | donor_id == "", ID, donor_id),
    outcome  = trimws(outcome),
    response_group = case_when(
      outcome %in% c("R","PR","CR","OR","ICB-PR","YES","Yes") ~ "Responder",
      outcome %in% c("NR","SD","ICB-SD","POST-ICI (RESISTANT)","Post-ICI (resistant)","NO","No") ~ "Non-responder",
      outcome %in% c("UT","UNTREATED","NE","N/A","NA","BASELINE","UNKNOWN") ~ NA_character_,
      TRUE ~ NA_character_
    )
  )

# --- DROP studies where response_group is entirely NA (your local behavior)
study_keep <- meta_all %>%
  dplyr::group_by(study) %>%
  dplyr::summarise(any_labeled = any(!is.na(response_group)), .groups = "drop") %>%
  dplyr::filter(any_labeled) %>%
  dplyr::pull(study)
meta_all <- meta_all %>% dplyr::filter(study %in% study_keep)

#write.csv(meta_all, "G:/Ranran/ICB_Response/meta_all.csv", row.names = FALSE)
expr_mat_all <- load_expr(expr_path)

all_genes <- rownames(expr_mat_all) %>% sort()
all_cell_types <- meta_all$cell_type %>% unique() %>% sort()

## ----- UI -----
ui <- fluidPage(
  titlePanel(app_title),
  sidebarLayout(
    sidebarPanel(
      selectizeInput(
        "gene", "Gene:", choices = all_genes, selected = "CD8A",
        options = list(placeholder = "Type a gene symbol…", create = FALSE)
      ),
      radioButtons(
        "pre_post_choice", "Timepoint:",
        choices = c("Pre", "Post", "Both"), selected = "Post", inline = TRUE
      ),
      selectizeInput(
        "cell_type_filter", "Cell types (optional):",
        choices = all_cell_types, multiple = TRUE
      ),
      uiOutput("plot_study_picker"),
      actionButton("run", "Update", class = "btn-primary")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel(
          "Heatmaps",
          br(),
          fluidRow(
            column(12, h4("Score"),     plotOutput("heat_signed", height = "600px")),
            column(12, h4("P-value"),   plotOutput("heat_pval",   height = "600px")),
            column(12, h4("Direction"), plotOutput("heat_dir",    height = "600px"))
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
          plotOutput("facet_plot", height = "700px"),
          br(),
          h4("Per-cell-type statistics"),
          DTOutput("stats_table")
        )
      )
    )
  )
)

## ----- Server -----
server <- function(input, output, session) {
  
  # Parameter state
  params <- reactiveVal(list(
    gene = "CD8A",
    prepost = "Post",
    celltypes = NULL
  ))
  
  observeEvent(input$run, {
    params(list(
      gene = input$gene,
      prepost = if (input$pre_post_choice == "Both") NULL else input$pre_post_choice,
      celltypes = input$cell_type_filter
    ))
  }, ignoreInit = FALSE)
  
  # Filtered meta (keeps 'study')
  meta_filt <- reactive({
    par <- params(); req(par)
    filter_meta(meta_all, par$prepost %||% "Both", par$celltypes)
  })
  
  # Which studies are available after filters?
  studies_available <- reactive({
    mf <- meta_filt(); req(nrow(mf) > 0)
    sort(unique(mf$study))
  })
  
  output$plot_study_picker <- renderUI({
    ch <- studies_available()
    selectInput("plot_study", "Dataset for Per-cell-type plot:", choices = ch, selected = ch[1])
  })
  
  # Combined stats across all studies (for heatmaps)
  combined_stats <- reactive({
    par <- params(); req(par)
    mf <- meta_filt(); req(nrow(mf) > 0)
    map_dfr(studies_available(), function(st) {
      res <- compute_one_study(st, par$gene, par$prepost, par$celltypes, mf, expr_mat_all)
      if (is.null(res)) return(NULL)
      res$stats
    })
  })
  
  # Patient df + stats for the selected study (per-cell-type tab)
  chosen_patient_df <- reactive({
    par <- params(); req(par, input$plot_study)
    mf <- meta_filt()
    res <- compute_one_study(input$plot_study, par$gene, par$prepost, par$celltypes, mf, expr_mat_all)
    validate(need(!is.null(res), "No data for this study with current filters."))
    res$patient
  })
  
  chosen_stats_tbl <- reactive({
    par <- params(); req(par, input$plot_study)
    mf <- meta_filt()
    res <- compute_one_study(input$plot_study, par$gene, par$prepost, par$celltypes, mf, expr_mat_all)
    validate(need(!is.null(res), "No data for this study with current filters."))
    res$stats
  })
  
  ## ---- Heatmaps ----
  output$heat_signed <- renderPlot({
    df <- combined_stats() %>% select(CellType = cell_type, Study, ThisStudy = signed_score)
    req(nrow(df) > 0)
    ggplot(df, aes(x = Study, y = fct_rev(factor(CellType)), fill = ThisStudy)) +
      geom_tile(color = "white") +
      geom_text(aes(label = ifelse(is.na(ThisStudy), "NA", sprintf("%.2f", ThisStudy))), size = 3) +
      scale_fill_gradient2(low = muted("blue"), mid = "white", high = muted("red"),
                           midpoint = 0, na.value = "grey90") +
      labs(x = NULL, y = NULL, fill = "Signed\nscore") +
      theme_minimal(base_size = 13) +
      theme(panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  })
  
  output$heat_pval <- renderPlot({
    df <- combined_stats() %>%
      dplyr::select(CellType = cell_type, Study, p_raw = p_value)
    
    req(nrow(df) > 0)
    
    df <- df %>%
      dplyr::mutate(
        # force numeric; anything non-numeric becomes NA
        p_raw = suppressWarnings(as.numeric(p_raw)),
        # label text (keep as character for printing on tiles)
        lab = ifelse(is.na(p_raw), "NA", format(p_raw, digits = 2, scientific = TRUE))
      )
    
    ggplot(df, aes(x = Study, y = forcats::fct_rev(factor(CellType)), fill = p_raw)) +
      geom_tile(color = "white") +
      geom_text(aes(label = lab), size = 3) +
      scale_fill_gradient(
        low = scales::muted("blue"),   # small p -> light
        high = scales::muted("red"), # large p -> dark
        limits = c(0, 1),
        oob = scales::squish,
        na.value = "grey90",
        guide = guide_colorbar(reverse = TRUE)  # reverse legend so 0 is at the top
      ) +
      labs(x = NULL, y = NULL, fill = "p-value") +
      theme_minimal(base_size = 13) +
      theme(
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)
      )
  })
  
  output$heat_dir <- renderPlot({
    df <- combined_stats() %>%
      select(CellType = cell_type, Study, ThisStudy = direction) %>%
      mutate(Direction = factor(ThisStudy, levels = c(-1,0,1), labels = c("-1","0","+1")))
    req(nrow(df) > 0)
    ggplot(df, aes(x = Study, y = fct_rev(factor(CellType)), fill = Direction)) +
      geom_tile(color = "white") +
      geom_text(aes(label = as.character(Direction)), size = 3) +
      scale_fill_manual(values = c("-1" = muted("blue"), "0" = "grey85", "+1" = muted("red")),
                        na.value = "grey90") +
      labs(x = NULL, y = NULL, fill = "Direction") +
      theme_minimal(base_size = 13) +
      theme(panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  })
  
  ## ---- Per-cell-type plot + table (selected study) ----
  output$facet_plot <- renderPlot({
    pec <- chosen_patient_df()
    par <- params(); req(nrow(pec) > 0, par$gene)
    
    ggplot(pec, aes(x = response_group, y = expr, fill = response_group)) +
      geom_boxplot(width = 0.35, outlier.shape = NA, alpha = 0.6) +
      geom_jitter(width = 0.08, size = 1.6, alpha = 0.7) +
      stat_summary(fun = median, geom = "point", size = 5, colour = "black", shape = 95) +
      facet_wrap(~ cell_type, scales = "free_y") +
      labs(
        title = paste(par$gene, "Expression (patient-level)",
                      ifelse(is.null(par$prepost), "Pre+Post", par$prepost),
                      input$plot_study),
        x = NULL, y = "Expression level"
      ) +
      theme_minimal(base_size = 14) +
      theme(legend.position = "none") +
      stat_compare_means(
        comparisons = list(c("Responder","Non-responder")),
        method = "wilcox.test", label = "p.signif"
      )
  }, res = 110)
  
  output$stats_table <- renderDT({
    st <- chosen_stats_tbl()
    req(nrow(st) > 0)
    st %>%
      mutate(
        study = Study,
        p_value = signif(p_value, 3),
        signed_score = round(signed_score, 3)
      ) %>%
      select(study, cell_type, n_R, n_NR, median_R, median_NR, direction, p_value, signed_score) %>%
      datatable(options = list(pageLength = 10), rownames = FALSE)
  })
}

shinyApp(ui, server)
