/*
 * popGenExp_doPairwiseComparisonOnHapsSharing.cpp
 *
 *  Created on: Jun 11, 2021
 *      Author: nick
 */




#include "popGenExp.hpp"
#include "elucidator/objects/BioDataObject.h"
#include <njhseq/objects/counters/DNABaseCounter.hpp>
#include <njhseq/PopulationGenetics.h>
#include <njhseq/objects/dataContainers/BasicPointMatrix.hpp>



namespace njhseq {




int popGenExpRunner::calc_pairwise_ccc_on_haps_sharing(const njh::progutils::CmdArgs & inputCommands){
	double minimumLociCoverageToKeepSamples = 0.90;
	HapsEncodedMatrix::SetWithExternalPars pars;
  uint32_t pairwise_factor_bin_size = 1000;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
  setUp.setOption(pairwise_factor_bin_size, "--pairwise_factor_bin_size", "pairwise_factor_bin_size");
	setUp.setOption(minimumLociCoverageToKeepSamples, "--minimumLociCoverageToKeepSamples", "minimum Loci Coverage To Keep Samples in post analysis steps, must have reads for at least this frction of the total loci");
  pars.setDefaults(setUp);

  setUp.processDirectoryOutputName(bfs::path(bfs::basename(pars.tableFnp)).string() + "_ccc_rmse_TODAY", true);
	setUp.finishSetUp(std::cout);

	setUp.startARunLog(setUp.pars_.directoryName_);


	setUp.timer_.setLapName("initial");
	setUp.timer_.startNewLap("encode haplotypes");
  HapsEncodedMatrix haps(pars);
	setUp.timer_.startNewLap("get hap probabilities");
	haps.calcHapProbs();
	setUp.timer_.startNewLap("add relative abundances");
  haps.add_relative_abundance();
  setUp.timer_.startNewLap("calc rmse and ccc");
	auto measures = haps.calc_ccc_rmse_measures(pairwise_factor_bin_size, setUp.pars_.verbose_);
  setUp.timer_.startNewLap("writing output matrices");
	OutputStream outSampNamesOut(njh::files::make_path(setUp.pars_.directoryName_, "sampleNames.tab.txt"));
	outSampNamesOut << njh::conToStr(haps.sampNamesVec_, "\n") << std::endl;
  {
	  auto ccc_out_fnp = njh::files::make_path(setUp.pars_.directoryName_, "ccc_on_targets_shared.tab.txt.gz");
	  OutputStream ccc_out(ccc_out_fnp);
	  for(const auto & ccc_row : measures.ccc){
	    ccc_out << njh::conToStr(ccc_row, "\t") << std::endl;
	  }
  }
  {
	  auto rmse_out_fnp = njh::files::make_path(setUp.pars_.directoryName_, "rmse_on_targets_shared.tab.txt.gz");
	  OutputStream rmse_out(rmse_out_fnp);
	  for(const auto & rmse_row : measures.rmse){
	    rmse_out << njh::conToStr(rmse_row, "\t") << std::endl;
	  }
  }

  {
	  auto targets_shared_out_fnp = njh::files::make_path(setUp.pars_.directoryName_, "targets_shared.tab.txt.gz");
	  OutputStream targets_shared_out(targets_shared_out_fnp);
	  for(const auto & targets_shared_row : measures.targets_shared){
	    targets_shared_out << njh::conToStr(targets_shared_row, "\t") << std::endl;
	  }
  }
  setUp.timer_.startNewLap("getting loci coverage info");

	std::unordered_map<std::string, double> lociCoveragePerSample = haps.getTargetCoveragePerSample();
	{
		table numTargetsPerSample = haps.getTableNumberTargetsPerSample(minimumLociCoverageToKeepSamples);
		OutputStream lociCoverageOut(njh::files::make_path(setUp.pars_.directoryName_, "loci_coverage_per_sample_info.tsv"));
		numTargetsPerSample.outPutContents(lociCoverageOut, "\t");
	}
	setUp.timer_.logLapTimes(setUp.rLog_.runLogFile_, true, 6, true);
	return 0;
}


static const std::string static_qmd_on_clustering = R"QMD(---
title: Processing groups
---


```{r setup}
library(HaplotypeRainbows)
library(tidyverse)
library(DT)

create_dt <- function(x) {
  DT::datatable(
    x,
    extensions = 'Buttons',
    options = list(
      dom = 'Blfrtip',
      buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
      lengthMenu = list(c(10, 25, 50,-1),
                        c(10, 25, 50, "All"))
    ),
    filter = "top"
  )
}
```

## Creating a haplotype rainbow

Create a haplotype rainbow and sort and cluster by the groups determined

```{r}


clusters = readr::read_tsv("clusters.tsv")
haps = readr::read_tsv("haps.tsv.gz")

haps_hr = haplotype_rainbow(haps,
                                               sample_col = "library_sample_name",
                                               target_col = "target_name",
                                               popuid_col = "seq",
                                               rel_abund_col = "within_sample_freq")$prep()

haps_hr$set_sample_meta(clusters, "sample")
haps_hr$sort_samples_by_clustering(abundance_weighted = T)$sort_samples_by_meta("group")$add_sample_cluster_gaps()

haps_hr$save_pdf(
  haps_hr$add_sample_annotation_to_plot(haps_hr$plot()), "haps_hr.pdf")

```


## Looking up two samples for their concordance, jaccard index and root mean square error

```{r, eval = F}
# read in the data and convert into matrixes for easy look up
sampleNames = readr::read_tsv("measures/sampleNames.tab.txt", col_names = "sample")

ccc = bind_cols(
  sampleNames,
  readr::read_tsv("measures/ccc_on_targets_shared.tab.txt.gz", col_names = sampleNames$sample)
)
ccc_mat <- ccc |>
  tibble::column_to_rownames("sample") |>
  as.matrix()

jaccard = bind_cols(
  sampleNames,
  readr::read_tsv("measures/jaccard_on_targets_shared.tab.txt.gz", col_names = sampleNames$sample)
)
jaccard_mat <- jaccard |>
  tibble::column_to_rownames("sample") |>
  as.matrix()

rmse = bind_cols(
  sampleNames,
  readr::read_tsv("measures/rmse_on_targets_shared.tab.txt.gz", col_names = sampleNames$sample)
)
rmse_mat <- rmse |>
  tibble::column_to_rownames("sample") |>
  as.matrix()



lookup_sample1_name = "samp1"
lookup_sample2_name = "samp2"


create_dt(tibble(
  sample1 = lookup_sample1_name,
  sample2 = lookup_sample2_name,
  ccc = ccc_mat[lookup_sample1_name, lookup_sample2_name]
  jaccard = jaccard_mat[lookup_sample1_name, lookup_sample2_name]
  rmse = rmse_mat[lookup_sample1_name, lookup_sample2_name]
))


```

## Looking up the connections for a group

```{r, eval = F}
adjacency_list = readr::read_tsv("adjacency_list.tsv.gz")

lookup_group = 0

adjacency_list_group0 = adjacency_list %>%
  filter(group == lookup_group)

create_dt(adjacency_list_group0)
```
)QMD";

static const std::string static_run_app_on_clustering = R"RUNAPP(#!/usr/bin/env Rscript

# =============================================================================
# Haplotype clustering explorer - launcher
# =============================================================================
#
# This is a small Shiny app that mirrors the `static_qmd_on_clustering` report
# (haplotype rainbow + pairwise sample look-ups + per-group connections) but
# adds interactive controls: uploading extra sample meta, sorting by it,
# choosing which samples / group to look up, and setting the PDF output path.
#
# ---- One-time package installation ------------------------------------------
#   install.packages(c("optparse", "shiny", "DT", "readr", "dplyr", "tibble"))
#   # HaplotypeRainbows (GitHub):
#   # remotes::install_github("nickjhathaway/HaplotypeRainbows")
#
# ---- Running the app --------------------------------------------------------
#   # from this directory, pointing at an output dir produced by
#   # `elucidator cluster_samples_using_ccc_of_microhaps`:
#   ./run_app.R --data-dir /path/to/<input>_cluster_TODAY
#
#   # serve on all interfaces (e.g. on a remote server) on port 8080:
#   ./run_app.R --data-dir /path/to/results --host 0.0.0.0 --port 8080
#
#   # then browse to http://<server-ip>:8080
#
# If `--data-dir` is omitted it defaults to the current directory and can also
# be changed from inside the app.
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(shiny)
})

option_list <- list(
  make_option(
    c("-d", "--data-dir"),
    type = "character",
    default = ".",
    help = "Directory containing clusters.tsv, haps.tsv.gz, adjacency_list.tsv.gz and measures/ [default: %default]"
  ),
  make_option(
    c("-p", "--port"),
    type = "integer",
    default = 3838,
    help = "Port to run the Shiny app on [default: %default]"
  ),
  make_option(
    c("-H", "--host"),
    type = "character",
    default = "127.0.0.1",
    help = "Host address to serve on (use 0.0.0.0 to expose to the network) [default: %default]"
  ),
  make_option(
    c("-b", "--browser"),
    action = "store_true",
    default = FALSE,
    help = "Automatically open a browser"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

# app.R reads this to pre-populate the data directory input
Sys.setenv(HR_CLUSTERING_DATA_DIR = normalizePath(opt$`data-dir`, mustWork = FALSE))

app_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
if (length(app_dir) == 0 || app_dir == "") app_dir <- "."

shiny::runApp(
  appDir = file.path(app_dir, "app.R"),
  host = opt$host,
  port = opt$port,
  launch.browser = opt$browser
)
)RUNAPP";

static const std::string static_app_on_clustering = R"APPR(################################################################################
# Haplotype clustering explorer - minimal Shiny app
#
# Mirrors `static_qmd_on_clustering`:
#   1. Creating a haplotype rainbow (sorted/clustered by group), with optional
#      extra sample meta (uploaded tsv/csv -> update_sample_meta) and sorting
#      by that meta, plus PDF export to a user-set path.
#   2. Looking up two samples for CCC / Jaccard / RMSE.
#   3. Looking up the connections for a group in the adjacency list.
#
# Packages: shiny, HaplotypeRainbows, readr, dplyr, tibble, DT
################################################################################

suppressPackageStartupMessages({
  library(shiny)
  library(HaplotypeRainbows)
  library(readr)
  library(dplyr)
  library(tibble)
  library(DT)
})

# -- helpers -------------------------------------------------------------------

create_dt <- function(x) {
  DT::datatable(
    x,
    extensions = "Buttons",
    options = list(
      dom = "Blfrtip",
      buttons = c("copy", "csv", "excel", "pdf", "print"),
      lengthMenu = list(c(10, 25, 50, -1), c(10, 25, 50, "All")),
      scrollX = TRUE
    ),
    filter = "top",
    rownames = FALSE
  )
}

quick_summary <- function(x, probs = c(0.25, 0.5, 0.75)) {
  q <- quantile(x, probs = probs, na.rm = TRUE)
  tibble(
    n = length(x[!is.na(x)]),
    sd = sd(x, na.rm = TRUE),
    mean = mean(x, na.rm = TRUE),
    min  = min(x, na.rm = TRUE),
    max  = max(x, na.rm = TRUE),
    !!!setNames(as.list(q), paste0("q", probs*100))
  )
}

# original_name drives extension detection; path is what actually gets read
read_delim_auto <- function(path, original_name = path) {
  base_ext  <- tolower(tools::file_ext(original_name))
  inner_ext <- if (base_ext == "gz") {
    tolower(tools::file_ext(tools::file_path_sans_ext(original_name)))
  } else {
    base_ext
  }
  if (inner_ext %in% c("tsv", "tab", "txt")) {
    read_tsv(path, show_col_types = FALSE)
  } else {
    read_csv(path, show_col_types = FALSE)
  }
}

# read a square measure matrix written alongside measures/sampleNames.tab.txt
read_measure_matrix <- function(measures_dir, fnp) {
  sample_names <- read_tsv(
    file.path(measures_dir, "sampleNames.tab.txt"),
    col_names = "sample", show_col_types = FALSE
  )
  mat <- read_tsv(
    file.path(measures_dir, fnp),
    col_names = sample_names$sample, show_col_types = FALSE
  )
  bind_cols(sample_names, mat) %>%
    column_to_rownames("sample") %>%
    as.matrix()
}

default_data_dir <- {
  d <- Sys.getenv("HR_CLUSTERING_DATA_DIR", unset = ".")
  if (nzchar(d)) d else "."
}

# -- UI ------------------------------------------------------------------------

ui <- fluidPage(
  titlePanel("Haplotype clustering explorer"),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      textInput("data_dir", "Data directory", value = default_data_dir),
      helpText("Directory with clusters.tsv, haps.tsv.gz,",
               "adjacency_list.tsv.gz and measures/."),
      actionButton("load", "Load data", class = "btn-primary"),
      tags$hr(),
      uiOutput("load_status")
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",

        # ---- Tab 1: haplotype rainbow ------------------------------------------
        tabPanel(
          "Haplotype rainbow",
          br(),
          fluidRow(
            column(
              4,
              h4("Column mapping"),
              selectInput("sample_col", "Sample column", choices = NULL),
              selectInput("target_col", "Target column", choices = NULL),
              selectInput("popuid_col", "Pop UID (seq) column", choices = NULL),
              selectInput("rel_abund_col", "Relative abundance column", choices = NULL)
            ),
            column(
              4,
              h4("Extra sample meta (optional)"),
              fileInput("meta_file", "Upload tsv / csv", accept = c(".tsv", ".txt", ".tab", ".csv")),
              selectInput("meta_match_col", "Meta match column (sample id)", choices = NULL),
              helpText("Added via update_sample_meta().")
            ),
            column(
              4,
              h4("Sorting & clustering"),
              checkboxInput("abund_weighted", "Cluster abundance weighted", value = TRUE),
              selectizeInput(
                "sort_meta_cols", "Sort samples by meta (in order)",
                choices = "group", selected = "group", multiple = TRUE
              ),
              checkboxInput("sort_desc", "Sort descending", value = FALSE),
              checkboxInput("cluster_gaps", "Add sample cluster gaps", value = TRUE)
            )
          ),
          fluidRow(
            column(
              6,
              selectizeInput(
                "annot_cols", "Annotation columns (blank = all meta)",
                choices = NULL, selected = NULL, multiple = TRUE
              )
            ),
            column(
              6,
              textInput("pdf_path", "Output PDF path", value = "haps_hr.pdf"),
              numericInput("pdf_width", "PDF width (in, blank = auto)", value = NA),
              numericInput("pdf_height", "PDF height (in, blank = auto)", value = NA)
            )
          ),
          fluidRow(
            column(
              12,
              actionButton("render", "Render preview", class = "btn-primary"),
              actionButton("export_pdf", "Export PDF"),
              uiOutput("pdf_status")
            )
          ),
          br(),
          plotOutput("rainbow_plot", height = "700px")
        ),

        # ---- Tab 2: pairwise sample look-up ------------------------------------
        tabPanel(
          "Pairwise look-up",
          br(),
          fluidRow(
            column(6, selectInput("lookup_s1", "Sample 1", choices = NULL)),
            column(6, selectInput("lookup_s2", "Sample 2", choices = NULL))
          ),
          helpText("Concordance (CCC), Jaccard index and RMSE for the two samples."),
          DTOutput("pairwise_tbl")
        ),

        # ---- Tab 3: group connections ------------------------------------------
        tabPanel(
          "Group connections",
          br(),
          selectInput("lookup_group", "Group", choices = NULL),
          DTOutput("group_tbl"),
          br(),
          h4("Group summary"),
          selectInput("summary_measure", "Summarize measure",
                      choices = c("ccc", "jaccard", "rmse"), selected = "ccc"),
          DTOutput("group_summary_tbl")
        )
      )
    )
  )
)

# -- server --------------------------------------------------------------------

server <- function(input, output, session) {

  data_store <- reactiveValues(
    clusters = NULL, haps = NULL, adjacency = NULL,
    ccc = NULL, jaccard = NULL, rmse = NULL,
    sample_names = NULL, error = NULL
  )

  meta_store <- reactiveVal(NULL)

  observeEvent(input$load, {
    dir <- input$data_dir
    data_store$error <- NULL
    tryCatch({
      req_file <- function(f) {
        p <- file.path(dir, f)
        if (!file.exists(p)) stop(sprintf("Missing file: %s", p))
        p
      }
      clusters  <- read_tsv(req_file("clusters.tsv"), show_col_types = FALSE)
      haps      <- read_tsv(req_file("haps.tsv.gz"), show_col_types = FALSE)
      adjacency <- read_tsv(req_file("adjacency_list.tsv.gz"), show_col_types = FALSE)

      measures_dir <- file.path(dir, "measures")
      ccc     <- read_measure_matrix(measures_dir, "ccc_on_targets_shared.tab.txt.gz")
      jaccard <- read_measure_matrix(measures_dir, "jaccard_on_targets_shared.tab.txt.gz")
      rmse    <- read_measure_matrix(measures_dir, "rmse_on_targets_shared.tab.txt.gz")

      data_store$clusters     <- clusters
      data_store$haps         <- haps
      data_store$adjacency    <- adjacency
      data_store$ccc          <- ccc
      data_store$jaccard      <- jaccard
      data_store$rmse         <- rmse
      data_store$sample_names <- rownames(ccc)

      # column mapping defaults, matching the qmd where present
      hap_cols <- names(haps)
      pick <- function(preferred) if (preferred %in% hap_cols) preferred else hap_cols[1]
      updateSelectInput(session, "sample_col", choices = hap_cols, selected = pick("library_sample_name"))
      updateSelectInput(session, "target_col", choices = hap_cols, selected = pick("target_name"))
      updateSelectInput(session, "popuid_col", choices = hap_cols, selected = pick("seq"))
      updateSelectInput(session, "rel_abund_col", choices = hap_cols, selected = pick("within_sample_freq"))

      # meta sort / annotation start with the clustering group column
      cluster_cols <- setdiff(names(clusters), "sample")
      updateSelectizeInput(session, "sort_meta_cols", choices = cluster_cols,
                           selected = intersect("group", cluster_cols))
      updateSelectizeInput(session, "annot_cols", choices = cluster_cols, selected = character(0))

      updateSelectInput(session, "lookup_s1", choices = data_store$sample_names)
      updateSelectInput(session, "lookup_s2", choices = data_store$sample_names,
                        selected = if (length(data_store$sample_names) > 1) data_store$sample_names[2] else NULL)

      groups <- sort(unique(adjacency$group))
      updateSelectInput(session, "lookup_group", choices = groups)

      measure_cols <- intersect(c("ccc", "jaccard", "jaccard", "rmse", "targets_shared"),
                                names(adjacency))
      updateSelectInput(session, "summary_measure", choices = measure_cols,
                        selected = if ("ccc" %in% measure_cols) "ccc" else measure_cols[1])
    }, error = function(e) {
      data_store$error <- conditionMessage(e)
    })
  })

  output$load_status <- renderUI({
    if (!is.null(data_store$error)) {
      div(style = "color:#b00;", strong("Error: "), data_store$error)
    } else if (!is.null(data_store$clusters)) {
      div(style = "color:#080;",
          sprintf("Loaded %d samples, %d adjacency rows.",
                  length(data_store$sample_names), nrow(data_store$adjacency)))
    } else {
      helpText("No data loaded yet.")
    }
  })

  # uploaded extra meta -> refresh match column + sort/annotation choices
  observeEvent(input$meta_file, {
    meta <- read_delim_auto(input$meta_file$datapath, input$meta_file$name)
    meta_store(meta)
    meta_cols <- names(meta)
    guess_match <- if ("sample" %in% meta_cols) "sample" else meta_cols[1]
    updateSelectInput(session, "meta_match_col", choices = meta_cols, selected = guess_match)

    value_cols <- setdiff(meta_cols, guess_match)
    cluster_cols <- if (!is.null(data_store$clusters)) setdiff(names(data_store$clusters), "sample") else character(0)
    all_meta <- union(cluster_cols, value_cols)
    updateSelectizeInput(session, "sort_meta_cols", choices = all_meta,
                         selected = isolate(input$sort_meta_cols))
    updateSelectizeInput(session, "annot_cols", choices = all_meta,
                         selected = isolate(input$annot_cols))
  })

  # build the prepped/sorted haplotype rainbow object
  build_rainbow <- reactive({
    req(data_store$haps, data_store$clusters)

    hr <- haplotype_rainbow(
      data_store$haps,
      sample_col    = input$sample_col,
      target_col    = input$target_col,
      popuid_col    = input$popuid_col,
      rel_abund_col = input$rel_abund_col
    )$prep()

    hr$set_sample_meta(data_store$clusters, "sample")

    meta <- meta_store()
    if (!is.null(meta) && !is.null(input$meta_match_col) && nzchar(input$meta_match_col)) {
      hr$update_sample_meta(meta, input$meta_match_col)
    }

    hr$sort_samples_by_clustering(abundance_weighted = isTRUE(input$abund_weighted))

    sort_cols <- input$sort_meta_cols
    if (length(sort_cols) > 0) {
      hr$sort_samples_by_meta(sort_cols, desc = isTRUE(input$sort_desc))
    }

    if (isTRUE(input$cluster_gaps)) {
      hr$add_sample_cluster_gaps()
    }
    hr
  })

  annotated_plot <- reactive({
    hr <- build_rainbow()
    annot_cols <- input$annot_cols
    if (length(annot_cols) == 0) annot_cols <- NULL
    hr$add_sample_annotation_to_plot(hr$plot(), cols = annot_cols)
  })

  rainbow_plot_rv <- eventReactive(input$render, {
    annotated_plot()
  })

  output$rainbow_plot <- renderPlot({
    rainbow_plot_rv()
  })

  observeEvent(input$export_pdf, {
    output$pdf_status <- renderUI(helpText("Rendering PDF..."))
    tryCatch({
      hr <- build_rainbow()
      w <- if (is.na(input$pdf_width)) NULL else input$pdf_width
      h <- if (is.na(input$pdf_height)) NULL else input$pdf_height
      hr$save_pdf(annotated_plot(), input$pdf_path, width = w, height = h)
      output$pdf_status <- renderUI(
        div(style = "color:#080;", sprintf("Saved: %s", normalizePath(input$pdf_path, mustWork = FALSE)))
      )
    }, error = function(e) {
      output$pdf_status <- renderUI(div(style = "color:#b00;", conditionMessage(e)))
    })
  })

  # ---- pairwise look-up --------------------------------------------------------
  output$pairwise_tbl <- renderDT({
    req(data_store$ccc, input$lookup_s1, input$lookup_s2)
    s1 <- input$lookup_s1
    s2 <- input$lookup_s2
    create_dt(tibble(
      sample1 = s1,
      sample2 = s2,
      ccc     = data_store$ccc[s1, s2],
      jaccard = data_store$jaccard[s1, s2],
      rmse    = data_store$rmse[s1, s2]
    ))
  })

  # ---- group connections -------------------------------------------------------
  group_rows <- reactive({
    req(data_store$adjacency, input$lookup_group)
    grp <- input$lookup_group
    # match the group column's type when filtering
    if (is.numeric(data_store$adjacency$group)) grp <- as.numeric(grp)
    dplyr::filter(data_store$adjacency, group == grp)
  })

  output$group_tbl <- renderDT({
    create_dt(group_rows())
  })

  output$group_summary_tbl <- renderDT({
    req(input$summary_measure)
    create_dt(dplyr::summarise(group_rows(), quick_summary(.data[[input$summary_measure]])))
  })
}

shinyApp(ui, server)
)APPR";

int popGenExpRunner::cluster_samples_dist_of_microhaps_sharing(const njh::progutils::CmdArgs & inputCommands){
	double minimumLociCoverageToKeepSamples = 0.90;
  double concordance_cut_off = 0.95;
  double rmse_cut_off = 0.15;
  bool cluster_on_rmse = false;
  njhUndirWeightedGraph<double, std::vector<double>>::dbscanPars dbscanPars;
  // dbscanPars.eps_ = 0.50;
  dbscanPars.minEpNeighbors_ = 5;
	HapsEncodedMatrix::SetWithExternalPars pars;
  uint32_t pairwise_factor_bin_size = 1000;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
  setUp.setOption(pairwise_factor_bin_size, "--pairwise_factor_bin_size", "pairwise_factor_bin_size");
	setUp.setOption(minimumLociCoverageToKeepSamples, "--minimumLociCoverageToKeepSamples", "minimum Loci Coverage To Keep Samples in post analysis steps, must have reads for at least this frction of the total loci");
  setUp.setOption(concordance_cut_off, "--concordance_cut_off", "concordance cut off");
  dbscanPars.eps_ = 1 - concordance_cut_off;
  setUp.setOption(dbscanPars.minEpNeighbors_, "--min_neighbors", "The minimum neighbors for the DB scan clustering");
  setUp.setOption(rmse_cut_off, "--rmse_cut_off", "RMSE cut off");
  setUp.setOption(cluster_on_rmse, "--cluster_on_rmse", "cluster on RMSE");
  if (cluster_on_rmse) {
    dbscanPars.eps_ = rmse_cut_off;
  }

  pars.setDefaults(setUp);

  setUp.processDirectoryOutputName(bfs::path(bfs::basename(pars.tableFnp)).string() + "_cluster_TODAY", true);
	setUp.finishSetUp(std::cout);

	setUp.startARunLog(setUp.pars_.directoryName_);

	setUp.timer_.setLapName("initial");
	setUp.timer_.startNewLap("encode haplotypes");

  HapsEncodedMatrix haps(pars);
	setUp.timer_.startNewLap("get hap probabilities");
	haps.calcHapProbs();
	setUp.timer_.startNewLap("add relative abundances");
  haps.add_relative_abundance();
  setUp.timer_.startNewLap("get index measures");
  auto indexRes = haps.genIndexMeasures(setUp.pars_.verbose_);
  setUp.timer_.startNewLap("calc rmse and ccc");
	auto original_measures = haps.calc_ccc_rmse_measures(pairwise_factor_bin_size, setUp.pars_.verbose_);
  setUp.timer_.startNewLap("writing output matrices");
  auto measures_dir = njh::files::make_path(setUp.pars_.directoryName_, "measures/");
  njh::files::makeDir(measures_dir);
  //write out measures
	OutputStream outSampNamesOut(njh::files::make_path(measures_dir, "sampleNames.tab.txt"));
	outSampNamesOut << njh::conToStr(haps.sampNamesVec_, "\n") << std::endl;
  {
	  auto ccc_out_fnp = njh::files::make_path(measures_dir, "ccc_on_targets_shared.tab.txt.gz");
	  OutputStream ccc_out(ccc_out_fnp);
	  for(const auto & ccc_row : original_measures.ccc){
	    ccc_out << njh::conToStr(ccc_row, "\t") << std::endl;
	  }
  }
  {
	  auto rmse_out_fnp = njh::files::make_path(measures_dir, "rmse_on_targets_shared.tab.txt.gz");
	  OutputStream rmse_out(rmse_out_fnp);
	  for(const auto & rmse_row : original_measures.rmse){
	    rmse_out << njh::conToStr(rmse_row, "\t") << std::endl;
	  }
  }
  {
	  auto jaccard_out_fnp = njh::files::make_path(measures_dir, "jaccard_on_targets_shared.tab.txt.gz");
	  OutputStream jaccard_out(jaccard_out_fnp);
	  for(const auto & jaccard_row : indexRes.byHapsTarShared){
	    jaccard_out << njh::conToStr(jaccard_row, "\t") << std::endl;
	  }
  }
  {
	  auto targets_shared_out_fnp = njh::files::make_path(measures_dir, "targets_shared.tab.txt.gz");
	  OutputStream targets_shared_out(targets_shared_out_fnp);
	  for(const auto & targets_shared_row : original_measures.targets_shared){
	    targets_shared_out << njh::conToStr(targets_shared_row, "\t") << std::endl;
	  }
  }
  setUp.timer_.startNewLap("getting loci coverage info");

	std::unordered_map<std::string, double> lociCoveragePerSample = haps.getTargetCoveragePerSample();
	{
		table numTargetsPerSample = haps.getTableNumberTargetsPerSample(minimumLociCoverageToKeepSamples);
		OutputStream lociCoverageOut(njh::files::make_path(setUp.pars_.directoryName_, "loci_coverage_per_sample_info.tsv"));
		numTargetsPerSample.outPutContents(lociCoverageOut, "\t");
	}

  setUp.timer_.setLapName("transforming matrix");
  auto inverse_ccc = original_measures.ccc;
	{
	  //for the distance functions below to work, have to transform CCC so that the lower the better, CCC runs from -1 to 1, so below will transform it so it runs from 0 to 2 with 0 being CCC of 1, 1 being CCC 0, and 2 being CCC -2
	  PairwisePairFactory pairFactory(inverse_ccc.size());
	  uint32_t pairBatchCount = 100000;
	  std::function<void()> transform_ccc =
    [&pairFactory,
      &pairBatchCount,
      &inverse_ccc]() {
      PairwisePairFactory::PairwisePairVec pairs;
      while (pairFactory.setNextPairs(pairs, pairBatchCount)) {
        for (const auto & pair : pairs.pairs_) {
          inverse_ccc[pair.row_][pair.col_] = -1 * (inverse_ccc[pair.row_][pair.col_] - 1);
          inverse_ccc[pair.col_][pair.row_] = inverse_ccc[pair.row_][pair.col_];
        }
      }
    };
	  njh::concurrent::runVoidFunctionThreaded(transform_ccc, pars.numThreads);
    // fill the diagonal
	  for (uint32_t pos = 0; pos < inverse_ccc.size(); ++pos) {
	    inverse_ccc[pos][pos] = 0;
	  }
	}

  const std::vector<std::vector<double> > &distance_matrix = cluster_on_rmse ? original_measures.rmse : inverse_ccc;

  setUp.timer_.setLapName("building matrix");
  auto dist_graph = std::make_unique<njhUndirWeightedGraph<double, std::vector<double> > > ();
  for (const auto & pos : iter::range(distance_matrix.size())) {
    dist_graph->addNode(estd::to_string(pos), distance_matrix[pos]);
  }
	{
	  uint32_t belowEp = 0;
	  PairwisePairFactory pairFactory(distance_matrix.size());
	  uint32_t pairBatchCount = 100000;
	  std::mutex graphMut;
	  struct PairDist {
	    PairDist(const PairwisePairFactory::PairwisePair & pair, double dist) :
          pair_(pair), dist_(dist) {
	    }
	    PairwisePairFactory::PairwisePair pair_;
	    double dist_;
	  };

    std::function<void()> addToGraph =
        [&graphMut, &pairFactory,&pairBatchCount,&belowEp,
          &distance_matrix, &dbscanPars,
          &haps, &dist_graph,
          &lociCoveragePerSample, &minimumLociCoverageToKeepSamples]() {
      PairwisePairFactory::PairwisePairVec pairs;
      std::vector<PairDist> belowEps;
      while (pairFactory.setNextPairs(pairs, pairBatchCount)) {
        for (const auto &pair: pairs.pairs_) {
          if (lociCoveragePerSample[haps.sampNamesVec_[pair.row_]] < minimumLociCoverageToKeepSamples ||
              lociCoveragePerSample[haps.sampNamesVec_[pair.col_]] < minimumLociCoverageToKeepSamples) {
            continue;
          }
          auto dist = distance_matrix[pair.row_][pair.col_];
          if (dist <= dbscanPars.eps_) {
            belowEps.emplace_back(PairDist{pair, dist});
          }
        }
      }
      if (!belowEps.empty()) {
        std::lock_guard<std::mutex> lock(graphMut);
        belowEp += belowEps.size();
        for (const auto &bEps: belowEps) {
          dist_graph->addEdge(estd::to_string(bEps.pair_.row_),
                              estd::to_string(bEps.pair_.col_),
                              bEps.dist_);
        }
      }
    };
	  njh::concurrent::runVoidFunctionThreaded(addToGraph, pars.numThreads);
	}

  setUp.timer_.startNewLap("dbscan");
  dist_graph->dbscan(dbscanPars);
  if (dbscanPars.minEpNeighbors_ > 2) {
    dist_graph->add_small_groups_in_off_nodes(dbscanPars);
  }
  if (setUp.pars_.verbose_) {
    std::cout << "Determined " << dist_graph->numberOfGroups_ << " groups" << std::endl;
  }


  setUp.timer_.startNewLap("output");
  OutputStream outFile(OutOptions(njh::files::make_path(setUp.pars_.directoryName_, "clusters.tsv")));
  outFile << "sample\tgroup";
  outFile << std::endl;

  std::map<uint32_t, std::vector<uint32_t>> groupIndexes;
  std::vector<uint32_t> allGroupedIndices;
  for(const auto & n : iter::enumerate(dist_graph->nodes_)) {
    groupIndexes[n.element->group_].emplace_back(n.index);
    if (std::numeric_limits<uint32_t>::max() != n.element->group_) {
      allGroupedIndices.emplace_back(n.index);
    }
  }
  std::unordered_map<uint32_t, std::map<std::string, double> > groups_ccc_stats;
  std::unordered_map<uint32_t, std::map<std::string, double> > groups_rmse_stats;
  std::unordered_map<uint32_t, std::map<std::string, double> > groups_jaccard_stats;
  for (const auto &group: groupIndexes) {
    std::vector<double> cccsWithinGroup;
    std::vector<double> rmsesWithinGroup;
    std::vector<double> jaccardsWithinGroup;
    {
      PairwisePairFactory pfac(group.second.size());
      PairwisePairFactory::PairwisePair pair;
      while (pfac.setNextPair(pair)) {
        //have to transform back due to the previous transform
        //cccsWithinGroup.emplace_back(measures.ccc[group.second[pair.col_]][group.second[pair.row_]] * -1 + 1);
        cccsWithinGroup.emplace_back(original_measures.ccc[group.second[pair.col_]][group.second[pair.row_]]);
        rmsesWithinGroup.emplace_back(original_measures.rmse[group.second[pair.col_]][group.second[pair.row_]]);
        jaccardsWithinGroup.emplace_back(indexRes.byHapsTarShared[group.second[pair.col_]][group.second[pair.row_]]);
      }
    }
    groups_ccc_stats[group.first] = getStatsOnVec(cccsWithinGroup);
    groups_rmse_stats[group.first] = getStatsOnVec(rmsesWithinGroup);
    groups_jaccard_stats[group.first] = getStatsOnVec(jaccardsWithinGroup);

    for (const auto &idx: group.second) {
      if (lociCoveragePerSample[haps.sampNamesVec_[idx]] < minimumLociCoverageToKeepSamples) {
        outFile << haps.sampNamesVec_[idx] << "\t" << "low_coverage_not_clustered";
      } else if (group.first == std::numeric_limits<uint32_t>::max()) {
        outFile << haps.sampNamesVec_[idx] << "\t" << "nogroup";
      } else {
        outFile << haps.sampNamesVec_[idx] << "\t" << group.first;
      }
      outFile << std::endl;
    }
  }

  OutputStream outGroupCountsFile(
    OutOptions(njh::files::make_path(setUp.pars_.directoryName_, "clusters_groupCounts.tsv")));
  outGroupCountsFile << "group\tsample_count";
  outGroupCountsFile << "\tmin_ccc\tmedian_ccc\tmean_ccc\tmax_ccc";
  outGroupCountsFile << "\tmin_rmse\tmedian_rmse\tmean_rmse\tmax_rmse";
  outGroupCountsFile << "\tmin_jaccard\tmedian_jaccard\tmean_jaccard\tmax_jaccard";
  outGroupCountsFile << std::endl;
  for (const auto &group: groupIndexes) {
    if (group.first == std::numeric_limits<uint32_t>::max()) {
      uint32_t low_coverage_not_clustered_cnt = 0;
      for (const auto &samp_idx: group.second) {
        if (lociCoveragePerSample[haps.sampNamesVec_[samp_idx]] < minimumLociCoverageToKeepSamples) {
          low_coverage_not_clustered_cnt++;
        }
      }
      outGroupCountsFile << "nogroup" << "\t" << group.second.size() - low_coverage_not_clustered_cnt;
      outGroupCountsFile << "\t" << "NA"
          << "\t" << "NA"
          << "\t" << "NA"
          << "\t" << "NA";
      outGroupCountsFile << "\t" << "NA"
          << "\t" << "NA"
          << "\t" << "NA"
          << "\t" << "NA";
      outGroupCountsFile << "\t" << "NA"
          << "\t" << "NA"
          << "\t" << "NA"
          << "\t" << "NA";
      outGroupCountsFile << std::endl;

      outGroupCountsFile << "low_coverage_not_clustered" << "\t" << low_coverage_not_clustered_cnt;
      outGroupCountsFile << "\t" << "NA"
          << "\t" << "NA"
          << "\t" << "NA"
          << "\t" << "NA";
      outGroupCountsFile << "\t" << "NA"
          << "\t" << "NA"
          << "\t" << "NA"
          << "\t" << "NA";
      outGroupCountsFile << "\t" << "NA"
          << "\t" << "NA"
          << "\t" << "NA"
          << "\t" << "NA";
      outGroupCountsFile << std::endl;
    } else {
      outGroupCountsFile << group.first << "\t" << group.second.size();
      outGroupCountsFile << "\t" << groups_ccc_stats[group.first]["min"]
          << "\t" << groups_ccc_stats[group.first]["median"]
          << "\t" << groups_ccc_stats[group.first]["mean"]
          << "\t" << groups_ccc_stats[group.first]["max"];
      outGroupCountsFile << "\t" << groups_rmse_stats[group.first]["min"]
          << "\t" << groups_rmse_stats[group.first]["median"]
          << "\t" << groups_rmse_stats[group.first]["mean"]
          << "\t" << groups_rmse_stats[group.first]["max"];
      outGroupCountsFile << "\t" << groups_jaccard_stats[group.first]["min"]
          << "\t" << groups_jaccard_stats[group.first]["median"]
          << "\t" << groups_jaccard_stats[group.first]["mean"]
          << "\t" << groups_jaccard_stats[group.first]["max"];
      outGroupCountsFile << std::endl;
    }
  }

  OutputStream adjacency_list_out(njh::files::make_path(setUp.pars_.directoryName_, "adjacency_list.tsv.gz"));
  adjacency_list_out << "node1\tnode2\tgroup\ttargets_shared\tccc\trmse\tjaccard" << std::endl;
  for (const auto & e : dist_graph->edges_) {
    if (e->on_) {
      auto node1 = e->nodeToNode_.begin()->second.lock();
      auto node2 = e->nodeToNode_.rbegin()->second.lock();
      auto node1_idx = njh::StrToNumConverter::stoToNum<uint32_t>(node1->name_);
      auto node2_idx = njh::StrToNumConverter::stoToNum<uint32_t>(node2->name_);
      adjacency_list_out << haps.sampNamesVec_[node1_idx]
        << "\t" << haps.sampNamesVec_[node2_idx]
        << "\t" << node1->group_
        << "\t" << original_measures.targets_shared[node1_idx][node2_idx]
        << "\t" << original_measures.ccc[node1_idx][node2_idx]
        << "\t" << original_measures.rmse[node1_idx][node2_idx]
        << "\t" << indexRes.byHapsTarShared[node1_idx][node2_idx]
        << std::endl;
    }
  }

	haps.exportEncodedTable().outPutContents(TableIOOpts::genTabFileOut(njh::files::make_path(setUp.pars_.directoryName_, "haps.tsv.gz")));
	OutputStream static_qmd_out(njh::files::make_path(setUp.pars_.directoryName_, "process.qmd"));
	static_qmd_out << static_qmd_on_clustering << std::endl;
	//write out a minimal shiny app (run_app.R + app.R) alongside the outputs so it can be run in place
	{
		auto run_app_fnp = njh::files::make_path(setUp.pars_.directoryName_, "run_app.R");
		OutputStream run_app_out(run_app_fnp);
		run_app_out << static_run_app_on_clustering << std::endl;
		::chmod(run_app_fnp.c_str(),
		S_IWUSR | S_IRUSR | S_IXUSR | S_IRGRP | S_IWGRP | S_IXGRP | S_IROTH);
	}
	{
		OutputStream app_out(njh::files::make_path(setUp.pars_.directoryName_, "app.R"));
		app_out << static_app_on_clustering << std::endl;
	}
	setUp.timer_.logLapTimes(setUp.rLog_.runLogFile_, true, 6, true);
	return 0;

}




int popGenExpRunner::doPairwiseComparisonOnHapsSharing(const njh::progutils::CmdArgs & inputCommands){
	bool writeOutDistMatrices = false;
	bool clusterOnJaccardIndexShared = false;
	njhUndirWeightedGraph<double, std::shared_ptr<BasicPointMatrix<double>::BasicPoint>>::dbscanPars dbscanPars;
	// dbscanPars.eps_ = 0.50;
	dbscanPars.eps_ = 0.10;
	dbscanPars.minEpNeighbors_ = 2;
	double minimumLociCoverageToKeepSamples = 0.90;
	bfs::path metaFnp;
	VecStr metaFieldsToCalcPopDiffs{};
	HapsEncodedMatrix::SetWithExternalPars pars;
	bool writeOutTarsAbsoluteShared = false;
	bool doNotBreakWithRmse = false;
	double rmseCutOffToBreak = 0.10;
	bool doNotWriteOutGroupedRMSEs = false;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(minimumLociCoverageToKeepSamples, "--minimumLociCoverageToKeepSamples", "minimum Loci Coverage To Keep Samples in post analysis steps, must have reads for at least this frction of the total loci");

	setUp.setOption(clusterOnJaccardIndexShared, "--clusterOnJaccardIndexShared", "cluster On Jaccard Index Shared");
	setUp.setOption(doNotBreakWithRmse, "--doNotBreakWithRmse", "do Not Break With Rmse");
	setUp.setOption(rmseCutOffToBreak, "--rmseCutOffToBreak", "rmse Cut Off To Break");
	setUp.setOption(doNotWriteOutGroupedRMSEs, "--doNotWriteOutGroupedRMSEs", "wriet Out Grouped RMSEs");
	bool writeOutGroupedRMSEs = !doNotWriteOutGroupedRMSEs;

	setUp.setOption(dbscanPars.eps_, "--eps", "Epsilon (distance sensitivity of algorithm)");
	setUp.setOption(dbscanPars.minEpNeighbors_, "--minpts", "The minimum number of epsilon neighbors");
	setUp.setOption(writeOutDistMatrices, "--writeOutDistMatrices", "write Out Dist Matrices");
	setUp.setOption(metaFnp, "--metaFnp", "Table of meta data for samples, needs a column named sample and each additional column will be the meta data associated with that sample");
	setUp.setOption(metaFieldsToCalcPopDiffs, "--metaFieldsToCalcPopDiffs", "Meta Fields To Calc Pop Diffs");
	setUp.setOption(writeOutTarsAbsoluteShared, "--writeOutTarsAbsoluteShared", "write Out Tars Absolute Shared");

  pars.setDefaults(setUp);

  setUp.processDirectoryOutputName(bfs::path(bfs::basename(pars.tableFnp)).string() + "_doPairwiseComparisonOnHapsSharing_TODAY", true);
	setUp.finishSetUp(std::cout);

	setUp.startARunLog(setUp.pars_.directoryName_);


	setUp.timer_.setLapName("initial");
	setUp.timer_.startNewLap("encode haplotypes");
  HapsEncodedMatrix haps(pars);

  if (!metaFnp.empty()) {
    haps.addMeta(metaFnp);
    if (!metaFieldsToCalcPopDiffs.empty()) {
      haps.meta_->checkForFieldsThrow(metaFieldsToCalcPopDiffs);
    }
  } else if (!metaFieldsToCalcPopDiffs.empty()) {
    haps.addMetaWithInputTab(njh::vecToSet(metaFieldsToCalcPopDiffs));
  }
	setUp.timer_.startNewLap("get hap probabilities");
	/**@todo look into whether or not this is being actually used */
	haps.calcHapProbs();

	haps.add_relative_abundance();
	setUp.timer_.startNewLap("writing sample info");


	setUp.timer_.startNewLap("get index measures");
	auto indexRes = haps.genIndexMeasures(setUp.pars_.verbose_);
	OutputStream outSampNamesOut(njh::files::make_path(setUp.pars_.directoryName_, "sampleNames.tab.txt"));
	outSampNamesOut << njh::conToStr(haps.sampNamesVec_, "\n") << std::endl;

	if(writeOutDistMatrices){
		OutputStream byTargetOut(njh::files::make_path(setUp.pars_.directoryName_, "percOfTarSharingAtLeastOneHap.tab.txt.gz"));
		OutputStream byHapOut(njh::files::make_path(setUp.pars_.directoryName_, "jaccardByAllHap.tab.txt.gz"));
		OutputStream byHapTarSharedOut(njh::files::make_path(setUp.pars_.directoryName_, "jaccardByHapsTarShared.tab.txt.gz"));
		OutputStream avgHapOut(njh::files::make_path(setUp.pars_.directoryName_, "avgJaccardPerTarget.tab.txt.gz"));

		OutputStream byHapTarSharedWeightedOut(njh::files::make_path(setUp.pars_.directoryName_, "jaccardByHapsTarSharedWeighted.tab.txt.gz"));
		OutputStream avgHapWeightedOut(njh::files::make_path(setUp.pars_.directoryName_, "avgJaccardPerTargetWeighted.tab.txt.gz"));
		OutputStream targetsSharedBetweenSampsOut(njh::files::make_path(setUp.pars_.directoryName_, "targetsSharedBetweenSamps.tab.txt.gz"));



		for(const auto & it : indexRes.byTarget){
			byTargetOut << njh::conToStr(it, "\t") << std::endl;
		}
		for(const auto & ih : indexRes.byHapsTarShared){
			byHapTarSharedOut << njh::conToStr(ih, "\t") << std::endl;
		}
		for(const auto & ih : indexRes.byAllHaps){
			byHapOut << njh::conToStr(ih, "\t") << std::endl;
		}
		for(const auto & ih : indexRes.avgJacard){
			avgHapOut << njh::conToStr(ih, "\t") << std::endl;
		}
		for(const auto & ih : indexRes.byHapsTarSharedWeighted){
			byHapTarSharedWeightedOut << njh::conToStr(ih, "\t") << std::endl;
		}
		for(const auto & ih : indexRes.avgJacardWeighted){
			avgHapWeightedOut << njh::conToStr(ih, "\t") << std::endl;
		}

		for(const auto & ih : indexRes.targetsShared){
			targetsSharedBetweenSampsOut << njh::conToStr(ih, "\t") << std::endl;
		}
	}

	std::unordered_map<std::string, double> lociCoveragePerSample = haps.getTargetCoveragePerSample();

	{
		table numTargetsPerSample = haps.getTableNumberTargetsPerSample(minimumLociCoverageToKeepSamples);
		OutputStream lociCoverageOut(njh::files::make_path(setUp.pars_.directoryName_, "loci_coverage_per_sample_info.tsv"));
		numTargetsPerSample.outPutContents(lociCoverageOut, "\t");
	}

	if(clusterOnJaccardIndexShared) {
		auto distFnp = njh::files::make_path(setUp.pars_.directoryName_, "1MinusjaccardByHapsTarShared.tab.txt.gz");
		{
			OutputStream byHapTarSharedOut(distFnp);
			for(const auto & ih : indexRes.byHapsTarShared){
				std::vector<double> outRow;
				outRow.reserve(ih.size());
				for(const auto & i : ih) {
					outRow.emplace_back(1 - i);
				}
				byHapTarSharedOut << njh::conToStr(outRow, "\t") << std::endl;
			}
		}



		OutputStream outFile(OutOptions(njh::files::make_path(setUp.pars_.directoryName_, "clusters_by_jaccardTargetsShared.tsv")));
		njh::stopWatch watch;
		watch.setLapName("Reading in");
		auto mat = BasicPointMatrix<double>::readInBasicMatrix(distFnp, dbscanPars);
		watch.startNewLap("Adding nodes");
		std::vector<std::vector<double>> pairwiseRMSEs;
	  std::vector<std::vector<double>> pairwise_cccs;
		if(doNotBreakWithRmse) {
			// mat.setGraph(pars.numThreads, setUp.pars_.verbose_);
			mat.graph_ = std::make_unique<
					njhUndirWeightedGraph<double,
							std::shared_ptr<BasicPointMatrix<double>::BasicPoint>>>();

			for (const auto & pos : iter::range(mat.points_.size())) {
				mat.graph_->addNode(estd::to_string(pos), mat.points_[pos]);
			}

			uint32_t belowEp = 0;
			/**@todo this appears to be actually fairly slow, i think it's mostly because the eu calculations is so fast, perhaps a better way of multithreading this can be done
			 *
			 */
			PairwisePairFactory pairFactory(mat.points_.size());
			uint32_t pairBatchCount = 100000;
			std::mutex graphMut;
			struct PairDist {
				PairDist(const PairwisePairFactory::PairwisePair & pair, double dist) :
						pair_(pair), dist_(dist) {
				}
				PairwisePairFactory::PairwisePair pair_;
				double dist_;
			};

			std::function<void()> addToGraph =
					[&graphMut, &pairFactory,&pairBatchCount,&belowEp,
						&mat, &dbscanPars,
						&haps,
						&lociCoveragePerSample, &minimumLociCoverageToKeepSamples]() {
						PairwisePairFactory::PairwisePairVec pairs;
						std::vector<PairDist> belowEps;
						while(pairFactory.setNextPairs(pairs, pairBatchCount)) {
							for(const auto & pair : pairs.pairs_) {
								if (lociCoveragePerSample[haps.sampNamesVec_[pair.row_]] < minimumLociCoverageToKeepSamples ||
									lociCoveragePerSample[haps.sampNamesVec_[pair.col_]] < minimumLociCoverageToKeepSamples) {
									continue;
								}
								auto dist = mat.points_[pair.row_]->vals_[pair.col_];
								if (dist < dbscanPars.eps_) {
									belowEps.emplace_back(PairDist{pair, dist});
								}
							}
						}
						if(!belowEps.empty()) {
							std::lock_guard<std::mutex> lock(graphMut);
							belowEp += belowEps.size();
							for(const auto & bEps : belowEps) {
								mat.graph_->addEdge(estd::to_string(bEps.pair_.row_), estd::to_string(bEps.pair_.col_),
										bEps.dist_);
							}
						}
			};
			njh::concurrent::runVoidFunctionThreaded(addToGraph, pars.numThreads);
		} else {
			pairwiseRMSEs = std::vector<std::vector<double>>(haps.sampNames_.size(), std::vector<double>(haps.sampNames_.size(),1.0));
		  pairwise_cccs = std::vector<std::vector<double>>(haps.sampNames_.size(), std::vector<double>(haps.sampNames_.size(),0.0));
		  for(size_t pos = 0; pos < haps.sampNames_.size(); ++pos){
				//set diagonal
				pairwiseRMSEs[pos][pos] = 0;
		    pairwise_cccs[pos][pos] = 1.0;
			}

			if (0 == pars.numThreads) {
				pars.numThreads = 1;
			}

			mat.graph_ = std::make_unique<
					njhUndirWeightedGraph<double,
							std::shared_ptr<BasicPointMatrix<double>::BasicPoint>>>();

			for (const auto & pos : iter::range(mat.points_.size())) {
				mat.graph_->addNode(estd::to_string(pos), mat.points_[pos]);
			}

			uint32_t belowEp = 0;
			/**@todo this appears to be actually fairly slow, i think it's mostly because the eu calculations is so fast, perhaps a better way of multithreading this can be done
			 *
			 */
			PairwisePairFactory pairFactory(mat.points_.size());
			uint32_t pairBatchCount = 100000;
			std::mutex graphMut;
			struct PairDist {
				PairDist(const PairwisePairFactory::PairwisePair & pair, double dist) :
						pair_(pair), dist_(dist) {
				}
				PairwisePairFactory::PairwisePair pair_;
				double dist_;
			};

			std::function<void()> addToGraph =
					[&graphMut, &pairFactory,&pairBatchCount,&belowEp,&mat, &haps,&rmseCutOffToBreak,
						&pairwiseRMSEs, &pairwise_cccs, &lociCoveragePerSample,
						&minimumLociCoverageToKeepSamples]() {
						PairwisePairFactory::PairwisePairVec pairs;
						std::vector<PairDist> belowEps;
						while(pairFactory.setNextPairs(pairs, pairBatchCount)) {
							for (const auto &pair: pairs.pairs_) {
								if (lociCoveragePerSample[haps.sampNamesVec_[pair.row_]] < minimumLociCoverageToKeepSamples ||
								    lociCoveragePerSample[haps.sampNamesVec_[pair.col_]] < minimumLociCoverageToKeepSamples) {
									continue;
								}
								//auto dist = mat.points_[pair.row_]->euDist(*mat.points_[pair.col_]);
								auto dist = mat.points_[pair.row_]->vals_[pair.col_];
								if (dist < mat.dbscanPars_.eps_) {
								  std::vector<double> row_values;
								  std::vector<double> col_values;
									std::vector<double> rmses;
									double sum = 0;
									for(const auto tpos : iter::range(haps.tarNamesVec_.size())) {
										if(haps.targetsEncodeBySamp_[pair.col_][tpos]  + haps.targetsEncodeBySamp_[pair.row_][tpos] == 2) {
											double current_sum = 0;
											for(const auto hapPos : iter::range(haps.numberOfHapsPerTarget_[tpos])) {
												current_sum += std::pow(haps.hapsEncodeBySampRelAbund_[pair.col_][haps.tarStart_[tpos] + hapPos] - haps.hapsEncodeBySampRelAbund_[pair.row_][haps.tarStart_[tpos] + hapPos],2);
												sum +=         std::pow(haps.hapsEncodeBySampRelAbund_[pair.col_][haps.tarStart_[tpos] + hapPos] - haps.hapsEncodeBySampRelAbund_[pair.row_][haps.tarStart_[tpos] + hapPos],2);
											  row_values.emplace_back(haps.hapsEncodeBySampRelAbund_[pair.row_][haps.tarStart_[tpos] + hapPos]);
											  col_values.emplace_back(haps.hapsEncodeBySampRelAbund_[pair.col_][haps.tarStart_[tpos] + hapPos]);
											}
											rmses.emplace_back(std::sqrt(current_sum));
										}
									}
									//only add if mean RMSE is less than the cut off
									//if(vectorMean(rmses) < rmseCutOffToBreak) {
								  auto ccc_calc = ConcordanceCalculator::lins_ccc_with_ci(row_values, col_values);
									pairwiseRMSEs[pair.col_][pair.row_] = std::sqrt(sum/rmses.size());
									pairwiseRMSEs[pair.row_][pair.col_] = std::sqrt(sum/rmses.size());
								  pairwise_cccs[pair.col_][pair.row_] = ccc_calc.ccc;
								  pairwise_cccs[pair.row_][pair.col_] = ccc_calc.ccc;
									// std::cout << __FILE__ << " : " << __LINE__ << std::endl;
									// std::cout << "lociCoveragePerSample[haps.sampNamesVec_[pair.row_]]: " << lociCoveragePerSample[haps.sampNamesVec_[pair.row_]] << std::endl;
									// std::cout << "lociCoveragePerSample[haps.sampNamesVec_[pair.col_]]: " << lociCoveragePerSample[haps.sampNamesVec_[pair.col_]] << std::endl;
									// std::cout << "haps.sampNamesVec_[pair.row_]: " << haps.sampNamesVec_[pair.row_] << std::endl;
									// std::cout << "haps.sampNamesVec_[pair.col_]: " << haps.sampNamesVec_[pair.col_] << std::endl;
									// std::cout << "rmses.size(): " << rmses.size() << std::endl;
									// std::cout << "std::sqrt(sum/rmses.size()): " << std::sqrt(sum/rmses.size()) << std::endl;
									// std::cout << "rmseCutOffToBreak          : " << rmseCutOffToBreak << std::endl << std::endl;

									if(std::sqrt(sum/rmses.size()) < rmseCutOffToBreak){
										belowEps.emplace_back(PairDist{pair, dist});
									}
								}
							}
						}
						if(!belowEps.empty()) {
							std::lock_guard<std::mutex> lock(graphMut);
							belowEp += belowEps.size();
							for(const auto & bEps : belowEps) {
								mat.graph_->addEdge(estd::to_string(bEps.pair_.row_), estd::to_string(bEps.pair_.col_),
										bEps.dist_);
							}
						}
			};

			njh::concurrent::runVoidFunctionThreaded(addToGraph, pars.numThreads);

			if (setUp.pars_.verbose_) {
				std::cout << std::endl;
				std::cout << "below: " << belowEp << "/" << pairFactory.totalCompares_ << std::endl;
			}


		}


		watch.startNewLap("dbscan");
		mat.graph_->dbscan(dbscanPars);
		// if(!doNotBreakWithRmse){
		// 						//first fill the relative abundance vector with the input relative abundance
		// 	std::vector<std::vector<double>> hapsEncodeBySampRelAbund = std::vector<std::vector<double>> (haps.sampNames_.size());
		// 	for(const auto & samp : haps.sampNames_){
		// 		hapsEncodeBySampRelAbund[haps.sampNamesKey_[samp]] = std::vector<double>(haps.totalHaps_, 0);
		// 	}
		// 	TableReader reReadHapTab(TableIOOpts(InOptions(haps.pars_.tableFnp), "\t", true));
		// 	VecStr row;
		// 	while(reReadHapTab.getNextRow(row)){
		// 		const auto& samp = row[reReadHapTab.header_.getColPos(haps.pars_.sampleCol)];
		// 		const auto& tar = row[reReadHapTab.header_.getColPos(haps.pars_.targetNameCol)];
		// 		if(!haps.pars_.selectSamples.empty() && !njh::in(samp, haps.pars_.selectSamples)){
		// 			continue;
		// 		}
		// 		if(!haps.pars_.selectTargets.empty() && !njh::in(tar, haps.pars_.selectTargets)){
		// 			continue;
		// 		}
		// 		const auto& hapName = row[reReadHapTab.header_.getColPos(haps.pars_.popIDCol)];
		// 		auto rBund = njh::StrToNumConverter::stoToNum<double>(row[reReadHapTab.header_.getColPos(haps.pars_.relAbundCol)]);
		// 		auto tKey = haps.tarNameKey_[tar];
		// 		auto hKey = haps.hapNamesKey_[tar][hapName];
		// 		//				if(rBund > 0 && rBund < 1){
		// 		//					std::cout << rBund << std::endl;
		// 		//				}
		// 		hapsEncodeBySampRelAbund[haps.sampNamesKey_[samp]][haps.tarStart_[tKey] + hKey] = rBund;
		// 	}
		// 	//now recalculate the relative abundance to be 0-1
		// 	for(const auto pos : iter::range(haps.sampNames_.size())) {
		// 		for(const auto tpos : iter::range(haps.tarNamesVec_.size())) {
		// 			if(haps.targetsEncodeBySamp_[pos][tpos] == 1) {
		// 				double sum = 0;
		// 				for(const auto hapPos : iter::range(haps.numberOfHapsPerTarget_[tpos])) {
		// 					sum += hapsEncodeBySampRelAbund[pos][haps.tarStart_[tpos] + hapPos];
		// 				}
		// 				for(const auto hapPos : iter::range(haps.numberOfHapsPerTarget_[tpos])) {
		// 					hapsEncodeBySampRelAbund[pos][haps.tarStart_[tpos] + hapPos] = hapsEncodeBySampRelAbund[pos][haps.tarStart_[tpos] + hapPos]/sum;
		// 				}
		// 			}
		// 		}
		// 	}
		// 	// for(const auto pos : iter::range(haps.sampNamesKey_.size())) {
		// 	// 	std::cout << njh::conToStr( hapsEncodeBySampRelAbund[pos], "\t") << std::endl;
		// 	// }
		// 	std::map<uint32_t, std::vector<uint32_t>> groupIndexes;
		// 	for(const auto & n : iter::enumerate(mat.graph_->nodes_)) {
		// 		groupIndexes[n.element->group_].emplace_back(n.index);
		// 	}
		// 	for(const auto & group : groupIndexes) {
		// 		PairwisePairFactory pair_factory(group.second.size());
		// 		PairwisePairFactory::PairwisePair pair;
		// 		while(pair_factory.setNextPair(pair)) {
		// 			std::vector<double> rmses;
		// 			double sum = 0;
		// 			for(const auto tpos : iter::range(haps.tarNamesVec_.size())) {
		// 				if(haps.targetsEncodeBySamp_[group.second[pair.col_]][tpos]  + haps.targetsEncodeBySamp_[group.second[pair.row_]][tpos] == 2) {
		// 					double current_sum = 0;
		// 					for(const auto hapPos : iter::range(haps.numberOfHapsPerTarget_[tpos])) {
		// 						current_sum += std::pow(hapsEncodeBySampRelAbund[group.second[pair.col_]][haps.tarStart_[tpos] + hapPos] - hapsEncodeBySampRelAbund[group.second[pair.row_]][haps.tarStart_[tpos] + hapPos],2);
		// 						sum += std::pow(hapsEncodeBySampRelAbund[group.second[pair.col_]][haps.tarStart_[tpos] + hapPos] - hapsEncodeBySampRelAbund[group.second[pair.row_]][haps.tarStart_[tpos] + hapPos],2);
		// 					}
		// 					rmses.emplace_back(std::sqrt(current_sum));
		// 				}
		// 			}
		// 			std::cout << haps.sampNamesVec_[group.second[pair.col_]] << "\t" <<  haps.sampNamesVec_[group.second[pair.row_]] << "\t" << std::sqrt(sum) << "\t" << std::sqrt(sum)/rmses.size() << "\t" << std::sqrt(sum/rmses.size()) << "\t" << vectorMean(rmses) << std::endl;
		// 		}
		// 	}
		// }

		watch.startNewLap("output");
		//mat.writeGraph(outFile);
		outFile << "sample\tgroup";
		outFile << std::endl;
		std::map<uint32_t, std::vector<uint32_t>> groupIndexes;
		std::vector<uint32_t> allGroupedIndices;
		for(const auto & n : iter::enumerate(mat.graph_->nodes_)) {
			groupIndexes[n.element->group_].emplace_back(n.index);
			if (std::numeric_limits<uint32_t>::max() != n.element->group_) {
				allGroupedIndices.emplace_back(n.index);
			}
		}
		if (!doNotBreakWithRmse) {
			//fill in all RMSEs now for the grouped data
			watch.startNewLap("calculate all RMSEs");
			uint32_t pairBatchCount = 100;
			PairwisePairFactory allGroupedSamplesFactory(allGroupedIndices.size());
			std::function<void()> calcRMSEs =
					[&allGroupedSamplesFactory,
					  &pairwiseRMSEs,
					  &pairwise_cccs,
					  &pairBatchCount,
						&haps,&allGroupedIndices]() {
						PairwisePairFactory::PairwisePairVec pairs;
						while(allGroupedSamplesFactory.setNextPairs(pairs, pairBatchCount)) {
							for(const auto & groupedPair : pairs.pairs_) {
								auto colSamplePos = allGroupedIndices[groupedPair.col_];
								auto rowSamplePos = allGroupedIndices[groupedPair.row_];
								std::vector<double> rmses;
							  std::vector<double> row_values;
							  std::vector<double> col_values;
								double sum = 0;
								for(const auto tpos : iter::range(haps.tarNamesVec_.size())) {
									if(haps.targetsEncodeBySamp_[colSamplePos][tpos]  + haps.targetsEncodeBySamp_[rowSamplePos][tpos] == 2) {
										double current_sum = 0;
										for(const auto hapPos : iter::range(haps.numberOfHapsPerTarget_[tpos])) {
											current_sum += std::pow(haps.hapsEncodeBySampRelAbund_[colSamplePos][haps.tarStart_[tpos] + hapPos] - haps.hapsEncodeBySampRelAbund_[rowSamplePos][haps.tarStart_[tpos] + hapPos],2);
											sum +=         std::pow(haps.hapsEncodeBySampRelAbund_[colSamplePos][haps.tarStart_[tpos] + hapPos] - haps.hapsEncodeBySampRelAbund_[rowSamplePos][haps.tarStart_[tpos] + hapPos],2);
										  row_values.emplace_back(haps.hapsEncodeBySampRelAbund_[rowSamplePos][haps.tarStart_[tpos] + hapPos]);
										  col_values.emplace_back(haps.hapsEncodeBySampRelAbund_[colSamplePos][haps.tarStart_[tpos] + hapPos]);
										}
										rmses.emplace_back(std::sqrt(current_sum));
									}
								}
								//only add if mean RMSE is less than the cut off
								//if(vectorMean(rmses) < rmseCutOffToBreak) {
							  auto ccc_calc = ConcordanceCalculator::lins_ccc_with_ci(row_values, col_values);
								pairwiseRMSEs[colSamplePos][rowSamplePos] = std::sqrt(sum/rmses.size());
								pairwiseRMSEs[rowSamplePos][colSamplePos] = std::sqrt(sum/rmses.size());
							  pairwise_cccs[colSamplePos][rowSamplePos] = ccc_calc.ccc;
							  pairwise_cccs[rowSamplePos][colSamplePos] = ccc_calc.ccc;
							}
						}
			};

			njh::concurrent::runVoidFunctionThreaded(calcRMSEs, pars.numThreads);

		}
		std::unordered_map<uint32_t, std::map<std::string, double>> groups_jaccard_stats;
		std::unordered_map<uint32_t, std::map<std::string, double>> groups_rmse_stats;
	  std::unordered_map<uint32_t, std::map<std::string, double>> groups_ccc_stats;
		for(const auto & group : groupIndexes) {

			if (!doNotBreakWithRmse && std::numeric_limits<uint32_t>::max() != group.first) {
				std::vector<double> rmsesWithinGroup;
			  std::vector<double> cccsWithinGroup;
				PairwisePairFactory pfac(group.second.size());
				PairwisePairFactory::PairwisePair pair;
				// std::cout << "group: " << group.first << std::endl;
				while (pfac.setNextPair(pair)) {
					// std::cout << haps.sampNamesVec_[group.second[pair.col_]] << " vs " << haps.sampNamesVec_[group.second[pair.row_]] << " rmse: " << pairwiseRMSEs[group.second[pair.col_]][group.second[pair.row_]] << std::endl;
					rmsesWithinGroup.emplace_back(pairwiseRMSEs[group.second[pair.col_]][group.second[pair.row_]]);
				  cccsWithinGroup.emplace_back(pairwise_cccs[group.second[pair.col_]][group.second[pair.row_]]);
				}
				groups_rmse_stats[group.first] = getStatsOnVec(rmsesWithinGroup);
			  groups_ccc_stats[group.first]  = getStatsOnVec(cccsWithinGroup);
				// std::cout << njh::conToStr(rmsesWithinGroup, ",") << std::endl;
				// std::cout << "stats: " << njh::json::toJson(stats) << std::endl;
				// std::cout << std::endl;
			}


			if (std::numeric_limits<uint32_t>::max() != group.first) {
				std::vector<double> jaccardWithinGroup;
				PairwisePairFactory pfac(group.second.size());
				PairwisePairFactory::PairwisePair pair;
				// std::cout << "group: " << group.first << std::endl;
				while (pfac.setNextPair(pair)) {
					// std::cout << haps.sampNamesVec_[group.second[pair.col_]] << " vs " << haps.sampNamesVec_[group.second[pair.row_]] << " rmse: " << pairwiseRMSEs[group.second[pair.col_]][group.second[pair.row_]] << std::endl;
					jaccardWithinGroup.emplace_back(1 - mat.points_[group.second[pair.col_]]->vals_[group.second[pair.row_]]);
				}
				groups_jaccard_stats[group.first] = getStatsOnVec(jaccardWithinGroup);
			}
			for (const auto &idx: group.second) {
				if (lociCoveragePerSample[haps.sampNamesVec_[idx]] < minimumLociCoverageToKeepSamples) {
					outFile << haps.sampNamesVec_[idx] << "\t" << "low_coverage_not_clustered";
				} else if (group.first == std::numeric_limits<uint32_t>::max()) {
					outFile << haps.sampNamesVec_[idx] << "\t" << "nogroup";
				} else {
					outFile << haps.sampNamesVec_[idx] << "\t" << group.first;
				}
				outFile << std::endl;
			}
		}
		if (writeOutGroupedRMSEs) {
			OutputStream outGroupCountsFile(OutOptions(njh::files::make_path(setUp.pars_.directoryName_, "rmses_of_grouped_samples.tsv.gz")));
			outGroupCountsFile << "sample";
			for (const auto & sample : allGroupedIndices) {
				outGroupCountsFile << "\t" << haps.sampNamesVec_[sample];
			}
			outGroupCountsFile << std::endl;
			for (const auto & sampleCol : allGroupedIndices) {
				outGroupCountsFile << haps.sampNamesVec_[sampleCol];
				for (const auto & sampleRow : allGroupedIndices) {
					outGroupCountsFile << "\t" << pairwiseRMSEs[sampleRow][sampleCol];
				}
				outGroupCountsFile << std::endl;
			}
		}

	  if (writeOutGroupedRMSEs) {
	    OutputStream outGroupCountsFile(OutOptions(njh::files::make_path(setUp.pars_.directoryName_, "ccc_of_grouped_samples.tsv.gz")));
	    outGroupCountsFile << "sample";
	    for (const auto & sample : allGroupedIndices) {
	      outGroupCountsFile << "\t" << haps.sampNamesVec_[sample];
	    }
	    outGroupCountsFile << std::endl;
	    for (const auto & sampleCol : allGroupedIndices) {
	      outGroupCountsFile << haps.sampNamesVec_[sampleCol];
	      for (const auto & sampleRow : allGroupedIndices) {
	        outGroupCountsFile << "\t" << pairwise_cccs[sampleRow][sampleCol];
	      }
	      outGroupCountsFile << std::endl;
	    }
	  }


		OutputStream outGroupCountsFile(OutOptions(njh::files::make_path(setUp.pars_.directoryName_, "clusters_by_jaccardTargetsShared_groupCounts.tsv")));
		outGroupCountsFile << "group\tsampleCount";
		outGroupCountsFile << "\tmin_jaccard\tmedian_jaccard\tmean_jaccard\tmax_jaccard";
		if (!doNotBreakWithRmse) {
			outGroupCountsFile << "\tmin_rmse\tmedian_rmse\tmean_rmse\tmax_rmse";
		}
	  outGroupCountsFile << "\tmin_ccc\tmedian_ccc\tmean_ccc\tmax_ccc";
		outGroupCountsFile << std::endl;
		for(const auto & group : groupIndexes) {
			if(group.first == std::numeric_limits<uint32_t>::max()) {
				uint32_t low_coverage_not_clustered_cnt = 0;
				for(const auto & samp_idx : group.second) {
					if (lociCoveragePerSample[haps.sampNamesVec_[samp_idx]] < minimumLociCoverageToKeepSamples) {
						low_coverage_not_clustered_cnt++;
					}
				}
				outGroupCountsFile << "nogroup" << "\t" << group.second.size() - low_coverage_not_clustered_cnt;
				outGroupCountsFile << "\t" << "NA"
						<< "\t" << "NA"
						<< "\t" << "NA"
						<< "\t" << "NA";
        if (!doNotBreakWithRmse) {
          outGroupCountsFile << "\t" << "NA"
              << "\t" << "NA"
              << "\t" << "NA"
              << "\t" << "NA";
        }
        outGroupCountsFile << "\t" << "NA"
            << "\t" << "NA"
            << "\t" << "NA"
            << "\t" << "NA";
				outGroupCountsFile << std::endl;

				outGroupCountsFile << "low_coverage_not_clustered" << "\t" << low_coverage_not_clustered_cnt;
				outGroupCountsFile << "\t" << "NA"
						<< "\t" << "NA"
						<< "\t" << "NA"
						<< "\t" << "NA";
				if (!doNotBreakWithRmse) {
					outGroupCountsFile << "\t" << "NA"
							<< "\t" << "NA"
							<< "\t" << "NA"
							<< "\t" << "NA";
				}
        outGroupCountsFile << "\t" << "NA"
            << "\t" << "NA"
            << "\t" << "NA"
            << "\t" << "NA";
				outGroupCountsFile << std::endl;
			} else {
				outGroupCountsFile << group.first << "\t" << group.second.size();
				outGroupCountsFile << "\t" << groups_jaccard_stats[group.first]["min"]
						<< "\t" << groups_jaccard_stats[group.first]["median"]
						<< "\t" << groups_jaccard_stats[group.first]["mean"]
						<< "\t" << groups_jaccard_stats[group.first]["max"];
				if (!doNotBreakWithRmse) {
					outGroupCountsFile << "\t" << groups_rmse_stats[group.first]["min"]
							<< "\t" << groups_rmse_stats[group.first]["median"]
							<< "\t" << groups_rmse_stats[group.first]["mean"]
							<< "\t" << groups_rmse_stats[group.first]["max"];
				}
        outGroupCountsFile << "\t" << groups_ccc_stats[group.first]["min"]
            << "\t" << groups_ccc_stats[group.first]["median"]
            << "\t" << groups_ccc_stats[group.first]["mean"]
            << "\t" << groups_ccc_stats[group.first]["max"];
				outGroupCountsFile << std::endl;
			}
		}
		if(setUp.pars_.verbose_){
			watch.logLapTimes(std::cout, true, 6, true);
		}
	}

	if(writeOutTarsAbsoluteShared){
		OutOptions hapsAbsoluteSharedBetweenSampsBetweenTargetsOutopts(njh::files::make_path(setUp.pars_.directoryName_, "hapsAbsoluteSharedBetweenSampsBetweenTargetsOut.tab.txt.gz"));
		haps.writeAbsoluteHapSharedPerSamplePerTar(hapsAbsoluteSharedBetweenSampsBetweenTargetsOutopts, setUp.pars_.verbose_);
	}

	{
		setUp.timer_.startNewLap("get population pairwise measures");
		std::vector<uint32_t> tarKeys(haps.tarNamesVec_.size());
		njh::iota<uint32_t>(tarKeys, 0);
		njh::concurrent::LockableQueue<uint32_t> tarQueue(tarKeys);
		OutputStream diversityMeasuresOut(njh::files::make_path(setUp.pars_.directoryName_, "diversityMeasuresPerTarget.tab.txt"));
		diversityMeasuresOut << "loci\tsampCount\ttotalHaps\tuniqueHaps\tSimpsonI\the\tExpP3\tExpP4\tExpP5\tsinglets\tdoublets\teffectiveNumOfAlleles\tShannonEntropyE" << '\n';
		std::mutex divOutMut;
		std::function<void()> getTargetInfo = [&tarQueue,&haps,&diversityMeasuresOut,&divOutMut](){
			uint32_t tarKey = std::numeric_limits<uint32_t>::max();
			while(tarQueue.getVal(tarKey)){


				std::vector<PopGenCalculator::PopHapInfo> hapsForTarget;
				for(const auto tarpos : iter::range(haps.numberOfHapsPerTarget_[tarKey])){
					hapsForTarget.emplace_back(PopGenCalculator::PopHapInfo(tarpos, 0));
				}
				uint32_t sampleCount = 0;
				for(const auto sampPos : iter::range(haps.hapsEncodeBySamp_.size())){
					if(haps.targetsEncodeBySamp_[sampPos][tarKey] > 0){
						++sampleCount;
					}
					for(const auto tarpos : iter::range(haps.numberOfHapsPerTarget_[tarKey])){
						if(haps.hapsEncodeBySamp_[sampPos][haps.tarStart_[tarKey] + tarpos] > 0){
							hapsForTarget[tarpos].unweighted_count_ += 1;
							hapsForTarget[tarpos].weighted_count_ += haps.hapsEncodeBySampRelAbund_[sampPos][haps.tarStart_[tarKey] + tarpos];
						}
					}
				}
				auto diversityForTar = PopGenCalculator::getGeneralMeasuresOfDiversity(hapsForTarget);
				auto totalHaps = PopGenCalculator::PopHapInfo::getTotalPopCount(hapsForTarget);
				{
					std::lock_guard<std::mutex> lock(divOutMut);
					diversityMeasuresOut << haps.tarNamesVec_[tarKey]
															<< "\t" << sampleCount
															<< "\t" << totalHaps
															<< "\t" << diversityForTar.alleleNumber_
															<< "\t" << diversityForTar.simpsonIndex_
															<< "\t" << diversityForTar.heterozygostiy_
															<< "\t" << (std::numeric_limits<long double>::max() == diversityForTar.expected_k_heterozygosities.at(3).k_heterozygosity_ ? "NA": estd::to_string(diversityForTar.expected_k_heterozygosities.at(3).k_heterozygosity_))
															<< "\t" << (std::numeric_limits<long double>::max() == diversityForTar.expected_k_heterozygosities.at(4).k_heterozygosity_ ? "NA": estd::to_string(diversityForTar.expected_k_heterozygosities.at(4).k_heterozygosity_))
															<< "\t" << (std::numeric_limits<long double>::max() == diversityForTar.expected_k_heterozygosities.at(5).k_heterozygosity_ ? "NA": estd::to_string(diversityForTar.expected_k_heterozygosities.at(5).k_heterozygosity_))
															<< "\t" << diversityForTar.singlets_
															<< "\t" << diversityForTar.doublets_
															<< "\t" << diversityForTar.effectiveNumOfAlleles_
															<< "\t" << diversityForTar.ShannonEntropyE_
															<< '\n';
				}
			}


		};
		njh::concurrent::runVoidFunctionThreaded(getTargetInfo, pars.numThreads);
	}



	if(!metaFieldsToCalcPopDiffs.empty()){
		setUp.timer_.startNewLap("get population pairwise measures");


		auto popMeasuresDir = njh::files::makeDir(setUp.pars_.directoryName_, njh::files::MkdirPar{"popDiffMeasures"});



		std::vector<uint32_t> tarKeys(haps.tarNamesVec_.size());
		njh::iota<uint32_t>(tarKeys, 0);
		for(const auto & field : metaFieldsToCalcPopDiffs){
			njh::concurrent::LockableQueue<uint32_t> tarQueue(tarKeys);
			OutputStream diversityMeasuresOut(njh::files::make_path(popMeasuresDir, njh::pasteAsStr(field, "_diversityMeasures.tab.txt.gz")));
			diversityMeasuresOut << field << "\tloci\tsampCount\ttotalHaps\tuniqueHaps\tSimpsonI\the\tExpP3\tExpP4\tExpP5\tsinglets\tdoublets\teffectiveNumOfAlleles\tShannonEntropyE\tIn" << '\n';
			OutputStream diffMeasuresOut(njh::files::make_path(popMeasuresDir, njh::pasteAsStr(field, "_diffMeasures.tab.txt.gz")));
			OutputStream pairwiseDiffMeasuresOut(njh::files::make_path(popMeasuresDir, njh::pasteAsStr(field, "_pairwiseDiffMeasures.tab.txt.gz")));
			diffMeasuresOut << "meta" << "\t"<< "loci"
					<<"\t"<<"totalHaps"
					<<"\t"<<"uniqueHaps"
					<<"\t"<<"nsamples"
					<<"\t"<<"HsSample"
					<<"\t"<<"HsEst"
					<<"\t"<<"HtSample"
					<<"\t"<<"HtEst"
					<<"\t"<<"Gst"
					<<"\t"<<"GstEst"
					<<"\t"<<"JostD"
					<<"\t"<<"JostDEst"
					<<"\t"<<"ChaoA"
					<<"\t"<<"ChaoB"
					<<"\t"<<"JostDChaoEst"
					<<"\t"<<"In"<< std::endl;

			pairwiseDiffMeasuresOut << "loci"
					<< "\t" << field << "1"
					<< "\t" << "popMeta" << "1_totalHaps"
					<< "\t" << "popMeta" << "1_uniqueHaps"
					<< "\t" << "popMeta" << "1_samples"
					<< "\t" << "hapsOnlyIn_popMeta" << "1"
					<< "\t" << "hapsOnlyIn_popMeta" << "1CumFreq"
					<< "\t" << field << "2"
					<< "\t" << "popMeta" << "2_totalHaps"
					<< "\t" << "popMeta" << "2_uniqueHaps"
					<< "\t" << "popMeta" << "2_samples"
					<< "\t" << "hapsOnlyIn_popMeta" << "2"
					<< "\t" << "hapsOnlyIn_popMeta" << "2CumFreq"
					<< "\t" << "uniqHapsCombinedPops"
					<< "\t" << "uniqHapsSharedInPops"
					<< "\t" << "HsSample"
									<< "\t" << "HsEst"
									<< "\t" << "HtSample"
									<< "\t" << "HtEst"
									<< "\t" << "Gst"
									<< "\t" << "GstEst"
									<< "\t" << "JostD"
									<< "\t" << "JostDEst"
									<< "\t" << "ChaoA"
									<< "\t" << "ChaoB"
									<< "\t" << "JostDChaoEst"
									<< "\t" << "In"

									<< "\t" << "brayCurtisDissim"
									<< "\t" << "brayCurtisRelativeDissim"
									<< "\t" << "jaccardIndexDissim"
									<< "\t" << "sorensenDistance"
									<< "\t" << "RMSE"
									<< "\t" << "correlationDissim"
									<< "\t" << "matchingCoefficientDistance"
									<< "\t" << "plainAvalance"
									<< std::endl;

			std::mutex divOutMut;
			std::vector<std::string> sampleToMeta;
			std::unordered_set<std::string> subFields;
			for(const auto sampPos : iter::range(haps.sampNamesVec_.size())){

				sampleToMeta.emplace_back(haps.meta_->groupData_[field]->getGroupForSample(haps.sampNamesVec_[sampPos]));
				subFields.emplace(sampleToMeta.back());
			}
			std::function<void()> getPopDiffMeasures = [&tarQueue,&haps, &diversityMeasuresOut,&diffMeasuresOut,&pairwiseDiffMeasuresOut,&divOutMut,&sampleToMeta,&subFields,&field](){

				uint32_t tarKey = std::numeric_limits<uint32_t>::max();
				while(tarQueue.getVal(tarKey)){

					std::unordered_map<std::string, std::vector<PopGenCalculator::PopHapInfo>> hapsForTargetPerPopulationRaw;
					for(const auto & subField : subFields){
						for(const auto tarpos : iter::range(haps.numberOfHapsPerTarget_[tarKey])){
							hapsForTargetPerPopulationRaw[subField].emplace_back(PopGenCalculator::PopHapInfo(tarpos, 0));
						}
					}
					std::unordered_map<std::string, uint32_t> sampleCount;
					for(const auto sampPos : iter::range(haps.hapsEncodeBySamp_.size())){
						if(haps.targetsEncodeBySamp_[sampPos][tarKey] > 0){
							++sampleCount[sampleToMeta[sampPos]];
						}
						for(const auto tarpos : iter::range(haps.numberOfHapsPerTarget_[tarKey])){
							if(haps.hapsEncodeBySamp_[sampPos][haps.tarStart_[tarKey] + tarpos] > 0){
								hapsForTargetPerPopulationRaw[sampleToMeta[sampPos]][tarpos].unweighted_count_ +=1;
								hapsForTargetPerPopulationRaw[sampleToMeta[sampPos]][tarpos].weighted_count_ += haps.hapsEncodeBySampRelAbund_[sampPos][haps.tarStart_[tarKey] + tarpos];

							}
						}
					}
					std::unordered_map<std::string, std::vector<PopGenCalculator::PopHapInfo>> hapsForTargetPerPopulation;
					for(const auto & pop : hapsForTargetPerPopulationRaw){
						for(const auto & hap : pop.second){
							if(hap.unweighted_count_ > 0){
								hapsForTargetPerPopulation[pop.first].emplace_back(hap);
							}
						}
					}

					PopGenCalculator::PopDifferentiationMeasures generalDiff;
					if(hapsForTargetPerPopulation.size() > 1){
						generalDiff = PopGenCalculator::getOverallPopDiffWeighted(hapsForTargetPerPopulation);
					}
					std::unordered_map<std::string, std::unordered_map<std::string, PopGenCalculator::PopDifferentiationMeasuresPairWise>> pairwiseDiffs;

					if(hapsForTargetPerPopulation.size() > 1){
						pairwiseDiffs = PopGenCalculator::getPairwisePopDiffWeighted(hapsForTargetPerPopulation);
					}
					std::unordered_map<std::string, PopGenCalculator::DiversityMeasures> divMeausresPerPop;
					for(const auto & hapsForPop : hapsForTargetPerPopulation){
						divMeausresPerPop[hapsForPop.first] = PopGenCalculator::getGeneralMeasuresOfDiversity(hapsForPop.second);
					}
					{
						std::lock_guard<std::mutex> lock(divOutMut);
						uint32_t grandTotalHaps = 0;
						uint32_t grandTotalSamples = 0;
						std::unordered_map<std::string, uint32_t> totalHapsPerPop;
						for(const auto & popDiv : divMeausresPerPop){
							auto totalHaps = PopGenCalculator::PopHapInfo::getTotalPopCount(hapsForTargetPerPopulation[popDiv.first]);
							totalHapsPerPop[popDiv.first] = totalHaps;
							grandTotalHaps += totalHaps;
							grandTotalSamples += sampleCount[popDiv.first];
							diversityMeasuresOut
							<< popDiv.first
							<< "\t" << haps.tarNamesVec_[tarKey]
																	<< "\t" << sampleCount[popDiv.first]
																	<< "\t" << totalHaps
																	<< "\t" << popDiv.second.alleleNumber_
																	<< "\t" << popDiv.second.simpsonIndex_
																	<< "\t" << popDiv.second.heterozygostiy_
							<< "\t" << (std::numeric_limits<long double>::max() == popDiv.second.expected_k_heterozygosities.at(3).k_heterozygosity_ ? "NA": estd::to_string(popDiv.second.expected_k_heterozygosities.at(3).k_heterozygosity_))
							<< "\t" << (std::numeric_limits<long double>::max() == popDiv.second.expected_k_heterozygosities.at(4).k_heterozygosity_ ? "NA": estd::to_string(popDiv.second.expected_k_heterozygosities.at(4).k_heterozygosity_))
							<< "\t" << (std::numeric_limits<long double>::max() == popDiv.second.expected_k_heterozygosities.at(5).k_heterozygosity_ ? "NA": estd::to_string(popDiv.second.expected_k_heterozygosities.at(5).k_heterozygosity_))
																	<< "\t" << popDiv.second.singlets_
																	<< "\t" << popDiv.second.doublets_
																	<< "\t" << popDiv.second.effectiveNumOfAlleles_
																	<< "\t" << popDiv.second.ShannonEntropyE_
																	<< "\t" << (hapsForTargetPerPopulation.size() > 1 ? generalDiff.informativenessForAssignPerPopulation_[popDiv.first] : 0)
																	<< '\n';
						}

						if(hapsForTargetPerPopulation.size() > 1){
							diffMeasuresOut << field << "\t"
									<< haps.tarNamesVec_[tarKey]
									<<"\t"<< grandTotalHaps
									<<"\t"<< haps.numberOfHapsPerTarget_[tarKey]
									<<"\t"<< grandTotalSamples
									<<"\t"<< generalDiff.hsSample_
									<<"\t"<< generalDiff.hsEst_
									<<"\t"<< generalDiff.htSample_
									<<"\t"<< generalDiff.htEst_
									<<"\t"<< generalDiff.gst_
									<<"\t"<< generalDiff.gstEst_
									<<"\t"<< generalDiff.jostD_
									<<"\t"<< generalDiff.jostDEst_
									<<"\t"<< generalDiff.chaoA_
									<<"\t"<< generalDiff.chaoB_
									<<"\t"<< generalDiff.jostDChaoEst_
									<<"\t"<< generalDiff.informativenessForAssign_<< std::endl;
							auto keys = getVectorOfMapKeys(pairwiseDiffs);
							njh::sort(keys);
							for(const auto & key : keys){
								auto subKeys = getVectorOfMapKeys(pairwiseDiffs.at(key));
								njh::sort(subKeys);
								for(const auto & subKey : subKeys){
									pairwiseDiffMeasuresOut << haps.tarNamesVec_[tarKey]
											<< "\t" << key
											<< "\t" << totalHapsPerPop[key]
											<< "\t" << divMeausresPerPop[key].alleleNumber_
											<< "\t" << sampleCount[key]
											<< "\t" << pairwiseDiffs.at(key).at(subKey).uniqueHapsInPop1_
											<< "\t" << pairwiseDiffs.at(key).at(subKey).uniqueHapsInPop1CumFreq_
											<< "\t" << subKey
											<< "\t" << totalHapsPerPop[subKey]
											<< "\t" << divMeausresPerPop[subKey].alleleNumber_
											<< "\t" << sampleCount[subKey]
											<< "\t" << pairwiseDiffs.at(key).at(subKey).uniqueHapsInPop2_
											<< "\t" << pairwiseDiffs.at(key).at(subKey).uniqueHapsInPop2CumFreq_

											<< "\t" << pairwiseDiffs.at(key).at(subKey).uniqueHapsAll_
											<< "\t" << pairwiseDiffs.at(key).at(subKey).uniqueHapsShared_
																<<"\t"<< pairwiseDiffs.at(key).at(subKey).genDiffMeasures_.hsSample_
																<<"\t"<< pairwiseDiffs.at(key).at(subKey).genDiffMeasures_.hsEst_
																<<"\t"<< pairwiseDiffs.at(key).at(subKey).genDiffMeasures_.htSample_
																<<"\t"<< pairwiseDiffs.at(key).at(subKey).genDiffMeasures_.htEst_
																<<"\t"<< pairwiseDiffs.at(key).at(subKey).genDiffMeasures_.gst_
																<<"\t"<< pairwiseDiffs.at(key).at(subKey).genDiffMeasures_.gstEst_
																<<"\t"<< pairwiseDiffs.at(key).at(subKey).genDiffMeasures_.jostD_
																<<"\t"<< pairwiseDiffs.at(key).at(subKey).genDiffMeasures_.jostDEst_
																<<"\t"<< pairwiseDiffs.at(key).at(subKey).genDiffMeasures_.chaoA_
																<<"\t"<< pairwiseDiffs.at(key).at(subKey).genDiffMeasures_.chaoB_
																<<"\t"<< pairwiseDiffs.at(key).at(subKey).genDiffMeasures_.jostDChaoEst_
																<<"\t"<< pairwiseDiffs.at(key).at(subKey).genDiffMeasures_.informativenessForAssign_


																<< "\t" << pairwiseDiffs.at(key).at(subKey).brayCurtisDissim_
																<< "\t" << pairwiseDiffs.at(key).at(subKey).brayCurtisRelativeDissim_
																<< "\t" << pairwiseDiffs.at(key).at(subKey).jaccardIndexDissim_
																<< "\t" << pairwiseDiffs.at(key).at(subKey).sorensenDistance_
																<< "\t" << pairwiseDiffs.at(key).at(subKey).RMSE_
																<< "\t" << pairwiseDiffs.at(key).at(subKey).halfR_
																<< "\t" << pairwiseDiffs.at(key).at(subKey).matchingCoefficientDistance_
																<< "\t" << pairwiseDiffs.at(key).at(subKey).plainAvalance_

																<< std::endl;
								}
							}
						}
					}
				}
			};


			njh::concurrent::runVoidFunctionThreaded(getPopDiffMeasures, pars.numThreads);


		}
	}




	setUp.timer_.logLapTimes(setUp.rLog_.runLogFile_, true, 6, true);
	return 0;
}


} //namespace njhseq
