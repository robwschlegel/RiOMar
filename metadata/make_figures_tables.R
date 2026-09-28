# metadata/make_figures_tables.R
#
# A CHECKLIST, not a pipeline: verifies every figure/table manuscript.tex
# expects has actually been generated, and reports where it lives and
# whether it's stale (older than its own data source). It does NOT
# generate, copy, or assemble any figure itself.
#
# All real figure/table generation lives in the numbered code/ pipeline
# (func/figure.R + func/figure.py, called from code/5_figures.py per
# CLAUDE.md: `python code/5_figures.py`; func/driver_interactions.R, called
# from code/4_time_series.py). Run those first. This script just tells you
# what's missing so nothing silently falls back to \figplaceholder's grey box
# in the compiled PDF.
#
# Driven by two CSV registries, both the single source of truth for their
# own concern:
#   - metadata/figure_table_registry.csv: "what manuscript slot is this,
#     what number does it currently have, which R function/script renders
#     it, where does its output live." Renumbering a figure/table (moving it
#     in manuscript.tex, or reassigning its number) is a one-row edit to
#     that CSV; this script and every figure.R/figure.py generator look up
#     their current_number/output_subdir from the same row instead of
#     hardcoding it.
#   - metadata/paragraph_source_registry.csv: "which script/data file
#     produced the numbers cited in this manuscript paragraph." One row per
#     Results/Discussion/Appendix paragraph containing a quantitative claim,
#     anchored by a verbatim text snippet (most of these paragraphs have no
#     \label of their own to key off).
#
# Usage:
#   Rscript metadata/make_figures_tables.R   (from the repo root)
# or, interactively, source() this file and call make_all_figures_tables().


# Setup -----------------------------------------------------------------

# This script lives in metadata/ (tracked) while manuscript.tex and
# references.bib live in the gitignored manuscript/ -- so locate the repo root
# from this script's own folder, then point at manuscript/ explicitly.
script_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)))
if (length(script_dir) == 1 && script_dir != "") {
  proj_dir <- normalizePath(file.path(script_dir, ".."))
} else {
  proj_dir <- getwd()  # fallback when source()-d interactively (RStudio project = repo root)
}
manuscript_dir <- file.path(proj_dir, "manuscript")

func_dir <- file.path(proj_dir, "func")
# figure.R loads tidyverse (for read_csv() etc.) then sources util.R itself
# (figure_table_registry, get_registry_row(), registry_filename()) -- sourcing
# util.R separately first would break on read_csv() before tidyverse loads.
source(file.path(func_dir, "figure.R"))

manuscript_tex_path <- file.path(manuscript_dir, "manuscript.tex")
references_bib_path <- file.path(manuscript_dir, "references.bib")

# Where the real pipeline (code/5_figures.py / code/4_time_series.py ->
# func/figure.R, func/driver_interactions.R) writes its output, per the
# figures/ARTICLE/FIGURE_X/ convention used for every manuscript figure.
riomar_figure_root <- file.path(proj_dir, "figures")


# Helpers -----------------------------------------------------------------

# Checks that a figure/table source file the real pipeline should have
# produced actually exists, reporting its path and whether it's older than a
# reference source file (a cheap staleness signal -- doesn't recompute
# anything, just compares mtimes). Never copies, generates, or assembles
# anything.
check_figure_exists <- function(label, path, newer_than = NULL) {
  if (!file.exists(path)) {
    message("[MISSING] ", label, ": ", path,
           "\n           -> run the pipeline stage that generates it (see metadata/figure_table_registry.csv).")
    return(invisible(FALSE))
  }
  stale_note <- ""
  if (!is.null(newer_than) && file.exists(newer_than) &&
      file.info(path)$mtime < file.info(newer_than)$mtime) {
    stale_note <- paste0(" [STALE -- older than ", newer_than, "]")
  }
  message("[ok] ", label, ": ", path, stale_note)
  invisible(TRUE)
}

# Cross-checks a figure row's registry-expected output path against the
# literal \includegraphics{...}/\figplaceholder{...} path manuscript.tex
# actually uses for that row's tex_label, catching the exact class of drift
# that let the validation figure/table sit at FIGURE_2/"Table 4" in code
# while already rendering as Figure S1/a main-sequence table in the PDF.
# Purely a text grep -- does not parse LaTeX, so a label spread across
# multiple lines or an unusual macro won't be caught.
check_manuscript_tex_path <- function(tex_label, expected_relpath) {
  if (!file.exists(manuscript_tex_path)) return(invisible(NA))
  tex_lines <- readLines(manuscript_tex_path, warn = FALSE)
  label_line <- grep(paste0("\\\\label\\{", tex_label, "\\}"), tex_lines, fixed = FALSE)
  if (length(label_line) == 0) {
    message("  [tex] no \\label{", tex_label, "} found in manuscript.tex")
    return(invisible(FALSE))
  }
  # \label{} normally follows its \includegraphics/\figplaceholder within the
  # same \begin{figure}...\end{figure} block -- search a small window around it.
  window <- tex_lines[max(1, label_line[1] - 6):label_line[1]]
  if (!any(grepl(expected_relpath, window, fixed = TRUE))) {
    message("  [tex MISMATCH] \\label{", tex_label, "}: manuscript.tex does not reference '",
           expected_relpath, "' near this label -- path/number may be stale.")
    return(invisible(FALSE))
  }
  invisible(TRUE)
}


# Per-row check ---------------------------------------------------------

check_registry_row <- function(row) {
  slot_label <- paste0(row$current_number, if (row$kind == "figure") " (fig)" else " (tab)",
                       " ", row$slot_key)

  if (row$kind == "figure") {
    if (is.na(row$output_subdir) || row$output_subdir == "") {
      message("[ok] ", slot_label, ": no numbered output folder for this slot (see notes).")
      return(invisible(TRUE))
    }
    filename <- if (!is.na(row$filename_override) && row$filename_override != "") {
      row$filename_override
    } else {
      registry_filename(row$output_subdir)
    }
    expected_path <- file.path(riomar_figure_root, "ARTICLE", row$output_subdir, filename)
    ok <- check_figure_exists(slot_label, expected_path)
    expected_relpath <- file.path("..", "figures", "ARTICLE", row$output_subdir, filename)
    check_manuscript_tex_path(row$tex_label, expected_relpath)
    return(ok)
  }

  # kind == "table"
  if (is.na(row$check_path) || row$check_path == "" || row$check_path == "hardcoded") {
    message("[ok] ", slot_label, ": hand-transcribed in manuscript.tex, no single source file to check.")
    return(invisible(TRUE))
  }
  # check_path can itself be a ";"-separated list (e.g. missing_days_table's
  # two source CSVs) -- same convention as source_files below, and the same
  # bug class if skipped: passing the raw joined string to file.exists()
  # checks for one file literally named "a.csv;b.csv", which never exists,
  # so a table with multiple real source files always reported [MISSING].
  paths <- trimws(strsplit(row$check_path, ";")[[1]])
  ok <- purrr::map_lgl(paths, ~ check_figure_exists(paste0(slot_label, " source data"), .x))
  invisible(all(ok))
}


# Paragraph numeric-source check -----------------------------------------
#
# Same checklist spirit as check_registry_row() above, one level deeper:
# metadata/paragraph_source_registry.csv maps each manuscript paragraph
# that states a quantitative claim (Results/Discussion/Appendix only -- see
# that CSV's own notes) to the script(s)/data file(s) that produced its
# numbers, anchored by a verbatim text snippet rather than a \label (most of
# these paragraphs have no \label of their own).

# Confirms row$paragraph_snippet still appears verbatim somewhere in
# manuscript.tex -- if not, the paragraph was edited, reworded, or removed
# since the row was written, and the row needs re-anchoring against the
# current text (or removal, if the claim itself is gone).
check_paragraph_snippet <- function(snippet) {
  if (!file.exists(manuscript_tex_path)) return(invisible(NA))
  tex_text <- paste(readLines(manuscript_tex_path, warn = FALSE), collapse = "\n")
  grepl(snippet, tex_text, fixed = TRUE)
}

# One semicolon-separated entry from source_files. Three special prefixes
# get a distinct check/message; anything else is a plain repo-relative path
# checked the same way check_figure_exists() checks a figure/table source.
check_source_entry <- function(entry) {
  entry <- trimws(entry)
  if (startsWith(entry, "PLOT-ONLY:")) {
    path <- sub("^PLOT-ONLY:", "", entry)
    if (file.exists(path)) {
      message("  [ok] (plot-only, no separate numeric source) ", path)
    } else {
      message("  [MISSING] (plot-only source) ", path)
    }
  } else if (startsWith(entry, "EXTERNAL-CITATION:")) {
    key <- sub("^EXTERNAL-CITATION:", "", entry)
    found <- file.exists(references_bib_path) &&
      any(grepl(paste0("@.*\\{", key, ","), readLines(references_bib_path, warn = FALSE)))
    if (found) {
      message("  [ok] (external citation) ", key)
    } else {
      message("  [MISSING] citekey '", key, "' not found in references.bib")
    }
  } else if (startsWith(entry, "UNVERIFIED:")) {
    message("  [WARN] no source script found: ", sub("^UNVERIFIED:", "", entry))
  } else {
    check_figure_exists(paste0("  source"), entry)
  }
}

check_paragraph_source_row <- function(row) {
  where <- paste(c(row$section, row$subsection,
                   if (!is.na(row$subsubsection) && row$subsubsection != "") row$subsubsection),
                 collapse = " > ")
  message("[", where, "] \"", row$paragraph_snippet, "...\"")

  if (isTRUE(check_paragraph_snippet(row$paragraph_snippet))) {
    message("  [ok] snippet found in manuscript.tex")
  } else {
    message("  [MISSING SNIPPET] not found in manuscript.tex -- paragraph may have been edited, reworded, moved, or removed; re-anchor this row")
  }

  purrr::walk(strsplit(row$source_files, ";")[[1]], check_source_entry)
  invisible(TRUE)
}


# Orchestrator ----------------------------------------------------------

make_all_figures_tables <- function() {

  message("== Figures ==")
  figures <- dplyr::filter(figure_table_registry, kind == "figure") |>
    dplyr::arrange(section == "supplement", current_number)
  invisible(purrr::pwalk(figures, function(...) check_registry_row(tibble::tibble(...))))

  message("\n== Tables ==")
  tables <- dplyr::filter(figure_table_registry, kind == "table") |>
    dplyr::arrange(section == "supplement", current_number)
  invisible(purrr::pwalk(tables, function(...) check_registry_row(tibble::tibble(...))))

  message("\n== Paragraph numeric sources ==")
  invisible(purrr::pwalk(paragraph_source_registry, function(...) check_paragraph_source_row(tibble::tibble(...))))

  message("\nDone. Re-run `pdflatex manuscript.tex` (x2) + `bibtex manuscript` to refresh the compiled PDF.")
}

# Only auto-run when invoked via `Rscript make_figures_tables.R`; sourcing
# this file interactively (e.g. in RStudio) will not trigger a full run.
if (identical(environment(), globalenv()) &&
    length(grep("--file=", commandArgs(trailingOnly = FALSE))) > 0) {
  make_all_figures_tables()
}
