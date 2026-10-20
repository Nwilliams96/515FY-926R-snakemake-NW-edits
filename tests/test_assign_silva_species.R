suppressPackageStartupMessages(library(readr))

script_argument <- commandArgs(trailingOnly = FALSE)
test_path <- sub("^--file=", "", script_argument[grepl("^--file=", script_argument)])
script_path <- normalizePath(file.path(
  dirname(test_path), "..", "workflow", "scripts", "assign_silva_species.R"
))
source(script_path)

test_skip_older_release <- function() {
  root <- tempfile("species-release-skip-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  taxonomy_path <- file.path(root, "input.tsv")
  output_path <- file.path(root, "output.tsv")
  summary_path <- file.path(root, "summary.tsv")
  writeLines(c(
    "#OTUID\ttaxonomy\tconfidence",
    "001\td__Bacteria; p__Test; g__Test\t0.95",
    "002\tEukaryota;Test:plas\t0.8"
  ), taxonomy_path)

  # Missing sequence/reference paths deliberately prove this branch cannot
  # consume SILVA 144 sequences when an older classifier was selected.
  assign_silva_species(
    taxonomy_path, file.path(root, "missing.fasta"), NULL,
    output_path, summary_path, enabled = FALSE, silva_version = "138.2"
  )
  stopifnot(identical(
    readBin(taxonomy_path, "raw", file.info(taxonomy_path)$size),
    readBin(output_path, "raw", file.info(output_path)$size)
  ))
  summary <- read_tsv(summary_path, col_types = cols(.default = col_character()))
  stopifnot(identical(summary$status, "skipped"))
  stopifnot(identical(summary$silva_version, "138.2"))
  stopifnot(grepl("only for SILVA 144", summary$reason))
  stopifnot(!"species_labels_added" %in% names(summary))
}

test_skip_older_release()
cat("Older SILVA species matching skip test passed\n")
