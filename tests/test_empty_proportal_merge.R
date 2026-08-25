suppressPackageStartupMessages(library(tidyverse))

script_path <- normalizePath(file.path(
  dirname(sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))])),
  "..", "workflow", "scripts", "correct_16S_18S_ASV-snakemake-v1.R"
))

setClass(
  "MockMergeSnakemake",
  slots = c(input = "list", output = "list")
)

root <- tempfile("empty-proportal-merge-")
dir.create(root, recursive = TRUE)
on.exit(unlink(root, recursive = TRUE), add = TRUE)

writeLines(
  c(
    "# Constructed from biom file",
    "#OTU ID\tSample-1\ttaxonomy",
    "prok-asv\t100\td__Bacteria;p__Test;c__Test;o__Test;f__Test;g__Test;s__Test"
  ),
  file.path(root, "16S.tsv")
)
writeLines(
  c(
    "# Constructed from biom file",
    "#OTU ID\tSample-1\ttaxonomy",
    "euk-asv\t50\tEukaryota;Test;Test;Test;Test;Test;Test;Test;Test"
  ),
  file.path(root, "18S.tsv")
)
write_tsv(
  tibble(PROK_reads = 100, EUK_reads = 50),
  file.path(root, "read-summary.tsv")
)
write_tsv(
  tibble(sample_type = c("16S", "18S"), amount_pM = c(1, 1)),
  file.path(root, "amplicon-molarities.tsv")
)

stats <- tibble(
  `sample-id` = c("type", "Sample-1"),
  input = c("numeric", "100"),
  `non-chimeric` = c("numeric", "80")
)
write_tsv(stats, file.path(root, "16S-stats.tsv"))
stats$input <- c("numeric", "50")
stats$`non-chimeric` <- c("numeric", "40")
write_tsv(stats, file.path(root, "18S-stats.tsv"))

# A valid ProPortal run may produce a completely empty file when no ASVs match.
file.create(file.path(root, "proportal.tsv"))

snakemake <- new(
  "MockMergeSnakemake",
  input = list(
    raw16S = file.path(root, "16S.tsv"),
    raw18S = file.path(root, "18S.tsv"),
    read_summary = file.path(root, "read-summary.tsv"),
    amplicon_concentrations = file.path(root, "amplicon-molarities.tsv"),
    stats16S = file.path(root, "16S-stats.tsv"),
    stats18S = file.path(root, "18S-stats.tsv"),
    proportalclassification = file.path(root, "proportal.tsv")
  ),
  output = list(
    mergedtabledada218Scorrected = file.path(root, "merged-corrected.tsv"),
    mergedtabledada2 = file.path(root, "merged-dada2.tsv"),
    mergedtableuncorrected = file.path(root, "merged-uncorrected.tsv")
  )
)

suppressMessages(source(script_path, local = globalenv()))

for (output_path in unname(snakemake@output)) {
  stopifnot(file.exists(output_path))
  merged <- read_tsv(output_path, show_col_types = FALSE)
  stopifnot("ProPortal_ASV_Ecotype" %in% names(merged))
  stopifnot(all(is.na(merged$ProPortal_ASV_Ecotype)))
  stopifnot(setequal(merged$ASV_hash, c("prok-asv", "euk-asv")))
}

cat("Empty ProPortal merge test passed\n")
