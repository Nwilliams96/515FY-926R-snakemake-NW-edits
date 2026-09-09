#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Biostrings)
  library(dada2)
  library(dplyr)
  library(readr)
  library(stringr)
})

append_species_labels <- function(taxonomy, species_matches) {
  taxonomy %>%
    left_join(species_matches, by = "ASV_hash") %>%
    mutate(
      Classified_Genus = str_match(
        Taxonomy,
        "(?:^|;\\s*)g__([^;]*)"
      )[, 2] %>% trimws() %>% na_if(""),
      Is_SILVA = str_detect(Taxonomy, "(?:^|;\\s*)d__"),
      Genus_Agrees =
        Is_SILVA &
        !is.na(Exact_Match_Genus) &
        !is.na(Exact_Match_Species) &
        !is.na(Classified_Genus) &
        str_to_lower(Classified_Genus) == str_to_lower(Exact_Match_Genus),
      Species = if_else(
        Genus_Agrees,
        paste(Exact_Match_Genus, Exact_Match_Species),
        NA_character_
      ),
      Taxonomy = if_else(
        !is.na(Species) & !str_detect(Taxonomy, "(?:^|;\\s*)s__"),
        paste0(str_remove(Taxonomy, "\\s*;\\s*$"), "; s__", Species),
        Taxonomy
      )
    )
}

assign_silva_species <- function(
  taxonomy_path,
  sequence_path,
  reference_path,
  output_path,
  summary_path
) {
  taxonomy <- read_tsv(
    taxonomy_path,
    show_col_types = FALSE,
    name_repair = "minimal"
  )
  if (ncol(taxonomy) < 3 || !"taxonomy" %in% names(taxonomy)) {
    stop("The exported taxonomy table must contain ASV ID, taxonomy, and confidence columns.")
  }
  names(taxonomy)[1] <- "ASV_hash"
  taxonomy <- taxonomy %>% rename(Taxonomy = taxonomy)

  sequences <- readDNAStringSet(sequence_path)
  query_ids <- names(sequences)
  if (is.null(query_ids) || any(query_ids == "") || anyDuplicated(query_ids)) {
    stop("The representative-sequence FASTA must contain unique ASV IDs.")
  }

  classified_genera <- str_match(
    taxonomy$Taxonomy,
    "(?:^|;\\s*)g__([^;]*)"
  )[, 2] %>% trimws() %>% na_if("")
  eligible_ids <- taxonomy$ASV_hash[
    str_detect(taxonomy$Taxonomy, "(?:^|;\\s*)d__") &
      !is.na(classified_genera)
  ]
  eligible_sequences <- sequences[names(sequences) %in% eligible_ids]

  if (length(eligible_sequences) == 0) {
    species_matches <- tibble(
      ASV_hash = character(),
      Exact_Match_Genus = character(),
      Exact_Match_Species = character()
    )
  } else {
    assignments <- dada2::assignSpecies(
      as.character(eligible_sequences),
      reference_path,
      allowMultiple = FALSE,
      tryRC = TRUE,
      verbose = TRUE
    )
    species_matches <- tibble(
      ASV_hash = names(eligible_sequences),
      Exact_Match_Genus = assignments[, "Genus"],
      Exact_Match_Species = assignments[, "Species"]
    )
  }

  augmented <- append_species_labels(taxonomy, species_matches)
  summary <- augmented %>%
    summarise(
      total_16S_ASVs = n(),
      eligible_SILVA_ASVs = sum(Is_SILVA & !is.na(Classified_Genus)),
      exact_unambiguous_species_matches = sum(!is.na(Exact_Match_Species)),
      species_labels_added = sum(Genus_Agrees, na.rm = TRUE),
      exact_matches_withheld_by_genus_check = sum(
        !is.na(Exact_Match_Species) & !Genus_Agrees,
        na.rm = TRUE
      ),
      unmatched_or_ambiguous =
        eligible_SILVA_ASVs - exact_unambiguous_species_matches
    )

  output <- augmented %>%
    select(ASV_hash, taxonomy = Taxonomy, confidence)
  names(output)[1] <- "#OTUID"

  write_tsv(output, output_path, na = "")
  write_tsv(summary, summary_path, na = "")
}

if (exists("snakemake")) {
  sink(snakemake@log[[1]], append = TRUE, split = TRUE)
  on.exit(sink(), add = TRUE)
  assign_silva_species(
    snakemake@input[["taxonomy"]],
    snakemake@input[["sequences"]],
    snakemake@input[["species_reference"]],
    snakemake@output[["taxonomy"]],
    snakemake@output[["summary"]]
  )
}
