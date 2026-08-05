#Purpose: normalize raw counts of 16S and 18S ASVs by dividing the 16S reads and 18S reads by the % passing DADA2
#then converting to proportions and multiplying counts of 18S sequences by the bias against them and adjusting the 16S sequences to match that bias (as user specifies).
#Written by Nathan Lloyd Robert Williams 17 March 2026 (snakemake version)

#Load dependencies
library(tidyverse)
library(data.table)

############ SNAKEMAKE DATA IMPORT ############

# Get input paths from Snakemake
raw16s_files     <- snakemake@input[["raw16S"]]
raw18s_files     <- snakemake@input[["raw18S"]]
readsum_files    <- snakemake@input[["read_summary"]]
bioanalyzer_path <- snakemake@input[["bioanalyzer"]]
stats16s_files   <- snakemake@input[["stats16S"]]
stats18s_files   <- snakemake@input[["stats18S"]]

# Import 16S
raw_16S <- purrr::map_dfr(raw16s_files, ~
  readr::read_delim(.x, delim="\t", skip=1, col_names=TRUE, show_col_types=FALSE)
) %>%
  rename(
    ASV_hash = `#OTU ID`,
    Taxonomy = taxonomy
  ) %>%
  as.data.table()

# Import 18S
raw_18S <- purrr::map_dfr(raw18s_files, ~
  readr::read_delim(.x, delim="\t", skip=1, col_names=TRUE, show_col_types=FALSE)
) %>%
  rename(
    ASV_hash = `#OTU ID`,
    Taxonomy = taxonomy
  ) %>%
  as.data.table()

# Import read_summary
read_summary <- readr::read_tsv(readsum_files, show_col_types = FALSE) %>% as.data.table()

# Import bioanalyzer. Force sample_type to character because readr otherwise
# parses values such as "16S" and "18S" as the numbers 16 and 18.
bioanalyzer_results <- readr::read_tsv(
  bioanalyzer_path,
  col_types = readr::cols(
    sample_type = readr::col_character(),
    amount_pM = readr::col_double()
  ),
  show_col_types = FALSE
) %>% as.data.table()

required_bioanalyzer_columns <- c("sample_type", "amount_pM")
missing_bioanalyzer_columns <- setdiff(required_bioanalyzer_columns, names(bioanalyzer_results))
if (length(missing_bioanalyzer_columns) > 0) {
  stop(
    "config/bioanalyzer.tsv is missing required column(s): ",
    paste(missing_bioanalyzer_columns, collapse = ", "),
    ". Expected columns: sample_type and amount_pM."
  )
}

bioanalyzer_results[, sample_type := toupper(trimws(as.character(sample_type)))]
bioanalyzer_results[, amount_pM := as.numeric(amount_pM)]

# Import Stats
statistics_18S <- purrr::map_dfr(stats18s_files, ~ readr::read_delim(.x, delim = "\t", show_col_types = FALSE)) %>% slice(-1)
statistics_16S <- purrr::map_dfr(stats16s_files, ~ readr::read_delim(.x, delim = "\t", show_col_types = FALSE)) %>% slice(-1)

#0. - Calculate the correction factor

#Turn read counts into objects
PROK_reads <-  read_summary$PROK_reads[1] #replace with total PROK reads if this was part of a run with multiple projects. (i.e For AMT30 G4, and mocks were run so total 16S reads was actually 13651543)
EUK_reads  <-  read_summary$EUK_reads[1] #replace with total EUK reads if this was part of a run with multiple projects. (i.e For AMT30 G4, and mocks were run so total 18S reads was actually 370459)

#Calculate PROK ratio
PROK_frac <- PROK_reads / (PROK_reads + EUK_reads)

#Calculate EUK ratio
EUK_frac <- EUK_reads / (PROK_reads + EUK_reads)

#Turn concentrations into objects
Concentration_16S <- bioanalyzer_results[sample_type == "16S", amount_pM]
Concentration_18S <- bioanalyzer_results[sample_type == "18S", amount_pM]

if (
  length(Concentration_16S) != 1 ||
  length(Concentration_18S) != 1 ||
  anyNA(c(Concentration_16S, Concentration_18S)) ||
  any(c(Concentration_16S, Concentration_18S) < 0) ||
  (Concentration_16S + Concentration_18S) <= 0
) {
  stop(
    "config/bioanalyzer.tsv must contain exactly one 16S row and one 18S row, ",
    "with non-negative numeric amount_pM values whose sum is greater than zero."
  )
}

#Calculate 16S ratio
bioanalyzer_frac_16S <- Concentration_16S / (Concentration_16S + Concentration_18S)

#Calculate 18S ratio
bioanalyzer_frac_18S <- Concentration_18S / (Concentration_16S + Concentration_18S)

#Set the correction factors for each file
#Correction factor for 18S and 16S, as determined empirically - please see notes for calculation instructions.
#Change this to reflect your dataset! Placeholder value of 2 from the mock communities
correction_factor_16S = bioanalyzer_frac_16S / PROK_frac
correction_factor_18S = bioanalyzer_frac_18S / EUK_frac

#Import .16S.all-16S-seqs.with-tax.tsv

#1. Calculate percent of reads that passed DADA2 denoising for both proks and euks
statistics_16S <- statistics_16S %>%
  rename(SampleID = `sample-id`) %>% 
  rename(non_chimeric = `non-chimeric`) %>% 
  mutate(non_chimeric = as.numeric(non_chimeric), input = as.numeric(input)) %>%
  as.data.table() %>% 
  select(SampleID, non_chimeric, input) %>%
  mutate(percentage_passed_final = non_chimeric/input) %>% #Calculate ratio of reads that passed DADA2 denoising 16S
  distinct(SampleID, .keep_all = TRUE)

statistics_18S <- statistics_18S %>% 
  rename(SampleID = `sample-id`) %>% 
  rename(non_chimeric = `non-chimeric`) %>%
  mutate(non_chimeric = as.numeric(non_chimeric), input = as.numeric(input)) %>%
  as.data.table() %>% 
  select(SampleID, non_chimeric, input) %>%
  mutate(percentage_passed_final = non_chimeric/input) %>% #Calculate ratio of reads that passed DADA2 denoising 18S
  distinct(SampleID, .keep_all = TRUE)

#2. Normalize ASV counts (divide counts of ASVs/ percent passed for each sample, multiply EUK ASV counts by the bias you specified)
#Arrange and connect 16S data and then make the DADA2 statistics calculation
raw_16S_long <- raw_16S %>%
  pivot_longer(cols = where(is.numeric),  
               names_to = "SampleID",   
               values_to = "read_abundance") %>%  
  mutate(read_abundance = as.numeric(read_abundance)) %>%
  left_join(statistics_16S, by = "SampleID") %>%
  mutate(corrected_read_abundance = (read_abundance * correction_factor_16S) / percentage_passed_final) %>%
  mutate(corrected_read_abundance_dada2 = read_abundance / percentage_passed_final)

#Arrange and connect 18S data and then make the DADA2 statistics calculation
raw_18S_long <- raw_18S %>%
  pivot_longer(cols = where(is.numeric),
               names_to = "SampleID",
               values_to = "read_abundance") %>%
  mutate(read_abundance = as.numeric(read_abundance)) %>%
  left_join(statistics_18S, by = "SampleID") %>%
  mutate(corrected_read_abundance = (read_abundance * correction_factor_18S) / percentage_passed_final) %>%
  mutate(corrected_read_abundance_dada2 = read_abundance / percentage_passed_final)

#Join the two data frames together, pivot wider and then export
combined_asv_long_corrected <- rbind(raw_18S_long, raw_16S_long) %>%
  select(SampleID, ASV_hash, Taxonomy, corrected_read_abundance)

combined_asv_long_corrected_dada2 <- rbind(raw_18S_long, raw_16S_long) %>%
  select(SampleID, ASV_hash, Taxonomy, corrected_read_abundance_dada2)

combined_asv_long_no_correction <-  rbind(raw_18S_long, raw_16S_long) %>%
  select(SampleID, ASV_hash, Taxonomy, read_abundance)

#Convert to wide format
combined_asv_corrected <- combined_asv_long_corrected %>%
  pivot_wider(names_from = SampleID, values_from = corrected_read_abundance, values_fill = 0)

#Convert to wide format
combined_asv_corrected_dada2 <- combined_asv_long_corrected_dada2 %>%
  pivot_wider(names_from = SampleID, values_from = corrected_read_abundance_dada2, values_fill = 0)

#Convert to wide format
combined_asv_no_correction <- combined_asv_long_no_correction %>%
  pivot_wider(names_from = SampleID, values_from = read_abundance, values_fill = 0)

ProPortal_Assigned_Cyanobacteria <- read_tsv(snakemake@input[["proportalclassification"]], col_names = FALSE)
colnames(ProPortal_Assigned_Cyanobacteria) <- c("ASV_hash", "ProPortal_ASV_Ecotype")

combined_asv_corrected <- combined_asv_corrected %>% left_join(ProPortal_Assigned_Cyanobacteria) %>%
  select(1:2, ProPortal_ASV_Ecotype, everything())

combined_asv_corrected_dada2 <- combined_asv_corrected_dada2 %>% left_join(ProPortal_Assigned_Cyanobacteria) %>%
  select(1:2, ProPortal_ASV_Ecotype, everything())

combined_asv_no_correction <- combined_asv_no_correction %>% left_join(ProPortal_Assigned_Cyanobacteria) %>%
  select(1:2, ProPortal_ASV_Ecotype, everything())

#Write out file
write_tsv(combined_asv_corrected, snakemake@output[["mergedtabledada218Scorrected"]])
write_tsv(combined_asv_corrected_dada2, snakemake@output[["mergedtabledada2"]])
write_tsv(combined_asv_no_correction, snakemake@output[["mergedtableuncorrected"]])
