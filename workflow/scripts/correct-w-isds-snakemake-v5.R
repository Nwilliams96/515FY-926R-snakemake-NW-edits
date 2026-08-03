#This script was written to convert ASV copies into ASV absolute copies by correcting data with internal standards. 
#It can only be used if genomic internal standards were added at the time of DNA extraction.
#Written by Nathan Williams 18/02/2026.

suppressPackageStartupMessages({
library(tidyverse)
})

#Set paths
isd_path  <- snakemake@input[["isd"]]
isd_added_path <- snakemake@input[["isd_added"]]
cruise_tag <- snakemake@config[["studyName"]]
isd_1_id <- snakemake@params[["intstd1name"]]
isd_2_id <- snakemake@params[["intstd2name"]]
isd_3_id <- snakemake@params[["intstd3name"]]
outdir <- "results/05-internal-std-corrected"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
figure_dir <- "results/06-figures"
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

write_isd_table <- function(df, stem) {
  readr::write_tsv(df, file.path(outdir, paste0(cruise_tag, ".", stem, ".tsv")))
}

#Import Data
asv_table <- read_tsv(snakemake@input[["asv_table"]])
isd <- read_tsv(isd_path , show_col_types = FALSE)
isd_added <- read_tsv(isd_added_path, show_col_types = FALSE) %>%
  rename(
    SampleID = sample,
    isd_1_ng = all_of(paste0(isd_1_id, "_ng")),
    isd_2_ng = all_of(paste0(isd_2_id, "_ng")),
    isd_3_ng = all_of(paste0(isd_3_id, "_ng"))
  ) %>%
  mutate(SampleID = str_replace_all(SampleID, "_", "-"))
samples <- isd_added

#Import Data local
isd_1_asvs <- read_delim(snakemake@input[[2]], delim = "\n", col_names = FALSE)
isd_2_asvs <- read_delim(snakemake@input[[3]], delim = "\n", col_names = FALSE)
isd_3_asvs <- read_delim(snakemake@input[[4]], delim = "\n", col_names = FALSE)

#Make the ISD dataframe lookup vectors
genome_len  <- setNames(isd$genome_len_bp, isd$internal_std_ID)
copy_number <- setNames(isd$rRNA_copy_number, isd$internal_std_ID)

#Calculate copies of each ISD added
#Set parameters
bp_weight <- 617.9
avogadro <- 6.022 * 1e23
# 1e9 is to convert to copies added per L.

isd_1len=isd$genome_len_bp[isd$internal_std_ID == isd_1_id]
isd_2len=isd$genome_len_bp[isd$internal_std_ID == isd_2_id]
isd_3len=isd$genome_len_bp[isd$internal_std_ID == isd_3_id]

isd_1copynum=isd$rRNA_copy_number[isd$internal_std_ID == isd_1_id]
isd_2copynum=isd$rRNA_copy_number[isd$internal_std_ID == isd_2_id]
isd_3copynum=isd$rRNA_copy_number[isd$internal_std_ID == isd_3_id]

# Do calculation
isd_copies_added <- isd_added %>% 
  select(SampleID, isd_1_ng, isd_2_ng, isd_3_ng) %>%
  mutate(isd_3_copies = ((((isd_3_ng/1e9) / (bp_weight * isd_3len)) * avogadro) * isd_3copynum)) %>%
  mutate(isd_2_copies = ((((isd_2_ng/1e9) / (bp_weight * isd_2len)) * avogadro) * isd_2copynum)) %>%
  mutate(isd_1_copies = ((((isd_1_ng/1e9) / (bp_weight * isd_1len)) * avogadro) * isd_1copynum))
          
#Pull internal standard copies out of ASV table frame
isd_1_ids <- pull(isd_1_asvs) %>% as.character()
isd_2_ids <- pull(isd_2_asvs) %>% as.character()
isd_3_ids <- pull(isd_3_asvs) %>% as.character()

isd_1_by_sample <- asv_table %>%
  filter(ASV_hash %in% isd_1_ids) %>%
  group_by(SampleID) %>%
  summarize(isd_1_copies_recovered = sum(Corrected_Sequence_Counts, na.rm = TRUE), .groups = "drop")

isd_2_by_sample <- asv_table %>%
  filter(ASV_hash %in% isd_2_ids) %>%
  group_by(SampleID) %>%
  summarize(isd_2_copies_recovered = sum(Corrected_Sequence_Counts, na.rm = TRUE), .groups = "drop")

isd_3_by_sample <- asv_table %>%
  filter(ASV_hash %in% isd_3_ids) %>%
  group_by(SampleID) %>%
  summarize(isd_3_copies_recovered = sum(Corrected_Sequence_Counts, na.rm = TRUE), .groups = "drop")

#Now calculate recovery ratio
merged_isd_data <- isd_copies_added %>%
  left_join(isd_1_by_sample) %>%
  left_join(isd_2_by_sample) %>%
  left_join(isd_3_by_sample)

#Calculate per-ISD recovery ratio columns (recovered / added)
merged_isd_data <- merged_isd_data %>%
  mutate(isd_3_recovery_ratio = isd_3_copies_recovered / isd_3_copies) %>%
  mutate(isd_1_recovery_ratio = isd_1_copies_recovered / isd_1_copies) %>%
  mutate(isd_2_recovery_ratio = isd_2_copies_recovered / isd_2_copies)

#Calculate the mean and median of all three combinations, as well as each combination of two.
merged_isd_data <- merged_isd_data %>%
  rowwise() %>%
  mutate(recovery_mean = mean( c(isd_1_recovery_ratio, isd_2_recovery_ratio, isd_3_recovery_ratio), na.rm = TRUE),
  recovery_median = median(c(isd_1_recovery_ratio, isd_2_recovery_ratio, isd_3_recovery_ratio),na.rm = TRUE),
  isd_1_isd_2_mean_recovery_ratio = mean( c(isd_1_recovery_ratio, isd_2_recovery_ratio), na.rm = TRUE),
  isd_1_isd_3_mean_recovery_ratio = mean( c(isd_1_recovery_ratio, isd_3_recovery_ratio), na.rm = TRUE),
  isd_2_isd_3_mean_recovery_ratio = mean( c(isd_2_recovery_ratio, isd_3_recovery_ratio), na.rm = TRUE)
  ) %>%
  ungroup()

#Make a subset of data for the plots
plot_data <- merged_isd_data %>%
  select(SampleID, isd_3_recovery_ratio, isd_1_recovery_ratio, isd_2_recovery_ratio, recovery_mean, 
         isd_1_isd_2_mean_recovery_ratio, isd_1_isd_3_mean_recovery_ratio, isd_2_isd_3_mean_recovery_ratio) %>%
  pivot_longer(
    cols = c(isd_3_recovery_ratio, isd_1_recovery_ratio, isd_2_recovery_ratio, recovery_mean, 
             isd_1_isd_2_mean_recovery_ratio, isd_1_isd_3_mean_recovery_ratio, isd_2_isd_3_mean_recovery_ratio),
    names_to = "Method",
    values_to = "Recovery"
  ) %>%
  filter(!is.na(Recovery))

#means dataframe
means_df <- plot_data %>% filter(Method %in% c("recovery_mean")) %>%
  distinct(SampleID, Recovery)

#Generate a plot to compare ratios
plot1 <- ggplot(plot_data, aes(x = SampleID, y = Recovery, colour = Method)) +
  geom_point(position = position_jitter(width = 0.15, height = 0), size = 2.5) +
  geom_line(data = means_df, aes(x = SampleID, y = Recovery, group = 1),
            inherit.aes = FALSE, colour = "black", linewidth = 1) +
  theme_minimal() +
  scale_colour_manual(values = c("Bacteria" = "#1F77B4", "Archaea" = "#2CA02C", "Eukaryota" = "#D62728", "Unassigned" = "#7F7F7F")) +
  labs(title = "Recovery Ratios per Sample",
       subtitle = "Dots = internal standards, Black line = per-sample mean",
       y = "Recovery Ratio",
       x = "Sample") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), legend.position = "bottom")

#Choose your method
recovery_ratios <- merged_isd_data %>% select(c("SampleID", "isd_3_recovery_ratio", "isd_1_recovery_ratio", "isd_2_recovery_ratio" ,
                                                "recovery_mean", "recovery_median", "isd_1_isd_2_mean_recovery_ratio", 
                                                "isd_1_isd_3_mean_recovery_ratio", "isd_2_isd_3_mean_recovery_ratio"))

#Join in recovery ratio
asv_table <- asv_table %>%
  left_join(recovery_ratios)


#add in unit for normalisation
isd_norm_fact <- samples %>% 
  select(SampleID,internal_std_normalization_factor,units)
  
asv_table <- asv_table %>% left_join(isd_norm_fact)
  
#Calculate recovery ratio
asv_table <- asv_table %>%
  mutate(Copies_isd_1_recovery_ratio = (Corrected_Sequence_Counts / isd_1_recovery_ratio)/internal_std_normalization_factor) %>%
  mutate(Copies_isd_2_recovery_ratio = (Corrected_Sequence_Counts / isd_2_recovery_ratio)/internal_std_normalization_factor) %>%
  mutate(Copies_isd_3_recovery_ratio = (Corrected_Sequence_Counts / isd_3_recovery_ratio)/internal_std_normalization_factor) %>%
  mutate(Copies_mean_recovery_ratio = (Corrected_Sequence_Counts / recovery_mean)/internal_std_normalization_factor) %>%
  mutate(Copies_median_recovery_ratio = (Corrected_Sequence_Counts / recovery_median)/internal_std_normalization_factor) %>%
  mutate(Copies_isd_1_isd_2_mean_recovery_ratio = (Corrected_Sequence_Counts /isd_1_isd_2_mean_recovery_ratio)/internal_std_normalization_factor) %>%
  mutate(Copies_isd_1_isd_3_mean_recovery_ratio = (Corrected_Sequence_Counts /isd_1_isd_3_mean_recovery_ratio)/internal_std_normalization_factor) %>%
  mutate(Copies_isd_2_isd_3_mean_recovery_ratio = (Corrected_Sequence_Counts /isd_2_isd_3_mean_recovery_ratio)/internal_std_normalization_factor)

#Remove ISDs from the data
asv_table <- asv_table %>%
  filter(!ASV_hash %in% isd_1_ids) %>%
  filter(!ASV_hash %in% isd_2_ids) %>%
  filter(!ASV_hash %in% isd_3_ids)

Domain_Totals <- asv_table %>%
  filter(!Domain %in% c("Unassigned")) %>%
  group_by(SampleID, Domain) %>%
  summarise(
    Total_isd_3_recovery_ratio     = sum(Copies_isd_3_recovery_ratio, na.rm = TRUE),
    Total_isd_1_recovery_ratio     = sum(Copies_isd_1_recovery_ratio, na.rm = TRUE),
    Total_isd_2_recovery_ratio     = sum(Copies_isd_2_recovery_ratio, na.rm = TRUE),
    Total_mean_recovery_ratio   = sum(Copies_mean_recovery_ratio, na.rm = TRUE),
    Total_median_recovery_ratio = sum(Copies_median_recovery_ratio, na.rm = TRUE),
    Total_isd_1_isd_2_mean_recovery_ratio = sum(Copies_isd_1_isd_2_mean_recovery_ratio, na.rm = TRUE),
    Total_isd_1_isd_3_mean_recovery_ratio = sum(Copies_isd_1_isd_3_mean_recovery_ratio, na.rm = TRUE),
    Total_isd_2_isd_3_mean_recovery_ratio = sum(Copies_isd_2_isd_3_mean_recovery_ratio, na.rm = TRUE),
    .groups = "drop"
  )
  
Domain_Totals_long <- Domain_Totals %>%
  pivot_longer(
    cols = starts_with("Total_"),
    names_to = "Method",
    values_to = "Total_Copies"
  ) %>%
  mutate(Method = recode(Method,
                         "Total_isd_3_recovery_ratio" = "isd 3 recovery ratio",
                         "Total_isd_1_recovery_ratio" = "isd 1 recovery ratio",
                         "Total_isd_2_recovery_ratio" = "isd 2 recovery ratio",
                         "Total_mean_recovery_ratio" = "mean recovery ratio",
                         "Total_median_recovery_ratio" = "median recovery ratio",
                         "Total_isd_1_isd_3_mean_recovery_ratio" = "isd 1 / isd 3 mean recovery ratio",
                         "Total_isd_1_isd_2_mean_recovery_ratio" = "isd 1 / isd 2 mean recovery ratio",
                         "Total_isd_2_isd_3_mean_recovery_ratio" = "isd 2 / isd 3 mean recovery ratio"))

plot2 <- ggplot(
  Domain_Totals_long,
  aes(x = SampleID, y = Total_Copies, colour = Domain, group = Domain)
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_y_log10() +
  facet_wrap(~Method, ncol = 1) +
  theme_minimal() +
  scale_colour_manual(values = c("Bacteria" = "#1F77B4", "Archaea" = "#2CA02C", "Eukaryota" = "#D62728", "Unassigned" = "#7F7F7F")) +
  labs(
    title = "Total Copies per Unit by Domain and Correction Method",
    y = "Total Copies per unit (log10 scale)",
    x = "Sample"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

# isd_3 recovery ratio wide table
asv_table_isd_3_recovery_ratio <- asv_table %>% select(SampleID, ASV_hash, Copies_isd_3_recovery_ratio) %>%
  group_by(ASV_hash, SampleID) %>%
  summarise(Abundance = sum(Copies_isd_3_recovery_ratio, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = SampleID, values_from = Abundance, values_fill = 0)

# isd_1 recovery ratio wide table
asv_table_isd_1_recovery_ratio <- asv_table %>% select(SampleID, ASV_hash, Copies_isd_1_recovery_ratio) %>% 
  group_by(ASV_hash, SampleID) %>%
  summarise(Abundance = sum(Copies_isd_1_recovery_ratio, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = SampleID, values_from = Abundance, values_fill = 0)

# isd_2 recovery ratio wide table
asv_table_isd_2_recovery_ratio <- asv_table %>% select(SampleID, ASV_hash, Copies_isd_2_recovery_ratio) %>% 
  group_by(ASV_hash, SampleID) %>%
  summarise(Abundance = sum(Copies_isd_2_recovery_ratio, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = SampleID, values_from = Abundance, values_fill = 0)

# mean_RR wide table
asv_table_mean_recovery_ratio <- asv_table %>%
  select(SampleID, ASV_hash, Copies_mean_recovery_ratio) %>% group_by(ASV_hash, SampleID) %>%
  summarise(Abundance = sum(Copies_mean_recovery_ratio, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = SampleID, values_from = Abundance, values_fill = 0)

# median_RR wide table
asv_table_median_recovery_ratio <- asv_table %>%
  select(SampleID, ASV_hash, Copies_median_recovery_ratio) %>%
  group_by(ASV_hash, SampleID) %>%
  summarise(Abundance = sum(Copies_median_recovery_ratio, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = SampleID, values_from = Abundance, values_fill = 0)

# mean_isd_1_isd_2 recovery_ratio wide table
asv_table_isd_1_isd_2_mean_recovery_ratio <- asv_table %>%
  select(SampleID, ASV_hash, Copies_isd_1_isd_2_mean_recovery_ratio) %>% group_by(ASV_hash, SampleID) %>%
  summarise(Abundance = sum(Copies_isd_1_isd_2_mean_recovery_ratio, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = SampleID, values_from = Abundance, values_fill = 0)

# mean_isd_1_isd_3 recovery_ratio wide table
asv_table_isd_1_isd_3_mean_recovery_ratio <- asv_table %>%
  select(SampleID, ASV_hash, Copies_isd_1_isd_3_mean_recovery_ratio) %>% group_by(ASV_hash, SampleID) %>%
  summarise(Abundance = sum(Copies_isd_1_isd_3_mean_recovery_ratio, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = SampleID, values_from = Abundance, values_fill = 0)

# mean_isd_2_isd_3 recovery_ratio wide table
asv_table_isd_2_isd_3_mean_recovery_ratio <- asv_table %>%
  select(SampleID, ASV_hash, Copies_isd_2_isd_3_mean_recovery_ratio) %>% group_by(ASV_hash, SampleID) %>%
  summarise(Abundance = sum(Copies_isd_2_isd_3_mean_recovery_ratio, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = SampleID, values_from = Abundance, values_fill = 0)

#Write tsv files - one for each correction method.
write_isd_table(asv_table, "ISD_corrected_asv_table")
write_isd_table(asv_table_isd_3_recovery_ratio, "asv_table_isd_3_recovery_ratio")
write_isd_table(asv_table_isd_1_recovery_ratio, "asv_table_isd_1_recovery_ratio")
write_isd_table(asv_table_isd_2_recovery_ratio, "asv_table_isd_2_recovery_ratio")
write_isd_table(asv_table_mean_recovery_ratio, "asv_table_mean_recovery_ratio")
write_isd_table(asv_table_median_recovery_ratio, "asv_table_median_recovery_ratio")
write_isd_table(asv_table_isd_1_isd_2_mean_recovery_ratio, "asv_table_isd_1_isd_2_mean_recovery_ratio")
write_isd_table(asv_table_isd_1_isd_3_mean_recovery_ratio, "asv_table_isd_1_isd_3_mean_recovery_ratio")
write_isd_table(asv_table_isd_2_isd_3_mean_recovery_ratio, "asv_table_isd_2_isd_3_mean_recovery_ratio")

#Write plots
pdf(snakemake@output[["recovery_plot"]], width = 12, height = 6)
print(plot1)
dev.off()
pdf(snakemake@output[["domain_plot"]], width = 12, height = 8)
print(plot2)
dev.off()

# PNG copies are embedded directly in the final self-contained HTML report.
png(snakemake@output[["recovery_plot_png"]], width = 1800, height = 900, res = 150)
print(plot1)
dev.off()
png(snakemake@output[["domain_plot_png"]], width = 1800, height = 1200, res = 150)
print(plot2)
dev.off()
