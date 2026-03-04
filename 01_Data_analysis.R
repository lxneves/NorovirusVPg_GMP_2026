# Proteomics analysis of VPg nucleotidyation in Murine norovirus #
# Leandro Xavier Neves, PhD (orcid.org/0000-0002-6074-1025 - github.com/lxneves) #

# Load packages
library(fs)
library(readr)
library(purrr)
library(dplyr)
library(stringr)
library(ggplot2)
library(tidyr)

# Define directory
path_fig <- "Y:/Leandro/2025/Nucleotidylation/"
path_exports <- "Y:/Leandro/2025/Nucleotidylation/Synthetic_peptides/AnnotatedSpectra/exports"


# HCD, ETD & EThcD optimisation----

# Find all .txt files recursively. These will contain intensity values for each fragment type (b, y, c, z, neutral losses, etc)
txt_files <- dir_ls(path_exports, recurse = TRUE, glob = "*.txt")

# Read .txt files into a list
msms_list <- txt_files %>%
  set_names(path_ext_remove(path_file(.))) %>%  # name list elements by filename (without .txt)
  map(read_tsv)

# Rename dfs with folder's name
msms_list <- imap(msms_list, ~ .x %>% 
                    mutate(source = .y))

# Remove rows of unidentified fragments
msms_list <- lapply(msms_list, function(df) {
  df[!is.na(df$`ion type`), ]
})

# Merge as a single df
desired_cols <- c("intensity", "ion type", "ion loss", "source", "signal to noise")

merged_msms <- msms_list %>%
  lapply(function(df) df[, desired_cols, drop = FALSE]) %>%
  bind_rows()

# Load up the Spectra_meta.tsv. This contains spectrum scores and summary of dissociation conditions
msms_tsv <- read.table("Y:/Leandro/2025/Nucleotidylation/Synthetic_peptides/AnnotatedSpectra/exports/Spectra_meta.tsv",
                       quote = NULL, check.names = FALSE, sep = "\t", header = TRUE)

# Create column source to match merged_msms
msms_tsv$source <- paste0(msms_tsv$`raw file`,"_",msms_tsv$scan)

# Copy information fragmentation details to merged_msms. Select columns of interest
desired_cols <- c("spectrum fragmentation", "fragmentation energy (ETD)",
                  "fragmentation energy (HCD)", "precursor m/z", "precursor charge",
                  "PSM score", "peptide sequence")

# Join into merged_msms by matching source == raw
merged_msms <- merged_msms %>%
  left_join(
    msms_tsv %>%
      select(any_of(c("source", desired_cols))) %>%
      distinct(`source`, .keep_all = TRUE),   
    by = c("source" = "source")
  )

# Summarise intensity
summary_msms <- merged_msms %>%
  group_by(
    `spectrum fragmentation`,
    `fragmentation energy (ETD)`,
    `fragmentation energy (HCD)`,
    `ion type`,
    `precursor charge`,
  ) %>%
  summarise(
    sum_intensity    = sum(intensity, na.rm = TRUE),
    mean_intensity   = mean(intensity, na.rm = TRUE),
    median_intensity = median(intensity, na.rm = TRUE),
    sd_intensity     = sd(intensity, na.rm = TRUE),
    n                = sum(!is.na(intensity)),
    mean_score       = mean(`PSM score`, na.rm = TRUE),
    sd_score         = sd(`PSM score`, na.rm = TRUE),
    .groups = "drop"
  )

# Summarise score
summary_score <- merged_msms %>%
  group_by(
    `spectrum fragmentation`,
    `fragmentation energy (ETD)`,
    `fragmentation energy (HCD)`,
    `precursor charge`,
  ) %>%
  summarise(
    n                = sum(!is.na(intensity)),
    mean_score       = mean(`PSM score`, na.rm = TRUE),
    sd_score         = sd(`PSM score`, na.rm = TRUE),
    .groups = "drop"
  )

## Plots----

# Define colours for each ion type
ion_palette <- c(
  "b" = "darkblue",
  "y" = "red",
  "c" = "black",
  "z" = "darkgrey",
  "Precursor" = "darkgreen",
  "Diagnostic" = "yellow3"
)

# Define ion types per fragmentation method
ETD <- c("c", "z", "Precursor", "Diagnostic")
HCD <- c("b", "y", "Precursor", "Diagnostic")
EThcD <- c("c", "z", "b", "y", "Precursor", "Diagnostic")

summary_msms$`precursor charge` <- gsub("2","Charge 2+",summary_msms$`precursor charge`)
summary_msms$`precursor charge` <- gsub("3","Charge 3+",summary_msms$`precursor charge`)

summary_score$`precursor charge` <- gsub("2","Charge 2+",summary_score$`precursor charge`)
summary_score$`precursor charge` <- gsub("3","Charge 3+",summary_score$`precursor charge`)

### ETD----
# subset ETD data
temp <- summary_msms %>%
  filter(`spectrum fragmentation` == "ETD") %>%
  mutate(
    `fragmentation energy (ETD)` = suppressWarnings(as.numeric(`fragmentation energy (ETD)`))
  ) %>%
  filter(!is.na(`fragmentation energy (ETD)`), mean_intensity > 0)
         
# Keep only ions from ETD series
temp <- temp %>%
  dplyr::filter(`ion type` %in% ETD)

# Add standard deviation
temp <- temp %>%
  dplyr::mutate(
    ymin = log10(mean_intensity - sd_intensity),
    ymax = log10(mean_intensity + sd_intensity),
    y    = log10(mean_intensity)
  )

# ggplot
etd <- ggplot(
  temp,
  aes(x = `fragmentation energy (ETD)`,
      y = log10(mean_intensity),
      color = `ion type`,
      group = `ion type`)
) +
  geom_line() +
  geom_point(size = 1.8) +
  facet_wrap(~ `precursor charge`, scales = NULL, ncol = 1) +
  scale_color_manual(
    values = ion_palette,
    limits = ETD   # only use colors for the selected ions, in this order
  ) +
  labs(
    title = "MS/MS ion series versus ETD reaction time",
    x = "ETD reaction time (ms)",
    y = "Log10-average intensity",
    color = "Ion series"
  ) +
  theme_classic() +
  theme(panel.grid.minor = element_blank()) + scale_x_continuous(breaks = c(50, 100, 150, 200, 250, 300))

# SD bars
etd +
  geom_errorbar(
    aes(ymin = ymin, ymax = ymax),
    width = 0.15,        # adjust the bar width
    linewidth = 0.5,
    alpha = 0.3,
    show.legend = FALSE
  )

ggsave(plot = last_plot(), filename = paste0(path_fig,"ETD_sd.svg"),
       dpi = 300, width = 110, height = 100, units = "mm", bg = "white")

### HCD----
# subset data
temp <- summary_msms %>%
  filter(`spectrum fragmentation` == "HCD") %>%
  mutate(
    `fragmentation energy (HCD)` = suppressWarnings(as.numeric(`fragmentation energy (HCD)`))
  ) %>%
  filter(!is.na(`fragmentation energy (HCD)`), mean_intensity > 0)

# Keep only ions from HCD series
temp <- temp %>%
  dplyr::filter(`ion type` %in% HCD)

temp <- temp %>%
  dplyr::filter(`fragmentation energy (HCD)` > 0)

# Add standard deviation
temp <- temp %>%
  dplyr::mutate(
    ymin = log10(mean_intensity - sd_intensity),
    ymax = log10(mean_intensity + sd_intensity),
    y    = log10(mean_intensity)
  )

# ggplot
hcd <- ggplot(
  temp,
  aes(x = `fragmentation energy (HCD)`,
      y = log10(mean_intensity),
      color = `ion type`,
      group = `ion type`)
) +
  geom_line() +
  geom_point(size = 1.8) +
  facet_wrap(~ `precursor charge`, scales = NULL, ncol = 1) +
  scale_color_manual(
    values = ion_palette,
    limits = HCD   # only use colors for the selected ions, in this order
  ) +
  labs(
    title = "MS/MS ion series per NCE",
    x = "Normalised Collision Energy (NCE)",
    y = "Log10-average intensity",
    color = "Ion series"
  ) +
  theme_classic() +
  theme(panel.grid.minor = element_blank(),)
        
hcd

# SD bars
hcd +
  geom_errorbar(
    aes(ymin = ymin, ymax = ymax),
    width = 0.15,        # adjust the bar width
    linewidth = 0.5,
    alpha = 0.3,
    show.legend = FALSE
  )

ggsave(plot = last_plot(), filename = paste0(path_fig,"HCD_sd.svg"),
       dpi = 300, width = 110, height = 100, units = "mm", bg = "white")


### EThcD----
# subset data
temp <- summary_msms %>%
  filter(`spectrum fragmentation` == "EThcD") %>%
  mutate(
    `fragmentation energy (ETD)` = suppressWarnings(as.numeric(`fragmentation energy (ETD)`)),
    `fragmentation energy (HCD)` = suppressWarnings(as.numeric(`fragmentation energy (HCD)`))
  ) %>%
  filter(!is.na(`fragmentation energy (ETD)`), mean_intensity > 0)

# Keep only ions from EThcD series
temp <- temp %>%
  dplyr::filter(`ion type` %in% EThcD)

# Add standard deviation
temp <- temp %>%
  dplyr::mutate(
    ymin = log10(mean_intensity - sd_intensity),
    ymax = log10(mean_intensity + sd_intensity),
    y    = log10(mean_intensity)
  )

temp$EThcD <- paste0(temp$`fragmentation energy (HCD)`,"_",temp$`fragmentation energy (ETD)`)

# Define X axis group order
temp <- temp %>%
  mutate(EThcD = factor(EThcD, levels = c("15_50","25_50","35_50","15_100","25_100",
                                          "35_100","15_150","25_150","35_150","15_200",
                                          "25_200","35_200","15_250","25_250","35_250",
                                          "15_300","25_300","35_300")))
# ggplot
ethcd <- ggplot(
  temp,
  aes(x = `EThcD`,
      y = log10(mean_intensity),
      color = `ion type`,
      group = `ion type`)
) +
  geom_line() +
  geom_point(size = 1.8) +
  facet_wrap(~ `precursor charge`, scales = NULL, ncol = 1) +
  scale_color_manual(
    values = ion_palette,
    limits = EThcD   # only use colors for the selected ions, in this order
  ) +
  labs(
    title = "MS/MS ion series per ETD reaction time and HCD NCE",
    x = "ETD reaction time (ms) & % NCE",
    y = "Log10-average intensity",
    color = "Ion series"
  ) +
  theme_classic() +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

ethcd

# SD bars
ethcd +
  geom_errorbar(
    aes(ymin = ymin, ymax = ymax),
    width = 0.15,        # adjust the bar width
    linewidth = 0.5,
    alpha = 0.3,
    show.legend = FALSE
  )

ggsave(plot = last_plot(), filename = paste0(path_fig,"EThcD_sd2.svg"),
       dpi = 300, width = 200, height = 100, units = "mm", bg = "white")

## Write .csv----
write.table(summary_score, "./summary_score.csv",
            sep = ",", quote = FALSE, row.names = FALSE)

write.table(msms_tsv, "./RAW_files.csv",
            sep = ",", quote = FALSE, row.names = FALSE)

write.table(merged_msms, "./summary_msms.csv",
            sep = ",", quote = FALSE, row.names = FALSE)


# PRM in  MNV-infected BV-2 cells----

# Load Skyline results
prm <- read.csv("Y:/Leandro/2025/Nucleotidylation/Skyline/LN_Precursor_Quant_MS1_MS2.csv",
                check.names = FALSE)

# Remove rows of test injection with synthetic peptide standards
test_inj <- unique(prm$`File Name`)[1] # .raw file name of the test injection 

prm <- prm[prm$`File Name` != test_inj,]

# Define conditions and bioreplicates
prm$Bioreplicate <- str_sub(prm$`Replicate Name`,-4,-4)

prm$Timepoint <- str_extract(prm$`Replicate Name`, "\\d+(?=h_)")

prm <- prm %>%
  mutate(Treatment = case_when(
    str_detect(`Replicate Name`, regex("mock", ignore_case = TRUE)) ~ "Mock",
    TRUE ~ "MNV"
  ))

prm$Group <- paste0(prm$Treatment," ",prm$Timepoint," hpi")

## Box plot----

# Prepare numeric columns
prm$`Total Area MS1` <- as.numeric(prm$`Total Area MS1`)
prm$`Total Area Fragment` <- as.numeric(prm$`Total Area Fragment`)

prm$`Total Area MS1`[prm$`Total Area MS1` == 0 ] <-NA
prm$`Total Area Fragment`[prm$`Total Area Fragment` == 0 ] <-NA

# Calculate LOD as 3*SD of mock samples
lod_stats <- prm %>%
  filter(Treatment == "Mock") %>%
  group_by(`Peptide Modified Sequence`) %>%
  summarise(
    mean_fragment = mean(`Total Area Fragment`, na.rm = TRUE),
    sd_fragment   = sd(`Total Area Fragment`, na.rm = TRUE),
    lod           = mean_fragment + 3 * sd_fragment
  )

# Add LOD to main df
prm <- prm %>%
  left_join(lod_stats %>% select(`Peptide Modified Sequence`, lod),
            by = "Peptide Modified Sequence")

# Define X axis group order
prm <- prm %>%
  mutate(Group = factor(Group, levels = c("Mock 0 hpi",
                                          "MNV 0 hpi",
                                          "MNV 4 hpi",
                                          "MNV 8 hpi",
                                          "MNV 12 hpi",
                                          "Mock 12 hpi")))

# plot
ggplot(prm, aes(x = Group, y = `Total Area Fragment`,
                color = `Peptide Modified Sequence`)) +
                #fill = `Peptide Modified Sequence`)) +
  geom_hline(data = lod_stats, aes(yintercept = lod, colour = `Peptide Modified Sequence`),
    linetype = "dashed", linewidth = 0.3)  +
  #geom_text(data = lod_stats, aes(x = Inf, y = lod,
  #                                label = "LOD",
  #                                colour = `Peptide Modified Sequence`),
  #          hjust = 1.1, vjust = -0.5, size = 2) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA, linewidth = 0.1,
               position = position_dodge(width = 0.75)) +
  geom_jitter(alpha = 0.4, size = 0.7, position = position_jitterdodge(
    jitter.width = 0.15, dodge.width = 0.75))+
  scale_y_log10() +
  scale_colour_manual(values = c("darkred","black")) +
    theme_minimal(base_size = 12) +
  labs(y = "Log10 Peptide Intensity", colour = "Peptide") +
  theme(axis.text.x = element_text(angle = 45, size = 6, hjust = 1, vjust = 1, face = "bold"),
        axis.line = element_line(color = "darkgrey"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 5),
        axis.title.y = element_text(size = 8),
        legend.position = "top",
        legend.text = element_text(size = 6),
        legend.title = element_blank(),
        panel.grid.major.x = element_blank()) 
  
ggsave(plot = last_plot(), filename = paste0(path_fig,"PRM_boxplot.svg"),
       dpi = 300, width = 70, height = 70, units = "mm", bg = "white")

## Occupancy----

# Filter data points below LOD
prm_filtered <- prm %>%
  filter(`Total Area Fragment` >= lod)

# Select data for occupancy calculation
occupancy_df <- prm_filtered %>%
  select(`Replicate Name`,
         `Peptide Modified Sequence`,
         `Total Area Fragment`,
         `Group`) %>%
  pivot_wider(
    names_from = `Peptide Modified Sequence`,
    values_from = `Total Area Fragment`)

# Calculate occupancy
occupancy_df <- occupancy_df %>%
  mutate(
    occupancy = `GLTDEEY[+345]DEFKK` /
      (`GLTDEEY[+345]DEFKK` + `GLTDEEYDEFKK`))

# Descriptive statistics
occupancy_summary <- occupancy_df %>%
  group_by(Group) %>%
  summarise(
    mean_occupancy = round(mean(occupancy, na.rm = TRUE),2),
    sd_occupancy   = round(sd(occupancy, na.rm = TRUE),2),
    n              = sum(!is.na(occupancy)))
