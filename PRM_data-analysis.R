# Proteomics analysis of VPg nucleotidyation in Murine norovirus #
# Leandro Xavier Neves, PhD (orcid.org/0000-0002-6074-1025 - github.com/lxneves) #


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



# Define ion series to be extracted from dfs
series <- c("b", "y", "z", "c", "Precursor", "Diagnostic")




# PRM-based quant in infected BV2 cells----

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
