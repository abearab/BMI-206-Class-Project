library(tidyverse)
library(dplyr)
library(tidyr)

Tcell <- read.table("500kb_Tcell_allcvar.txt", sep=' ', header=TRUE)
fibroblast <- read.table("500kb_fibroblast_allcvar_part.txt", sep=' ', header=TRUE)
Tcell_nocov <- read.table("500kb_Tcell_nocvar_parts.txt", sep=' ', header=TRUE)


Tcell_track <- Tcell %>%
  filter(p < 0.05 & boot_basic_p < 0.05) %>%
  # Separate the peak column into chr, start, and end
  separate(peak, into = c("chr", "start", "end"), sep = "[:-]") %>%
  # Convert start and end columns to numeric (optional)
  mutate(start = as.numeric(start), end = as.numeric(end))

Tcell_nocov_track <- Tcell_nocov %>%
  filter(p < 0.05 & boot_basic_p < 0.05) %>%
  # Separate the peak column into chr, start, and end
  separate(peak, into = c("chr", "start", "end"), sep = "[:-]") %>%
  # Convert start and end columns to numeric (optional)
  mutate(start = as.numeric(start), end = as.numeric(end))

fibroblast_track <- fibroblast %>%
  filter(p < 0.05 & boot_basic_p < 0.05) %>%
  # Separate the peak column into chr, start, and end
  separate(peak, into = c("chr", "start", "end"), sep = "[:-]") %>%
  # Convert start and end columns to numeric (optional)
  mutate(start = as.numeric(start), end = as.numeric(end))

bed_file <- Tcell_nocov_track %>%
  dplyr::select(chr, start, end, gene)  # Select relevant columns for the BED file

# Write the data frame to a BED file
write.table(
  bed_file,
  file = "Tcell_nocov.bed",         # Output file name
  quote = FALSE,              # No quotes around values
  sep = "\t",                 # Tab-delimited format
  row.names = FALSE,          # No row names
  col.names = FALSE           # No column headers (BED files don't use headers)
) 

