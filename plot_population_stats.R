# Load required packages
library(dplyr)
library(ggplot2)
library(patchwork)  # to combine plots
library(tidyr)
library(readr)
library(tidyverse)

# Function to read and reshape pixy output files
pixy_to_long <- function(pixy_files) {
  pixy_df <- list()

  for(i in seq_along(pixy_files)) {
    stat_file_type <- gsub(".*_|.txt", "", pixy_files[i])
    
    df <- read_delim(pixy_files[i], delim = "\t")

    if (stat_file_type == "pi") {
      df_long <- df %>%
        pivot_longer(cols = -c(pop, window_pos_1, window_pos_2, chromosome),
                     names_to = "statistic", values_to = "value") %>%
        rename(pop1 = pop) %>%
        mutate(pop2 = NA_character_)
    } else if (stat_file_type == "fst") {
      df_long <- df %>%
        pivot_longer(cols = -c(pop1, pop2, window_pos_1, window_pos_2, chromosome),
                     names_to = "statistic", values_to = "value")
    } else {
      stop(paste("Unexpected stat file type:", stat_file_type))
    }

    pixy_df[[i]] <- df_long
  }

  bind_rows(pixy_df) %>%
    arrange(pop1, pop2, chromosome, window_pos_1, statistic)
}

# pixy code: pixy --vcf vcfFile.vcf.gz --populations populations.txt --output_folder pixy_output --output_prefix sexFst  --stats fst pi  --bypass_invariant_check  --window_size 1000
# File setup
# folder containing files for pixy results for nucleotide diversity and Fixation index (FST) 
pixy_folder <- "pixy_output"
# pi.txt --> nucleotide diversity
# fst.txt --> Fst (fixation index) 
pixy_files <- list.files(pixy_folder, full.names = TRUE)
pixy_df <- pixy_to_long(pixy_files)


# Labeler for stats
pixy_labeller <- as_labeller(c(
  avg_pi = "pi",          # Greek symbol for pi
  avg_wc_fst = "F[ST]"
), default = label_parsed)


library(tidyverse)
library(patchwork)  # for stacking plots easily

# Assuming pixy_df is your long data frame with columns:
# pop1, pop2, chromosome, window_pos_1, window_pos_2, statistic, value

# Clean chromosome name
pixy_df <- pixy_df %>%
  mutate(
    chrom_clean = ifelse(grepl("^patChr", chromosome),
                        sub(".*r", "", chromosome),
                        chromosome),
    chrom_clean = factor(chrom_clean, levels = c(as.character(1:9), "X", "Y")),
    midpoint = (window_pos_1 + window_pos_2) / 2,
    bin = floor(midpoint / 1e6)
  )

# Male and Female pi separately (pop1 is male or female, statistic == avg_pi)
pi_separate <- pixy_df %>%
  filter(statistic == "avg_pi", pop2 %>% is.na() | pop2 == "") %>%
  group_by(pop1, chrom_clean, bin) %>%
  summarise(
    mean_value = mean(value, na.rm = TRUE),
    start_pos = bin * 1e6,
    end_pos = (bin + 1) * 1e6,
    .groups = "drop"
  ) %>%
  mutate(mid_bin_pos = (start_pos + end_pos) / 2 / 1e6)

# Total pi (average of male & female)
pi_total <- pi_separate %>%
  group_by(chrom_clean, bin) %>%
  summarise(
    mean_value = mean(mean_value, na.rm = TRUE),
    start_pos = first(start_pos),
    end_pos = first(end_pos),
    .groups = "drop"
  ) %>%
  mutate(mid_bin_pos = (start_pos + end_pos) / 2 / 1e6)

# Fst data (pop1 == female, pop2 == male, statistic == avg_wc_fst)
fst <- pixy_df %>%
  filter(statistic == "avg_wc_fst", pop1 == "female", pop2 == "male") %>%
  group_by(chrom_clean, bin) %>%
  summarise(
    mean_value = mean(value, na.rm = TRUE),
    start_pos = bin * 1e6,
    end_pos = (bin + 1) * 1e6,
    .groups = "drop"
  ) %>%
  mutate(mid_bin_pos = (start_pos + end_pos) / 2 / 1e6)

# Plot 1: Male and Female pi lines together
p1 <- ggplot(pi_separate, aes(x = mid_bin_pos, y = mean_value, color = pop1)) +
  geom_line(size = 1) +
  facet_wrap(~chrom_clean, scales = "free_x", nrow = 1) +
  labs(title = NULL, x =NULL, y = expression(pi)) +
  theme_classic() +
  scale_color_manual(values = c("female" = "red", "male" = "blue")) +
  theme(legend.position = "top")

# Plot 2: Total pi
p2 <- ggplot(pi_total, aes(x = mid_bin_pos, y = mean_value)) +
  geom_line(color = "skyblue", size = 1) +
  facet_wrap(~chrom_clean, scales = "free_x", nrow = 1) +
  labs(title = NULL, x = NULL, y = expression(pi)) +
  theme_classic()

# Plot 3: FST
p3 <- ggplot(fst, aes(x = mid_bin_pos, y = mean_value)) +
  geom_line(color = "darkgreen", size = 1) +
  facet_wrap(~chrom_clean, scales = "free_x", nrow = 1) +
  labs(title = NULL, x = "Position (Mb)", y = expression(F[ST])) +
  theme_classic()

p1 <- p1 + theme(
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank()
)

p2 <- p2 + theme(
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank()
)
p2 <- p2 + theme(strip.background = element_blank(), strip.text.x = element_blank())
p3 <- p3 + theme(strip.background = element_blank(), strip.text.x = element_blank())
final_plot <- p1 / p2 / p3 + plot_layout(heights = c(1, 1, 1))

# Display plot
print(final_plot)

#outputFile <- 'outputName'
ggsave("output_file_name.pdf", final_plot, width = 14, height = 10)
ggsave("output_file_name.svg", final_plot, width = 14, height = 10)

options(cli.num_colors = 0)
options(pillar.bold = FALSE)
options(pillar.subtle = FALSE)
options(pillar.sigfig = 4) 


pixy_df <- pixy_df %>%
  mutate(
    region = case_when(
      chrom_clean == "Y" & midpoint >= 1 & midpoint <= 94e6 ~ "SDR",
      chrom_clean == "Y" & midpoint > 94e6 ~ "YPAR",
      chrom_clean == "X" & midpoint >= 1 & midpoint <= 205e6 ~ "X-NR",
      chrom_clean == "X" & midpoint > 205e6 ~ "XPAR",
      chrom_clean %in% as.character(1:9) ~ "Autosome",
      TRUE ~ "Other"
    )
  )

  
library(dplyr)

# Pi by sex and region
pi_mean_region_sex <- pixy_df %>%
  filter(statistic == "avg_pi", pop1 %in% c("male", "female")) %>%
  group_by(pop1, region) %>%
  summarise(
    mean_pi   = mean(value, na.rm = TRUE),
    median_pi = median(value, na.rm = TRUE)
  ) %>%
  arrange(pop1, region)

print("Mean and median nucleotide diversity (pi) by region for males and females:")
print(pi_mean_region_sex)

# Total pi by region (all samples)
pi_mean_region_total <- pixy_df %>%
  filter(statistic == "avg_pi") %>%
  group_by(region) %>%
  summarise(
    mean_pi   = mean(value, na.rm = TRUE),
    median_pi = median(value, na.rm = TRUE)
  ) %>%
  arrange(region)

print("Mean and median nucleotide diversity (pi) by region (all samples):")
print(pi_mean_region_total)


# Filter FST data for regions of interest
fst_regions <- pixy_df %>%
  filter(statistic == "avg_wc_fst", region %in% c("SDR", "Autosome", "X-NR", "XPAR", "YPAR"))

# Calculate mean and median FST by region
fst_stats_region <- fst_regions %>%
  group_by(region) %>%
  summarise(
    mean_fst   = mean(value, na.rm = TRUE),
    median_fst = median(value, na.rm = TRUE)
  ) %>%
  arrange(region)

# Print results
print("Mean and median FST by region:")
print(fst_stats_region)



fst_subset_autosome <- fst_regions %>% filter(region %in% c("SDR", "Autosome"))
fst_subset_xnr      <- fst_regions %>% filter(region %in% c("SDR", "X-NR"))
fst_subset_xpar     <- fst_regions %>% filter(region %in% c("SDR", "XPAR"))
fst_subset_ypar     <- fst_regions %>% filter(region %in% c("SDR", "YPAR"))

fst_sdr_vs_autosome <- wilcox.test(value ~ region, data = fst_subset_autosome)
fst_sdr_vs_xnr      <- wilcox.test(value ~ region, data = fst_subset_xnr)
fst_sdr_vs_xpar     <- wilcox.test(value ~ region, data = fst_subset_xpar)
fst_sdr_vs_ypar     <- wilcox.test(value ~ region, data = fst_subset_ypar)


print("Mean total nucleotide diversity (pi) by region (all samples):")
#print(as.data.frame(pi_mean_region_total))

print(paste("SDR vs Autosome p-value:", fst_sdr_vs_autosome$p.value))
print(paste("SDR vs X-NR p-value:    ", fst_sdr_vs_xnr$p.value))
print(paste("SDR vs XPAR p-value:    ", fst_sdr_vs_xpar$p.value))
print(paste("SDR vs YPAR p-value:    ", fst_sdr_vs_ypar$p.value))
