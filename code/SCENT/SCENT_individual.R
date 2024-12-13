library(tidyverse)
library(dplyr)
library(tidyr)
Tcell <- read.table("/Users/zoechen/Documents/Study/UCSF/Fall_2024/BMI_206_Biostats/Projects/500kb_Tcell_allcvar.txt", sep=' ', header=TRUE)
fibroblast <- read.table("/Users/zoechen/Documents/Study/UCSF/Fall_2024/BMI_206_Biostats/Projects/500kb_fibroblast_allcvar_part.txt", sep=' ', header=TRUE)
Tcell_nocov <- read.table("/Users/zoechen/Documents/Study/UCSF/Fall_2024/BMI_206_Biostats/Projects/500kb_Tcell_nocvar_parts.txt", sep=' ', header=TRUE)

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
  file = "/Users/zoechen/Documents/Study/UCSF/Fall_2024/BMI_206_Biostats/Projects/Tcell_nocov.bed",         # Output file name
  quote = FALSE,              # No quotes around values
  sep = "\t",                 # Tab-delimited format
  row.names = FALSE,          # No row names
  col.names = FALSE           # No column headers (BED files don't use headers)
) 

#-------------------------------------------------------------------------------------
#transfer to pair bedpe
library(biomaRt)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
Tcell_genes <- Tcell$gene
gene_locations <- getBM(
  attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position", "strand"),
  filters = "hgnc_symbol",
  values = Tcell_genes,
  mart = ensembl
)
colnames(gene_locations) <- c("gene", "gene_chr", "gene_start", "gene_end", "strand")

Tcell_full <- merge(Tcell, gene_locations, by = "gene")
Tcell_full$gene_chr <- paste0("chr", Tcell_full$gene_chr)
Tcell_full <- Tcell_full[!grepl("chrHSCHR|chrHG", Tcell_full$gene_chr), ]

Tcell_full <- tidyr::separate(Tcell_full, peak, into = c("chromosome", "start", "end"), sep = "[:-]")

Tcell_full$start <- as.numeric(Tcell_full$start)
Tcell_full$end <- as.numeric(Tcell_full$end)
Tcell_full$gene_start <- as.numeric(Tcell_full$gene_start)
Tcell_full$gene_end <- as.numeric(Tcell_full$gene_end)

Tcell_full$distanceToTSS <- round((Tcell_full$start + (Tcell_full$end - Tcell_full$start)/2) - 
  (Tcell_full$gene_start + (Tcell_full$gene_end - Tcell_full$gene_start)/2), 0)

Tcell_full_filter <- Tcell_full %>%
  filter(p < 0.05)
# update ggplot
ggplot(Tcell_full_filter, aes(x = distanceToTSS, y = beta)) +
  # Scatter plot with density-based coloring
  geom_bin2d(bins = 150) +  # Adjust `bins` for resolution
  scale_fill_gradient(low = "lightblue", high = "darkblue") +  # Density color scale
  
  # Horizontal line at y = 0
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  
  # Classic theme
  theme_classic() +
  theme(
    axis.line = element_line(),  # Add axis lines
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5)  # Adjust border line width
  ) +
  
  # X-axis customization
  scale_x_continuous(
    limits = c(-700000, 700000),  # Set x-axis range
    breaks = seq(-500000, 500000, by = 250000),  # Tick marks every 250k
    labels = c("-500K", "-250K", "0", "250K", "500K")  # Custom labels
  ) +
  
  # Y-axis customization
  scale_y_continuous(
    limits = c(-0.5, 1.2),  # Set y-axis range
    breaks = seq(-1, 1.5, by = 0.5)  # Tick marks every 0.5
  ) +
  
  # Labels and title
  labs(
    x = "Distance to TSS (bp)",
    y = "Effect Estimate (Beta)",
    title = "Effect Estimate vs. Distance to TSS"
  )




# Add the significant column
Tcell_annotated <- Tcell_full %>%
  mutate(significant = ifelse(-log10(p) >= -log10(0.05) & boot_basic_p < 0.05, "#008080", "gray"))  # 同时满足两个条件


# plot
ggplot(Tcell_annotated, aes(x = distanceToTSS, y = -log10(p), color = significant)) +
  geom_point(alpha = 0.6) +  # 使用显著性颜色绘制散点
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", size = 1) +  # 添加水平线
  geom_vline(xintercept = c(-50000, 50000), linetype = "dashed", color = "blue", size = 1) +  # 添加垂直线
  theme_classic() +
  theme(
    axis.line = element_line(),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5)
  ) +
  scale_color_identity() +  # 保持指定颜色
  scale_x_continuous(
    limits = c(-500000, 500000),
    breaks = seq(-500000, 500000, by = 250000),
    labels = c("-500K", "-250K", "0", "250K", "500K")
  ) +
  labs(
    x = "Distance to TSS (bp)",
    y = "-Log10(p_val)",
    title = "P-val vs. Distance to TSS"
  )



ggplot(Tcell_annotated, aes(x = distanceToTSS, y = -log10(boot_basic_p), color = significant)) +
  geom_point(alpha = 0.6) +  
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", size = 1) +  # 添加水平线
  geom_vline(xintercept = c(-50000, 50000), linetype = "dashed", color = "blue", size = 1) +  # 添加垂直线
  theme_classic() +
  theme(
    axis.line = element_line(),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5)
  ) +
  scale_color_identity() +  
  scale_x_continuous(
    limits = c(-500000, 500000),
    breaks = seq(-500000, 500000, by = 250000),
    labels = c("-500K", "-250K", "0", "250K", "500K")
  ) +
  labs(
    x = "Distance to TSS (bp)",
    y = "-Log10(boot_basic_p)",
    title = "boot_basic_p vs. Distance to TSS"
  )

#---------------------------------------------------------------------------
#Add a distance column
Tcell_annotated <- Tcell_annotated %>%
  mutate(
    distance_group = case_when(
      abs(distanceToTSS) <= 3000 ~ "≤3kb",
      abs(distanceToTSS) > 3000 & abs(distanceToTSS) <= 50000 ~ "3kb-50kb",
      abs(distanceToTSS) > 50000 ~ ">50kb"
    ),
    significance_label = ifelse(p <= 0.05, "sig", "not sig")
    )

# plot Volcano Plot
ggplot(Tcell_annotated, aes(x = beta, y = -log10(p), color = distance_group)) +
  geom_point(alpha = 0.6, size = 1.5) +  # 设置散点透明度和大小
  scale_color_manual(
    values = c("≤3kb" = "green", "3kb-50kb" = "blue", ">50kb" = "red")  # 自定义颜色
  ) +
  theme_classic() +
  labs(
    x = "Effect Size (Beta)",
    y = "-Log10(P-value)",
    color = "Distance from TSS",
    title = "Volcano Plot with Distance Grouping"
  ) +
  theme(
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.position = "right"
  )


# Make sure the levels are correct for plotting
Tcell_annotated$significance_label <- factor(
  Tcell_annotated$significance_label, 
  levels = c( "sig","not sig")
)

# summary data for number of `sig` and `not sig` 
summary_data <- Tcell_annotated %>%
  group_by(distance_group, significance_label) %>%
  summarize(count = n(), .groups = "drop") %>%
  mutate(distance_group_total = ave(count, distance_group, FUN = sum))  
# plot barplot
ggplot(summary_data, aes(x = distance_group, y = count, fill = significance_label)) +
  geom_bar(stat = "identity", position = "stack") +  
  scale_fill_manual(
    values = c("not sig" = "gray", "sig" = "#008080"), 
    labels = c("sig" = "Significant", "not sig" = "Not Significant")
  ) +
  coord_flip() +
  theme_classic() +
  labs(
    x = "Distance from TSS",
    y = "Count",
    fill = "Significance",
    title = "Barplot of Significance by Distance Group"
  ) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )
#------------------------------------------------------------------
# 假设 bin 大小为 5000 bp
bin_size <- 100

# 创建区间范围，确保正确覆盖正负值
breaks <- seq(-600000, 600000, by = bin_size)

# 创建格式化的标签
labels <- paste0(
  format(breaks[-length(breaks)], scientific = FALSE),  # 区间开始值
  " to ",
  format(breaks[-1], scientific = FALSE)               # 区间结束值
)

# 生成 bin，带格式化标签
Tcell_annotated2 <- Tcell_full_filter %>%
  mutate(
    bin = cut(distanceToTSS, 
              breaks = breaks, 
              include.lowest = TRUE, 
              right = FALSE, 
              labels = labels)  # 使用手动标签
  )

bin_centers <- sapply(levels(Tcell_annotated2$bin), function(b) {
  range <- as.numeric(unlist(strsplit(gsub(" to ", ",", b), ",")))  # 转换区间为数值
  mean(range)  # 计算区间中心点
})

# 创建 bin 对应的中心点映射表
bin_mapping <- data.frame(
  bin = levels(Tcell_annotated2$bin), 
  bin_center = bin_centers
)

# 将中心点映射回原始数据框
Tcell_annotated2 <- Tcell_annotated2 %>%
  left_join(bin_mapping, by = "bin")

# 计算每个 bin 的平均 beta 和 95% 置信区间
summary_data <- Tcell_annotated2 %>%
  group_by(bin, bin_center) %>%
  summarize(
    mean_beta = max(beta, na.rm = TRUE),
    ci_lower = max(beta, na.rm = TRUE) - 1.96 * sd(beta, na.rm = TRUE) / sqrt(n()),
    ci_upper = max(beta, na.rm = TRUE) + 1.96 * sd(beta, na.rm = TRUE) / sqrt(n())
  ) %>%
  ungroup()

# 绘制图形
ggplot(summary_data, aes(x = bin_center, y = mean_beta)) +
  geom_point(color = "#008080", alpha = 0.6) +  # 添加散点图
  geom_smooth(method = "loess", color = "blue", fill = "lightblue", alpha = 0.4) +  # 添加平滑线 + 
  scale_x_continuous(
    limits = c(-500000, 500000),  # 设置 x 轴范围
    breaks = seq(-500000, 500000, by = 250000),  # 以 250k 为间隔设置刻度
    labels = c("-500K", "-250K", "0", "250K", "500K")  # 标签显示
  )  +
  theme_classic() +
  labs(
    x = "Distance to TSS (bp)", 
    y = "Mean Effect Size (Beta)", 
    title = "Effect Size vs. Distance to TSS"
  )

#---------------------------------------------------------------------------
# read HiChIP Loop
Naive_H3K27ac_Loops_lifted <- process_loop_file(Naive_H3K27ac_Loops, chain)
Th17_H3K27ac_Loops_lifted <- process_loop_file(Th17_H3K27ac_Loops, chain)
Treg_H3K27ac_Loops_lifted <- process_loop_file(Treg_H3K27ac_Loops, chain)

Naive_H3K27ac_Loops_hg38 <- Naive_H3K27ac_Loops_lifted %>% na.omit()
Th17_H3K27ac_Loops_hg38 <- Th17_H3K27ac_Loops_lifted %>% na.omit()
Treg_H3K27ac_Loops_hg38 <- Treg_H3K27ac_Loops_lifted %>% na.omit()

# save as bedpe
write.table(Naive_H3K27ac_Loops_hg38, file = "/Users/zoechen/Documents/Study/UCSF/Fall_2024/BMI_206_Biostats/Projects/Naive_H3K27ac_Loops_hg38.bedpe", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(Th17_H3K27ac_Loops_hg38, file = "/Users/zoechen/Documents/Study/UCSF/Fall_2024/BMI_206_Biostats/Projects/Th17_H3K27ac_Loops_hg38.bedpe", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(Treg_H3K27ac_Loops_hg38, file = "/Users/zoechen/Documents/Study/UCSF/Fall_2024/BMI_206_Biostats/Projects/Treg_H3K27ac_Loops_hg38.bedpe", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
#---------------------------------------------------------------------------
library(GenomicRanges)
library(regioneR)
Tcell_distal <- Tcell_annotated %>%
  filter(distance_group == ">50kb" & significance_label == "sig")
Tcell_short <- Tcell_annotated %>%
  filter(distance_group == "3kb-50kb" & significance_label == "sig")

tcell_distal <- Tcell_distal[,2:4] %>% distinct()
tcell_short <- Tcell_short[,2:4] %>% distinct()


ENCODE_Tcell=toGRanges(read.table("/Users/zoechen/Documents/Study/UCSF/Fall_2024/BMI_206_Biostats/Projects/ENCODE_Tcell/recurrent_peaks.merge.bed", sep='\t'))
tcell_distal=toGRanges(tcell_distal)
tcell_short=toGRanges(tcell_short)


numOverlaps(tcell_distal, ENCODE_Tcell,count.once=TRUE)
numOverlaps(tcell_short, ENCODE_Tcell,count.once=TRUE)

tcell_distal <- tcell_distal[width(tcell_distal) > 0]
ENCODE_Tcell <- ENCODE_Tcell[width(ENCODE_Tcell) > 0]
# Perform enrichment analysis for Tcell_distal
set.seed(10)
library(BSgenome.Hsapiens.UCSC.hg38)
#Generate a set of random genomic regions of the same size as hits1
genome(tcell_distal) <- "hg38"  # Replace "hg19" with the correct genome
genome(ENCODE_Tcell_GR) <- "hg38"
genome(tcell_short) <- "hg38"

options(mc.cores = 1)

distal_enrichment <- overlapPermTest(A = tcell_distal, B = ENCODE_Tcell, ntimes = 1000, genome = "hg38")

# Perform enrichment analysis for Tcell_short
short_enrichment <- overlapPermTest(A = tcell_short, B = ENCODE_Tcell, ntimes = 1000, genome = "hg38")
# Set up a multi-panel layout (1 row, 2 columns)
par(mfrow = c(1, 2))

# Plot the distal enrichment
plot(distal_enrichment, main = "Distal Enrichment")

# Plot the short enrichment
plot(short_enrichment, main = "Short Enrichment")

write.table(tcell_distal, file = "/Users/zoechen/Documents/Study/UCSF/Fall_2024/BMI_206_Biostats/Projects/tcell_distal.bed", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(tcell_short, file = "/Users/zoechen/Documents/Study/UCSF/Fall_2024/BMI_206_Biostats/Projects/tcell_short.bed", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

# Define the contingency table
contingency_table <- matrix(c(4953, 1998, 13627, 15948), nrow = 2, byrow = TRUE)
colnames(contingency_table) <- c("Overlap", "Not Overlap")
rownames(contingency_table) <- c("3kb - 50kb", ">50kb")

# Perform Fisher's exact test
fisher.test(contingency_table)
#------------------------------------------------------------------------------------
#transfer to pair bedpe
Tcell_promo <- Tcell_annotated %>%
  filter(distance_group == "≤3kb" & significance_label == "sig")
Tcell_promo_bedpe <- Tcell_promo[,c(10:12, 2:4)] %>% distinct()
Tcell_short_bedpe <- Tcell_short[,c(10:12, 2:4)] %>% distinct()
Tcell_distal_bedpe <- Tcell_distal[,c(10:12, 2:4)] %>% distinct()

Tcell_promo_bedpe_adj <- Tcell_promo_bedpe %>%
  mutate(gene_start = gene_start - 500,
         gene_end = gene_start + 500,
         mid = start+round((end - start)/2,0),
         start = mid - 500, 
         end = mid + 500) %>% dplyr::select(!mid)

Tcell_short_bedpe_adj <- Tcell_short_bedpe %>%
  mutate(gene_start = gene_start - 500,
         gene_end = gene_start + 500,
         mid = start+round((end - start)/2,0),
         start = mid - 500, 
         end = mid + 500) %>% dplyr::select(!mid)

Tcell_distal_bedpe_adj <- Tcell_distal_bedpe %>%
  mutate(gene_start = gene_start - 500,
         gene_end = gene_start + 500,
         mid = start+round((end - start)/2,0),
         start = mid - 500, 
         end = mid + 500) %>% dplyr::select(!mid)

write.table(Tcell_distal_bedpe_adj, file = "/Users/zoechen/Documents/Study/UCSF/Fall_2024/BMI_206_Biostats/Projects/tcell_distal_adj.bedpe", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(Tcell_short_bedpe_adj, file = "/Users/zoechen/Documents/Study/UCSF/Fall_2024/BMI_206_Biostats/Projects/tcell_short_adj.bedpe", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(Tcell_promo_bedpe_adj, file = "/Users/zoechen/Documents/Study/UCSF/Fall_2024/BMI_206_Biostats/Projects/tcell_promo_adj.bedpe", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)


T_Cell_HiChIP <- read.table("/Users/zoechen/Documents/Study/UCSF/Fall_2024/BMI_206_Biostats/Projects/enriched_pixels_10000.bedpe", sep='\t', header=TRUE)
T_Cell_HiChIP_f <- T_Cell_HiChIP %>%
  filter(fdrBL < 0.1 & fdrDonut < 0.1 & fdrH < 0.1 & fdrV < 0.1)


# Define the contingency table
contingency_table <- matrix(c(420, 5854, 3578, 22945), 
                            nrow = 2, 
                            dimnames = list(c("3kb - 50kb", ">50kb"), c("Overlap", "Not Overlap")))

# Perform Fisher's Exact Test
fisher.test(contingency_table)



#--------------------------------------------------------------
# Create a confusion matrix
confusion_matrix <- matrix(
  c(420, 6274, 5854, 22945), 
  nrow = 2, 
  byrow = TRUE,
  dimnames = list(
    "Distance Group" = c("3kb - 50kb", ">50kb"),
    "Overlap Status" = c("Overlap", "Not Overlap")
  )
)

# Print the confusion matrix
print(confusion_matrix)

library(RColorBrewer)
library(plotrix)

heatmap_colors <- colorRampPalette(brewer.pal(9, "Blues"))(256)

heatmap_obj <- heatmap(confusion_matrix, Rowv = NA, Colv = NA, 
                       col = heatmap_colors, scale = "none", 
                       main = "Confusion Matrix Heatmap", 
                       xlab = "Overlap Status", ylab = "Distance Group",
                       keep.dendro = FALSE)

# Add annotations (numbers)
text_labels <- as.character(confusion_matrix)
x_coords <- rep(1:ncol(confusion_matrix), each = nrow(confusion_matrix))
y_coords <- rep(1:nrow(confusion_matrix), times = ncol(confusion_matrix))

# Add numbers to the heatmap
text(x_coords, rev(y_coords), labels = text_labels, col = "white", font = 2)

#--------------------------------------------------------------

# Create a new confusion matrix
confusion_matrix <- matrix(
  c(4953, 13627, 1998, 15948), 
  nrow = 2, 
  byrow = TRUE,
  dimnames = list(
    "Distance Group" = c("3kb - 50kb", ">50kb"),
    "Overlap Status" = c("Overlap", "Not Overlap")
  )
)

# Print the confusion matrix
print(confusion_matrix)

# Optional: Visualize the confusion matrix using a blue gradient heatmap with annotations
library(RColorBrewer)

# Calculate percentages
percent_matrix <- 100 * confusion_matrix / sum(confusion_matrix)


heatmap_colors <- colorRampPalette(brewer.pal(9, "Blues"))(256)

# Adjust margins for axis labels
par(mar = c(5, 6, 4, 2) + 0.1)  # Bottom, left, top, right margins

# Draw heatmap with percentage-based colors
heatmap_obj <- heatmap(percent_matrix, Rowv = NA, Colv = NA, 
                       col = heatmap_colors, scale = "none", 
                       main = "Confusion Matrix Heatmap (Percentage)", 
                       xlab = "Overlap Status", ylab = "Distance Group",
                       cex.main = 1.5, cex.lab = 1.2, cex.axis = 1.2,
                       keep.dendro = FALSE)

