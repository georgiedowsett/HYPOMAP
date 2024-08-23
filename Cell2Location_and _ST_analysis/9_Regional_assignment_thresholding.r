#### Regional annotation of nucseq clusters based off C3 cell abundance mappings 
 
library(Seurat)
library(ggplot2)
library(tidyr)
library(dplyr)

# Define paths as variables
base_path <- "~/path/to/your/base/directory/"
tables_path <- paste0(base_path, "tables/")
nucseq_objects_path <- paste0(base_path, "data/nucseq-objects/")
output_file <- paste0(tables_path, "nucseq_clusters_region_assignment.txt")

# Load C3 ungrouped median table 
avg_abundance_C3_ungrouped <- data.table::fread(paste0(tables_path, "all_clusters_avg_abundance_C3_regions_NN27_res_0.5_named_median_240621.txt"), data.table = F)
avg_abundance_C3_grouped <- data.table::fread(paste0(tables_path, "all_clusters_avg_abundance_C3_regions_NN27_res_0.5_grouped_median_240621.txt"), data.table = F)

# Make a grouped to ungrouped df
region_grouping <- data.frame(region_name_ungrouped = avg_abundance_C3_ungrouped[, 1])
region_grouping$region_name_grouped <- stringr::str_remove(region_grouping$region_name_ungrouped, pattern = " [0-9]+")

# Load the color scheme for the grouped simplified anno
clusters_grouped_colors <- data.table::fread(paste0(tables_path, "spatial_clusters_colours.txt"), data.table = F)
region_color_mapping <- as.character(clusters_grouped_colors$col)
names(region_color_mapping) <- as.character(clusters_grouped_colors$V1)

# Add to above also
region_grouping <- dplyr::left_join(region_grouping, clusters_grouped_colors, by = c("region_name_grouped" = "V1"))
region_color_mapping_ungrp <- as.character(region_grouping$col)
names(region_color_mapping_ungrp) <- as.character(region_grouping$region_name_ungrouped)

# Load the avg c4 abundances per spatial cluster
avg_abundances_c4 <- data.table::fread(paste0(tables_path, "all_clusters_avg_abundance_C4_regions_NN27_res_0.5_named_mean_240621.txt"), data.table = F)
avg_abundances_c3 <- data.table::fread(paste0(tables_path, "all_clusters_avg_abundance_C3_regions_NN27_res_0.5_named_mean_240621.txt"), data.table = F)
avg_abundances_c2 <- data.table::fread(paste0(tables_path, "all_clusters_avg_abundance_C2_regions_NN27_res_0.5_named_mean_240621.txt"), data.table = F)

# Make long format
avg_abundances_c4_long <- avg_abundances_c4 %>% tidyr::gather(key = "cluster", value = "abundance", -V1) %>% dplyr::rename(region = V1) %>% dplyr::mutate(abundance = as.numeric(abundance))
avg_abundances_c3_long <- avg_abundances_c3 %>% tidyr::gather(key = "cluster", value = "abundance", -V1) %>% dplyr::rename(region = V1) %>% dplyr::mutate(abundance = as.numeric(abundance))
avg_abundances_c2_long <- avg_abundances_c2 %>% tidyr::gather(key = "cluster", value = "abundance", -V1) %>% dplyr::rename(region = V1) %>% dplyr::mutate(abundance = as.numeric(abundance))

all_abundances <- list(
  C2 = avg_abundances_c2_long,
  C3 = avg_abundances_c3_long,
  C4 = avg_abundances_c4_long
)

human_hypo_combined <- readRDS(paste0(nucseq_objects_path, "human_hypo_combined.rds"))
human_hypo_combined@meta.data <- human_hypo_combined@meta.data[, 1:27] 

##########
### Check assigned region cluster by taking max after adjusting + filter for MAD
##########

all_abundances_adj_max_mad_stats <- list()
for (i in 1:length(all_abundances)) {
  target_level <- names(all_abundances)[i]
  print(target_level)
  # with adjust
  avg_abundances_long_stats <- all_abundances[[i]] %>%
    dplyr::group_by(region) %>% 
    dplyr::mutate(abundance_adj = abundance - median(abundance)) %>%
    dplyr::group_by(cluster) %>% 
    dplyr::mutate(
      median_ab = median(abundance_adj),
      max_ab = max(abundance_adj),
      max_median_ratio = max_ab / median_ab,
      current_max_ratio = abundance_adj / max_ab,
      abundance_mad = mad(abundance_adj),
      mad_x = abundance_adj / abundance_mad
    ) %>%
    dplyr::slice_max(order_by = abundance_adj, n = 1)
  
  all_abundances_adj_max_mad_stats[[target_level]] <- avg_abundances_long_stats
  
  add_regions <- avg_abundances_long_stats %>%
    dplyr::left_join(region_grouping[, 1:2], by = c("region" = "region_name_ungrouped")) %>%
    dplyr::select(cluster, region_name_grouped) %>% 
    dplyr::distinct() %>% 
    dplyr::mutate(cluster = stringr::str_replace(cluster, "\\.", "-"))
  colnames(add_regions) <- c(target_level, paste0(target_level, "_region_adj_max"))
  temp_meta <- human_hypo_combined@meta.data %>% dplyr::left_join(add_regions) %>% as.data.frame()
  rownames(temp_meta) <- temp_meta$Cell_ID
  human_hypo_combined@meta.data <- temp_meta
}

# Create a table with the max region assignment for each C3 neuronal and C2 non-neuronal cluster 
c3 <- all_abundances_adj_max_mad_stats$C3
c2 <- all_abundances_adj_max_mad_stats$C2

# Read in cluster levels from snRNAseq
nucseq_clusters <- read.delim(paste0(nucseq_objects_path, "human_hypo_combined_cluster_annotation.txt"))

glia_c3 <- c('C3-1', 'C3-12', 'C3-2', 'C3-10', 'C3-3', 'C3-4', 'C3-5', 'C3-13', 'C3-6', 'C3-11', 'C3-7', 'C3-8', 'C3-9', 'C3-14', 'C3-15', 'C3-18', 'C3-21', 'C3-16', 'C3-19', 'C3-17', 'C3-20', 'C3-141', 'C3-146', 'C3-154', 'C3-155', 'C3-147', 'C3-142', 'C3-143', 'C3-148', 'C3-149', 'C3-156', 'C3-150', 'C3-144', 'C3-151', 'C3-152', 'C3-145', 'C3-153')
glia_c3 <- stringr::str_replace(glia_c3, "-", "\\.")

glia_c2 <- c('C2-1', 'C2-2', 'C2-3', 'C2-4', 'C2-5', 'C2-6', 'C2-7', 'C2-8', 'C2-48', 'C2-49', 'C2-50', 'C2-51', 'C2-52')
glia_c2 <- stringr::str_replace(glia_c2, "-", "\\.")

region_assign <- c3[!c3$cluster %in% glia_c3, ]
c2_glia <- c2[c2$cluster %in% glia_c2, ]

region_assign <- rbind(region_assign, c2_glia)

region_assign <- region_assign[, c(1:4, 9:10)]

nas <- c('C3.50', 'C3.53', 'C3.97', 'C3.99', 'C3.23', 'C3.26', 'C3.30', 'C3.33', 'C3.39', 'C3.40', 'C3.110')

# Creating a table with all the correct manually assigned and auto-assigned region assignments for the C3.C2 levels. 
region_assign1 <- region_assign %>%
  dplyr::left_join(region_grouping[, 1:2], by = c("region" = "region_name_ungrouped")) %>%
  dplyr::distinct() %>%
  mutate(region_threshold = ifelse(mad_x > 10, region_name_grouped, NA)) %>%
  mutate(region_threshold_manual = case_when(
    cluster %in% nas ~ NA,
    cluster == "C3.48" ~ "LH",
    cluster == "C3.34" ~ "LPOA",
    cluster == "C3.90" ~ "MAM",
    cluster == "C3.83" ~ "MAM",
    cluster == "C3.79" ~ "MPOA",
    cluster == "C3.80" ~ "POA",
    cluster == "C3.100" ~ "Thalamaus",
    TRUE ~ region_threshold
  ))

# Function to get region_assign1 data based on cluster value
get_region_data <- function(cluster_value) {
  region_assign1 %>% filter(cluster == cluster_value)
}

# Initialize an empty list to store the results
results <- list()

# Loop through each row of nucseq_clusters
for (i in 1:nrow(nucseq_clusters)) {
  C2_value <- gsub("-", ".", nucseq_clusters$C2[i])
  C3_value <- gsub("-", ".", nucseq_clusters$C3[i])
  C4_value <- nucseq_clusters$C4[i]
  
  # Check if C2_value is present in region_assign1$cluster
  region_data <- get_region_data(C2_value)
  
  if (nrow(region_data) == 0) {
    # If not found, check if C3_value is present
    region_data <- get_region_data(C3_value)
  }
  
  if (nrow(region_data) > 0) {
    # If found, add C2, C3, and C4 columns to region_data
    region_data <- region_data %>%
      mutate(C2 = nucseq_clusters$C2[i], C3 = nucseq_clusters$C3[i], C4 = C4_value)
    
    # Append the result to the list
    results[[length(results) + 1]] <- region_data
  }
}

# Combine all results into a single data frame
final_result <- bind_rows(results)

# Select the desired columns
final_result <- final_result %>%
  select(C2, C3, C4, abundance, abundance_adj, abundance_mad, mad_x, region_name_grouped, region_threshold, region_threshold_manual)

# View the result
head(final_result)

Idents(human_hypo_combined) <- 'C4'
current_idents <- Idents(human_hypo_combined)

# Create a data frame to map the cluster identities to region_threshold_manual values
cluster_mapping <- final_result %>%
  select(C4, region_threshold_manual) %>%
  distinct()
metadata_vector <- setNames(cluster_mapping$region_threshold_manual, cluster_mapping$C4)
mapped_metadata <- metadata_vector[as.character(current_idents)]
human_hypo_combined <- AddMetaData(human_hypo_combined, metadata = mapped_metadata, col.name = "region_threshold_manual")
head(human_hypo_combined@meta.data)

for (i in 1:length(all_abundances)) {
  target_level <- names(all_abundances)[i]
  print(target_level)
  # with adjust
  avg_abundances_long_stats <- all_abundances[[i]] %>%
    dplyr::group_by(region) %>%
    dplyr::mutate(abundance_adj = abundance - median(abundance)) %>%
    dplyr::group_by(cluster) %>%
    dplyr::mutate(median_ab = median(abundance_adj), max_ab = max(abundance_adj), max_median_ratio = max_ab / median_ab, current_max_ratio = abundance_adj / max_ab, abundance_mad = mad(abundance_adj), mad_x = abundance_adj / abundance_mad) %>%
    dplyr::slice_max(order_by = abundance_adj, n = 1) %>%
    dplyr::filter(mad_x >= 10)
  
  all_abundances_adj_max_mad_stats[[target_level]] <- avg_abundances_long_stats
  
  add_regions <- avg_abundances_long_stats %>% 
    dplyr::left_join(region_grouping[, 1:2], by = c("region" = "region_name_ungrouped")) %>%
    dplyr::select(cluster, region_name_grouped) %>% 
    dplyr::distinct() %>% 
    dplyr::mutate(cluster = stringr::str_replace(cluster, "\\.", "-"))
  colnames(add_regions) <- c(target_level, paste0(target_level, "_region_adj_max_cleaned"))
  temp_meta <- human_hypo_combined@meta.data %>% dplyr::left_join(add_regions) %>% as.data.frame()
  rownames(temp_meta) <- temp_meta$Cell_ID
  human_hypo_combined@meta.data <- temp_meta
}

# Make plots
colnames(human_hypo_combined@meta.data)
sapply(all_abundances_adj_max_mad_stats, nrow)
p2 <- DimPlot(human_hypo_combined, group.by = c("C2_region_adj_max_cleaned"), reduction = "umap_scvi_hypo", label = F, shuffle = TRUE, raster.dpi = c(1536, 1536), pt.size = 1.3) + scale_color_manual(values = region_color_mapping)
p3 <- DimPlot(human_hypo_combined, group.by = c("C3_region_adj_max_cleaned"), reduction = "umap_scvi_hypo", label = F, shuffle = TRUE, raster.dpi = c(1536, 1536), pt.size = 1.3) + scale_color_manual(values = region_color_mapping)
p4 <- DimPlot(human_hypo_combined, group.by = c("C4_region_adj_max_cleaned"), reduction = "umap_scvi_hypo", label = F, shuffle = TRUE, raster.dpi = c(1536, 1536), pt.size = 1.3) + scale_color_manual(values = region_color_mapping)
cowplot::plot_grid(plotlist = list(p2, p3, p4), nrow = 2)

# Write the final result to a file
final_result1 <- final_result[, -1]
write.table(final_result1, output_file, sep = '\t', row.names = F)
