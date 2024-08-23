# Code to create the mean and median abundance values for each cluster at each clustering level.

# Load packages
library(Seurat)

# Define paths as variables
data_path <- "data/humanHYPOMAP_spatial_c2l_240621.RDS"
output_dir <- "tables/"

# Load dataset
merged <- readRDS(data_path)

# Create function to calculate average abundance per cluster
cell2loc_avg_abundance <- function(data, assay, clusters, average_function) {
  x <- data.frame(t(data@assays[[assay]]@data))
  y <- data@meta.data[[clusters]]
  print('Calculating average cell abundance per region..')
  x$cluster <- y
  avex <- aggregate(.~cluster, x, FUN = average_function)
  rownames(avex) <- avex$cluster
  avex <- avex[, -c(1)]
  return(avex)
}

# Calculate and save average (mean) abundance per cluster
for (i in c('C1', 'C2', 'C3', 'C4')) {
  for (j in c('region_cluster_NN_27_res_0.5', 'NN27_res_0.5_grouped', 'NN27_res_0.5_named')) {
    print(i)
    print(j)
    x <- cell2loc_avg_abundance(data = merged, assay = i, clusters = j, average_function = mean)
    write.table(x, paste0(output_dir, "all_clusters_avg_abundance_", i, "_regions_", j, "_mean_240621.txt"), sep = '\t')
    print('saved')
  }
}

# Calculate and save median abundance per cluster
for (i in c('C1', 'C2', 'C3', 'C4')) {
  for (j in c('region_cluster_NN_27_res_0.5', 'NN27_res_0.5_grouped', 'NN27_res_0.5_named')) {
    print(i)
    print(j)
    x <- cell2loc_avg_abundance(data = merged, assay = i, clusters = j, average_function = median)
    write.table(x, paste0(output_dir, "all_clusters_avg_abundance_", i, "_regions_", j, "_median_240621.txt"), sep = '\t')
    print('saved')
  }
}
