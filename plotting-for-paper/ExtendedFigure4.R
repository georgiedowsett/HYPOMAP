library(Seurat)
library(ggplot2)
library(cowplot)
library(patchwork)
library(viridis)
library(viridisLite)
library(scales)

setwd('/path/to/humanHYPOMAP/')
# Load data
merged <- readRDS('data/humanHYPOMAP_spatial_c2l_240621.RDS')

# Set image order
ordered_images <- c('slice5A', 'slice6B', 'sliceC1A', 'sliceB1B', 'slice2B', 'slice4B', 'slice7A', 'slice3A', 'slice8B')

# Set output directory
output_dir <- 'figs_241004/'

# Cell2location plotting function
my_cell2loc_palette <- colorRampPalette(c(magma(7, begin = 0.2)), bias = 1)(n = 256)
cell2loc_plots <- function(data, cluster, scale_cutoff = 0.992, assay = data@active.assay, image_order = ordered_images, slide_metadata = 'captureArea', stroke = NA, pt.size.factor = 1.3) {
  gene_expression <- t(data@assays[[assay]]@data[cluster, , drop = FALSE])
  slide_info <- as.matrix(data@meta.data[[slide_metadata]])
  gene_slide <- as.data.frame(cbind(gene_expression, slide_info))
  quantiles <- aggregate(gene_expression ~ V2, data = gene_slide, FUN = function(x) quantile(x, probs = scale_cutoff))
  x <- max(quantiles[,2])
  if (x < 0.05) {
    scale <- signif(x, digits = 2)
  } else {
    scale <- round(x, digits = 1)
  }
  SpatialFeaturePlot(data, cluster, image.alpha = 0, images = image_order, crop = F, pt.size.factor = pt.size.factor, stroke = stroke) & 
    plot_layout(nrow = 1) & 
    scale_fill_gradientn(limits = c(0, scale), 
                         breaks = c(0, scale), 
                         colors = my_cell2loc_palette, 
                         oob = squish)
}

# Read in clusters and identify astrocyte cluster numbers
nucseq_clusters <- read.delim('data/nucseq-objects/human_hypo_combined_cluster_annotation.txt')
astro <- nucseq_clusters[nucseq_clusters$C1 == 'C1-1',]

# Let's look at the C2 spatial clusters 
DefaultAssay(merged) <- 'C2'
plots <- list()

for (i in 1:length(unique(astro$C2))) {
  # Replace '-' with '.' in the cluster name
  cluster <- gsub("-", ".", astro$C2[i])
  title <- astro$C2_named[i]
  p <- cell2loc_plots(merged, cluster = cluster, image_order = 'slice3A', assay = 'C2', stroke = NA, pt.size.factor = 1.44) &
    theme(legend.position = "right", 
          legend.box = "vertical", 
          legend.margin = margin(t = 10, r = 5, b = -20, l = -5), 
          legend.text = element_text(size = 15), 
          legend.title = element_blank(), 
          legend.key.height = unit(0.5, "cm"), 
          plot.title = element_text(size = 10))
  p <- p + ggtitle(title)
  plots <- append(plots, list(p))
}

pdf(paste0(output_dir, 'humanHYPOMAP_extfig4a.pdf'), width = 10, height = 10)
plot_grid(plotlist = plots, ncol = 3, nrow = 1, rel_widths = c(1.3, 1.3, 1.3, 1.3))
dev.off()

png(paste0(output_dir, 'humanHYPOMAP_extfig4a.png'), width = 10, height = 5, res = 400, units = 'in')
plot_grid(plotlist = plots, ncol = 3, nrow = 1, rel_widths = c(1.3, 1.3, 1.3, 1.3))
dev.off()


# Spatial feature plots showing Tanycyte Markers 
my_spatial_palette <- colorRampPalette(c(viridis(7, begin = 0.2)), bias = 2)(n = 256)
plot_list <- list()
for (i in c('CRYM', 'FRZB', 'DIO2', 'FZD5', 'STOML3', 'LPAR3')) {
  plot <- SpatialFeaturePlot(merged, 
                             features = i, 
                             images = 'slice2B', 
                             crop = F, 
                             image.alpha = 0, 
                             pt.size.factor = 1.3,
                             stroke = NA) & 
    scale_fill_gradientn(colors = my_spatial_palette, 
                         oob = squish) &
    theme(legend.position = "right", 
          legend.box = "vertical", 
          legend.margin = margin(t = 0, 
                                 r = 20, 
                                 b = -100, 
                                 l = 0), 
          legend.text = element_text(size = 30), 
          legend.title = element_blank(), 
          legend.key.height = unit(1, "cm"), 
          plot.title = element_text(face = "italic", 
                                    size = 40, 
                                    hjust = 0.5, 
                                    margin = margin(b = 10))) & 
    labs(title = i)
  plot_list[[i]] <- plot
}
plots <- wrap_plots(plot_list, ncol = 6)

png(paste0(output_dir, 'humanHYPOMAP_extfig4d.png'), width = 30, height = 5, units = 'in', res = 400)
plots
dev.off()
