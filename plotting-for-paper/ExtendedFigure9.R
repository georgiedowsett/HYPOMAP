library(Seurat)
library(ggplot2)
library(cowplot)
library(viridis)
library(viridisLite)
library(scales)
library(patchwork)

setwd('/path/to/humanHYPOMAP/')

# Load in data
merged <- readRDS('data/humanHYPOMAP_spatial_c2l_240621.RDS')

# Set image order
ordered_images <- c('slice5A', 'sliceC1A', 'slice6B',  'sliceB1B', 'slice2B', 'slice4B', 'slice7A', 'slice3A', 'slice8B' )

# Set output directory
output_dir <- 'figs_241004/'

my_spatial_palette <- colorRampPalette(c(viridis(7, begin = 0.2)), bias = 2)(n = 256)

huhy_plots <- function(data, gene, percentile = NA, image.alpha = 0, assay.use = data@active.assay, ncol = NULL, image = NULL, pt.size.factor = 1.4, crop = FALSE, stroke = NA) {
  if (is.na(percentile)){
    gene_expression <- data@assays[[assay.use]]@data[gene,]
    p <- SpatialFeaturePlot(data, gene, image = image, image.alpha = image.alpha, pt.size.factor = pt.size.factor, crop = crop, ncol = ncol, stroke = stroke) & 
      scale_fill_gradientn(colors = my_spatial_palette, oob=squish) &
      theme(plot.title = element_blank())
    p <- wrap_plots(p)
    
  } else {
    gene_expression <- data@assays[[assay.use]]@data[gene,]
    scale <- quantile(gene_expression, percentile)
    scale <- as.numeric(round(scale, 1))
    p <- SpatialFeaturePlot(data, gene, image = image, image.alpha = image.alpha, pt.size.factor = pt.size.factor, crop = crop, ncol = ncol, stroke = stroke) & 
      scale_fill_gradientn(limits = c(0, scale), breaks = c(0,scale), colors = my_spatial_palette, oob=squish) &
      theme(plot.title = element_blank())
    p <- wrap_plots(p) + plot_layout(guide= 'collect') &
      theme(legend.position = 'top', legend.direction = 'horizontal', legend.key.size = unit(2, 'cm'), legend.text = element_text(size = 40), legend.title = element_text(size = 40))
  }
  print(p)
}

# Plotting of cell2location cell abundances
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
                         breaks = c(0,scale), 
                         colors = my_cell2loc_palette, 
                         oob=squish)
}

# Use the modified function
png(paste0(output_dir, 'humanHYPOMAP_extendedfig9_GLP1R.png'), width = 35, height = 10, units = 'in', res = 400)
huhy_plots(data = merged, gene = 'GLP1R', 
           percentile = 1, 
           image.alpha = 0, 
           assay.use = 'Spatial', 
           ncol = 9, 
           image = ordered_images, 
           pt.size.factor = 1.45, 
           crop = FALSE)
dev.off()

png(paste0(output_dir, 'humanHYPOMAP_extendedfig9_GIPR.png'), width = 35, height = 10, units = 'in', res = 400)
huhy_plots(data = merged, gene = 'GIPR', 
           percentile = 1, 
           image.alpha = 0, 
           assay.use = 'Spatial', 
           ncol = 9, 
           image = ordered_images, 
           pt.size.factor = 1.45, 
           crop = FALSE)
dev.off()

# Spatial mapping plots 
# Read in clusters
nucseq_clusters <- read.delim('data/nucseq-objects/human_hypo_combined_cluster_annotation.txt')
incretin <- nucseq_clusters[nucseq_clusters$C4 %in% c('C4-158', 'C4-373', 'C4-312', 'C3-12'),]

slide <- 'slice2B'
DefaultAssay(merged) <- 'C4'
title <- incretin$C4_named[4]
png(paste0(output_dir, 'humanHYPOMAP_extendedfig9_C4-158.png'), width = 15, height = 10, units = 'in', res = 400)
cell2loc_plots(merged, 
               cluster = 'C4.158', 
               scale_cutoff = 0.995, 
               assay = 'C4', 
               image_order = slide, 
               pt.size.factor = 1.5) &
  theme(legend.position = 'right', 
        legend.key.size = unit(1.5, 'cm'), 
        legend.text = element_text(size = 30), 
        legend.title = element_blank(),
        title = element_text(size = 30)) &
  ggtitle(title)
dev.off()

title <- incretin$C4_named[3]
png(paste0(output_dir, 'humanHYPOMAP_extendedfig9_C4-373.png'), width = 15, height = 10, units = 'in', res = 400)
cell2loc_plots(merged, 
               cluster = 'C4.373', 
               scale_cutoff = 0.992, 
               assay = 'C4', 
               image_order = slide, 
               pt.size.factor = 1.5) &
  theme(legend.position = 'right', 
        legend.key.size = unit(1.5, 'cm'), 
        legend.text = element_text(size = 30), 
        legend.title = element_blank(),
        title = element_text(size = 30)) &
  ggtitle(title)
dev.off()

title <- incretin$C4_named[2]
png(paste0(output_dir, 'humanHYPOMAP_extendedfig9_C4-312.png'), width = 15, height = 10, units = 'in', res = 400)
cell2loc_plots(merged, 
               cluster = 'C4.312', 
               scale_cutoff = 0.9, 
               assay = 'C4', 
               image_order = slide, 
               pt.size.factor = 1.5) &
  theme(legend.position = 'right', 
        legend.key.size = unit(1.5, 'cm'), 
        legend.text = element_text(size = 30), 
        legend.title = element_blank(),
        title = element_text(size = 30)) &
  ggtitle(title)
dev.off()

DefaultAssay(merged) <- 'C3'
title <- incretin$C4_named[1]
png(paste0(output_dir, 'humanHYPOMAP_extendedfig9_C3-12.png'), width = 15, height = 10, units = 'in', res = 400)
cell2loc_plots(merged, 
               cluster = 'C3.12', 
               scale_cutoff = 0.992, 
               assay = 'C3', 
               image_order = slide, 
               pt.size.factor = 1.5) &
  theme(legend.position = 'right', 
        legend.key.size = unit(1.5, 'cm'), 
        legend.text = element_text(size = 30), 
        legend.title = element_blank(),
        title = element_text(size = 30)) &
  ggtitle(title)
dev.off()
