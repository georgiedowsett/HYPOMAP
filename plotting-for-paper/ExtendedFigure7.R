library(Seurat)
library(ggplot2)
library(cowplot)
library(scales)
library(viridis)
library(viridisLite)
library(patchwork)

setwd('/path/to/humanHYPOMAP/')

# Load data
merged <- readRDS('data/humanHYPOMAP_spatial_c2l_240621.RDS')

# Set image order
ordered_images <- c('slice5A', 'sliceC1A', 'slice6B', 'sliceB1B', 'slice2B', 'slice4B', 'slice7A', 'slice3A', 'slice8B')

# Set output directory
output_dir <- 'figs_241004/'

my_spatial_palette <- colorRampPalette(c(viridis(7, begin = 0.2)), bias = 2)(n = 256)

huhy_plots <- function(data, gene, percentile = NA, image.alpha = 0, assay.use = data@active.assay, ncol = NULL, image = NULL, pt.size.factor = 1.4, crop = FALSE, stroke = NA){
  if (is.na(percentile)){
    gene_expression <- data@assays[[assay.use]]@data[gene,]
    p <- SpatialFeaturePlot(data, gene, image = image, image.alpha = image.alpha, pt.size.factor = pt.size.factor, crop = crop, ncol = ncol, stroke = stroke) & 
      scale_fill_gradientn(colors = my_spatial_palette, oob = squish) &
      theme(plot.title = element_blank())
    p <- wrap_plots(p)
    
  } else {
    gene_expression <- data@assays[[assay.use]]@data[gene,]
    scale <- quantile(gene_expression, percentile)
    scale <- as.numeric(round(scale, 1))
    p <- SpatialFeaturePlot(data, gene, image = image, image.alpha = image.alpha, pt.size.factor = pt.size.factor, crop = crop, ncol = ncol, stroke = stroke) & 
      scale_fill_gradientn(limits = c(0, scale), breaks = c(0, scale), colors = my_spatial_palette, oob = squish) &
      theme(plot.title = element_blank())
    p <- wrap_plots(p) + plot_layout(guide = 'collect') &
      theme(legend.position = 'top', legend.direction = 'horizontal', legend.key.size = unit(2, 'cm'), legend.text = element_text(size = 40), legend.title = element_text(size = 40))
  }
  print(p)
}

pdf(paste0(output_dir, 'humanHYPOMAP_extendedfig7_POMC.pdf'), width = 10, height = 10)
huhy_plots(data = merged, gene = 'POMC', 
           percentile = 1, 
           image.alpha = 0, 
           assay.use = 'Spatial', 
           ncol = 3, 
           image = ordered_images, 
           pt.size.factor = 1.5, 
           crop = FALSE) 
dev.off()

pdf(paste0(output_dir, 'humanHYPOMAP_extendedfig7_AGRP.pdf'), width = 10, height = 10)
huhy_plots(data = merged, gene = 'AGRP', 
           percentile = 1, 
           image.alpha = 0, 
           assay.use = 'Spatial', 
           ncol = 3, 
           image = ordered_images, 
           pt.size.factor = 1.5, 
           crop = FALSE) 
dev.off()

png(paste0(output_dir, 'humanHYPOMAP_extendedfig7_POMC.png'), width = 10, height = 10, units = 'in', res = 400)
huhy_plots(data = merged, gene = 'POMC', 
           percentile = 1, 
           image.alpha = 0, 
           assay.use = 'Spatial', 
           ncol = 3, 
           image = ordered_images, 
           pt.size.factor = 1.45, 
           crop = FALSE)  
dev.off()

png(paste0(output_dir, 'humanHYPOMAP_extendedfig7_AGRP.png'), width = 10, height = 10, units = 'in', res = 400)
huhy_plots(data = merged, gene = 'AGRP', 
           percentile = 1, 
           image.alpha = 0, 
           assay.use = 'Spatial', 
           ncol = 3, 
           image = ordered_images, 
           pt.size.factor = 1.45, 
           crop = FALSE) 
dev.off()
