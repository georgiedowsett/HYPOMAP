library(Seurat)
library(ggplot2)
library(cowplot)
library(viridis)
library(viridisLite)
library(scales)
library(patchwork)

setwd('/path/to/humanHYPOMAP/')
# Load data
merged <- readRDS('data/humanHYPOMAP_spatial_c2l_240621.RDS')

# Set image order
ordered_images <- c('slice5A', 'sliceC1A', 'slice6B', 'sliceB1B', 'slice2B', 'slice4B', 'slice7A', 'slice3A', 'slice8B')

# Set output directory
output_dir <- 'figs_241004/'

my_spatial_palette <- colorRampPalette(c(viridis(7, begin = 0.2)), bias = 2)(n = 256)
DefaultAssay(merged) <- 'Spatial'

# Transcription factor plotting
png(paste0(output_dir, 'humanHYPOMAP_extfig3_spatial.png'), width = 35, height = 40, res = 400, units = 'in')
wrap_plots(SpatialFeaturePlot(merged, 
                              c('MEIS2', 'FEZF1', 'SIX3', 'SIM1', 'LHX6', 'TBX3', 'OTP', 'FOXB1'),
                              image = ordered_images, 
                              image.alpha = 0, 
                              pt.size.factor = 1.5, 
                              crop = FALSE, 
                              ncol = 9, 
                              stroke = FALSE) & 
             scale_fill_gradientn(limits = c(0, 3), 
                                  breaks = c(0,1,2,3), 
                                  colors = my_spatial_palette, 
                                  oob = squish) &
             theme(plot.title = element_blank())) +
  plot_layout(guide = 'collect') &
  theme(legend.position = 'right', 
        legend.key.size = unit(2, 'cm'), 
        legend.text = element_text(size = 40), 
        legend.title = element_blank())
dev.off()
