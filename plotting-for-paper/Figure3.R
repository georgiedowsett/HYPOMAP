# Load packages
library(Seurat)
library(ggplot2)
library(patchwork)
library(viridis)
library(viridisLite)
library(scales)
library(grid)
source("/path/to/scripts/utility_functions.R")
source("/path/to/scripts/plot_functions.R")
setwd('/path/to/humanHYPOMAP/')
# Load data
merged <- readRDS('data/humanHYPOMAP_spatial_c2l_240621.RDS')

load_plot_params()
load_colors()
getLongPalette = colorRampPalette(long_palette_strong)
getOkabeItoPalette = colorRampPalette(short_palette)
output_dir <- 'figs_241004/'

# Set image order
ordered_images <- c('slice5A', 'sliceC1A', 'slice6B', 'sliceB1B', 'slice2B', 'slice4B', 'slice7A', 'slice3A', 'slice8B')

# Figure 3B
my_spatial_palette <- colorRampPalette(c(viridis(7, begin = 0.2)), bias = 2)(n = 256)
plot_list <- list()
for (i in c('SLC17A6', 'SLC32A1', 'TBX3', 'FEZF1', 'SIM1')) {
  plot <- SpatialFeaturePlot(merged, 
                            features = i, 
                            images = 'slice2B', 
                            crop = F, 
                            image.alpha = 0, 
                            pt.size.factor = 1.3,
                            stroke = NA) & 
    scale_fill_gradientn(colors = my_spatial_palette, 
                         oob=squish) &
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
plots <- wrap_plots(plot_list, ncol = 5)

pdf(paste0(output_dir, 'humanHYPOMAP_figure3B.pdf'), width = 25, height = 5)
plots
dev.off()

png(paste0(output_dir, 'humanHYPOMAP_figure3B.png'), width = 25, height = 5, units = 'in', res = 400)
plots
dev.off()


# Figure 3G spatial plots colored by regions 
# Set color scheme
col <- getOkabeItoPalette(23)
names(col) <- unique(merged$NN27_res_0.5_grouped)
col <- data.frame(col)
write.table(col, 'tables/spatial_clusters_colours_OkabeItoPalette.txt', sep = '\t', row.names = T)
cols <- as.list(col$col)
names(cols) <- rownames(col)

pdf(paste0(output_dir, 'humanHYPOMAP_figure3G_nolegend.pdf'), width = 35, height = 5)
SpatialDimPlot(merged,
               group.by = 'NN27_res_0.5_grouped', 
               ncol = 9, 
               crop = F, 
               image.alpha = 0, 
               images = ordered_images, 
               pt.size.factor = 1.45, 
               cols = cols,
               stroke = NA) & 
  theme(legend.title = element_blank(), plot.title = element_blank()) & NoLegend() 
dev.off()

png(paste0(output_dir, 'humanHYPOMAP_figure3G_nolegend.png'), width = 35, height = 5, units = 'in', res = 400)
SpatialDimPlot(merged,
               group.by = 'NN27_res_0.5_grouped', 
               ncol = 9, 
               crop = F, 
               image.alpha = 0, 
               images = ordered_images, 
               pt.size.factor = 1.45, 
               cols = cols,
               stroke = NA) & 
  theme(legend.title = element_blank(), plot.title = element_blank()) & NoLegend() 
dev.off()

# Printing the legends
legend <- cowplot::get_legend(DimPlot(merged, group.by = 'NN27_res_0.5_grouped', cols = cols))
# Create a blank plot with just the legend
legend_plot <- ggplot() + 
  theme_void() + 
  annotation_custom(legend)
# Display the legend plot
grid.newpage()

pdf(paste0(output_dir, 'humanHYPOMAP_figure3G_legend.pdf'))
grid.draw(legend)
dev.off()

rm(merged)
gc()

#### UMAP plot of the NucSeq data colored by regional assignment
human_hypo_combined <- readRDS('data/nucseq-objects/human_hypo_combined.rds')

umap <- DimPlot(human_hypo_combined,
        group.by = "region_threshold_manual",
        raster = FALSE,
        cols = cols,
        reduction = "umap_scvi_hypo",
        label = F,
        shuffle = TRUE,
        raster.dpi = c(1536,1536)) + 
  NoAxes() +
  ggtitle(element_blank())

ggsave(filename = paste0(output_dir, 'humanHYPOMAP_figure3f.png'),
       plot = umap, "png", dpi = 600, width = 350, height = 300, units = "mm")
