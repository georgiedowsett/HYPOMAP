# Load libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

setwd('/path/to/humanHYPOMAP/')

# Load Lukas plotting functions
source("/path/to/scripts/utility_functions.R")
source("/path/to/scripts/plot_functions.R")

load_plot_params()

# Load data 
data <- readRDS('/path/to/humanHYPOMAP/data/nucseq-objects/human_hypo_combined.rds')

# Set output directory
output_dir <- 'figs_241004/' 

#######
# POMC/AGRP plot
#######
Idents(data) <- data@meta.data$C4
pomc <- WhichCells(data, idents = c('C4-373', 'C4-374', 'C4-375'))
agrp <- WhichCells(data, idents = c('C4-161', 'C4-293', 'C4-355'))
data$pomcagrp <- ifelse(colnames(data) %in% pomc, 'POMC',
                        ifelse(colnames(data) %in% agrp, 'AGRP', 'NO'))
cols <- c('POMC' = '#325A9BFF', 'AGRP' = '#B10DA1FF', 'NO' = 'lightgrey')

plot <- DimPlot(data, 
                group.by = 'pomcagrp', 
                reduction = 'umap_scvi_hypo', 
                raster = TRUE, 
                cols = cols,
                raster.dpi = c(rasterize_px, rasterize_px),
                pt.size = 2,
                order = TRUE) &
  NoAxes() &
  theme(legend.text = element_text(size = 20)) &
  ggtitle(element_blank())

ggsave(filename = paste0(output_dir, 'nucseq_plot_pomcagrp.png'),
       plot = plot, "png", dpi = 600, width = 350, height = 300, units = "mm")

####### 
# MC3R and MC4R feature plots for figure 5 
######

mcrs <- FeaturePlot(data, 
                    features = c('MC3R', 'MC4R'), 
                    cols = c('lightgrey', 'darkblue'), 
                    pt.size = 1, 
                    raster = FALSE,
                    order = TRUE,
                    label = FALSE,
                    reduction = 'umap_scvi_hypo') &
  NoAxes() & 
  theme(legend.text = element_text(size = 30), title = element_text(size = 40))

ggsave(filename = paste0(output_dir, 'Ext_fig_MCRs_feature.png'),
       plot = mcrs, "png", dpi = 600, width = 1200, height = 600, units = "mm")


mcrs <- FeaturePlot(data, 
                    features = c('GLP1R', 'GIPR'), 
                    cols = c('lightgrey', 'darkblue'), 
                    pt.size = 1, 
                    raster = FALSE,
                    order = TRUE,
                    label = FALSE,
                    reduction = 'umap_scvi_hypo') &
  NoAxes() & 
  theme(legend.text = element_text(size = 30), title = element_text(size = 40))

ggsave(filename = paste0(output_dir, 'Ext_fig_incretins_feature.png'),
       plot = mcrs, "png", dpi = 600, width = 1200, height = 600, units = "mm")
