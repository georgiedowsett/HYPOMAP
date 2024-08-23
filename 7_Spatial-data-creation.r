##### This script will add the cell2location output data to the spatial R object #####

# First, we load packages 
library(Seurat)
library(ggplot2)
library(tidyverse)
library(scCustomize)

# Define paths as variables
base_path <- "~/path/to/your/base/directory/"
barcode_path <- "data/barcodes-to-remove/"
cell2loc_path <- "data/c2l-outs/240424/"
output_path <- "data/humanHYPOMAP_spatial_c2l_240510.RDS"
leiden_clustering_path <- "tables/leiden_regional_clustering.csv"

# Load raw data
b2 <- Load10X_Spatial(paste0(base_path, "new-spatial/2B/"), filename = "filtered_feature_bc_matrix.h5", slice = "slice2B")
a3 <- Load10X_Spatial(paste0(base_path, "new-spatial/3A/"), filename = "filtered_feature_bc_matrix.h5", slice = "slice3A")
b4 <- Load10X_Spatial(paste0(base_path, "new-spatial/4B/"), filename = "filtered_feature_bc_matrix.h5", slice = "slice4B")
a5 <- Load10X_Spatial(paste0(base_path, "new-spatial/5A/"), filename = "filtered_feature_bc_matrix.h5", slice = "slice5A")
b6 <- Load10X_Spatial(paste0(base_path, "new-spatial/6B/"), filename = "filtered_feature_bc_matrix.h5", slice = "slice6B")
a7 <- Load10X_Spatial(paste0(base_path, "new-spatial/7A/"), filename = "filtered_feature_bc_matrix.h5", slice = "slice7A")
b8 <- Load10X_Spatial(paste0(base_path, "new-spatial/8B/"), filename = "filtered_feature_bc_matrix.h5", slice = "slice8B")
b1_b <- Load10X_Spatial(paste0(base_path, "new-spatial/B1_B/"), filename = "filtered_feature_bc_matrix.h5", slice = "sliceB1B")
c1_a <- Load10X_Spatial(paste0(base_path, "new-spatial/C1_A/"), filename = "filtered_feature_bc_matrix.h5", slice = "sliceC1A")

# Remove unwanted barcodes
barcodes <- list.files(path = barcode_path, pattern = ".csv")
remove <- list()
for (i in barcodes) {
  remove[[i]] <- read.csv(paste0(barcode_path, i))
}
remove <- remove[-3]

slides <- c(b1_b, c1_a, b2, a3, b4, a5, b6, a7, b8)
slidenames <- c('b1_b', 'c1_a', 'b2', 'a3', 'b4', 'a5', 'b6', 'a7', 'b8')
for (i in 1:length(slides)) {
  print(dim(slides[[i]]))
  keep <- setdiff(colnames(slides[[i]]), remove[[i]]$Barcode)
  slides[[i]] <- subset(slides[[i]], cells = keep)
  print(dim(slides[[i]]))
}

# Add metadata
# Experiment location
for (i in 1:9) {
  slides[[i]] <- AddMetaData(slides[[i]], metadata = slidenames[i], col.name = "captureArea")
}
for (i in 1:2) {
  slides[[i]] <- AddMetaData(slides[[i]], metadata = "Cambridge", col.name = "Experiment")
}
for (i in 3:9) {
  slides[[i]] <- AddMetaData(slides[[i]], metadata = "Copenhagen", col.name = "Experiment")
}
# Slide / batch
for (i in 1:2) {
  slides[[i]] <- AddMetaData(slides[[i]], metadata = "slide5", col.name = "slideNumber")
}
slides[[3]] <- AddMetaData(slides[[3]], metadata = "slide1", col.name = "slideNumber")
for (i in 4:5) {
  slides[[i]] <- AddMetaData(slides[[i]], metadata = "slide2", col.name = "slideNumber")
}
for (i in 6:7) {
  slides[[i]] <- AddMetaData(slides[[i]], metadata = "slide3", col.name = "slideNumber")
}
for (i in 8:9) {
  slides[[i]] <- AddMetaData(slides[[i]], metadata = "slide4", col.name = "slideNumber")
}

# Add other metadata: donorID, Sex, Age, BMI
slides[[1]] <- AddMetaData(slides[[1]], metadata = "IbV8d", col.name = "donorID")
slides[[2]] <- AddMetaData(slides[[2]], metadata = "q65KC", col.name = "donorID")
for (i in c(3, 5)) {
  slides[[i]] <- AddMetaData(slides[[i]], metadata = "oYtkk", col.name = "donorID")
}
slides[[4]] <- AddMetaData(slides[[4]], metadata = "hp9Yr", col.name = "donorID")
slides[[7]] <- AddMetaData(slides[[7]], metadata = "OgTiF", col.name = "donorID")
for (i in c(6, 8)) {
  slides[[i]] <- AddMetaData(slides[[i]], metadata = "JpLuT", col.name = "donorID")
}
slides[[9]] <- AddMetaData(slides[[9]], metadata = "bT5r4", col.name = "donorID")

# Add Sex
for (i in c(1:3, 5, 6, 8, 9)) {
  slides[[i]] <- AddMetaData(slides[[i]], metadata = "Female", col.name = "Sex")
}
for (i in c(4, 7)) {
  slides[[i]] <- AddMetaData(slides[[i]], metadata = "Male", col.name = "Sex")
}

# Add age 
slides[[1]] <- AddMetaData(slides[[1]], metadata = "95", col.name = "Age")
slides[[2]] <- AddMetaData(slides[[2]], metadata = "86", col.name = "Age")
for (i in c(3, 5)) {
  slides[[i]] <- AddMetaData(slides[[i]], metadata = "79", col.name = "Age")
}
slides[[4]] <- AddMetaData(slides[[4]], metadata = "75", col.name = "Age")
slides[[7]] <- AddMetaData(slides[[7]], metadata = "73", col.name = "Age")
for (i in c(6, 8)) {
  slides[[i]] <- AddMetaData(slides[[i]], metadata = "51", col.name = "Age")
}
slides[[9]] <- AddMetaData(slides[[9]], metadata = "54", col.name = "Age")

# Add BMI
for (i in c(1, 2, 7)) {
  slides[[i]] <- AddMetaData(slides[[i]], metadata = "unknown", col.name = "BMI")
}
for (i in c(3, 5)) {
  slides[[i]] <- AddMetaData(slides[[i]], metadata = "29.5", col.name = "BMI")
}
slides[[4]] <- AddMetaData(slides[[4]], metadata = "16.6", col.name = "BMI")
for (i in c(6, 8)) {
  slides[[i]] <- AddMetaData(slides[[i]], metadata = "16.5", col.name = "BMI")
}
slides[[9]] <- AddMetaData(slides[[9]], metadata = "41", col.name = "BMI")

##### Merge and Normalize data ######
merged <- merge(slides[[3]], c(slides[[4]], slides[[5]], slides[[6]], slides[[7]], slides[[8]], slides[[9]], slides[[1]], slides[[2]]))
# Clear workspace 
rm(slides)
rm(b2, b4, b6, b8, a3, a5, a7, b1_b, c1_a)
gc()

merged <- NormalizeData(merged)

###### Adding cell2location data as separate assays #######
cell2loc <- read.delim(paste0(cell2loc_path, "C4/C4_clusters_obs_table.csv"), sep = ",")

table(substring(colnames(merged), 18, last = 20))
cell2loc$barcodes <- substr(cell2loc$spot_id, nchar(cell2loc$spot_id) - 17, nchar(cell2loc$spot_id))

cell2loc$slide <- NA
cell2loc$slide[substr(cell2loc$spot_id, 1, 3) == "2B_"] <- "_1"
cell2loc$slide[substr(cell2loc$spot_id, 1, 3) == "3A_"] <- "_2"
cell2loc$slide[substr(cell2loc$spot_id, 1, 3) == "4B_"] <- "_3"
cell2loc$slide[substr(cell2loc$spot_id, 1, 3) == "5A_"] <- "_4"
cell2loc$slide[substr(cell2loc$spot_id, 1, 3) == "6B_"] <- "_5"
cell2loc$slide[substr(cell2loc$spot_id, 1, 3) == "7A_"] <- "_6"
cell2loc$slide[substr(cell2loc$spot_id, 1, 3) == "8B_"] <- "_7"
cell2loc$slide[substr(cell2loc$spot_id, 1, 5) == "B1_B_"] <- "_8"
cell2loc$slide[substr(cell2loc$spot_id, 1, 5) == "C1_A_"] <- "_9"

cell2loc$spot_ID <- paste0(cell2loc$barcodes, cell2loc$slide)

cell2loc <- subset(cell2loc, select = -c(barcodes, slide, spot_id))
rownames(cell2loc) <- cell2loc$spot_ID
cell2loc <- subset(cell2loc, select = -c(spot_ID))

# Remove other metadata
cell2loc <- cell2loc[, 22:473]

# Add as a new assay 
c4 <- t(cell2loc)
merged[["C4"]] <- CreateAssayObject(data = c4)

# c3
cell2loc <- read.delim(paste0(cell2loc_path, "C3/C3_clusters_obs_table.csv"), sep = ",")

table(substring(colnames(merged), 18, last = 20))
cell2loc$barcodes <- substr(cell2loc$spot_id, nchar(cell2loc$spot_id) - 17, nchar(cell2loc$spot_id))

cell2loc$slide <- NA
cell2loc$slide[substr(cell2loc$spot_id, 1, 3) == "2B_"] <- "_1"
cell2loc$slide[substr(cell2loc$spot_id, 1, 3) == "3A_"] <- "_2"
cell2loc$slide[substr(cell2loc$spot_id, 1, 3) == "4B_"] <- "_3"
cell2loc$slide[substr(cell2loc$spot_id, 1, 3) == "5A_"] <- "_4"
cell2loc$slide[substr(cell2loc$spot_id, 1, 3) == "6B_"] <- "_5"
cell2loc$slide[substr(cell2loc$spot_id, 1, 3) == "7A_"] <- "_6"
cell2loc$slide[substr(cell2loc$spot_id, 1, 3) == "8B_"] <- "_7"
cell2loc$slide[substr(cell2loc$spot_id, 1, 5) == "B1_B_"] <- "_8"
cell2loc$slide[substr(cell2loc$spot_id, 1, 5) == "C1_A_"] <- "_9"

cell2loc$spot_ID <- paste0(cell2loc$barcodes, cell2loc$slide)

cell2loc <- subset(cell2loc, select = -c(barcodes, slide, spot_id))
rownames(cell2loc) <- cell2loc$spot_ID
cell2loc <- subset(cell2loc, select = -c(spot_ID))

# Remove other metadata
cell2loc <- cell2loc[, 22:ncol(cell2loc)]

# Add as a new assay 
c3 <- t(cell2loc)
merged[["C3"]] <- CreateAssayObject(data = c3)

# c2
cell2loc <- read.delim(paste0(cell2loc_path, "C2/C2_clusters_obs_table.csv"), sep = ",")

table(substring(colnames(merged), 18, last = 20))
cell2loc$barcodes <- substr(cell2loc$spot_id, nchar(cell2loc$spot_id) - 17, nchar(cell2loc$spot_id))

cell2loc$slide <- NA
cell2loc$slide[substr(cell2loc$spot_id, 1, 3) == "2B_"] <- "_1"
cell2loc$slide[substr(cell2loc$spot_id, 1, 3) == "3A_"] <- "_2"
cell2loc$slide[substr(cell2loc$spot_id, 1, 3) == "4B_"] <- "_3"
cell2loc$slide[substr(cell2loc$spot_id, 1, 3) == "5A_"] <- "_4"
cell2loc$slide[substr(cell2loc$spot_id, 1, 3) == "6B_"] <- "_5"
cell2loc$slide[substr(cell2loc$spot_id, 1, 3) == "7A_"] <- "_6"
cell2loc$slide[substr(cell2loc$spot_id, 1, 3) == "8B_"] <- "_7"
cell2loc$slide[substr(cell2loc$spot_id, 1, 5) == "B1_B_"] <- "_8"
cell2loc$slide[substr(cell2loc$spot_id, 1, 5) == "C1_A_"] <- "_9"

cell2loc$spot_ID <- paste0(cell2loc$barcodes, cell2loc$slide)

cell2loc <- subset(cell2loc, select = -c(barcodes, slide, spot_id))
rownames(cell2loc) <- cell2loc$spot_ID
cell2loc <- subset(cell2loc, select = -c(spot_ID))

# Remove other metadata
cell2loc <- cell2loc[, 22:ncol(cell2loc)]

# Add as a new assay 
c2 <- t(cell2loc)
merged[["C2"]] <- CreateAssayObject(data = c2)

# c1
cell2loc <- read.delim(paste0(cell2loc_path, "C1/C1_clusters_obs_table.csv"), sep = ",")

table(substring(colnames(merged), 18, last = 20))
cell2loc$barcodes <- substr(cell2loc$spot_id, nchar(cell2loc$spot_id) - 17, nchar(cell2loc$spot_id))

cell2loc$slide <- NA
cell2loc$slide[substr(cell2loc$spot_id, 1, 3) == "2B_"] <- "_1"
cell2loc$slide[substr(cell2loc$spot_id, 1, 3) == "3A_"] <- "_2"
cell2loc$slide[substr(cell2loc$spot_id, 1, 3) == "4B_"] <- "_3"
cell2loc$slide[substr(cell2loc$spot_id, 1, 3) == "5A_"] <- "_4"
cell2loc$slide[substr(cell2loc$spot_id, 1, 3) == "6B_"] <- "_5"
cell2loc$slide[substr(cell2loc$spot_id, 1, 3) == "7A_"] <- "_6"
cell2loc$slide[substr(cell2loc$spot_id, 1, 3) == "8B_"] <- "_7"
cell2loc$slide[substr(cell2loc$spot_id, 1, 5) == "B1_B_"] <- "_8"
cell2loc$slide[substr(cell2loc$spot_id, 1, 5) == "C1_A_"] <- "_9"

cell2loc$spot_ID <- paste0(cell2loc$barcodes, cell2loc$slide)

cell2loc <- subset(cell2loc, select = -c(barcodes, slide, spot_id))
rownames(cell2loc) <- cell2loc$spot_ID
cell2loc <- subset(cell2loc, select = -c(spot_ID))

# Remove other metadata
cell2loc <- cell2loc[, 22:ncol(cell2loc)]

# Add as a new assay 
c1 <- t(cell2loc)
merged[["C1"]] <- CreateAssayObject(data = c1)

#functions for rotating seurat objects taken from: https://github.com/satijalab/seurat/issues/2702
rotimat=function(foo,rotation){
  if(!is.matrix(foo)){
    cat("Input is not a matrix")
    return(foo)
  }
  if(!(rotation %in% c("180","Hf","Vf", "R90", "L90"))){
    cat("Rotation should be either L90, R90, 180, Hf or Vf\n")
    return(foo)
  }
  if(rotation == "180"){
    foo <- foo %>% 
      .[, dim(.)[2]:1] %>%
      .[dim(.)[1]:1, ]
  }
  if(rotation == "Hf"){
    foo <- foo %>%
      .[, dim(.)[2]:1]
  }
  
  if(rotation == "Vf"){
    foo <- foo %>%
      .[dim(.)[1]:1, ]
  }
  if(rotation == "L90"){
    foo = t(foo)
    foo <- foo %>%
      .[dim(.)[1]:1, ]
  }
  if(rotation == "R90"){
    foo = t(foo)
    foo <- foo %>%
      .[, dim(.)[2]:1]
  }
  return(foo)
}

rotateSeuratImage = function(seuratVisumObject, slide = "slice1", rotation="Vf"){
  if(!(rotation %in% c("180","Hf","Vf", "L90", "R90"))){
    cat("Rotation should be either 180, L90, R90, Hf or Vf\n")
    return(NULL)
  }else{
    seurat.visium = seuratVisumObject
    ori.array = (seurat.visium@images)[[slide]]@image
    img.dim = dim(ori.array)[1:2]/(seurat.visium@images)[[slide]]@scale.factors$lowres
    new.mx <- c()  
    # transform the image array
    for (rgb_idx in 1:3){
      each.mx <- ori.array[,,rgb_idx]
      each.mx.trans <- rotimat(each.mx, rotation)
      new.mx <- c(new.mx, list(each.mx.trans))
    }
    
    # construct new rgb image array
    new.X.dim <- dim(each.mx.trans)[1]
    new.Y.dim <- dim(each.mx.trans)[2]
    new.array <- array(c(new.mx[[1]],
                         new.mx[[2]],
                         new.mx[[3]]), 
                       dim = c(new.X.dim, new.Y.dim, 3))
    
    #swap old image with new image
    seurat.visium@images[[slide]]@image <- new.array
    
    ## step4: change the tissue pixel-spot index
    img.index <- (seurat.visium@images)[[slide]]@coordinates
    
    #swap index
    if(rotation == "Hf"){
      seurat.visium@images[[slide]]@coordinates$imagecol <- img.dim[2]-img.index$imagecol
    }
    
    if(rotation == "Vf"){
      seurat.visium@images[[slide]]@coordinates$imagerow <- img.dim[1]-img.index$imagerow
    }
    
    if(rotation == "180"){
      seurat.visium@images[[slide]]@coordinates$imagerow <- img.dim[1]-img.index$imagerow
      seurat.visium@images[[slide]]@coordinates$imagecol <- img.dim[2]-img.index$imagecol
    }
    
    if(rotation == "L90"){
      seurat.visium@images[[slide]]@coordinates$imagerow <- img.dim[2]-img.index$imagecol
      seurat.visium@images[[slide]]@coordinates$imagecol <- img.index$imagerow
    }
    
    if(rotation == "R90"){
      seurat.visium@images[[slide]]@coordinates$imagerow <- img.index$imagecol
      seurat.visium@images[[slide]]@coordinates$imagecol <- img.dim[1]-img.index$imagerow
    }
    
    return(seurat.visium)
  }  
}

#rotate images and adding C2L clusters 
merged<-rotateSeuratImage(merged, slide = 'slice2B', rotation = 'R90')
merged<-rotateSeuratImage(merged, slide = 'slice3A', rotation = 'L90')
merged<-rotateSeuratImage(merged, slide = 'slice4B', rotation = 'L90')
merged<-rotateSeuratImage(merged, slide = 'slice5A', rotation = 'L90')
merged<-rotateSeuratImage(merged, slide = 'slice6B', rotation = 'L90')
merged<-rotateSeuratImage(merged, slide = 'slice7A', rotation = 'L90')
merged<-rotateSeuratImage(merged, slide = 'slice8B', rotation = 'R90')
merged<-rotateSeuratImage(merged, slide = 'sliceB1B', rotation = 'R90')
merged<-rotateSeuratImage(merged, slide = 'sliceC1A', rotation = '180') 

# Load leiden clustering information
cell2loc <- read.delim(leiden_clustering_path, sep = ",")

cell2loc<-cell2loc[,c(1, 175:179)]

table(substring(colnames(merged), 18, last = 20))
cell2loc$barcodes <- substr(cell2loc$spot_id, nchar(cell2loc$spot_id) - 17, nchar(cell2loc$spot_id))

cell2loc$slide <- NA  
cell2loc$slide[substr(cell2loc$spot_id, 1, 3) == "2B_"] <- "_1"
cell2loc$slide[substr(cell2loc$spot_id, 1, 3) == "3A_"] <- "_2"
cell2loc$slide[substr(cell2loc$spot_id, 1, 3) == "4B_"] <- "_3"
cell2loc$slide[substr(cell2loc$spot_id, 1, 3) == "5A_"] <- "_4"
cell2loc$slide[substr(cell2loc$spot_id, 1, 3) == "6B_"] <- "_5"
cell2loc$slide[substr(cell2loc$spot_id, 1, 3) == "7A_"] <- "_6"
cell2loc$slide[substr(cell2loc$spot_id, 1, 3) == "8B_"] <- "_7"
cell2loc$slide[substr(cell2loc$spot_id, 1, 5) == "B1_B_"] <- "_8"
cell2loc$slide[substr(cell2loc$spot_id, 1, 5) == "C1_A_"] <- "_9"

cell2loc$spot_ID<-paste0(cell2loc$barcodes, cell2loc$slide)
clusters<-cell2loc[,c(2:6, 9)]

for (i in colnames(clusters[,1:5])){
  leiden_vector<-setNames(clusters[[i]], clusters$spot_ID)
  merged<-AddMetaData(merged, metadata = leiden_vector, col.name = print(i))
}

#Creating order of slides 
ordered_images<-c('slice5A', 'sliceC1A', 'slice6B', 'sliceB1B', 'slice2B', 'slice4B', 'slice7A', 'slice3A', 'slice8B' )

rm(remove, c1, c2, c3, c4)
gc()

#choosing NN27 res 0.5 as the best clustering 
Idents(merged)<-'region_cluster_NN_27_res_0.5'

markers<-FindAllMarkers(merged, only.pos = T, logfc.threshold = 0.4)

#Look at the C3 clusters markers 
DefaultAssay(merged)<-'C3'
c3_markers<-FindAllMarkers(merged, only.pos = T)


#Adding names to NN12 res 0.5 clusters 
NN27_names<- c('0' = 'LH/RCN',
              '1' = 'Fx/LH',
              '2' = 'Unassigned',
              '3' = 'LH 1',
              '4' = 'VMH',
              '5' = 'OT',
              '6' = 'Fx/OT/ac',
              '7' = 'PVN',
              '8' = 'ARC',
              '9' = 'TMN',
              '10' = 'LPOA',
              '11' = 'Vascular',
              '12' = 'Ventricular',
              '13' = 'MPOA',
              '14' = 'SON',
              '15' = 'LH 2',
              '16' = 'Periventricular',
              '17' = 'POA',
              '18' = 'Thalamus 1',
              '19' = 'DMH',
              '20' = 'LTN',
              '21' = 'ME',
              '22' = 'SCN',
              '23' = 'LH 3',
              '24' = 'Thalamus 2',
              '25' = 'MAM 1',
              '26' = 'MAM 2')
merged$NN27_res_0.5_named <- NN27_names[as.character(merged@meta.data$region_cluster_NN_27_res_0.5)]

#Adding names to NN12 res 0.5 clusters 
NN27_grouped<- c('0' = 'LH/RCN',
                 '1' = 'Fx/LH',
                 '2' = 'Unassigned',
                 '3' = 'LH',
                 '4' = 'VMH',
                 '5' = 'OT',
                 '6' = 'Fx/OT/ac',
                 '7' = 'PVN',
                 '8' = 'ARC',
                 '9' = 'TMN',
                 '10' = 'LPOA',
                 '11' = 'Vascular',
                 '12' = 'Ventricular',
                 '13' = 'MPOA',
                 '14' = 'SON',
                 '15' = 'LH',
                 '16' = 'Periventricular',
                 '17' = 'POA',
                 '18' = 'Thalamus',
                 '19' = 'DMH',
                 '20' = 'LTN',
                 '21' = 'ME',
                 '22' = 'SCN',
                 '23' = 'LH',
                 '24' = 'Thalamus',
                 '25' = 'MAM',
                 '26' = 'MAM')
merged$NN27_res_0.5_grouped <- NN27_grouped[as.character(merged@meta.data$region_cluster_NN_27_res_0.5)]

#Creating PCA data for figure legends
DefaultAssay(merged)<-'Spatial'
merged<-FindVariableFeatures(merged)
merged<-ScaleData(merged)
merged<-RunPCA(merged)

saveRDS(merged, output_path)

