
# HYPOMAP 

This repository contains the analysis scripts for the spatial transcriptomics aspects of the human HYPOMAP dataset and paper. 


## Installation

The following docker image was used to perform the cell2location pipeline 
```bash
docker pull brianlamx/cell2location
```
The following docker image used to perform analysis using R
```bash 
docker pull brianlamx/seurat_rstudio
```
    
## Documentation

Each analysis script found in [Cell2lLocation_and_ST_analysis](https://github.com/georgiedowsett/HYPOMAP/tree/main/Cell2Location_and%20_ST_analysis) performs the following: 

[1_checkNucSeqObjectforC2L.py](https://github.com/georgiedowsett/HYPOMAP/blob/main/Cell2Location_and%20_ST_analysis/1_checkNucSeqObjectforC2L.py) looks at the snRNAseq h5ad object to see if it is suitable for cell2location analysis pipeline 

[2_cell2location_nucseq_model_refsigs.py](https://github.com/georgiedowsett/HYPOMAP/blob/main/Cell2Location_and%20_ST_analysis/2_cell2location_nucseq_model_refsigs.py) Uses the Cell2location pipeline (part 1) to create reference signatures for each cluster at each cluster level in the snRNAseq object 

 3-6 Uses the cell2location pipeline to estimate snRNAseq cluster abundance at each spot in the cell2location object for each level of clustering in the mrtree. Leiden clustering is performed to create regional clustering in the spatial datasets

 [7_Spatial-data-creation.r](https://github.com/georgiedowsett/HYPOMAP/blob/main/Cell2Location_and%20_ST_analysis/7_Spatial-data-creation.r) Creates a spatial transcriptomics seurat object containing normalised count data, and the output from all 4 cell2location runs. It adds the leiden clustering, and labels the clusters based on regions. 

 [8_average_abundance_spatial_tables.r](https://github.com/georgiedowsett/HYPOMAP/blob/main/Cell2Location_and%20_ST_analysis/8_average_abundance_spatial_tables.r) Creates mean and median cell abundance scores for each snRNAseq cluster in each regional cluster (from the leiden clustering in 5)

 [9_Regional_assignment_thresholding.r](https://github.com/georgiedowsett/HYPOMAP/blob/main/Cell2Location_and%20_ST_analysis/9_Regional_assignment_thresholding.r) Assigns a region to each snRNAseq cluster at the C3 level. 



 


## Used By

This analysis is used in the following bioRxiv paper: [Human HYPOMAP: A comprehensive spatio-cellular map of the human hypothalamus](https://www.biorxiv.org/content/10.1101/2023.09.15.557967v1)


## Authors

- [@georgiedowsett](https://github.com/georgiedowsett)


## License

For the purpose of open access, the author has applied a Creative Commons Attribution (CC BY) licence to any Author Accepted Manuscript version arising from this submission

