###### This script is for checking the snRNAseq h5ad object is appropriately formatted ######

#import packages
import anndata as ad
import scanpy as sc
import pandas as pd
from typing import Optional

# Define file paths as variables
snRNAseq_file_path: str = "/path/to/human_hypo_combined.h5ad"
output_file_path: str = "/path/to/output.txt"

def main() -> None:
    """
    Main function to check the formatting of the snRNAseq h5ad object and write metadata and data samples to a text file.
    """
    # Open a text file for writing
    with open(output_file_path, "w") as file:
        # Load snRNAseq file
        adata = sc.read_h5ad(snRNAseq_file_path)

        # Write metadata columns to file
        file.write("Metadata columns:\n")
        file.write(adata.obs.head(10).to_string())
        file.write("\n\n")

        # Look at first 10 rows and columns of raw data
        if adata.raw:
            file.write("Head of raw data (dense format with labels):\n")
            raw_data_dense = adata.raw[:, :].X[:10, :10].toarray()
            raw_data_df = pd.DataFrame(raw_data_dense, index=adata.obs_names[:10], columns=adata.var_names[:10])
            file.write(raw_data_df.to_string())
            file.write("\n\n")

        # Look at normalized/main data (first 10x10 elements)
        file.write("Head of normalized/main data (dense format with labels):\n")
        main_data_dense = adata.X[:10, :10].toarray()  # Convert sparse matrix to dense
        main_data_df = pd.DataFrame(main_data_dense, index=adata.obs_names[:10], columns=adata.var_names[:10])
        file.write(main_data_df.to_string())
        file.write("\n\n")

        # Look at genes (variable names)
        file.write("Variable names (genes/features):\n")
        file.write(str(adata.var_names[:10]))
        file.write("\n\n")

if __name__ == "__main__":
    main()
