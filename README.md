# Gene Expression Network Analysis

This repository contains R scripts for processing and analyzing gene expression data, along with network information. The main steps of the analysis include:

1. Loading and preprocessing data
2. Loading and preprocessing metadata
3. Creating a Seurat object and normalizing data
4. Selecting gene names based on network data
5. Filtering and processing expression set
6. Processing data and creating output files
7. Performing additional analysis (optional)

Here is a detailed explanation of each topic:

1. **Loading and preprocessing data:** The script begins by loading essential libraries and importing the gene expression data from the file ("data/GSE84465_GBM_All_data.csv") [(NCBI GEO, GSE84465)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84465). The data is converted into a dataframe, and unnecessary rows and columns are removed. Gene names are then assigned as row names to facilitate further analysis.

2. **Loading and preprocessing metadata:** The script loads and preprocesses the metadata by fetching the GEO series matrix file ("data/GSE84465_series_matrix.txt") [(NCBI GEO, GSE84465)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84465) using the `getGEO()` function from the GEOquery library. The metadata is then cleaned and reformatted by removing row names, setting column names, and extracting relevant information from the characteristics columns. This processed metadata is later used to create the expression set object and for filtering the expression data based on specific criteria.

3. **Creating a Seurat object and normalizing data:** In this step, the script utilizes the `CreateSeuratObject()` function from the Seurat library to create a Seurat object, which is a specialized data structure for storing and analyzing single-cell RNA-seq data. The Seurat object is created using the preprocessed gene expression data (`counts`) and the corresponding metadata (`p`). Next, the data is normalized using the `SCTransform()` function, which applies a variance-stabilizing transformation to the data [(article)](https://doi.org/10.1186/s13059-019-1874-1). This transformation helps to account for technical variations, and makes the data more suitable for downstream analyses like clustering and dimensionality reduction. The normalization process results in a Seurat object (`pbmc`) with normalized expression values, which is then used for further analysis. The normalized expression values are also extracted and saved in the `fullcounts` dataframe.

4. **Selecting gene names based on network data:** The script reads the network information from the "Network_v3.xls" file, which contains information about the interactions between genes (nodes) in the form of activation and inhibition. The network was constructed using MetaCore [(site)](https://clarivate.com/products/biopharma/research-development/early-research-intelligence-solutions/). The network data is processed to filter and extract relevant gene names (nodes) that are related to the transcription regulation mechanism. This step allows us to focus on the genes of interest in the context of the biological system being studied. After extracting the relevant gene names, the script then filters the gene expression data to include only the expression values corresponding to these selected genes, creating a new dataset (`newcounts`) for further analysis. By narrowing down the gene list based on network data, we can concentrate on the most important genes and their interactions.

5. **Filtering and processing expression set:** During this step, the script constructs an expression set object by combining the selected gene expression data (`newcounts`) and the associated metadata (`p`). The expression set is then subjected to filtering based on predefined criteria, such as cell type ("Neoplastic") and tissue type ("Tumor"), resulting in a refined set of data. In addition, genes that exhibit only null values across all samples are removed, retaining only those genes with non-zero expression levels. This filtered expression set (`eset`) is subsequently used in the following analysis steps.

6. **Processing data and creating output files:** In this stage, the script works with the filtered expression set (`eset`) and the network data (`nodes`) to generate adjacency matrices representing activation and inhibition interactions within the network [(see matrixfromtable repository)](https://github.com/marcosgvjunior/graph-matrix-and-combinatorics). These matrices are subsequently saved as "act_BT_ALL_seurat.xlsx" and "sup_BT_ALL_seurat.xlsx" files. In addition, the refined gene expression data, which has been filtered based on the network and specific criteria, is stored as "datapoints_seurat_BT_ALL.xlsx" for further analysis or visualization.

7. **Optional - Performing additional analysis:** The script offers an optional function, additional_analysis(), which is a preliminary implementation that offers a sample of possible further analysis. It is not optimized and should be improved before using it. It is divided into four main steps:

    1. Data filtering and preparation: Filtering samples based on is_filtered and is_filtered_genes, selecting genes of interest, and computing the distance matrix based on the scaled expression data.

    2. Clustering methods: Performing various clustering techniques, including DBSCAN, hierarchical clustering, K-means clustering, and Louvain clustering using Seurat library functions.

    3. Dimensionality reduction and visualization techniques: Applying Principal Component Analysis (PCA), t-Distributed Stochastic Neighbor Embedding (t-SNE), and Uniform Manifold Approximation and Projection (UMAP) for dimensionality reduction and visualization. The generated plots can be saved and utilized for more in-depth exploration.

    4. Differential expression analysis: Conducting differential expression analysis for two clusters using the MAST test.

These steps serve as a framework for the analysis of gene expression data, interconnecting it with a network and searching for patterns within the gene expression profiles.

## Usage

To use the script, run the R code provided in the file `gene_expression_network_analysis.R`. Ensure that the necessary input files are available in the specified paths, and the required R packages are installed.

## Input Files

1. Gene expression data: "data/GSE84465_GBM_All_data.csv" [(download page)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84465)
2. Annotation information: "data/GSE84465_series_matrix.txt" [(download page)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84465)
3. Network information: "Network_v3.xls" (not included)

## Output Files

1. Combined node information: "rede_BT_ALL_seurat.xlsx"
2. Activation adjacency matrix: "act_BT_ALL_seurat.xlsx"
3. Inhibition adjacency matrix: "sup_BT_ALL_seurat.xlsx"
4. Filtered gene expression data: "datapoints_seurat_BT_ALL.xlsx"

## Applications
These functions can be applied to the analysis of gene expression data and the construction of potential gene regulatory networks. This can aid in the understanding of cellular processes, disease mechanisms, and therapeutic targets.

## Dependencies
To use the functions in this repository, you will need the following R packages:
- Seurat
- dplyr
- ggplot2
- Matrix
- readxl
- openxlsx
- ... 

## Contributing

We welcome contributions to this repository. If you have any suggestions, improvements, or bug fixes, please feel free to submit a pull request or open an issue.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.