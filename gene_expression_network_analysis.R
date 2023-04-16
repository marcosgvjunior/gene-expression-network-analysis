#Instituto Oswaldo Cruz (PGBCS/IOC/Fiocruz)
#Ph.D. program on Computational and Systems Biology
#Student: Marcos Vieira
#September, 2022

# scRNA-Seq analysis

rm(list=ls())
gc()

# Load required libraries
library(stats4)
library(tidyverse)
library(data.table)
library(readxl)
library(writexl)
library(pryr)
library(GEOquery)
library(Seurat)
library(MAST)
library(dbscan)
library(dendextend)

source("matrixfromtable.R")

# Function to load and preprocess data
load_and_preprocess_data <- function() {

  xdata <- data.table::fread("data/GSE84465_GBM_All_data.csv", sep = " ")
  data  <- as.data.frame(xdata)
  
  mem_change(rm(xdata))  
  counts <- head(column_to_rownames(data, "V1"), -5)  
  mem_change(rm(data))
  
  return(counts)
}

# Function to load and preprocess metadata
load_and_preprocess_metadata <- function() {
  gse <- getGEO(GEO = NULL, filename = "data/GSE84465_series_matrix.txt", destdir = tempdir(),
                GSElimits = NULL, GSEMatrix = TRUE, AnnotGPL = FALSE, getGPL = FALSE,
                parseCharacteristics = TRUE)
    
  # Check dimensions of the metadata object
  cat("Dimensions of metadata object: ", dim(gse), "\n")

  p <- remove_rownames(pData(gse))
  
  mem_change(rm(gse))
  
  p <- column_to_rownames(p, "description.1")
  p <- p %>% select("title", "geo_accession","characteristics_ch1.1","characteristics_ch1.2","characteristics_ch1.3","characteristics_ch1.4",
                    "characteristics_ch1.5", "characteristics_ch1.6", "characteristics_ch1.7", "description" )
  p <- p %>% mutate( "characteristics_ch1.6" = gsub( "cell type: ", "", characteristics_ch1.6 ))
  p <- p %>% mutate( "characteristics_ch1.3" = gsub( "tissue: ", "", characteristics_ch1.3 ))
  p <- p %>% mutate( "characteristics_ch1.1" = gsub( "plate id: ", "", characteristics_ch1.1 ))
  p <- p %>% mutate( "characteristics_ch1.2" = gsub( "well: ", "", characteristics_ch1.2 ))
  p <- p %>% mutate( "characteristics_ch1.4" = gsub( "patient id: ", "", characteristics_ch1.4 ))
  p <- p %>% mutate( "characteristics_ch1.7" = gsub( "neoplastic: ", "", characteristics_ch1.7 ))
  p <- p %>% mutate( "characteristics_ch1.5" = gsub( "tsne cluster: ", "", characteristics_ch1.5 ))
  
  return(p)
}

# Function to create Seurat object and perform normalization
create_seurat_object_and_normalize <- function(counts, p) {

  seuset     <- CreateSeuratObject(as.matrix(counts), meta.data = as.data.frame(p))
  pbmc       <- SCTransform(seuset)
  fullcounts <- as.data.frame( GetAssayData(object = pbmc, slot = "data") )# ('counts', 'data', and 'scale.data')

  return(list(fullcounts = fullcounts, pbmc = pbmc))
}

selection_gene_names <- function(fullcounts){

  network               <- data.frame( read_excel( "Network_v3.xls" ) )
  nodes                 <- data.frame( network %>% select( nodeout, nodein, action, mechanism ) %>% filter( mechanism == "Transcription regulation" ) )
  names( nodes )        <- c( "nodeout", "nodein", "action", "mechanism" )
  nnames                <- dplyr::union( nodes[["nodeout"]], nodes[["nodein"]] )
  nnamessort            <- sort( nnames )
  totalnodes            <- length( nnames )
  newcounts             <- fullcounts[rownames(fullcounts) %in% nnamessort,]

  return(list(newcounts=newcounts, nodes=nodes))
}

# Function to filter and process expression set
filter_expression_set <- function(newcounts, p) {
  # Filter expression set object based on metadata
  eset_raw <- ExpressionSet( assayData = as.matrix(newcounts), phenoData = AnnotatedDataFrame(p)) 
  eset_raw <- eset_raw[ , eset_raw$characteristics_ch1.6 == c("Neoplastic")] 
  eset_raw <- eset_raw[ , eset_raw$characteristics_ch1.3 == c("Tumor")]
  
  # Check dimensions of the expression set object
  cat("Dimensions of expression set object: ", dim(eset_raw), "\n")
  
  # Print the row names
  cat("Row names of eset object: ", rownames(eset_raw), "\n")

  # Remove genes with only null values
  keep <- rowMeans(exprs(eset_raw)) > 0

  # Check number of genes to keep (non-zero row means)
  cat("Number of genes to keep: ", sum(keep), "\n")

  eset <- eset_raw[keep, ]

  # Check row names of the final expression set object
  cat("Sorted row names of expression set object: ", sort(rownames(eset)), "\n")
  
  return(eset)
}

# Function to process the data and create output files
process_data_and_create_output_files <- function(eset, nodes) {

  nodesfilter  <- nodes %>% filter( nodeout %in% rownames(eset) ) %>% filter( nodein %in% rownames(eset) )
  nodesfilterA <- nodesfilter %>% filter( action == "Activation" ) %>% mutate( action = 1 )
  nodesfilterI <- nodesfilter %>% filter( action == "Inhibition" ) %>% mutate( action = 2 )
  newnodes     <- bind_rows( nodesfilterA, nodesfilterI )

  write_xlsx( newnodes, "rede_BT_ALL_seurat.xlsx" )
  redeact   <- matrixfromtable( "rede_BT_ALL_seurat.xlsx", 1, "act_BT_ALL_seurat.xlsx" )
  redeinb   <- matrixfromtable( "rede_BT_ALL_seurat.xlsx", 2, "sup_BT_ALL_seurat.xlsx" )

  selecionarnos <- rownames(eset) %in% dplyr::union( newnodes[["nodeout"]], newnodes[["nodein"]] )
  datapoints    <- as.data.frame(eset) %>% rownames_to_column(var = "Cell")
  swha_seurat   <- as.data.frame(datapoints[,2:(length(selecionarnos)+1)][selecionarnos])
  swha_seurat   <- swha_seurat %>% select(sort(names(.)))
  swha_seurat   <- swha_seurat + abs(min(swha_seurat)) 

  # Check row names of the final expression set object
  cat("Column names of the final gene list: ", colnames(swha_seurat), "\n")
  
  write_xlsx( swha_seurat, "datapoints_seurat_BT_ALL.xlsx")
}

additional_analysis <- function(pbmc, is_filtered = FALSE,  is_filtered_genes = FALSE) {
  
  if (is_filtered) {
    # Filter data for Neoplastic and Tumor
    pbmc <- pbmc[, pbmc@meta.data$characteristics_ch1.6 == "Neoplastic" & pbmc@meta.data$characteristics_ch1.3 == "Tumor"]
  }

  if(is_filtered_genes){
    # Select the genes of interest
    selected_genes <- c("Gene1", "Gene2", "etc")
    pbmc <- pbmc[selected_genes, ]
  }

  # Compute distance matrix
  dist_matrix <- dist(t(pbmc@assays$SCT@scale.data), method = "euclidean")

  # Perform DBSCAN clustering
  dbscan_result <- dbscan::dbscan(dist_matrix, eps = 50, minPts = 3) # You may need to adjust the eps and minPts parameters
  pbmc@meta.data$DBSCAN_clusters <- as.factor(dbscan_result$cluster)
  
  # Perform hierarchical clustering
  hclust_result   <- hclust(dist_matrix, method = "complete") # methods: "single", "average", "complete", etc.
  # hclust_dendrogram <- as.dendrogram(hclust_result)
  # plot(hclust_dendrogram, main = "Hierarchical Clustering Dendrogram", xlab = "Samples", ylab = "Distance")
  cutree_clusters <- cutree(hclust_result, k = 7) # Choose the number of clusters after visualizing the dendrogram
  pbmc@meta.data$HClust_clusters <- as.factor(cutree_clusters)

  ### For the k-means clustering ###
  ### Set the range of possible cluster numbers
  cluster_range <- 3:12

  # Calculate the WCSS for each number of clusters
  WCSS_values <- sapply(cluster_range, function(k) {
  kmeans_result <- kmeans(t(pbmc@assays$SCT@scale.data), centers = k)
  kmeans_result$tot.withinss
  })

  # Plot the WCSS against the number of clusters
  elbow_plot_kmeans <- ggplot(data.frame(Clusters = cluster_range, WCSS = WCSS_values),
                            aes(x = Clusters, y = WCSS)) +
  geom_point() +
  geom_line() +
  xlab("Number of Clusters") +
  ylab("Within-Cluster Sum of Squares (WCSS)") +
  ggtitle("Elbow Plot for K-means Clustering")

  # Save the elbow plot
  ggsave("elbow_plot_kmeans.png", elbow_plot_kmeans, width = 6, height = 4, dpi = 300, bg = "white")
  ###
  ###

  # Perform k-means clustering
  kmeans_clusters <- kmeans(t(pbmc@assays$SCT@scale.data), centers = 10)
  pbmc@meta.data$KMeans_clusters <- as.factor(kmeans_clusters$cluster)

  # Perform PCA
  pbmc <- RunPCA(pbmc, npcs = 30, verbose = FALSE)
  elbow_plot <- ElbowPlot(pbmc, ndims = 30)
  ggsave("elbow_plot.png", elbow_plot, width = 6, height = 4, dpi = 300, bg="white")

  # Perform clustering
  pbmc <- FindNeighbors(pbmc, dims = 1:10)
  pbmc <- FindClusters(pbmc, resolution = 0.2)

  # Perform t-SNE
  pbmc <- RunTSNE(pbmc, dims = 1:10, perplexity = 80)
  
  tsne_plot_dbscan <- DimPlot(pbmc, group.by = "DBSCAN_clusters", label = TRUE, repel = TRUE, pt.size = 1)
  tsne_plot_hclust <- DimPlot(pbmc, group.by = "HClust_clusters", label = TRUE, repel = TRUE, pt.size = 1)
  tsne_plot_kmeans <- DimPlot(pbmc, group.by = "KMeans_clusters", label = TRUE, repel = TRUE, pt.size = 1)  
  tsne_plot_seurat <- DimPlot(pbmc, group.by = "seurat_clusters", label = TRUE, repel = TRUE, pt.size = 1)

  ggsave("tsne_plot_dbscan.png", tsne_plot_dbscan, width = 6, height = 4, dpi = 300)
  ggsave("tsne_plot_hclust.png", tsne_plot_hclust, width = 6, height = 4, dpi = 300) 
  ggsave("tsne_plot_kmeans.png", tsne_plot_kmeans, width = 6, height = 4, dpi = 300)    
  ggsave("tsne_plot_seurat.png", tsne_plot_seurat, width = 6, height = 4, dpi = 300)

  # Perform UMAP
  pbmc <- RunUMAP(pbmc, dims = 1:10)
  dim_plot <- DimPlot(pbmc, group.by = "characteristics_ch1.6", label = TRUE, repel = TRUE)
  ggsave("dim_plot.png", dim_plot, width = 6, height = 4, dpi = 300)

  # Perform differential expression analysis
  cluster1_markers <- FindMarkers(pbmc, ident.1 = 0, test.use = "MAST")
  cluster2_markers <- FindMarkers(pbmc, ident.1 = 1, test.use = "MAST")

  # Save output files
  write.table(cluster1_markers, file = "cluster1_markers.txt", sep = "\t", quote = FALSE, row.names = TRUE)
  write.table(cluster2_markers, file = "cluster2_markers.txt", sep = "\t", quote = FALSE, row.names = TRUE)

  return(pbmc)
}

# Main execution
counts <- load_and_preprocess_data()
p      <- load_and_preprocess_metadata()

# Check if row names of the metadata match the column names of the counts data
cat("Matching column names: ", all(rownames(p) == colnames(counts)), "\n")

fullcounts_pmbc <- create_seurat_object_and_normalize(counts, p)

newcounts_nodes <- selection_gene_names(fullcounts_pmbc$fullcounts)

# Filter the expression set
eset <- filter_expression_set(newcounts_nodes$newcounts, p)

# Call the process_data_and_create_output_files function
process_data_and_create_output_files(eset, newcounts_nodes$nodes)

# Optional - Performing additional analysis
# additional_analysis(fullcounts_pmbc$pbmc, FALSE, FALSE)

