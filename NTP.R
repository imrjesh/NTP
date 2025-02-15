### Load Libraries ####
library(tidyverse)
library(pheatmap)
library(Rtsne)
library(umap)
library(ggplot2)
library(readxl)
#### Load the excel sheet dataframe #### 
berzin <- read_xlsx("/Path to file/Cell_Line_DTP.xlsx")
### Load signature 
clstr1 <- read.table("/Gene signature/clst1_sig.tsv", sep="\t", header=T, check.names=F)
clstr2 <- read.table("/Gene Signature/clst2_sig.tsv", sep="\t", header=T, check.names=F)
clstr3 <- read.table("/Gene Signature/clst3_sig.tsv", sep="\t", header=T, check.names=F)
# Convert Berzin data to matrix format
berzin_matrix <- as.matrix(berzin[, -1])
rownames(berzin_matrix) <- berzin$Gene


###### Second methodology based on the science paper #####
library(tidyverse)
library(pheatmap)

# Assuming berzin, clstr1, clstr2, and clstr3 are already loaded

# Function to get signature genes from Berzin dataset
get_signature_genes <- function(cluster_data, berzin_data, top_n = 50) {
  genes <- cluster_data$Gene
  berzin_subset <- berzin_data[berzin_data$Gene %in% genes, ]
  
  # Calculate mean expression and variance for each gene
  gene_stats <- berzin_subset %>%
    pivot_longer(cols = -Gene, names_to = "Sample", values_to = "Expression") %>%
    group_by(Gene) %>%
    summarise(
      mean_expr = mean(Expression),
      var_expr = var(Expression)
    ) %>%
    arrange(desc(mean_expr), desc(var_expr))
  
  # Select top genes based on mean expression and variance
  top_genes <- head(gene_stats, top_n)$Gene
  
  return(top_genes)
}
# Get signature genes for each cluster
cluster1_genes <- get_signature_genes(clstr1, berzin)
cluster2_genes <- get_signature_genes(clstr2, berzin)
cluster3_genes <- get_signature_genes(clstr3, berzin)

# Combine all signature genes
all_signature_genes <- unique(c(cluster1_genes, cluster2_genes, cluster3_genes))

# Extract expression data for signature genes
signature_expr <- berzin[berzin$Gene %in% all_signature_genes, ]
rownames(signature_expr) <- signature_expr$Gene
signature_expr_2 <- signature_expr[, -1]  # Remove Gene column
#### Setting rownames 
rownames(signature_expr_2) <- signature_expr$Gene
# Calculate z-scores
signature_expr_z <- t(scale(t(signature_expr_2)))

# Create annotation for genes
gene_annotation <- data.frame(
  Cluster = case_when(
    rownames(signature_expr_z) %in% cluster1_genes ~ "Cluster1",
    rownames(signature_expr_z) %in% cluster2_genes ~ "Cluster2",
    rownames(signature_expr_z) %in% cluster3_genes ~ "Cluster3"
  )
)

# Plot heatmap
#pheatmap(signature_expr_z,
#         show_rownames = FALSE,
#         show_colnames = FALSE,
#        annotation_row = gene_annotation,
#         main = "Expression of Signature Genes in Berzin Dataset",
#         color = colorRampPalette(c("blue", "white", "red"))(100))

# Function to calculate cluster scores
calculate_cluster_score <- function(expr_data, gene_list) {
  subset_data <- expr_data[rownames(expr_data) %in% gene_list, ]
  colMeans(subset_data)
}

# Calculate cluster scores
cluster1_scores <- calculate_cluster_score(signature_expr_z, cluster1_genes)
cluster2_scores <- calculate_cluster_score(signature_expr_z, cluster2_genes)
cluster3_scores <- calculate_cluster_score(signature_expr_z, cluster3_genes)

# Combine scores
cluster_scores <- data.frame(
  Sample = colnames(signature_expr_z),
  Cluster1 = cluster1_scores,
  Cluster2 = cluster2_scores,
  Cluster3 = cluster3_scores
)

# Assign samples to clusters
cluster_scores$Assigned_Cluster <- apply(cluster_scores[, c("Cluster1", "Cluster2", "Cluster3")], 1, which.max)
cluster_scores$Assigned_Cluster <- paste0("Cluster", cluster_scores$Assigned_Cluster)

# Print summary of cluster assignments
table(cluster_scores$Assigned_Cluster)
write.table(cluster_scores, "/Path to NTP classification/cluster_assignment_cell_line_DTP.tsv",
            row.names = FALSE, sep = "\t")

# Plot cluster scores
cluster_scores_long <- pivot_longer(cluster_scores, cols = c(Cluster1, Cluster2, Cluster3), 
                                    names_to = "Cluster", values_to = "Score")

ggplot(cluster_scores_long, aes(x = Cluster, y = Score, fill = Cluster)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Distribution of Cluster Scores", x = "Cluster", y = "Score")

#### saving gene expression matrix as well ####
# Add cluster information to the z-scored expression matrix
signature_expr_z_with_cluster <- as.data.frame(signature_expr_z)
signature_expr_z_with_cluster$Gene <- rownames(signature_expr_z_with_cluster)  # Add gene names as a column
signature_expr_z_with_cluster$Cluster <- gene_annotation$Cluster  # Add cluster information

# Reorder columns to make it more readable (Gene, Cluster, followed by expression values)
signature_expr_z_with_cluster <- signature_expr_z_with_cluster %>%
  relocate(Gene, Cluster)

# Save the dataset to a CSV file
write.csv(signature_expr_z_with_cluster, "/Path to saved file/signature_expr_z_with_clusters_cell_line_DTP.csv", row.names = FALSE)

# Print a message to confirm saving
cat("The z-scored expression matrix with cluster information has been saved as 'signature_expr_z_with_clusters.csv'.\n")

