# Load required libraries
library(Rtsne)
library(ggplot2)
library(ConsensusClusterPlus)

# 1. Load your ATRT methylation beta values from your file.
# Adjust this if your RDS file contains a list rather than a simple matrix.
CRINET_betas <- readRDS("/Users/mariadanielahernandez/Downloads/metadata_betas.rds")

# 2. Filter for the 5000 most variable probes based on standard deviation
betassd <- apply(CRINET_betas, 1, sd)
betas_redmat <- CRINET_betas[rev(order(betassd))[1:5000],]

# Check for any NA values in the column names (sample names)
na_columns <- which(is.na(colnames(betas_redmat)))
if(length(na_columns) > 0){
  cat("NA found in sample names at positions:", na_columns, "\n")
}

# 3. tSNE Analysis
# tSNE expects samples as rows, so we transpose the matrix

# -----------------------------
# 1. Remove duplicates and run tSNE
# -----------------------------

# Transpose the reduced beta matrix so that samples are rows
tsne_input <- t(betas_redmat)

# Remove duplicate rows (i.e. duplicate samples)
unique_tsne_input <- tsne_input[!duplicated(tsne_input), ]

# Print number of unique samples
cat("Number of unique samples for tSNE:", nrow(unique_tsne_input), "\n")

# Run tSNE on the unique samples
set.seed(42)  # For reproducibility
result_tSNE <- Rtsne(unique_tsne_input, perplexity = 12)

# Create a data frame with tSNE coordinates and sample names
tsne_df <- data.frame(
  Sample = rownames(unique_tsne_input),
  coord1 = result_tSNE$Y[, 1],
  coord2 = result_tSNE$Y[, 2]
)

# Plot tSNE result
library(ggplot2)
tsne_plot <- ggplot(tsne_df, aes(x = coord1, y = coord2)) +
  geom_point(size = 2, color = "blue") +
  theme_minimal() +
  labs(title = "tSNE Plot of ATRT Methylation Data",
       x = "tSNE 1", y = "tSNE 2")
print(tsne_plot)





# Load required libraries
library(Rtsne)
library(ggplot2)
library(ConsensusClusterPlus)
library(RColorBrewer)  # For color palettes

# -----------------------------
# 1. Load Data and Filter Probes
# -----------------------------

# Load your ATRT methylation beta values from your file.
CRINET_betas <- readRDS("/Users/mariadanielahernandez/Downloads/metadata_betas.rds")

# Calculate standard deviation for each probe and select top 5000 most variable probes
betassd <- apply(CRINET_betas, 1, sd)
betas_redmat <- CRINET_betas[rev(order(betassd))[1:5000],]

# Check for NA in column names (sample names)
na_columns <- which(is.na(colnames(betas_redmat)))
if(length(na_columns) > 0){
  cat("NA found in sample names at positions:", na_columns, "\n")
}

# -----------------------------
# 2. tSNE Analysis
# -----------------------------

# Transpose so that samples are rows for tSNE
tsne_input <- t(betas_redmat)

# Remove duplicate samples (rows)
unique_tsne_input <- tsne_input[!duplicated(tsne_input), ]
cat("Number of unique samples for tSNE:", nrow(unique_tsne_input), "\n")

# Run tSNE on the unique samples
set.seed(42)  # For reproducibility
result_tSNE <- Rtsne(unique_tsne_input, perplexity = 12)

# Create a data frame with tSNE coordinates and sample names
tsne_df <- data.frame(
  Sample = rownames(unique_tsne_input),
  coord1 = result_tSNE$Y[, 1],
  coord2 = result_tSNE$Y[, 2]
)

# -----------------------------
# 3. Consensus Clustering
# -----------------------------
# Remove duplicate columns (samples) to match with tSNE unique samples.
unique_columns <- !duplicated(t(betas_redmat))
betas_for_consensus <- betas_redmat[, unique_columns]

# Convert to matrix explicitly
betas_for_consensus <- as.matrix(betas_for_consensus)

# Run Consensus Clustering using Pearson distance
ATRT_consensus <- ConsensusClusterPlus(betas_for_consensus,
                                       maxK = 6,
                                       reps = 1000,
                                       pItem = 0.8,
                                       pFeature = 1,
                                       title = "ATRT",
                                       clusterAlg = "pam", 
                                       distance = "pearson",
                                       writeTable = TRUE,
                                       seed = 1262118388.71279)

# Choose a consensus result, e.g., for k = 4 clusters
consensus_clusters <- ATRT_consensus[[4]]$consensusClass

# Check that the sample names in consensus clustering match those in tSNE.
consensus_sample_names <- colnames(betas_for_consensus)
cat("Consensus clustering on", length(consensus_sample_names), "samples.\n")

# Merge consensus clusters with tSNE results (assumes sample names are identical)
tsne_df$Cluster <- consensus_clusters[tsne_df$Sample]

# -----------------------------
# 4. Customized tSNE Plot (No Sample Names)
# -----------------------------
tsne_cluster_plot <- ggplot(tsne_df, aes(x = coord1, y = coord2, color = as.factor(Cluster))) +
  geom_point(size = 3) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_color_brewer(palette = "Set1") +
  labs(title = "tSNE Plot with Consensus Clusters (ATRT)",
       x = "tSNE 1",
       y = "tSNE 2",
       color = "Cluster")

# Display the plot
print(tsne_cluster_plot)

# -----------------------------
# 5. Save Plot to PDF
# -----------------------------
pdf("ATRT_tSNE_Consensus_Clusters.pdf", width = 7, height = 7)
print(tsne_cluster_plot)
dev.off()



# Load required libraries
library(readxl)    # For reading Excel files
library(dplyr)     # For data manipulation

# -----------------------------
# 1. Read the Excel File
# -----------------------------
# Use the updated file path for the Excel file.
idat_info <- read_excel("/Users/mariadanielahernandez/Downloads/IDATs_ATRT.xlsx")

# Inspect the first few rows and the column names to verify the correct column is used.
head(idat_info)
colnames(idat_info)

# -----------------------------
# 2. Extract Sample IDs from the Beta Matrix
# -----------------------------
# Assuming that 'betas_for_consensus' is the beta matrix used in consensus clustering
# and that its column names are the sample IDs.
beta_sample_ids <- colnames(betas_for_consensus)

# -----------------------------
# 3. Check for Matching Sample IDs
# -----------------------------
# Here we assume the Excel file contains a column named "IDAT Methylierung" with sample IDs.
matching_samples <- intersect(beta_sample_ids, idat_info$`IDAT Methylierung`)
cat("Number of matching samples:", length(matching_samples), "\n")

# -----------------------------
# 4. Ensure Consensus Clusters Have Sample IDs as Names
# -----------------------------
# If not already done, assign names to the consensus_clusters vector using beta sample IDs.
# This example assumes that the consensus_clusters vector came from:
# consensus_clusters <- ATRT_consensus[[4]]$consensusClass
# and that its names attribute is currently NULL.
if (is.null(names(consensus_clusters))) {
  names(consensus_clusters) <- beta_sample_ids
}

# -----------------------------
# 5. Assign Consensus Clusters to the Matching Samples
# -----------------------------
# Create a new column in the Excel data to store the cluster assignments.
idat_info <- idat_info %>%
  mutate(ConsensusCluster = consensus_clusters[match(`IDAT Methylierung`, names(consensus_clusters))])

# Check the result
head(idat_info)

# -----------------------------
# 6. (Optional) Save the Updated Data
# -----------------------------
# Save the updated Excel data as a CSV file for further analysis.
write.csv(idat_info, "IDATs_ATRT_with_Clusters.csv", row.names = FALSE)






# Load ggrepel for non-overlapping text labels
library(ggrepel)

# Create a new column "Label" to include the sample name only for matching samples
tsne_df <- tsne_df %>%
  mutate(Label = ifelse(Sample %in% matching_samples, Sample, ""))

# Create the tSNE plot with labels for matching samples using ggrepel
tsne_cluster_plot_labeled <- ggplot(tsne_df, aes(x = coord1, y = coord2, color = as.factor(Cluster))) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = Label), size = 3, 
                  force = 1,         # adjust force if needed for spacing
                  max.overlaps = Inf) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_color_brewer(palette = "Set1") +
  labs(title = "tSNE Plot with Consensus Clusters (ATRT)\n(Matching Samples Labeled)",
       x = "tSNE 1", y = "tSNE 2", color = "Cluster")

# Display the plot
print(tsne_cluster_plot_labeled)

# Optionally, save the plot as a PDF
pdf("ATRT_tSNE_Consensus_Clusters_MatchingLabeled.pdf", width = 7, height = 7)
print(tsne_cluster_plot_labeled)
dev.off()




# Filter to only include matching samples
result_table_matched <- idat_info %>% 
  filter(`IDAT Methylierung` %in% matching_samples) %>%
  select(`IDAT Methylierung`, Subgruppe, ConsensusCluster)

# Print the resulting table
print(result_table_matched)

# Optionally, display the table in a formatted way (e.g., using knitr)
library(knitr)
kable(result_table_matched, caption = "Matched Samples: Subgroups and Consensus Clusters")


library(dplyr)
library(ggplot2)
library(ggrepel)

# Assuming tsne_df contains the tSNE results with columns:
#   Sample, coord1, coord2, Cluster
# and result_table_matched contains:
#   `IDAT Methylierung`, Subgruppe, ConsensusCluster

# Merge the tSNE results with the Excel matched info based on sample IDs
tsne_matched <- tsne_df %>%
  inner_join(result_table_matched, by = c("Sample" = "IDAT Methylierung"))

# Create a new column "Label" that will hold the subgroup value for the matched samples.
# For samples not in the matched set, you can leave the label empty.
tsne_df <- tsne_df %>%
  mutate(Label = ifelse(Sample %in% tsne_matched$Sample, 
                        tsne_matched$Subgruppe[match(Sample, tsne_matched$Sample)],
                        ""))

# Now, plot the full tSNE with cluster colors and add labels (subgroup) for matched samples.
tsne_cluster_plot_with_subgroup <- ggplot(tsne_df, aes(x = coord1, y = coord2, color = as.factor(Cluster))) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = Label), size = 3, max.overlaps = Inf) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_color_brewer(palette = "Set1") +
  labs(title = "tSNE Plot with Consensus Clusters and Subgroups (Matched Samples)",
       x = "tSNE 1", y = "tSNE 2", color = "Cluster")

# Display the plot
print(tsne_cluster_plot_with_subgroup)

# Optionally, save the plot as a PDF:
pdf("ATRT_tSNE_Consensus_Clusters_With_Subgroups.pdf", width = 7, height = 7)
print(tsne_cluster_plot_with_subgroup)
dev.off()



# Load required libraries
library(uwot)
library(ggplot2)
library(ggrepel)
library(dplyr)

# -----------------------------
# 1. Run UMAP Analysis
# -----------------------------
# Use the same unique_tsne_input from before (samples as rows) for consistency.
# unique_tsne_input <- tsne_input[!duplicated(tsne_input), ]
# Run UMAP on unique samples; adjust n_neighbors and min_dist as desired.
set.seed(42)  # For reproducibility
umap_results <- umap(unique_tsne_input, n_neighbors = 15, min_dist = 0.1)

# Create a UMAP results data frame
umap_df <- data.frame(
  Sample = rownames(unique_tsne_input),
  UMAP1 = umap_results[, 1],
  UMAP2 = umap_results[, 2]
)

# -----------------------------
# 2. Add Consensus Cluster Assignments
# -----------------------------
# Here we assume that the consensus_clusters vector is already named with sample IDs.
# If not, assign them as shown earlier.
umap_df$Cluster <- consensus_clusters[umap_df$Sample]

# -----------------------------
# 3. Merge with Excel Matched Data for Subgroups
# -----------------------------
# result_table_matched is assumed to contain the matched samples from the Excel file:
# columns: `IDAT Methylierung`, Subgruppe, ConsensusCluster.
# Merge UMAP results with the matched Excel data by sample ID.
umap_matched <- umap_df %>%
  inner_join(result_table_matched, by = c("Sample" = "IDAT Methylierung"))

# Add a Label column: for samples in the matched set, use the "Subgruppe" value; otherwise, leave blank.
umap_df <- umap_df %>%
  mutate(Label = ifelse(Sample %in% umap_matched$Sample,
                        umap_matched$Subgruppe[match(Sample, umap_matched$Sample)],
                        ""))

# -----------------------------
# 4. Create the UMAP Plot with Subgroup Labels for Matched Samples
# -----------------------------
umap_plot <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = as.factor(Cluster))) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = Label), size = 3, max.overlaps = Inf) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_color_brewer(palette = "Set1") +
  labs(title = "UMAP Plot with Consensus Clusters and Subgroups (Matched Samples)",
       x = "UMAP1", y = "UMAP2", color = "Cluster")

# Display the plot
print(umap_plot)

# -----------------------------
# 5. Save the UMAP Plot to PDF
# -----------------------------
pdf("ATRT_UMAP_Consensus_Clusters_With_Subgroups.pdf", width = 7, height = 7)
print(umap_plot)
dev.off()
















