library(Seurat)
library(harmony)
library(ggplot2)
library(enrichR)
library(RColorBrewer)
library(dplyr)


path <- "/Users/mariadanielahernandez/Desktop/SeuratATRT"

# Load all samples
I007_017 <- readRDS(paste0(path, "/I007_017.Rds"))
I010_019 <- readRDS(paste0(path, "/I010_019.Rds"))
I014_002 <- readRDS(paste0(path, "/I014_002.Rds"))
I018_063 <- readRDS(paste0(path, "/I018_063.Rds"))
I022_056 <- readRDS(paste0(path, "/I022_056.Rds"))
I034_060 <- readRDS(paste0(path, "/I034_060.Rds"))
I036_058 <- readRDS(paste0(path, "/I036_058.Rds"))
I044_012 <- readRDS(paste0(path, "/I044_012.Rds"))
I056_040 <- readRDS(paste0(path, "/I056_040.Rds"))
I074_008 <- readRDS(paste0(path, "/I074_008.Rds"))
I074_013 <- readRDS(paste0(path, "/I074_013.Rds"))
I084_007 <- readRDS(paste0(path, "/I084_007.Rds"))
I040_006 <- readRDS(paste0(path, "/I040_006.Rds"))
I040_006P <- readRDS(paste0(path, "/I040_006P.Rds"))
I056_040P <- readRDS(paste0(path, "/I056_040P.Rds"))





all_samples <- merge(I007_017, y = c(I010_019, I014_002, I018_063, I022_056,I034_060,I036_058, I044_012,I056_040, I074_008,I074_013,I084_007, I040_006,I040_006P,I056_040P),
                     add.cell.ids = c("I007_017","I010_019"," I014_002", "I018_063", "I022_056","I034_060","I036_058", "I044_012","I056_040","I074_008","I074_013","I084_007","I040_006","I040_006P","I056_040P"),
                     project = "ATRT_all_samples")


# Add DNA methylation subgroup annotations
all_samples$methylation_subgroup <- NA
all_samples$methylation_subgroup[all_samples$orig.ident == "I007_017"] <- "ATRT-SHH"
all_samples$methylation_subgroup[all_samples$orig.ident == "I010_019"] <- "ATRT-TYR"
all_samples$methylation_subgroup[all_samples$orig.ident == "I014_002"] <- "ATRT-SHH"
all_samples$methylation_subgroup[all_samples$orig.ident == "I018_063"] <- "ATRT-SHH"
all_samples$methylation_subgroup[all_samples$orig.ident == "I022_056"] <- "ATRT-MYC"
all_samples$methylation_subgroup[all_samples$orig.ident == "I034_060"] <- "ATRT-SHH"
all_samples$methylation_subgroup[all_samples$orig.ident == "I036_058"] <- "ATRT-MYC"
all_samples$methylation_subgroup[all_samples$orig.ident == "I044_012"] <- "ATRT-SHH"
all_samples$methylation_subgroup[all_samples$orig.ident == "I056_040"] <- "ATRT-SHH"
all_samples$methylation_subgroup[all_samples$orig.ident == "I074_008"] <- "ATRT-TYR"
all_samples$methylation_subgroup[all_samples$orig.ident == "I074_013"] <- "ATRT-TYR"
all_samples$methylation_subgroup[all_samples$orig.ident == "I084_007"] <- "NA"
all_samples$methylation_subgroup[all_samples$orig.ident == "I040_006"] <- "ATRT-TYR"
all_samples$methylation_subgroup[all_samples$orig.ident == "I040_006P"] <-"ATRT-TYR"
all_samples$methylation_subgroup[all_samples$orig.ident == "I056_040P"] <-"ATRT-SHH"



# Display the head of methylation subgroup column to verify

head(all_samples[["methylation_subgroup"]])








all_samples$annot <- NA
all_samples$annot[all_samples$orig.ident == "I007_017"] <- "Relapse"
all_samples$annot[all_samples$orig.ident == "I010_019"] <- "Relapse"
all_samples$annot[all_samples$orig.ident == "I014_002"] <- "Relapse"
all_samples$annot[all_samples$orig.ident == "I018_063"] <- "Relapse"
all_samples$annot[all_samples$orig.ident == "I022_056"] <- "Metastasis"
all_samples$annot[all_samples$orig.ident == "I034_060"] <- "Primary"
all_samples$annot[all_samples$orig.ident == "I036_058"] <- "Relapse"
all_samples$annot[all_samples$orig.ident == "I044_012"] <- "Relapse"
all_samples$annot[all_samples$orig.ident == "I056_040"] <- "Relapse"
all_samples$annot[all_samples$orig.ident == "I074_008"] <- "Relapse"
all_samples$annot[all_samples$orig.ident == "I074_013"] <- "Relapse"
all_samples$annot[all_samples$orig.ident == "I084_007"] <- "Relapse"
all_samples$annot[all_samples$orig.ident == "I040_006"] <- "Relapse"
all_samples$annot[all_samples$orig.ident == "I040_006P"] <- "Primary"
all_samples$annot[all_samples$orig.ident == "I056_040P"] <- "Primary"


all_samples <- NormalizeData(all_samples)
all_samples <- FindVariableFeatures(all_samples)
all_samples <- ScaleData(all_samples)


#PCA
all_samples <- RunPCA(all_samples)
#Harmony for batch correction based on 'orig.ident' 
all_samples <- RunHarmony(all_samples, group.by.vars = "orig.ident", plot_convergence = TRUE)
all_samples <- RunUMAP(all_samples, reduction = "harmony", dims = 1:20)
all_samples <- FindNeighbors(all_samples, reduction = "harmony", dims = 1:20)
all_samples <- FindClusters(all_samples, resolution = 0.5)


#UMAP Plot seurat clusters
DimPlot(all_samples, reduction = "umap", group.by = "orig.ident", label = FALSE)


#UMAP Plot with condition annotations
DimPlot(all_samples, reduction = "umap", group.by = "annot", label = TRUE)

#UMAP Plot with methylation subgroups annotations
#DimPlot(all_samples, reduction = "umap", group.by = "methylation_subgroup", label = TRUE)

#UMAP Plot seurat clusters
DimPlot(all_samples, reduction = "umap", group.by = "seurat_clusters", label = FALSE)

#UMAP Plot with methylation subgroups annotations
DimPlot(all_samples, reduction = "umap", group.by = "methylation_subgroup", label = TRUE)



# Define the PDF file path
pdf("/Users/mariadanielahernandez/ATRT/umap_plots.pdf", width = 10, height = 8)

# Generate and save UMAP plots
print(DimPlot(all_samples, reduction = "umap", group.by = "orig.ident", label = FALSE))
print(DimPlot(all_samples, reduction = "umap", group.by = "annot", label = TRUE))
print(DimPlot(all_samples, reduction = "umap", group.by = "seurat_clusters", label = FALSE))
print(DimPlot(all_samples, reduction = "umap", group.by = "methylation_subgroup", label = TRUE))

# Close the PDF device
dev.off()

# Confirm that the file was saved
cat("UMAP plots saved successfully at /Users/mariadanielahernandez/ATRT/umap_plots.pdf\n")




#join layers
all_samples <- JoinLayers(all_samples, assay = "RNA")

all_markers <- FindAllMarkers(all_samples, 
                              assay = "RNA", 
                              only.pos = TRUE,  
                              min.pct = 0.25,    
                              logfc.threshold = 0.25)

head(all_markers)





# Assuming all_markers has been generated by FindAllMarkers as specified in your script.
# Sort and select the top 10 markers for each cluster based on adjusted p-value
top_markers_per_cluster <- all_markers %>%
  group_by(cluster) %>%
  arrange(p_val_adj) %>%
  slice_head(n = 10)

# Print the top 10 markers for each cluster
print(top_markers_per_cluster)

# Print more rows, e.g., the first 20 rows
print(top_markers_per_cluster, n = 90)

# Optionally, you can output this to a more readable format or save it for further use:
write.csv(top_markers_per_cluster, "/Users/mariadanielahernandez/ATRT/Top10MarkersPerCluster.csv")

# Select top markers (e.g., 5 per cluster)
top_markers <- all_markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) %>%
  pull(gene) %>%
  unique()

# Define the output file path
output_file <- "top_markers_dotplot.pdf"

# Save the Dot Plot to a PDF
pdf(output_file, width = 10, height = 8)  # Set PDF dimensions

DotPlot(all_samples, features = top_markers, group.by = "seurat_clusters") +
  RotatedAxis() + 
  ggtitle("Top Marker Genes Across Clusters") +
  theme(axis.text.x = element_text(size = 6),   # Reduce font size for x-axis
        axis.text.y = element_text(size = 6),   # Reduce font size for y-axis
        plot.title = element_text(size = 10))   # Reduce title font size

dev.off()

# Confirmation
cat("Dot Plot saved successfully as:", output_file, "\n")


#DEG SHH TYR
# Ensure all necessary libraries are loaded
library(Seurat)

comparison_object <- all_samples


# Join necessary layers. Adjust the assay and layers as per your specific setup.
comparison_object <- JoinLayers(comparison_object, assay = "RNA", layers = c("data", "scale.data"))

# Proceed with normalization and feature selection if not done
comparison_object <- NormalizeData(comparison_object)
comparison_object <- FindVariableFeatures(comparison_object)
comparison_object <- ScaleData(comparison_object)

# Perform PCA and clustering as required
comparison_object <- RunPCA(comparison_object)
comparison_object <- RunUMAP(comparison_object, dims = 1:20)

####

unique(comparison_object$orig.ident)  # or replace 'ident' with the actual column name you are using
unique(comparison_object$methylation_subgroup)
deg_results <- FindMarkers(comparison_object, 
                           ident.1 = "ATRT-SHH", ident.2 = "ATRT-TYR", 
                           min.pct = 0.25, logfc.threshold = 0.25, 
                           group.by = "methylation_subgroup")

comparison_object$ident <- comparison_object$methylation_subgroup  # Set identities based on methylation subgroup
# Check if the methylation_subgroup column is assigned correctly
table(comparison_object$methylation_subgroup)
# Set the identities using the methylation_subgroup column
comparison_object <- SetIdent(comparison_object, value = "methylation_subgroup")

# Now run the FindMarkers function
deg_results <- FindMarkers(comparison_object, ident.1 = "ATRT-SHH", ident.2 = "ATRT-TYR", 
                           min.pct = 0.25, logfc.threshold = 0.25)
# Clean any extra spaces in the methylation_subgroup column
comparison_object$methylation_subgroup <- trimws(comparison_object$methylation_subgroup)

# Now try the FindMarkers function again
deg_results <- FindMarkers(comparison_object, ident.1 = "ATRT-SHH", ident.2 = "ATRT-TYR", 
                           min.pct = 0.25, logfc.threshold = 0.25)
# Check the current identities
table(Idents(comparison_object))
deg_results <- FindMarkers(comparison_object, ident.1 = "ATRT-SHH", ident.2 = "ATRT-TYR", 
                           min.pct = 0.25, logfc.threshold = 0.25)


# Display the results
head(deg_results[order(deg_results$p_val_adj), ])

# Optionally, save the results
write.csv(deg_results, "/Users/mariadanielahernandez/ATRT/SHH_vs_TYR_DEGs.csv")

print("Differential expression analysis results saved successfully.")


# Assuming you have already loaded ggplot2 and have the deg_results dataframe ready
library(ggplot2)
library(ggplot2)

# Assuming deg_df is your data frame with DEG results
deg_results$neg_log10_p_val <- -log10(deg_results$p_val_adj)  # Transform p-values for plotting
deg_results$gene <- rownames(deg_results)

# Create the volcano plot with adjusted aesthetics
volcano_plot <- ggplot(deg_results, aes(x = avg_log2FC, y = neg_log10_p_val)) +
  geom_point(aes(color = p_val_adj < 0.05), alpha = 0.5) +  # Color code significant points
  scale_color_manual(values = c("grey", "red")) +
  labs(x = "Average Log2 Fold Change", y = "-Log10 Adjusted P-value", title = "Volcano Plot of DEGs between SHH and TYR") +
  theme_minimal() +
  geom_text(data = subset(deg_results, p_val_adj < 0.05 & abs(avg_log2FC) > 1),  # Label highly significant and largely changed genes
            aes(label = gene), vjust = 1.5, hjust = 1.5, check_overlap = TRUE, size = 3)

# Save the plot to a PDF with adjusted dimensions
ggsave("/Users/mariadanielahernandez/ATRT/volcano_plot_SHH_TYR.pdf", plot = volcano_plot, device = "pdf", width = 12, height = 8)  # Adjust dimensions as needed

print("Volcano plot saved as PDF successfully.")


library(ggplot2)
library(Seurat)
library(gridExtra)

# Define genes of interest
genes_of_interest <- c("MYCN", "GLI2", "ASCL1", "HES5", "HES6", "DLL1", "DLL3", 
                       "EZH2", "SUZ12", "EED", "AURKA", "HDAC1", "HDAC2", "ERBB2", "VEGFA")

# Set identity to methylation subgroup
comparison_object <- SetIdent(comparison_object, value = "methylation_subgroup")

# Subset only ATRT-SHH and ATRT-TYR
comparison_object_subset <- subset(comparison_object, idents = c("ATRT-SHH", "ATRT-TYR"))

# Define colors for SHH and TYR only
color_palette <- c("ATRT-SHH" = "red", "ATRT-TYR" = "blue")

# Save multiple violin plots per page in a PDF
pdf("/Users/mariadanielahernandez/ATRT/violin_plots_SHH_vs_TYR.pdf", width = 10, height = 12)

plot_list <- list()  # Store plots for grid arrangement
counter <- 1  # Track the number of plots per page

for (gene in genes_of_interest) {
  p <- VlnPlot(comparison_object_subset, features = gene, group.by = "methylation_subgroup", 
               pt.size = 1,  # Add dots for individual cells
               cols = color_palette) + 
    ggtitle(paste("Expression of", gene)) + 
    theme_minimal(base_size = 14) + 
    theme(panel.grid.major = element_blank(),  # Remove grid lines
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"))
  
  plot_list[[counter]] <- p
  counter <- counter + 1
  
  # If 4 plots per page, arrange them in a grid and reset
  if (counter > 4) {
    grid.arrange(grobs = plot_list, ncol = 2, nrow = 2)
    plot_list <- list()  # Reset list
    counter <- 1
  }
}

# Print remaining plots (if any left)
if (length(plot_list) > 0) {
  grid.arrange(grobs = plot_list, ncol = 2, nrow = ceiling(length(plot_list)/2))
}

dev.off()

cat("Violin plots saved successfully at /Users/mariadanielahernandez/ATRT/violin_plots_SHH_vs_TYR.pdf\n")


library(ggplot2)
library(Seurat)
library(gridExtra)
library(cowplot)  # For better figure arrangement

# Define genes of interest
genes_of_interest <- c("MYCN", "GLI2", "ASCL1", "HES5", "HES6", "DLL1", "DLL3", 
                       "EZH2", "SUZ12", "EED", "AURKA", "HDAC1", "HDAC2", "ERBB2", "VEGFA")

# Set identity to methylation subgroup
comparison_object <- SetIdent(comparison_object, value = "methylation_subgroup")

# Subset only ATRT-SHH and ATRT-TYR
comparison_object_subset <- subset(comparison_object, idents = c("ATRT-SHH", "ATRT-TYR"))

# Define professional color palette
color_palette <- c("ATRT-SHH" = "#D55E00", "ATRT-TYR" = "#0072B2")  # Orange & Blue for contrast

# Store individual plots
plot_list <- list()
counter <- 1

for (gene in genes_of_interest) {
  p <- VlnPlot(comparison_object_subset, features = gene, group.by = "methylation_subgroup", 
               pt.size = 0.5,  # Keep dots small for readability
               cols = color_palette) + 
    ggtitle(paste("Expression of", gene)) + 
    theme_classic(base_size = 16) +  # Use classic theme for cleaner plots
    theme(
      axis.title = element_text(size = 18, face = "bold"),
      axis.text = element_text(size = 14),
      legend.title = element_blank(),  # Remove redundant legend title
      legend.text = element_text(size = 14),
      panel.grid = element_blank()  # Remove background grid
    )
  
  plot_list[[counter]] <- p
  counter <- counter + 1
}

# Arrange plots in 2x2 grids per page and save to PDF & PNG
output_pdf <- "/Users/mariadanielahernandez/ATRT/violin_plots_SHH_vs_TYR.pdf"
output_png <- "/Users/mariadanielahernandez/ATRT/violin_plots_SHH_vs_TYR.png"

pdf(output_pdf, width = 14, height = 10)
for (i in seq(1, length(plot_list), by = 4)) {
  print(plot_grid(plotlist = plot_list[i:min(i+3, length(plot_list))], ncol = 2, nrow = 2))
}
dev.off()

# Save as high-resolution PNG for preview
png(output_png, width = 14, height = 10, units = "in", res = 600)
for (i in seq(1, length(plot_list), by = 4)) {
  print(plot_grid(plotlist = plot_list[i:min(i+3, length(plot_list))], ncol = 2, nrow = 2))
}
dev.off()

cat("Publication-ready violin plots saved at:\n", output_pdf, "\n", output_png, "\n")

library(ggplot2)

# Assuming deg_results is your data frame with DEG results
deg_results$neg_log10_p_val <- -log10(deg_results$p_val_adj)  # Transform p-values
deg_results$gene <- rownames(deg_results)  # Extract gene names

# Create the volcano plot with same aesthetics but without the grid
volcano_plot <- ggplot(deg_results, aes(x = avg_log2FC, y = neg_log10_p_val)) +
  geom_point(aes(color = p_val_adj < 0.05), alpha = 0.5) +  # Color code significant points
  scale_color_manual(values = c("grey", "red")) +
  labs(x = "Average Log2 Fold Change", 
       y = "-Log10 Adjusted P-value", 
       title = "Volcano Plot of DEGs between SHH and TYR") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank()) +  # Remove minor grid lines
  geom_text(data = subset(deg_results, p_val_adj < 0.05 & abs(avg_log2FC) > 1),  
            aes(label = gene), vjust = 1.5, hjust = 1.5, check_overlap = TRUE, size = 3) +
  annotate("text", x = max(deg_results$avg_log2FC) * 0.8, y = -1, 
           label = "ATRT-TYR", color = "black", size = 6, fontface = "bold") + 
  annotate("text", x = min(deg_results$avg_log2FC) * 0.8, y = -1, 
           label = "ATRT-SHH", color = "black", size = 6, fontface = "bold")

# Save the plot to a PDF with adjusted dimensions
ggsave("/Users/mariadanielahernandez/ATRT/volcano_plot_SHH_TYR.pdf", 
       plot = volcano_plot, device = "pdf", width = 12, height = 8)  

print("Volcano plot saved as PDF successfully.")




# Define marker genes for ATRT and other cell types
marker_genes <- list(
  Tumor = c("SOX2", "MYC", "MYCN", "OLIG2", "NES", "VIM", "LIN28A", "CCND1", "EZH2", "MKI67", "SMARCB1",  # Core ATRT markers
            "GLI2", "ASCL1", "HES5", "HES6", "DLL1", "DLL3",  # SHH pathway markers
            "MITF", "TYR", "MLANA", "SOX10", "GPNMB",  # TYR markers
            "N-MYC", "E2F3", "AURKA", "TOP2A",  # MYC markers
            "SUZ12", "EED", "AURKA", "HDAC1", "HDAC2", "ERBB2", "VEGFA"),  # Epigenetic regulators & growth factors
  Immune = c("CD3D", "CD3E", "CD8A", "CD19", "CD68", "AIF1", "PTPRC"),
  Stromal = c("COL1A1", "COL1A2", "DCN", "LUM", "PDGFRA", "ACTA2"),
  Neuronal = c("SYT1", "SNAP25", "MAP2", "RBFOX3", "DCX", "NEUROD1"),
  Glial = c("GFAP", "S100B", "AQP4", "OLIG1", "MBP", "SOX10")
)



library(Seurat)
library(dplyr)
library(ggplot2)

# Initialize a column for cell type annotations
all_samples$cell_type <- "Unknown"

# Annotate clusters based on marker gene expression
for (cluster in unique(all_samples$seurat_clusters)) {
  cluster_markers <- all_markers %>% filter(cluster == !!cluster)
  top_cluster_markers <- cluster_markers$gene[1:5]  # Top 5 markers for the cluster
  
  # Check which cell type markers are most represented
  cell_type_scores <- sapply(names(marker_genes), function(cell_type) {
    sum(top_cluster_markers %in% marker_genes[[cell_type]])
  })
  
  # Assign the cell type with the highest score
  if (max(cell_type_scores) > 0) {
    all_samples$cell_type[all_samples$seurat_clusters == cluster] <- names(which.max(cell_type_scores))
  }
}

# Check the distribution of annotated cell types
print(table(all_samples$cell_type))


# Identify tumor clusters
tumor_clusters <- unique(all_samples$seurat_clusters[all_samples$cell_type == "Tumor"])

# Generate UMAP plot with tumor clusters highlighted in red
umap_plot <- DimPlot(all_samples, group.by = "seurat_clusters", 
                     cells.highlight = WhichCells(all_samples, idents = tumor_clusters), 
                     cols.highlight = "red", cols = "grey") +
  ggtitle("Tumor Clusters Highlighted in Red") +
  theme_minimal()

# Add cluster labels to the plot
umap_plot <- LabelClusters(plot = umap_plot, id = "seurat_clusters", size = 6, repel = TRUE)

# Display the plot
print(umap_plot)



library(Seurat)
library(ggplot2)

# Generate UMAP plot with cell type annotations
umap_plot <- DimPlot(all_samples, group.by = "cell_type") +
  ggtitle("UMAP of Annotated Cell Types") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.background = element_blank(),  # Remove background
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12)
  )

# Add cluster labels to the plot
umap_plot <- LabelClusters(plot = umap_plot, id = "cell_type", size = 6, repel = TRUE)

# Display the plot
print(umap_plot)


# Identify tumor cells
tumor_cells <- WhichCells(all_samples, idents = unique(all_samples$seurat_clusters[all_samples$cell_type == "Tumor"]))

# Create a subset with only tumor cells
tumor_subset <- subset(all_samples, cells = tumor_cells)

# Run UMAP for the tumor subset (if not already computed)
tumor_subset <- RunUMAP(tumor_subset, dims = 1:20)  

# Generate UMAP plot for tumor cells only
tumor_umap <- DimPlot(tumor_subset, group.by = "seurat_clusters") +
  ggtitle("UMAP of ATRT Tumor Cells") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12)
  )

# Add cluster labels to the plot
tumor_umap <- LabelClusters(plot = tumor_umap, id = "seurat_clusters", size = 6, repel = TRUE)

# Display the plot
print(tumor_umap)


#pairs psudotime analysis
library(monocle3)
if (!requireNamespace("SeuratWrappers", quietly = TRUE)) {
  remotes::install_github("satijalab/seurat-wrappers")
}
library(SeuratWrappers)


# Extract tumor cells only from the primary-relapse pairs
tumor_pairs <- subset(tumor_subset, subset = orig.ident %in% c("I040_006P", "I040_006", "I056_040P", "I056_040"))

# Check the number of tumor cells per sample
table(tumor_pairs$orig.ident)

# Visualize UMAP with the primary-relapse pairs
DimPlot(tumor_pairs, group.by = "orig.ident", label = TRUE)

# Ensure Monocle3 is installed and loaded
library(monocle3)

# Set orig.ident as the identity
tumor_pairs <- SetIdent(tumor_pairs, value = "orig.ident")

# Check if identities are correctly assigned
table(Idents(tumor_pairs))
# Select root cells from primary tumors
root_cells <- WhichCells(tumor_pairs, idents = c("I040_006P", "I056_040P"))

# Verify root cells
length(root_cells)

# Order cells in pseudotime using primary tumors as the root
tumor_cds <- order_cells(tumor_cds, root_cells = root_cells)

# Plot pseudotime trajectory
plot_cells(tumor_cds, color_cells_by = "pseudotime")



library(ggplot2)
library(dplyr)

# Extract UMAP coordinates and metadata for annotation
umap_df <- as.data.frame(reducedDims(tumor_cds)$UMAP)  # Get UMAP coordinates

# Rename columns to standard names (if needed)
colnames(umap_df) <- c("UMAP1", "UMAP2")  # Adjust column names as needed

# Add `orig.ident` metadata
umap_df$orig.ident <- colData(tumor_cds)$orig.ident  

# Compute centroid positions for each `orig.ident`
label_positions <- umap_df %>%
  group_by(orig.ident) %>%
  summarize(UMAP1 = median(UMAP1), UMAP2 = median(UMAP2))  # Compute median UMAP position

# Plot pseudotime and add custom labels
plot_cells(tumor_cds, 
           color_cells_by = "pseudotime",  # Keep pseudotime colors
           show_trajectory_graph = TRUE, 
           label_groups_by_cluster = FALSE,  
           label_leaves = FALSE,  
           label_branch_points = FALSE,  
           graph_label_size = 3) +  
  geom_text(data = label_positions, aes(x = UMAP1, y = UMAP2, label = orig.ident),  
            color = "black", size = 7, fontface = "bold")  # Adjust label size and font


# Save Monocle3 pseudotime object
saveRDS(tumor_cds, "/Users/mariadanielahernandez/ATRT/tumor_pseudotime_pairs.rds")

# Extract pseudotime values
pseudotime_values <- data.frame(
  Cell = colnames(tumor_cds),
  Pseudotime = pseudotime(tumor_cds)
)

# Save pseudotime values to CSV
write.csv(pseudotime_values, "/Users/mariadanielahernandez/ATRT/pseudotime_values.csv", row.names = FALSE)

cat("Pseudotime results saved successfully!")


# Check the names of the UMAP dimensions
print(colnames(reducedDims(tumor_cds)$UMAP))

# Assuming the names might be V1, V2, rename them for clarity
umap_df <- as.data.frame(reducedDims(tumor_cds)$UMAP)
colnames(umap_df) <- c("UMAP1", "UMAP2")

# Add pseudotime and orig.ident to the dataframe
umap_df$pseudotime <- pseudotime(tumor_cds)
umap_df$orig.ident <- colData(tumor_cds)$orig.ident

# Check the first few rows to confirm all is well
head(umap_df)


# Assuming 'umap_df$is_primary' distinguishes between primary and relapse
umap_df$is_primary <- ifelse(umap_df$orig.ident %in% c("I040_006P", "I056_040P"), "Primary", "Relapse")

library(ggplot2)

ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = pseudotime)) +
  geom_point(alpha = 0.8) +
  scale_color_viridis_c(option = "magma", end = 0.9) +  # Use viridis color scale for clear visualization
  facet_wrap(~orig.ident) +  # Facet by original identifier to differentiate primary from relapse
  labs(title = "Pseudotime Trajectory with Primary and Relapse Samples",
       x = "UMAP 1", y = "UMAP 2") +
  theme_minimal() +
  theme(
    legend.position = "right",
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.background = element_rect(fill = "white", colour = "white")  # Set background color
  )



# Summary statistics
library(dplyr)
stats <- umap_df %>%
  group_by(is_primary) %>%
  summarize(mean_pseudotime = mean(pseudotime, na.rm = TRUE),
            median_pseudotime = median(pseudotime, na.rm = TRUE),
            sd_pseudotime = sd(pseudotime, na.rm = TRUE))

# Output the stats
print(stats)

# Wilcoxon test to compare distributions
test_result <- wilcox.test(pseudotime ~ is_primary, data = umap_df)
print(test_result)

#FIXED PSEUDOTIME LATEST!

###############################################################################
# 1. Load required libraries
###############################################################################
library(Seurat)
library(SeuratWrappers)  # if needed: devtools::install_github("satijalab/seurat-wrappers")
library(monocle3)        # if needed: devtools::install_github("cole-trapnell-lab/monocle3")
library(dplyr)
library(ggplot2)


###############################################################################
# 2. Subset the Seurat object to the primary & relapse samples of interest
###############################################################################
# Assume you already have a Seurat object "tumor_subset" containing your tumor cells.
# We'll define which samples correspond to primary (marked with "P") vs. relapse.

pairs_of_interest <- c("I040_006P", "I040_006", "I056_040P", "I056_040")
tumor_pairs <- subset(tumor_subset, subset = orig.ident %in% pairs_of_interest)


###############################################################################
# 3. Convert Seurat object to Monocle 3 cell_data_set
###############################################################################
cds <- as.cell_data_set(tumor_pairs)


###############################################################################
# 4. Monocle 3 Trajectory Inference
###############################################################################
# 4A. Preprocess (select top PCs, etc.)
cds <- preprocess_cds(cds, num_dim = 20)  # Adjust num_dim as appropriate

# 4B. Reduce dimensions using UMAP
cds <- reduce_dimension(cds, reduction_method = "UMAP")

# 4C. Optionally cluster cells with Monocle 3
cds <- cluster_cells(cds)

# 4D. Learn the principal graph
cds <- learn_graph(cds)


###############################################################################
# 5. Order cells in pseudotime
###############################################################################
# Here we define the root cells as those from primary samples ("P" in their orig.ident).
root_cells <- colnames(cds)[cds@colData$orig.ident %in% c("I040_006P", "I056_040P")]
cds <- order_cells(cds, root_cells = root_cells)


###############################################################################
# 6. Branch Analysis: Identify Branch-Dependent Genes
###############################################################################
# "graph_test()" checks how genes vary along branches of the trajectory graph.
branch_de <- graph_test(cds, neighbor_graph = "principal_graph", cores = 4)

# Filter for statistically significant genes (adjust threshold if needed)
sig_genes <- rownames(subset(branch_de, q_value < 0.05))
head(sig_genes)

# Monocle 3 requires 'gene_short_name' in rowData for gene labeling
if (!"gene_short_name" %in% colnames(rowData(cds))) {
  rowData(cds)$gene_short_name <- rownames(rowData(cds))
}

# Plot expression in pseudotime for a few significant genes
plot_genes_in_pseudotime(
  cds[sig_genes[1:6], ], 
  color_cells_by = "pseudotime",
  label_by_short_name = TRUE
)


###############################################################################
# 7. Create a Data Frame for Custom Pseudotime + UMAP Plot
###############################################################################
# Extract UMAP embeddings (Monocle 3's dimensional reduction)
umap_coords <- reducedDims(cds)$UMAP  # matrix with two columns: UMAP1, UMAP2

# Extract pseudotime values (numeric vector, one value per cell)
pt <- pseudotime(cds)

# Extract sample IDs (orig.ident) or other metadata from colData(cds)
sample_ids <- colData(cds)$orig.ident

# Construct a data frame for ggplot
umap_df <- data.frame(
  UMAP1 = umap_coords[, 1],
  UMAP2 = umap_coords[, 2],
  pseudotime = pt,
  orig.ident = sample_ids
)

# Create a column indicating whether sample is "Primary" or "Relapse"
umap_df$is_primary <- ifelse(
  umap_df$orig.ident %in% c("I040_006P", "I056_040P"),
  "Primary",
  "Relapse"
)


###############################################################################
# 8. Plot UMAP Faceted by Sample ID, Colored by Pseudotime
###############################################################################
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = pseudotime)) +
  geom_point(alpha = 0.8) +
  scale_color_viridis_c(option = "magma", end = 0.9) +
  facet_wrap(~ orig.ident) +
  labs(
    title = "Pseudotime Trajectory with Primary and Relapse Samples",
    x = "UMAP 1", 
    y = "UMAP 2"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", colour = "white")
  )


###############################################################################
# 9. Print Statistics for Primary vs. Relapse (Confirm Root Placement)
###############################################################################
# Summarize pseudotime distribution in each group
pseudotime_stats <- umap_df %>%
  group_by(is_primary) %>%
  summarize(
    mean_pt   = mean(pseudotime, na.rm = TRUE),
    median_pt = median(pseudotime, na.rm = TRUE),
    min_pt    = min(pseudotime, na.rm = TRUE),
    max_pt    = max(pseudotime, na.rm = TRUE),
    n_cells   = n()
  )

print(pseudotime_stats)

# Statistical test (Wilcoxon rank-sum) to see if distributions differ significantly
primary_pt <- umap_df$pseudotime[umap_df$is_primary == "Primary"]
relapse_pt <- umap_df$pseudotime[umap_df$is_primary == "Relapse"]
wilcox_res <- wilcox.test(primary_pt, relapse_pt)
wilcox_res



  
  
branch_de <- graph_test(cds, neighbor_graph = "principal_graph")
sig_genes <- rownames(subset(branch_de, q_value < 0.05))
# 1. Create a metadata column ("condition" in this example)
tumor_pairs$condition <- ifelse(tumor_pairs$orig.ident %in% c("I040_006P","I056_040P"), 
                                "Primary", 
                                "Relapse")

# 2. Set this new column as the active identity
tumor_pairs <- SetIdent(tumor_pairs, value = "condition")

# 3. Verify identities
table(Idents(tumor_pairs))  # Should show "Primary" and "Relapse"

# 4. Now run FindMarkers comparing ident.1 = "Primary" vs. ident.2 = "Relapse"
DE_results <- FindMarkers(tumor_pairs, ident.1 = "Primary", ident.2 = "Relapse")

DE_results <- FindMarkers(tumor_pairs, ident.1 = "Primary", ident.2 = "Relapse")

common_genes <- intersect(sig_genes, rownames(DE_results))
branch_relevant_DE <- DE_results[common_genes, ]
sig_common_genes_DE <- subset(branch_relevant_DE, p_val_adj < 0.05)
non_ribo_df <- subset(sig_common_genes_DE, !grepl("^RPL|^RPS", rownames(sig_common_genes_DE)))
non_ribo_df_sorted <- non_ribo_df[order(non_ribo_df$p_val_adj), ]
head(non_ribo_df_sorted)
genes_to_plot <- rownames(non_ribo_df_sorted)


# 'non_ribo_df' is the data frame of branch-dependent + DE genes
# with ribosomal (RPL/RPS) genes removed.
genes_to_plot <- rownames(non_ribo_df)
genes_to_plot
# e.g. c("HES4", "AURKAIP1", "PARK7", "ENO1", "PIK3CD", "SRM", ...)


library(monocle3)
library(dplyr)

# A. Pseudotime for each cell
pt <- pseudotime(cds)  # numeric vector named by cell IDs

# B. Metadata (including 'is_primary' or similar)
meta <- colData(cds) %>% as.data.frame()

# C. Expression matrix
expr_mat <- exprs(cds)  # or 'counts(cds)' / 'logcounts(cds)', depending on your pipeline
df_list <- lapply(genes_to_plot, function(g) {
  data.frame(
    cell       = colnames(cds),
    gene       = g,
    pseudotime = pt[colnames(cds)],
    is_primary = meta[colnames(cds), "is_primary"],  # or your grouping column
    expression = as.numeric(expr_mat[g, colnames(cds)])
  )
})

df <- do.call(rbind, df_list)
head(df)
library(ggplot2)

p <- ggplot(df, aes(x = pseudotime, y = expression, color = is_primary)) +
  geom_point(alpha = 0.6) +                           # Points for each cell
  geom_smooth(method = "loess", color = "black") +    # Black smoothing curve
  facet_wrap(~ gene, scales = "free_y") +             # Separate facet per gene
  scale_y_log10() +                                   # Log10 y-axis
  labs(
    title = "Pseudotime Plot of Non-Ribo Genes",
    x = "Pseudotime",
    y = "Expression (log10 scale)",
    color = NULL       # Removes legend title
  ) +
  theme_bw() + 
  theme(
    panel.grid.major = element_blank(),  # Remove major grid
    panel.grid.minor = element_blank(),  # Remove minor grid
    legend.title = element_blank()       # Also remove legend title
  )

print(p)

# Required library
library(pheatmap)
library(RColorBrewer)  # For more color options

# Prepare the data and annotations as previously described
genes_subset <- head(genes_to_plot, 20)  # Adjust this as needed
expr_mat <- GetAssayData(tumor_pairs, assay = "RNA", slot = "data")
heatmap_mat <- expr_mat[genes_subset, ]
pt_vals <- pseudotime(cds)
ordered_cells <- names(sort(pt_vals))
heatmap_mat <- heatmap_mat[, ordered_cells]
heatmap_mat_scaled <- t(scale(t(heatmap_mat)))

# Create a dataframe for annotations including pseudotime
annotation_col <- data.frame(
  Condition = ifelse(grepl("P$", cds@colData$orig.ident[ordered_cells]), "Primary", "Relapse"),
  Pseudotime = pt_vals[ordered_cells]
)
rownames(annotation_col) <- ordered_cells

# Color settings for annotations
annotation_colors <- list(
  Condition = c("Primary" = "green", "Relapse" = "red"),
  Pseudotime = colorRampPalette(c("blue", "yellow"))(100)  # Pseudotime color gradient
)

# --- Generate Heatmap for Primary ---
primary_cells <- ordered_cells[annotation_col$Condition == "Primary"]
heatmap_mat_primary <- heatmap_mat_scaled[, primary_cells]

# Create annotation specifically for Primary
primary_annotation <- annotation_col[primary_cells, , drop = FALSE]

pdf("heatmap_primary_with_pseudotime.pdf", width = 10, height = 8)
pheatmap(heatmap_mat_primary,
         annotation_col = primary_annotation,
         annotation_colors = annotation_colors,
         cluster_cols = FALSE,
         show_colnames = FALSE,
         main = "Heatmap of Top 20 Branch-Dependent Genes - Primary")
dev.off()
cat("Primary heatmap with pseudotime saved successfully.\n")

# --- Generate Heatmap for Relapse ---
relapse_cells <- ordered_cells[annotation_col$Condition == "Relapse"]
heatmap_mat_relapse <- heatmap_mat_scaled[, relapse_cells]

# Create annotation specifically for Relapse
relapse_annotation <- annotation_col[relapse_cells, , drop = FALSE]

pdf("heatmap_relapse_with_pseudotime.pdf", width = 10, height = 8)
pheatmap(heatmap_mat_relapse,
         annotation_col = relapse_annotation,
         annotation_colors = annotation_colors,
         cluster_cols = FALSE,
         show_colnames = FALSE,
         main = "Heatmap of Top 20 Branch-Dependent Genes - Relapse")
dev.off()
cat("Relapse heatmap with pseudotime saved successfully.\n")





#subgroup relapses key genes
library(Seurat)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(cowplot)

# Subset for relapses in SHH and TYR
relapse_subset <- subset(all_samples, subset = methylation_subgroup %in% c("ATRT-SHH", "ATRT-TYR") & annot == "Relapse")

# Set identity to methylation subgroup for the relapse subset
relapse_subset <- SetIdent(relapse_subset, value = "methylation_subgroup")

# Define the genes of interest
genes_of_interest <- c("MYCN", "GLI2", "ASCL1", "HES5", "HES6", "DLL1", "DLL3", 
                       "EZH2", "SUZ12", "EED", "AURKA", "HDAC1", "HDAC2", "ERBB2", "VEGFA")

# Define color palette
color_palette <- c("ATRT-SHH" = "#D55E00", "ATRT-TYR" = "#0072B2")

# Save multiple violin plots per page in a PDF
pdf("/Users/mariadanielahernandez/ATRT/relapse_violin_plots_SHH_vs_TYR.pdf", width = 14, height = 10)
plot_list <- list()
counter <- 1  # Initialize counter to manage plot grid

for (gene in genes_of_interest) {
  p <- VlnPlot(relapse_subset, features = gene, group.by = "methylation_subgroup", 
               pt.size = 1, cols = color_palette) + 
    ggtitle(paste("Expression of", gene)) + 
    theme_classic(base_size = 16) +
    theme(
      axis.title = element_text(size = 18, face = "bold"),
      axis.text = element_text(size = 14),
      legend.title = element_blank(),
      legend.text = element_text(size = 14),
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 0, hjust = .5)) +
    scale_x_discrete(labels = c("ATRT-SHH" = "ATRT-SHH relapses", "ATRT-TYR" = "ATRT-TYR relapses"))
  
  plot_list[[counter]] <- p
  counter <- counter + 1
  
  # Check if the plot list reaches four plots or is at the end of the gene list
  if (counter > 4 || gene == tail(genes_of_interest, n = 1)) {
    print(plot_grid(plotlist = plot_list, ncol = 2, nrow = 2))
    plot_list <- list()  # Reset the list
    counter <- 1  # Reset counter
  }
}

dev.off()

# Confirmation message
cat("Violin plots for relapse comparison between SHH and TYR saved successfully at /Users/mariadanielahernandez/ATRT/relapse_violin_plots_SHH_vs_TYR.pdf\n")

#pseudotime MONOCLE PRIMARY

###############################################################################
# 1. Load required libraries
###############################################################################
library(Seurat)
library(SeuratWrappers)
library(monocle3)
library(dplyr)
library(ggplot2)

###############################################################################
# 2. Subset the Seurat object to the primary samples of interest
###############################################################################
# Define primary sample IDs, assuming they are marked with "P" or similarly.
#primary_samples <- c("I040_006P", "I056_040P")  
primary_samples <- c("I040_006") #only one subgroup 
primary_subset <- subset(tumor_subset, subset = orig.ident %in% primary_samples)

###############################################################################
# 3. Convert Seurat object to Monocle 3 cell_data_set
###############################################################################
cds <- as.cell_data_set(primary_subset)

###############################################################################
# 4. Monocle 3 Trajectory Inference
###############################################################################
# 4A. Preprocess (select top PCs, etc.)
cds <- preprocess_cds(cds, num_dim = 20)  # Adjust num_dim as appropriate

# 4B. Reduce dimensions using UMAP
cds <- reduce_dimension(cds, reduction_method = "UMAP")

# 4C. Optionally cluster cells with Monocle 3
cds <- cluster_cells(cds)

# 4D. Learn the principal graph
cds <- learn_graph(cds)

###############################################################################
# 5. Order cells in pseudotime
###############################################################################
# Define the root cells as those from primary samples (if specific cells are known, use them)
root_cells <- colnames(cds)[cds@colData$orig.ident %in% primary_samples]
cds <- order_cells(cds, root_cells = root_cells)

###############################################################################
# 6. Branch Analysis: Identify Branch-Dependent Genes
###############################################################################
branch_de <- graph_test(cds, neighbor_graph = "principal_graph", cores = 4)

# Filter for statistically significant genes (adjust threshold if needed)
sig_genes <- rownames(subset(branch_de, q_value < 0.05))

# Monocle 3 requires 'gene_short_name' in rowData for gene labeling
if (!"gene_short_name" %in% colnames(rowData(cds))) {
  rowData(cds)$gene_short_name <- rownames(rowData(cds))
}

# Plot expression in pseudotime for a few significant genes
plot_genes_in_pseudotime(
  cds[sig_genes[1:6], ], 
  color_cells_by = "pseudotime",
  label_by_short_name = TRUE
)

###############################################################################
# 7. Create a Data Frame for Custom Pseudotime + UMAP Plot
###############################################################################
# Extract UMAP embeddings and pseudotime values
umap_coords <- reducedDims(cds)$UMAP  # matrix with two columns: UMAP1, UMAP2
pt <- pseudotime(cds)
sample_ids <- colData(cds)$orig.ident

# Construct a data frame for ggplot
umap_df <- data.frame(
  UMAP1 = umap_coords[, 1],
  UMAP2 = umap_coords[, 2],
  pseudotime = pt,
  orig.ident = sample_ids,
  is_primary = ifelse(sample_ids %in% primary_samples, "Primary", "Other")
)

###############################################################################
# 8. Plot UMAP Colored by Pseudotime, Highlighting Primary Samples
###############################################################################
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = pseudotime)) +
  geom_point(alpha = 0.8) +
  scale_color_viridis_c(option = "magma", end = 0.9) +
  labs(title = "Pseudotime Trajectory of Primary Samples", x = "UMAP 1", y = "UMAP 2") +
  theme_minimal() +
  theme(legend.position = "right",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", colour = "white"))

###############################################################################
# 9. Summary and Statistics for Primary Samples Only
###############################################################################
pseudotime_stats <- umap_df %>%
  filter(is_primary == "Relapse") %>%
  summarize(
    mean_pt   = mean(pseudotime, na.rm = TRUE),
    median_pt = median(pseudotime, na.rm = TRUE),
    min_pt    = min(pseudotime, na.rm = TRUE),
    max_pt    = max(pseudotime, na.rm = TRUE),
    n_cells   = n()
  )

print(pseudotime_stats)

# Display dimensions and a snippet of the expression matrix
expr_matrix <- exprs(cds)
print(dim(expr_matrix))  # Should show the number of genes and number of cells

# Display the first few rows and columns of the expression matrix
print(head(expr_matrix[, 1:5]))



# Calculate mean expression for genes that are expressed in at least some cells
mean_expressions <- apply(expr_matrix, 1, function(x) {
  expressed_values <- x[x > 0]  # Filter to include only expressed values
  if (length(expressed_values) > 0) mean(expressed_values, na.rm = TRUE) else 0
})

# Create a data frame and exclude genes with a mean expression of zero
gene_expression_df <- data.frame(
  gene = rownames(expr_matrix),
  expression = mean_expressions,
  stringsAsFactors = FALSE
) %>%
  dplyr::filter(expression > 0)  # Include only genes with non-zero expression

# Identify the top 10 expressed genes
top_genes <- gene_expression_df %>%
  dplyr::arrange(desc(expression)) %>%
  dplyr::slice_head(n = 10) %>%
  dplyr::pull(gene)

print(top_genes)

# Install gridExtra if you haven't already
if (!requireNamespace("gridExtra", quietly = TRUE)) {
  install.packages("gridExtra")
}

library(gridExtra)
library(ggplot2)
plot_gene_expression_pseudotime <- function(cds, gene_name) {
  df <- data.frame(
    pseudotime = pseudotime(cds),
    expression = as.numeric(exprs(cds)[gene_name, ]),
    sample_type = "Relapse"
  )
  
  df <- na.omit(df)
  df$expression_log = log10(df$expression + 1)
  
  p <- ggplot(df, aes(x = pseudotime, y = expression_log, color = sample_type)) +
    geom_point(alpha = 0.4, size = 1) +
    geom_smooth(method = 'loess', color = "black", se = FALSE) +
    labs(
      title = paste("Expression of", gene_name, "across Pseudotime"),
      x = "Pseudotime", y = "Expression (log scale)"
    ) +
    theme_minimal()
  
  return(p)
}
# Apply the plotting function to each of the top genes and collect all plots
plots <- lapply(top_genes, function(gene) plot_gene_expression_pseudotime(primary_cds, gene))
#HERE
# Save all plots in a single PDF file
pdf("/Users/mariadanielahernandez/ATRT/rcombined_expression_plots.pdf", width = 11, height = 8.5)
grid.arrange(grobs = plots, ncol = 2, nrow = 5)  # Adjust ncol and nrow as needed
dev.off()
# Assuming gene_expression_df is already created with all genes and their expressions
# Exclude "HBA2" and "HB1" and select the next top genes

excluded_genes <- c("HBA2", "HB1")

# Filter out the excluded genes and select the top 10 after exclusion
top_genes <- gene_expression_df %>%
  filter(!(gene %in% excluded_genes)) %>%
  arrange(desc(expression)) %>%
  slice_head(n = 12) %>%  # Select 12 assuming "HBA2" and "HB1" were in the initial top 12
  pull(gene)

# Check the selected top genes
print(top_genes)

# Apply the plotting function to each of the top genes and save plots
plots <- lapply(top_genes, function(gene) plot_gene_expression_pseudotime(cds, gene))

# Save all plots in a single PDF file
pdf("/Users/mariadanielahernandez/ATRT/combined_expression_plots.pdf", width = 11, height = 8.5)
grid.arrange(grobs = plots, ncol = 2, nrow = 5)  # Adjust ncol and nrow to fit your layout preference
dev.off()

# Calculate the appropriate number of rows and columns
num_plots <- length(plots)
ncol <- 2  # Preset number of columns
nrow <- ceiling(num_plots / ncol)  # Calculate the number of rows needed

# Use grid.arrange with the calculated number of rows and columns
pdf("/Users/mariadanielahernandez/ATRT/1R_combined_expression_plots.pdf", width = 11, height = 8.5)
if (num_plots > 0) {
  gridExtra::grid.arrange(grobs = plots, ncol = ncol, nrow = nrow)
} else {
  cat("No plots to arrange.")
}
dev.off()
# Assuming gene_expression_df is already created with all genes and their expressions
# Filter out genes starting with "HB" and select the top expressed genes afterwards

# Filter out the excluded genes that start with "HB"
top_genes <- gene_expression_df %>%
  filter(!grepl("^HB", gene)) %>%  # Exclude genes starting with "HB"
  arrange(desc(expression)) %>%
  slice_head(n = 10) %>%  # Adjust the number as needed to get the top N genes after filtering
  pull(gene)

# Check the selected top genes to ensure "HB" genes are excluded
print(top_genes)
# Apply the plotting function to each of the top genes and save plots
plots <- lapply(top_genes, function(gene) plot_gene_expression_pseudotime(cds, gene))

# Assuming you have a dynamic setup for grid dimensions as previously described
num_plots <- length(plots)
ncol <- 2  # Preset number of columns
nrow <- ceiling(num_plots / ncol)  # Calculate the number of rows needed

# Save all plots in a single PDF file with appropriate grid settings
pdf("/Users/mariadanielahernandez/ATRT/1Rcombined_expression_plots.pdf", width = 11, height = 8.5)
if (num_plots > 0) {
  gridExtra::grid.arrange(grobs = plots, ncol = ncol, nrow = nrow)
} else {
  cat("No plots to arrange.")
}
dev.off()

#primary relapse peseuodtime 

# Define primary sample IDs based on your dataset
#primary_samples <- c("I040_006P", "I056_040P")  # Adjust these identifiers as needed
primary_samples <- c("I056_040P") 
# Subset the Seurat object for primary samples
primary_subset <- subset(tumor_subset, subset = orig.ident %in% primary_samples)

# Convert to Monocle 3 cell_data_set
cds_primary <- as.cell_data_set(primary_subset)
# Define relapse sample IDs
#relapse_samples <- c("I040_006", "I056_040")  # Adjust these identifiers as needed
relapse_samples <- c("I056_040")
# Subset the Seurat object for relapse samples
relapse_subset <- subset(tumor_subset, subset = orig.ident %in% relapse_samples)

# Convert to Monocle 3 cell_data_set
cds_relapse <- as.cell_data_set(relapse_subset)

# Ensure UMAP dimension reduction
cds_primary <- reduce_dimension(cds_primary, reduction_method = "UMAP")
cds_relapse <- reduce_dimension(cds_relapse, reduction_method = "UMAP")

# Cluster cells if not already done
cds_primary <- cluster_cells(cds_primary)
cds_relapse <- cluster_cells(cds_relapse)

# Learn the trajectory graph
cds_primary <- learn_graph(cds_primary)
cds_relapse <- learn_graph(cds_relapse)

# Order cells to calculate pseudotime
cds_primary <- order_cells(cds_primary, reduction_method = "UMAP")
cds_relapse <- order_cells(cds_relapse, reduction_method = "UMAP")

# Function to extract expression data
extract_gene_expression <- function(cds, gene_name, condition) {
  if (!(gene_name %in% rownames(rowData(cds)))) {
    stop("Gene not found: ", gene_name)
  }
  
  # Extract pseudotime and expression
  df <- data.frame(
    pseudotime = pseudotime(cds),
    expression = as.numeric(exprs(cds)[gene_name, ]),
    condition = condition  # "Primary" or "Relapse"
  )
  
  df <- na.omit(df)  # Remove missing values
  return(df)
}

# Define genes of interest
genes_to_compare <- c("HES4", "PIK3CD")

# Extract expression data for both conditions
expression_data_primary <- bind_rows(lapply(genes_to_compare, function(gene) extract_gene_expression(cds_primary, gene, "Primary")))
expression_data_relapse <- bind_rows(lapply(genes_to_compare, function(gene) extract_gene_expression(cds_relapse, gene, "Relapse")))

# Combine both datasets
expression_data <- bind_rows(expression_data_primary, expression_data_relapse)

# Print first few rows of expression data
print(head(expression_data))

# Function to extract HES4 expression values
extract_hes4_expression <- function(cds, condition) {
  gene_name <- "HES4"  # Target gene
  
  # Check if gene exists in dataset
  if (!(gene_name %in% rownames(rowData(cds)))) {
    stop("Gene not found: ", gene_name)
  }
  
  # Extract pseudotime and expression
  df <- data.frame(
    pseudotime = pseudotime(cds),
    expression = as.numeric(exprs(cds)[gene_name, ]),
    condition = condition  # "Primary" or "Relapse"
  )
  
  df <- na.omit(df)  # Remove missing values
  return(df)
}

# Extract HES4 expression for PRIMARY and RELAPSE
hes4_primary <- extract_hes4_expression(cds_primary, "Primary")
hes4_relapse <- extract_hes4_expression(cds_relapse, "Relapse")

# Combine both datasets
hes4_expression_data <- bind_rows(hes4_primary, hes4_relapse)

# Print the first few rows of expression data
print(head(hes4_expression_data))
library(ggplot2)

ggplot(hes4_expression_data, aes(x = pseudotime, y = expression, color = condition)) +
  geom_point(alpha = 0.4, size = 1) +
  geom_smooth(method = "loess", se = FALSE) +
  labs(
    title = "HES4 Expression in Primary vs Relapse",
    x = "Pseudotime",
    y = "Expression"
  ) +
  theme_minimal() +
  scale_color_manual(values = c("Primary" = "blue", "Relapse" = "red"))



# Load required libraries
library(Seurat)
library(ggplot2)
library(dplyr)

# Define mesenchymal markers
mesenchymal_markers <- c("VIM", "CD44", "SPP1", "MSLN", "MMP14", "LAMA1", 
                         "CALD1", "S100A4", "SNAI2", "SNAI1", "TWIST1", 
                         "ZEB1", "FSP1", "COL1A1")

# Step 1: UMAP FeaturePlot
feature_plot <- FeaturePlot(all_samples, features = mesenchymal_markers, reduction = "umap", ncol = 4)

# Step 2: Violin Plot
vln_plot <- VlnPlot(all_samples, features = mesenchymal_markers, group.by = "seurat_clusters", pt.size = 0) +
  theme(text = element_text(size = 12))

# Step 3: DotPlot
dot_plot <- DotPlot(all_samples, features = mesenchymal_markers, group.by = "seurat_clusters") + 
  RotatedAxis() +
  ggtitle("Mesenchymal Marker Expression Across Clusters")

# Step 4: Save all plots in a PDF
pdf("/Users/mariadanielahernandez/ATRT/mesenchymal_markers.pdf", width = 14, height = 10)
print(feature_plot)
print(vln_plot)
print(dot_plot)
dev.off()

# Step 5: Display success message
cat("Mesenchymal marker plots saved successfully at /Users/mariadanielahernandez/ATRT/mesenchymal_markers.pdf\n")

#SUBSET TUMOR CELLS FROM THE HARMONIZED ALL SAMPLES AND SAVE IT AS RDS 