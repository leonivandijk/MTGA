# --------------------------------  preparation --------------------------------
library(pheatmap)
library(DESeq2)
library(apeglm)
library(WGCNA)
library(dplyr)
library(limma)
library(flashClust)
library(biomaRt)
library(ggplot2)
library(ggrepel)
library(gridExtra)
setwd("/Users/leonivandijk/Desktop/thesis/")

# read files
norm_counts <- read.table(
  "data/syn3388564/ROSMAP_RNAseq_FPKM_gene.tsv",
  header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE,
  row.names = 1, as.is = TRUE)
raw_counts <- read.table(
  "data/syn3388564/ROSMAP_all_counts_matrix.txt.gz",
  header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE,
  row.names = 1, as.is = TRUE)
meta1 <- read.table(
  "data/syn3388564/ROSMAP_clinical.csv",
  header = TRUE, sep = ",", check.names = FALSE, stringsAsFactors = FALSE,
  as.is = TRUE)
meta2 <- read.table(
  "data/syn3388564/ROSMAP_biospecimen_metadata.csv",
  header = TRUE, sep = ",", check.names = FALSE, stringsAsFactors = FALSE,
  as.is = TRUE)
row.names(meta1) <- meta1$individualID
metadata <- merge(meta1, meta2, by="individualID", all.x=TRUE)
metadata <- subset(subset(metadata, assay=='rnaSeq'), specimenID %in% colnames(raw_counts))
rm(meta1, meta2)


# prepare metadata
#   - extract batch number
row_index <- substring(colnames(norm_counts)[-1], 1, nchar(colnames(norm_counts)[-1])-2)
batch <- substring(colnames(norm_counts)[-1], nchar(colnames(norm_counts)[-1]), nchar(colnames(norm_counts)[-1]))
batch_mapping <- data.frame(specimenID = row_index, batch = batch)
#   - remove duplicate person
metadata <- merge(metadata, batch_mapping, by = "specimenID", all.x = TRUE)
metadata <- subset(metadata, !(metadata$specimenID == '492_120515' & metadata$batch %in% c(6,7)))
metadata$batch <- factor(metadata$batch)
#   - complete metadata
metadata <- metadata[complete.cases(metadata$pmi),]

#   - remove outlier persons (from later analysis)
#metadata <- subset(metadata, !(metadata$specimenID %in% c('380_120503', '367_120502', '500_120515', '507_120515')))

#   - factorize gender and filter on males
metadata$msex <- ifelse(metadata$msex == 0, "male", "female")
metadata <- subset(metadata, msex == 'male')
metadata$msex <- factor(metadata$msex)

#   - set diagnosis
metadata$disease <- ifelse((metadata$cogdx == 1), "control", 
                           ifelse((metadata$cogdx == 4), "alzheimer", NA))
metadata$disease <- factor(metadata$disease, levels = c("control","alzheimer"))
metadata <- subset(metadata, disease %in% c("control", "alzheimer"))
rownames(metadata) <- metadata$specimenID

rm(norm_counts, batch_mapping)


#   - print sample groups
length(rownames(subset(metadata, disease=='alzheimer')))
length(rownames(subset(metadata, disease=='control')))


# prepare expression data
raw_counts <- raw_counts[5:nrow(raw_counts), ]
raw_counts <- raw_counts[, rownames(metadata)]
#    - keep only protein_coding RNA
rownames(raw_counts) <- sub("\\..*", "", rownames(raw_counts))
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
annot <- getBM(attributes=c("ensembl_gene_id", "gene_biotype"),
               values = rownames(raw_counts),
               mart = ensembl)
rownames(annot) <- annot$ensembl_gene_id
annot <- annot[rownames(raw_counts),]
annot <- subset(annot, gene_biotype == 'protein_coding')
raw_counts <- raw_counts[rownames(annot), ]
#    - pre-filter out low count genes: keep genes with at least 10 counts > 100
#       - reduces the memory of dds, count of 10 is reasonable for bulk RNA-seq
raw_counts <- raw_counts[rowSums(raw_counts >= 10) >= 100, ]
#    - build DESeqDataSet
dds <- DESeqDataSetFromMatrix(
  countData = raw_counts,
  colData = metadata,
  design = ~ batch + disease
)

rm(ensembl, annot)

# ----------------------  differential expression analysis ---------------------
#   - normalize
dds <- DESeq(dds) # original counts preserved in counts(dds)
norm_counts <- counts(dds, normalized=TRUE)
norm_counts <- norm_counts[goodGenes(t(norm_counts)),]
#   - variance stabilizing transformation
rlog <- vst(dds, blind=FALSE)

# save dataset
mat <- assay(rlog)
mm <- model.matrix(~disease, colData(rlog))
mat <- limma::removeBatchEffect(mat, batch=rlog$batch, design=mm)
write.table(mat, "/Users/leonivandijk/Desktop/thesis/data/syn3388564/normalized/275male_vst_highfiltered_prtncoding.csv",
            row.names = TRUE, col.names = TRUE, sep = "\t")
write.table(metadata, "/Users/leonivandijk/Desktop/thesis/data/syn3388564/normalized/pheno_275male.csv",
            row.names = TRUE, col.names = TRUE, sep = "\t")

# results differential expression
res05 <- results(dds, alpha=0.05)
summary(res05)

sum(res05$padj < 0.05, na.rm=TRUE)

plotCounts(dds, gene=which.min(res05$padj), intgroup="disease")

# ------------------------------ data visualizations -------------------------
# Perform PCA
pca_result <- prcomp(t(mat), scale. = FALSE, rank. = 5)
# Define custom colors
custom_colors <- c("control" = "#6BAED6", "alzheimer" = "#FFD700")
# Extract PCA data
pca_data <- data.frame(pca_result$x[, 1:5])
pca_data[] <- lapply(pca_data, function(x) x / sqrt(sum((x - mean(x))^2)))
pca_data$Sample <- rownames(pca_data)
pca_data$disease <- metadata$disease
# Calculate the outliers
outlier_indices <- which(pca_data$PC1^2 + pca_data$PC2^2 > quantile(pca_data$PC1^2 + pca_data$PC2^2, 0.95))
outliers <- pca_data[outlier_indices, ]

pca_plot <- autoplot(pca_result, data = metadata, colour = "disease", frame = TRUE, frame.type = 'norm', label.size = 4, 
                     x = 1, y = 2) +
  ggtitle("2D PCA-plot from AD Dataset: PC1 & PC2") +
  scale_colour_manual(values = custom_colors) +
  geom_text_repel(data = outliers, aes(x = PC1, y = PC2, label = Sample), size = 3) +
  scale_fill_manual(values = custom_colors) +
  geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed")+
  geom_point(size = 2, shape = 21, colour = "black")+
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01)) + scale_y_continuous(labels = scales::number_format(accuracy = 0.01))+
  theme_minimal() +
  theme(
    plot.title = element_text(size = 15, face = "bold"),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 15)
  ) +
  labs(x = paste0("Dim1 (", round(summary(pca_result)$importance[2,1] * 100, 1), "% of Variance)"),
       y = paste0("Dim2 (", round(summary(pca_result)$importance[2,2] * 100, 1), "% of Variance)"))

# Display the plot
print(pca_plot)
ggsave("/Users/leonivandijk/Desktop/thesis/data/syn3388564/PCA_1_2.png", plot = pca_plot)

# Calculate the outliers
outlier_indices <- which(pca_data$PC4^2 + pca_data$PC5^2 > quantile(pca_data$PC4^2 + pca_data$PC5^2, 0.95))
outliers <- pca_data[outlier_indices, ]

pca_plot <- autoplot(pca_result, data = metadata, colour = "disease", frame = TRUE, frame.type = 'norm', label.size = 4, 
                     x = 4, y = 5) +
  ggtitle("2D PCA-plot from AD Dataset: PC4 & PC5") +
  scale_colour_manual(values = custom_colors) +
  geom_text_repel(data = outliers, aes(x = PC4, y = PC5, label = Sample), size = 3) +
  scale_fill_manual(values = custom_colors) +
  geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed")+
  geom_point(size = 2, shape = 21, colour = "black")+
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01)) + scale_y_continuous(labels = scales::number_format(accuracy = 0.01))+
  theme_minimal() +
  theme(
    plot.title = element_text(size = 15, face = "bold"),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 15)
  ) +
  labs(x = paste0("Dim4 (", round(summary(pca_result)$importance[2,4] * 100, 1), "% of Variance)"),
       y = paste0("Dim5 (", round(summary(pca_result)$importance[2,5] * 100, 1), "% of Variance)"))

# Display the plot
print(pca_plot)
ggsave("/Users/leonivandijk/Desktop/thesis/data/syn3388564/PCA_3_4.png", plot = pca_plot)

# check for outlier samples
sampleTree = flashClust(dist(t(mat)), method = "average"); #do clustering
traitColors = numbers2colors(metadata$cogdx, signed = FALSE);
par(cex = 0.6);
par(mar = c(0,4,2,0))
pdf("/Users/leonivandijk/Desktop/thesis/data/syn3388564/sampledendo.pdf", width=30)
plotDendroAndColors(sampleTree, traitColors, groupLabels = "disease", rowText = NULL)
plot(sampleTree, main = "Sample clustering to detect outliers",  labels = metadata$cogdx, sub="", xlab="", cex.lab = 1.5,cex.axis = 1.5, cex.main = 2)
dev.off()

labels = cutreeStatic(sampleTree, cutHeight = 100,minSize=5)
labels <- as.data.frame(labels)
rownames(labels) <- rownames(metadata)
rownames(subset(labels, labels == 0))

library(pheatmap)
cors <- cor(mat, method="pearson")
pdf("/Users/leonivandijk/Desktop/thesis/data/syn3388564/cor_heatmap.pdf", width=12, height=12)
plot_pheatmap = pheatmap(cors)
dev.off()


sampleDists <- dist(t(assay(rlog)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- rlog$disease
colnames(sampleDistMatrix) <- rlog$disease
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

cors <- cor(mat, method="pearson")
pheatmap(cors)


