# --------------------------------  preparation --------------------------------
library(DESeq2)
library(apeglm)
library(pheatmap)
library(biomaRt)  
library(WGCNA)
library(ggrepel)
library(gridExtra)
library(ggfortify)
setwd("/Users/leonivandijk/Desktop/thesis/")
urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
path <- paste(urld, "acc=GSE64810", "file=GSE64810_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep="&");

# read files
raw_counts <- as.matrix(data.table::fread(path, header=T, colClasses="integer"), rownames=1)
metadata   <- read.table(
              "rfiles/HD_sample_info_as_numeric__only_condition_death.txt", header = TRUE,
              sep = "\t", row.names = 1, check.names = FALSE, dec = ".")
genes_ad <-
  rownames(read.table(
    "data/syn3388564/normalized/275male_vst_highfiltered_prtncoding.csv", header = TRUE,
    sep = "\t", check.names = FALSE, stringsAsFactors = FALSE, row.names = 1,
    as.is = TRUE , dec = "."))
genes_ad <- sub("\\..*", "", genes_ad)

# prepare metadata
metadata$disease <- ifelse(metadata$Condition == 1, "huntington", 
                         ifelse(metadata$Condition == 0, "control", NA))
metadata$disease <- factor(metadata$disease, levels = c("control","huntington"))


# prepare expression data
#    - keep only protein_coding RNA (switch to ensembl ids)
#    - pre-filter out low count genes: keep genes with at least 10 counts > 25 samples
#    - build DESeqDataSet

# pre-filter low count genes
raw_counts <- raw_counts[rowSums(raw_counts >= 10) >= 25, ]

ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")  
entrezgene = rownames(raw_counts)
genes <- getBM(filters="entrezgene_id", attributes=c("ensembl_gene_id","entrezgene_id"), values=entrezgene, mart=ensembl)

entrez_ids = data.frame(rownames(raw_counts))
entrez_ids['ensembl_mapping'] = NaN
rownames(entrez_ids) = entrez_ids$rownames.raw_counts.

j = 0
for (id in entrez_ids$rownames.raw_counts.){
  j = j+1
  if (j%%1000 == 0){
    print(j)
  }
  mapping = subset(genes, genes$entrezgene_id==id)
  # no ensemble id mapping found
  if (length(rownames(mapping)) == 0){
    entrez_ids[id,2] = "no mapping"
  }
  # unique ensembl id mapping found
  else if (length(rownames(mapping)) == 1){
    entrez_ids[id,2] = mapping[1,1]
  }
  # if we find multiple mappings, we try to align with AD dataset
  else{
    found_alignment = 0 
    for (i in c(1:length(mapping[,1]))){
      if (mapping[i,1] %in% genes_ad){
        entrez_ids[id,2] = mapping[i,1]
        found_alignment = found_alignment+1
      }
    }
    # no alignment found just take a random one
    if (found_alignment == 0) {
      i = sample(1:length(mapping[,1]),1)
      entrez_ids[id,2] = mapping[i,1]
    } 
  }
  }


# remove entrez ids without mapping
entrez_ids = subset(entrez_ids, entrez_ids$ensembl_mapping != "no mapping")
# remove duplications in emsembl mapping
doublemapping = subset(entrez_ids, duplicated(entrez_ids['ensembl_mapping']))
entrez_ids = subset(entrez_ids, !(entrez_ids$rownames.raw_counts. %in% doublemapping$rownames.raw_counts.))

# keep only protein coding RNA
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
annot <- getBM(attributes=c("ensembl_gene_id", "gene_biotype"),
               values = entrez_ids$ensembl_mapping,
               mart = ensembl)
rownames(annot) <- annot$ensembl_gene_id
annot <- annot[entrez_ids$ensembl_mapping,]
annot <- subset(annot, gene_biotype == 'protein_coding')

raw_counts <- raw_counts[rownames(entrez_ids), ]
rownames(raw_counts) <- entrez_ids$ensembl_mapping
raw_counts <- raw_counts[rownames(annot), ]

# remove outliers (based on inspection)
colnames(raw_counts) <- rownames(metadata)
metadata <- subset(metadata, !(rownames(metadata) %in% c('C_0018', 'C_0009')))
raw_counts <- raw_counts[,rownames(metadata)]

#   - print sample groups
length(rownames(subset(metadata, disease=='huntington')))
length(rownames(subset(metadata, disease=='control')))

dds <- DESeqDataSetFromMatrix(
  countData = raw_counts,
  colData = metadata,
  design = ~disease
)

rm(ensembl, annot, doublemapping, entrez_ids, mapping, genes)

# ----------------------  differential expression analysis ---------------------
#   - operate on raw counts
#   - for later clustering (WGCNA) we will add additional transformations
dds <- DESeq(dds) # original counts preserved in counts(dds)
norm_counts <- counts(dds, normalized=TRUE)
norm_counts <- norm_counts[goodGenes(t(norm_counts)),]
rlog <- vst(dds, blind=FALSE) #actually do rlog later at night

# save dataset
mat <- assay(rlog)
mm <- model.matrix(~disease, colData(rlog))

write.table(mat, "/Users/leonivandijk/Desktop/thesis/data/GSE64810/normalized/67_vst_highfiltered_prtncoding.csv",
            row.names = TRUE, col.names = TRUE, sep = "\t"
)
write.table(metadata, "/Users/leonivandijk/Desktop/thesis/data/GSE64810/normalized/pheno_67male.csv",
            row.names = TRUE, col.names = TRUE, sep = "\t"
)


res <- results(dds)

resLFC <- lfcShrink(dds, coef="disease_huntington_vs_control", type="apeglm")
resLFC

resOrdered <- res[order(res$pvalue),]

res05 <- results(dds, alpha=0.05)
summary(res05)

sum(res$padj < 0.05, na.rm=TRUE)
plotMA(resLFC, ylim=c(-2,2))
plotCounts(dds, gene=which.min(res$padj), intgroup="disease")


# Perform PCA
pca_result <- prcomp(t(mat), scale. = FALSE, rank. = 5)
# Define custom colors
custom_colors <- c("control" = "#6BAED6", "huntington" = "#FFD700")
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
  ggtitle("2D PCA-plot from HD Dataset: PC1 & PC2") +
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
ggsave("/Users/leonivandijk/Desktop/thesis/data/GSE64810/PCA_3_4.png", plot = pca_plot)


# sample clustering
sampleTree = flashClust(dist(t(mat)), method = "average"); #do clustering
traitColors = numbers2colors(metadata$Condition, signed = FALSE);
par(cex = 0.6);
par(mar = c(0,4,2,0))
#pdf("/Users/leonivandijk/Desktop/thesis/data/GSE64810/dendo_sample.pdf", width=15)
plot_Dendro = plotDendroAndColors(sampleTree, traitColors, groupLabels = "disease", rowText = NULL)
dev.off()

# correlation between samples
library(pheatmap)
cors <- cor(mat, method="pearson")
#pdf("/Users/leonivandijk/Desktop/thesis/data/GSE64810/cor_heatmap.pdf", width=12, height=12)
plot_pheatmap = pheatmap(cors)
dev.off()


#combine figures

# Create individual plots
p1 <- PCA12

p2 <- PCA34

p3 <- grid::grid.grabExpr(plot_Dendro)

p4 <- grid::grid.grabExpr(plot_pheatmap)
# Arrange plots in a 2x2 grid
grid.arrange(p1, p2, p3, p4, ncol=2)
