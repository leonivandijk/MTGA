# --------------------------------  preparation --------------------------------
library("WGCNA")
library(flashClust)
options(digits = 16)
sizeGrWindow(12,9)
setwd("/Users/leonivandijk/Desktop/thesis/")
output_path = 'rfiles/output/july/HD'

# read files
phenoHD <-
  read.table(
    "data/GSE64810/normalized/pheno_67male.csv", header = TRUE,
    sep = "\t", row.names = 1, check.names = FALSE, dec = ".")

genecounts_HD <-
  read.table(
    "data/GSE64810/normalized/67_vst_highfiltered_prtncoding.csv", header = TRUE,
    sep = "\t", check.names = FALSE, stringsAsFactors = FALSE, row.names = 1,
    as.is = TRUE)

pheno <-
  read.table(
    "data/syn3388564/normalized/pheno_275male.csv", header = TRUE,
    sep = "\t", row.names = 1, check.names = FALSE, dec = ".")
pheno <- pheno[c('disease', 'braaksc', 'ceradsc')]
pheno$disease <- ifelse((pheno$disease == "control"),0, 
                                         ifelse((pheno$disease == "alzheimer"), 1, NA))
# outliers
pheno <- subset(pheno, !(rownames(pheno) %in% c('380_120503', '367_120502', '500_120515', '507_120515')))

genecounts <-
  read.table(
    "data/syn3388564/normalized/275male_vst_highfiltered_prtncoding.csv", header = TRUE,
    sep = "\t", check.names = FALSE, stringsAsFactors = FALSE, row.names = 1,
    as.is = TRUE , dec = ".")
rownames(genecounts) <- sub("\\..*", "", rownames(genecounts))
genecounts <- t(genecounts)
genecounts <-genecounts[rownames(pheno),]

#if we run for hd
pheno <- phenoHD
genecounts <- t(genecounts_HD)


# ----------------------- WGCNA: determine soft power --------------------------
powers <- c(c(1:10), seq(from = 12, to=30, by=2))
cex1 <- 0.9
opath_powers = paste(output_path, "powers", sep="/")


# compute network topology for the different power options and save results
sft <- pickSoftThreshold(genecounts, powerVector = powers, verbose = 5,
                         networkType="signed hybrid", corFnc = "bicor", corOptions = list(use = "p"))

pdf(paste(opath_powers, "sft.pdf", sep="/"))

# scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1],
     -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,unsigned R^2",type="n",
     main = paste("Scale independence"))

text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")

abline(h=0.80,col="red")

# mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))

text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

write.table(sft$fitIndices,
            paste(opath_powers,"sft_table_signed_hybrid.txt", sep="/"),
            col.names=T, row.names=F,sep="\t")


# ------------------------ WGCNA: calculate matrices ---------------------------
#clear space in environment
rm(sft, cex1, powers)

soft_power <- 10

# compute adjacancy matrix based on expression correlations
adjacency <-
  adjacency(genecounts, power = soft_power, type = "signed hybrid",
            corFnc = "bicor", corOptions = list(use = "p"))

# turn adjacency matrix into topological overlap
TOM <- TOMsimilarity(adjacency, TOMType = "signed")
dissTOM <- 1-TOM



# -------------------------- WGCNA: detect modules -----------------------------
# clear space in environment
rm(adjacency, TOM)

min_module_size <- 25
opath_modules = paste(output_path, "modules", sep="/")


# cluster the genes according to the dissimilarity matrix and plot the dendogram
gene_tree <- flashClust(as.dist(dissTOM), method = "average")

# identify modules using dynamic tree cut
dynamic_modules <-
  cutreeDynamic(dendro = gene_tree, distM = dissTOM, deepSplit = 4,
                pamRespectsDendro = FALSE, minClusterSize = min_module_size)

# covert numeric labels into colors and plot the dendogram
colors <- labels2colors(dynamic_modules)
pdf(paste(opath_modules, "gene_tree_colors_filtered_hyb10_dp4_size25.pdf", sep = "/"), width=10)
plotDendroAndColors(gene_tree, cbind(colors), 
                    cbind("Dynamic Tree Cut"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()


# -------------------------- WGCNA: refine modules -----------------------------
MEList <- moduleEigengenes(genecounts, colors = colors)

# calculate dissimilarity of module eigengenes
MEs <- MEList$eigengenes
MEDiss <- 1-cor(MEs)

# cluster module eigengenes that represent similar expression and plot modules
ME_diss_thres <- 0.1
METree <- flashClust(as.dist(MEDiss), method = "average")

# merge similar modules
merge <- mergeCloseModules(genecounts, colors, cutHeight = ME_diss_thres,
                           verbose = 3)
merged_colors <- merge$colors

# compute new eigengenes and plot the new modules
merged_MEs <- merge$newMEs
pdf(paste(opath_modules, "mergME_AVGclust_ME_hyb8_thres01_size25_dspl4.pdf", sep = "/"), width=15)
plotDendroAndColors(gene_tree, cbind(colors, merged_colors),
                    c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels = FALSE,
                    hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()

# overwrite previous modules
moduleColors <- merged_colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = merged_MEs


# ---- relating modules to clinical traits and identifying important genes -----
n_genes <- ncol(genecounts)
n_samples <- nrow(genecounts)
opath_traits = opath_modules
MEs = orderMEs(moduleEigengenes(genecounts, moduleColors)$eigengenes)

module_trait_correlation <- cor(MEs, pheno$Condition, use = "p")
module_trait_pvalue = corPvalueStudent(module_trait_correlation, n_samples)
write.table(module_trait_correlation,
            paste(opath_traits, "module_trait_correlation.csv", sep = "/"),
            col.names = TRUE, row.names = TRUE, sep = "\t")
write.table(module_trait_pvalue,
            paste(opath_traits, "module_trait_pvalue.csv", sep = "/"),
            col.names = TRUE, row.names = TRUE, sep = "\t")

# display p-values in a heat-map
pdf(paste(opath_traits, "modHeatmap_pvalue.pdf", sep = "/"))
textMatrix = paste(signif(module_trait_correlation, 2), "\n(",
                   signif(module_trait_pvalue, 1), ")", sep = "")
dim(textMatrix) = dim(module_trait_correlation)
par(mar = c(6, 8.5, 3, 3))

labeledHeatmap(Matrix = module_trait_pvalue,
               xLabels = names(pheno)[1],
               yLabels = names(MEs),
               ySymbols = names(MEs),
               cex.lab.y = 0.3,
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.15,
               zlim = c(0,0.1),
               main = paste("Module-trait relationships"))
dev.off()

# display cor-values in a heat-map
pdf(paste(opath_traits, "modHeatmap_corvalue.pdf", sep = "/"))
textMatrix = paste(signif(module_trait_correlation, 2), " (",
                   signif(module_trait_pvalue, 1), ")", sep = "")
dim(textMatrix) = dim(module_trait_correlation)
par(mar = c(6, 8.5, 3, 3))

labeledHeatmap(Matrix = module_trait_correlation,
               xLabels = names(pheno)[1],
               yLabels = names(MEs),
               ySymbols = names(MEs),
               cex.lab.y = 0.3,
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.3,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

module_counts <- data.frame(table(moduleColors))
length(rownames(module_counts))
length(module_trait_pvalue[,1][abs(module_trait_pvalue[,1]) < 0.05])
mean(module_counts$Freq)
median(module_counts$Freq)
max(abs(module_trait_correlation[,1]))
mean(abs(module_trait_correlation[,1]))
median(abs(module_trait_correlation[,1]))

moduleColors <- data.frame(moduleColors)
rownames(moduleColors) <- colnames(genecounts)
for (color in module_counts$moduleColors){
  if (module_trait_pvalue[paste("ME",color, sep=""), "disease"] < 0.05){
    genes_module <- rownames(moduleColors)[moduleColors == color]
    write.table(genes_module, paste(opath_traits, paste(color, 'csv', sep='.'), sep = "/"),
                row.names = TRUE, col.names = TRUE, sep = "\t"
    )}
}
  




