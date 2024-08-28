library("ggVennDiagram")
library("ggvenn")
# this script is used to analyse commonalities and differences of the differential gene expression analysis of AD and HD data.
# saving the DEGs is required before performing this analysis.


# load data
setwd("/Users/leonivandijk/Desktop/thesis/")

hd <-
  read.table(
    "data/GSE64810/normalized/67_vst_highfiltered_prtncoding.csv", header = TRUE,
    sep = "\t", check.names = FALSE, stringsAsFactors = FALSE, row.names = 1,
    as.is = TRUE)

hd_deseq <- 
  read.table(
    "pyfiles/MCGA/data/result_DESeq15989.csv", header = TRUE,
    sep = "\t", check.names = FALSE, stringsAsFactors = FALSE, row.names = 1,
    as.is = TRUE)
hd_deseq <- subset(hd_deseq, hd_deseq$padj < 0.05)
hd_up <- subset(hd_deseq, hd_deseq$log2FoldChange > 0)
hd_down <- subset(hd_deseq, hd_deseq$log2FoldChange < 0)

pheno_ad <-
  read.table(
    "data/syn3388564/normalized/pheno_275male.csv", header = TRUE,
    sep = "\t", row.names = 1, check.names = FALSE, dec = ".")
pheno_ad <- subset(pheno_ad, !(rownames(pheno_ad) %in% c('380_120503', '367_120502', '500_120515', '507_120515')))

ad <-
  read.table(
    "data/syn3388564/normalized/275male_vst_highfiltered_prtncoding.csv", header = TRUE,
    sep = "\t", check.names = FALSE, stringsAsFactors = FALSE, row.names = 1,
    as.is = TRUE , dec = ".")
rownames(ad) <- sub("\\..*", "", rownames(ad))
ad <- t(ad)
ad <-ad[rownames(pheno_ad),]
ad <- data.frame(t(ad))
rm(pheno_ad)

ad_deseq <- 
  read.table(
    "pyfiles/MCGA/data/result_DESeq15380.csv", header = TRUE,
    sep = "\t", check.names = FALSE, stringsAsFactors = FALSE, row.names = 1,
    as.is = TRUE)
ad_deseq <- subset(ad_deseq, ad_deseq$padj < 0.05)
ad_up <- subset(ad_deseq, ad_deseq$log2FoldChange > 0)
ad_down <- subset(ad_deseq, ad_deseq$log2FoldChange < 0)

# make lists of the genes
gene_lists <- list()
gene_lists[["Genes in AD dataset"]] <- rownames(ad)
gene_lists[["AD DEGs"]] <- rownames(ad_deseq)
gene_lists[["Genes in HD dataset"]] <- rownames(hd)
gene_lists[["HD DEGs"]] <- rownames(hd_deseq)


gene_lists_1 <-list()
gene_lists_1[["Genes in AD dataset"]] <- rownames(ad)
gene_lists_1[["Genes in HD dataset"]] <- rownames(hd)

gene_lists_2 <-list()
gene_lists_2[["AD DEGs"]] <- rownames(ad_deseq)
gene_lists_2[["HD DEGs"]] <- rownames(hd_deseq)

updown_list <- list()
updown_list[["AD UP"]] <- rownames(ad_up)
updown_list[["AD DOWN"]] <- rownames(ad_down)
updown_list[["HD UP"]] <- rownames(hd_up)
updown_list[["HD DOWN"]] <- rownames(hd_down)



# make plot
p <- ggvenn(
  gene_lists_1, 
  fill_color = c("#ACECF7", "#F1EAC8"),
  stroke_size = 0.5, set_name_size = 5, text_size=5.5
)
ggsave("data/venn_diagram_genes.png", plot = p, width = 5, height = 5, units = "in", dpi = 300)

p <- ggvenn(
  gene_lists_2, 
  fill_color = c("#ACECF8", "#F1EAC9"),
  stroke_size = 0.5, set_name_size = 5, text_size=5.5
)
ggsave("data/venn_diagram_degs.png", plot = p, width = 5, height = 5, units = "in", dpi = 300)



data <- data.frame(
  Category = c("up-up", 
               "down-down", 
               "up(AD) - down(HD)", 
               "down(AD) - up(HD)"),
  Count = c(749, 697, 74, 48)
)

custom_colors <- c("up-up" = "#ACECF7",
                   "down-down" = "#F1EAC8",
                   "up(AD) - down(HD)" = "#6BAED6",
                   "down(AD) - up(HD)" = "#F3D694")

# Create the bar chart
 p <- ggplot(data, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity", color = "black", size = 0.2) +  # Add border to bars
  theme_minimal() +
  ggtitle("Differentially Expressed Genes in AD and HD") +
  scale_fill_manual(values = custom_colors) +
  theme(
    plot.title = element_text(size = 15),
    axis.title = element_text(size = 10),
    axis.text.x = element_blank(),  # Remove x-axis text
    axis.ticks.x = element_blank(), # Remove x-axis ticks
    axis.text.y = element_text(size = 10),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 15),
  ) +
  labs(x = "Category", y = "Number of Genes")

ggsave("data/updown-degs.png", plot = p, width = 5, height = 5, units = "in", dpi = 300)


# export background list for genes
combined_genes <- c(gene_lists$`Genes in AD dataset`, gene_lists$`Genes in HD dataset`)
background <- unique(combined_genes)
write.table(background, file = "data/background.txt", col.names = FALSE, row.names=FALSE)

upup <- c(updown_list$`AD UP`, updown_list$`HD UP`)
upup <- upup[duplicated(upup)]
write.table(upup, file = "data/upup", col.names = FALSE, row.names=FALSE)

downdown <- c(updown_list$`AD DOWN`, updown_list$`HD DOWN`)
downdown <- downdown[duplicated(downdown)]
write.table(downdown, file = "data/downdown", col.names = FALSE, row.names=FALSE)

updown <- c(updown_list$`AD UP`, updown_list$`HD DOWN`)
updown <- updown[duplicated(updown)]
write.table(updown, file = "data/updown", col.names = FALSE, row.names=FALSE)

downup <- c(updown_list$`AD DOWN`, updown_list$`HD UP`)
downup <- downup[duplicated(downup)]
write.table(downup, file = "data/downup", col.names = FALSE, row.names=FALSE)

