## Analysis after alignment


###---------Normalisation-------###
setwd("~/Folder1") #set working directory here to where appropriate 
library(SCnorm)
mouse = read.table('count_matrix.txt', header= F, quote="", comment.char="", sep="\t",
row.names = 1)
mouse = mouse[-1:-4,] # keep everything but the first to fourth columns as these are unmapped genes
conditions = rep(c(1), each=917) #set each cell as its own condition
mouse_normalised = SCnorm(Data=mouse, Conditions = conditions, PrintProgressPlots =
F, reportSF = F, FilterCellNum = 10, FilterExpression = 0,Thresh = 0.1, K=NULL, NCores =
NULL, ditherCounts = F, PropToUse = 0.25, Tau = 0.5, withinSample = NULL, useSpikes =
F, useZerosToScale = F) # kept default parameters


###--------------Select for insulin positive 2 cells-------------###
Ins2 <- normalised_mouse[grepl("ENSMUSG00000000215",
rownames(normalised_mouse)),] # find row name which contains the Ensembl ID for Ins2
gene
Ins2_sorted <- Ins2[,order(Ins2[1,])] # sort Ins2 expression values in ascending order
ins2_positive = Ins2_sorted[, -1:-307] # filter out ins2 cells with read count less than 500
ins2_pos_cells = colnames(ins2_positive) # find the column names for insulin positive cells
library(dplyr)
ins2_pos_all = select(normalised_mouse, one_of(ins2_pos_cells)) # select the insulin positive 2 columns from the normalised count matrix
transposed_ins2_pos = as.data.frame(t(ins2_pos_all)) # convert to dataframe and transpose
library(ggplot2)
ggplot(transposed_positive, aes(ENSMUSG00000029644.7)) +
 geom_histogram(fill="#43a2ca") +
  labs(title = "Frequency of pdx1 in Ins2 positive cells", x = "Cells", y = "Frequency") +
  theme(panel.background = element_blank())


###-------------PDX1 extraction from Insulin positive cells-----------####
PDX1 <- ins2_pos_all[grepl("ENSMUSG00000029644", rownames(ins2_pos_all)),] #select PDX1 gene from count matrix
PDX1_sorted <- PDX1[,order(PDX1[1,])] #order the cells from low to high expressing PDX1
all_PDX1_cells = colnames(PDX1_sorted) #select cells 
PDX1_all = select(normalised_mouse, one_of(all_PDX1_cells)) #now add other cells so it’s the count matrix but cells ordered based on PDX1 expression
PDX1_lowest_third <- PDX1_all[,1:203]
PDX1_mid_third <- PDX1_all[,204:407]
PDX1_top_third <- PDX1_all[408:610]
low_high_pdx1 = cbind(PDX1_lowest_third, PDX1_top_third) #count matrix with low and high PDX1 levels
write.table(low_high_pdx1, "~/Folder/low_high_pdx1.txt", sep="\t")



###------Clustering, PCA----------####
#Clustering
library(mclust)
ins2_positive_pca <- prcomp(ins2_pos_all, center = T, scale = T)
PC12 <- data.frame(PC1=ins2_positive_pca$x[,1],PC2=ins2_positive_pca$x[,2])
plot(PC12)
clust <- Mclust(PC12)
summary(clust)
plot(clust, what= "BIC")
all_metadata <- c(rep(c("lowest_third"), each=203), rep(c("middle_third"),
                                                        each=204),rep(c("top_third"), each=203))
cell_names <- colnames(PDX1_all)
metadata_matrix <- cbind(cell_names, all_metadata)
PDX1 <- c(rep(c("lowest"), each=203), rep(c("middle"), each=204),rep(c("highest"),
                                                                     each=203))
cell_names <- colnames(PDX1_all)
metadata_matrix <- cbind(PDX1, cell_names)
library(SingleCellExperiment)
library(SC3)
library(scater)
sce <- SingleCellExperiment(
  assays = list(
    counts = as.matrix(PDX1_all),
    logcounts = log2(as.matrix(PDX1_all) + 1)
  ),
  colData = metadata_matrix
)
plotPCA(sce, colour_by = "PDX1", shape_by = "PDX1")
#expression plot
sce <- SingleCellExperiment(
  assays = list(
    counts = as.matrix(PDX1_all),
    logcounts = log2(as.matrix(PDX1_all) + 1)
  ),
  colData = metadata_matrix
)
plotExpression(sce, rownames(sce)[16129], x="cell_names", exprs_values = "counts",
               colour_by = "PDX1")



###-------Differential Gene analysis-------###
rounded_all_PDX1 <- round(PDX1_all, digits=0) # round to one digit
PDX1_low_and_high <- rounded_all_PDX1[,-204:-407]
metadata_low_and_high <- c(rep(c("lowest_third"), each=203), rep(c("top_third"),
                                                                 each=203))
library(DESeq2)
cell_names <- colnames(PDX1_low_and_high)
metadata_matrix_low_and_high <- cbind(cell_names, metadata_low_and_high)
dds <- DESeqDataSetFromMatrix(countData = PDX1_low_and_high, colData =
                                metadata_matrix_low_and_high, design = ~ metadata_low_and_high)
ones <- c(rep(1, 406))
sizeFactors(dds) <- ones # divide by scale factor one as data is already normalised
dds2 <- DESeq(dds)
BM <- sapply( levels(dds2$metadata_low_and_high), function(lvl) rowMeans(
  counts(dds2,normalized=TRUE)[,dds2$metadata_low_and_high == lvl] ) )
BM <- as.data.frame(BM)
low_vs_high_results <- results(dds2, name =
                                 "metadata_low_and_high_top_third_vs_lowest_third")
low_vs_high_results2 <- na.omit(low_vs_high_results) # get rid of missing data
low_v_high_significant <- low_vs_high_results2[low_vs_high_results2$padj<0.05,] # select for signficant DEGs with p value <0.05
write.table(low_vs_high_results, "DEG_low_v_high_DESeq_nofilter2", append = FALSE,
            sep = "\t", dec = ".", row.names = TRUE, col.names = FALSE)
write.table(low_v_high_significant, "DEseq_low_vs_high_significant", append = FALSE, sep
            = "\t", dec = ".", row.names = TRUE, col.names = TRUE)

#change name for all names after deseq
DEG_low_v_high_DESeq_nofilter <- read.delim('DEG_low_v_high_DESeq_nofilter2',
                                            header = F)
library('biomaRt')
mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
DEG_low_v_high_DESeq_nofilter$V1 <- gsub("\\.[0-9]*$", "",
                                         DEG_low_v_high_DESeq_nofilter$V1) #removes decimals in the names
Genes <- DEG_low_v_high_DESeq_nofilter$V1
G_list <- getBM(filters= "ensembl_gene_id",
                attributes= c("ensembl_gene_id" , "mgi_symbol"),
                values= Genes,
                mart= mart)
write.table(G_list, "~/Folder/G_list.txt", sep="\t")
mouse_genes = rep(NA, nrow(DEG_low_v_high_DESeq_nofilter))
DEG_low_v_high_DESeq_nofilter = cbind(mouse_genes,
                                      DEG_low_v_high_DESeq_nofilter)
DEG_low_v_high_DESeq_nofilter$mouse_genes <-
  G_list$mgi_symbol[match(DEG_low_v_high_DESeq_nofilter$V1, G_list$ensembl_gene_id)]
write.table(DEG_low_v_high_DESeq_nofilter,
            "~/Folder/DEG_low_v_high_DESeq_nofilter_named.txt", sep="\t")



###---------Plot log fold change for heatmap---------###
all_genes <- read.delim("DEG_low_v_high_DESeq_nofilter_named.txt", row.names = 1)
all_genes <- all_genes[, -1:-2] #remove the gene names in the data frame
all_genes2 <- na.omit(all_genes) #rows which contains na's
all_genes3 <- subset(all_genes2, padj<0.05 ) #select significant genes
all_genes3 <- all_genes3[order(all_genes3$padj),] #order the significant genes from most to
least significant
ordered_log <- all_genes3[order(-all_genes3$log2foldchange), ] #order log change from pos
to neg lfc
top_20_pos_lfc <- ordered_log[1:50, ] #top 20 highest log fold change
top_20_gene_name <- as.vector(top_20_pos_lfc$ensemble_ID) #take ensemble ID/names
of these genes
low_high_pdx1 <- read.delim("low_high_pdx1.txt") #load in the low and high pdx1 count
matrix
G_list <- read.delim("G_list.txt", row.names = 1)#load in the conversion table for ensemble
ID to gene names
top_20_gene_count = low_high_pdx1[c(top_20_gene_name),]#select the genes which have
the highest lfc
mouse_gene_names = row.names(top_20_gene_count) #filter out the gene ensemble IDs
mouse_gene_names2 = gsub("\\.[0-9]*$", "", mouse_gene_names) #remove the dots or any
other symbol to allow name conversion
mouse_gene_symbol = rep(NA, nrow(top_20_gene_count)) # make an empty vector for the
real names to be filled in
mouse_gene_names3 = as.data.frame(cbind(mouse_gene_names2, mouse_gene_symbol))
#bind the empty vector and the matrix
mouse_gene_names3$mouse_gene_symbol <-
  G_list$mgi_symbol[match(mouse_gene_names3$mouse_gene_names2,
                          G_list$ensembl_gene_id)] #match the ensemble ID and the mouse symbols using the G_list
conversion table
rownames(top_20_gene_count) = mouse_gene_names3$mouse_gene_symbol #replace the
row names with the real gene names
log_top_20_gene_count2 = sign(top_20_gene_count) * log(abs(top_20_gene_count))
#absolute log transform the counts so they are relative to each other
library(pheatmap)
pheatmap(log_top_20_gene_count2, main="Log Fold Change", sep = "",
         color = colorRampPalette(c("white", "turquoise4","black"))(10),
         cluster_cols = FALSE, cluster_rows = FALSE, show_colnames = F)



###---------plot padj heatmap--------###
all_genes <- read.delim("DEG_low_v_high_DESeq_nofilter_named.txt", row.names = 1)
all_genes <- all_genes[, -1:-2] #remove gene names
all_genes2 <- na.omit(all_genes) #remove rows containing nas from data frame
ordered_padj <- all_genes2[order(all_genes2$padj), ] #order the padj from lowest to highest
top_20_pos_padj <- ordered_padj[1:20, ] #top 20 highest adj
top_20_gene_name <- as.vector(top_20_pos_padj$ensemble_ID) #take ensemble ID of
these
low_high_pdx1 <- read.delim("low_high_pdx1.txt") #load in the low and high pdx1 count
matrix
G_list <- read.delim("G_list.txt", row.names = 1)#load in the conversion table for ensemble
ID to gene names
top_20_gene_count = low_high_pdx1[c(top_20_gene_name),]#select the genes which have
the highest lfc
mouse_gene_names = row.names(top_20_gene_count) #filter out the gene ensemble IDs
mouse_gene_names2 = gsub("\\.[0-9]*$", "", mouse_gene_names) #remove the dots or any
other symbol to allow name conversion
mouse_gene_symbol = rep(NA, nrow(top_20_gene_count)) # make an empty vector for the
real names to be filled in
mouse_gene_names3 = as.data.frame(cbind(mouse_gene_names2, mouse_gene_symbol))
#bind the empty vector and the matrix
mouse_gene_names3$mouse_gene_symbol <-
  G_list$mgi_symbol[match(mouse_gene_names3$mouse_gene_names2,
                          G_list$ensembl_gene_id)] #match the ensemble ID and the mouse symbols using the G_list
conversion table
rownames(top_20_gene_count) = mouse_gene_names3$mouse_gene_symbol #replace the
row names with the real gene names
log_top_20_gene_count2 = sign(top_20_gene_count) * log(abs(top_20_gene_count))


#absolute log transform the counts so they are relative to each other
library(pheatmap)
pheatmap(log_top_20_gene_count2, main="Padj", sep = "",
         color = colorRampPalette(c("white", "turquoise4","black"))(10),
         cluster_cols = FALSE, cluster_rows = FALSE, show_colnames = F)


###---------volcano plot ggplot2----------###
library("ggplot2")
library("ggrepel") #Avoid overlapping labels
library(dplyr)
mutateddf <- mutate(all_genes2, sig=ifelse(all_genes2$padj<0.05, "Significant genes", "Not
Significant Genes")) #Will have different colors depending on significance
input <- cbind(gene=rownames(mutateddf ), mutateddf ) #convert the rownames to a column
volc = ggplot(input, aes(log2foldchange, -log10(padj))) +
  geom_point(aes(col=sig)) + #add points colored by significance
  scale_color_manual(values=c("black", "red")) +
  ggtitle("Volcano plot") + #e.g. 'Volcanoplot DESeq2'
  ylim(0, 40) +
  xlim(-5, 5) +
  theme(plot.title = element_text(hjust = 0.5)) #volcanoplot with log2Foldchange versus
pvalue
volc+geom_text_repel(data = head(input, 0), aes(label= ensemble_ID)) #option can include
adding name text for the top 20 genes in the plot but we decided not to



###-----Tables for DAVID analysis---###
named_significant_high_low <- read.delim("named_high_low_significant.txt")
low_v_high_upregulated <-
  named_significant_high_low[named_significant_high_low$log2FoldChange>0,]
low_v_high_downregulated <-
  named_significant_high_low[named_significant_high_low$log2FoldChange<0,]




###-----Enrichment plots---------###
#Enrichment scores -- upregulated (Functional Annotation)
upregulated_genes <- read.delim('Low_v_high_upregulated_functional_annotation.txt',
                                header = F)
library(dplyr)
x <- c("Annotation Cluster 1", "Annotation Cluster 2",
       "Annotation Cluster 3", "Annotation Cluster 4", "Annotation Cluster 5") #Extract rows that
have enrichment scores and cluster.
upregulated_genes_metadat = matrix(nrow = 5, ncol = 2)
upregulated_genes_metadat[1:5] = x
upregulated_genes_metadat = as.data.frame(upregulated_genes_metadat)
Enrichment_Score = rep(NA, nrow(upregulated_genes_metadat))
upregulated_genes_metadat3 = cbind(upregulated_genes_metadat, Enrichment_Score)
upregulated_genes_metadat3$Enrichment_Score =
  upregulated_genes$V2[match(upregulated_genes_metadat3$V1, upregulated_genes$V1)]
upregulated_genes_metadat4 = upregulated_genes_metadat3[, -2]
upregulated_genes_metadat4$Enrichment_Score = gsub("Enrichment Score: ", "",
                                                   upregulated_genes_metadat4$Enrichment_Score)
row.names(upregulated_genes_metadat4) = upregulated_genes_metadat4$V1
Clusters = c("Splicing", "Nucleoside biosynthetic process", "Regulation of peptide hormone
secretion",
             "Nicotinamide nucleotide metabolic process", "Positive regulation of RNA
biosynthetic process")
upregulated_genes_metadat5 = cbind(Clusters, upregulated_genes_metadat4)
write.table(upregulated_genes_metadat4, "~/Folder/upregulated_genes_metadat4.txt",
            sep="\t")
counts <- as.numeric(upregulated_genes_metadat4$Enrichment_Score)
dev.new(width = 25, height = 4)
par(mar= c(5, 16, 5, 1))
barplot(counts, main= "Functional annotation on upregulated genes", horiz=TRUE,
        names.arg= Clusters, cex.names=0.8, las=2, xlim=c(0,10), col = "Turquoise", xlab =
          "Enrichment Score")


#Enrichment plot -- downregulated (Functional Annotation)
downregulated_genes <- read.delim('annotation_downregulated.txt', header = F)
library(dplyr)
x <- c("Annotation Cluster 1", "Annotation Cluster 2",
       "Annotation Cluster 3", "Annotation Cluster 4", "Annotation Cluster 5") #Extract rows that
have enrichment scores and cluster.
downregulated_genes_metadat = matrix(nrow = 5, ncol = 2)
downregulated_genes_metadat[1:5] = x
downregulated_genes_metadat = as.data.frame(downregulated_genes_metadat)
Enrichment_Score = rep(NA, nrow(downregulated_genes_metadat))
downregulated_genes_metadat3 = cbind(downregulated_genes_metadat,
                                     Enrichment_Score)
downregulated_genes_metadat3$Enrichment_Score =
  downregulated_genes$V2[match(downregulated_genes_metadat3$V1,
                               downregulated_genes$V1)]
downregulated_genes_metadat4 = downregulated_genes_metadat3[, -2]
downregulated_genes_metadat4$Enrichment_Score = gsub("Enrichment Score: ", "",
                                                     downregulated_genes_metadat4$Enrichment_Score)
row.names(downregulated_genes_metadat4) = downregulated_genes_metadat4$V1
Clusters = c("Regulation of protein
 kinase activity", "Synaptic signaling", "Positive regulation of
protein kinase activity",
             "Positive regulation of
 protein phosphorylation", "Focal adhesion
 assembly")
downregulated_genes_metadat5 = cbind(Clusters, downregulated_genes_metadat4)
write.table(upregulated_genes_metadat4, "~/Module 6/upregulated_genes_metadat4.txt",
            sep="\t")
counts <- as.numeric(downregulated_genes_metadat4$Enrichment_Score)
dev.new(width = 25, height = 4)
par(mar= c(5, 16, 5, 1))
barplot(counts, main="Functional annotation on downregulated genes", horiz=TRUE,
        names.arg= Clusters, cex.names=1.2, cex.main=1.5, las=2, xlim=c(0,6), col =
          "#c994c7", xlab="Enrichment score" )


#Enrichment plot -- upregulated (Functiona; Classification)
upregulated_genes <- read.delim('classification_upregulated.txt', header = F)
library(dplyr)
x <- c("Gene Group 1", "Gene Group 2",
       "Gene Group 3","Gene Group 4", "Gene Group 5") #Extract rows that have enrichment
scores and cluster.
upregulated_genes_metadat = matrix(nrow = 3, ncol = 2)
upregulated_genes_metadat[1:3] = x
upregulated_genes_metadat = as.data.frame(upregulated_genes_metadat)
Enrichment_Score = rep(NA, nrow(upregulated_genes_metadat))
upregulated_genes_metadat3 = cbind(upregulated_genes_metadat, Enrichment_Score)
upregulated_genes_metadat3$Enrichment_Score =
  upregulated_genes$V2[match(upregulated_genes_metadat3$V1, upregulated_genes$V1)]
upregulated_genes_metadat4 = upregulated_genes_metadat3[, -2]
upregulated_genes_metadat4$Enrichment_Score = gsub("Enrichment Score: ", "",
                                                   upregulated_genes_metadat4$Enrichment_Score)
row.names(upregulated_genes_metadat4) = upregulated_genes_metadat4$V1
counts <- as.numeric(upregulated_genes_metadat4$Enrichment_Score)
dev.new(width = 25, height = 4)
par(mar= c(5, 8, 5, 1))
barplot(counts, main="Functional classification on upregulated genes", horiz=TRUE,
        names.arg= x, cex.names=1.2, cex.main=1.5, las=2, xlim=c(0,4), col = "#c994c7", xlab
        = "Enrichment score" )


#Enrichment plot -- downregulated (Functional Classification)
downregulated_genes <- read.delim('classification_downregulated.txt', header = F)
library(dplyr)
x <- c("Gene Group 1", "Gene Group 2",
       "Gene Group 3") #Extract rows that have enrichment scores and cluster.
downregulated_genes_metadat = matrix(nrow = 3, ncol = 2)
downregulated_genes_metadat[1:3] = x
downregulated_genes_metadat = as.data.frame(downregulated_genes_metadat)
Enrichment_Score = rep(NA, nrow(downregulated_genes_metadat))
downregulated_genes_metadat3 = cbind(downregulated_genes_metadat,
                                     Enrichment_Score)
downregulated_genes_metadat3$Enrichment_Score =
  downregulated_genes$V2[match(downregulated_genes_metadat3$V1,
                               downregulated_genes$V1)]
downregulated_genes_metadat4 = downregulated_genes_metadat3[, -2]
downregulated_genes_metadat4$Enrichment_Score = gsub("Enrichment Score: ", "",
                                                     downregulated_genes_metadat4$Enrichment_Score)
row.names(downregulated_genes_metadat4) = downregulated_genes_metadat4$V1
counts <- as.numeric(downregulated_genes_metadat4$Enrichment_Score)
dev.new(width = 25, height = 4)
par(mar= c(5, 8, 5, 1))
barplot(counts, main="Functional classification on downregulated genes", horiz=TRUE,
        names.arg= x, cex.names=1.2, cex.main=1.5, las=2, xlim=c(0,4), col = "#c994c7", xlab
        = "Enrichment score" )