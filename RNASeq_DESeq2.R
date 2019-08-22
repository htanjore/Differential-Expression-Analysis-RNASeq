setwd("~/RNAseq/AnalysisR")


library(tidyverse)
library(DESeq2)
library(ggplot2)
#library(dbplyr)
#library(biomaRt)
library(plotly)
#library(readxl)
library(ggrepel)
library(goseq)
library(clusterProfiler)
library(DOSE)
library(circlize)
library(readxl)
#library(ComplexHeatmap)
#library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(vsn)
library(annotables)
library(grid)
require(scales)
library(goseq)
library(org.Mm.eg.db)
library(pathview)

# Read counts from raw data and create DDS matrix for further DE analysis. 
# Give genotyep and create a condition column with number of replicates
# point the condition columns to the columns in the raw data 

counts <- read.table('/Users/hari/RNAseq/AnalysisR/Repcounts.txt', header = TRUE, row.names = 1)
counts <- counts[6:ncol(counts)]
genotype <- c("wt","wt","wt","wt","wt",'HSP1', 'HSP1','HSP1','HSP1','HSP1',"wt","wt","wt","wt","wt")
condition <- c(rep("Untreated",5),rep("HPS_untreated",5),rep("REP",5))
colData <- data.frame(genotype, condition)
rownames(colData) <- colnames(counts)
match(rownames(colData), colnames(counts))
reordered_counts <- counts[, match(rownames(colData), colnames(counts))]
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = colData,
                              design = ~condition)
dds <- dds[rowSums(counts(dds)) > 10, ]
dds <- estimateSizeFactors(dds)
dds_nor <- counts(dds, normalized = TRUE)
dds_nor <-as.data.frame(dds_nor) %>% rownames_to_column('ensgene')


# unsupervised clustering analysis using varianceStabilizingTransformation in R and 
# extract the matrix of transformed counts
#compute correlation between samples
vsd_dds <- vst(dds, blind = TRUE)
vsd_dds <- assay(vsd_dds)
vsd_dds_cor <- cor(vsd_dds)
head(vsd_dds_cor)
# plot the heatmap of VST
annotation <- colData %>% dplyr::select(condition)
pheatmap(vsd_dds_cor, annotation = annotation, fontsize = 20)

# Unsupervised clustering analysis using PCA analysis 
# the below code commented out didnot work for transformed values
#plotPCA(vsd_dds, intgroup = 'condition') so used rlogtransformation which is same
rld <- rlogTransformation(dds, blind = TRUE)
plotPCA(rld, intgroup = 'condition')+
  geom_point(size=10)+
  theme(text = element_text(size=20))
dds <- DESeq(dds)
plotDispEsts(dds)

# Comparing between repetitive and untreated group. 
# results(dds,contrast = c('condition_factor','level_to_compare','base_level, alpha = 0.05)
# here base level is untreated and level_to_compare is REP(repetitive group) 
#and significance value by alpha = standard alpha is 0.05 and if we want to increase biological genes 
# we can use lfcThreshold = 0.32 after alpha. 1.25 fold change is  = 0.32
res_REP <- results(dds,contrast = c('condition','REP','Untreated'), alpha = 0.05,  lfcThreshold = 0.32)
# MA plot shows the mean of the normalized counts vs log2fold change for all genes tested
plotMA(res_REP, ylim = c(-8,8))

# to improve the estimated fold changes we can use log2foldchange. 
res_REPLFC <- lfcShrink(dds, contrast = c('condition','REP','Untreated'), res = res_REP)
plotMA(res_REPLFC, ylim = c(-8,8))
mcols(res_REPLFC)
#uncomment to run
print(summary(res_REP))
sink('/Users/hari/RNAseq/AnalysisR/output_data/rep_wt_summary.txt')
print(summary(res_REPLFC))
sink()


# results for HPS vs Untreated
res_HPS <- results(dds, contrast = c('condition','HPS_untreated','Untreated'), alpha = 0.05)
plotMA(res_HPS, ylim = c(-8,8))
res_HPSLFC <- lfcShrink(dds, contrast = c('condition','HPS_untreated','Untreated'), res = res_HPS)
# cooksCuttoff use this to check outliers 
#res_HPS <- results(dds, cooksCutoff=FALSE)
plotMA(res_HPSLFC, ylim = c(-8,8))
mcols(res_HPSLFC)
# uncomment to run
print(summary(res_HPSLFC))
sink( '/Users/hari/RNAseq/AnalysisR/output_data/HPS_wt_summary.txt')
print(summary(res_HPSLFC))
sink()

# Annotating genes in repetive vs Untreated data with annotables and 
# merge the dataframes from DESeq2 analysis and normalized counts for visualizations and heatmap. 

res_REP_genes <- data.frame(res_REPLFC) %>% 
  rownames_to_column(var = "ensgene")
res_REP_genes <- left_join(res_REP_genes, dds_nor, by='ensgene')
all_genes <-  left_join(x = res_REP_genes, y = grcm38[, c('ensgene', 'symbol','description')], by ='ensgene')
Gm <- all_genes[!grepl("Gm", all_genes$symbol),]
Ig<- all_genes[!grepl("Ig", all_genes$symbol),]
all_genes<- merge(Gm,Ig)
all_genes <- all_genes[!duplicated(all_genes$ensgene),] 

#Found some pseudo genes and ig genes and filtering for pseudogenes and ig. it is not necessary but just to see what happens
#write.csv(as.data.frame(all_genes), file = '/Users/hari/RNAseq/AnalysisR/output_data/all_genes.csv', row.names = FALSE)
# Annotating genes in HPS vs Untreateds data 

res_HPS_genes <- data.frame(res_HPSLFC) %>% 
  rownames_to_column(var = "ensgene")
res_HPS_genes <- left_join(res_HPS_genes, dds_nor, by='ensgene')
all_HPS_genes <-  left_join(x = res_HPS_genes, y = grcm38[, c('ensgene', 'symbol','description')], by ='ensgene')

#BH Method of P value
sum(all_genes$padj < 0.05, na.rm = TRUE)
sum(all_HPS_genes$padj < 0.05, na.rm = TRUE)
REP_genes <- all_genes %>% filter(padj < 0.05) %>% arrange(padj)
HPS_genes <- all_HPS_genes %>% filter(padj < 0.05) %>% arrange(padj)
# uncomment to write to CSV file 
#write.csv(as.data.frame(REP_genes), file = '/Users/hari/RNAseq/AnalysisR/output_data/Rep_vs_wt_DE_genes.csv', row.names = FALSE)
#write.csv(as.data.frame(all_HPS_genes), file = '/Users/hari/RNAseq/AnalysisR/output_data/HPS_genes.csv', row.names = FALSE)
# HPS vs Untreated
# reading files gave 5655 genes
DE_rep_genes  <-  read.csv('/Users/hari/RNAseq/AnalysisR/output_data/Rep_vs_wt_DE_genes.csv', stringsAsFactors = FALSE)
DE_rep_genes <- DE_rep_genes[!duplicated(DE_rep_genes$symbol),] 
DE_rep_genes_up <- DE_rep_genes %>%  
                  filter(log2FoldChange >0) %>% 
                  arrange(padj)
DE_rep_genes_down <- DE_rep_genes %>%  
  filter(log2FoldChange < -1) %>% 
  arrange(padj)

# Writing up and down genes for analysis
#write.csv(as.data.frame(DE_rep_genes_up), file = '/Users/hari/RNAseq/AnalysisR/output_data//DE_rep_genes_up.csv', row.names = FALSE)
#write.csv(as.data.frame(DE_rep_genes_down), file = '/Users/hari/RNAseq/AnalysisR/output_data/DE_rep_genes_down.csv', row.names = FALSE)
#results_Rep_Untreated <- results_Rep_Untreated[!duplicated(results_Rep_Untreated$ensgene),] 
# drop NAs
DE_rep_genes <-DE_rep_genes  %>% drop_na()  %>% arrange(padj)
Rep_Untreated_DE_genes<-subset(DE_rep_genes[8:23])
Rep_Untreated_DE_genes <- Rep_Untreated_DE_genes %>% remove_rownames %>% column_to_rownames(var="symbol")
sum(is.na(DE_rep_genes))
str(Rep_Untreated_DE_genes)


# Choose heatmap color palette
heat_colors <- brewer.pal(n = 11, name = "RdBu")
annotation <- colData %>% dplyr::select(condition)
# Plot heatmap for 5435 genes
pheatmap(Rep_Untreated_DE_genes,
         color = heat_colors, 
         cluster_rows = T, 
         show_rownames = F,
         annotation = annotation,  fontsize = 15,
         scale = 'row')


top_50_Rep <-DE_rep_genes %>%  
  filter(log2FoldChange >1 | log2FoldChange < -1) %>% 
  arrange(padj) %>% 
  dplyr::slice(1:50)
Heatmap_top_50_Rep <- subset(top_50_Rep[8:23])
Heatmap_top_50_Rep<- Heatmap_top_50_Rep %>% remove_rownames %>% column_to_rownames(var="symbol")

heat_colors <- brewer.pal(n = 11, name = "RdBu")
pheatmap(Heatmap_top_50_Rep,
         color = heat_colors, 
         cluster_rows = T, 
         show_rownames = T,
         annotation = annotation, 
         scale = 'row', fontsize = 10,cex =  1.25,
         clustering_distance_rows = "correlation",#Pearson's
         # clustering_method = "average",
         fontsize_row=7)

# heat map for HPSvs Untreated

results_HPS_Untreated <-  read.csv('/Users/hari/RNAseq/AnalysisR/output_data/HPS_genes.csv', stringsAsFactors = FALSE)
results_HPS_Untreated <- results_HPS_Untreated[!duplicated(results_HPS_Untreated$ensgene),] 
HPS_Untreated_DE_genes <- results_HPS_Untreated  %>% 
                          drop_na() %>%  
                          filter(log2FoldChange >1 | log2FoldChange < -1) %>%  
                          arrange(padj)


# filter genes from normalized counts with significant genes 
HPS_Untreated_DE_genes<-subset(HPS_Untreated_DE_genes, select = c('symbol', "WT2", "WT3", "WT4", 
                                                                  "WT5", "WT6", "HPS1.9", "HPS1.10", 'HPS1.11', "HPS1.12", "HPS1.13"))
HPS_Untreated_DE_genes<- HPS_Untreated_DE_genes %>% remove_rownames %>% column_to_rownames(var="symbol")

str(HPS_Untreated_DE_genes)
# Choose heatmap color palette
heat_colors <- brewer.pal(n = 11, name = "RdBu")

# Plot heatmap
pheatmap(HPS_Untreated_DE_genes,
         color = heat_colors, 
         cluster_rows = T, 
         show_rownames = F,
         annotation = annotation, fontsize = 16,
         scale = 'row')
# slicing top 50 genes for heatmap. 
Heatmap_top_50_HPS <-results_HPS_Untreated  %>% 
                     drop_na() %>%  
                     filter(log2FoldChange >1 | log2FoldChange < -1) %>%  
                    arrange(padj) %>% 
                    dplyr::slice(1:50)
Heatmap_top_50_HPS <- subset(Heatmap_top_50_HPS, 
                             select = c('symbol', "HPS1.9", "HPS1.10", 
                                        'HPS1.11', "HPS1.12", "HPS1.13","WT2", "WT3", "WT4", 
                                        "WT5", "WT6")) %>% 
                        remove_rownames %>% column_to_rownames(var="symbol")

heat_colors <- brewer.pal(n = 11, name = "RdBu")
pheatmap(Heatmap_top_50_HPS,
         color = heat_colors, 
         cluster_rows = T, 
         show_rownames = T,
         annotation = annotation, fontsize = 10,cex =  1.25,
         scale = 'row',
         fontsize_row=7)


# plotting expression values 
colData2 <- colData %>% rownames_to_column(var = 'samplename')
top_25 <- DE_rep_genes[1:25, ] %>%  dplyr::select('ensgene','symbol', "WT2", "WT3", "WT4", "WT5", "WT6","REP16","REP17","REP18","REP19","REP20") %>% 
  gather(key = 'samplename', value = 'normalized_counts', WT2:REP20)

top_25 <- inner_join(top_25, colData2, by = 'samplename')

ggplot(top_25, aes(x=symbol, y= normalized_counts, color = condition))+
  geom_point(size = 3)+
  scale_colour_manual(values = c("Blue", "Red"))+
  scale_y_log10(labels=comma)+
  xlab('DE Genes')+
  ylab('Normalized Counts(log transformed)')+
  theme(text = element_text(size=20))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#### plotting top  25 upregulted genes
# plotting expression values 

colData2 <- colData %>% rownames_to_column(var = 'samplename')
top_25_upregulated <- DE_rep_genes %>% filter(log2FoldChange >1) 
top_25_upregulated <- top_25_upregulated[1:25, ]%>%  dplyr::select('ensgene','symbol', "WT2", "WT3", "WT4", "WT5", "WT6","REP16","REP17","REP18","REP19","REP20") %>% 
  gather(key = 'samplename', value = 'normalized_counts', WT2:REP20)


top_25_upregulated <- inner_join(top_25_upregulated, colData2, by = 'samplename')

ggplot(top_25_upregulated, aes(x=symbol, y= normalized_counts, color = condition))+
  geom_point(size = 5)+
  scale_colour_manual(values = c("Blue", "Red"))+
  scale_y_log10(labels=comma)+
  xlab('Top 25 DE Genes > 1 logfold change')+
  ylab('Normalized Counts(log transformed)')+
  theme(text = element_text(size=20))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#### plotting top  25 downregulated genes
# plotting expression values 
colData2 <- colData %>% rownames_to_column(var = 'samplename')
top_25_downregulated <- DE_rep_genes %>% filter(log2FoldChange < -1) 
top_25_downregulated <- top_25_downregulated[1:25, ]%>%  dplyr::select('ensgene','symbol', "WT2", "WT3", "WT4", "WT5", "WT6","REP16","REP17","REP18","REP19","REP20") %>% 
  gather(key = 'samplename', value = 'normalized_counts', WT2:REP20)


top_25_downregulated <- inner_join(top_25_downregulated, colData2, by = 'samplename')
ggplot(top_25_downregulated, aes(x=symbol, y= normalized_counts, color = condition))+
  geom_point(size = 5)+
  scale_colour_manual(values = c("Blue", "Red"))+
  scale_y_log10(labels=comma)+
  xlab('Top 25 DE Genes < -1 logfold change')+
  ylab('Normalized Counts(log transformed)')+
  theme(text = element_text(size=20))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#Volcono plot for Rep vs WT 

DE_rep_genes$up <- DE_rep_genes$log2FoldChange >  1  & DE_rep_genes$padj < 0.05
DE_rep_genes$down  <- DE_rep_genes$log2FoldChange < -1  & DE_rep_genes$padj < 0.05
DE_rep_genes$threshold <- as.factor(abs(DE_rep_genes$log2FoldChange) > 1 & DE_rep_genes$padj < 0.05)
sum(DE_rep_genes$up, na.rm = TRUE)
sum(DE_rep_genes$down , na.rm = TRUE)

# Slicing to annotate text in the volcano plot
top_25_upreg_genes <-  DE_rep_genes  %>%
  filter(log2FoldChange >1, padj< 0.05) %>% 
  arrange(padj) %>% dplyr::slice(1:25)
top_25_down_genes = DE_rep_genes %>%
  filter(log2FoldChange < -1, padj< 0.05) %>% 
  arrange(padj)%>% dplyr::slice(1:25)


#write.csv(as.data.frame(top_down_genes), file = '/Users/hari/RNAseq/AnalysisR/Output_files/downregulatedinrep.csv')
#write.csv(as.data.frame(top_down_genes), file = '/Users/hari/RNAseq/AnalysisR/Output_files/top50downregulatedinrep.csv')

volcano_plot <- ggplot(data=DE_rep_genes, aes(x=log2FoldChange, y=-log10(padj))) +
  geom_point(data=DE_rep_genes, size=1, colour="gray", aes(text = symbol)) +
  geom_point(data=DE_rep_genes[DE_rep_genes$up==TRUE, ], size=2, colour="Red",aes(text = symbol)) +
  geom_point(data=DE_rep_genes[DE_rep_genes$down  ==TRUE, ], size=2, colour="Blue",aes(text = symbol)) +
  #geom_point(aes(text = Genename))+
  xlab("log2 fold change") +
  ylab("-log10 p-value adjusted") +
  ggtitle("WT Rep Bleo vs WT Untreated")+
  geom_vline(xintercept = 1, colour="red", linetype="dashed")+
  geom_vline(xintercept = -1, colour="blue", linetype="dashed")+
  geom_text_repel(data=top_25_upreg_genes , aes(label=symbol), size = 4)+
  geom_text_repel(data=top_25_down_genes, aes(label=symbol), size = 4)+
  #scale_x_continuous() +
  xlim(-10,10)+
  scale_y_continuous() +
  annotate('text', x = 7, y = 150, label = '660 genes \n Upregulated \n log2foldchange >1 \n padj <0.05', color = 'Red', size = 6)+
  annotate('text', x = -7, y = 150, label = '1318 genes \n downregulated \n log2foldchange < -1 \n padj <0.05', color = 'Blue',size = 6)+
  theme(axis.title.y = element_text(face="bold", size=16),
        axis.title.x = element_text(face="bold", size=16, colour="black"),
        axis.text = element_text(size=20),
        legend.title =element_blank() ,
        legend.text = element_text(size = 12)) +
  theme(plot.title = element_text(hjust = 0.5))

# to convert plot to plotly uncomment the below code
#volcano_plot <- ggplotly(volcano_plot, hoverinfo = 'text')
volcano_plot

# HPS vs Untreated

results_HPS_Untreated$up <- results_HPS_Untreated$log2FoldChange >  1  & results_HPS_Untreated$padj < 0.05
results_HPS_Untreated$down   <- results_HPS_Untreated$log2FoldChange < -1  & results_HPS_Untreated$padj < 0.05
results_HPS_Untreated$threshold <- as.factor(abs(results_HPS_Untreated$log2FoldChange) > 1 & results_HPS_Untreated$padj < 0.05)

sum(results_HPS_Untreated$up, na.rm = TRUE)
sum(results_HPS_Untreated$down, na.rm = TRUE)

top_20_upreg_genes_HPS1 = results_HPS_Untreated %>%
  filter(log2FoldChange >1, padj< 0.05) %>% 
  arrange(padj)%>% dplyr::slice(1:20)


#write.csv(as.data.frame(top_upreg_genes_HPS1), file = '/Users/hari/RNAseq/AnalysisR/Output_files/upregulatedinHPS1.csv')

top_20_down_genes_HPS1 <-  results_HPS_Untreated%>%
  filter(log2FoldChange < -1, padj< 0.05) %>% 
  arrange(padj)%>% dplyr::slice(1:20)

#write.csv(as.data.frame(top_down_genes_HPS1), file = '/Users/hari/RNAseq/AnalysisR/Output_files/DownregulatedinHPS1.csv')
# y axis negative logarithm to the base 10 of the t-test p values
volcano_plot2 <- ggplot(data=results_HPS_Untreated, aes(x=log2FoldChange, y=-log10(padj))) +
  geom_point(data=results_HPS_Untreated, size=1, colour="gray",aes(text = symbol)) + 
  geom_point(data=results_HPS_Untreated[results_HPS_Untreated$up==TRUE, ], size=2, colour="Red",aes(text = symbol)) +
  geom_point(data=results_HPS_Untreated[results_HPS_Untreated$down  ==TRUE, ], size=2, colour="Blue",aes(text = symbol)) +
  geom_vline(xintercept = 1, colour="red", linetype="dashed")+
  geom_vline(xintercept = -1, colour="blue", linetype="dashed")+
  geom_text_repel(data=top_20_upreg_genes_HPS1, aes(label=symbol), size = 4)+
  geom_text_repel(data=top_20_down_genes_HPS1, aes(label=symbol), size = 4)+
  xlab("log2 fold change") +
  ylab("-log10 p-value adjusted") +
  ggtitle("HPS untreated vs WT_Untreated")+
  #scale_x_continuous() +
  scale_y_continuous() +
  xlim(-10,10)+
  annotate('text', x = 6, y = 50, label = '36 genes \n Upregulated \n log2foldchange >1 \n padj <0.05', color = 'Red',size = 6)+
  annotate('text', x = -6, y = 50, label = '92 genes \n downregulated \n log2foldchange < -1 \n padj <0.05', color = 'Blue',size = 6)+
  theme_bw() +
  theme(axis.title.y = element_text(face="bold", size=16),
        axis.title.x = element_text(face="bold", size=16, colour="black"),
        axis.text = element_text(size=12),
        legend.title =element_blank() ,
        legend.text = element_text(size = 12)) +
  theme(plot.title = element_text(hjust = 0.5))
#volcano_plot2 <- ggplotly(volcano_plot2)
volcano_plot2

# Generating rank file

all_genes_rnk <-  read.csv('/Users/hari/RNAseq/AnalysisR/output_data/all_genes.csv')

all_genes_rnk <- all_genes_rnk %>% arrange(padj)
all_genes_rnk$symbol <- toupper(all_genes_rnk$symbol)

#x<- x %>% filter(padj < 0.25)
all_genes_rnk$fcsign <- sign(all_genes_rnk$log2FoldChange)
all_genes_rnk$logP=-log10(all_genes_rnk$pvalue)
all_genes_rnk$metric= all_genes_rnk$logP/all_genes_rnk$fcsign
y<-all_genes_rnk[,c("symbol", "metric")]
sum(is.na(y$symbol))
sum(is.na(y$metric))
y <- y %>% drop_na()
y<- y[!duplicated(y$symbol),] 
head(y)


#write.table(y,file='/Users/hari/RNAseq/AnalysisR/output_data/hari_de.rnk',quote=F,sep="\t",row.names=F)

Reactome_up <- read.csv('/Users/hari/RNAseq/AnalysisR/Output_files/output datacamp DESeq2/gsea_report_for_na_pos_1551207478364.csv')
Reactome_down <- read.csv('/Users/hari/RNAseq/AnalysisR/Output_files/output datacamp DESeq2/gsea_report_for_na_neg_1551207478364.csv')
Reactome_up$NAME <- gsub("_", " ",  Reactome_up$NAME)
Reactome_up$NAME <-gsub("REACTOME", " ",  Reactome_up$NAME) 
Reactome_up <- Reactome_up %>% filter(FDR.q.val < 0.05) %>% dplyr::slice(1:40)
Reactome_down$NAME <- gsub("_", " ",  Reactome_down$NAME)
Reactome_down$NAME <-gsub("REACTOME", " ",  Reactome_down$NAME) 
Reactome_down <- Reactome_down %>% filter(FDR.q.val < 0.05) %>% dplyr::slice(1:40)


ggplot(Reactome_up, aes(reorder(NAME, NES), NES)) +
  geom_point(stat = 'identity',aes(size = NES, color = -NES)) + 
  coord_flip()+
  xlab('') +
  ylab('Normalized Enrichment Score(NES)') +
  theme(text = element_text(size=15))+
  ggtitle('Gene Set Enrichment Analysis: \n Upregulated Reactome Pathways') +
  theme(axis.title.y = element_text(size = 4))+
  theme(legend.position="none")+
  geom_text(aes(label=SIZE), size=4,colour = 'black',position = position_stack(vjust = 1.08))


ggplot(Reactome_down, aes(reorder(NAME, NES), NES)) +
  geom_point(stat = 'identity',aes(size = NES, color = -NES)) + 
  coord_flip()+
  xlab('') +
  scale_x_discrete(position = 'left')+
  ylab('Normalized Enrichment Score(NES)') +
  theme(text = element_text(size=15))+
  ggtitle('Gene Set Enrichment Analysis: \n Top 25 Downregulated Reactome Pathways') +
  theme(axis.title.y = element_text(size = 4))+
  theme(legend.position="none")+
  geom_text(aes(label=SIZE), size=4,colour = 'black',position = position_stack(vjust = -0.40))

KEGG_up <- read.csv('/Users/hari/RNAseq/AnalysisR/Output_files/output datacamp DESeq2/gsea_report_for_na_pos_KEGG_1551210816292.csv')
KEGG_down <- read.csv('/Users/hari/RNAseq/AnalysisR/Output_files/output datacamp DESeq2/gsea_report_for_na_neg_KEGG_1551210816292.csv')
KEGG_up $NAME <- gsub("_", " ",  KEGG_up $NAME)
KEGG_up $NAME <-gsub("KEGG", " ",  KEGG_up$NAME) 
KEGG_up  <-   KEGG_up %>% filter(FDR.q.val < 0.05) %>% dplyr::slice(1:35)
KEGG_down$NAME <- gsub("_", " ",  KEGG_down$NAME)
KEGG_down$NAME <-gsub("KEGG", " ",  KEGG_down$NAME) 
KEGG_down <- KEGG_down%>% filter(FDR.q.val < 0.05) %>% dplyr::slice(1:35)


ggplot(KEGG_up, aes(reorder(NAME, NES), NES)) +
  geom_point(stat = 'identity',aes(size = NES, color = -NES)) + 
    coord_flip()+
  xlab('') +
  ylab('Normalized Enrichment Score(NES)') +
  theme(text = element_text(size=15))+
  ggtitle('Gene Set Enrichment Analysis FDR < 0.05: \n Upregulated KEGG Pathways') +
  theme(axis.title.y = element_text(size = 4))+
  theme(legend.position="none")+
  geom_text(aes(label=SIZE), size=4,colour = 'black',position = position_stack(vjust = 1.08))
# geom_text(aes(label=RANK.AT.MAX), size=4,colour = 'yellow',position = position_stack(vjust = .8))


ggplot(KEGG_down, aes(reorder(NAME, NES), NES)) +
  geom_point(stat = 'identity',aes(size = NES, color = -NES)) + 
  coord_flip()+
  xlab('') +
  scale_x_discrete(position = 'left')+
  ylab('Normalized Enrichment Score(NES)') +
  theme(text = element_text(size=15))+
  ggtitle('Gene Set Enrichment Analysis FDR < 0.05: \n Downregulated KEGG Pathways') +
  theme(axis.title.y = element_text(size = 4))+
  theme(legend.position="none")+
  geom_text(aes(label=SIZE), size=4,colour = 'black',position = position_stack(vjust = -0.25))

# Hall mark Pathways

Hallmark_up <- read.csv('/Users/hari/RNAseq/AnalysisR/Output_files/output datacamp DESeq2/gsea_report_for_na_pos_hallmark_1551282126528.csv')
Hallmark_down <- read.csv('/Users/hari/RNAseq/AnalysisR/Output_files/output datacamp DESeq2/gsea_report_for_na_neg_hallmark_551282126528.csv')
Hallmark_up$NAME <- gsub("_", " ",  Hallmark_up$NAME)
Hallmark_up$NAME <-gsub("HALLMARK", " ",  Hallmark_up$NAME) 
Hallmark_up <-   Hallmark_up %>% filter(FDR.q.val < 0.05) %>% dplyr::slice(1:40)
Hallmark_down$NAME <- gsub("_", " ",  Hallmark_down$NAME)
Hallmark_down$NAME <-gsub("HALLMARK", " ",  Hallmark_down$NAME) 
Hallmark_down <- Hallmark_down %>% filter(FDR.q.val < 0.05) %>% dplyr::slice(1:25)


ggplot(Hallmark_up, aes(reorder(NAME, NES), NES)) +
  geom_point(stat = 'identity',aes(size = NES, color = -NES)) + 
  coord_flip()+
  xlab('') +
  ylab('Normalized Enrichment Score(NES)') +
  theme(text = element_text(size=15))+
  ggtitle('Gene Set Enrichment Analysis FDR < 0.05: \n HallMark Upregulated Pathways') +
  theme(axis.title.y = element_text(size = 4))+
  theme(legend.position="none")+
  geom_text(aes(label=SIZE), size=4,colour = 'black',position = position_stack(vjust = 1.08))
# geom_text(aes(label=RANK.AT.MAX), size=4,colour = 'yellow',position = position_stack(vjust = .8))


ggplot(Hallmark_down, aes(reorder(NAME, NES), NES)) +
  geom_point(stat = 'identity',aes(size = NES, color = -NES)) + 
  coord_flip()+
  xlab('') +
  scale_x_discrete(position = 'left')+
  ylab('Normalized Enrichment Score(NES)') +
  theme(text = element_text(size=15))+
  ggtitle('Gene Set Enrichment Analysis FDR < 0.05: \n  HallMark Downregulated Pathways') +
  theme(axis.title.y = element_text(size = 4))+
  theme(legend.position="none")+
  geom_text(aes(label=SIZE), size=4,colour = 'black',position = position_stack(vjust = -0.09))



## pathway exploration

pathway_genes <- read.csv('/Users/hari/RNAseq/AnalysisR/Output_files/output datacamp DESeq2/pathway_data/Short_listed_expression/total_genes.csv')
#pathway_genes <- pathway_genes %>% mutate(simple_ttest = as.numeric(simple_ttest))
str(pathway_genes)
# run from these for other pathways
#countmatrix_symbol <- read.csv('./data/foldchange_raw_values_rep_wt.csv')
test <- read.csv('/Users/hari/RNAseq/AnalysisR/Output_files/output datacamp DESeq2/pathway_data/HALLMARK_HYPOXIA.csv')
#count <- count %>% mutate(Genename = tolower(Genename)) %>% mutate(Genename = capitalize(Genename))
filtered_genes<- pathway_genes[pathway_genes[,9] %in% test[,2], ]
filtered_genes <- filtered_genes %>%  filter(padj< 0.05)  %>% arrange(padj)
#write.csv(as.data.frame(filtered_genes), file = '/Users/hari/RNAseq/AnalysisR/Output_files/output datacamp DESeq2/pathway_data/Short_listed_expression/hypoxia.csv')




### BAR PLOT 

genes_x <- read.csv('/Users/hari/RNAseq/AnalysisR/Output_files/output datacamp DESeq2/pathway_data/finalpathwaygenes.csv')
genes_x$pathway <- gsub("_", " ",  genes_x$pathway)
#genes_x<- genes_x[!duplicated(genes_x$Genename),]
c <- genes_x %>%  mutate(downregulated = log2FoldChange < 0) %>%  filter(downregulated == TRUE) %>%
  group_by(pathway) %>% dplyr::summarise(Downregulated =n()) 

d <- genes_x %>%  mutate(upregulated = log2FoldChange > 0) %>%  filter(upregulated == TRUE)%>%
  group_by(pathway) %>% dplyr::summarise(Upregulated =n())

cd <- full_join(x=c, y=d, by.x = "pathway", by.y = 'pathway')
cd[is.na(cd)] <- 0


#cd <- cd %>%  mutate(Total = DownregulatedGenes+UpregulatedGenes) %>% arrange(desc(Total))
cd.long <- gather(cd, Genes, Value, Downregulated:Upregulated)

cd.long <- cd.long %>% group_by(pathway,Genes) %>% 
  summarise(count = sum(Value,na.rm = T))

ggplot(cd.long, aes(reorder(pathway, count), count, fill = Genes)) +
  geom_bar(stat = "identity") + coord_flip()+
  geom_text(data=subset(cd.long, count>0), aes(y=count, label=count), size=4,position = position_stack(vjust = 0.5))+
  #geom_text(size = 3, position = position_stack(vjust = 0.5))+
  #geom_text(aes(label=Genes), vjust=-0.3, size=3.5)+
  scale_fill_manual(values = c("orange","green"))+
  xlab('Pathway') +
  theme(text = element_text(size=15))+
  theme(axis.title.y = element_text(size = 1))+
  ylab('Number of Genes') +
  ggtitle('        Gene Ontology pathway of DEG') 



# pathview

## showing all genes in pathway map 

pathway_genes2 <- read.csv('/Users/hari/RNAseq/AnalysisR/Output_files/output datacamp DESeq2/pathway_data/Short_listed_expression/total_genes.csv', stringsAsFactors = FALSE)
#pathway_genes2 <- pathway_genes2 %>% mutate(simple_ttest = as.numeric(simple_ttest))
test2 <- read.csv('/Users/hari/RNAseq/data/transfac_hif1a.csv')
filtered_genes2<- pathway_genes2[pathway_genes2[,9] %in% test2[,2], ] 
#write.csv(as.data.frame(filtered_genes2), file = '/Users/hari/RNAseq/AnalysisR/Output_files/output datacamp DESeq2/pathway_data/tgfb1_unfiltered.csv')

pathway_genes_filtered <- pathway_genes2 %>% filter(padj< 0.05)

gene.data <- as.numeric(filtered_genes2$log2FoldChange)
names(gene.data) <- filtered_genes2$symbol
pathID <- "04066"
pview <- pathview(gene.data=gene.data,
                  gene.idtype="symbol",
                  pathway.id=pathID,
                  species="mmu",
                  out.suffix="kegg_HIF1a",
                  kegg.native=T,
                  same.layer=T)


gene.data <- as.numeric(pathway_genes_filtered$log2FoldChange)
names(gene.data) <- pathway_genes_filtered$symbol
pathID <- "04350"
pview <- pathview(gene.data=gene.data,
                  gene.idtype="symbol",
                  pathway.id=pathID,
                  species="mmu",
                  out.suffix="kegg_tgftest_filtered",
                  kegg.native=T,
                  same.layer=T)

## pathview KEGG pathways
data(korg)
head(korg, 100)

##### expression plot for hypoxia genes 
colData2 <- colData %>% rownames_to_column(var = 'samplename')
Hypoxia_genes <- filtered_genes[1:33, ] %>%  dplyr::select('ensgene','symbol', "WT2", "WT3", "WT4", "WT5", "WT6","REP16","REP17","REP18","REP19","REP20") %>% 
  gather(key = 'samplename', value = 'normalized_counts', WT2:REP20)


# plotting expression values 
Hypoxia_genes <- inner_join(Hypoxia_genes, colData2, by = 'samplename')
ggplot(Hypoxia_genes, aes(x=symbol, y= normalized_counts, color = condition))+
  geom_point(size = 3)+
  scale_colour_manual(values = c("Blue", "Red"))+
  scale_y_log10(labels=comma)+
  xlab('Hypoxia DE Genes ')+
  ylab('Normalized Counts(log transformed)')+
  theme(text = element_text(size=20))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

