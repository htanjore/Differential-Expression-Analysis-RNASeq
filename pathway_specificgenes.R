setwd("~/RNAseq/AnalysisR/Differential-Expression-Analysis-RNASeq/")


library(tidyverse)
library(ComplexHeatmap)
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(tidyr)
library(ggpubr)
library(R.utils)
library(janitor)
library(jsonlite)


# reading total genes

total_genes <- read.csv('./data/total_genes_wt_rep.csv')

# reading json file from website and save it as csv and comment the write csv after saving file
# HIF1a CHIPSEQ Predicted
chea_hifa_json <- fromJSON("https://amp.pharm.mssm.edu/Harmonizome/api/1.0/gene_set/HIF1A/CHEA+Transcription+Factor+Targets")
chea_hifa_json <- chea_hifa_json$associations
#write.csv(as.data.frame(chea_hifa_json), file = './data/chea_hifa_json.csv')
# HIF1a Tranfec predicted
transfec_hifa_curated_json  <- fromJSON("https://amp.pharm.mssm.edu/Harmonizome/api/1.0/gene_set/HIF1A/TRANSFAC+Curated+Transcription+Factor+Targets")
transfec_hifa_curated_json <- transfec_hifa_curated_json$associations
#write.csv(as.data.frame(transfec_hifa_curated_json), file = './data/transfec_hifa_curated_json.csv')
transfec_hifa_predicted_json  <- fromJSON("https://amp.pharm.mssm.edu/Harmonizome/api/1.0/gene_set/HIF1A/TRANSFAC+Predicted+Transcription+Factor+Targets")
transfec_hifa_predicted_json <- transfec_hifa_predicted_json$associations
#write.csv(as.data.frame(transfec_hifa_predicted_json), file = './data/transfec_hifa_predicted_json.csv')
spib_Jaspar_predicted_json  <- fromJSON("https://amp.pharm.mssm.edu/Harmonizome/api/1.0/gene_set/SPIB/JASPAR+Predicted+Transcription+Factor+Targets")
spib_Jaspar_predicted_json <- spib_Jaspar_predicted_json$associations
#write.csv(as.data.frame(spib_Jaspar_predicted_json), file = './data/spib_Jaspar_predicted_json.csv')
DDIT3_transfec_predicted_json  <- fromJSON("https://amp.pharm.mssm.edu/Harmonizome/api/1.0/gene_set/DDIT3/TRANSFAC+Curated+Transcription+Factor+Targets")
DDIT3_transfec_predicted_json <- DDIT3_transfec_predicted_json$associations
#write.csv(as.data.frame(DDIT3_transfec_predicted_json), file = './data/DDIT3_transfec_predicted_json.csv')

pulmonary_Fibrosis_CTD  <- fromJSON("https://amp.pharm.mssm.edu/Harmonizome/api/1.0/gene_set/Pulmonary+Fibrosis/CTD+Gene-Disease+Associations")
pulmonary_Fibrosis_CTD <- pulmonary_Fibrosis_CTD$associations
pulmonary_Fibrosis_CTD 

chea_hifa <- read.csv('./data/chea_hifa_json.csv')
chea_hifa <- chea_hifa %>% dplyr::select(gene.symbol)
chea_hifa_list<- total_genes[total_genes[,8] %in% chea_hifa[,1 ], ]
chea_hifa_list <- chea_hifa_list %>% filter(padj<0.05)# %>% filter(log2FoldChange >1)%>% arrange(padj)
#write.csv(as.data.frame(chea_hifa_list), file = './data/chea_hifa_list.csv')



transfec_hifa_curated <- read.csv('./data/transfec_hifa_curated_json.csv')
transfec_hifa_curated <- transfec_hifa_curated %>% dplyr::select(gene.symbol)
transfec_hifa_curated_list<- total_genes[total_genes[,8] %in% transfec_hifa_curated[,1 ], ]
transfec_hifa_curated_list <- transfec_hifa_curated_list %>% filter(padj<0.05)# %>% filter(log2FoldChange >1) %>% arrange(padj)
#write.csv(as.data.frame(transfec_hifa_curated_list), file = './data/transfec_hifa_curated_list.csv')



transfec_hifa_predicted <- read.csv('./data/transfec_hifa_predicted_json.csv')
transfec_hifa_predicted <- transfec_hifa_predicted %>% dplyr::select(gene.symbol)
transfec_hifa_predicted_list<- total_genes[total_genes[,8] %in% transfec_hifa_predicted[,1 ], ]
transfec_hifa_predicted_list<- transfec_hifa_predicted_list %>% filter(padj<0.05) #%>% filter(log2FoldChange >1) %>% arrange(padj)
#write.csv(as.data.frame(transfec_hifa_predicted_list), file = './data/transfec_hifa_predicted_list.csv')


combined <- rbind(chea_hifa_list, transfec_hifa_curated_list, transfec_hifa_predicted_list)
combined<- combined[!duplicated(combined$ensgene),] %>% drop_na() %>% 
          arrange(padj) #%>% 
          #filter(log2FoldChange >1 | log2FoldChange < -1)
combined_hif1<-subset(combined , select = c('symbol', "WT2", "WT3", "WT4",
                                                                  "WT5", "WT6", "REP16", "REP17", 'REP18', "REP19", "REP20"))

combined_hif1<- combined_hif1 %>% remove_rownames %>% column_to_rownames(var="symbol")

colnames(combined_hif1)

genotype <- c("wt","wt","wt","wt","wt","wt","wt","wt","wt","wt")
Bleomycin<- c(rep("Untreated",5),rep("Repetitive",5))
colData <- data.frame(genotype, Bleomycin)
rownames(colData) <- colnames(combined_hif1)
annotation <- colData %>% dplyr::select(Bleomycin)
#sum(is.na(combined_hif1))
#str(combined_hif1)


heat_colors <- brewer.pal(n = 11, name = "RdBu")
pheatmap(combined_hif1,
         color = heat_colors, 
         cluster_rows = T, 
         show_rownames = F,
         cellwidth = 15,
         annotation = annotation, 
         scale = 'row', fontsize = 10,cex =  1.25,
         clustering_distance_rows = "correlation",#Pearson's
         # clustering_method = "average",
         fontsize_row=6)

top_75_Rep <-combined %>% drop_na() %>% 
  arrange(padj) %>%
  filter(log2FoldChange >1 | log2FoldChange < -1) %>%
  dplyr::slice(1:75) 

top_75_Rep<-subset(top_75_Rep , select = c('symbol',"WT2", "WT3", "WT4",
                                           "WT5", "WT6","REP16", "REP17", 'REP18', "REP19", "REP20"))

top_75_Rep<- top_75_Rep %>% remove_rownames %>% column_to_rownames(var="symbol")

# Choose heatmap color palette
heat_colors <- rev(brewer.pal(n = 11, name = "RdBu"))
annotation <- colData %>% dplyr::select(Bleomycin)
# Plot heatmap for 5435 genes
pheatmap(top_75_Rep,
         color = inferno(256, direction = -1), 
         #color =heat_colors, 
         cluster_rows = T, 
         show_rownames = T,
         cellwidth=15,
         annotation = annotation,  
         scale = 'row', fontsize = 8,cex =  1.25,
         #clustering_distance_rows = "correlation",#Pearson's
         clustering_method = "average",
         fontsize_row=5)





