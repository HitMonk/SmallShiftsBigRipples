library(tidyverse)
library(tximport)
library('DESeq2')
library(EnhancedVolcano)
library("pheatmap")
library(ggpubr)
library("ComplexHeatmap")
library(tximportData)
library("topGO")
library(biomaRt)
library(VennDiagram)
library(dplyr)
library(Rhdf5lib)
library(rhdf5)
library(data.table)
library(AnnotationForge)
library(ggpattern)
library(readxl)
library(clusterProfiler)

setwd("/home/pyre/Project_data/mittag_group/temperature_Transcriptome/")

annotD <- read.csv("Gene2GeneDesc.txt",
            sep="\t", head = T)
keggs <- read.csv("Gene2annotation.txt",
            sep="\t", head = T)
t2g <- read.csv("Trans2Gene.csv")
geneID2GO <- readMappings("~/Project_data/mittag_group/temperature_Transcriptome/Gofile.txt")
GOTab <- read.csv("~/Project_data/mittag_group/temperature_Transcriptome/Gofile.txt", sep="\t")
GOTab_long <- GOTab %>%
  separate_rows(GO_term, sep = ",\\s*")

head(geneID2GO)
#make GO ontology universe file
#GOfile <- GO %>% group_by(gene_name) %>% reframe(GO_term =str_c(sort(go_id), collapse = ", "))
#write.table(x = GOfile, "Gofile.txt", sep="\t", quote = F, row.names = F)

#setup paths
dir <- file.path("~/Project_data/mittag_group/temperature_Transcriptome/Salmon_Counts_MA_tempTranscriptome_allTrans/")
list.files(dir)
s2c <- read.table(file.path(dir,"cre_new_mono.csv"), header=TRUE,
                  stringsAsFactors=FALSE, sep=",")
files <- file.path(dir,s2c$Filename,"quant.sf")                           
names(files) <- paste0(s2c$Abv)

#import using TXimport
txi.kallisto <- tximport(files, type = "salmon",tx2gene = t2g, 
                         txIn = TRUE, txOut = FALSE, countsFromAbundance = "no")
lapply(txi.kallisto, head)
counts_comp <- as.data.frame(txi.kallisto$counts)
counts_comp$name <- rownames(txi.kallisto$counts)
abund_comp <- as.data.frame(txi.kallisto$abundance)
abund_comp$name <-rownames(txi.kallisto$abundance)

ddsComp <- DESeqDataSetFromTximport(txi.kallisto, s2c, design=~Condition)
keepComp <- rowSums(counts(ddsComp) >= 10) >= 6
dds <- ddsComp[keepComp,]
dds <- DESeq(dds)
res <- results(dds)

# summary(res)
rld <- rlog(dds, blind = FALSE)
plotPCA(rld, intgroup="Condition")+geom_point(size=8, pch=0.2)+
  #geom_text_repel(aes(label = name))+ ggtitle("PCA of Cre6.1 Monoculture")+ 
  theme_bw()+
  scale_color_manual(values = c("T18" = "#258CCC", "T23" = "black", 
                                "T28" = "#428033", "T33" = "#FD3DB5"))+
  theme(plot.title = element_text(color="Black", size=22, face="italic",hjust = 0.5),
        legend.text = element_text(size=18),
        legend.title = element_text(size=22),
        axis.title =element_text(size=22),
        axis.text=element_text(size=18)) 


mapped_counts <- as.data.frame(round(colSums(as.data.frame(txi.kallisto$counts))))
colnames(mapped_counts) <- "mapped"
mapped_counts$sample <- row.names(mapped_counts)


#RNAseq filter plots function
rna_pipe <- function(desObj=dds,txi_mat=txi.kallisto$abundance,
                     cond1, cond2, samplefile=s2c, pvalAdj=0.01, FC=1){
  resCont <- results(desObj, contrast=c("Condition", cond1, cond2))
  resDF <- as.data.frame(resCont[order(resCont$padj),])
  resDF$name <- row.names(resDF)
  
  ## filtered DF
  filterDF <- resDF %>%
    filter(padj <= pvalAdj & 
             (log2FoldChange >= FC | log2FoldChange <= -(FC)))
  
  
  ## count and store degs in variable deg
  deg <- resDF %>% 
    filter(padj <= pvalAdj) %>%
    summarize(up_deg=sum(log2FoldChange >FC),
              down_deg=sum(log2FoldChange < -(FC)))
  keyvals <- ifelse(
    resDF$log2FoldChange < -FC & resDF$padj < 0.01,'royalblue',
    ifelse(resDF$log2FoldChange > FC & resDF$padj < 0.01,'darkorange2',
           'black'))
  keyvals[is.na(keyvals)] <- 'black'
  names(keyvals)[keyvals == 'darkorange2'] <- 'Up regulated'
  names(keyvals)[keyvals == 'black'] <- 'Not significant'
  names(keyvals)[keyvals == 'royalblue'] <- 'Down regulated'
  
  title.text <- paste(cond1,'vs', cond2,
                      'filtered at Padj:',pvalAdj,
                      'FC:',FC,
                      '\nDEG:',sum(deg[1],deg[2]),
                      '; UP:',deg[1],
                      'Down:',deg[2],
                      sep=" "
  )
  
  vol.p <- EnhancedVolcano(resDF,
                           lab = "",
                           colCustom = keyvals,
                           x = 'log2FoldChange',
                           y = 'pvalue',
                           title = title.text,
                           pCutoff = pvalAdj,
                           pCutoffCol="padj",
                           FCcutoff = FC,
                           pointSize = 1.5,
                           labSize = 9.0,
                           subtitle = "",
                           legendPosition = 'right',
                           legendLabSize = 0,
                           legendIconSize = 5.0,
                           drawConnectors = TRUE,
                           widthConnectors = 0.5,
                           colConnectors = 'grey50',
                           gridlines.major = TRUE,
                           gridlines.minor = FALSE,
                           border = 'full',
                           borderWidth = 1.0,
                           borderColour = 'black')
  return(list(filterDF, resDF, vol.p))  
}  

T18v23 <- rna_pipe(cond1="T23", cond2="T18",
                   desObj=dds,txi_mat=txi.kallisto$abundance)
filterT18v23 <- as.data.frame(T18v23[1])
resT18v23 <- as.data.frame(T18v23[2])
T18v23[3]
#ggsave("Ppr_Tri.tiff",units="in", width = 10 , height = 10 ,dpi = 300,compression="lzw")
filterT18v23_desc <- merge(merge(filterT18v23, annotD), counts_comp, by="name", all=F) %>%
  dplyr::select(-lfcSE, -baseMean, -stat, -pvalue)
resT18v23_desc <- merge(merge(resT18v23, annotD), counts_comp, by="name", all=F) %>%
  dplyr::select(-lfcSE, -baseMean, -stat, -pvalue)

T18v28 <- rna_pipe(cond1="T28", cond2="T18",
                   desObj=dds,txi_mat=txi.kallisto$abundance)
filterT18v28 <- as.data.frame(T18v28[1])
resT18v28 <- as.data.frame(T18v28[2])
T18v28[3]
#ggsave("Ppr_Tri.tiff",units="in", width = 10 , height = 10 ,dpi = 300,compression="lzw")
filterT18v28_desc <- merge(merge(filterT18v28, annotD), counts_comp, by="name", all=F) %>%
  dplyr::select(-lfcSE, -baseMean, -stat, -pvalue)
resT18v28_desc <- merge(merge(resT18v28, annotD), counts_comp, by="name", all=F) %>%
  dplyr::select(-lfcSE, -baseMean, -stat, -pvalue)

T18v33 <- rna_pipe(cond1="T33", cond2="T18",
                   desObj=dds,txi_mat=txi.kallisto$abundance)
filterT18v33 <- as.data.frame(T18v33[1])
resT18v33 <- as.data.frame(T18v33[2])
T18v33[3]
#ggsave("Ppr_Tri.tiff",units="in", width = 10 , height = 10 ,dpi = 300,compression="lzw")
filterT18v33_desc <- merge(merge(filterT18v33, annotD), counts_comp, by="name", all=F) %>%
  dplyr::select(-lfcSE, -baseMean, -stat, -pvalue)
resT18v33_desc <- merge(merge(resT18v33, annotD), counts_comp, by="name", all=F) %>%
  dplyr::select(-lfcSE, -baseMean, -stat, -pvalue)


###Heatmap
top_genes <- filterT18v28_desc %>% arrange(desc(log2FoldChange)) %>% head(500)
bottom_genes <- filterT18v28_desc %>% arrange(log2FoldChange) %>% head(500)
selected_genes <- bind_rows(top_genes, bottom_genes) %>% distinct(name, .keep_all=TRUE)

rld_sign <- assay(rld)[selected_genes$name,]

# Create column annotations
annotation_col <- data.frame(
  Group = sapply(colnames(rld_sign), function(col) {
    if(grepl("^Cre_18", col)) return("T18")
    if(grepl("^Cre_23", col)) return("T23")
    if(grepl("^Cre_28", col)) return("T28")
    if(grepl("^Cre_33", col)) return("T33")
    return(NA)
  })
)
rownames(annotation_col) <- colnames(rld_sign)
annotation_colors <- list(
  Group = c("T18" = "blue", "T23" = "black", "T28" = "green", "T33" = "red")
)


library(factoextra)

# Function to determine optimal k using the elbow method
find_optimal_k <- function(data, max_k = 20) {
  # Scale the data
  data_scaled <- t(scale(t(data)))
  
  # Calculate within-cluster sum of squares for different k values
  wss <- numeric(max_k)
  for (k in 1:max_k) {
    km <- kmeans(data_scaled, centers = k, nstart = 25)
    wss[k] <- km$tot.withinss
  }
  
  # Plot elbow curve
  plot(1:max_k, wss, type="b", pch=19,
       xlab="Number of clusters K",
       ylab="Total within-clusters sum of squares")
  
  # Find elbow point (simplified approach)
  diffs <- diff(wss)
  elbow_point <- which.min(abs(diffs - mean(diffs))) + 1
  abline(v = elbow_point, col = "red", lty = 2)
  
  return(elbow_point)
}

set.seed(1234)  # For reproducibility
optimal_k <- find_optimal_k(rld_sign)
cat("Optimal number of clusters:", optimal_k, "\n")

set.seed(123)
k <-   ComplexHeatmap::pheatmap(rld_sign, scale="row",
                                row_km=8,
                                show_rownames = F,cluster_cols=F,cluster_rows=T,
                                clustering_distance_rows = "euclidean", 
                                clustering_distance_cols = "euclidean", 
                                clustering_method = "complete",border_color = FALSE)
k_drawn <- ComplexHeatmap::draw(k)

# 3. Extract the row order from the *drawn object*
#    This will no longer produce a warning.
cluster_indices_list <- ComplexHeatmap::row_order(k_drawn)

# 4. Get the gene names from your original data
all_gene_names <- rownames(rld_sign)

# 5. Map the indices to the gene names
cluster_gene_list <- lapply(cluster_indices_list, function(indices) {
  all_gene_names[indices]
})


####GeneOntology Enrichment
#begin GO enrichment analysis
#Enrichment function
goFunc <- function(onto="BP", opt="up",fTable){
  if(opt=="up"){
    deg <- fTable %>% filter(log2FoldChange >0)
  }
  else if(opt=="down"){
    deg <- fTable %>% filter(log2FoldChange <0)
  }
  else if(opt=="all"){
    deg <- fTable
  }
  geneList <- factor(as.integer(annotD$name %in% deg$name))
  names(geneList) <- annotD$name
  length(geneList)
  
  GO <- new("topGOdata", ontology = onto, allGenes = geneList,
            annot = annFUN.gene2GO, gene2GO = geneID2GO)
  resultweight01 <- runTest(GO, algorithm = "weight01", statistic = "fisher")
  #create datatable:
  allRes <- GenTable(GO, weight01FS = resultweight01,
                     topNodes = 20, numChar=120)
  p.adj=round(p.adjust(allRes$weight01FS,method="BH"),digits = 4)
  all_res_final <- cbind(allRes,p.adj, Cat=onto) %>% subset(p.adj < 0.05)
  return(all_res_final)
}

#enrichment 
B_up <- goFunc(fTable = filterT18v28,onto = "BP", opt="up")
M_up <- goFunc(fTable = filterT18v28,onto = "MF", opt="up")
C_up <- goFunc(fTable = filterT18v28,onto = "CC", opt="up")
GO_up<- rbind(B_up, M_up, C_up)

ggplot(GO_up, aes(y = reorder(Term, as.numeric(as.factor(Cat))), 
                  x = Significant, size = Annotated)) + 
  geom_point(aes(color = p.adj), alpha = 10.0) + 
  geom_tile(aes(width = Inf, fill = Cat), alpha = 0.1) + 
  scale_fill_manual(values = c("green", "red", "blue"))+
  xlab("Number of genes in genelist")+ylab("GO term")+ 
  ggtitle("Upegulated enriched terms in \n T18 vs T23")+ 
  theme_bw()+ theme(plot.title = element_text(color="Black", size=22, face="italic",hjust = 0.5),,
                    legend.text = element_text(size=18),
                    legend.title = element_text(size=22),
                    axis.title =element_text(size=22),
                    axis.text=element_text(size=20, color="Black")) 
#ggsave("upreg_Ppr_Pfl.tiff",units="in", width = 24 , height = 12 ,dpi = 300,
#       compression="lzw")


B_down <- goFunc(fTable = filterT18v28,onto = "BP", opt="down")
M_down <- goFunc(fTable = filterT18v28,onto = "MF", opt="down")
C_down <- goFunc(fTable = filterT18v28,onto = "CC", opt="down")
GO_down<- rbind(B_down, M_down, C_down)

ggplot(GO_down, aes(y = reorder(Term, as.numeric(as.factor(Cat))), 
                    x = Significant, size = Annotated)) + 
  geom_point(aes(color = p.adj), alpha = 10.0) + 
  geom_tile(aes(width = Inf, fill = Cat), alpha = 0.1) + 
  ggtitle("Downregulated enriched terms in \n Cre_Pfl vs Cre_Ppr") + 
  scale_fill_manual(values = c("green", "red", "blue"))+
  xlab("Number of genes in genelist")+ylab("GO term")+
  theme_bw()+ theme(plot.title = element_text(color="Black", size=22, face="italic",hjust = 0.5),
                    legend.text = element_text(size=18),
                    legend.title = element_text(size=22),
                    axis.title =element_text(size=22),
                    axis.text=element_text(size=20, color="Black")) 



#####Output merge transcriptome results
# First: Rename log2FC and padj columns before joining for clarity
resT18v23_sub <- resT18v23_desc %>%
  dplyr::select(
    name, 
    Log2FC_18v23 = log2FoldChange, 
    padj_18v23 = padj,
    gene_description,
    Cre_18.1, Cre_18.2, Cre_18.3,
    Cre_23.1, Cre_23.2, Cre_23.3,
    Cre_28.1, Cre_28.2, Cre_28.3,
    Cre_33.1, Cre_33.2, Cre_33.3
  )

resT18v28_sub <- resT18v28_desc %>%
  dplyr::select(
    name, 
    Log2FC_18v28 = log2FoldChange, 
    padj_18v28 = padj
  )

resT18v33_sub <- resT18v33_desc %>%
  dplyr::select(
    name, 
    Log2FC_18v33 = log2FoldChange, 
    padj_18v33 = padj
  )

# Join all dataframes by 'name'
summary_df <- resT18v23_sub %>%
  left_join(resT18v28_sub, by = "name") %>%
  left_join(resT18v33_sub, by = "name")

# Keep only requested columns and optionally reorder
summary_df <- summary_df %>%
  dplyr::select(
    name,
    Log2FC_18v23, padj_18v23,
    Log2FC_18v28, padj_18v28,
    Log2FC_18v33, padj_18v33,
    gene_description,
    Cre_18.1, Cre_18.2, Cre_18.3,
    Cre_23.1, Cre_23.2, Cre_23.3,
    Cre_28.1, Cre_28.2, Cre_28.3,
    Cre_33.1, Cre_33.2, Cre_33.3
  )


KeggP <- read.csv("2025FinalizedPaperJulyAugust/KeggP/NewRun/KeggPath_SeparatedComplete.csv")
GeneSymb <- read.csv("References/algal_references/geneName_Symbol.csv")

KeggP_annotated <- KeggP %>%
  # 1. Join the tables matching Cre.ID to ens_gene
  left_join(GeneSymb, by = c("Cre.ID" = "ens_gene"))

KeggP_annotated <- KeggP %>%
  left_join(GeneSymb, by = c("Cre.ID" = "ens_gene")) %>%
  mutate(
    # 1. Convert empty strings to NA
    geneSymbol = na_if(geneSymbol, ""),
    # 2. Fill NAs with Cre.ID
    geneSymbol = coalesce(geneSymbol, Cre.ID)
  )

# View result
head(KeggP_annotated$geneSymbol)

#write.table(summary_df,"2025FinalizedPaperJulyAugust/Transcriptome_Table.txt",
#           sep="\t", quote = F, row.names = F)

##Get imp ids only
ImpIds <- read_xls("2025FinalizedPaperJulyAugust/Tab5.xls", sheet=3)
head(ImpIds)
colnames(ImpIds)[1] <- "name"
tab5 <- merge(ImpIds, summary_df, by="name", all.x=T)
#write.table(tab5, "2025FinalizedPaperJulyAugust/selected_ids.txt",
#          sep="\t", quote = F, row.names = F)


#####KEGG Path
Keggs_clean <-keggs %>%
  select(name, kegg_id) %>%
  distinct() %>%
  filter(kegg_id != "" & !is.na(kegg_id)) # Remove empty mappings


sig_gene_list <- filterT18v33$name
universe_list <- unique(Keggs_clean$name)

kegg_data <- download_KEGG(species = "ko", keggType = "KEGG", keyType = "kegg")
kegg_pathway_to_ko <- kegg_data$KEGGPATHID2EXTID

term2gene_df <- kegg_pathway_to_ko %>%
  rename(kegg_id = to, pathway_id = from) %>%
  # Remove the "ko" prefix from pathway IDs for cleaner merging if needed
  mutate(pathway_id = str_replace(pathway_id, "ko", "map")) %>% 
  inner_join(Keggs_clean, by = "kegg_id", relationship = "many-to-many") %>%
  select(pathway_id, name) # Format must be: Column 1 = Term (Pathway), Column 2 = Gene

term2name_df <- kegg_data$KEGGPATHID2NAME %>%
  mutate(from = str_replace(from, "ko", "map"))

kegg_enrichment <- enricher(
  gene = sig_gene_list,         # The significant genes
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",         # Benjamini-Hochberg correction
  universe = universe_list,     # Background genes
  minGSSize = 10,               # Exclude small pathways (<10 genes annotated)
  maxGSSize = 500,              # Exclude huge, vague pathways
  TERM2GENE = term2gene_df,     # The map we created: Pathway -> Gene
  TERM2NAME = term2name_df      # The map: Pathway ID -> Pathway Name
)

Keggenrichr <- (kegg_enrichment@result)

# Dotplot (The standard visualization)
dotplot(kegg_enrichment, showCategory = 25, title = "KEGG Pathway Enrichment")
results_df <- as.data.frame(kegg_enrichment)

