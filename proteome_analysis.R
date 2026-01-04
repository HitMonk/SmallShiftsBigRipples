###Proteome data
protDat <- readxl::read_xlsx("/home/pyre/Project_data/mittag_group/temperature_Transcriptome/Finalized_Complete_results/20250704_Flagella_Supernatant/20250702_Flagella/FragPipeAnalysis/Flagella_Full_dataset.xlsx")
annotD <- read.csv("/home/pyre/Project_data/mittag_group/temperature_Transcriptome/annotD.txt",
                     header = T)
colnames(annotD) <- c("GeneID","GeneDescription")

cilDB_clean_annot <- read.table("~/Project_data/mittag_group/temperature_Transcriptome/cilDB_annotations_compartments.txt",
                                header = T)
cilDB_clean_annot <- cilDB_clean_annot %>%
  mutate_at(vars(Axo, `M.M`,
                 KC1,Terg ), ~replace(., is.na(.), 0))

# Using tidyverse/dplyr
# Using stringr
protDat$GeneID <- stringr::str_extract(protDat$Protein.ID, "^[^.]+\\.+[^_]+_4532")
protDat$X18_vs_X28_log2.fold.change <- protDat$`T18_vs_T28_log2-fold difference`*-1
protDat_annot <- merge(protDat, annotD, by.x="GeneID", by.y="GeneID", all.x=TRUE)[,c(1,2,57,53,58)]
write.table(protDat_annot, "/home/pyre/Project_data/mittag_group/temperature_Transcriptome/2025FinalizedPaperJulyAugust/Supplemental_tables+updatedFigures/ProtDatComplete3271.txt",
            sep="\t", quote = F, row.names = F)

#protDatF <- protDat[,-c(1,2,10)]
protDatF_id <- protDat_annot[protDat_annot$GeneID %in% cilDB_clean_annot$GID,]
protDatF_mis <- protDat_annot[!(protDat_annot$GeneID %in% cilDB_clean_annot$GID),]
write.table(protDatF_id,"/home/pyre/Project_data/mittag_group/temperature_Transcriptome/Finalized_Complete_results/CilDBProt_newChiaChi.txt",
            sep="\t", quote = F,
            row.names = F)
protDatF_mis_merged_data <- merge(protDatF_mis, 
                                  filterT18v28_desc, 
                                  by.x = "GeneID",
                                  by.y = "name",
                                  all = FALSE
)

merged_data <- merge(protDatF_id, 
                     filterT18v28_desc, 
                     by.x = "GeneID", 
                     by.y = "name", 
                     all = FALSE
)

IFTFlaProt <- merge(protDat_annot, impIdsComp, by="GeneID")

protDat_annotSig <- protDat_annot[
  (protDat_annot$X18_vs_X28_log2.fold.change > 0.5 | protDat_annot$X18_vs_X28_log2.fold.change < -0.5) &
    (protDat_annot$T18_vs_T28_p.adj < 0.05),]

protDatF_idSig <- protDat_annotSig[protDat_annotSig$GeneID %in% cilDB_clean_annot$GID,]
protDatF_misSig <- protDat_annotSig[!(protDat_annotSig$GeneID %in% cilDB_clean_annot$GID),]
protDatF_mis_merged_dataSig <- merge(protDatF_misSig, 
                                  resT18v28_desc, 
                                  by.x = "GeneID",
                                  by.y = "name",
                                  all.x = T
)

merged_dataSig <- merge(protDatF_idSig, 
                        resT18v28_desc, 
                     by.x = "GeneID", 
                     by.y = "name", 
                     all.x = T
)


write.table(protDatF_mis_merged_dataSig, file = "Sig_Non-cilia proteins.txt",
            sep="\t", quote = F, row.names = F)
write.table(merged_dataSig, file = "Sig_Cilia proteins.txt",
            sep="\t", quote = F, row.names = F)


# Define thresholds
FC <- 0.5  # Fold-change threshold
pvalAdj <- 0.05  # Adjusted p-value threshold

# Calculate the number of upregulated and downregulated genes
deg <- protDatF_id %>% 
  filter(T18_vs_T28_p.adj < pvalAdj) %>%
  summarize(
    up_deg = sum(X18_vs_X28_log2.fold.change > FC),
    down_deg = sum(X18_vs_X28_log2.fold.change < -FC)
  )

# Define custom colors for the points
keyvals <- ifelse(
  protDatF_id$X18_vs_X28_log2.fold.change < -FC & protDatF_id$T18_vs_T28_p.adj < 0.05, 'royalblue',
  ifelse(protDatF_id$X18_vs_X28_log2.fold.change > FC & protDatF_id$T18_vs_T28_p.adj < 0.05, 'darkorange2',
         'black')
)
keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'darkorange2'] <- 'Up regulated'
names(keyvals)[keyvals == 'black'] <- 'Not significant'
names(keyvals)[keyvals == 'royalblue'] <- 'Down regulated'


EnhancedVolcano(
  protDatF_id,
  lab = "",  # Labels for points
  x = 'X18_vs_X28_log2.fold.change',  # Log2 fold-change column
  y = 'T18_vs_T28_p.adj',  # Adjusted p-value column
  xlab = "Log2 Fold Change",
  ylab = "-Log10 Adjusted P-value",
  title = "Volcano Plot",
  pCutoff = pvalAdj,  # Adjusted p-value threshold
  FCcutoff = FC,  # Fold-change threshold
  colCustom = keyvals,  # Custom colors
  pointSize = 1.5,  # Size of points
  labSize = 9.0,  # Size of labels
  subtitle = paste0("Up: ", deg$up_deg, ", Down: ", deg$down_deg),  # Subtitle with counts
  legendPosition = 'right',  # Position of the legend
  legendLabSize = 10,  # Size of legend labels
  legendIconSize = 5.0,  # Size of legend icons
  drawConnectors = TRUE,  # Draw connectors for labels
  widthConnectors = 0.5,  # Width of connectors
  colConnectors = 'grey50',  # Color of connectors
  gridlines.major = TRUE,  # Show major gridlines
  gridlines.minor = FALSE,  # Hide minor gridlines
  border = 'full',  # Full border
  borderWidth = 1.0,  # Width of the border
  borderColour = 'black',  # Color of the border
  ylim=c(0,4.5)
)


###Look if any photoreceptors are identified
protDat_annot[protDat_annot$GeneID %in% photo$Cre.ID,]

#####Other plots
protplotIds <- read.csv("~/Project_data/mittag_group/temperature_Transcriptome/Finalized_Complete_results/Final_tables/proteome_tubulin_vals.csv")


##Pull out FAPs
faps_id <- protDat_annot[protDat_annot$GeneDescription %like% "Flagellar"& 
                protDat_annot$X18_vs_X28_p.adj < 0.05,]
write.table(faps_id, "faps_id.txt",
            sep="\t", quote = F, row.names = F)


#####Other plots
protplotIds <- read.csv("~/Project_data/mittag_group/temperature_Transcriptome/Finalized_Complete_results/Final_tables/protein_sub.csv")

protImp <- merge(protplotIds, protDat_annot, by.x="GeneID", by.y="GeneID")
colnames(protImp)[5] <- "log2FoldChange"
colnames(protImp)[6] <- "padj"

protImp_new <- protImp %>%
  mutate(
    label = case_when(   # Customize significance labels based on p-value
      padj < 0.0001 ~ "****",  
      padj < 0.001  ~ "***",    
      padj < 0.01   ~ "**",     
      padj < 0.05   ~ "*", 
      padj >= 0.05 ~ "",
    )
  )

protImp_new$Direction <- ifelse(as.numeric(protImp_new$log2FoldChange) > 0, "Positive", "Negative")
protImp_new$log2FoldChange <- as.numeric(abs(protImp_new$log2FoldChange))
# Arrange by Order column and create factor levels accordingly
protImp_new <- protImp_new %>%
  arrange(Order) %>%  # Sort by your Order column (1-24)
  mutate(Def = factor(Def, levels = unique(Def)))  # Preserve the Order-based sequence

ggplot(protImp_new, aes(x = log2FoldChange, y = reorder(Def, -Order), fill=Direction)) +
  geom_bar(stat = "identity", width=0.75,    
           position=position_dodge(0.75), colour="black") +
  scale_fill_manual(values = c("Negative" = "royalblue", "Positive" = "darkorange2")) +
  geom_text(aes(label = label), vjust = 0.8, hjust=-0.1,size = 10, position = position_dodge(width = 0.9)) +
  #facet_grid(Function ~ ., scales = "free_y", space = "free_y", switch = "y")+
  theme_minimal()+ xlim(0,6.5)+
  labs(y= "Gene names", x="Log2 Foldchange")+
  theme(axis.title.x = element_text(face = "bold", size=20), 
        axis.title.y = element_text(face = "bold", size=20),
        axis.text.x = element_text(face="bold", color="#000000", size=20),
        axis.text.y = element_text(face="bold", color="#000000", size=20),
        legend.title=element_text(face="bold", size=5),
        panel.grid.major = element_line(color = "grey", size=1), 
        strip.text.x = element_text(size = 18),
        strip.text = element_text(size=10),
        strip.placement = "outside")
