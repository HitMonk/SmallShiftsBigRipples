load("/home/pyre/Project_data/mittag_group/temperature_Transcriptome/Finalized_Complete_results/datasets_Phytozome.Rdata")

library(ggpattern)
###Load secretome data
secDat <- readxl::read_xlsx("/home/pyre/Project_data/mittag_group/temperature_Transcriptome/Finalized_Complete_results/20250704_Flagella_Supernatant/20250702_Supernatant/FargPipeAnalysis/Supernatant_Full_dataset.xlsx")
#secID <- read.csv("/home/pyre/Project_data/mittag_group/temperature_Transcriptome/Finalized_Complete_results/IdeentifiedSecProt.csv", head=F)
#colnames(secID) <- "GeneID"

deepLoc <- readxl::read_xlsx("/home/pyre/Project_data/mittag_group/temperature_Transcriptome/Finalized_Complete_results/deepLoc.xlsx")
unique(deepLoc$name)

secDat$GeneID <- stringr::str_extract(secDat$Protein.ID, "^[^.]+\\.+[^_]+_4532")
secDat$X18_vs_X28_log2.fold.change <- secDat$`T18_vs_T28_log2-fold difference`*-1
secDat_annot <- merge(secDat, annotD, by.x="GeneID", by.y="name", all.x=TRUE)[,c(1,2,57,53,58)]

#secProt <-secDat_annot[secDat_annot$GeneID %in%  secID$GeneID,]
deepProt <-secDat_annot[secDat_annot$GeneID %in%  deepLoc$name,]

# Find shared IDs (intersection)
shared_ids <- intersect(deepLoc_ids, deepProt_ids)
length(shared_ids)
# Find IDs unique to each dataframe
unique_to_deepLoc <- setdiff(deepLoc_ids, secProt_ids)
unique_to_deepProt <- setdiff(deepProt_ids, deepLoc_ids)
#write.table(secProt,"/home/pyre/Project_data/mittag_group/temperature_Transcriptome/2025FinalizedPaperJulyAugust/Supplemental_tables+updatedFigures/secDat539.txt",
#            sep="\t", quote = F, row.names = F)

#write.table(secDat_annot, "/home/pyre/Project_data/mittag_group/temperature_Transcriptome/2025FinalizedPaperJulyAugust/Supplemental_tables+updatedFigures/secDatComplete2763.txt",
#            sep="\t", quote = F, row.names = F)



##Get significant ids
secDat_annotSig <- secDat_annot[
  (secDat_annot$X18_vs_X28_log2.fold.change > 0.5 | 
     secDat_annot$X18_vs_X28_log2.fold.change < -0.5) &
    (secDat_annot$T18_vs_T28_p.adj < 0.05),]

secProtSig <- secProt[
  (secProt$X18_vs_X28_log2.fold.change > 0.5 | secProt$X18_vs_X28_log2.fold.change < -0.5) &
    (secProt$T18_vs_T28_p.adj < 0.05),]

write.table(secProtSig,"2025FinalizedPaperJulyAugust/secretedProts.csv",
           quote = F, row.names = F, sep="\t")
##Now prep the data to feed into enrichment analysis.
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
  geneList <- factor(as.integer(annotD$GeneID %in% deg$GeneID))
  names(geneList) <- annotD$GeneID
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

colnames(secDat_annotSig) <- c("GeneID",
                           "log2FoldChange",
                           "padj", "GeneDescription")
colnames(secProtSig) <- c("GeneID",
                               "log2FoldChange",
                               "padj", "GeneDescription")
#enrichment 
B_up <- goFunc(fTable = secProtSig,onto = "BP", opt="up")
M_up <- goFunc(fTable = secProtSig,onto = "MF", opt="up")
C_up <- goFunc(fTable = secProtSig,onto = "CC", opt="up")
GO_up<- rbind(B_up, M_up, C_up)
write.table(GO_up, "Up_GOenrich_18v28.txt",
            quote = F, sep="\t",
            row.names = F)

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

###enrichment TopGO for downregulated terms
B_down <- goFunc(fTable = secProtSig,onto = "BP", opt="down")
M_down <- goFunc(fTable = secProtSig,onto = "MF", opt="down")
C_down <- goFunc(fTable = secProtSig,onto = "CC", opt="down")
GO_down<- rbind(B_down, M_down, C_down)
write.table(GO_down, "Down_GOenrich_18v28.txt",
            quote = F, sep="\t",
            row.names = F)

###Get importantIds 
impIdsComp <- read.csv("~/Project_data/mittag_group/temperature_Transcriptome/Finalized_Complete_results/SecImpIds.csv", head=T)
MatingSec <- readxl::read_xlsx("~/Downloads/proteomes-06-00036-s001/SupplementalTable3_SolubleMatingSecretome_MassSpecData_Rev.xlsx")
head(impIdsComp)
head(MatingSec)

impIdsComp$GeneID_clean <- gsub("_4532$", "", impIdsComp$GeneID)

# Compare: which GeneIDs from impIdsComp exist in MatingSec
impIdsComp$In_MatingSec <- impIdsComp$GeneID_clean %in% MatingSec$GeneID

# To view the result:
mating <- impIdsComp[impIdsComp$In_MatingSec=="TRUE",]
mating[mating$GeneID%in%deepProt$GeneID,]


MatingSecretome <- merge(deepProt,impIdsComp, by="GeneID",all.y=T)


ECP33 <- resT18v33[resT18v33$name %in% impIdsComp$GeneID,c(2,6,7)]
ECP33$Temp <- "T33"

ECP28 <- resT18v28[resT18v28$name %in% impIdsComp$GeneID,c(2,6,7)]
ECP28$Temp <- "T28"


combined_df <- rbind(ECP33,ECP28)
combined_df <- merge(combined_df,impIdsComp,by.x="name", by.y="GeneID")

combined_df <- combined_df %>%
  mutate(
    label = case_when(   # Customize significance labels based on p-value
      padj < 0.00001 ~ "****",
      padj < 0.0001 ~ "***",  
      padj < 0.001  ~ "**",    
      padj < 0.01   ~ "*",  
    )
  )
# Convert Temp to factor after modification
combined_df$Temp <- as.factor(combined_df$Temp)
combined_df$Direction <- ifelse(as.numeric(combined_df$log2FoldChange) > 0, "Positive", "Negative")
combined_df$log2FoldChange <- as.numeric(abs(combined_df$log2FoldChange))


combined_df2 <- transform(combined_df, Gene = paste(name, Function, sep = " "))
range(combined_df$log2FoldChange)
table(combined_df2$Function)

symbols_to_keep <- combined_df2 %>%
  filter(log2FoldChange > 1) %>%
  filter(padj <0.01) %>%
  pull(name) %>%
  unique()

combined_df2_filtered <- combined_df2 %>%
  filter(name %in% symbols_to_keep)


ecpP <- ggplot(combined_df2_filtered, aes(x = log2FoldChange, y = reorder(Symbol,-Order), group = Temp, fill=Direction,
                                       pattern = Temp)) +
  geom_bar_pattern(stat = "identity", 
                   width = 0.75,    
                   position = position_dodge(0.75), 
                   colour = "black",
                   pattern_fill = "black",           # Color of the pattern
                   pattern_density = 0.2,           # Density of the pattern
                   pattern_spacing = 0.02,          # Spacing between pattern elements
                   pattern_angle = 45) +
  scale_fill_manual(values = c("Negative" = "royalblue", "Positive" = "darkorange2")) +
  geom_text(aes(label = label), vjust = 0.8, hjust=-0.1,size = 20, position = position_dodge(width = 0.9)) +
  theme_minimal()+ xlim(0,8.5)+
  labs(y= "Gene names", x="Log2 Foldchange")+
  facet_grid(Function ~ ., scales = "free_y", space = "free_y", switch = "y")+
  scale_pattern_manual(values = c(T28 = "none", T33 = "stripe")) +
  theme(axis.title.x = element_text(face = "bold", size=20), 
        axis.title.y = element_text(face = "bold", size=20),
        axis.text.x = element_text(face="bold", color="#000000", size=20),
        axis.text.y = element_text(face="bold", color="#000000", size=20),
        legend.title=element_text(face="bold", size=5),
        panel.grid.major = element_line(color = "grey", size=1), 
        strip.text.x = element_text(size = 18),
        strip.text = element_text(size=10),
        strip.placement = "outside")
ecpP


####
ecp_df <- summary_df[summary_df$name %in% impIdsComp$GeneID,]
ecp_dfAnnot <- merge(ecp_df, GeneSymb, by.x="name", by.y="ens_gene")
write.table(ecp_dfAnnot, "ecp_annotated_table.txt",
            sep="\t", quote = )

####
prot_Matches <- merge(secDat, protplotIds, by.x="GeneID", by.y="Name")


colnames(protplotIds)
combined_df <- prot_Matches %>%
  mutate(
    label = case_when(   # Customize significance labels based on p-value
      X18_vs_X28_p.adj < 0.0001 ~ "****",  
      X18_vs_X28_p.adj < 0.001  ~ "***",    
      X18_vs_X28_p.adj < 0.01   ~ "**",     
      X18_vs_X28_p.adj < 0.05   ~ "*", 
      X18_vs_X28_p.adj >= 0.05 ~ "",
    )
  )
# Convert Temp to factor after modification
combined_df <- combined_df[!(duplicated(combined_df$GeneID) & combined_df$label == ""), ]

combined_df$Direction <- ifelse(as.numeric(combined_df$X18_vs_X28_log2.fold.change) > 0, "Positive", "Negative")
combined_df$X18_vs_X28_log2.fold.change <- as.numeric(abs(combined_df$X18_vs_X28_log2.fold.change))

ggplot(combined_df, aes(x = X18_vs_X28_log2.fold.change, y = reorder(Symbol, -Order), fill = Direction)) +
  geom_bar(stat = "identity", width = 0.75, position = position_dodge(0.75), colour = "black") +
  scale_fill_manual(values = c("Negative" = "royalblue", "Positive" = "darkorange2")) +
  geom_text(aes(label = label), vjust = 0.8, hjust = -0.1, size = 10, position = position_dodge(width = 0.9)) +
  theme_minimal() +
  xlim(0, 7) +
  labs(y = "Gene names", x = "Log2 Foldchange") +
  theme(
    axis.title.x = element_text(face = "bold", size = 20),
    axis.title.y = element_text(face = "bold", size = 20),
    axis.text.x = element_text(face = "bold", color = "#000000", size = 20),
    axis.text.y = element_text(face = "bold", color = "#000000", size = 20),
    legend.title = element_text(face = "bold", size = 12),
    panel.grid.major = element_line(color = "grey", size = 1),
    strip.text.x = element_text(size = 18),
    strip.text = element_text(size = 10),
    strip.placement = "outside"
  )


prot_Matches <- merge(secDat, impIdsComp, by.x="GeneID", by.y="GeneID")
combined_df <- prot_Matches %>%
  mutate(
    label = case_when(   # Customize significance labels based on p-value
      T18_vs_T28_p.adj < 0.0001 ~ "****",  
      T18_vs_T28_p.adj < 0.001  ~ "***",    
      T18_vs_T28_p.adj < 0.01   ~ "**",     
      T18_vs_T28_p.adj < 0.05   ~ "*", 
      T18_vs_T28_p.adj >= 0.05 ~ "",
    )
  )
# Convert Temp to factor after modification
combined_df <- combined_df[!(duplicated(combined_df$GeneID) & combined_df$label == ""), ]

combined_df$Direction <- ifelse(as.numeric(combined_df$X18_vs_X28_log2.fold.change) > 0, "Positive", "Negative")
combined_df$X18_vs_X28_log2.fold.change <- as.numeric(abs(combined_df$X18_vs_X28_log2.fold.change))


ggplot(combined_df, aes(x = X18_vs_X28_log2.fold.change, y = reorder(Symbol, -Order), fill = Direction)) +
  geom_bar(stat = "identity", width = 0.75, position = position_dodge(0.75), colour = "black") +
  scale_fill_manual(values = c("Negative" = "royalblue", "Positive" = "darkorange2")) +
  geom_text(aes(label = label), vjust = 0.8, hjust = -0.1, size = 10, position = position_dodge(width = 0.9)) +
  theme_minimal() +
  xlim(0, 6) +
  labs(y = "Gene names", x = "Log2 Foldchange") +
  theme(
    axis.title.x = element_text(face = "bold", size = 20),
    axis.title.y = element_text(face = "bold", size = 20),
    axis.text.x = element_text(face = "bold", color = "#000000", size = 20),
    axis.text.y = element_text(face = "bold", color = "#000000", size = 20),
    legend.title = element_text(face = "bold", size = 12),
    panel.grid.major = element_line(color = "grey", size = 1),
    strip.text.x = element_text(size = 18),
    strip.text = element_text(size = 10),
    strip.placement = "outside"
  )



prot_Matches <- merge(secDat, impIdsComp, by.x="GeneID", by.y="GeneID")
combined_df <- prot_Matches %>%
  mutate(
    label = case_when(   # Customize significance labels based on p-value
      T18_vs_T28_p.adj < 0.0001 ~ "****",  
      T18_vs_T28_p.adj < 0.001  ~ "***",    
      T18_vs_T28_p.adj < 0.01   ~ "**",     
      T18_vs_T28_p.adj < 0.05   ~ "*", 
      T18_vs_T28_p.adj >= 0.05 ~ "",
    )
  )
# Convert Temp to factor after modification
combined_df <- combined_df[!(duplicated(combined_df$GeneID) & combined_df$label == ""), ]

combined_df$Direction <- ifelse(as.numeric(combined_df$X18_vs_X28_log2.fold.change) > 0, "Positive", "Negative")
combined_df$X18_vs_X28_log2.fold.change <- as.numeric(abs(combined_df$X18_vs_X28_log2.fold.change))

ggplot(combined_df, aes(x = X18_vs_X28_log2.fold.change, y = reorder(Symbol, -Order), fill = Direction)) +
  geom_bar(stat = "identity", width = 0.75, position = position_dodge(0.75), colour = "black") +
  scale_fill_manual(values = c("Negative" = "royalblue", "Positive" = "darkorange2")) +
  geom_text(aes(label = label), vjust = 0.8, hjust = -0.1, size = 10, position = position_dodge(width = 0.9)) +
  theme_minimal() +
  xlim(0, 12) +
  labs(y = "Gene names", x = "Log2 Foldchange") +
  theme(
    axis.title.x = element_text(face = "bold", size = 20),
    axis.title.y = element_text(face = "bold", size = 20),
    axis.text.x = element_text(face = "bold", color = "#000000", size = 20),
    axis.text.y = element_text(face = "bold", color = "#000000", size = 20),
    legend.title = element_text(face = "bold", size = 12),
    panel.grid.major = element_line(color = "grey", size = 1),
    strip.text.x = element_text(size = 18),
    strip.text = element_text(size = 10),
    strip.placement = "outside"
  )


# Define thresholds
FC <- 0.5  # Fold-change threshold
pvalAdj <- 0.05  # Adjusted p-value threshold

# Calculate the number of upregulated and downregulated genes
deg <- deepProt %>% 
  filter(T18_vs_T28_p.adj < pvalAdj) %>%
  summarize(
    up_deg = sum(X18_vs_X28_log2.fold.change > FC),
    down_deg = sum(X18_vs_X28_log2.fold.change < -FC)
  )

# Define custom colors for the points
keyvals <- ifelse(
  deepProt$X18_vs_X28_log2.fold.change < -FC & deepProt$T18_vs_T28_p.adj < 0.05, 'royalblue',
  ifelse(deepProt$X18_vs_X28_log2.fold.change > FC & deepProt$T18_vs_T28_p.adj < 0.05, 'darkorange2',
         'black')
)
keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'darkorange2'] <- 'Up regulated'
names(keyvals)[keyvals == 'black'] <- 'Not significant'
names(keyvals)[keyvals == 'royalblue'] <- 'Down regulated'


EnhancedVolcano(
  deepProt,
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
  ylim=c(0,6.5)
)


saltSec <- read.csv("SaltSecID.csv", head=T)
saltSecMatch <- secDat[secDat$GeneID %in% saltSec$CreSaltID,c(11,3,5)]
saltSecMatch <- merge(saltSecMatch, annotD, by="GeneID")
write.table(secProt, "/home/pyre/Project_data/mittag_group/temperature_Transcriptome/Finalized_Complete_results/SecretedProteome.txt",
            sep="\t", quote = F, row.names = F)

##Make table of extracellular identified proteins
ECPidSec <- read.csv("/home/pyre/Project_data/mittag_group/temperature_Transcriptome/Finalized_Complete_results/SecImpIdsExtended.csv")

ECP33 <- resT18v33[resT18v33$name %in% ECPidSec$GeneID,c(2,6,7)]
ECP33$Temp <- "T33"

ECP28 <- resT18v28[resT18v28$name %in% ECPidSec$GeneID,c(2,6,7)]
ECP28$Temp <- "T28"


combined_df <- rbind(ECP33,ECP28)
combined_df <- merge(combined_df,ECPidSec,by.x="name", by.y="GeneID")

combined_df <- combined_df %>%
  mutate(
    label = case_when(   # Customize significance labels based on p-value
      padj < 0.0001 ~ "****",  
      padj < 0.001  ~ "***",    
      padj < 0.01   ~ "**",     
      padj < 0.05   ~ "*", 
      padj >= 0.05 ~ "",
    )
  )
# Convert Temp to factor after modification
combined_df$Temp <- as.factor(combined_df$Temp)
combined_df$Direction <- ifelse(as.numeric(combined_df$log2FoldChange) > 0, "Positive", "Negative")
combined_df$log2FoldChange <- as.numeric(abs(combined_df$log2FoldChange))


combined_df2 <- transform(combined_df, Gene = paste(name, Function, sep = " "))
range(combined_df$log2FoldChange)
table(combined_df2$Function)

#Drop non sig prots
combined_df2Sig <- combined_df2[!combined_df2$Symbol == "MMP1",]

ecpP <- ggplot(combined_df2Sig, aes(x = log2FoldChange, y = reorder(Symbol,-Order), group = Temp, fill=Direction,
                                 pattern = Temp)) +
  geom_bar_pattern(stat = "identity", 
                   width = 0.75,    
                   position = position_dodge(0.75), 
                   colour = "black",
                   pattern_fill = "black",           # Color of the pattern
                   pattern_density = 0.2,           # Density of the pattern
                   pattern_spacing = 0.02,          # Spacing between pattern elements
                   pattern_angle = 45) +
  scale_fill_manual(values = c("Negative" = "royalblue", "Positive" = "darkorange2")) +
  geom_text(aes(label = label), vjust = 0.8, hjust=-0.1,size = 10, position = position_dodge(width = 0.9)) +
  theme_minimal()+xlim(0,8.5)+
  labs(y= "Gene names", x="Log2 Foldchange")+
  #facet_grid(Function ~ ., scales = "free_y", space = "free_y", switch = "y")+
  scale_pattern_manual(values = c(T28 = "none", T33 = "stripe")) +
  theme(axis.title.x = element_text(face = "bold", size=20), 
        axis.title.y = element_text(face = "bold", size=20),
        axis.text.x = element_text(face="bold", color="#000000", size=20),
        axis.text.y = element_text(face="bold", color="#000000", size=20),
        legend.title=element_text(face="bold", size=5),
        panel.grid.major = element_line(color = "grey", size=1), 
        strip.text.x = element_text(size = 18),
        strip.text = element_text(size=10),
        strip.placement = "outside")
ecpP


ECPidSec <- read.csv("/home/pyre/Project_data/mittag_group/temperature_Transcriptome/Finalized_Complete_results/SecImpIdsExtended.csv")
prot_Matches <- merge(deepProt, ECPidSec, by.x="GeneID", by.y="GeneID")
merge(prot_Matches, deepLoc, by.x="GeneID", by.y="name")

combined_df <- prot_Matches %>%
  mutate(
    label = case_when(   # Customize significance labels based on p-value
      T18_vs_T28_p.adj < 0.0001 ~ "****",  
      T18_vs_T28_p.adj < 0.001  ~ "***",    
      T18_vs_T28_p.adj < 0.01   ~ "**",     
      T18_vs_T28_p.adj < 0.05   ~ "*", 
      T18_vs_T28_p.adj >= 0.05 ~ "",
    )
  )
# Convert Temp to factor after modification
combined_df <- combined_df[!(duplicated(combined_df$GeneID) & combined_df$label == ""), ]

combined_df$Direction <- ifelse(as.numeric(combined_df$X18_vs_X28_log2.fold.change) > 0, "Positive", "Negative")
combined_df$X18_vs_X28_log2.fold.change <- as.numeric(abs(combined_df$X18_vs_X28_log2.fold.change))
combined_dfSub <- combined_df[combined_df$T18_vs_T28_p.adj < 0.05,]
combined_dfSub <- combined_dfSub[!duplicated(combined_dfSub$Symbol), ]

combined_df2_sub_filtered <- combined_dfSub %>%
  group_by(Symbol) %>%
  filter(any(X18_vs_X28_log2.fold.change > 0.5)) %>%
  ungroup()

ggplot(combined_dfSub, aes(x = X18_vs_X28_log2.fold.change, y = reorder(Symbol, -Order), fill = Direction)) +
  geom_bar(stat = "identity", width = 0.75, position = position_dodge(0.75), colour = "black") +
  scale_fill_manual(values = c("Negative" = "royalblue", "Positive" = "darkorange2")) +
  geom_text(aes(label = label), vjust = 0.8, hjust = -0.1, size = 10, position = position_dodge(width = 0.9)) +
  theme_minimal() +
  xlim(0, 12) +
  labs(y = "Gene names", x = "Log2 Foldchange") +
  theme(
    axis.title.x = element_text(face = "bold", size = 20),
    axis.title.y = element_text(face = "bold", size = 20),
    axis.text.x = element_text(face = "bold", color = "#000000", size = 20),
    axis.text.y = element_text(face = "bold", color = "#000000", size = 20),
    legend.title = element_text(face = "bold", size = 12),
    panel.grid.major = element_line(color = "grey", size = 1),
    strip.text.x = element_text(size = 18),
    strip.text = element_text(size = 10),
    strip.placement = "outside"
  )

####Correlation analysis
# Rename transcriptome columns for consistency
colnames(resT18v28_desc)[1] <- "GeneID"
colnames(resT18v28_desc)[2] <- "X18_vs_X28_log2.fold.change_trans"

# Extract relevant columns
prot <- protDat[, c("GeneID", "X18_vs_X28_log2.fold.change")]
sec <- secProt[, c("GeneID", "X18_vs_X28_log2.fold.change")]
trans <- resT18v28_desc[, c("GeneID", "X18_vs_X28_log2.fold.change_trans")]

# Merge dataframes
mergedSecFla <- merge(prot, sec, by = "GeneID", suffixes = c("_prot", "_sec"))
mergedSecTrans <- merge(sec, trans, by = "GeneID")
mergedSecFla <- merge(mergedSecFla, annotD, by="GeneID")

cor_prot_sec <- cor(mergedSecFla$X18_vs_X28_log2.fold.change_sec, 
                    mergedSecFla$X18_vs_X28_log2.fold.change_prot, 
                    use = "complete.obs")

cor_prot_trans <- cor(merged$X18_vs_X28_log2.fold.change_sec, 
                      merged$X18_vs_X28_log2.fold.change_trans, 
                      use = "complete.obs")

# Proteome vs. Transcriptome
ggplot(mergedSecFla, aes(x = X18_vs_X28_log2.fold.change_sec, y = X18_vs_X28_log2.fold.change_prot)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Proteome vs. Transcriptome log2FC", x = "Proteome log2FC", y = "Transcriptome log2FC")

