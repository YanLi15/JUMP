rm(list = ls())
load("Data analysis/Mouse cerebellum (Slide-seq + Slide-seqV2)/MC_crop50.RData")
load("Data analysis/Mouse cerebellum (Slide-seq + Slide-seqV2)/MC_Results_crop50.RData")

##----------------------------------------------------------
## Moran's I
##----------------------------------------------------------
source("R/ValidationFuncs.R")
library(spdep)

genes.overlap <- intersect(intersect(genes_rep_bh, genes_rep_maxp), genes_rep_jump)
genes.bh.only  <- genes_rep_bh[!genes_rep_bh%in%genes.overlap]
genes.jump.only <- genes_rep_jump[!genes_rep_jump%in%genes.overlap]

# MI based on Slide-seq data
stat.res1 <- corrValue(counts1, location1)
MI1 <- stat.res1$MI
mi_all1 <- MI1[overlap]

mi_bh_only1 <- MI1[genes.bh.only]
mi_jump_only1 <- MI1[genes.jump.only]

# Moran's I for SVGs additionally identified by different methods
mean(mi_all1); mean(mi_jump_only1); median(mi_all1); median(mi_jump_only1)
# 0.0011; 0.0079; 1.77e-4; 0.0061
mean(mi_bh_only1); median(mi_bh_only1)
# 0.0083; 0.0068
factor1 <- factor(rep(c("JUMP only", "All (Slide-seq)"), 
                      times = c(length(genes.jump.only), length(overlap))))
dataset1 <- data.frame(value = c(mi_jump_only1, mi_all1), group = factor1)
par(mar = c(3, 4.5, 1, 1))
boxplot(value ~ group, dataset1, col = c("antiquewhite1", "lightSalmon"), 
        xlab = NULL, ylab = "Moran's I statistics",outline = FALSE,
        cex.axis = 1.4, cex.lab = 1.6) # 4.5 * 8 inches

# MI based on Slide-seqV2 data
stat.res2 <- corrValue(counts2, location2)
MI2 <- stat.res2$MI
mi_all2 <- MI2[overlap]

mi_bh_only2 <- MI2[genes.bh.only]
mi_jump_only2 <- MI2[genes.jump.only]

# Moran's I for different methods
mean(mi_all2); mean(mi_maxp2); mean(mi_bh2); mean(mi_jump2)
median(mi_all2); median(mi_maxp2); median(mi_bh2); median(mi_jump2)

factor2 <- factor(rep(c("All (Slide-seqV2)", "BH", "MaxP-BH", "MaxP-ST"),
                     times = c(length(overlap),length(genes_rep_bh),
                               length(genes_rep_maxp),
                               length(genes_rep_jump))))
dataset2 <- data.frame(value = c(mi_all2, mi_bh2, mi_maxp2, mi_jump2),
                       group = factor2)
par(mar = c(3, 4.5, 1, 1))
boxplot(value ~ group, dataset2,
        col = c("antiquewhite1", "#56B4E9", "#8FBC8F", "#F4B183"),
        xlab = NULL, ylab = "Moran's I statistics",outline = FALSE,
        cex.axis = 1.4, cex.lab = 1.6)


# Moran's I for SVGs additionally identified by different methods
mean(mi_all2); mean(mi_jump_only2); median(mi_all2); median(mi_jump_only2)
# 0.0031; 0.0176; 0.0011; 0.0143
mean(mi_bh_only2); median(mi_bh_only2)
# 0.0185; 0.0152
factor2 <- factor(rep(c("JUMP only", "All (Slide-seqV2)"), 
                      times = c(length(genes.jump.only), length(overlap))))
dataset2 <- data.frame(value = c(mi_jump_only2, mi_all2), group = factor2)
par(mar = c(3, 4.5, 1, 1))
boxplot(value ~ group, dataset2, col = c("antiquewhite1", "lightSalmon"), 
        xlab = NULL, ylab = "Moran's I statistics",outline = FALSE,
        cex.axis = 1.4, cex.lab = 1.7) # 4.8 * 8.4

save(stat.res1, stat.res2, file = "Data analysis/Mouse cerebellum (Slide-seq + Slide-seqV2)/MC_MI.RData")

##--------------------------------------------------------------------
## Summarized spatial expression patterns
##--------------------------------------------------------------------
library(amap)
source("R/PlotFuncs.R")

## Summarized patterns in Slide-seq data
vst_count1 <- var_stabilize(counts1) # R function in funcs.R
sig_vst_count1 <- vst_count1[genes_rep_jump, ]
sig_vst_res1 <- t(apply(sig_vst_count1, 1, LMReg, T = log(apply(counts1, 2, sum))))

hc1 <- hcluster(sig_vst_res1, method = "euc", link = "ward", nbproc = 1, 
                doubleprecision = TRUE)
numC <- 3
memb1 <- cutree(hc1, k = numC)
cent1 <- NULL
for (k in 1:numC) {
  cent1 <- cbind(cent1, colMeans(sig_vst_res1[memb1 == k, , drop = FALSE]))
}
position_cord1 <- location1
rownames(position_cord1) <- rownames(cent1)

rel_cent1 <- t(apply(cent1, 1, relative_func))
# rel_cent1 <- rel_cent1[,c(2,1,3)]
pd1 <- setNames(cbind.data.frame(position_cord1, rel_cent1), 
                c("x", "y", paste0("Pattern ", c("III", "II", "I"))))
MBP1 <- lapply(1:numC, function(x) {
  pattern_plot3(pd1, x, xy = T, main = T, titlesize = 9.5)
})
grid.arrange(grobs = MBP1[numC:1], nrow = 3) # 12.5 * 37.5

## Summarized patterns in Slide-seqV2 data
vst_count2 <- var_stabilize(counts2) 
sig_vst_count2 <- vst_count2[genes_rep_jump, ]
sig_vst_res2 <- t(apply(sig_vst_count2, 1, LMReg, T = log(apply(counts2, 2, sum))))

hc2 <- hcluster(sig_vst_res2, method = "euc", link = "ward", nbproc = 1, 
                doubleprecision = TRUE)
numC <- 3
memb2 <- cutree(hc2, k = numC)
cent2 <- NULL
for (k in 1:numC) {
  cent2 <- cbind(cent2, colMeans(sig_vst_res2[memb2 == k, , drop = FALSE]))
}
position_cord2 <- location2
rownames(position_cord2) <- rownames(cent2)

rel_cent2 <- t(apply(cent2, 1, relative_func))
rel_cent2 <- rel_cent2[,c(3,2,1)]
pd2 <- setNames(cbind.data.frame(position_cord2, rel_cent2), 
                c("x", "y", paste0("Pattern ", c("III", "II", "I"))))
MBP2 <- lapply(1:numC, function(x) {
  pattern_plot3(pd2, x, xy = T, main = T, titlesize = 9.5)
})
grid.arrange(grobs = MBP2[numC:1], nrow = 3) # 12.5 * 37.5

##--------------------------------------------------------------------
## Spatial expression patterns of randomly selected SVGs
##--------------------------------------------------------------------
source("R/PlotFuncs.R")

# All genes identified by MaxP_S
randi <- sample(length(genes_rep_jump), 24, replace = FALSE)
gene_plot <- genes_rep_jump[randi]

# three representative genes corresponding to the main patterns
gene_plot <- c("Pcp2", "Mbp", "Snap25")
# three representative genes uniquely identified by MaxP-ST corresponding to the main patterns
# gene_plot <- c("Ppp1r17", "Grb14", "Ndufb10")
randi <- sample(length(genes.jump.only), 24, replace = FALSE)
gene_plot <- genes.jump.only[randi]
gene_plot <- c('Mrps12', 'Ndufb11', 'Gm14033', 'Ndufa3', 'Golga7b', 'Gsn', 'Celf2',
               'Banp', 'Aplp2', 'Rpl35a', 'Atrx', 'Gria2', 'Tuba4a', 'Homer3',
               'Cox5a', 'Atp5f1', 'Npy', 'Slc6a1', 'Pllp', 'Snrpn',
               'Pebp1', 'Dnm1', 'Psd3', 'Trim9')

genes.jump.unique <- genes_rep_jump[!genes_rep_jump%in%genes_rep_bh]
randi <- sample(length(genes.jump.unique), 10, replace = FALSE)
gene_plot <- genes.jump.unique[randi]
# 'Max', 'Map4', 'Ndufb10', 'Gdf10', 'Bod1l', 'Tuba1b', 'Thy1', 'Ndufb11', 'Trim9',
# 'Txndc15', 'Nfix', 'Uqcrh', 'Dclk1', 'Eef1b2', 'Hsp90aa1', 'Atpif1', '2010107E04Rik'

# Based on Slide-seq data
vst_count1 <- var_stabilize(counts1) # R function in funcs.R
sig_vst_ct1 <- vst_count1[gene_plot, ]
rel_vst_ct1 <- apply(sig_vst_ct1, 1, relative_func)
pltdat1 <- cbind.data.frame(location1,rel_vst_ct1)
genetitle <- gene_plot
pp1 <- lapply(1:(ncol(pltdat1)-2),
              function(x){pattern_plot3(pltdat1,x,main=T,titlesize=9.5,title=genetitle[x])})
grid.arrange(grobs=pp1, nrow=3) # 37.5 * 100 inches

##--------------------------------------------------------------------
vst_count2 <- var_stabilize(counts2) # R function in funcs.R
sig_vst_ct2 <- vst_count2[gene_plot, ]
rel_vst_ct2 <- apply(sig_vst_ct2, 1, relative_func)
pltdat2 <- cbind.data.frame(location2,rel_vst_ct2)
genetitle <- gene_plot
pp2 <- lapply(1:(ncol(pltdat2)-2),
              function(x){pattern_plot3(pltdat2,x,main=T,titlesize=9.5,title=genetitle[x])})
grid.arrange(grobs=pp2, nrow=3) # 37.5 * 100 inches

##--------------------------------------------------------------------
## Mouse cerebellum related studies validation
##--------------------------------------------------------------------
# library(readxl)
# # Cell type marker genes (4152) (Wizeman et al., 2019)
# cellmarker_genes <- read_xlsx("Data analysis/Mouse cerebellum (Slide-seq + Slide-seqV2)/Validation sets/elife-42388-supp1-v2.xlsx",
#                               col_names = TRUE)
# 
# cellmarker_genes <- unique(cellmarker_genes$gene) # 4152 genes
# 
# length(intersect(genes_rep_maxp, cellmarker_genes)) # crop50: 152/279; crop25: 172/316
# length(intersect(genes.bh.only, cellmarker_genes)) # 54/115; 55/128
# length(intersect(genes.jump.only, cellmarker_genes)) # 78/169; 71/169

# Genes related to the cerebellum in Harmonizome database (1867)
library(rjson)
library(stringr)
## Allen Brain Atlas & BioGPS mouse cell type
Harmonizome_Allen <- fromJSON(paste(readLines("Data analysis/Mouse cerebellum (Slide-seq + Slide-seqV2)/Validation sets/Harmonizome_Allen Atlas.json")))
Harmonizome_Allen_cortex <- fromJSON(paste(readLines("Data analysis/Mouse cerebellum (Slide-seq + Slide-seqV2)/Validation sets/Harmonizome_Allen_cortex.json")))
Harmonizome_Allen_hemisphere <- fromJSON(paste(readLines("Data analysis/Mouse cerebellum (Slide-seq + Slide-seqV2)/Validation sets/Harmonizome_Allen_hemisphere.json")))
# Harmonizome_Allen_nuclei <- fromJSON(paste(readLines("Data analysis/Mouse cerebellum (Slide-seq + Slide-seqV2)/Validation sets/Harmonizome_Allen_nuclei.json")))
# Harmonizome_Allen_vermis <- fromJSON(paste(readLines("Data analysis/Mouse cerebellum (Slide-seq + Slide-seqV2)/Validation sets/Harmonizome_Allen_vermis.json")))
# Harmonizome_Allen_white  <- fromJSON(paste(readLines("Data analysis/Mouse cerebellum (Slide-seq + Slide-seqV2)/Validation sets/Harmonizome_Allen_white matter.json")))
# Harmonizome_BioGPS <- fromJSON(paste(readLines("Data analysis/Mouse cerebellum (Slide-seq + Slide-seqV2)/Validation sets/Harmonizome_BioGPS.json")))
Harmonizome_Allen <- Harmonizome_Allen[["associations"]]
Harmonizome_Allen_cortex <- Harmonizome_Allen_cortex[["associations"]]
Harmonizome_Allen_hemisphere <- Harmonizome_Allen_hemisphere[["associations"]]
# Harmonizome_Allen_nuclei <- Harmonizome_Allen_nuclei[["associations"]]
# Harmonizome_Allen_vermis <- Harmonizome_Allen_vermis[["associations"]]
# Harmonizome_Allen_white <- Harmonizome_Allen_white[["associations"]]

# Harmonizome_BioGPS <- Harmonizome_BioGPS[["associations"]]
Harmonizome.genes <- NULL
for (i in 1:length(Harmonizome_Allen)){
  Harmonizome.genes   <- c(Harmonizome.genes, Harmonizome_Allen[[i]][["gene"]][["symbol"]])
}
for (i in 1:length(Harmonizome_Allen_cortex)){
  Harmonizome.genes   <- c(Harmonizome.genes, Harmonizome_Allen_cortex[[i]][["gene"]][["symbol"]])
}
for (i in 1:length(Harmonizome_Allen_hemisphere)){
  Harmonizome.genes   <- c(Harmonizome.genes, Harmonizome_Allen_hemisphere[[i]][["gene"]][["symbol"]])
}
# for (i in 1:length(Harmonizome_Allen_nuclei)){
#   Harmonizome.genes   <- c(Harmonizome.genes, Harmonizome_Allen_nuclei[[i]][["gene"]][["symbol"]])
# }
# for (i in 1:length(Harmonizome_Allen_vermis)){
#   Harmonizome.genes   <- c(Harmonizome.genes, Harmonizome_Allen_vermis[[i]][["gene"]][["symbol"]])
# }
# for (i in 1:length(Harmonizome_Allen_white)){
#   Harmonizome.genes   <- c(Harmonizome.genes, Harmonizome_Allen_white[[i]][["gene"]][["symbol"]])
# }
# for (i in 1:length(Harmonizome_BioGPS)){
#   Harmonizome.genes   <- c(Harmonizome.genes, Harmonizome_BioGPS[[i]][["gene"]][["symbol"]])
# }

Harmonizome.genes <- Harmonizome.genes[!duplicated(Harmonizome.genes)]
Harmonizome.genes <- tolower(Harmonizome.genes)
Harmonizome.genes <- str_to_title(Harmonizome.genes) # 3000 genes

length(intersect(genes_rep_maxp, Harmonizome.genes)) # 84/279; 76/316
length(intersect(genes.bh.only, Harmonizome.genes)) # 28/115; 23/128
length(intersect(genes.jump.only, Harmonizome.genes)) # 41/169; 30/169

# # spatially non-random genes in the Purkinje layer identified from Slide-seq data (669)
# SlideSeq.genes <- read.table("Data analysis/Mouse cerebellum (Slide-seq + Slide-seqV2)/Validation sets/Slide-seq.txt")
# SlideSeq.genes <- t(SlideSeq.genes) # 669
# 
# length(intersect(genes_rep_maxp, SlideSeq.genes)) # 91/279; 101/308
# length(intersect(genes.bh.only, SlideSeq.genes)) # 28/115; 23/118
# length(intersect(genes.jump.only, SlideSeq.genes)) # 39/169; 29/149

# Kozareva et al.
library(readxl)

cluster_genes <- read_xlsx("Data analysis/Mouse cerebellum (Slide-seq + Slide-seqV2)/Validation sets/Kozareva et al.xlsx",
                              col_names = TRUE)
cluster_genes = cluster_genes[which(abs(cluster_genes$logFC)>=0.5),] 
cluster_genes = as.vector(unique(cluster_genes$gene)) # 3976

length(intersect(genes_rep_maxp, cluster_genes)) # crop50: 214/279; 166/316
length(intersect(genes.bh.only, cluster_genes)) # 73/115; 49/128
length(intersect(genes.jump.only, cluster_genes)) # 107/169; 58/169

##--------------------------------------------------------------------
## GO enrichment
##--------------------------------------------------------------------
library(clusterProfiler)
library(org.Mm.eg.db)
library(dplyr)
library(ggplot2)
library(ggrepel)
# options(connectionObserver = NULL) # run if library(org.Mm.eg.db) failed

genes.all <- bitr(overlap, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Mm.eg.db)
# jump.only <- genes_rep_jump[!genes_rep_jump%in%genes_rep_bh]
genes_jump_only <- bitr(genes.jump.only, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Mm.eg.db)
go.jump.only <- enrichGO(
  gene          = genes_jump_only$ENTREZID,
  universe      = genes.all$ENTREZID,
  OrgDb         = org.Mm.eg.db,
  ont           = "All" ,
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.99,
  qvalueCutoff  = 0.99,
  readable      = TRUE,
  pool = TRUE
)
sum(go.jump.only$p.adjust<0.05) # 179

library(clusterProfiler)
library(org.Mm.eg.db)
library(dplyr)
library(ggplot2)
library(ggrepel)
# options(connectionObserver = NULL) # run if library(org.Mm.eg.db) failed

genes.jump <- bitr(genes_rep_jump, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Mm.eg.db)
go.jump <- enrichGO(
  gene          = genes.jump$ENTREZID,
  universe      = genes.all$ENTREZID,
  OrgDb         = org.Mm.eg.db,
  ont           = "All" ,
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.99,
  qvalueCutoff  = 0.99,
  readable      = TRUE,
  pool = TRUE
)
sum(go.jump$p.adjust<0.05) # 725(50); 696(30); 723(25); 708(40)
sum(go.jump$p.adjust<0.01) # 452(50); 431; 472; 427

genes.bh <- bitr(genes_rep_bh, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Mm.eg.db)
go.bh <- enrichGO(
  gene          = genes.bh$ENTREZID,
  universe      = genes.all$ENTREZID,
  OrgDb         = org.Mm.eg.db,
  ont           = "All" ,
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.99,
  qvalueCutoff  = 0.99,
  readable      = TRUE,
  pool = TRUE
)
sum(go.bh$p.adjust<0.05) # 688(50); 683; 727; 672
sum(go.bh$p.adjust<0.01) # 418(50); 421; 447; 403

genes.maxp <- bitr(genes_rep_maxp, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Mm.eg.db)
go.maxp <- enrichGO(
  gene          = genes.maxp$ENTREZID,
  universe      = genes.all$ENTREZID,
  OrgDb         = org.Mm.eg.db,
  ont           = "All" ,
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.99,
  qvalueCutoff  = 0.99,
  readable      = TRUE,
  pool = TRUE
)
sum(go.maxp$p.adjust<0.05) # 517(50); 613; 663; 622
sum(go.maxp$p.adjust<0.01) # 282(50); 359; 419; 356

sig.go.jump <- filter(go.jump, p.adjust<.01)
sig.go.bh <- filter(go.bh, p.adjust<.01)
sig.go.maxp <- filter(go.maxp, p.adjust<.01)

library(VennDiagram)
venn.diagram(
  x = list(sig.go.jump@result$ID, sig.go.bh@result$ID, sig.go.maxp@result$ID),
  category.names = c("JUMP", "BH", "MaxP"),
  fill = c("#F4B183", "#56B4E9", "#8FBC8F"),
  filename = 'Data analysis/Mouse cerebellum (Slide-seq + Slide-seqV2)/venn_go_MC_crop50.tiff',
  output = TRUE,
  margin = 0.02,
  cex = 1.5,
  cat.cex = 1.5,
  cat.dist = c(0.07, 0.05, 0.05),
)

sig.go.overlap <- intersect(sig.go.bh@result$ID, sig.go.jump@result$ID) # 394
go.jump.only <- sig.go.jump[!sig.go.jump@result$ID%in%sig.go.overlap] # 58

write.csv(go.jump,file = "./Data analysis/Mouse cerebellum (Slide-seq + Slide-seqV2)/MC_go_jump.csv", quote = FALSE)
write.csv(go.jump.only,file = "./Data analysis/Mouse cerebellum (Slide-seq + Slide-seqV2)/MC_go_jump_only.csv", quote = FALSE)

results <- go.jump@result
results <- results[sample(nrow(results)),]
results <- results[,c("ID","ONTOLOGY", "pvalue", "Count", "Description")]
results[,"Category"] = NA
results$Category[which(results$ONTOLOGY == "BP")] = 1
results$Category[which(results$ONTOLOGY == "CC")] = 2
results$Category[which(results$ONTOLOGY == "MF")] = 3
results <- results[order(results$Category),]

results[,"BP"] <- NA
results$BP[which(results$Category == 1)] = 1:sum(results$Category == 1)
results$BP[which(results$Category == 2)] = 1:sum(results$Category == 2)
results$BP[which(results$Category == 3)] = 1:sum(results$Category == 3)

don <- results %>% 
  
  # Compute group size
  group_by(ONTOLOGY) %>% 
  summarise(len=max(BP)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(len)-len) %>%
  dplyr::select(-len) %>%
  
  # Add this info to the initial dataset
  left_join(results, ., by=c("ONTOLOGY"="ONTOLOGY")) %>%
  
  # Add a cumulative position of group
  arrange(ONTOLOGY, BP) %>%
  mutate(BPcum=BP+tot)

axisdf = don %>% group_by(ONTOLOGY) %>% summarize(center=(max(BPcum) + min(BPcum))/2)
# MaxP_S
annotated = c("glutamatergic synapse", "neurotransmitter transport", "synapse organization",
              "exocytic process", "parallel fiber to Purkinje cell synapse",
              "GABA-ergic synapse", "regulation of neuron projection development")
#,
#              "regulation of nervous system process")
# # MaxP_S only
# annotated = c("regulation of neurotransmitter levels", "neuron projection terminus",
#               "regulation of trans-synaptic signaling", "locomotory behavior",
#               "intrinsic component of synaptic vesicle membrane", "myelin sheath")

# q-value 1e-6 = p-value 4e-8; q-value 0.05 = p-value 0.0183
ggplot(don, aes(x = BPcum, y=-log10(pvalue))) +
  geom_point(aes(color = as.factor(ONTOLOGY), size = Count), alpha=0.8) +
  scale_colour_manual(name="",  
                      values = c("BP"="Tan", "CC"="DarkSeaGreen", "MF"="#6496D2")) +
  scale_x_continuous(label = axisdf$ONTOLOGY, breaks= axisdf$center) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 40)) +     # remove space between plot area and x axis
  # geom_hline(yintercept = -log10(0.0097), color = '#545454', size= 1.2, linetype = "dashed") + 
  geom_hline(yintercept = -log10(0.0012), color = '#545454', size= 1.2, linetype = "dashed") + 
  guides(color = "none") + 
  theme_bw() +
  theme( 
    legend.position = c(0.02,0.98),
    legend.justification = c(0.02,0.98),
    panel.border = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 13),
    axis.text.x = element_text(size = 12.5),
    axis.text.y = element_text(size = 11),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(size = 0.5, linetype = "solid")
  ) + 
  geom_text_repel(
    data = don[don$Description %in% annotated,],
    aes(label = Description),
    size = 4,
    segment.color = "black", show.legend = FALSE) # 5 * 7.5 inches

##--------------------------------------------------------------------
## KEGG enrichment analysis
##--------------------------------------------------------------------
# install.packages('R.utils')
R.utils::setOption( "clusterProfiler.download.method",'auto')
library(cowplot)

# kegg.jump.only <- enrichKEGG(
#   gene         = genes_jump_only$ENTREZID,
#   organism     = 'mmu',
#   universe     = genes.all$ENTREZID,
#   pvalueCutoff = 0.9,
#   qvalueCutoff = 0.9
# )
# sum(kegg.jump.only$p.adjust<0.05) # 43
# sig.jump.only <- filter(kegg.jump.only, p.adjust<.05)

kegg.jump <- enrichKEGG(
  gene         = genes.jump$ENTREZID,
  organism     = 'mmu',
  universe     = genes.all$ENTREZID,
  pvalueCutoff = 0.9,
  qvalueCutoff = 0.9
)
sum(kegg.jump$p.adjust<0.05) # 61


kegg.bh <- enrichKEGG(
  gene         = genes.bh$ENTREZID,
  organism     = 'mmu',
  universe     = genes.all$ENTREZID,
  pvalueCutoff = 0.9,
  qvalueCutoff = 0.9
)
sum(kegg.bh$p.adjust<0.05) # 58


kegg.maxp <- enrichKEGG(
  gene         = genes.maxp$ENTREZID,
  organism     = 'mmu',
  universe     = genes.all$ENTREZID,
  pvalueCutoff = 0.9,
  qvalueCutoff = 0.9
)
sum(kegg.maxp$p.adjust<0.05) # 48


sig.kegg.jump <- filter(kegg.jump, p.adjust<.05)
sig.kegg.bh <- filter(kegg.bh, p.adjust<.05)
kegg.overlap    <- intersect(sig.kegg.jump@result$ID, sig.kegg.bh@result$ID)
kegg.jump.only  <- sig.kegg.jump[!sig.kegg.jump@result$ID%in%kegg.overlap]

sig.jump.only <- sig.kegg.jump %>% filter(ID %in% kegg.jump.only$ID)

# kegg.bar <- barplot(sig.jump.only,showCategory=20,color = "pvalue") +
#   theme(
#     legend.position = c(0.98,0.02),
#     legend.justification = c(0.98,0.02),
#     axis.text.y = element_text(size = 15),
#     axis.text.x = element_text(size = 12), # face = "bold",
#     axis.title = element_text(size = 15),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank()
#   ) + scale_color_continuous()

# par(mar = c(4, 4.5, 0.5, 1))
# kegg.dot <- dotplot(sig.jump.only,showCategory=20,color = "pvalue") +
#   theme(
#     legend.position = c(0.98,0.02),
#     legend.justification = c(0.98,0.02),
#     axis.text.y = element_text(size = 15),
#     axis.text.x = element_text(size = 12), # face = "bold",
#     axis.title = element_text(size = 15),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank()
#   ) + scale_color_continuous()

par(mar = c(4, 4.5, 0.5, 1))
dotplot(sig.kegg.jump,showCategory=10,color = "pvalue") +
  theme(
    legend.position = c(0.85,0.32),
    axis.text.y = element_text(size = 15),
    axis.text.x = element_text(size = 12), # face = "bold",
    axis.title = element_text(size = 15),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) + scale_color_continuous() # 6 * 9 inches


