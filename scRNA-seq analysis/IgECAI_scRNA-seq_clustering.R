library(Seurat)
library(patchwork)
library(ggplot2)
library(dplyr)
library(cowplot)
library(stringr)
library(enrichplot)
library(data.table)
library(clusterProfiler)
library(org.Mm.eg.db)
library(ReactomePA)
library(slingshot)
library(circlize)
library(ComplexHeatmap)
library(tradeSeq)
library(BiocParallel)
library(magrittr)

# installation of rDBEC package from Github: https://github.com/s-shichino1989/TASSeq-paper)
# is required for the reproduction of results.
library(rDBEC)

#Distibution-based background subtration,hashtag/sampletag demultiplexing, and Seurat Object creation for 
# BD Rhapsody-derived WTA data was conducted by using rDBEC package as previously described 
# Shichino, S. et al. Commun Biol. (2022). doi:10.1038/s42003-022-03536-0
# Github: https://github.com/s-shichino1989/TASSeq-paper)

#=================Background subtraction by DEBC filtering and create Seurat object===============================#
#conducted under Seurat version 2 workflow
setwd("/media/owner/Data/GEO submission fiiles/geo_submission_IgECAI/processed_file")
# read processed count data matrix ("matrix_inflection_MiyakeWTA2.txt.gz")
fnames = dir(pattern = "matrix_inflection_")
tablelist = lapply(fnames, tableread_fast_sparse)

fnames = gsub("\\..+$", "", fnames) # remove ".txt" from file names
fnames = gsub("matrix_inflection_", "", fnames)
gc()

names(tablelist)=fnames
DBEC_filter = background_subtraction_Biex(tablelist, min.event=100, minimum.max.expr=7,
                                          species="mmu", min.ave=6.0, AutoThreshold = FALSE,
                                          min.diff=5.5, modelnames="E",
                                          uncert.thre=1, nthreads=12, sample.name=fnames)

names(DBEC_filter) = fnames

DBEC_res = apply_DBEC_filter(tablelist, DBEC_filter=DBEC_filter, nthreads=24, sample.name = fnames)
names(DBEC_res) = fnames

save(DBEC_res, file="DBEC_Miyake.rda")
save(DBEC_filter, file="DBEC_filter_Miyake.rda")

##create Seurat object (Seurat v2 workflow) and annotate by hashtag
for (i in 1:length(DBEC_res)){
  colnames(DBEC_res[[i]]) = paste(names(DBEC_res)[i], colnames(DBEC_res[[i]]), sep='_')
}

seu = lapply(DBEC_res[[1]], CreateSeuratObject, min.cells = 5, min.genes = 500)

#Add hashtag annotation metadata
hashtag_fnames = dir(pattern = "Hashtag_top1M_")
hashtag_data=lapply(hashtag_fnames, read.table, header=TRUE, sep="\t", quote="", stringsAsFactors=F)
for (i in 1:length(hashtag_data)){
  hashtag_data[[i]][,1] = paste(names(tablelist)[i], hashtag_data[[i]][,1], sep='_')
}

for (i in 1:length(hashtag_data)){
  seu[[i]] = Demultiplex_DNAtags(hashtag_data[[i]], seu[[i]],
                                 nTags=4, nCells_in_cartridge=20000, 
                                 sample.name=fnames[[i]], scale.factor = 4001613,
                                 dir.name="./")
}


dir.name=("./")
sample.name=("Miyake_")

mBC = seu[[1]]

#add metadata
colnames(mBC@meta.data)[2]="nReads"

mito.genes = grep(pattern = "^mt.", x = rownames(x = mBC@data), value = TRUE)
percent.mito = Matrix::colSums(mBC@raw.data[mito.genes, ])/Matrix::colSums(mBC@raw.data)
mBC = AddMetaData(object = mBC, metadata = percent.mito, col.name = "percent.mito")

ribo.genes = grep(pattern = "^Rp[sl][[:digit:]]", x = rownames(x = mBC@data), value = TRUE)
percent.ribo = Matrix::colSums(mBC@raw.data[ribo.genes, ])/Matrix::colSums(mBC@raw.data)
mBC = AddMetaData(object = mBC, metadata = percent.ribo, col.name = "percent.ribo")

nReads_log = log10(Matrix::colSums(mBC@raw.data))
mBC = AddMetaData(object = mBC, metadata = nReads_log, col.name = "nReads.log")


#plot statistics
file.name=paste(dir.name, sample.name, "stats.png", sep='')
png(file.name, width = 1024, height = 800)
par(mfrow = c(2, 2))
GenePlot(object = mBC, gene1 = "nGene", gene2 = "percent.mito", cex.use = 0.3)
GenePlot(object = mBC, gene1 = "nGene", gene2 = "percent.ribo", cex.use = 0.3)
GenePlot(object = mBC, gene1 = "nReads", gene2 = "nGene", cex.use = 0.3)
dev.off()

#filter out doublets
mBC = FilterCells(object = mBC, subset.names = c("percent.mito"),
                  low.thresholds = c(-Inf), high.thresholds = c(0.2))

doublets = mBC@meta.data$TagIDs %in% c("doublet", "not_detected")
mBC = SubsetData(object = mBC, cells.use=rownames(mBC@meta.data[!doublets,]),
                 subset.raw = TRUE)

raw.data = mBC@assays$RNA@counts
raw.data = mBC@raw.data
tmp = paste(mBC@meta.data$TagIDs, rownames(mBC@meta.data), sep="_")
colnames(raw.data)=tmp
colnames(raw.data)=gsub("not_detected", "not-detected", colnames(raw.data))
raw.data=as.matrix(raw.data)
fwrite(as.data.frame(raw.data), file = "matrix_inflection_demulti_DBEC_MiyakeWTA1.txt.gz", 
       row.names=T, col.names=T, sep="\t", eol="\n", quote=F, compress="gzip")

#===============================Seurat Object Updating, Normalization, Variable Feature Detection===============================#
#Following scripts were conducted under Seurat ver 4
mBC <- UpdateSeuratObject(mBC)
View(mBC@meta.data)
metadata<- mBC@meta.data
metadata$Mouse_origin <- "WT"
metadata$Mouse_origin[which(str_detect(metadata$TagIDs, paste(c("A0313","A0314"), collapse = "|")))] <- "CCR2"
metadata$sample <- "WT Day3"
metadata$sample[which(str_detect(metadata$TagIDs, "A0312"))] <- "WT Day5"
metadata$sample[which(str_detect(metadata$TagIDs, "A0313"))] <- "CCR2 Day3"
metadata$sample[which(str_detect(metadata$TagIDs, "A0314"))] <- "CCR2 Day5"
metadata$exp_day <- "Day3"
metadata$exp_day[which(str_detect(metadata$TagIDs, paste(c("A0312","A0314"), collapse = "|")))] <- "Day5"
View(metadata)
mouse_levels = c("WT","CCR2")
metadata$Mouse_origin= factor(x=metadata$Mouse_origin, levels = mouse_levels)
sample_levels = c("WT Day3","WT Day5","CCR2 Day3","CCR2 Day5")
metadata$sample= factor(x=metadata$sample, levels = sample_levels)
View(metadata)
metadata$mito_plus_ribo = metadata$percent.mito + metadata$percent.ribo
metadata$log10GenesPerUMI <- log10(metadata$nFeature_RNA) / log10(metadata$nCount_RNA)
mBC@meta.data <- metadata

# Seurat clustering 
mBC = NormalizeData(object = mBC, scale.factor=1000000)
all.genes <- rownames(mBC)
mBC = ScaleData(object = mBC, vars.to.regress = c("nReads"), features = all.genes, verbose = TRUE)
mBC = FindVariableFeatures(mBC, selection.method = "vst", nfeatures = 2000)
mBC<- RunPCA(mBC, npcs = 100, verbose = FALSE)

mBC = JackStraw(mBC, num.replicate = 100, dims = 100) 
mBC = ScoreJackStraw(mBC, dims = 1:100)
JackStrawPlot(mBC, dims = 1:100)
ElbowPlot(object = mBC)

#Clustering
tmp = as.data.frame(mBC@reductions$pca@jackstraw@overall.p.values)
tmp1 = tmp[tmp$Score>0.05, 1]
dims = c(1:(min(tmp1)-1))

#UMA
mBC <- RunUMAP(mBC, reduction = "pca", dims = dims)
mBC <- FindNeighbors(mBC, reduction = "pca", dims = dims)
mBC <- FindClusters(mBC, resolution = 1.0)

DimPlot(mBC, reduction = "umap", label = T)+coord_fixed()
DimPlot(mBC, reduction = "umap", split.by = "sample", ncol=2, label=T)+coord_fixed()
DimPlot(mBC, reduction = "umap", split.by = "Mouse_origin", ncol=2, label=T)+coord_fixed()

Idents(mBC) <- "Mouse_origin"
mBC_WT = subset(mBC, idents = "WT")
mBC_KO = subset(mBC, idents = "CCR2")
Idents(mBC) <- "seurat_clusters"
Idents(mBC_WT) <- "seurat_clusters"
Idents(mBC_KO) <- "seurat_clusters"

FeaturePlot(mBC_KO, features = c("Il1r1","Pdgfra","Cxcl5"), coord.fixed = T, ncol=2)

# find markers for every cluster compared to all remaining cells, report only the positive ones
mBC.markers <- FindAllMarkers(mBC, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
mBC.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(mBC_new, features = top10$gene) + NoLegend()
write.csv(mBC.markers, file= "filtered_UMAP_markergenes.csv")

# MoMac re-clustering 1st round
MoMac = subset(mBC, idents = c(3,9,12,13,14,20))
MoMac = NormalizeData(MoMac, scale.factor = 1000000)
MoMac = FindVariableFeatures(MoMac, mean.function = ExpMean, 
                             dispersion.function = LogVMR, 
                             mean.cutoff = c(0.1,Inf),
                             dispersion.cutoff = c(0.5,Inf))
all.genes <- rownames(MoMac)

MoMac = ScaleData(MoMac, vars.to.regress = c("nReads"), features = all.genes)
MoMac = RunPCA(MoMac, features = VariableFeatures(MoMac),
               npcs = 100)
MoMac = ProjectDim(MoMac)
MoMac = JackStraw(MoMac, num.replicate = 100, dims = 100)
MoMac = ScoreJackStraw(MoMac, dims = 1:100)
JackStrawPlot(MoMac, dims = 1:40)
ElbowPlot(object = MoMac)
tmp = as.data.frame(MoMac@reductions$pca@jackstraw@overall.p.values)
tmp1 = tmp[tmp$Score>0.05, 1]
dims = c(1:(min(tmp1)-1))
MoMac = FindNeighbors(MoMac, dims = dims)
MoMac = FindClusters(MoMac, resolution = 0.5)
MoMac = RunUMAP(MoMac, dims = dims)
DimPlot(MoMac, reduction = "umap", pt.size = 1, label=T)+coord_fixed()
DimPlot(WTMac, reduction = "umap", split.by="sample", pt.size = 1, ncol=2)+coord_fixed()
FeaturePlot(MoMac, features = c("Treml4","Ly6c1","Ccr2"),pt.size = 0.5, coord.fixed = T, ncol=1)

DimPlot(MoMac, reduction = "umap", split.by="Mouse_origin", ncol=1)+coord_fixed()

cols = c('0'="#00B0F0", "1"="#00B050", "2"="#FF0000","3"="#7030A0",  "4"="#EF05DE")
FeaturePlot(MoMac,features = c("Ly6c2","Adgre1","Gas6","Folr2","H2-Ab1","Dpt"),pt.size = 0.5, coord.fixed = T)


save(MoMac,dims,file="230204_MoMac_reclustering_withDCs_plus13mac.rda")
load("/media/owner/Data/Single cell RNA-seq 202011/to三宅先生(IgECAI_single_cell)/230204_MoMac_reclustering_withDCs.rda")

# MoMac re-clustering 2nd round
Mac_lin = subset(MoMac, idents = c(0,1,2,3,4))

Mac_lin = NormalizeData(Mac_lin, scale.factor = 1000000)
Mac_lin = FindVariableFeatures(Mac_lin, mean.function = ExpMean, 
                               dispersion.function = LogVMR, 
                               mean.cutoff = c(0.1,Inf),
                               dispersion.cutoff = c(0.5,Inf))
all.genes <- rownames(Mac_lin)

Mac_lin = ScaleData(Mac_lin, vars.to.regress = c("nReads"), features = all.genes)
Mac_lin = RunPCA(Mac_lin, features = VariableFeatures(Mac_lin),
                 npcs = 100)
Mac_lin = ProjectDim(Mac_lin)
Mac_lin = JackStraw(Mac_lin, num.replicate = 100, dims = 100)
Mac_lin = ScoreJackStraw(Mac_lin, dims = 1:100)
JackStrawPlot(Mac_lin, dims = 1:40)
ElbowPlot(object = Mac_lin)
tmp = as.data.frame(Mac_lin@reductions$pca@jackstraw@overall.p.values)
tmp1 = tmp[tmp$Score>0.05, 1]
dims = c(1:(min(tmp1)-1))
Mac_lin = FindNeighbors(Mac_lin, dims = dims)
Mac_lin = FindClusters(Mac_lin, resolution = 0.5)
Mac_lin = RunUMAP(Mac_lin, dims = dims)
#00B050
cols = c('0'="#FF0000", "1"="#00B0F0", '2'="#00B050", '3'="#FFC000",'4'="#002060", '5'="#7030A0")
DimPlot(Mac_lin, reduction = "umap", pt.size = 1, cols=cols, label=T)+coord_fixed()
my_levels = c(5,1,0,3,4,2)
Mac_lin@active.ident = factor(x=Mac_lin@active.ident, levels = my_levels)

Idents(Mac_lin) = "TagIDs"
WTMac_lin = subset(Mac_lin, idents=c("A0311", "A0312"))
CCR2Mac_lin = subset(Mac_lin, idents=c("A0313", "A0314"))
Idents(WTMac_lin) = "seurat_clusters"
Idents(CCR2Mac_lin) = "seurat_clusters"
Idents(Mac_lin) = "seurat_clusters"

VlnPlot(WTMac_lin,feature = c("H2-Ab1","H2-DMb1","Ciita"), pt.size = 0, ncol=3, cols = cols)
VlnPlot(WTMac_lin,feature = c("Ly6c2","Ccr2","Csf1r"), pt.size = 0, ncol=3, cols = cols)
VlnPlot(WTMac_lin,feature = c("Flt3","Cd209a","Tmem123"), pt.size = 0, ncol=3, cols = cols)
save(Mac_lin,dims,file="230205_Mac_lin_re-reclustering_withDCs-2.rda")

MacwoDC =subset(Mac_lin, idents =c(0,1,3,5))
Idents(MacwoDC) = "TagIDs"
WTMacwoDC = subset(MacwoDC, idents=c("A0311", "A0312"))
CCR2MacwoDC = subset(MacwoDC, idents=c("A0313", "A0314"))
Idents(WTMacwoDC) = "seurat_clusters"
Idents(CCR2MacwoDC) = "seurat_clusters"
Idents(MacwoDC) = "seurat_clusters"

DimPlot(MacwoDC, reduction = "umap", pt.size = 1, cols=cols, label=T)+coord_fixed()
DimPlot(MacwoDC, reduction = "umap", split.by="Mouse_origin", pt.size = 1, ncol=1,cols=cols, label=T)+coord_fixed()
DimPlot(WTMacwoDC, reduction = "umap", split.by="exp_day", ncol=1,cols=cols, label=T)+coord_fixed()
FeaturePlot(WTMacwoDC,feature = c("Ly6c2","Adgre1","Ccr2"), pt.size = 0.5, ncol=1, coord.fixed = T)

my_levels = c(5,1,0,3)
WTMacwoDC@active.ident = factor(x=WTMacwoDC@active.ident, levels = my_levels)

VlnPlot(WTMacwoDC, features = c("Arg1","Mrc1","Cd163", "Pdpn","Retnla","Il10"), pt.size = 0, cols=cols)
VlnPlot(WTMacwoDC, features = c("Ly6c1", "Retnla","Chil3","Mrc1","Vegfa", "Ccl2"), cols=cols, pt.size = 0)
VlnPlot(WTMacwoDC, features = c("Apoe", "C1qa","Folr2","Abca1","Abcg1", "Trem2"), cols=cols, pt.size = 0)
VlnPlot(WTMacwoDC, features = c("Msx3", "Tmem26","Gja1","Il1rap","Il10","Il4ra"), cols=cols, pt.size = 0)

#DotPlot
my_levels2 = c(3,0,1,5)
WTMacwoDC@active.ident = factor(x=WTMacwoDC@active.ident, levels = my_levels2)
features <- c( "Mertk", "Itgav","Itgb5","Trem2", "Pros1","Gas6","Mfge8","Stab1","C1qa")
DotPlot(WTMacwoDC, features = features, dot.scale=12) + RotatedAxis()

################################################################################
# Differentail expression analysis between classical mono and early CMDMs
################################################################################
DEGs <- FindMarkers(WTMacwoDC, ident.1 = "5", ident.2 = "1", 
                    logfc.threshold = -Inf, min.pct = -Inf, min.diff.pct = -Inf)
DEGs$q_val <- p.adjust(DEGs$p_val, method = "BH")
DEGs$gene <- rownames(DEGs)

head(DEGs, n = 15)

#Gene Set Enrichment Analysis by using clusterProfiler package
DEGlist = DEGs[,c("gene", "avg_log2FC")] %>% arrange(desc(avg_log2FC))
DEGgenes = DEGlist$gene
DEGgenes = bitr(DEGgenes, 
                fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = org.Mm.eg.db)
colnames(DEGgenes)[1]="gene" 
DEGlist= merge(DEGlist, DEGgenes, by="gene")
DEG = as.vector(DEGlist$avg_log2FC)
names(DEG) <- as.vector(DEGlist$ENTREZID)
DEG = sort(DEG, decreasing = T)
gseGO_diff <- gseGO(geneList     = DEG,
                    OrgDb        = org.Mm.eg.db,
                    ont          = "BP",
                    pvalueCutoff = 0.5,
                    pAdjustMethod = "BH",
                    verbose      = FALSE)
res_GO = as.data.frame(gseGO_diff)

# GSEA visualization
GsID = "GO:0006119"
GsName = res_GO$Description[res_GO$ID==GsID]
Gspval=res_GO$qvalue[res_GO$ID==GsID]
GsNES=res_GO$NES[res_GO$ID==GsID]
p<-gseaplot(gseGO_diff, geneSetID = GsID, color.line = "red", color.vline = "red", title = paste0(GsID,": ",GsName))
p[[2]] <-p[[2]]+
  annotate("text",x=2500,y=-0.4, label=paste0("NES: ",round(GsNES,3)),size=6,hjust=0, vjust=0)+
  annotate("text",x=6000,y=-0.3, label=paste0("qvalue: ",formatC(Gspval, digits = 2, format = "E")),size=6)
p[[1]]+p[[2]]+plot_layout(ncol = 1)

GsID = "GO:0034341"
GsName = res_GO$Description[res_GO$ID==GsID]
Gspval=res_GO$qvalue[res_GO$ID==GsID]
GsNES=res_GO$NES[res_GO$ID==GsID]
p<-gseaplot(gseGO_diff, geneSetID = GsID, color.line = "red", color.vline = "red", title = paste0(GsID,": ",GsName))
p[[2]] <-p[[2]]+
  annotate("text",x=15000,y=0.3, label=paste0("NES: ",round(GsNES,3)),size=6,hjust=0, vjust=0)+
  annotate("text",x=18000,y=0.4, label=paste0("qvalue: ",formatC(Gspval, digits = 2, format = "E")),size=6)
p[[1]]+p[[2]]+plot_layout(ncol = 1)

################################################################################
#Differentail expression analysis between early CMDMs and late CMDMs
################################################################################
DEGs <- FindMarkers(WTMacwoDC, ident.1 = "1", ident.2 = "0", 
                    logfc.threshold = -Inf, min.pct = -Inf, min.diff.pct = -Inf)
DEGs$q_val <- p.adjust(DEGs$p_val, method = "BH")
DEGs$gene <- rownames(DEGs)

head(DEGs, n = 15)
DEGs_latevsresident = DEGs

#Gene Set Enrichment Analysis by using clusterProfiler package
DEGlist = DEGs[,c("gene", "avg_log2FC")] %>% arrange(desc(avg_log2FC))
DEGgenes = DEGlist$gene
DEGgenes = bitr(DEGgenes, 
                fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = org.Mm.eg.db)
colnames(DEGgenes)[1]="gene" 
DEGlist= merge(DEGlist, DEGgenes, by="gene")
DEG = as.vector(DEGlist$avg_log2FC)
names(DEG) <- as.vector(DEGlist$ENTREZID)
DEG = sort(DEG, decreasing = T)
DEG_latevsresident = DEG
gseGO_diff <- gseGO(geneList     = DEG,
                    OrgDb        = org.Mm.eg.db,
                    ont          = "BP",
                    pvalueCutoff = 0.5,
                    pAdjustMethod = "BH",
                    verbose      = FALSE)

res_GO = as.data.frame(gseGO_diff)
gseGO_diff_simplified = simplify(gseGO_diff)


# GSEA visualization
dotplot(gseGO_diff_simplified, showCategory=40, x="NES", font.size = 10, label_format=100)

GsID = c("GO:0006909","GO:0016042","GO:0006910")
gg <- NULL
for (i in 1:length(GsID)){
GsName = res_GO$Description[res_GO$ID==GsID[i]]
Gspval=res_GO$qvalue[res_GO$ID==GsID[i]]
GsNES=res_GO$NES[res_GO$ID==GsID[i]]
p<-gseaplot(gseGO_diff, geneSetID = GsID[i], color.line = "red", color.vline = "red", title = paste0(GsID[i],": ",GsName))
p[[2]] <-p[[2]]+
  annotate("text",x=2500,y=-0.4, label=paste0("NES: ",round(GsNES,3)),size=6,hjust=0, vjust=0)+
  annotate("text",x=6000,y=-0.3, label=paste0("qvalue: ",formatC(Gspval, digits = 2, format = "E")),size=6)
gg[[i]] <- p[[1]]+p[[2]]+plot_layout(ncol = 1)
}
gg[[1]]
gg[[2]]
gg[[3]]

#KEGG GSEA
gseGO_KEGG <- gseKEGG(geneList     = DEG,
                      organism        = "mmu",
                      keyType = "ncbi-geneid",
                      pvalueCutoff = 0.05,
                      use_internal_data = F, 
                      verbose      = FALSE)
res_KEGG = as.data.frame(gseGO_KEGG)
View(res_KEGG)

GsID = c("mmu04142","mmu04979")
gg <- NULL
for (i in 1:length(GsID)){
GsName = res_KEGG$Description[res_KEGG$ID==GsID[i]]
Gspval=res_KEGG$qvalue[res_KEGG$ID==GsID[i]]
GsNES=res_KEGG$NES[res_KEGG$ID==GsID[i]]
p <-gseaplot(gseGO_KEGG, geneSetID = GsID[i], color.line = "red", color.vline = "red", title = paste0(GsID[i],": ",GsName))
p[[2]]<-p[[2]]+
  annotate("text",x=2500,y=-0.4, label=paste0("NES: ",round(GsNES,3)),size=6,hjust=0, vjust=0)+
  annotate("text",x=6000,y=-0.3, label=paste0("qvalue: ",formatC(Gspval, digits = 2, format = "E")),size=6)
gg[[i]] <- p[[1]]+p[[2]]+plot_layout(ncol = 1)
}
gg[[1]]
gg[[2]]

# Reactome GSEA
gseGO_Reactome <- gsePathway(geneList= DEG,
                             organism = "mouse",
                             minGSSize=120, 
                             pvalueCutoff=0.05,
                             pAdjustMethod="BH", 
                             verbose=FALSE)
res_Rea = as.data.frame(gseGO_Reactome)

GsID = "R-MMU-9013149"
GsName = res_Rea$Description[res_Rea$ID==GsID]
Gspval=res_Rea$qvalue[res_Rea$ID==GsID]
GsNES=res_Rea$NES[res_Rea$ID==GsID]
p<-gseaplot(gseGO_Reactome, geneSetID = GsID, color.line = "red", color.vline = "red", title = paste0(GsID,": ",GsName))
p[[2]] <-p[[2]]+
  annotate("text",x=2500,y=-0.4, label=paste0("NES: ",round(GsNES,3)),size=6,hjust=0, vjust=0)+
  annotate("text",x=6000,y=-0.3, label=paste0("qvalue: ",formatC(Gspval, digits = 2, format = "E")),size=6)
p[[1]]+p[[2]]+plot_layout(ncol = 1)


###################################################################################
#Differentail expression analysis between late CMDMs and resident-like Macs
###################################################################################
DEGs <- FindMarkers(WTMacwoDC, ident.1 = "0", ident.2 = "3", 
                    logfc.threshold = -Inf, min.pct = -Inf, min.diff.pct = -Inf)
DEGs$q_val <- p.adjust(DEGs$p_val, method = "BH")
DEGs$gene <- rownames(DEGs)

head(DEGs, n = 15)

#Gene Set Enrichment Analysis by using clusterProfiler package
DEGlist = DEGs[,c("gene", "avg_log2FC")] %>% arrange(desc(avg_log2FC))
DEGgenes = DEGlist$gene
DEGgenes = bitr(DEGgenes, 
                fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = org.Mm.eg.db)
colnames(DEGgenes)[1]="gene" 
DEGlist= merge(DEGlist, DEGgenes, by="gene")
DEG = as.vector(DEGlist$avg_log2FC)
names(DEG) <- as.vector(DEGlist$ENTREZID)
DEG = sort(DEG, decreasing = T)

gseGO_diff <- gseGO(geneList     = DEG,
                    OrgDb        = org.Mm.eg.db,
                    ont          = "BP",
                    pvalueCutoff = 0.5,
                    pAdjustMethod = "BH",
                    verbose      = FALSE)

res_GO = as.data.frame(gseGO_diff)

# GSEA visualization
GsID = "GO:0006909"
GsName = res_GO$Description[res_GO$ID==GsID]
Gspval=res_GO$qvalue[res_GO$ID==GsID]
GsNES=res_GO$NES[res_GO$ID==GsID]
p<-gseaplot(gseGO_diff, geneSetID = GsID, color.line = "red", color.vline = "red", title = paste0(GsID,": ",GsName))
p[[2]] <-p[[2]]+
  annotate("text",x=15000,y=0.2, label=paste0("NES: ",round(GsNES,3)),size=6,hjust=0, vjust=0)+
  annotate("text",x=18000,y=0.3, label=paste0("qvalue: ",formatC(Gspval, digits = 2, format = "E")),size=6)
p[[1]]+p[[2]]+plot_layout(ncol = 1)

#KEGG GSEA
gseGO_KEGG <- gseKEGG(geneList     = DEG,
                      organism        = "mmu",
                      keyType = "ncbi-geneid",
                      pvalueCutoff = 0.05,
                      use_internal_data = F, 
                      verbose      = FALSE)
res_KEGG = as.data.frame(gseGO_KEGG)

GsID = "mmu04151"
GsName = res_KEGG$Description[res_KEGG$ID==GsID]
Gspval=res_KEGG$qvalue[res_KEGG$ID==GsID]
GsNES=res_KEGG$NES[res_KEGG$ID==GsID]
p <-gseaplot(gseGO_KEGG, geneSetID = GsID, color.line = "red", color.vline = "red", title = paste0(GsID,": ",GsName))
p[[2]]<-p[[2]]+
  annotate("text",x=15000,y=0.2, label=paste0("NES: ",round(GsNES,3)),size=6,hjust=0, vjust=0)+
  annotate("text",x=18000,y=0.3, label=paste0("qvalue: ",formatC(Gspval, digits = 2, format = "E")),size=6)
p[[1]]+p[[2]]+plot_layout(ncol = 1)

##################################################
#Pseudotime analysis by slingshot
##################################################
sce <- as.SingleCellExperiment(WTMacwoDC)
reducedDim(sce) <- reducedDim(sce)[, dims]
sce <- suppressWarnings(slingshot(
  sce,
  reducedDim = "UMAP",
  clusterLabels = WTMacwoDC@active.ident,
  start.clus = "5"
))
sds <- SlingshotDataSet(sce)

require(scales)

cell_pal <- function(cell_vars, pal_fun,...) {
  if (is.numeric(cell_vars)) {
    pal <- pal_fun(100, ...)
    return(pal[cut(cell_vars, breaks = 100)])
  } else {
    categories <- sort(unique(cell_vars))
    pal <- setNames(pal_fun(length(categories), ...), categories)
    return(pal[cell_vars])
  }
}

#
par(xpd=TRUE)
par(mar=c(3,3,3,3))
cell_colors_clust <- cell_pal(WTMacwoDC@active.ident, hue_pal())
plot(reducedDim(sds), col = cell_colors_clust, pch=16, cex = 0.5,asp=1)
lines(sds, lwd = 2, type = 'lineages', col = 'black')


#pseudotime color
nc <- 1  
pt <- slingPseudotime(sds)
nms <- colnames(pt)
nr <- ceiling(length(nms)/nc)


col_fun = colorRamp2(c(0, 7, 14), c("#120789", "#CD4A76", "#F0F921"))
par(mfrow = c(nr, nc))
for (i in nms) {
  colors <- col_fun(pt[,i])
  plot(reducedDim(sds), col = colors, pch = 16, cex = 0.8, main = i,  asp=1, bty="n")
  lines(sds, lwd = 2, col = 'black', type = 'lineages')
}

lgd = Legend(col_fun = col_fun, title = "pseudotime")
draw(lgd, x = unit(1, "npc"), y = unit(1, "npc"), just = c("right", "top"))


# fit negative binomial GAM
BPPARAM <- BiocParallel::bpparam()
sce <- fitGAM(sce, BPPARAM = BPPARAM)

Cluster = as.vector(WTMacwoDC@meta.data$seurat_clusters)
head(Cluster)
cols = c('0'="#FF0000", "1"="#00B0F0", '2'="#00B050", '3'="#FFC000",'4'="#002060", '5'="#7030A0")
counts <- counts(sce)
count_mtx = as.data.frame(counts)

save(sce,Cluster,file="230206_slingshot_pseudotime_recluster.rda")

Gene=c("Ly6c2","Adgre1","Pdcd1lg2")
plotSmoothers(sce, counts, gene = Gene[1], lwd=1,curvesCols = "black", border =F, pointCol=Cluster)+
  ggplot2::scale_x_continuous(limits = c(0, 15))+
  ggplot2::scale_color_manual(values = cols)+
  ggplot2::ggtitle(paste(Gene[1]))+
  theme(plot.title = element_text(size = 16, face = "italic"))+
  theme(plot.title = element_text(hjust = 0.5))+
plotSmoothers(sce, counts, gene = Gene[2], lwd=1,curvesCols = "black", border =F, pointCol=Cluster)+
  ggplot2::scale_x_continuous(limits = c(0, 15))+
  ggplot2::scale_color_manual(values = cols)+
  ggplot2::ggtitle(paste(Gene[2]))+
  theme(plot.title = element_text(size = 16, face = "italic"))+
  theme(plot.title = element_text(hjust = 0.5))+ 
plotSmoothers(sce, counts, gene = Gene[3], lwd=1,curvesCols = "black", border =F, pointCol=Cluster)+
  ggplot2::scale_x_continuous(limits = c(0, 15))+
  ggplot2::scale_color_manual(values = cols)+
  ggplot2::ggtitle(paste(Gene[3]))+
  theme(plot.title = element_text(size = 16, face = "italic"))+
  theme(plot.title = element_text(hjust = 0.5))+ 
  plot_layout(ncol = 1)