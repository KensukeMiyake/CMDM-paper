library(TCC)
library(dplyr)
library(ggplot2)
library(ggbiplot)
library(Seurat)

setwd("/media/owner/Data/bulk RNA-seq moMac sorting")
data = read.table("/media/owner/Data/bulk RNA-seq moMac sorting/Miyake18_RNASeq_expression_table.txt", header=T)
head(data)
rownames(data) = data[,c(1)]
data = data[, c(2:13)]
head(data)
data = data[, c(1,5,6,7,8,9,10,11,12,2,3,4)]
colnames(data)
data = data[, c(1,5,9,2,6,10,3,7,11,4,8,12)]
data = data[, c(10,11,12,7,8,9,4,5,6,1,2,3)]
colnames(data) =c("cMono_1","cMono_2","cMono_3","EarlyCMDM_1","EarlyCMDM_2","EarlyCMDM_3",
                  "LateCMDM_1","LateCMDM_2","LateCMDM_3","Resident_1","Resident_2","Resident_3")
head(data)
dim(data)

#construct TCC cluss object
group = c(rep(1,"cMono", 3), rep(2,"EarlyCMDM", 3), rep(3,"LateCMDM", 3),rep(4,"Resident", 3))
tcc = new("TCC", data, group)
head(tcc)

#filter low-count genes
dim(tcc$count)
tcc <- filterLowCountGenes(tcc)
dim(tcc$count)

#normalization of multi-group count data
#DEGES/DESeq2
tcc <- calcNormFactors(tcc, norm.method="tmm", test.method="edger",
                       iteration=3, FDR=0.1, floorPDEG=0.05)
tcc$norm.factors
normalized.count <- getNormalizedData(tcc)
write.csv(normalized.count, file = "230804_Bulk_MoMac_normalized_count.csv")
head (normalized.count)
dim(normalized.count)


#DE analysis among 3 groups
designall <- model.matrix(~as.factor(group))
designall
tcc <- estimateDE(tcc, test.method="edger", FDR=0.1,
                  design=designall,coef=c(2:4))
result <- getResult(tcc, sort=FALSE)
head(result)
sum(tcc$stat$q.value < 0.1) 
tmp <- cbind( normalized.count, result)
head(tmp)

#Fold Change
mean_G1 <- log2(apply(sweep(as.matrix(normalized.count[,group==1]), 1, 0.125, FUN="+"),1,mean))
mean_G2 <- log2(apply(sweep(as.matrix(normalized.count[,group==2]), 1, 0.125, FUN="+"),1,mean))
mean_G3 <- log2(apply(sweep(as.matrix(normalized.count[,group==3]), 1, 0.125, FUN="+"),1,mean))
mean_G4 <- log2(apply(sweep(as.matrix(normalized.count[,group==4]), 1, 0.125, FUN="+"),1,mean))

FoldChange <-cbind (mean_G1,mean_G2,mean_G3,mean_G4)
head(FoldChange)
MAX_FC<- apply(FoldChange,1,max)
MIN_FC<- apply(FoldChange,1,min)
FoldChange = cbind( FoldChange,MAX_FC,MIN_FC)
tmp<- cbind(tmp,FoldChange)
colnames(tmp)
tmp =tmp[,-c(14:15)]
head(tmp)
write.csv(tmp, file="230804_Monosorting_DEG_TCC_all.csv", quote=F, row.names=T)

#PCA
norm.count = as.matrix(normalized.count)
pca = prcomp(t(norm.count), scale = TRUE)
group = c(rep("cMono", 3), rep("EarlyCMDM", 3), rep("LateCMDM", 3),rep("Resident", 3))
ggbiplot(pca, groups=group, var.axes  =F,ellipse = T, ellipse.prob = 0.68)+
  scale_color_manual(values=c("#7030A0", "#00B0F0", "#FF0000", "#FFC000"))+
  theme_bw()+
  ggtitle("Principal Component Analysis")


###########################################################
# DE analysis between 2 specific groups
###########################################################
# Group2 vs Group3
design <- model.matrix(~0+as.factor(group))
design
design2<-matrix(c(design[,c(1)],rep(1,12),design[,c(3,4)]), nrow=12, ncol=4)
design2
tcc1 <- estimateDE(tcc, test.method="edger", FDR=0.1,
                   design=design2,coef=3)

result_CMDM <- getResult(tcc1, sort=FALSE)
sum(tcc1$stat$q.value < 0.1) 
m.value_calc<-mean_G3 - mean_G2
result_CMDM<-cbind(result_CMDM,m.value_calc)
head(result_CMDM)
result_CMDM = result_CMDM[,-c(2:3)]
tmp_CMDM <- cbind(normalized.count, result_CMDM)
head(tmp_CMDM)

# Group1 vs Group2
design <- model.matrix(~0+as.factor(group))
design
design1<-matrix(c(rep(1,12),design[,c(2:4)]), nrow=12, ncol=4)
design1
tcc2 <- estimateDE(tcc, test.method="edger", FDR=0.1,
                   design=design1,coef=2)
result_cMono <- getResult(tcc2, sort=FALSE)
sum(tcc2$stat$q.value < 0.1) 
m.value_calc<-mean_G2 - mean_G1
result_cMono<-cbind(result_cMono,m.value_calc)
head(result_cMono)
result_cMono = result_cMono[,-c(2:3)]
tmp_cMono <- cbind(normalized.count, result_cMono)
head(tmp_cMono)

# Group3 vs Group4
design <- model.matrix(~0+as.factor(group))
design
design1<-matrix(c(rep(1,12),design[,c(2:4)]), nrow=12, ncol=4)
design1
tcc5 <- estimateDE(tcc, test.method="edger", FDR=0.1,
                   design=design1,coef=4)
result_resMono <- getResult(tcc5, sort=FALSE)
sum(tcc5$stat$q.value < 0.1) 
m.value_calc<-mean_G4 - mean_G1
result_resMono<-cbind(result_resMono,m.value_calc)
head(result_resMono)
result_resMono = result_resMono[,-c(2:3)]
tmp_resMono <- cbind(normalized.count, result_resMono)
head(tmp_resMono)
write.csv(tmp, file="230404_BMBA_Difamilast.csv", quote=F, row.names=T)

# Group1 vs Group4
design <- model.matrix(~0+as.factor(group))
design
design3<-matrix(c(design[,c(1:2)],rep(1,12),design[,c(4)]), nrow=12, ncol=4)
design3
tcc3 <- estimateDE(tcc, test.method="edger", FDR=0.1,
                   design=design3,coef=4)

result_resident <- getResult(tcc3, sort=FALSE)
sum(tcc3$stat$q.value < 0.1) 
m.value_calc<-mean_G4 - mean_G3
result_resident<-cbind(result_resident,m.value_calc)
head(result_resident)
result_resident = result_resident[,-c(2:3)]
tmp_resident <- cbind(normalized.count, result_resident)
head(tmp_resident)
write.csv(tmp, file="230404_BMBA_Difamilast.csv", quote=F, row.names=T)


#enrichment analysis
DEGs_earlyCMDM1 = result_CMDM%>% filter(m.value_calc< 0) %>% filter(estimatedDEG==1)
DEGs_earlyCMDM2 = result_cMono %>% filter(m.value_calc > 0) %>% filter(estimatedDEG==1) 
DEGs_lateCMDM1 = result_CMDM%>% filter(m.value_calc> 0) %>% filter(estimatedDEG==1)
DEGs_lateCMDM2 = result_resident %>% filter(m.value_calc < 0) %>% filter(estimatedDEG==1) 
DEGs_cMono1 = result_cMono %>% filter(m.value_calc < 0) %>% filter(estimatedDEG==1) 
DEGs_cMono2 = result_resMono %>% filter(m.value_calc < 0) %>% filter(estimatedDEG==1) 
DEGs_res1 = result_resident %>% filter(m.value_calc > 0) %>% filter(estimatedDEG==1) 
DEGs_res2 = result_resMono %>% filter(m.value_calc > 0) %>% filter(estimatedDEG==1) 

earlyCMDMgenes <- intersect(rownames(DEGs_earlyCMDM1),rownames(DEGs_earlyCMDM2))
lateCMDMgenes <- intersect(rownames(DEGs_lateCMDM1),rownames(DEGs_lateCMDM2))
cMonogenes <- intersect(rownames(DEGs_cMono1),rownames(DEGs_cMono2))
Residentgenes <- intersect(rownames(DEGs_res1),rownames(DEGs_res2))

#load 2nd reclustered Mo-Mac Seurat object (named as Mac_lin)
load("Mac_lin_re-reclustering_withDCs.rda")
MacwoDC =subset(Mac_lin, idents =c(0,1,3,5))
Idents(MacwoDC) = "TagIDs"
WTMacwoDC = subset(MacwoDC, idents=c("A0311", "A0312"))
CCR2MacwoDC = subset(MacwoDC, idents=c("A0313", "A0314"))
Idents(WTMacwoDC) = "seurat_clusters"
Idents(CCR2MacwoDC) = "seurat_clusters"
Idents(MacwoDC) = "seurat_clusters"

earlyCMDM_list = list(intersect(earlyCMDMgenes, as.vector(rownames(WTMacwoDC))))
lateCMDM_list = list(intersect(lateCMDMgenes, as.vector(rownames(WTMacwoDC))))
cMono_list = list(intersect(cMonogenes, as.vector(rownames(WTMacwoDC))))
Resident_list = list(intersect(Residentgenes, as.vector(rownames(WTMacwoDC))))

#Add module scores to Seurat Object
WTMacwoDC = AddModuleScore(WTMacwoDC, features = earlyCMDM_list, name = "earlyCMDM")
WTMacwoDC = AddModuleScore(WTMacwoDC, features = lateCMDM_list, name = "lateCMDM")
WTMacwoDC = AddModuleScore(WTMacwoDC, features = cMono_list, name = "cMono")
WTMacwoDC = AddModuleScore(WTMacwoDC, features = Resident_list, name = "Resident")

FeaturePlot(WTMacwoDC, features = c("cMono1","earlyCMDM1","lateCMDM1","Resident1"), pt.size = 1, ncol=2, min.cutoff = 'q2',max.cutoff='q98',coord.fixed = T)
cols = c('0'="#FF0000", "1"="#00B0F0", '3'="#FFC000",'5'="#7030A0")
my_levels = c(5,1,0,3)
WTMacwoDC@active.ident = factor(x=WTMacwoDC@active.ident, levels = my_levels)
VlnPlot(WTMacwoDC, features = c("cMono1","earlyCMDM1","lateCMDM1","Resident1"), cols = cols, pt.size = 0, ncol=4)+NoLegend()
