setwd("/Users/houhouhou/opt/R/GSE128033")
folders=list.files('./',pattern="^IPF")
folders
library(Seurat)
scList = lapply(folders,function(folder){ CreateSeuratObject(counts = Read10X(folder),project = folder,min.cells = 3, min.features = 200)})
scList = lapply(folders,function(folder){ CreateSeuratObject(counts = Read10X_h5(folder),project = folder,min.cells = 3, min.features = 200)})
scList = lapply(folders,function(folder){readRDS(folder)})

HC<- merge(scList[[1]], 
           y = c(scList[[2]], scList[[3]],scList[[4]],scList[[5]]), 
           add.cell.ids = c("HC1","HC2", "HC3","HC4","HC5"), 
           project = "HC")
metadata=HC@meta.data
HC@meta.data$group="HC"

IPF<- merge(scList[[6]], 
            y = c(scList[[7]], scList[[8]], scList[[9]],scList[[10]]), 
            add.cell.ids = c("IPF1","IPF2", "IPF3", "IPF4","IPF5"), 
            project = "IPF")
metadata=IPF@meta.data
IPF@meta.data$group="IPF"

#计算线粒体基因比例
HC[["percent.mt"]] <- PercentageFeatureSet(HC,pattern = "^MT")
IPF[["percent.mt"]] <- PercentageFeatureSet(IPF,pattern = "^MT")
VlnPlot(HC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, 
        group.by = "orig.ident", 
        pt.size = 0)
VlnPlot(IPF, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, 
        group.by = "orig.ident", 
        pt.size = 0)

#进行质控
HC<- subset(HC, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 15)
IPF<- subset(IPF, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 15)
VlnPlot(HC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, 
        group.by = "orig.ident", 
        pt.size = 0)
VlnPlot(IPF, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, 
        group.by = "orig.ident", 
        pt.size = 0)

#使用FindIntegrationAnchors合并数据，IntegrateData去除批次效应
HC<- NormalizeData(HC)
HC<- FindVariableFeatures(HC, nfeatures = 2000)
IPF<- NormalizeData(IPF)
IPF<- FindVariableFeatures(IPF, nfeatures = 2000)
sampleList <- list(HC, IPF)
scedata <- FindIntegrationAnchors(object.list = sampleList, dims = 1:30)
IPF_combined<- IntegrateData(anchorset = scedata, dims = 1:30)
save(IPF_combined, file = "IPF_combined.RData")

#细胞周期
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
IPF_combined <- CellCycleScoring(IPF_combined, 
                                 s.features = s.genes, 
                                 g2m.features = g2m.genes,
                                 set.ident = TRUE)
VlnPlot(IPF_combined,features = c("S.Score","G2M.Score"),group.by = "orig.ident")

#降维聚类
IPF_combined <- ScaleData(IPF_combined, vars.to.regress = c("S.Score", "G2M.Score"), verbose = FALSE)
IPF_combined <- RunPCA(IPF_combined, npcs = 50, verbose = FALSE)
IPF_combined <- RunHarmony(IPF_combined, reduction="pca",group.by.vars="group",reduction.save="harmony")
DimPlot(IPF_combined,reduction="pca")
DimPlot(IPF_combined,reduction="harmony")
ElbowPlot(IPF_combined, ndims=50, reduction="pca")
Loadings(object = IPF_combined[["pca"]])   ###获取基因在pc轴上的投射值
Embeddings(object = IPF_combined[["pca"]])  ###获取各个细胞的pc值
Stdev(IPF_combined)   ###获取各pc轴解释量方差
DimHeatmap(IPF_combined, dims = 1:20, nfeatures=10, cells = 500, balanced = TRUE)###查看决定pc值的top10基因在500个细胞中的热图
IPF_combined <- FindNeighbors(IPF_combined, reduction = "pca", dims = 1:30)
IPF_combined <- FindClusters(IPF_combined, 
                             resolution=seq(from = 0.1, 
                                            to = 1.0, 
                                            by = 0.1))
IPF_combined <- FindClusters(IPF_combined, 
                             resolution=0.4)
IPF_combined <- RunUMAP(IPF_combined, reduction = "pca", dims = 1:30)
IPF_combined <- RunUMAP(IPF_combined, reduction = "harmony", dims = 1:30)
IPF_combined <- RunTSNE(IPF_combined, reduction = "pca", dims = 1:30)
IPF_combined <- RunTSNE(IPF_combined, reduction = "harmony", dims = 1:30)
library(clustree)
clustree(IPF_combined)
Idents(IPF_combined) <- "integrated_snn_res.0.3" 
DimPlot(IPF_combined)
DimPlot(IPF_combined,reduction="tsne")
IPF_combined$seurat_clusters <- IPF_combined@active.ident
DimPlot(IPF_combined,label = T,reduction="umap",split.by = "orig.ident",ncol = 3)
DimPlot(IPF_combined,label = T,reduction="umap",split.by = "group")
DimPlot(IPF_combined,label = T,reduction="tsne", split.by = "orig.ident",ncol = 3)
DimPlot(IPF_combined,label = T,reduction="tsne",split.by = "group")

#鉴定各个群的高表达基因
DefaultAssay(IPF_combined) <- "RNA"
all.markers  <- FindAllMarkers(IPF_combined, 
                               only.pos = TRUE, 
                               min.pct = 0.1, 
                               logfc.threshold = 0.25)
significant.markers  <- all.markers [all.markers $p_val_adj < 0.2, ]
write.csv(significant.markers, file = "significant.markers.csv")
library(tidyverse)
top_markers <- significant.markers%>%group_by(cluster)%>%top_n(10, wt = avg_log2FC)
DoHeatmap(IPF_combined, features=top_markers$gene, group.by ="seurat_clusters",slot= "data") 

#呼吸系统常规markers
markers <- c("PTPRC","EPCAM","PECAM1","LYVE1","VWF", "COL1A2")
DotPlot(IPF_combined,features = markers)+coord_flip()

FeaturePlot(IPF_combined,features = c("PTPRC","EPCAM","PECAM1","LYVE1","VWF","COL1A2"))

alldata <- ScaleData(IPF_combined, 
                     features = markers, 
                     assay = "RNA")
DoHeatmap(alldata, 
          features = markers,
          group.by = "seurat_clusters",
          assay = "RNA")

#上皮细胞
FeaturePlot(IPF_combined,features = c("AGER","SFTPC","SCGB1A1","KRT7","TPPP3"))

#添加细胞群注释
new.cluster.ids <- c("0"="SPP1+ Macrophage", 
                     "1"="FABP4+ Macrophage", 
                     "2"="T cell", 
                     "3"="Fibroblast", 
                     "4"="AT1/AT2", 
                     "5"="Endothelial cell", 
                     "6"="Monocyte/Macrophage", 
                     "7"="T cell", 
                     "8"="Monocyte", 
                     "9"="NK cell", 
                     "10"="B/T cell/Monocyte", 
                     "11"="Endothelial cell", 
                     "12"="HMGB1+ Macrophage", 
                     "13"="Epithelial cell", 
                     "14"="Ciliated cell", 
                     "15"="Fibroblast", 
                     "16"="B cell",
                     "17"="Endothelial cell", 
                     "18"="Monocyte", 
                     "19"="Endothelial cell",
                     "20"="AT1/AT2")
IPF_combined <- RenameIdents(IPF_combined, new.cluster.ids)                        
IPF_combined$celltype <- IPF_combined@active.ident
DimPlot(IPF_combined, group.by = "celltype", label = T, reduction = "tsne")
DimPlot(IPF_combined, reduction="tsne",group.by = "celltype")
save(IPF_combined, file = "IPF_combined.RData")
DimPlot(IPF_combined, split.by = "group",ncol = 3)
DimPlot(IPF_combined, reduction="tsne", split.by = "group",ncol = 3)


#新细胞群的markers
newall.markers  <- FindAllMarkers(IPF_combined, 
                                  only.pos = TRUE, 
                                  min.pct = 0.1, 
                                  logfc.threshold = 0.25)
newmarkers  <- newall.markers [newall.markers $p_val_adj < 0.2, ]

#每个样本细胞数
table(IPF_combined$orig.ident)
prop.table(table(Idents(IPF_combined)))
table(Idents(IPF_combined), IPF_combined$orig.ident)
Cellratio <- prop.table(table(Idents(IPF_combined), IPF_combined$orig.ident), margin = 2)#计算各组样本不同细胞群比例
Cellratio
Cellratio <- as.data.frame(Cellratio)
colourCount = length(unique(Cellratio$Var1))
library(ggplot2)
ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  coord_flip()+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))

#每组细胞数
table(IPF_combined$group)
prop.table(table(Idents(IPF_combined)))
table(Idents(IPF_combined), IPF_combined$group)
Cellratio <- prop.table(table(Idents(IPF_combined), IPF_combined$group), margin = 2)
Cellratio
Cellratio <- as.data.frame(Cellratio)
colourCount = length(unique(Cellratio$Var1))
library(ggplot2)
ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  coord_flip()+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))

#几个感兴趣的细胞群
VlnPlot(IPF_combined, features = c("Pecam1", "Cd53","Hopx", "Sftpc", "Cd74","Lyz2"),
        ncol = 3, 
        group.by = "orig.ident",
        pt.size = 0)
FeaturePlot(IPF_combined,features = c("Pecam1", "Cd53","Hopx", "Sftpc", "Cd74","Lyz2"),
            ncol = 3, 
            split.by = "orig.ident")

#整体的分组差异分析
DefaultAssay(IPF_combined) <- "RNA"
Idents(IPF_combined) = IPF_combined$group 
table(Idents(IPF_combined))
deg_all = FindMarkers(IPF_combined,ident.1 = 'SSc',
                      ident.2 = 'HC')
head(deg_all[order(deg_all$p_val),])
table(Idents(IPF_combined))

library(EnhancedVolcano)
head(deg_all)
EnhancedVolcano(deg_all,
                lab = rownames(deg_all),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pCutoff = 0.05,
                FCcutoff = 0.5,
                pointSize = 3.0,
                labSize = 6.0,
                title = "SSc vs HC")

#内皮细胞亚群
Endothelial=subset(IPF_combined,idents=c("5","11","17","19"))
Endothelial<- NormalizeData(Endothelial)
Endothelial<- FindVariableFeatures(Endothelial, nfeatures = 2000)
Endothelial <- ScaleData(Endothelial, vars.to.regress = c("S.Score", "G2M.Score"), verbose = FALSE)
Endothelial <- RunPCA(Endothelial, npcs = 50, verbose = FALSE)
head(Endothelial)
Endothelial= SCTransform(Endothelial)
Endothelial <- RunHarmony(Endothelial, group.by.vars="orig.ident", assay.use="SCT", reduction.save= "harmony" )
# DimPlot(Endothelial,reduction="pca")
ElbowPlot(Endothelial, ndims=50, reduction="harmony")
# ElbowPlot(Endothelial, ndims=50, reduction="harmony")
Loadings(object = Endothelial[["pca"]])   ###获取基因在pc轴上的投射值
Embeddings(object = Endothelial[["pca"]])  ###获取各个细胞的pc值
Stdev(IPF_combined)   ###获取各pc轴解释量方差
DimHeatmap(IPF_combined, dims = 1:20, nfeatures=10, cells = 500, balanced = TRUE)###查看决定pc值的top10基因在500个细胞中的热图
Endothelial <- FindNeighbors(Endothelial, reduction = "harmony", dims = 1:20)
Endothelial <- FindClusters(Endothelial,resolution=seq(from = 0.1, 
                                                       to = 1.0, 
                                                       by = 0.1))
Endothelial <- RunUMAP(Endothelial, reduction = "harmony", dims = 1:20)
Endothelial<- RunTSNE(Endothelial, reduction = "harmony", dims = 1:20)
library(clustree)
clustree(Endothelial)
Idents(Endothelial) <- "integrated_snn_res.0.2" 
DimPlot(Endothelial,label =T)
DimPlot(Endothelial,reduction="tsne",label=T)
Endothelial$seurat_clusters <- Endothelial@active.ident
DimPlot(Endothelial,label = T,split.by = "group",ncol = 3) 
DimPlot(Endothelial,label = T,reduction="tsne", split.by = "group",ncol = 3) 

#Endothelial分群的marker分析
DefaultAssay(Endothelial) <- "RNA"
endothelial_markers  <- FindAllMarkers(Endothelial, 
                                       only.pos = TRUE, 
                                       min.pct = 0.1, 
                                       logfc.threshold = 0.25)
write.csv(endothelial_markers, file = "endothelial_markers.csv")
sig_endothelialmarkers  <- endothelial_markers [endothelial_markers $p_val_adj < 0.2, ]
write.csv(sig_endothelialmarkers, file = "sig_endothelialmarkers.csv")

top_endothelial_markers <- endothelial_markers%>%group_by(cluster)%>%top_n(10, wt = avg_log2FC)
DoHeatmap(Endothelial, features=top_endothelial_markers$gene, group.by ="seurat_clusters") 
DoHeatmap(Endothelial,features = c("VWF","SLC6A2","CD200","VCAM1","BST1",
                                   "MPG","PLAT","GJA4","GJA5","EDN1","FBLN2","HEY1",
                                   "GPIHBP1","SEMA3C","CADM1","HILPDA",
                                   "CA4","FIBIN","CYP4B1","EDNRB","KDR",
                                   "MMM1","CCL21","FXYD6","FGL2","CP","CCL21A",
                                   "MKI67","TOP2A"),group.by = "seurat_clusters")
#endothelial的注释
Endothelial=subset(Endothelial, idents =c("0","1","2","4","5","6"))
new.cluster.ids <- c("0"="Vein",
                     "1"="Vein", 
                     "2"="aCap", 
                     "4"="Lymph",
                     "5"="Artery",
                     "6"="gCap")
Endothelial <- RenameIdents(Endothelial, new.cluster.ids)                        
Endothelial$celltype <- Endothelial@active.ident
DimPlot(Endothelial, group.by = "celltype")
DimPlot(Endothelial, reduction="tsne",group.by = "celltype")
save(Endothelial, file = "Endothelial.RData")
DimPlot(Endothelial, split.by = "group",ncol = 3)
DimPlot(Endothelial, reduction="tsne", split.by = "group",ncol = 3)
DimPlot(Endothelial, reduction="tsne", split.by = "orig.ident",ncol = 3)

#内皮细胞分群比例
table(Endothelial$group)
prop.table(table(Idents(Endothelial)))
table(Idents(Endothelial), Endothelial$group)
Cellratio <- prop.table(table(Idents(Endothelial), Endothelial$group), margin = 2)#计算各组样本不同细胞群比例
Cellratio
Cellratio <- as.data.frame(Cellratio)
colourCount = length(unique(Cellratio$Var1))
library(ggplot2)
ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  coord_flip()+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))



#分组分析出图
#vein
vein=subset(Endothelial,idents = "Vein")
Idents(vein) = vein$group 
table(Idents(vein))
deg_vein = FindMarkers(vein,ident.1 = 'IPF',
                       ident.2 = 'HC',min.pct = 0,
                       logfc.threshold = 0)
head(deg_vein[order(deg_vein$p_val),])
res1=deg_vein
sig_res1=subset(res1,p_val_adj<0.05&abs(avg_log2FC)>0.15)
library(EnhancedVolcano)
write.csv(sig_res1,file="res1.csv")
head(res1)
EnhancedVolcano(res1,
                lab = rownames(res1),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pCutoff = 0.05,
                FCcutoff = 0.5,
                pointSize = 3.0,
                labSize = 6.0,
                title = 'Vein:IPF vs HC')

#Artery
artery=subset(Endothelial,idents = "Artery")
Idents(artery) = artery$group 
table(Idents(artery))
deg_artery = FindMarkers(artery,ident.1 = 'IPF',
                         ident.2 = 'HC',min.pct = 0,
                         logfc.threshold = 0)
head(deg_artery[order(deg_artery$p_val),])
library(EnhancedVolcano)
res2=deg_artery
sig_res2=subset(res2,p_val_adj<0.05&abs(avg_log2FC)>0.15)
write.csv(sig_res2,file="res2.csv")
head(res2)
EnhancedVolcano(res2,
                lab = rownames(res2),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pCutoff = 0.05,
                FCcutoff = 0.5,
                pointSize = 3.0,
                labSize = 6.0,
                title = 'Artery:IPF vs HC')

#aCap
aCap=subset(Endothelial,idents = "aCap")
Idents(aCap) = aCap$group 
table(Idents(aCap))
deg_aCap = FindMarkers(aCap,ident.1 = 'IPF',
                       ident.2 = 'HC',min.pct = 0,
                       logfc.threshold = 0)
head(deg_aCap[order(deg_aCap$p_val),])
res3=deg_aCap
sig_res3=subset(res3,p_val_adj<0.05&abs(avg_log2FC)>0.15)
write.csv(sig_res3,file="res3.csv")
head(res3)
EnhancedVolcano(res3,
                lab = rownames(res3),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pCutoff = 0.05,
                FCcutoff = 0.5,
                pointSize = 3.0,
                labSize = 6.0,
                title = 'aCap:IPF vs HC')

#gCap
gCap=subset(Endothelial,idents = "gCap")
Idents(gCap) = gCap$group 
table(Idents(gCap))
deg_gCap = FindMarkers(gCap,ident.1 = 'IPF',
                       ident.2 = 'HC',min.pct = 0,
                       logfc.threshold = 0)
head(deg_gCap[order(deg_gCap$p_val),])
library(EnhancedVolcano)
res4=deg_gCap
sig_res4=subset(res4,p_val_adj<0.05&abs(avg_log2FC)>0.15)
write.csv(sig_res4,file="res4.csv")
head(res4)
EnhancedVolcano(res4,
                lab = rownames(res4),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pCutoff = 0.05,
                FCcutoff = 0.5,
                pointSize = 3.0,
                labSize = 6.0,
                title = 'gCap:IPF vs HC')

#Lymph
lym=subset(Endothelial,idents = "Lymph")
Idents(lym) = lym$group 
table(Idents(lym))
deg_lym = FindMarkers(lym,ident.1 = 'IPF',
                      ident.2 = 'HC',min.pct = 0,
                      logfc.threshold = 0)
head(deg_lym[order(deg_lym$p_val),])
library(EnhancedVolcano)
res5=deg_lym
sig_res5=subset(res5,p_val_adj<0.05&abs(avg_log2FC)>0.15)
write.csv(sig_res5,file="res5.csv")
head(res5)
EnhancedVolcano(res5,
                lab = rownames(res5),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pCutoff = 0.05,
                FCcutoff = 0.5,
                pointSize = 3.0,
                labSize = 6.0,
                title = 'Lymph:IPF vs HC')

#富集分析
diff1<-read.csv(file="/Users/houhouhou/opt/R/GSE128033/res1.csv") 
diff2<-read.csv(file="/Users/houhouhou/opt/R/GSE128033/res2.csv")
diff3<-read.csv(file="/Users/houhouhou/opt/R/GSE128033/res3.csv")
diff4<-read.csv(file="/Users/houhouhou/opt/R/GSE128033/res4.csv")
diff5<-read.csv(file="/Users/houhouhou/opt/R/GSE128033/res5.csv")

library(AnnotationDbi)
library(org.Hs.eg.db)#基因注释包
library(clusterProfiler)#富集包
library(dplyr)
library(ggplot2)#画图包

#2、基因id转换，kegg和go富集用的ID类型是ENTREZID）
#TCGA数据框如果没有进行基因注释，那么fromType应该是Ensembl，各种ID之间可以互相转换,toType可以是一个字符串，也可以是一个向量，看自己需求
gene.df1 <- bitr(diff1$SYMBOL,fromType="SYMBOL",toType="ENTREZID", OrgDb = org.Hs.eg.db)
gene1 <- gene.df1$ENTREZID

gene.df2 <- bitr(diff2$SYMBOL,fromType="SYMBOL",toType="ENTREZID", OrgDb = org.Hs.eg.db)
gene2 <- gene.df2$ENTREZID

gene.df3 <- bitr(diff3$SYMBOL,fromType="SYMBOL",toType="ENTREZID", OrgDb = org.Hs.eg.db)
gene3 <- gene.df3$ENTREZID

gene.df4 <- bitr(diff4$SYMBOL,fromType="SYMBOL",toType="ENTREZID", OrgDb = org.Hs.eg.db)
gene4 <- gene.df4$ENTREZID

gene.df5 <- bitr(diff5$SYMBOL,fromType="SYMBOL",toType="ENTREZID", OrgDb = org.Hs.eg.db)
gene5 <- gene.df5$ENTREZID
#3、GO富集
##CC表示细胞组分，MF表示分子功能，BP表示生物学过程，ALL表示同时富集三种过程，选自己需要的,我一般是做BP,MF,CC这3组再合并成一个数据框，方便后续摘取部分通路绘图。
ego_ALL <- enrichGO(gene = gene1,#我们上面定义了
                    OrgDb=org.Hs.eg.db,
                    keyType = "ENTREZID",
                    ont = "ALL",#富集的GO类型
                    pAdjustMethod = "BH",#这个不用管，一般都用的BH
                    minGSSize = 1,
                    pvalueCutoff = 0.01,#P值可以取0.05
                    qvalueCutoff = 0.05,
                    readable = TRUE)

ego_CC <- enrichGO(gene = gene1,
                   OrgDb=org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "CC",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.01,
                   qvalueCutoff = 0.05,
                   readable = TRUE)

ego_BP <- enrichGO(gene = gene1,
                   OrgDb=org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "BP",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.01,
                   qvalueCutoff = 0.05,
                   readable = TRUE)

ego_MF <- enrichGO(gene = gene1,
                   OrgDb=org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "MF",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.01,
                   qvalueCutoff = 0.05,
                   readable = TRUE)

#4、将结果保存到当前路径
ego_ALL <- as.data.frame(ego_ALL)
ego_result_BP <- as.data.frame(ego_BP)
ego_result_CC <- as.data.frame(ego_CC)
ego_result_MF <- as.data.frame(ego_MF)
ego <- rbind(ego_result_BP,ego_result_CC,ego_result_MF)#或者这样也能得到ego_ALL一样的结果
write.csv(ego_ALL,file = "ego_ALL.csv",row.names = T)
write.csv(ego_result_BP,file = "ego_result_BP.csv",row.names = T)
write.csv(ego_result_CC,file = "ego_result_CC.csv",row.names = T)
write.csv(ego_result_MF,file = "ego_result_MF.csv",row.names = T)
write.csv(ego,file = "ego.csv",row.names = T)

#5、但是有时候我们富集到的通路太多了，不方便绘制图片，可以选择部分绘制，这时候跳过第4步，来到第5步
display_number = c(10, 10, 10)#这三个数字分别代表选取的BP、CC、MF的通路条数，这个自己设置就行了
ego_result_BP <- as.data.frame(ego_BP)[1:display_number[1], ]
ego_result_CC <- as.data.frame(ego_CC)[1:display_number[2], ]
ego_result_MF <- as.data.frame(ego_MF)[1:display_number[3], ]

##将以上我们摘取的部分通路重新组合成数据框
go_enrich_df <- data.frame(
  ID=c(ego_result_BP$ID, ego_result_CC$ID, ego_result_MF$ID),
  Description=c(ego_result_BP$Description,ego_result_CC$Description,ego_result_MF$Description),
  GeneNumber=c(ego_result_BP$Count, ego_result_CC$Count, ego_result_MF$Count),
  type=factor(c(rep("biological process", display_number[1]), 
                rep("cellular component", display_number[2]),
                rep("molecular function", display_number[3])), 
              levels=c("biological process", "cellular component","molecular function" )))

##通路的名字太长了，我选取了通路的前五个单词作为通路的名字
for(i in 1:nrow(go_enrich_df)){
  description_splite=strsplit(go_enrich_df$Description[i],split = " ")
  description_collapse=paste(description_splite[[1]][1:5],collapse = " ") #这里的5就是指5个单词的意思，可以自己更改
  go_enrich_df$Description[i]=description_collapse
  go_enrich_df$Description=gsub(pattern = "NA","",go_enrich_df$Description)
}

##开始绘制GO柱状图
###横着的柱状图
go_enrich_df$type_order=factor(rev(as.integer(rownames(go_enrich_df))),labels=rev(go_enrich_df$Description))#这一步是必须的，为了让柱子按顺序显示，不至于很乱
COLS <- c("#66C3A5", "#8DA1CB", "#FD8D62")#设定颜色

ggplot(data=go_enrich_df, aes(x=type_order,y=GeneNumber, fill=type)) + #横纵轴取值
  geom_bar(stat="identity", width=0.8) + #柱状图的宽度，可以自己设置
  scale_fill_manual(values = COLS) + ###颜色
  coord_flip() + ##这一步是让柱状图横过来，不加的话柱状图是竖着的
  xlab("GO term") + 
  ylab("Gene_Number") + 
  labs(title = "The Most Enriched GO Terms")+
  theme_bw()

###竖着的柱状图 
go_enrich_df$type_order=factor(rev(as.integer(rownames(go_enrich_df))),labels=rev(go_enrich_df$Description))
COLS <- c("#66C3A5", "#8DA1CB", "#FD8D62")
ggplot(data=go_enrich_df, aes(x=type_order,y=GeneNumber, fill=type)) + 
  geom_bar(stat="identity", width=0.8) + 
  scale_fill_manual(values = COLS) + 
  theme_bw() + 
  xlab("GO term") + 
  ylab("Num of Genes") + 
  labs(title = "The Most Enriched GO Terms")+ 
  theme(axis.text.x=element_text(face = "bold", color="gray50",angle = 70,vjust = 1, hjust = 1 ))#angle是坐标轴字体倾斜的角度，可以自己设置

#内皮细胞KEGG
#1、KEGG富集
k1 <- enrichKEGG(gene = gene1,keyType = "kegg",organism= "human", qvalueCutoff = 0.01, pvalueCutoff=0.01,use_internal_data = TRUE)

#2、可视化
###柱状图
hh <- as.data.frame(k1)
rownames(hh) <- 1:nrow(hh)
hh$order=factor(rev(as.integer(rownames(hh))),labels = rev(hh$Description))
ggplot(hh,aes(y=order,x=Count,fill=p.adjust))+
  geom_bar(stat = "identity",width=0.7)+
  #coord_flip()+
  scale_fill_gradient(low = "darksalmon",high ="darkslateblue" )+
  labs(title = "KEGG Pathways Enrichment",
       x = "Gene numbers", 
       y = "Pathways")+
  theme(axis.title.x = element_text(face = "bold",size = 16),
        axis.title.y = element_text(face = "bold",size = 16),
        legend.title = element_text(face = "bold",size = 16))+
  theme_bw()

###气泡图
hh <- as.data.frame(kk)
rownames(hh) <- 1:nrow(hh)
hh$order=factor(rev(as.integer(rownames(hh))),labels = rev(hh$Description))
ggplot(hh,aes(y=order,x=Count))+
  geom_point(aes(size=Count,color=-1*p.adjust))+
  scale_color_gradient(low="green",high = "red")+
  labs(color=expression(p.adjust,size="Count"), 
       x="Gene Number",y="Pathways",title="KEGG Pathway Enrichment")+
  theme_bw()

#vein
vein<- NormalizeData(vein)
vein<- FindVariableFeatures(vein, nfeatures = 2000)
vein <- ScaleData(vein, vars.to.regress = c("S.Score", "G2M.Score"), verbose = FALSE)
vein <- RunPCA(vein, npcs = 50, verbose = FALSE)
vein <- RunHarmony(vein, reduction="pca",group.by.vars="group",reduction.save="harmony")
DimPlot(vein,reduction="pca")
ElbowPlot(vein, ndims=50, reduction="pca")
Loadings(object = vein[["pca"]])   ###获取基因在pc轴上的投射值
Embeddings(object = vein[["pca"]])  ###获取各个细胞的pc值
Stdev(IPF_combined)   ###获取各pc轴解释量方差
DimHeatmap(IPF_combined, dims = 1:20, nfeatures=10, cells = 500, balanced = TRUE)###查看决定pc值的top10基因在500个细胞中的热图
vein <- FindNeighbors(vein, reduction = "harmony", dims = 1:20)
vein <- FindClusters(vein,resolution=seq(from = 0.1, 
                                         to = 1.0, 
                                         by = 0.1))
vein <- RunUMAP(vein, reduction = "harmony", dims = 1:20)
vein<- RunTSNE(vein, reduction = "harmony", dims = 1:20)
library(clustree)
clustree(vein)
Idents(vein) <- "RNA_snn_res.0.2" 
DimPlot(vein,label =T)
DimPlot(vein,reduction="tsne",label=T)
vein$seurat_clusters <- vein@active.ident
DimPlot(vein,label = T,split.by = "group",ncol = 3) 
DimPlot(vein,label = T,reduction="tsne", split.by = "group",ncol = 3) 

DefaultAssay(vein) <- "RNA"
vein_markers  <- FindAllMarkers(vein, 
                                only.pos = TRUE,
                                min.pct = 0.1, 
                                logfc.threshold = 0.25)
sigvein_markers  <- vein_markers [vein_markers $p_val_adj < 0.2, ]
write.csv(sigvein_markers, file = "sigvein_markers.csv")

top_vein_markers <- vein_markers%>%group_by(cluster)%>%top_n(10, wt = avg_log2FC)
DoHeatmap(vein, features=top_vein_markers$gene, group.by ="seurat_clusters") 
