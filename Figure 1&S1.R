title: "Nociceptive interoception drives a self-reinforcing neurovascular circuit in aortic dissection"


# GSE254132 ----------------------------------------------------------
plan("multisession", workers = 20) #nbrOfWorkers();availableCores()
options(future.globals.maxSize = 20000 * 1024^2)
DefaultAssay(sc) <- "integrated"
sc <- FindVariableFeatures(sc, selection.method = "vst", nfeatures = 2000)
sc <- ScaleData(sc, features = VariableFeatures(sc))
sc <- RunPCA(sc, features = VariableFeatures(sc))
ElbowPlot(sc) 
sc <- RunUMAP(sc, dims = 1:15) 
sc <- FindNeighbors(sc, dims = 1:15)
sc <- FindClusters(sc, resolution = 0.1)
P1 <- DimPlot(sc, reduction = 'umap', raster = FALSE, label = T)
P2 <- DimPlot(sc, group.by = "orig.ident", reduction = "umap", raster = FALSE, label = T)
P1 + P2

P1|DotPlot(sc, feature = "MKI67") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

SMC = c("RAMP1","MYH11","TAGLN","ACTA2","MYL9","CSRP2","RGS5","TPM2","HSPB1","YNS1","NPNT","MYLK","TINAGL1","HSBP1","PDGFRA")
SynSMC = c("DCN","ACTA2","COL8A1")
ConSMC = c("MYH11","MYL9","ACTA2")
Fibroblast =  c("DCN","COL1A1","COL3A1","FGF7","MME","PDGFRA")

Idents(sc) <- sc$seurat_clusters
new.cluster.ids <- c("0"="Macrophage", 
                     "1"="Fibroblast", 
                     "2"="T", 
                     "4"="DC",
                     "3"="Neutrophil",
                     "5"="Endothelial",
                     "6"="Prolif",
                     "7"="Mast",
                     "8"="B",
                     "9"="SMC",
                     "10"="Erythrocyte"
)
sc <- RenameIdents(sc, new.cluster.ids) 
sc$major_cell_type <- sc@active.ident
table(sc$major_cell_type)
save(sc, file = "~/Aortic Dissection/GSE254132/AD_sc.Rdata")

sc_smc <- sc[,sc$major_cell_type %in% c("SMC")]
AD_GSE254132_smc <- sc_smc
save(AD_GSE254132_smc, file = "~/Aortic Dissection/GSE254132/AD_GSE254132_smc.Rdata")
dim(sc_smc)

cols = c("#FFC685","#99C5E3","#CC99BB","#638B66","#F5D8C1","#E57E7E","#A0A6A7","#CEB152","#26456E","#5586B3","#A3123A","#BB6449")
sc$major_cell_type <- factor(sc$major_cell_type, 
                             levels <- c("Macrophage", "Fibroblast", "T", "DC", "Neutrophil", "Endothelial",
                                         "Prolif", "Mast","B", "SMC", "Erythrocyte"))      
P1|DimPlot(sc, reduction = "umap", label = F, pt.size = 0.5, label.box = T, group.by = "major_cell_type", cols = cols)

DefaultAssay(sc) <- "RNA"
VlnPlot(sc, feature = cgrp_ramp1_genes, cols = color, group.by = "major_cell_type", pt.size = 0) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

VlnPlot(sc, feature = "RAMP1", cols = color, group.by = "major_cell_type", pt.size = 0) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

cgrp_ramp1_genes <- c("CALCA", "CALCB", "CALCRL", "RAMP1")
sc <- AddModuleScore(sc, features = list(cgrp_ramp1_genes), name = "CGRP_Pathway")
VlnPlot(sc, features = "CGRP_Pathway1", group.by = "major_cell_type", cols = cols, pt.size = 0)


library(AUCell)
receptor_geneset <- list("CGRP_Receptor" = c("CALCRL", "RAMP1"))
expr_matrix <- GetAssayData(sc, slot = "data")
rankings <- AUCell_buildRankings(expr_matrix, nCores = 10)
cells_AUC <- AUCell_calcAUC(receptor_geneset, rankings)
sc$CGRP_AUC <- as.vector(getAUC(cells_AUC)["CGRP_Receptor", ])
VlnPlot(sc, features = "CGRP_AUC", cols = color, group.by = "major_cell_type", pt.size = 0) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



# GSE222318 -----------------------------------------------------
library(Seurat)
library(data.table)
expr_mat <- fread("GSE222318_Raw_gene_counts_matrix.csv.gz", data.table = FALSE)
rownames(expr_mat) <- expr_mat[[1]]
expr_mat <- expr_mat[, -1]   
expr_mat <- as.matrix(expr_mat)
AD <- CreateSeuratObject(counts = expr_mat, project = "AD")

metadata <- read.csv("GSE222318_Cell_metadata.csv.gz", row.names = 1)
AD <- AddMetaData(AD, metadata = metadata)
head(AD@meta.data)
table(AD@meta.data$Position)

DefaultAssay(AD) <- "RNA"
Idents(AD2) <- "All"

VlnPlot(AD, c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0, group.by = NULL)
AD <- subset(AD, subset = nFeature_RNA < 5000 & nCount_RNA < 20000 & percent.mt < 15)
AD2 <- subset(AD, subset = nFeature_RNA < 4000 & nFeature_RNA > 200 & nCount_RNA > 200 & nCount_RNA < 15000 & percent.mt < 5)
dim(AD)
Assays(AD)

rb_genes <- rownames(AD)[grep("^RP[SL]",rownames(AD))]
counts <- GetAssayData(AD, assay = "RNA")
counts <- counts[-(which(rownames(AD) %in% rb_genes)),]
AD <- subset(AD, features = rownames(counts))

colnames(AD@meta.data)
head(AD@meta.data)
table(AD@meta.data$cellType1, AD@meta.data$subCluster)
table(AD@meta.data$cellType2, AD@meta.data$subCluster)
table(AD@meta.data$subCluster)

plan("multisession", workers = 10) #nbrOfWorkers();availableCores()
options(future.globals.maxSize = 20000 * 1024^2)

AD <- NormalizeData(AD, normalization.method = "LogNormalize", scale.factor = 10000)
AD <- FindVariableFeatures(AD, selection.method = "vst", nfeatures = 2000)
AD <- ScaleData(AD, features = VariableFeatures(AD))
AD <- RunPCA(AD, features = VariableFeatures(AD))
ElbowPlot(AD) 
AD <- RunUMAP(AD, dims = 1:10, 
              #n.neighbors = 15, 
              min.dist = 0.4,
              spread = 0.6
) 
AD <- FindNeighbors(AD, dims = 1:10)
AD <- FindClusters(AD, resolution = 0.1)
P1 <- DimPlot(AD, reduction = 'umap', raster = FALSE, label = T)
P2 <- DimPlot(AD, group.by = "orig.ident", reduction = "umap", raster = FALSE, label = T)
P1 + P2

AD_GSE222318 <- AD
save(AD_GSE222318, file = "~/Aortic Dissection/GSE222318/AD_GSE222318.Rdata")

p <- DimPlot(AD, reduction = 'umap', raster=FALSE, label= T)
color <- c("#FFC685","#8493CA","#E8CCC7", "#5586B3", 
           "#B6C796", "#C98C9A","#72A77A", "#F6A4A4", "#9F77CF",
           "#7D8BAE", "#638B66", 
           
           "#FFDD99",  "#D66F62")

DimPlot(AD, reduction = 'umap', cols = color, raster=FALSE, label= T) #+ NoLegend()

AD_markers <- FindAllMarkers(AD, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
AD_markers %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC) -> top20

p | DotPlot(AD, feature = Mast) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

Idents(AD) <- AD$seurat_clusters
interesting_groups <- c("0","1","2","3","5","6","7","8","9")
AD <- subset(x = AD, idents = interesting_groups)
new.cluster.ids <- c("0"="Macrophage", 
                     "1"="T", 
                     "2"="Fibroblast", 
                     "3"="SMC",
                     "5"="Neutrophil",
                     "6"="Endothelial",
                     "7"="Mast",
                     "8"="B",
                     "9"="Prolif"
)
AD <- RenameIdents(AD, new.cluster.ids) 
AD$major_cell_type <- AD@active.ident
table(AD$major_cell_type)
save(AD, file = "~/Aortic Dissection/AD_smc(GSE222318).Rdata")

AD_smc <- AD[,AD$major_cell_type %in% c("SMC")]
AD_GSE222318_smc <- AD_smc
save(AD_GSE222318_smc, file = "~/Aortic Dissection/GSE222318/AD_GSE222318_smc.Rdata")

SMC = c("RAMP1","MYH11","TAGLN","ACTA2","MYL9","CSRP2","RGS5","TPM2","HSPB1","YNS1","NPNT","MYLK","TINAGL1","HSBP1","PDGFRA")
SynSMC = c(AD_smc <- AD[,AD$major_cell_type %in% c("CD8_T","CD4_T","NK")]
           "DCN","ACTA2","COL8A1")
ConSMC = c("MYH11","MYL9","ACTA2")
Fibroblast =  c("DCN","COL1A1","COL3A1","FGF7","MME","PDGFRA")

Mph <- c("TNF","IL1B","IL6","IL12","LYVE1","CD68","CD80","CD86","SEPP1")
Mph <- c("MERTK","MRC1","IL10","F13A1","CD163","FOLR2","CX3CR1","SPP1","C1QA","C1QB","C1QC")
Mph <-c("C1QA","C1QB","C1QC")
Mph <- c("CD169","CD206","VCAM1","SPP1")
Mph <- c("SPI1","MYB","MAF","MAFB","ZEB2") # TRM
Monocyte <- c("LITAF","SOD2","IFITM2","S100A8","IFITM2","SAT1","NAMPT","MS3A4","CCR2","F13A1","FCGR4","CSF1R") # Monocyte
DC <- c("CCDC126","CLEC9A","IRF8","CPNE3","SRGN","FLT3","HLA-DOB") # DC
DC = c("CD209","ITGAM","ADGRE1","CCR7","FLT3","CCL22","IFITM1")
DC = c("CD123","CD304") # pDC in human
DC = c("IRF4","CD11C","IRF8","CD141","HLA-DRA","HLA-DPA1","HLA-DPB1") # mDC in human
DC = c("CLEC10A","CLEC10A","CD1C","HLA-DQA1","HLA-DPB1","HLA-DPA1","HLA-DRA") # DC in AD
Neutrophil = c("CSF3R","S100A8","FCGR3B","NCF1","CXCL8","CST3","NAMPT") # Neutrophil
CD4 = c("CD3D","CD3E")
CD4 = c("CREM","CXCR6","RGCC","MR3C1","GZMB") # CD4_active
CD4 = c("CCR7","IL7R","CCL20","KLRB1") # Resting CD4 T
CD4 = c("IL2RA","CTLA4","TNFRSF18","ID3","LTB") # CD4_Treg
CD8 = c("GZMK","CRTAM","CCL4","CMC1","CD8") # CD8_active
CD8 = c("GZMA",'FCGR3A','FGFBP2','TYROBP','CX3CR1','GNLY','KLRD1') # NKT
B = c('CD79A','CD79B',"CD19","MZB1") # B cell
Plasma = c("SSR4","SEC11C","PRDX4","FKBP11") # Plasma cell
Mast = c("MCPT4","GATA2","KIT") # Mast cell
Basophils = c("MCPT8","CPA3","CD200R3","FCER1") # Basophils
Endo = c('CLDN5','FLT1','RAMP2',"PECAM1") #Endo
ILC2 = c("KLRG1","AREG","RORA","GATA3") #ILC2
Prolif = c('MKI67','STMN1','TUBA1B','HIST1H4C')



p | DotPlot(AD, feature = "RAMP1", group.by = "Group2") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p | DotPlot(AD, feature = "RAMP1", group.by = "Group1") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

VlnPlot(AD, feature = cgrp_ramp1_genes, cols = color, group.by = "major_cell_type", pt.size = 0) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

cgrp_ramp1_genes <- c("CALCA", "CALCB", "CALCRL", "RAMP1",
                      "GNAS", "ADCY3", "ADCY5", "PRKACA", 
                      "CREB1", "PDE4D", "ARRB1", "ARRB2", 
                      "IL6", "TNF", "NR4A1", "PTGS2")
AD <- AddModuleScore(AD, features = list(cgrp_ramp1_genes), name = "CGRP_Pathway")
VlnPlot(AD, features = "CGRP_Pathway1", group.by = "major_cell_type", pt.size = 0)


library(AUCell)
receptor_geneset <- list("CGRP_Receptor" = c("CALCRL", "RAMP1"))
expr_matrix <- GetAssayData(AD, slot = "data")
rankings <- AUCell_buildRankings(expr_matrix, nCores = 1)
cells_AUC <- AUCell_calcAUC(receptor_geneset, rankings)
AD$CGRP_AUC <- as.vector(getAUC(cells_AUC)["CGRP_Receptor", ])
VlnPlot(AD, features = "CGRP_AUC", cols = color, group.by = "major_cell_type", pt.size = 0) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# GSE189795 ---------------------------------------------------------------
file_path <- path.expand("~/Aortic Dissection/GSE189795/GSE189795_RAW.tar")
cmd <- paste("tar -xvf", shQuote(file_path))
system(cmd)
data_dir <- "GSE189795_RAW"

library(data.table)
expr_file <- "~/Aortic Dissection/GSE189795/GSE189795_AllSample.NormalizedCounts.txt.gz"
expr_data <- fread(expr_file, data.table = FALSE)
head(expr_data)
head(expr_data[,1:5]) 
dim(expr_data)

rownames(expr_data) <- expr_data[,1]
expr_data <- expr_data[,-1]
expr_mat <- as.matrix(expr_data)
AD_GSE189795 <- CreateSeuratObject(counts = expr_mat)
print(AD_GSE189795)

plan("multisession", workers = 30) #nbrOfWorkers();availableCores()
options(future.globals.maxSize = 20000 * 1024^2)
AD_GSE189795 <- NormalizeData(AD_GSE189795)
AD_GSE189795 <- FindVariableFeatures(AD_GSE189795, selection.method = "vst", nfeatures = 2000)
AD_GSE189795 <- ScaleData(AD_GSE189795, features = VariableFeatures(AD_GSE189795))
AD_GSE189795 <- RunPCA(AD_GSE189795, features = VariableFeatures(AD_GSE189795))
ElbowPlot(AD_GSE189795) 
AD_GSE189795 <- RunUMAP(AD_GSE189795, dims = 1:15)
AD_GSE189795 <- FindNeighbors(AD_GSE189795, dims = 1:15)
AD_GSE189795 <- FindClusters(AD_GSE189795, resolution = 0.1)
P1 <- DimPlot(AD_GSE189795, reduction = 'umap', raster = FALSE, label = T)
P2 <- DimPlot(AD_GSE189795, group.by = "orig.ident", reduction = "umap", raster = FALSE, label = T)
P1 + P2
save(AD_GSE189795, file = "~/Aortic Dissection/GSE189795/AD_GSE189795.Rdata")


P1|DotPlot(AD_GSE189795, feature = "RAMP1") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

SMC = c("RAMP1","MYH11","TAGLN","ACTA2","MYL9","CSRP2","RGS5","TPM2","HSPB1","YNS1","NPNT","MYLK","TINAGL1","HSBP1","PDGFRA")
SynSMC = c("DCN","ACTA2","COL8A1")
ConSMC = c("MYH11","MYL9","ACTA2")
Fibroblast =  c("DCN","COL1A1","COL3A1","FGF7","MME","PDGFRA")

Idents(AD_GSE189795) <- AD_GSE189795$seurat_clusters
AD_GSE189795_smc <- AD_GSE189795[,AD_GSE189795$seurat_clusters %in% c("5")]
save(AD_GSE189795_smc, file = "~/Aortic Dissection/GSE189795/AD_GSE189795_smc.Rdata")
dim(AD_GSE189795_smc)



# GSE213740 ---------------------------------------------------------------
file_path <- path.expand("~/Aortic Dissection/GSE213740/GSE213740_RAW.tar")
cmd <- paste("tar -xvf", shQuote(file_path))
system(cmd)
data_dir <- "GSE213740_RAW"

library(Seurat)
library(Matrix)
library(dplyr)


samples <- list(
  "AD1" = "GSM6593315_AD_replicate_1",
  "AD2" = "GSM6593316_AD_replicate_2",
  "AD3" = "GSM6593317_AD_replicate_3",
  "AD4" = "GSM6593318_AD_replicate_4",
  "AD5" = "GSM6593319_AD_replicate_5",
  "AD6" = "GSM6593320_AD_replicate_6",
  "Normal1" = "GSM6593321_Normal_replicate_1",
  "Normal2" = "GSM6593322_Normal_replicate_2",
  "Normal3" = "GSM6593323_Normal_replicate_3"
)

seurat_list <- list()
for (i in names(samples)) {
  path_prefix <- samples[[i]]
  data_dir <- paste0("./", path_prefix)
  

  R.utils::gunzip(paste0(data_dir, "_barcodes.tsv.gz"), overwrite = TRUE)
  R.utils::gunzip(paste0(data_dir, "_features.tsv.gz"), overwrite = TRUE)
  R.utils::gunzip(paste0(data_dir, "_matrix.mtx.gz"), overwrite = TRUE)
  

  counts <- ReadMtx(
    mtx = paste0(data_dir, "_matrix.mtx"),
    features = paste0(data_dir, "_features.tsv"),
    cells = paste0(data_dir, "_barcodes.tsv"),
    feature.column = 2
  )
  
  seurat_obj <- CreateSeuratObject(counts = counts, project = i)
  seurat_obj$sample <- i
  seurat_list[[i]] <- seurat_obj
}

combined <- merge(seurat_list[[1]], y = seurat_list[-1], add.cell.ids = names(seurat_list), project = "GSE213740")
save(combined, file = "~/Aortic Dissection/GSE213740/AD_GSE213740_combined.Rdata")

plan("multisession", workers = 30) #nbrOfWorkers();availableCores()
options(future.globals.maxSize = 20000 * 1024^2)

combined[["percent.mt"]] <- PercentageFeatureSet(combined, pattern = "^MT-")
VlnPlot(combined, c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0)
combined <- subset(combined, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA > 200 & nCount_RNA < 20000 & percent.mt < 20)
dim(combined)
AD_GSE213740 <- combined
save(AD_GSE213740, file = "~/Aortic Dissection/GSE213740/AD_GSE213740.Rdata")

combined <- SCTransform(combined, method = "glmGamPoi",  verbose = TRUE)
combined <- RunPCA(combined)
ElbowPlot(combined)
combined <- RunUMAP(combined, dims = 1:20)
combined <- FindNeighbors(combined, dims = 1:20)
combined <- FindClusters(combined, resolution = 0.1)
P1 <- DimPlot(combined, reduction = 'umap', raster = FALSE, label = T)
P2 <- DimPlot(combined, group.by = "orig.ident", reduction = "umap", raster = FALSE, label = T)
P1 + P2

P1|DotPlot(combined, feature = SMC) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

SMC = c("RAMP1","MYH11","TAGLN","ACTA2","MYL9","CSRP2","RGS5","TPM2","HSPB1","NPNT","MYLK","TINAGL1","HSBP1","PDGFRA")
SynSMC = c("DCN","ACTA2","COL8A1")
ConSMC = c("MYH11","MYL9","ACTA2")
Fibroblast =  c("DCN","COL1A1","COL3A1","FGF7","MME","PDGFRA")

Idents(combined) <- combined$seurat_clusters
combined_smc <- combined[,combined$seurat_clusters %in% c("2","3")]
AD_GSE213740_smc <- combined_smc
save(AD_GSE213740_smc, file = "~/Aortic Dissection/GSE213740/AD_GSE213740_smc.Rdata")
dim(combined_smc)

# GSE224559 ---------------------------------------------------------------
setwd("~/Aortic Dissection/GSE224559")
file_path <- path.expand("~/Aortic Dissection/GSE224559/GSE224559_RAW.tar")
cmd <- paste("tar -xvf", shQuote(file_path))
system(cmd)

data_dir <- "~/Aortic Dissection/GSE224559"
files <- list.files(data_dir, full.names = TRUE)
mouse_files <- files[grepl("mouse", files)]
file.remove(mouse_files)

samples <- list(
  "AD3_Intima" = "GSM7025944_human_dissection_3_intima",
  "AD3_Adventitia" = "GSM7025945_human_dissection_3_adventitia",
  "AD6" = "GSM7025946_human_dissection_6"
)

seurat_list <- list()
for (i in names(samples)) {
  prefix <- samples[[i]]
  data_dir <- "~/Aortic Dissection/GSE224559"
  

  barcode_file <- file.path(data_dir, paste0(prefix, "_barcodes.tsv.gz"))
  feature_file <- file.path(data_dir, paste0(prefix, "_features.tsv.gz"))
  matrix_file <- file.path(data_dir, paste0(prefix, "_matrix.mtx.gz"))
  

  gunzip(barcode_file, overwrite = TRUE)
  gunzip(feature_file, overwrite = TRUE)
  gunzip(matrix_file, overwrite = TRUE)
  

  barcode_file_unzipped <- sub("\\.gz$", "", barcode_file)
  feature_file_unzipped <- sub("\\.gz$", "", feature_file)
  matrix_file_unzipped <- sub("\\.gz$", "", matrix_file)
  

  counts <- ReadMtx(
    mtx = matrix_file_unzipped,
    features = feature_file_unzipped,
    cells = barcode_file_unzipped,
    feature.column = 2
  )
  

  seurat_obj <- CreateSeuratObject(counts = counts, project = i)
  seurat_obj$sample <- i
  seurat_list[[i]] <- seurat_obj
}

combined <- merge(seurat_list[[1]], y = seurat_list[-1], add.cell.ids = names(seurat_list), project = "GSE224559_Human")
save(combined, file = "~/Aortic Dissection/GSE224559/AD_GSE224559_human_combined.Rdata")

combined[["percent.mt"]] <- PercentageFeatureSet(combined, pattern = "^MT-")
VlnPlot(combined, c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0)
combined <- subset(combined, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & nCount_RNA > 200 & nCount_RNA < 30000 & percent.mt < 10)
dim(combined)
AD_GSE224559 <- combined
save(AD_GSE224559, file = "~/Aortic Dissection/GSE224559/AD_GSE224559.Rdata")

plan("multisession", workers = 30) #nbrOfWorkers();availableCores()
options(future.globals.maxSize = 20000 * 1024^2)
combined <- SCTransform(combined, verbose = TRUE)
combined <- RunPCA(combined, verbose = FALSE)
ElbowPlot(combined)
combined <- RunUMAP(combined, dims = 1:15)
combined <- FindNeighbors(combined, dims = 1:15)
combined <- FindClusters(combined, resolution = 0.1)

P1 <- DimPlot(combined, reduction = 'umap', raster = FALSE, label = T)
P2 <- DimPlot(combined, group.by = "orig.ident", reduction = "umap", raster = FALSE, label = T)
P1 + P2

P1|DotPlot(combined, feature = Fibroblast) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

SMC = c("RAMP1","MYH11","TAGLN","ACTA2","MYL9","CSRP2","RGS5","TPM2","HSPB1","NPNT","MYLK","TINAGL1","HSBP1","PDGFRA")
SynSMC = c("DCN","ACTA2","COL8A1")
ConSMC = c("MYH11","MYL9","ACTA2")
Fibroblast =  c("DCN","COL1A1","COL3A1","FGF7","MME","PDGFRA")

Idents(combined) <- combined$seurat_clusters
combined_smc <- combined[,combined$seurat_clusters %in% c("7")]
AD_GSE224559_smc <- combined_smc
save(AD_GSE224559_smc, file = "~/Aortic Dissection/GSE224559/AAD_GSE224559_smc.Rdata")
dim(AD_GSE224559_smc)



# GSE226492 ---------------------------------------------------------------
setwd("~/Aortic Dissection/GSE226492")
file_path <- path.expand("~/Aortic Dissection/GSE226492/GSE226492_RAW.tar")
cmd <- paste("tar -xvf", shQuote(file_path))
system(cmd)

data_dir <- "~/Aortic Dissection/GSE226492"
files <- list.files(data_dir, full.names = TRUE)
aa_files <- files[grepl("AA", basename(files))]
file.remove(aa_files)
cat("\n", paste(basename(aa_files), collapse = "\n"))

data_dir <- "~/Aortic Dissection/GSE226492"
samples <- list(
  "Con1" = "GSM7077892_Con1",
  "Con2" = "GSM7077893_Con2",
  "Con3" = "GSM7077894_Con3",
  "AD1"  = "GSM7077898_AD1",
  "AD2"  = "GSM7077899_AD2",
  "AD3"  = "GSM7077900_AD3",
  "AD4_Marfan" = "GSM7077901_MF"
)

seurat_list <- list()
for (i in names(samples)) {
  prefix <- samples[[i]]
  

  barcode_file <- file.path(data_dir, paste0(prefix, "_barcodes.tsv.gz"))
  feature_file <- file.path(data_dir, paste0(prefix, "_features.tsv.gz"))
  matrix_file <- file.path(data_dir, paste0(prefix, "_matrix.mtx.gz"))
  

  if (file.exists(barcode_file)) gunzip(barcode_file, overwrite = TRUE)
  if (file.exists(feature_file)) gunzip(feature_file, overwrite = TRUE)
  if (file.exists(matrix_file)) gunzip(matrix_file, overwrite = TRUE)
  

  barcode_file_unzipped <- sub("\\.gz$", "", barcode_file)
  feature_file_unzipped <- sub("\\.gz$", "", feature_file)
  matrix_file_unzipped <- sub("\\.gz$", "", matrix_file)
  

  counts <- ReadMtx(
    mtx = matrix_file_unzipped,
    features = feature_file_unzipped,
    cells = barcode_file_unzipped,
    feature.column = 2
  )
  

  seurat_obj <- CreateSeuratObject(counts = counts, project = i)
  seurat_obj$sample <- i
  seurat_list[[i]] <- seurat_obj
}

combined <- merge(x = seurat_list[[1]], y = seurat_list[-1], add.cell.ids = names(seurat_list), project = "GSE226492")
combined$dataset <- "GSE226492"
save(combined, file = "~/Aortic Dissection/GSE226492/AD_GSE226492_combined_raw.Rdata")

combined[["percent.mt"]] <- PercentageFeatureSet(combined, pattern = "^MT-")
VlnPlot(combined, c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0)
combined <- subset(combined, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & nCount_RNA > 200 & nCount_RNA < 20000 & percent.mt < 10)
dim(combined)
AD_GSE226492 <- combined
save(AD_GSE226492, file = "~/Aortic Dissection/GSE226492/AD_GSE226492.Rdata")

plan("multisession", workers = 4)
options(future.globals.maxSize = 20000 * 1024^2)

combined <- SCTransform(combined, verbose = TRUE)
combined <- RunPCA(combined, verbose = FALSE)
ElbowPlot(combined)

combined <- RunUMAP(combined, dims = 1:15)
combined <- FindNeighbors(combined, dims = 1:15)
combined <- FindClusters(combined, resolution = 0.1)

P1 <- DimPlot(combined, reduction = 'umap', raster = FALSE, label = TRUE)
P2 <- DimPlot(combined, group.by = "orig.ident", reduction = "umap", raster = FALSE, label = TRUE)
P1 + P2

SMC <- c("RAMP1","MYH11","TAGLN","ACTA2","MYL9","CSRP2","RGS5","TPM2","HSPB1","NPNT","MYLK","TINAGL1","HSBP1","PDGFRA")
Fibroblast <-  c("DCN","COL1A1","COL3A1","FGF7","MME","PDGFRA")

DotPlot(combined, features = SMC) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

Idents(combined) <- combined$seurat_clusters
AD_GSE226492_smc <- combined[, combined$seurat_clusters %in% c("3")]
save(AD_GSE226492_smc, file = "~/Aortic Dissection/GSE226492/AD_GSE226492_smc.Rdata")
dim(AD_GSE226492_smc)



# All integration ---------------------------------------------------------
setwd("~/Aortic Dissection/All cluster")

load("~/Aortic Dissection/GSE254132/AD_GSE254132.Rdata")
load("~/Aortic Dissection/GSE222318/AD_GSE222318.Rdata")
load("~/Aortic Dissection/GSE189795/AD_GSE189795.Rdata")
load("~/Aortic Dissection/GSE213740/AD_GSE213740.Rdata")
load("~/Aortic Dissection/GSE226492/AD_GSE226492.Rdata")
load("~/Aortic Dissection/GSE224559/AD_GSE224559.Rdata")

AD_GSE254132$dataset <- "AD_GSE254132"
AD_GSE189795$dataset <- "AD_GSE189795"
AD_GSE213740$dataset <- "AD_GSE213740"
AD_GSE222318$dataset <- "AD_GSE222318"
AD_GSE224559$dataset <- "AD_GSE224559"
AD_GSE226492$dataset <- "AD_GSE226492"

seurat_list <- list(
  AD_GSE254132 = AD_GSE254132,
  AD_GSE189795 = AD_GSE189795,
  AD_GSE213740 = AD_GSE213740,
  AD_GSE222318 = AD_GSE222318,
  AD_GSE224559 = AD_GSE224559,
  AD_GSE226492 = AD_GSE226492
)

sapply(seurat_list, function(obj) ncol(obj))

plan("multisession", workers = 1)
options(future.globals.maxSize = 100 * 1024^3)

seurat_list <- lapply(seurat_list, SCTransform, verbose = FALSE)
features <- SelectIntegrationFeatures(object.list = seurat_list, nfeatures = 3000)
seurat_list <- PrepSCTIntegration(object.list = seurat_list, anchor.features = features)
seurat_list <- lapply(seurat_list,function(x){RunPCA(x,features = features,npcs = 50,verbose = FALSE)})
anchors <- FindIntegrationAnchors(object.list = seurat_list,
                                  normalization.method = "SCT",
                                  reduction = "rpca",
                                  anchor.features = features)

AD_all <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
save(AD_all, file = "~/Aortic Dissection/All cluster/AD_all.Rdata")

AD_all[["percent.mt"]] <- PercentageFeatureSet(AD_all, pattern = "^MT-")
VlnPlot(AD_all, c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0)
AD_all <- subset(AD_all, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & nCount_RNA > 200 & nCount_RNA < 20000)
dim(AD_all)
save(AD_all, file = "~/Aortic Dissection/All cluster/AD_all_260622.Rdata")

DefaultAssay(AD_all) <- "integrated"
AD_all <- FindVariableFeatures(AD_all,nfeatures = 3000)
AD_all <- ScaleData(AD_all,features = VariableFeatures(AD_all))
AD_all <- RunPCA(AD_all,npcs = 40,verbose = FALSE)
ElbowPlot(AD_all,ndims = 40)

library(ggplot2)
library(patchwork)
umap_list <- list()
for(d in 20:30){
  
  AD_tmp <- FindNeighbors(AD_all,dims = 1:d,verbose = FALSE)
  AD_tmp <- FindClusters(AD_tmp,resolution = 0.3,verbose = FALSE)
  AD_tmp <- RunUMAP(AD_tmp,dims = 1:d,verbose = FALSE)
  
  umap_list[[as.character(d)]] <- DimPlot(AD_tmp,group.by = "dataset") + ggtitle(paste0("dims = ",d))
}
wrap_plots(umap_list,ncol = 3)

p <- wrap_plots(umap_list, ncol = 3)
ggsave(filename = "~/Aortic Dissection/All cluster/UMAP_dims_20_30.pdf",plot = p,width = 16,height = 12)


AD_all <- FindNeighbors(AD_all,dims = 1:25)
AD_all <- FindClusters(AD_all,resolution = 0.2)
AD_all <- RunUMAP(AD_all,dims = 1:25)
save(AD_all, file = "~/Aortic Dissection/All cluster/AD_all_260623.Rdata")


p1 <- DimPlot(AD_all,cols = cluster_colors, label = TRUE, repel = TRUE, raster=FALSE) + NoLegend()
p1
p2 <- DimPlot(AD_all,group.by = "dataset")

DefaultAssay(AD_all) <- "RNA"

markers <- c(
  "FCN1","LST1","S100A8","S100A9","VCAN",
  "LYZ","C1QA","C1QB","C1QC","APOE","TREM2",
  "FCGR3B","CSF3R","MPO",
  "ITGAX","HLA-DRA","CLEC10A","LILRA4",
  "CD3D","CD3E","TRAC","CD8A","IL7R",
  "MS4A1","CD79A","CD74",
  "MZB1","JCHAIN","XBP1","IGKC","SDC1",
  "ACTA2","TAGLN","MYH11","CNN1","RAMP1",
  "COL1A1","COL1A2","DCN","LUM","COL3A1",
  "HBA1","HBA2","HBB"
)

cluster10_markers <- FindMarkers(
  AD_all,
  ident.1 = 10,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

head(cluster10_markers, 20)

DotPlot(AD_all, features = markers) + RotatedAxis()
VlnPlot(AD_all, features = markers) + RotatedAxis()

Idents(AD_all) <- "seurat_clusters"
new.cluster.ids <- c(
  "Monocyte",
  "T",
  "Fibroblast",
  "Macrophage",
  "Endothelial",
  "SMC",
  "SMC",
  "Neutrophil",
  "SMC",
  "Plasma",
  "Mast",
  "Fibroblast",
  "B",
  "Macrophage"
)
names(new.cluster.ids) <- levels(AD_all)
AD_all <- RenameIdents(AD_all,new.cluster.ids)
AD_all$major_cell_type <- Idents(AD_all)
table(AD_all$major_cell_type)

color = c("#5586B3","#7FB7BE","#9F77CF","#6FA8DC","#E8CCC7","#B6C796","#F6A4A4",
          "#FFDD99","#D66F62","#F2B36D", "#638B66","#72A77A", "#C8A2C8")

Idents(AD_all) <- "major_cell_type"
DimPlot(AD_all, cols = color, label=TRUE, raster=FALSE, repel=TRUE) + NoLegend()
DimPlot(AD_all, group.by="major_cell_type", label=TRUE, repel=TRUE)


VlnPlot(AD_all, feature = "RAMP1", cols = color, pt.size = 0.00000000000000001, raster=FALSE) + NoLegend()
VlnPlot(AD_all, feature = "RAMP1")

library(ggplot2)
df <- FetchData(AD_all,vars = c("RAMP1", "major_cell_type"))
ggplot(df, aes(x = major_cell_type, y = RAMP1, fill = major_cell_type)) +
  geom_boxplot(outlier.size = 0.2) +
  scale_fill_manual(values = color) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

library(scCustomize)
reds_dark <- c("#fee0d2", "#fc9272", "#de2d26", "#a50f15", "#67000d")
reds_continuous <- colorRampPalette(reds_dark)(100)
Plot_Density_Custom(seurat_object = AD_all, features = "RAMP1", reduction = "umap", custom_palette = reds_continuous, aspect_ratio = 1)



# SMC Integration -------------------------------------------------------------
AD_GSE254132_smc <- sc_smc
AD_GSE254132_smc$dataset <- "AD_GSE254132_smc"
AD_GSE189795_smc$dataset <- "AD_GSE189795_smc"
AD_GSE213740_smc$dataset <- "AD_GSE213740_smc"
AD_GSE222318_smc$dataset <- "AD_GSE222318_smc"
AD_GSE224559_smc$dataset <- "AD_GSE224559_smc"
AD_GSE226492_smc$dataset <- "AD_GSE226492_smc"

seurat_list <- list(
  AD_GSE254132_smc = AD_GSE254132_smc,
  AD_GSE189795_smc = AD_GSE189795_smc,
  AD_GSE213740_smc = AD_GSE213740_smc,
  AD_GSE222318_smc = AD_GSE222318_smc,
  AD_GSE224559_smc = AD_GSE224559_smc,
  AD_GSE226492_smc = AD_GSE226492_smc
)

sapply(seurat_list, function(obj) ncol(obj))

plan("multisession", workers = 1)
options(future.globals.maxSize = 20000 * 1024^2)


seurat_list <- lapply(seurat_list, SCTransform, verbose = FALSE)

features <- SelectIntegrationFeatures(object.list = seurat_list, nfeatures = 3000)
seurat_list <- PrepSCTIntegration(object.list = seurat_list, anchor.features = features)
anchors <- FindIntegrationAnchors(object.list = seurat_list,
                                  normalization.method = "SCT",
                                  anchor.features = features)

AD_SMC <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
save(AD_SMC, file = "~/Aortic Dissection/AD_SMC.Rdata")

AD_SMC[["percent.mt"]] <- PercentageFeatureSet(AD_SMC, pattern = "^MT-")
VlnPlot(AD_SMC, c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0)
AD_SMC <- subset(AD_SMC, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA > 200 & nCount_RNA < 20000 & percent.mt < 15)
dim(AD_SMC)
save(AD_SMC, file = "~/Aortic Dissection/AD_SMC_0721.Rdata")


DefaultAssay(AD_SMC) <- "integrated"



AD_SMC <- ScaleData(AD_SMC, verbose = FALSE)
AD_SMC <- RunPCA(AD_SMC, features = features, npcs = 50)
ElbowPlot(AD_SMC, ndims = 50)

plots <- list()
for (d in 10:15) {
  cat("Processing dims = 1:", d, "\n")
  
  AD_SMC <- RunUMAP(AD_SMC, dims = 1:d)
  AD_SMC <- FindNeighbors(AD_SMC, dims = 1:d)
  AD_SMC <- FindClusters(AD_SMC, resolution = 0.2)
  
  P1 <- DimPlot(AD_SMC, reduction = 'umap', raster = FALSE, label = TRUE) + ggtitle(paste("dims = 1:", d, sep = ""))
  P2 <- DimPlot(AD_SMC, group.by = "dataset", reduction = "umap", raster = FALSE)
  
  plots[[paste0("dims_", d)]] <- P1 + P2
}

wrap_plots(plots, ncol = 2)

DefaultAssay(AD_SMC) <- "integrated"
AD_SMC <- RunUMAP(AD_SMC, dims = 1:12)
AD_SMC <- FindNeighbors(AD_SMC, dims = 1:12)
AD_SMC <- FindClusters(AD_SMC, resolution = 0.3)


P1 <- DimPlot(AD_SMC, reduction = 'umap', raster = FALSE, label = TRUE)
P2 <- DimPlot(AD_SMC, group.by = "dataset", reduction = "umap", raster = FALSE)
P1 + P2

DefaultAssay(AD_SMC) <- "RNA"
SMC <- c("RAMP1","MYH11","TAGLN","ACTA2","MYL9","CSRP2","RGS5","TPM2","HSPB1","NPNT","MYLK","TINAGL1","HSBP1","PDGFRA")
Fibroblast <-  c("DCN","COL1A1","COL3A1","FGF7","MME","PDGFRA")
Neuropeptide_regulated_SMC <- c(
  "VIPR1", "VIPR2",
  "NPY1R", "NPY2R", "NPY5R",
  "SSTR2", "SSTR3", "SSTR5",
  "OXTR", "AVPR1A", "CHRM3",
  "CALCRL", "RAMP1",   #  CGRP receptor components
  "ADCY1", "PRKACA", "CREB1",
  "NOS1", "NOS3",
  "KCNMB1", "CACNA1C", "HCN4"
)

P1|DotPlot(AD_SMC, feature = Neuropeptide_regulated_SMC) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
P1|DotPlot(AD_SMC, feature = SMC, group.by = "dataset") + theme(axis.text.x = element_text(angle = 45, hjust = 1))


AD_SMC@meta.data$RNA_snn_res.0.5 <- NULL
AD_SMC@meta.data$RNA_snn_res.0.3 <- NULL
AD_SMC@meta.data$RNA_snn_res.0.2 <- NULL
AD_SMC@meta.data$Group1 <- NULL
AD_SMC@meta.data$Group2 <- NULL
AD_SMC@meta.data$major_cell_type <- NULL
AD_SMC@meta.data$integrated_snn_res.0.1 <- NULL
AD_SMC@meta.data$CGRP_Pathway1 <- NULL
AD_SMC@meta.data$CGRP_AUC <- NULL
AD_SMC@meta.data$scissor <- NULL
AD_SMC@meta.data$celltype <- NULL
AD_SMC@meta.data$integrated_snn_res.0.5 <- NULL
AD_SMC@meta.data$seurat_clusters <- NULL
AD_SMC@meta.data$celltype <- NULL
head(AD_SMC@meta.data)
table(AD_SMC@meta.data$dataset)
save(AD_SMC, file = "~/Aortic Dissection/AD_SMC_0721.Rdata")


# Cluster -----------------------------------------------------------------
head(AD_SMC@meta.data)
table(AD_SMC@meta.data$dataset)
table(AD_SMC@meta.data$Group, AD_SMC@meta.data$integrated_snn_res.0.3)
table(AD_SMC@meta.data$Subgroup, AD_SMC@meta.data$integrated_snn_res.0.3)

Assays(AD_SMC)
DefaultAssay(AD_SMC) <- "RNA"


color = c("#5586B3", "#638B66", "#B6C796", "#E8CCC7", "#FFC685", "#C98C9A", "#9F77CF", "#72A77A", "#FFDD99", "#8493CA", "#F6A4A4", "#D66F62", "#C8A2C8")

color = c("#C9A9A6", "#D7BEA8", "#B7B7A4", "#A5B463", "#94A7A1", "#81A29D", "#6B9080", "#778DA9", "#9290C3", "#A4907C", "#B98B82", "#C8A2C8", "#F6A4A4")

P1 <- DimPlot(AD_SMC, reduction = 'umap', cols = color, raster = FALSE, label = TRUE)
P2 <- DimPlot(AD_SMC, group.by = "dataset", reduction = "umap", raster = FALSE)
P1 + P2

DefaultAssay(AD_SMC) <- "RNA"
Idents(AD_SMC) <- AD_SMC$integrated_snn_res.0.3
AD_SMC_markers <- FindAllMarkers(AD_SMC, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
library(dplyr)
top30 <- AD_SMC_markers %>% group_by(cluster) %>% top_n(30, avg_log2FC)
save(top30, file = "~/Aortic Dissection/AD_SMC_allmarkers.csv")

Contractile_VSMC = c("ACTA2", "TAGLN", "MYH11", "CNN1", "TPM2", "MYL9")
Synthetic_VSMC = c("COL1A1", "COL1A2", "FN1", "SPP1", "VCAN", "TIMP1")
Secretory_VSMC <- c("FN1", "SPP1", "VCAN", "TIMP1", "TGFB1", "CCL2", "IL6")
Fibromyocyte <- c("ACTA2", "TAGLN", "DCN", "LUM", "COL3A1", "PDGFRA", "FBLN1", "FBLN2")
Mesenchymal_like_VSMC = c("PDGFRA", "DCN", "LUM", "THY1", "COL3A1")
Fibroblast_like_VSMC= c("FBLN1", "FBLN2", "COL6A1", "PDGFRA")

Osteogenic_VSMC = c("RUNX2", "BGLAP", "SPP1", "ALPL", "BMP2")
Chondrogenic_VSMC = c("SOX9", "ACAN", "COL2A1", "COMP")

Macrophage_like_VSMC = c("CD68", "CD14", "LYZ", "AIF1")
Proinflammatory_VSMC = c("CCL2", "CXCL2", "IL6", "TNF", "NFKBIA", "ICAM1")

Endothelial_like_VSMC = c("VWF", "PECAM1", "CDH5", "KDR","VCAM1")

Stem_like_VSMC = c("NES", "PROM1", "SOX2", "KLF4")
Pericytes = c("RGS5", "PDGFRB", "CSPG4", "ACTA2")

Neuroresponsive_VSMC <- c("RAMP1","PRKACA","CREB1","VIPR2","VIPR1","CALCRL","ADCYAP1R1","NPY1R")
IFN_SMC <- c(

  "ISG15", "IFI6", "IFI27", "IFIT1", "IFIT2", "IFIT3", "IRF7", "OAS1", "OAS2", "OAS3",
  "MX1", "MX2", "STAT1", "STAT2", "RSAD2", "XAF1", "BST2",

  "DDX58", "IFI44", "IFI44L", "TRIM22", "HERC5", "EPSTI1",

  "CXCL9", "CXCL10", "CXCL11", "HLA-A", "HLA-B", "HLA-C", "HLA-DRA"
)

P1|DotPlot(AD_SMC, feature = c("NES","IL6")) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
P1|DotPlot(AD_SMC, feature = Proinflammatory_VSMC) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

save(AD_SMC, file = "~/Aortic Dissection/AD_SMC_0721.Rdata")


interesting_groups <- c('0','1','2','3','4','5','6','7','8','9','10')
Idents(AD_SMC) <- AD_SMC@meta.data$integrated_snn_res.0.3
AD_smc <- subset(x = AD_SMC, idents = interesting_groups)
table(AD_smc$Group, AD_smc$integrated_snn_res.0.3)


Idents(AD_smc) <- AD_smc@meta.data$integrated_snn_res.0.3
new.cluster.ids <- c("0"="Synthetic SMC",
                     "1"="Mesenchymal-like SMC", 
                     "2"="Immune-regulatory SMC", 
                     "3"="IFN-enriched SMC",
                     "4"="Proinflammatory SMC",
                     "5"="Neuroresponsive SMC",
                     "6"="Contractile SMC",
                     "7"="Stressed SMC",
                     "8"="Secretory SMC",
                     "9"="Fibroblast-like SMC",
                     "10"="Macrophage-like SMC")
AD_smc <- RenameIdents(AD_smc, new.cluster.ids) 
AD_smc$major_cell_type <- AD_smc@active.ident
table(AD_smc$major_cell_type)
AD_SMC <- AD_smc
save(AD_SMC, file = "~/Aortic Dissection/AD_SMC_0723.Rdata")



# Plot - UMAP / Cell Ratio / Heatmap --------------------------------------------------------------
head(AD_SMC@meta.data)
table(AD_SMC@meta.data$dataset)
table(AD_SMC@meta.data$Group, AD_SMC@meta.data$major_cell_type)
table(AD_SMC@meta.data$Subgroup, AD_SMC@meta.data$major_cell_type)
P1 <- DimPlot(AD_SMC, group.by = "integrated_snn_res.0.3", reduction = "umap", label = T)

color1 = c("#5586B3", "#638B66", "#B6C796", "#9F77CF", "#ead4d4",  "#C98C9A", "#fad4a6","#72A77A", "#cee0e6", "#9ac1f0", "#f5c0bf", "#D66F62", "#C8A2C8")
color2 = c("#5586B3", "#81A29D", "#B6C796", "#9290C3", "#ead4d4",  "#C98C9A", "#fad4a6","#72A77A", "#cee0e6", "#9ac1f0", "#f5c0bf", "#D66F62", "#C8A2C8")
color3 = c("#5586B3", "#45958e", "#B6C796", "#9F77CF", "#ead4d4",  "#C98C9A", "#fad4a6","#72A77A", "#cee0e6", "#9ac1f0", "#f5c0bf", "#D66F62", "#C8A2C8")

P1|DimPlot(AD_SMC, reduction = 'umap', cols = color3, raster = FALSE, label = F)
p1 <- DimPlot(AD_SMC, reduction = 'umap', cols = color1, raster = FALSE, label = F)
p2 <- DimPlot(AD_SMC, reduction = 'umap', cols = color2, raster = FALSE, label = F)
p3 <- DimPlot(AD_SMC, reduction = 'umap', cols = color3, raster = FALSE, label = F)
p1|p2|p3


DimPlot(AD_SMC, reduction = 'umap', cols = color3, raster = FALSE, label = T)

DimPlot(AD_SMC, reduction = 'umap', group.by = "Subgroup", raster = FALSE, label = T)
DimPlot(AD_SMC, reduction = "umap", split.by = "Subgroup", cols = color3, label = TRUE, raster = FALSE)
AD_SMC@meta.data$Subgroup <- ordered(AD_SMC@meta.data$Subgroup,
                                     levels = c("Chronic","Subacute","Acute","CT"))
table(AD_SMC@meta.data$Subgroup)
Idents(AD_SMC) <- AD_SMC@meta.data$major_cell_type
prop.table(table(Idents(AD_SMC)))
table(Idents(AD_SMC), AD_SMC@meta.data$Subgroup)
Cellratio <- prop.table(table(Idents(AD_SMC), AD_SMC@meta.data$Subgroup), margin = 2)
Cellratio <- as.data.frame(Cellratio)

ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1), stat = "identity", width = 0.5, size = 0.4) + 
  theme_classic() +
  labs(x='', y = 'Ratio') +
  scale_fill_manual(values = color) +
  theme(
    panel.border = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_text(color = "black", size = 14),
    axis.title = element_text(color = "black", size = 14),
    legend.text = element_text(color = "black", size = 14),
    legend.title = element_text(color = "black", size = 14),
    plot.title = element_text(color = "black", size = 14, hjust = 0.5)
  ) +
  coord_flip()


DimPlot(AD_SMC, reduction = 'umap', group.by = "dataset", cols = color3, raster = FALSE, label = T)
DimPlot(AD_SMC, reduction = "umap", split.by = "dataset", cols = color3, label = TRUE, raster = FALSE)

table(AD_SMC@meta.data$dataset)
Idents(AD_SMC) <- AD_SMC@meta.data$major_cell_type
prop.table(table(Idents(AD_SMC)))
table(Idents(AD_SMC), AD_SMC@meta.data$dataset)
Cellratio <- prop.table(table(Idents(AD_SMC), AD_SMC@meta.data$dataset), margin = 2)
Cellratio <- as.data.frame(Cellratio)

ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1), stat = "identity", width = 0.5, size = 0.4) + 
  theme_classic() +
  labs(x='', y = 'Ratio') +
  scale_fill_manual(values = color3) +
  theme(
    panel.border = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_text(color = "black", size = 14),
    axis.title = element_text(color = "black", size = 14),
    legend.text = element_text(color = "black", size = 14),
    legend.title = element_text(color = "black", size = 14),
    plot.title = element_text(color = "black", size = 14, hjust = 0.5)
  ) +
  coord_flip()

## Gene expression Hheatmap
mito_genes <- rownames(AD_SMC)[grep("^MT-", rownames(AD_SMC))] 
rb_genes <- rownames(AD_SMC)[grep("^RP[SL]",rownames(AD_SMC))]
counts <- GetAssayData(AD_SMC, assay = "RNA")
counts <- counts[-(which(rownames(AD_SMC) %in% mito_genes)),]
counts <- counts[-(which(rownames(AD_SMC) %in% rb_genes)),]
AD_SMC <- subset(AD_SMC, features = rownames(counts))

AD_SMC <- ScaleData(AD_SMC)
table(AD_SMC@meta.data$major_cell_type)
Idents(AD_SMC) <- AD_SMC$major_cell_type
AD_SMC_markers <- FindAllMarkers(AD_SMC, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

AD_SMC_markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top5

DoHeatmap(AD_SMC, features = top20$gene, group.colors = color3, size = 2, draw.lines = F) +
  scale_fill_gradientn(colors = c("white","#ff9600"))

DoHeatmap(AD_SMC, features = top20$gene, group.colors = color3, size = 2, draw.lines = F) +
  scale_fill_gradientn(colors = c("white", "white", "#FFD580", "#FF9600"), 
                       values = c(0, 0.4, 0.6, 1))

DoHeatmap(AD_SMC, features = top20$gene, group.colors = color3,  size = 2, draw.lines = FALSE) +
  scale_fill_gradientn(colors = c("white", "white", "#A8E6A3"), 
                       values = c(0, 0.4, 0.6, 1))


top_genes <- unique(top5$gene)


DotPlot(AD_SMC,
        features = top_genes,
        group.by = "major_cell_type",
        cols = c("white", "#F08080"),
        assay = "SCT") +
  RotatedAxis() +
  theme_test() +
  theme(
    axis.text = element_text(size = 16, colour = "black"),
    axis.title = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

## Gene expression DotPlot
desired_order <- c("Contractile SMC", "Fibroblast-like SMC", "IFN-enriched SMC", 
                   "Immune-regulatory SMC", "Macrophage-like SMC", "Mesenchymal-like SMC",
                   "Neuroresponsive SMC", "Proinflammatory SMC", "Secretory SMC",
                   "Stressed SMC", "Synthetic SMC")
AD_SMC@meta.data$major_cell_type <- ordered(AD_SMC@meta.data$major_cell_type, levels = desired_order)

desired_order <- rev(c("Contractile SMC", "Fibroblast-like SMC", "IFN-enriched SMC", 
                       "Immune-regulatory SMC", "Macrophage-like SMC", "Mesenchymal-like SMC",
                       "Neuroresponsive SMC", "Proinflammatory SMC", "Secretory SMC",
                       "Stressed SMC", "Synthetic SMC"))
AD_SMC$major_cell_type <- factor(AD_SMC$major_cell_type, levels = desired_order)

color3 = c("#5586B3", "#45958e", "#B6C796", "#9F77CF", "#ead4d4",  "#C98C9A", "#fad4a6","#72A77A", "#cee0e6", "#9ac1f0", "#f5c0bf")
color_vec <- c("#e76e50", "#e18d64", "#e6b471", "#e4c57b", "#bcbc80","#a5b47e", "#92af81", "#70a68b", "#4d9a8e", "#3a6f6f", "#334a55")

color3_cont <- colorRampPalette(unname(rev(color3)))(100)



color3_cont <- colorRampPalette(unname(color3))(100)

DotPlot(AD_SMC,group.by = "major_cell_type",assay = "SCT",cols = color3_cont,
        features = c(
          "CNN1", "ACTA2", "TAGLN",
          "DCN", "FBLN1","COL1A1",
          "ISG15", "IFI6","MX1",
          "ZEB2","APOE","STEAP4","CCL21",
          "CD74", "HLA-DRA","SPP1",
          "TMSB4X","TNS1","NFIA",
          "ACTG1","TGM2","SLC39A14","CP",
          "EEF1G","SEPTIN7","CCL5", 
          "CCL2","ADAMTS4","XBP1",
          "MT1M","RGS16","DDIT3",
          "FN1","VCAN","TIMP1"
        )) + 
  RotatedAxis() +
  theme_test() +
  scale_x_discrete(position = "top") +
  scale_colour_gradientn(colours = color3_cont) +
  theme(
    axis.text = element_text(size = 16, colour = "black"),
    axis.title = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    axis.text.x = element_text(angle = 270, hjust = 1, vjust = 0)
  )

library(pheatmap)

desired_order <- c("Contractile SMC", "Fibroblast-like SMC", "IFN-enriched SMC", 
                   "Immune-regulatory SMC", "Macrophage-like SMC", "Mesenchymal-like SMC",
                   "Neuroresponsive SMC", "Proinflammatory SMC", "Secretory SMC",
                   "Stressed SMC", "Synthetic SMC")

avg_exp <- AverageExpression(
  AD_SMC,
  group.by = "major_cell_type",
  assays = "SCT",
  features = rev(c(
    "CNN1", "ACTA2", "TAGLN",
    "DCN", "FBLN1","COL1A1",
    "ISG15", "IFI6","MX1",
    "ZEB2","APOE","STEAP4","CCL21",
    "CD74", "HLA-DRA","SPP1",
    "TMSB4X","TNS1","NFIA",
    "ACTG1","TGM2","SLC39A14","CP",
    "EEF1G","SEPTIN7","CCL5", 
    "CCL2","ADAMTS4","XBP1",
    "MT1M","RGS16","DDIT3",
    "FN1","VCAN","TIMP1"
  )
  ))

mat <- avg_exp$SCT
mat_ordered <- mat[, desired_order]

library(pheatmap)
pheatmap(
  mat_ordered,
  color = colorRampPalette(rev(c( "#e18d64", "#ffffff", "#92af81")))(100),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  scale = "row",
  fontsize = 12,
  border_color = NA
)

library(pheatmap)

pheatmap(
  mat_ordered,
  color = colorRampPalette(
    c(
      "#6F79AE",
      "#A7D2DA",
      "#F3F1EE",
      "#E7B98A",
      "#E07B4F"
    )
  )(50),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  scale = "row",
  fontsize = 12,
  border_color = NA
)


## KEGG
mito_genes <- rownames(AD_SMC)[grep("^MT-", rownames(AD_SMC))] 
rb_genes <- rownames(AD_SMC)[grep("^RP[SL]",rownames(AD_SMC))]
counts <- GetAssayData(AD_SMC, assay = "RNA")
counts <- counts[-(which(rownames(AD_SMC) %in% mito_genes)),]
counts <- counts[-(which(rownames(AD_SMC) %in% rb_genes)),]
AD_SMC <- subset(AD_SMC, features = rownames(counts))

Idents(AD_SMC) <- AD_SMC$major_cell_type
AD_SMC_markers <- FindAllMarkers(AD_SMC, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

library(clusterProfiler)
library(org.Hs.eg.db)
ids = bitr(AD_SMC_markers$gene,'SYMBOL','ENTREZID','org.Hs.eg.db')
AD_SMC_markers = merge(AD_SMC_markers,ids,by.x='gene',by.y='SYMBOL')

AD_SMC_Synthetic <- subset(AD_SMC_markers, AD_SMC_markers$cluster %in% c("Synthetic SMC"))
AD_SMC_Mesenchymal <- subset(AD_SMC_markers, AD_SMC_markers$cluster %in% c("Mesenchymal-like SMC"))
AD_SMC_Immune <- subset(AD_SMC_markers, AD_SMC_markers$cluster %in% c("Immune-regulatory SMC"))
AD_SMC_IFN <- subset(AD_SMC_markers, AD_SMC_markers$cluster %in% c("IFN-enriched SMC"))
AD_SMC_Proinflammatory <- subset(AD_SMC_markers, AD_SMC_markers$cluster %in% c("Proinflammatory SMC"))
AD_SMC_Neuroresponsive <- subset(AD_SMC_markers, AD_SMC_markers$cluster %in% c("Neuroresponsive SMC"))
AD_SMC_Contractile <- subset(AD_SMC_markers, AD_SMC_markers$cluster %in% c("Contractile SMC"))
AD_SMC_Stressed <- subset(AD_SMC_markers, AD_SMC_markers$cluster %in% c("Stressed SMC"))
AD_SMC_Secretory <- subset(AD_SMC_markers, AD_SMC_markers$cluster %in% c("Secretory SMC"))
AD_SMC_Fibroblast <- subset(AD_SMC_markers, AD_SMC_markers$cluster %in% c("Fibroblast-like SMC"))
AD_SMC_Macrophage <- subset(AD_SMC_markers, AD_SMC_markers$cluster %in% c("Macrophage-like SMC"))

enrichKK = enrichKEGG(gene =  AD_SMC_Macrophage$ENTREZID,
                      organism     = 'hsa',
                      pvalueCutoff = 0.01,
                      qvalueCutoff = 0.01)
enrichKK = setReadable(enrichKK,OrgDb = 'org.Hs.eg.db',keyType = "ENTREZID")
barplot(enrichKK)
kegg_Macrophage = enrichKK@result
kegg_Macrophage$log10 = -log10(kegg_Macrophage$pvalue)
kegg_Macrophage = kegg_Macrophage[order(kegg_Macrophage$log10),]
kegg_Macrophage$group = 'kegg_Macrophage'
write.csv(kegg_Macrophage, file = "~/Aortic Dissection/KEGG/kegg_Macrophage.csv")


read_kegg_csv <- function(path, category) {
  df <- read.csv(path, stringsAsFactors = FALSE)
  df$category <- category
  return(df)
}

df_list <- list(
  read_kegg_csv("~/Aortic Dissection/KEGG/Revised/2pathway/kegg_Mesenchymal_inR.csv", "Mesenchymal"),
  read_kegg_csv("~/Aortic Dissection/KEGG/Revised/2pathway/kegg_Immune_inR.csv", "Immune"),
  read_kegg_csv("~/Aortic Dissection/KEGG/Revised/2pathway/kegg_IFN_inR.csv", "IFN"),
  read_kegg_csv("~/Aortic Dissection/KEGG/Revised/2pathway/kegg_Proinflammatory_inR.csv", "Proinflammatory"),
  read_kegg_csv("~/Aortic Dissection/KEGG/Revised/2pathway/kegg_Neuroresponsive_inR.csv", "Neuroresponsive"),
  read_kegg_csv("~/Aortic Dissection/KEGG/Revised/2pathway/kegg_Contractile_inR.csv", "Contractile"),
  read_kegg_csv("~/Aortic Dissection/KEGG/Revised/2pathway/kegg_Stressed_inR.csv", "Stressed"),
  read_kegg_csv("~/Aortic Dissection/KEGG/Revised/2pathway/kegg_Secretory_inR.csv", "Secretory"),
  read_kegg_csv("~/Aortic Dissection/KEGG/Revised/2pathway/kegg_Fibroblast_inR.csv", "Fibroblast"),
  read_kegg_csv("~/Aortic Dissection/KEGG/Revised/2pathway/kegg_Macrophage_inR.csv", "Macrophage"),
  read_kegg_csv("~/Aortic Dissection/KEGG/Revised/2pathway/kegg_Synthetic_inR.csv", "Synthetic")
  
)

kegg_combined <- do.call(rbind, df_list)
head(kegg_combined)
table(kegg_combined$category)

category_order <- c(
  "Contractile", "Fibroblast", "IFN", "Immune", "Macrophage",
  "Mesenchymal", "Neuroresponsive", "Proinflammatory", "Secretory", 
  "Stressed", "Synthetic"
)

kegg_combined$category <- factor(kegg_combined$category, levels = category_order)
kegg_combined <- kegg_combined %>% arrange(category)

ggplot(kegg_combined, aes(category, factor(Description, levels = unique(Description)))) +
  geom_point(aes(color = log10, size = GeneRatio)) + 
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.line = element_line(color = "black")
  ) +
  scale_color_gradient(low = "#f87b54", high = "#4a3374") +
  labs(x = NULL, y = NULL) + 
  guides(size = guide_legend(order = 1))


# SMC miloR - Group ---------------------------------------------------------------
## miloR
library(stringr)
library(Seurat)
library(miloR)
library(scater)
library(scran)
library(dplyr)
library(patchwork)
library(qs)
library(BiocParallel)



umap_coords <- Embeddings(AD_SMC, reduction = "umap")
colnames(umap_coords) <- c("UMAP1", "UMAP2")
AD_SMC_sce <- as.SingleCellExperiment(AD_SMC)
reducedDim(AD_SMC_sce, "UMAP") <- umap_coords


AD_SMC_milo <- Milo(AD_SMC_sce)
reducedDim(AD_SMC_milo,"UMAP") <- reducedDim(AD_SMC_sce,"UMAP")
AD_SMC_milo <- runPCA(AD_SMC_milo, ncomponents = 50)
AD_SMC_milo
AD_SMC_milo <- buildGraph(AD_SMC_milo, k = 30, d = 15,reduced.dim = "PCA")
AD_SMC_milo <- makeNhoods(AD_SMC_milo, prop = 0.2, k = 30, d = 15, refined = TRUE, reduced_dims = "PCA")
plotNhoodSizeHist(AD_SMC_milo)

AD_SMC_milo <- countCells(AD_SMC_milo, meta.data = as.data.frame(colData(AD_SMC_milo)), sample = "orig.ident")
head(nhoodCounts(AD_SMC_milo))

traj_design <- data.frame(colData(AD_SMC_milo))[,c("orig.ident", "Group")]
traj_design <- distinct(traj_design)
rownames(traj_design) <- traj_design$orig.ident
traj_design <- traj_design[colnames(nhoodCounts(AD_SMC_milo)), , drop=FALSE]
traj_design$group <- factor(traj_design$Group, level = c("CT","AD"))
traj_design

AD_SMC_milo <- calcNhoodDistance(AD_SMC_milo, d=15, reduced.dim = "PCA")
da_results <- testNhoods(AD_SMC_milo,design = ~ Group, design.df = traj_design)
da_results %>%
  arrange(SpatialFDR) %>%
  head()

ggplot(da_results, aes(PValue)) + geom_histogram(bins=50)
ggplot(da_results, aes(logFC, -log10(SpatialFDR))) +
  geom_point() +
  geom_hline(yintercept = 1)

AD_SMC_milo <- buildNhoodGraph(AD_SMC_milo)
umap_pl <- plotReducedDim(AD_SMC_milo, dimred = "UMAP",
                          colour_by="Group", text_by = "major_cell_type",
                          text_size=3, point_size=0.5) + guides(fill="none")
umap_pl

nh_graph_pl <- plotNhoodGraphDA(AD_SMC_milo, da_results, layout="UMAP",alpha=0.1)
umap_pl + nh_graph_pl +
  plot_layout(guides="collect")

nh_graph_pl +
  plot_layout(guides="collect")

da_results <- annotateNhoods(AD_SMC_milo, da_results, coldata_col = "major_cell_type")
head(da_results)
ggplot(da_results, aes(minor_DC_state_fraction)) + geom_histogram(bins=50)

str(da_results)
table(da_results$major_cell_type)
range(da_results$SpatialFDR)

desired_order <- c("Contractile SMC", "Fibroblast-like SMC", "IFN-enriched SMC", 
                   "Immune-regulatory SMC", "Macrophage-like SMC", "Mesenchymal-like SMC",
                   "Neuroresponsive SMC", "Proinflammatory SMC", "Secretory SMC",
                   "Stressed SMC", "Synthetic SMC")
da_results$major_cell_type <- factor(da_results$major_cell_type,
                                     levels = desired_order)
plotDAbeeswarm(da_results, group.by = "major_cell_type", alpha = 0.1)



# SMC Ro/e - Subgroup ----------------------------------------------------------------
# Ro/e
### +++, Ro/e > 1; ++, 0.8 < Ro/e  1; +, 0.2  Ro/e  0.8; +/, 0 < Ro/e < 0.2; , Ro/e = 0
do.tissueDist <- function(cellInfo.tb = cellInfo.tb,
                          meta.cluster = cellInfo.tb$meta.cluster,
                          colname.patient = "patient",
                          loc = cellInfo.tb$loc,
                          out.prefix,
                          pdf.width=3,
                          pdf.height=5,
                          verbose=0){

    library(data.table)
  dir.create(dirname(out.prefix),F,T)
  
  cellInfo.tb = data.table(cellInfo.tb)
  cellInfo.tb$meta.cluster = as.character(meta.cluster)
  
  if(is.factor(loc)){
    cellInfo.tb$loc = loc
  }else{cellInfo.tb$loc = as.factor(loc)}
  
  loc.avai.vec <- levels(cellInfo.tb[["loc"]])
  count.dist <- unclass(cellInfo.tb[,table(meta.cluster,loc)])[,loc.avai.vec]
  freq.dist <- sweep(count.dist,1,rowSums(count.dist),"/")
  freq.dist.bin <- floor(freq.dist * 100 / 10)
  print(freq.dist.bin)
  
  {
    count.dist.melt.ext.tb <- test.dist.table(count.dist)
    p.dist.tb <- dcast(count.dist.melt.ext.tb,rid~cid,value.var="p.value")
    OR.dist.tb <- dcast(count.dist.melt.ext.tb,rid~cid,value.var="OR")
    OR.dist.mtx <- as.matrix(OR.dist.tb[,-1])
    rownames(OR.dist.mtx) <- OR.dist.tb[[1]]
  }
  
  sscVis::plotMatrix.simple(OR.dist.mtx,
                            out.prefix=sprintf("%s.OR.dist",out.prefix),
                            show.number=F,
                            waterfall.row=T,par.warterfall = list(score.alpha = 2,do.norm=T),
                            exp.name=expression(italic(OR)),
                            z.hi=4,
                            palatte=viridis::viridis(7),
                            pdf.width = 4, pdf.height = pdf.height)
  if(verbose==1){
    return(list("count.dist.melt.ext.tb"=count.dist.melt.ext.tb,
                "p.dist.tb"=p.dist.tb,
                "OR.dist.tb"=OR.dist.tb,
                "OR.dist.mtx"=OR.dist.mtx))
  }else{
    return(OR.dist.mtx)
  }
}

test.dist.table <- function(count.dist,min.rowSum=0){
  count.dist <- count.dist[rowSums(count.dist)>=min.rowSum,,drop=F]
  sum.col <- colSums(count.dist)
  sum.row <- rowSums(count.dist)
  count.dist.tb <- as.data.frame(count.dist)
  setDT(count.dist.tb,keep.rownames=T)
  count.dist.melt.tb <- melt(count.dist.tb,id.vars="rn")
  colnames(count.dist.melt.tb) <- c("rid","cid","count")
  count.dist.melt.ext.tb <- as.data.table(ldply(seq_len(nrow(count.dist.melt.tb)), function(i){
    this.row <- count.dist.melt.tb$rid[i]
    this.col <- count.dist.melt.tb$cid[i]
    this.c <- count.dist.melt.tb$count[i]
    other.col.c <- sum.col[this.col]-this.c
    this.m <- matrix(c(this.c,
                       sum.row[this.row]-this.c,
                       other.col.c,
                       sum(sum.col)-sum.row[this.row]-other.col.c),
                     ncol=2)
    res.test <- fisher.test(this.m)
    data.frame(rid=this.row,
               cid=this.col,
               p.value=res.test$p.value,
               OR=res.test$estimate)
  }))
  count.dist.melt.ext.tb <- merge(count.dist.melt.tb,count.dist.melt.ext.tb,
                                  by=c("rid","cid"))
  count.dist.melt.ext.tb[,adj.p.value:=p.adjust(p.value,"BH")]
  return(count.dist.melt.ext.tb)
}

meta.tb = AD_SMC@meta.data
meta.tb$meta.cluster = meta.tb$major_cell_type
meta.tb$loc = meta.tb$Subgroup
table(meta.tb$loc)

OR.NEU.list <- do.tissueDist(cellInfo.tb = meta.tb,
                             out.prefix = "./ORplot",
                             pdf.width = 4, pdf.height = 6, verbose = 1)
OR.NEU.list_p = OR.NEU.list$p.dist.tb
OR.NEU.list_heat = OR.NEU.list$OR.dist.mtx
OR.NEU.list_bin = OR.NEU.list$freq.dist.bin

head(OR.NEU.list_heat)
OR.NEU.list_heat_t <- t(OR.NEU.list_heat)
OR.NEU.list_heat_t <- OR.NEU.list_heat_t[
  c("CT", "Acute", "Subacute", "Chronic"),
  , drop = FALSE
]

orderedcell = pheatmap::pheatmap(OR.NEU.list_heat_t,
                                 show_colnames = T,
                                 scale = "row", #'none', 'row' or 'column'
                                 cluster_rows = F, 
                                 cluster_cols = F, 
                                 border_color = "black", border_width = 2,
                                 cellwidth = 30, cellheight = 24,
                                 display_numbers = T, fontsize_number = 10, number_color = "black",
                                 color = colorRampPalette(c("#a6bddb", "#f0f0f0", "#fbb4ae"))(100))


# SMC bootstrap -----------------------------------------------------------
table(AD_SMC@meta.data$major_cell_type)
table(AD_SMC@meta.data$Subgroup)

meta <- AD_SMC@meta.data
table(meta$major_cell_type, meta$Subgroup)

library(dplyr)
library(boot)


bootstrap_proportion <- function(data, group_col, target_col, group_value, target_value, R = 1000) {

  df_sub <- data[data[[group_col]] == group_value, ]
  n <- nrow(df_sub)
  
  boot_fn <- function(data, indices) {
    sampled <- data[indices, ]
    mean(sampled[[target_col]] == target_value)
  }
  
  boot_obj <- boot(df_sub, boot_fn, R = R)
  mean_prop <- mean(boot_obj$t)
  ci <- quantile(boot_obj$t, probs = c(0.025, 0.975))
  
  return(data.frame(
    Subgroup = group_value,
    major_cell_type = target_value,
    proportion = mean_prop,
    lower = ci[1],
    upper = ci[2]
  ))
}


subgroups <- c("CT", "Acute", "Subacute", "Chronic")
cell_types <- unique(meta$major_cell_type)


bootstrap_results <- do.call(rbind, lapply(cell_types, function(cell_type) {
  do.call(rbind, lapply(subgroups, function(subgroup) {
    bootstrap_proportion(meta, "Subgroup", "major_cell_type", subgroup, cell_type, R = 500)
  }))
}))


table(AD_SMC@meta.data$major_cell_type)
target_cell_type <- "Mesenchymal-like SMC"
color3 = c("#5586B3", "#45958e", "#B6C796", "#9F77CF", "#ead4d4",  "#C98C9A", "#fad4a6","#72A77A", "#cee0e6", "#9ac1f0", "#f5c0bf")


target_df <- bootstrap_results %>%
  filter(major_cell_type == target_cell_type) %>%
  mutate(Subgroup = factor(Subgroup, levels = c("CT", "Acute", "Subacute", "Chronic")))


missing_subgroups <- setdiff(levels(target_df$Subgroup), target_df$Subgroup)
if (length(missing_subgroups) > 0) {
  fill_df <- data.frame(
    Subgroup = factor(missing_subgroups, levels = c("CT", "Acute", "Subacute", "Chronic")),
    major_cell_type = target_cell_type,
    proportion = 0,
    lower = 0,
    upper = 0
  )
  target_df <- bind_rows(target_df, fill_df)
}


target_df$Subgroup <- factor(target_df$Subgroup, levels = c("CT", "Acute", "Subacute", "Chronic"))


ggplot(target_df, aes(x = Subgroup, y = proportion)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, group = 1), fill = "#c6e0dc", alpha = 0.3) +
  geom_line(group = 1, color = "#45958e", size = 1.2) +
  geom_point(color = "#45958e", size = 3) +
  labs(
    y = "Proportion (bootstrap)",
    x = "",
    title = paste(target_cell_type)
  ) +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black", linewidth = 0.8),
    
    axis.text = element_text(color = "black", size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title = element_text(color = "black", size = 16),
    plot.title = element_text(color = "black", size = 18, hjust = 0.5)
  )


# SMC Neuro subtype -------------------------------------------------------

neuro_receptor_sets <- list(
  
  Sympathetic_Adrenergic = c("ADRA1A", "ADRA1B", "ADRA1D",
                             "ADRA2A", "ADRA2B", "ADRA2C",
                             "ADRB1", "ADRB2", "ADRB3"),
  
  
  Parasympathetic_Muscarinic = c("CHRM1", "CHRM2", "CHRM3", "CHRM4", "CHRM5"),
  
  
  Sensory_Neuropeptide = c("TACR1",
                           "CALCRL",
                           "RAMP1",
                           "NTRK1", "NTRK2", "NTRK3",
                           "TRPV1", "TRPA1")
)


library(UCell)
AD_SMC <- AddModuleScore_UCell(AD_SMC,features = neuro_receptor_sets,name = NULL)
colnames(AD_SMC@meta.data)

library(ggplot2)
VlnPlot(AD_SMC, features = "Sympathetic_Adrenergic", pt.size = 0)+ NoLegend()
VlnPlot(AD_SMC, features = "Parasympathetic_Muscarinic", pt.size = 0)+ NoLegend()
VlnPlot(AD_SMC, features = "Sensory_Neuropeptide", pt.size = 0)+ NoLegend()

library(ggplot2)
library(tidyr)
library(dplyr)

library(data.table)
heatmap_dt <- as.data.table(AD_SMC@meta.data)
heatmap_data <- heatmap_dt[, .(Sympathetic_Adrenergic = mean(Sympathetic_Adrenergic, na.rm = TRUE),
                               Parasympathetic_Muscarinic = mean(Parasympathetic_Muscarinic, na.rm = TRUE),
                               Sensory_Neuropeptide = mean(Sensory_Neuropeptide, na.rm = TRUE)),
                           by = major_cell_type] %>%
  melt(id.vars = "major_cell_type",
       variable.name = "Receptor_Type",
       value.name = "Score")


ggplot(heatmap_data, aes(x = Receptor_Type, y = major_cell_type, fill = Score)) +
  geom_tile(color = "white", size = 0.5) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = mean(heatmap_data$Score)) +
  labs(x = "Neuroreceptor Type", 
       y = "Cell Type", 
       title = "Average Neuroreceptor Scores by Cell Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5))

# SMC Nebula --------------------------------------------------------------
library(Seurat)
library(Nebulosa)


library(scCustomize)
reds_dark <- c("#fee0d2", "#fc9272", "#de2d26", "#a50f15", "#67000d")
reds_continuous <- colorRampPalette(reds_dark)(100)
Plot_Density_Custom(seurat_object = AD_SMC, features = "RAMP1", reduction = "umap", custom_palette = reds_continuous, aspect_ratio = 1)
Plot_Density_Custom(seurat_object = AD_SMC, features = "CALCRL", reduction = "umap", custom_palette = reds_continuous, aspect_ratio = 1)


