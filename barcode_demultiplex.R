
library(ggplot2)
library(Seurat)
library(cowplot)
library(harmony)
library(dplyr)


cc05c=Read10X_h5("/Volumes/singlecell1/hashing\ realign/cc05ashdep_p/outs/filtered_feature_bc_matrix.h5")
cc05=CreateSeuratObject(counts=cc05c$`Gene Expression`)

# Normalize RNA data with log normalization
cc05 <- NormalizeData(cc05)
# Find and scale variable features
cc05 <- FindVariableFeatures(cc05, selection.method = "mean.var.plot")
cc05 <- ScaleData(cc05, features = VariableFeatures(cc05))


cc05[["HTO"]] <- CreateAssayObject(counts = cc05c$`Antibody Capture`)
cc05 <- NormalizeData(cc05, assay = "HTO", normalization.method = "CLR")

cc05 <- HTODemux(cc05, assay = "HTO", positive.quantile = 0.99)

table(cc05$HTO_classification.global)


table(Idents(cc05))
cc05=RenameIdents(cc05, "HTO1"="HTII-280+", "HTO2"="HTII-280-")
cc05$ht2=Idents(cc05)
cc05$sample="cc05"
cc05$group="Donor"
saveRDS(cc05, file="cc05.rds")

Idents(cc05) <- "HTO_classification.global"
VlnPlot(cc05, features = "nCount_RNA", pt.size = 0.1, log = TRUE)


cc05.subset <- subset(cc05, idents = "Negative", invert = TRUE)

# Calculate a tSNE embedding of the HTO data
DefaultAssay(cc05.subset) <- "HTO"
cc05.subset <- ScaleData(cc05.subset, features = rownames(cc05.subset),
                         verbose = FALSE)
cc05.subset <- RunPCA(cc05.subset, features = rownames(cc05.subset), approx = FALSE)
cc05.subset <- RunTSNE(cc05.subset, dims = 1:8, perplexity = 100)
DimPlot(cc05.subset)
HTOHeatmap(cc05, assay = "HTO", ncells = 5000)

# Extract the singlets
cc05.singlet <- subset(cc05, idents = "Singlet")
HTOHeatmap(cc05.singlet, assay = "HTO", ncells = 5000)


saveRDS(cc05.singlet, file="cc05.singlet.rds")



library(ggplot2)
library(Seurat)
library(cowplot)
library(harmony)
library(dplyr)


cc06c=Read10X_h5("/Volumes/singlecell1/hashing\ realign/cc06hashdep_p/outs/filtered_feature_bc_matrix.h5")
cc06=CreateSeuratObject(counts=cc06c$`Gene Expression`)

# Normalize RNA data with log normalization
cc06 <- NormalizeData(cc06)
# Find and scale variable features
cc06 <- FindVariableFeatures(cc06, selection.method = "mean.var.plot")
cc06 <- ScaleData(cc06, features = VariableFeatures(cc06))


cc06[["HTO"]] <- CreateAssayObject(counts = cc06c$`Antibody Capture`)
cc06 <- NormalizeData(cc06, assay = "HTO", normalization.method = "CLR")

cc06 <- HTODemux(cc06, assay = "HTO", positive.quantile = 0.99)

table(cc06$HTO_classification.global)






###cc06 <- HTODemux(cc06, assay = "HTO", positive.quantile = 0.95)
###table(cc06$HTO_classification.global)
Idents(cc06) <- "HTO_maxID"
RidgePlot(cc06, assay = "HTO", features = c("HTO1", "HTO2"), ncol = 2)

FeatureScatter(cc06, feature1 = "HTO1", feature2 = "HTO2")





table(Idents(cc06))
cc06=RenameIdents(cc06, "HTO1"="HTII-280+", "HTO2"="HTII-280-")
cc06$ht2=Idents(cc06)
cc06$sample="cc06"
cc06$group="Donor"
saveRDS(cc06, file="cc06.rds")

Idents(cc06) <- "HTO_classification.global"
VlnPlot(cc06, features = "nCount_RNA", pt.size = 0.1, log = TRUE)


cc06.subset <- subset(cc06, idents = "Negative", invert = TRUE)

# Calculate a tSNE embedding of the HTO data
DefaultAssay(cc06.subset) <- "HTO"
cc06.subset <- ScaleData(cc06.subset, features = rownames(cc06.subset),
                         verbose = FALSE)
cc06.subset <- RunPCA(cc06.subset, features = rownames(cc06.subset), approx = FALSE)
cc06.subset <- RunTSNE(cc06.subset, dims = 1:8, perplexity = 100)
DimPlot(cc06.subset)
HTOHeatmap(cc06, assay = "HTO", ncells = 5000)

# Extract the singlets
cc06.singlet <- subset(cc06, idents = "Singlet")
HTOHeatmap(cc06.singlet, assay = "HTO", ncells = 5000)


saveRDS(cc06.singlet, file="cc06.singlet.rds")



library(ggplot2)
library(Seurat)
library(cowplot)
library(harmony)
library(dplyr)


ipf07fepc=Read10X_h5("/Volumes/singlecell1/hashing\ realign/ipf07fep/outs/filtered_feature_bc_matrix.h5")
ipf07fep=CreateSeuratObject(counts=ipf07fepc$`Gene Expression`)

# Normalize RNA data with log normalization
ipf07fep <- NormalizeData(ipf07fep)
# Find and scale variable features
ipf07fep <- FindVariableFeatures(ipf07fep, selection.method = "mean.var.plot")
ipf07fep <- ScaleData(ipf07fep, features = VariableFeatures(ipf07fep))


ipf07fep[["HTO"]] <- CreateAssayObject(counts = ipf07fepc$`Antibody Capture`)
ipf07fep <- NormalizeData(ipf07fep, assay = "HTO", normalization.method = "CLR")

ipf07fep <- HTODemux(ipf07fep, assay = "HTO", positive.quantile = 0.99)

table(ipf07fep$HTO_classification.global)






###ipf07fep <- HTODemux(ipf07fep, assay = "HTO", positive.quantile = 0.95)
###table(ipf07fep$HTO_classification.global)
Idents(ipf07fep) <- "HTO_maxID"
RidgePlot(ipf07fep, assay = "HTO", features = c("HTO1", "HTO2"), ncol = 2)

FeatureScatter(ipf07fep, feature1 = "HTO1", feature2 = "HTO2")





table(Idents(ipf07fep))
ipf07fep=RenameIdents(ipf07fep, "HTO1"="HTII-280+", "HTO2"="HTII-280-")
ipf07fep$ht2=Idents(ipf07fep)
ipf07fep$sample="ipf07fep"
ipf07fep$group="IPF"
saveRDS(ipf07fep, file="ipf07fep.rds")

Idents(ipf07fep) <- "HTO_classification.global"
VlnPlot(ipf07fep, features = "nCount_RNA", pt.size = 0.1, log = TRUE)


ipf07fep.subset <- subset(ipf07fep, idents = "Negative", invert = TRUE)

# Calculate a tSNE embedding of the HTO data
DefaultAssay(ipf07fep.subset) <- "HTO"
ipf07fep.subset <- ScaleData(ipf07fep.subset, features = rownames(ipf07fep.subset),
                             verbose = FALSE)
ipf07fep.subset <- RunPCA(ipf07fep.subset, features = rownames(ipf07fep.subset), approx = FALSE)
ipf07fep.subset <- RunTSNE(ipf07fep.subset, dims = 1:8, perplexity = 100)
DimPlot(ipf07fep.subset)
HTOHeatmap(ipf07fep, assay = "HTO", ncells = 5000)

# Extract the singlets
ipf07fep.singlet <- subset(ipf07fep, idents = "Singlet")
HTOHeatmap(ipf07fep.singlet, assay = "HTO", ncells = 5000)


saveRDS(ipf07fep.singlet, file="ipf07fep.singlet.rds")



library(ggplot2)
library(Seurat)
library(cowplot)
library(harmony)
library(dplyr)


ipf04fepc=Read10X_h5("/Volumes/singlecell1/hashing\ realign/ipf04fephash/outs/filtered_feature_bc_matrix.h5")
ipf04fep=CreateSeuratObject(counts=ipf04fepc$`Gene Expression`)

# Normalize RNA data with log normalization
ipf04fep <- NormalizeData(ipf04fep)
# Find and scale variable features
ipf04fep <- FindVariableFeatures(ipf04fep, selection.method = "mean.var.plot")
ipf04fep <- ScaleData(ipf04fep, features = VariableFeatures(ipf04fep))


ipf04fep[["HTO"]] <- CreateAssayObject(counts = ipf04fepc$`Antibody Capture`)
ipf04fep <- NormalizeData(ipf04fep, assay = "HTO", normalization.method = "CLR")

ipf04fep <- HTODemux(ipf04fep, assay = "HTO", positive.quantile = 0.99)

table(ipf04fep$HTO_classification.global)






###ipf04fep <- HTODemux(ipf04fep, assay = "HTO", positive.quantile = 0.95)
###table(ipf04fep$HTO_classification.global)
Idents(ipf04fep) <- "HTO_maxID"
RidgePlot(ipf04fep, assay = "HTO", features = c("HTO1", "HTO2"), ncol = 2)

FeatureScatter(ipf04fep, feature1 = "HTO1", feature2 = "HTO2")





table(Idents(ipf04fep))
ipf04fep=RenameIdents(ipf04fep, "HTO1"="HTII-280+", "HTO2"="HTII-280-")
ipf04fep$ht2=Idents(ipf04fep)
ipf04fep$sample="ipf04fep"
ipf04fep$group="IPF"
saveRDS(ipf04fep, file="ipf04fep.rds")

Idents(ipf04fep) <- "HTO_classification.global"
VlnPlot(ipf04fep, features = "nCount_RNA", pt.size = 0.1, log = TRUE)


ipf04fep.subset <- subset(ipf04fep, idents = "Negative", invert = TRUE)

# Calculate a tSNE embedding of the HTO data
DefaultAssay(ipf04fep.subset) <- "HTO"
ipf04fep.subset <- ScaleData(ipf04fep.subset, features = rownames(ipf04fep.subset),
                             verbose = FALSE)
ipf04fep.subset <- RunPCA(ipf04fep.subset, features = rownames(ipf04fep.subset), approx = FALSE)
ipf04fep.subset <- RunTSNE(ipf04fep.subset, dims = 1:8, perplexity = 100)
DimPlot(ipf04fep.subset)
HTOHeatmap(ipf04fep, assay = "HTO", ncells = 5000)

# Extract the singlets
ipf04fep.singlet <- subset(ipf04fep, idents = "Singlet")
HTOHeatmap(ipf04fep.singlet, assay = "HTO", ncells = 5000)


saveRDS(ipf04fep.singlet, file="ipf04fep.singlet.rds")





library(ggplot2)
library(Seurat)
library(cowplot)
library(harmony)
library(dplyr)


ipf04nepc=Read10X_h5("/Volumes/singlecell1/hashing\ realign/ipf04nephash/outs/filtered_feature_bc_matrix.h5")
ipf04nep=CreateSeuratObject(counts=ipf04nepc$`Gene Expression`)

# Normalize RNA data with log normalization
ipf04nep <- NormalizeData(ipf04nep)
# Find and scale variable features
ipf04nep <- FindVariableFeatures(ipf04nep, selection.method = "mean.var.plot")
ipf04nep <- ScaleData(ipf04nep, features = VariableFeatures(ipf04nep))


ipf04nep[["HTO"]] <- CreateAssayObject(counts = ipf04nepc$`Antibody Capture`)
ipf04nep <- NormalizeData(ipf04nep, assay = "HTO", normalization.method = "CLR")

ipf04nep <- HTODemux(ipf04nep, assay = "HTO", positive.quantile = 0.99)

table(ipf04nep$HTO_classification.global)






###ipf04nep <- HTODemux(ipf04nep, assay = "HTO", positive.quantile = 0.95)
###table(ipf04nep$HTO_classification.global)
Idents(ipf04nep) <- "HTO_maxID"
RidgePlot(ipf04nep, assay = "HTO", features = c("HTO1", "HTO2"), ncol = 2)

FeatureScatter(ipf04nep, feature1 = "HTO3", feature2 = "HTO4")





table(Idents(ipf04nep))
ipf04nep=RenameIdents(ipf04nep, "HTO3"="HTII-280+", "HTO4"="HTII-280-")
ipf04nep$ht2=Idents(ipf04nep)
ipf04nep$sample="ipf04nep"
ipf04nep$group="IPF"
saveRDS(ipf04nep, file="ipf04nep.rds")

Idents(ipf04nep) <- "HTO_classification.global"
VlnPlot(ipf04nep, features = "nCount_RNA", pt.size = 0.1, log = TRUE)


ipf04nep.subset <- subset(ipf04nep, idents = "Negative", invert = TRUE)

# Calculate a tSNE embedding of the HTO data
DefaultAssay(ipf04nep.subset) <- "HTO"
ipf04nep.subset <- ScaleData(ipf04nep.subset, features = rownames(ipf04nep.subset),
                             verbose = FALSE)
ipf04nep.subset <- RunPCA(ipf04nep.subset, features = rownames(ipf04nep.subset), approx = FALSE)
ipf04nep.subset <- RunTSNE(ipf04nep.subset, dims = 1:8, perplexity = 100)
DimPlot(ipf04nep.subset)
HTOHeatmap(ipf04nep, assay = "HTO", ncells = 5000)

# Extract the singlets
ipf04nep.singlet <- subset(ipf04nep, idents = "Singlet")
HTOHeatmap(ipf04nep.singlet, assay = "HTO", ncells = 5000)


saveRDS(ipf04nep.singlet, file="ipf04nep.singlet.rds")

levels(hash)
DimPlot(hash)

DimPlot(hash, group.by="ht2")


















