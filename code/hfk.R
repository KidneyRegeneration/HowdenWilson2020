library(Seurat)
library(tidyverse)
library(ggraph)
source(here::here("functions.R"))

# load samples after downloading from geo

#------------Generate Seurat Objects-----------------
options(future.globals.maxSize = Inf)
# sample1 <- Read10x()
# sample2 <- Read10x()

spence.male <- CreateSeuratObject(counts = sample1)
spence.male <- AddMetaData(object = spence.male, metadata = "SPENCE_male", col.name = "Sample")
spence.female <- CreateSeuratObject(sample2)
spence.female <- AddMetaData(object = spence.female, metadata = "SPENCE_female", col.name = "Sample")

# generate list for both samples
spence.list <- list(spence.male, spence.female)
spence.list <- list(spence.int[, spence.int$Sample == "SPENCE_male"],
                    spence.int[, spence.int$Sample == "SPENCE_female"])

#------------Quality Control and Dimensional Reduction-----------------
for (i in 1:length(spence.list)) {
  spence.list[[i]] <- CellCycleScoring(spence.list[[i]], s.features = cc.genes$s.genes,
                                       g2m.features = cc.genes$g2m.genes, set.ident = F)
  spence.list[[i]] <- PercentageFeatureSet(spence.list[[i]], pattern = "^MT-", col.name = "percent.mt")
  ribo <- c(grep(pattern = "^RPS", x = rownames(spence.list[[i]]), value = T),
            grep(pattern = "^RPL", rownames(spence.list[[i]]), value = T))
  spence.list[[i]] <- PercentageFeatureSet(spence.list[[i]], features = ribo, col.name = "percent.ribo")
  spence.list[[i]] <- SCTransform(spence.list[[i]], verbose = F, vars.to.regress = c("S.Score", "G2M.Score"))

}

features <- SelectIntegrationFeatures(spence.list, nfeatures = 3000)
spence.list <- PrepSCTIntegration(spence.list, anchor.features = features, verbose = F)
anchors <- FindIntegrationAnchors(spence.list, normalization.method = "SCT", anchor.features = features, verbose = F)
spence.int <- IntegrateData(anchorset = anchors,
                            normalization.method = "SCT", verbose = T)
spence.int <- RunPCA(spence.int)

# 3d plot
spence.int <- RunUMAP(spence.int, dims = 1:30, seed.use = 250395, n.components = 3)
spence.int@reductions$umap3d <- spence.int@reductions$umap
# 2d plot
spence.int <- RunUMAP(spence.int, dims = 1:30, seed.use = 250395, n.components = 2)

#------------Clustering------------------------------

DefaultAssay(spence.int) <- "integrated"
spence.int <- subset(spence.int, subset = nFeature_RNA > 1500 & nCount_RNA < 2e+05 & percent.mt < 40)
spence.int <- FindNeighbors(spence.int, dims = 1:30)
spence.int <- FindClusters(spence.int, resolution = seq(0, 0.5, 0.05))

spence.int$all_stringent <- spence.int$integrated_snn_res.0.4
spence.int$all_broad <- spence.int$integrated_snn_res.0.1

DefaultAssay(spence.int) <- "RNA"
spence.int <- NormalizeData(spence.int, verbose = FALSE)

spence.int$Broad_Annotation<- as.character(spence.int$all_broad)

spence.int$Broad_Annotation <- gsub(pattern = "1",
                                    replacement = "Nephron_Progenitors",
                                    x = spence.int$Broad_Annotation)
spence.int$Broad_Annotation <- gsub(pattern = "2",
                                    replacement = "Podocytes",
                                    x = spence.int$Broad_Annotation)
spence.int$Broad_Annotation <- gsub(pattern = "0",
                                    replacement = "Stroma",
                                    x = spence.int$Broad_Annotation)
spence.int$Broad_Annotation <- gsub(pattern = "5",
                                    replacement = "Endothelium",
                                    x = spence.int$Broad_Annotation)
spence.int$Broad_Annotation <- gsub(pattern = "3",
                                    replacement = "Distal_Nephron",
                                    x = spence.int$Broad_Annotation)
spence.int$Broad_Annotation <- gsub(pattern = "4",
                                    replacement = "Proximal_Nephron",
                                    x = spence.int$Broad_Annotation)
spence.int$Broad_Annotation <- gsub(pattern = "6",
                                    replacement = "Ureteric Epithelium",
                                    x = spence.int$Broad_Annotation)
spence.int$Broad_Annotation <- gsub(pattern = "7",
                                    replacement = "Immune",
                                    x = spence.int$Broad_Annotation)
spence.int$Broad_Annotation <- gsub(pattern = "8",
                                    replacement = "Immune",
                                    x = spence.int$Broad_Annotation)
spence.int <- SetIdent(spence.int, value = "all_broad")


#--------------------Sub Clustering------------------

others <- spence.int[, spence.int$all_broad %in% c(5,7,8)]
nephron <- spence.int[, spence.int$all_broad %in% c(1,2,3,4)]
ue <- spence.int[, spence.int$all_broad %in% c(6)]
stroma <- spence.int[, spence.int$all_broad %in% c(0)]

## Ureteric Epithelium

ue.list <- list(ue[, ue$Sample == "SPENCE_male"],
                ue[, ue$Sample == "SPENCE_female"])
options(future.globals.maxSize = 10000 * 1024^2)

for (i in 1:length(ue.list)) {
  ue.list[[i]] <- SCTransform(ue.list[[i]], verbose = F, vars.to.regress = c("S.Score", "G2M.Score"))
}

features <- SelectIntegrationFeatures(ue.list)
ue.list <- PrepSCTIntegration(ue.list, anchor.features = features, verbose = F)
anchors <- FindIntegrationAnchors(ue.list,
                                  normalization.method = "SCT",
                                  anchor.features = features, verbose = F)
ue <- IntegrateData(anchorset = anchors,
                    normalization.method = "SCT", verbose = T)

ue <- RunPCA(ue, verbose=F)
ue <- FindNeighbors(ue, dims = 1:30, verbose = F)
ue <- FindClusters(ue, resolution = seq(0,2,0.2))
# 3d plot
ue <- RunUMAP(ue, dims = 1:30, seed.use = 250395, n.components = 3, verbose = F)
ue@reductions$umap3d <- ue@reductions$umap

#2d plot
ue <- RunUMAP(ue, dims = 1:30, seed.use = 250395, n.components = 2, verbose = F)

ue$highres <- ue$integrated_snn_res.1.8
ue <- SetIdent(ue, value = "highres")

ue$Highres_Ann <- as.character(ue$highres)

ue$Highres_Ann <- gsub(pattern = "5",
                       replacement = "U.Tip",
                       x = ue$Highres_Ann)
ue$Highres_Ann <- gsub(pattern = "0",
                       replacement = "U.Cortical",
                       x = ue$Highres_Ann)
ue$Highres_Ann <- gsub(pattern = "1",
                       replacement = "U.Cortical",
                       x = ue$Highres_Ann)
ue$Highres_Ann <- gsub(pattern = "2",
                       replacement = "U.Cortical",
                       x = ue$Highres_Ann)
ue$Highres_Ann <- gsub(pattern = "6",
                       replacement = "U.Cortical",
                       x = ue$Highres_Ann)
ue$Highres_Ann <- gsub(pattern = "7",
                       replacement = "U.Cortical",
                       x = ue$Highres_Ann)
ue$Highres_Ann <- gsub(pattern = "4",
                       replacement = "U.Med_Outer",
                       x = ue$Highres_Ann)
ue$Highres_Ann <- gsub(pattern = "8",
                       replacement = "U.Med_Outer",
                       x = ue$Highres_Ann)
ue$Highres_Ann <- gsub(pattern = "3",
                       replacement = "U.Med_Inner",
                       x = ue$Highres_Ann)

## Nephron

nephron@meta.data <- nephron@meta.data[,c(1:12, 25:27)]
neph <- nephron
neph.list <- list(neph[, neph$Sample == "SPENCE_male"],
                  neph[, neph$Sample == "SPENCE_female"])
for (i in 1:length(neph.list)) {
  neph.list[[i]] <- SCTransform(neph.list[[i]], verbose = F, vars.to.regress = c("S.Score", "G2M.Score"))
}

features <- SelectIntegrationFeatures(neph.list)
neph.list <- PrepSCTIntegration(neph.list, anchor.features = features, verbose = F)
anchors <- FindIntegrationAnchors(neph.list,
                                  normalization.method = "SCT",
                                  anchor.features = features, verbose = F)
neph <- IntegrateData(anchorset = anchors,
                      normalization.method = "SCT", verbose = T)
neph <- RunPCA(neph, verbose = F)
neph <- FindNeighbors(neph, dims = 1:30)
neph <- FindClusters(neph, resolution = seq(0,2,0.1))

#3d plot
neph <- RunUMAP(neph, dims = 1:30, seed.use = 250395, n.components = 3, verbose = F)
neph@reductions$umap3d <- neph@reductions$umap
#2d plot
neph <- RunUMAP(neph, dims = 1:30, seed.use = 250395, n.components = 2, verbose = F)

nephron <- neph
DefaultAssay(nephron) <- "RNA"
nephron <- NormalizeData(nephron)
nephron$highres <- nephron$integrated_snn_res.1

dct <- nephron[, nephron$highres == 14]
DefaultAssay(dct) <- "RNA"
dct <- NormalizeData(dct)
DefaultAssay(dct) <- "integrated"
dct <- FindVariableFeatures(dct)
dct <- RunPCA(dct)
dct <- FindNeighbors(dct)
dct <- FindClusters(dct, resolution = 0.2)

nephron$Highres_Ann <- as.character(nephron$integrated_snn_res.1)
nephron$Highres_Ann <- gsub(pattern = "22",
                            replacement = "N.NPC_Primed",
                            x = nephron$Highres_Ann)
nephron$Highres_Ann <- gsub(pattern = "21",
                            replacement = "N.PT_Mat",
                            x = nephron$Highres_Ann)
nephron$Highres_Ann <- gsub(pattern = "20",
                            replacement = "N.LoH",
                            x = nephron$Highres_Ann)
nephron$Highres_Ann <- gsub(pattern = "19",
                            replacement = "N.NP_Str",
                            x = nephron$Highres_Ann)
nephron$Highres_Ann <- gsub(pattern = "18",
                            replacement = "N.PEC",
                            x = nephron$Highres_Ann)
nephron$Highres_Ann <- gsub(pattern = "17",
                            replacement = "N.NPC_CC",
                            x = nephron$Highres_Ann)
nephron$Highres_Ann <- gsub(pattern = "16",
                            replacement = "N.Distal_EN",
                            x = nephron$Highres_Ann)
nephron$Highres_Ann <- gsub(pattern = "15",
                            replacement = "N.DST",
                            x = nephron$Highres_Ann)
nephron$Highres_Ann <- gsub(pattern = "14",
                            replacement = "N.DCT_CS",
                            x = nephron$Highres_Ann)
nephron$Highres_Ann <- gsub(pattern = "13",
                            replacement = "N.PEC",
                            x = nephron$Highres_Ann)
nephron$Highres_Ann <- gsub(pattern = "12",
                            replacement = "N.NPC_PTA_CC",
                            x = nephron$Highres_Ann)
nephron$Highres_Ann <- gsub(pattern = "11",
                            replacement = "N.Distal_EN",
                            x = nephron$Highres_Ann)
nephron$Highres_Ann <- gsub(pattern = "10",
                            replacement = "N.LoH",
                            x = nephron$Highres_Ann)
nephron$Highres_Ann <- gsub(pattern = "9",
                            replacement = "N.PT_Mat",
                            x = nephron$Highres_Ann)
nephron$Highres_Ann <- gsub(pattern = "8",
                            replacement = "N.RV",
                            x = nephron$Highres_Ann)
nephron$Highres_Ann <- gsub(pattern = "7",
                            replacement = "N.Pod_Dev",
                            x = nephron$Highres_Ann)
nephron$Highres_Ann <- gsub(pattern = "6",
                            replacement = "N.Medial_EN",
                            x = nephron$Highres_Ann)
nephron$Highres_Ann <- gsub(pattern = "5",
                            replacement = "N.PT_Dev",
                            x = nephron$Highres_Ann)
nephron$Highres_Ann <- gsub(pattern = "4",
                            replacement = "N.NPC",
                            x = nephron$Highres_Ann)
nephron$Highres_Ann <- gsub(pattern = "3",
                            replacement = "N.NPC",
                            x = nephron$Highres_Ann)
nephron$Highres_Ann <- gsub(pattern = "2",
                            replacement = "N.NPC_Primed",
                            x = nephron$Highres_Ann)
nephron$Highres_Ann <- gsub(pattern = "1",
                            replacement = "N.Pod_Mat",
                            x = nephron$Highres_Ann)
nephron$Highres_Ann <- gsub(pattern = "0",
                            replacement = "N.NPC_PTA",
                            x = nephron$Highres_Ann)

dct$Highres_Ann <- as.character(dct$integrated_snn_res.0.2)

dct$Highres_Ann <- gsub(pattern = "0",
                        replacement = "N.CS",
                        x = dct$Highres_Ann)
dct$Highres_Ann <- gsub(pattern = "1",
                        replacement = "N.DCT",
                        x = dct$Highres_Ann)
nephron.data <- data.frame(cell = colnames(nephron), stringsAsFactors = F)
ann.data <- data.frame(cell = c(colnames(nephron[, nephron$Highres_Ann != "N.DCT_CS"]), colnames(dct)),
                       ann = c(nephron[, nephron$Highres_Ann != "N.DCT_CS"]$Highres_Ann, dct$Highres_Ann), stringsAsFactors = F)


ann.data <- ann.data %>% arrange(cell)
nephron.data <- left_join(nephron.data, ann.data, by = "cell")
head(nephron.data, 20)


nephron$Highres_Ann <- nephron.data$ann

## Stroma

stroma@meta.data <- stroma@meta.data[,c(1:12, 25:27)]
stroma.list <- list(stroma[, stroma$Sample == "SPENCE_male"],
                    stroma[, stroma$Sample == "SPENCE_female"])
options(future.globals.maxSize = 10000 * 1024^2)

for (i in 1:length(stroma.list)) {
  stroma.list[[i]] <- SCTransform(stroma.list[[i]], verbose = F, vars.to.regress = c("S.Score", "G2M.Score"))
}

features <- SelectIntegrationFeatures(stroma.list)
stroma.list <- PrepSCTIntegration(stroma.list, anchor.features = features, verbose = F)
anchors <- FindIntegrationAnchors(stroma.list,
                                  normalization.method = "SCT",
                                  anchor.features = features, verbose = F)
stroma <- IntegrateData(anchorset = anchors,
                        normalization.method = "SCT", verbose = T)
stroma <- RunPCA(stroma, verbose = F)
stroma <- FindNeighbors(stroma, dims = 1:30)
stroma <- FindClusters(stroma, resolution = seq(0,1,0.1))
#3d
stroma <- RunUMAP(stroma, dims = 1:30, seed.use = 250395, n.components = 3, verbose = F)
stroma@reductions$umap3d <- stroma@reductions$umap
#2d
stroma <- RunUMAP(stroma, dims = 1:30, seed.use = 250395, n.components = 2, verbose = F)

DefaultAssay(stroma) <- "RNA"
stroma <- NormalizeData(stroma, verbose = F)
stroma$highres <- stroma$integrated_snn_res.0.2
stroma <- SetIdent(stroma, value = "highres")

stroma$Highres_Ann <- as.character(stroma$highres)

stroma$Highres_Ann <- gsub(pattern = "1",
                           replacement = "S.IC",
                           x = stroma$Highres_Ann) # ALDH1A2 weaker MEIS2 PBX1
stroma$Highres_Ann <- gsub(pattern = "0",
                           replacement = "S.Med",
                           x = stroma$Highres_Ann) # ALX1 CLDN11 DCN ...
stroma$Highres_Ann <- gsub(pattern = "2",
                           replacement = "S.OC_NZ",
                           x = stroma$Highres_Ann) # NTN1 FOXD1 MEIS2 PBX1
stroma$Highres_Ann <- gsub(pattern = "3",
                           replacement = "S.Mesangial",
                           x = stroma$Highres_Ann)
stroma$Highres_Ann <- gsub(pattern = "4",
                           replacement = "S.NP_Str",
                           x = stroma$Highres_Ann)

hfk.data <- data.frame(cell = (colnames(spence.int)), stringsAsFactors = F)
ann.data <- data.frame(cell = (c(colnames(ue), colnames(nephron),
                                 colnames(others), colnames(stroma)
)),
ann = (c(ue$Highres_Ann,
         nephron$Highres_Ann,
         others$Broad_Annotation,
         stroma$Highres_Ann)), stringsAsFactors = F)
ann.data <- ann.data %>% arrange(cell)
hfk.data <- left_join(hfk.data, ann.data, by = "cell")
spence.int$Highres_Ann <- hfk.data$ann

spence.int$Identity <- spence.int$Highres_Ann
spence.int$Identity <- factor(spence.int$Identity,
                              levels = rev(c("S.OC_NZ", "S.IC", "S.Med", "S.Mesangial", "S.NP_Str","N.NP_Str",
                                             "N.NPC", "N.NPC_CC", "N.NPC_Primed", "N.NPC_PTA", "N.NPC_PTA_CC",
                                             "N.RV", "N.PEC", "N.Pod_Dev", "N.Pod_Mat",
                                             "N.Medial_EN", "N.PT_Dev", "N.PT_Mat", "N.Distal_EN",
                                             "N.LoH", "N.DST", "N.DCT", "N.CS",
                                             "U.Tip", "U.Cortical", "U.Med_Outer", "U.Med_Inner",
                                             "Endothelium", "Immune")))

spence.int$IdentityRandom <- factor(spence.int$Highres_Ann,
                                    levels = c(as.character(rev(unique(spence.int$Highres_Ann)[1:12])),
                                               as.character(unique(spence.int$Highres_Ann)[13:29])))

save(spence.int, file = here::here("data/Holloway2020.rda"))




