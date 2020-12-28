library(tidyverse)
library(Seurat)

seurat <- read_rds("../output/seurat/HumanFetalKidney_Reference_Spence.rds")
clusters <- unique(seurat$Identity)
clustlength <- length(clusters)
value.list <- list()
seurat <- SetIdent(seurat, value = "Identity")

for (i in 1:length(clusters)) {
  for (j in 1:length(clusters)) {
    if (i != j) {
      temp <- FindMarkers(object = seurat, group.by = "Identity", ident.1 = paste0(clusters[i]),
                          random.seed = 250395, 
                          ident.2 = paste0(clusters[j]), only.pos = T, test.use = "t")
      value.list[[paste0(clusters[i], ".", clusters[j])]] <- temp
      
    }
  }
}

write_rds(value.list, "HFK_Genesets_Identity.rds")

## end of job

