# Standard dotplot function
SWDP.bw <- function(x, features, group.by = "Identity", assay = "RNA", col.min = 0, col.max = 10, split.by = NULL) {
  DotPlot(object = x, features = features, col.min = col.min, col.max = col.max,
          cols = c("lightgrey", "black"),
          group.by = group.by, dot.min = 0.1,
          assay = assay) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0.5)) +
    theme(axis.text=element_text(size=rel(0.5))) +
    theme(panel.grid.major = element_line(colour = "lightgray"),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
}

SWDP.col <- function(x, features, group.by = "Identity", assay = "RNA", col.min = 0, col.max = 10, split.by = NULL, col = "red") {
  DotPlot(object = x, features = features, col.min = col.min, col.max = col.max,
          cols = c("lightgrey", col),
          group.by = group.by, dot.min = 0.1,
          assay = assay) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0.5)) +
    theme(axis.text=element_text(size=rel(1))) +
    theme(panel.grid.major = element_line(colour = "lightgray"),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
}

SWDP.vir <- function(x, features, group.by = "Identity", assay = "RNA", col.min = 0, col.max = 10, split.by = NULL) {
  DotPlot(object = x, features = features, col.min = col.min, col.max = col.max,
          cols = rev(c("#482677FF", "#95D840FF")),
          group.by = group.by, dot.min = 0.1,
          assay = assay) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0.5)) +
    theme(axis.text=element_text(size=rel(1))) +
    theme(panel.grid.major = element_line(colour = "lightgray"),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
}

# export DE genes as excel
exportDEsheet <- function(seurat.file, seurat.ident, markers, filename) {
  export <- lapply(0:(length(unique((seurat.file@active.ident)))-1),
                   function(x) {
                     markers %>%
                       dplyr::filter(cluster == unique(neph.merge@active.ident)[x]
                                     , p_val_adj < 0.05, avg_logFC > 0) %>%
                       dplyr::arrange(-avg_logFC) %>%
                       select(Gene = gene, LogFC = avg_logFC, pVal = p_val_adj)
                   })
  WriteXLS::WriteXLS(export,
                     ExcelFileName = filename,
                     SheetNames = paste0("Cluster ", 0:(length(unique(colnames(seurat.file@metadata$seurat.ident)))-1)))
}

'%!in%' <- function(x,y)!('%in%'(x,y))

ggplotColors <- function(g){
  d <- 360/g
  h <- cumsum(c(15, rep(d,g - 1)))
  hcl(h = h, c = 100, l = 65)
}

## add genes in set to results table
addGenes <- function(results, cp_index, geneSymbols, topTab=NULL, pval=NULL){

  results = as.data.frame(results)
  results$genes = NA_character_

  for(i in 1:nrow(results)){
    symbols = geneSymbols[unlist(cp_index[rownames(results)[i]], use.names=FALSE)]
    symbols[symbols == ""] <- NA
    results$genes[i] = paste(symbols, collapse=",")
    if(!is.null(topTab)){
      tmp = topTab[, grepl("symbol",colnames(topTab),ignore.case = TRUE)]

      results$up[i] = paste(tmp[topTab$adj.P.Val < pval & topTab$logFC > 0 &
                                  tmp %in% symbols],
                            collapse = ",")
      results$down[i] = paste(tmp[topTab$adj.P.Val < pval & topTab$logFC < 0 &
                                    tmp %in% symbols],
                              collapse = ",")
    }
  }

  results
}
