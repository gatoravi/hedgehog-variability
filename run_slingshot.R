library(slingshot)
library(Seurat)
library(grDevices)
library(RColorBrewer)
library(ggplot2)

load("all_samples.Robject")

all_samples_notg2 <- all_samples
all_samples_notg2 <- all_samples[, all_samples@meta.data$Phase != "G2M"]
DimPlot(all_samples_notg2, group.by = "seurat_clusters", label = T)


# Based on https://bioconductor.org/packages/release/bioc/vignettes/slingshot/inst/doc/vignette.html
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
#load(file = "sce.Robject")
sce <- as.SingleCellExperiment(all_samples_notg2)


## Run slingshot
sce <- slingshot(sce, clusterLabels = sce$seurat_clusters, reducedDim = "UMAP", extend = "y", stretch = 0)
sds <- SlingshotDataSet(sce)

cluster_colors_palette <- c('#8dd3c7','#fccde5','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69')
hedgehog_gradient_palette <- c('#f1eef6','#f1eef6','#f1eef6','#2b8cbe','#045a8d')
ident_colors_palette <- c('#66c2a5','#fc8d62','#8da0cb')
cells_clusters_colors <- cluster_colors_palette[as.numeric(all_samples_notg2$seurat_clusters)]
hedgehog_gradient_colors <- hedgehog_gradient_palette[all_samples_notg2$hedgehog + 1]
ident_colors <- ident_colors_palette[as.numeric(as.factor(all_samples_notg2$orig.ident), levels = c("untreated", "hr17", "hr30"))]
print(table(all_samples_notg2$orig.ident))
print(class(all_samples_notg2$orig.ident))

pdf("slingshot.pdf")
    DimPlot(all_samples_notg2, group.by = "orig.ident", cols = ident_colors_palette)
    DimPlot(all_samples_notg2, group.by = "seurat_clusters", cols = cluster_colors_palette)
dev.off()
#    plot(sds, show.constraints = T, type = "both")
pdf("slingshot_seuratclusters.pdf")
    plot(reducedDims(sce)$UMAP, asp = 1, col = cells_clusters_colors)
    lines(SlingshotDataSet(sce), lwd=2, col='black')
dev.off()
#    plot(reducedDims(sce)$UMAP, asp = 1, col = cells_clusters_colors )
#    lines(SlingshotDataSet(sce), type = 'lineages', lwd=2, col='black')
pdf("slingshot_hedgehogexpression.pdf")
    plot(reducedDims(sce)$UMAP, asp = 1, col = hedgehog_gradient_colors)
    lines(SlingshotDataSet(sce), lwd=2, col='black')
dev.off()
pdf("slingshot_ident.pdf")
    plot(reducedDims(sce)$UMAP, asp = 1, col = ident_colors)
    lines(SlingshotDataSet(sce), lwd=2, col='black')
dev.off()
# Comment out for now to save running time
#save(sce, file = "sce.Robject")
