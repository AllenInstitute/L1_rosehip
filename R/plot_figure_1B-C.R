library(pheatmap)
library(RColorBrewer)
library(feather)
library(igraph)


#### Figure 1B - PCA ####
pca1 <- read.csv("data/pca_sample_coords.csv", row.names = 1)
clid <- read.csv("data/clid.csv", row.names = 1)
pca.cl <- merge(pca1, clid, by = "row.names")

pdf(file = "output/Figure1B_pca.pdf", width = 5, height = 5, onefile = FALSE)
plot(PC2_rot ~ PC1_rot, col = x, data = pca.cl)
dev.off()



#### Figure 1B - heatmap ####
cl.cons.corr <- read.csv("data/cl.cons.csv.gz", row.names = 1)

hm.colors <- colorRampPalette(c("white", brewer.pal(9, "YlOrRd")))(100)

pdf(file = "output/Figure1B_cc_heatmap.pdf", width = 10, height = 10, 
    onefile = FALSE)
pheatmap(cl.cons.corr, clustering_method = "ward.D2", 
         clustering_distance_rows = as.dist(1 - cl.cons.corr),
         clustering_distance_cols = as.dist(1 - cl.cons.corr),
         color = hm.colors, border_color = NA, 
         fontsize_row = 2, fontsize_col = 2)
dev.off()



#### Figure 1C - constellation ####
cl.anno <- as.data.frame(read_feather("data/human/anno.feather"))
cl.anno$cluster_label <- sapply(cl.anno$cluster_label, 
                                    function(x) strsplit(x, "_")[[1]][1])
cl.type <- read.csv("data/cl.anno.csv", row.names = 1)
cl.cons.mean <- apply(cl.cons.corr, 1, function(x) tapply(x, cl.anno$cluster_label, mean))
cl.cons.mean <- apply(cl.cons.mean, 1, function(x) tapply(x, cl.anno$cluster_label, mean))
keep.cl <- grep("[ieg]", row.names(cl.cons.mean))
cocl <- cl.cons.mean[keep.cl, keep.cl]

# Plot network
ig <- graph.adjacency(cocl, mode="undirected", weighted=TRUE, diag=FALSE)
E(ig)$width <- 30*E(ig)$weight
V(ig)$size <- 3 * sqrt(table(cl.anno$cluster_label)[keep.cl])
V(ig)$label.color <- "black"
V(ig)$color <- unique(cl.anno$cluster_color[match(row.names(cocl), cl.anno$cluster_label)])
coords <- as.matrix(read.csv("data/cl_constellation_coords.csv", header = FALSE))
E(ig)$color <- ifelse(E(ig)$weight >= 0.05, "grey40", NA)

pdf("output/cl_constellation.pdf", width = 12, height = 12)
plot.igraph(ig, layout = coords)
legend("topright", col = "grey40", cex = 0.6,
       lwd = 30*c(0.05, 0.1, 0.2), 
       legend = c(0.05, 0.1, 0.2))
dev.off()



#### Figure 1C - pie chart ####
cl.anno.subset <- droplevels(subset(cl.anno, cluster_type_label %in% c("inh"),
                                 select = c("cluster_label", "cluster_color")))

pie.color <- cl.anno.subset$cluster_color
names(pie.color) <- cl.anno.subset$cluster_label
cl.prop <- sort(table(cl.anno.subset$cluster_label))
pdf("output/inh_prop_pie.pdf", width = 4, height = 4)
pie(cl.prop, col = as.character(pie.color[names(cl.prop)]))
dev.off()
