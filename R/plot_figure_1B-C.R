library(WGCNA)
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



#### Figure 1B - tree ####
cl.anno <- as.data.frame(read_feather("data/human/anno.feather"))
cl.anno$cluster_label <- sapply(cl.anno$cluster_label, 
                                function(x) strsplit(x, "_")[[1]][1])

cl.cons.corr <- read.csv("data/cl.cons.csv.gz", row.names = 1)
cl.cons.dist <- as.dist(1 - cl.cons.corr)
hc1 <- hclust(cl.cons.dist, method = "ward.D2")

pdf(file = "output/Figure1B_cc_dend.pdf", width = 12, height = 4)
pal1 <- rainbow(length(unique(cl.anno$cluster_label)))
plotDendroAndColors(hc1, colors = pal1[as.numeric(as.factor(cl.anno$cluster_label))], 
                    dendroLabels = FALSE, groupLabels = cl.anno$cluster_label)
dev.off()



#### Figure 1C - constellation ####
load("output/coreCells_clusterConfidences.RData", verbose = TRUE)
cl.type <- read.csv("data/cl.anno.csv", row.names = 1)
cocl <- cbind(clusterOverlap, e1 = 0, g1 = 0, g2 = 0, g3 = 0, g4 = 0)
cocl <- rbind(cocl, matrix(0, 5, ncol(cocl)))
short.cl.names <- sapply(colnames(cocl), function(x) strsplit(x, "_")[[1]][1])
row.names(cocl) <- short.cl.names
colnames(cocl) <- short.cl.names


# Plot network
ig <- graph.adjacency(cocl, mode="undirected", weighted=TRUE, diag=FALSE)
E(ig)$width <- E(ig)$weight
cl.size <- sqrt(table(cl.anno$cluster_label))
V(ig)$size <- 3 * cl.size[match(row.names(cocl), names(cl.size))]
V(ig)$label.color <- "black"
V(ig)$color <- unique(cl.anno$cluster_color[match(row.names(cocl), 
                                                  cl.anno$cluster_label)])

# Manual layout
coords.file <- "data/cl_constellation_coords.csv"
if (file.exists(coords.file)) {
  coords <- as.matrix(read.csv(coords.file, header = FALSE))
} else {
  id <- tkplot(ig,  layout = coords)
  coords <- tkplot.getcoords(id)
  write.table(coords, file = coords.file, 
              row.names = FALSE, col.names = FALSE, sep = ",")
}

pdf("output/cl_constellation.pdf", width = 8, height = 5)
plot.igraph(ig, layout = coords, asp = FALSE, vertex.label = NA)
# Number of intermediate nuclei
legend("topright", col = "grey40", cex = 0.6,
       lwd = c(2, 4, 8), 
       legend = c(2, 4, 8))
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
