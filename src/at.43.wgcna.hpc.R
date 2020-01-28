require(plyr)
require(ape)
require(ggplot2)
require(tidyr)
require(WGCNA)

options(stringsAsFactors = FALSE)

dirw = '/home/springer/zhoux379/scratch/briggs2'
diro = file.path(dirw, '52.wgcna')

fi = file.path(dirw, "36.long.filtered.tsv")
ti = read.table(fi, sep = "\t", header = T, as.is = T)

enableWGCNAThreads()
allowWGCNAThreads()

gts = c("B73", "Mo17", "B73xMo17")
softPowers = c(B73=8, Mo17=8, B73xMo17=10)

for (gt in gts) {
tiw = spread(ti[ti$Genotype == gt, -c(3,4)], Tissue, fpkm)

datExpr = t(as.matrix(tiw[,-1]))
colnames(datExpr) = tiw[,1]

### hclust and cut tree into modules
softPower = as.numeric(softPowers[gt])
adjacency = adjacency(datExpr, power = softPower, type = "unsigned", corFnc = "cor")
dim(adjacency)

TOM = TOMsimilarity(adjacency, TOMType = "unsigned")
dissTOM = 1-TOM
geneTree = hclust(as.dist(dissTOM), method = "average")

#as.dist(1-cor(x), method = cor_opt))

# Plot the resulting clustering tree (dendrogram)
fp1 = sprintf("%s/03.tree.%s.pdf", diro, gt)
pdf(fp1, width = 12, height = 9)
#sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
    labels = FALSE, hang = 0.04)
dev.off()

# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
              deepSplit = 2, pamRespectsDendro = FALSE,
              minClusterSize = minModuleSize);
table(dynamicMods)

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath

fp2 = sprintf("%s/04.tree.%s.pdf", diro, gt)
pdf(fp2, width = 12, height = 9)
#sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                  dendroLabels = FALSE, hang = 0.03,
                  addGuide = TRUE, guideHang = 0.05,
                  main = "Gene dendrogram and module colors")
dev.off()
wg_clus1 = dynamicMods; wg_cols1 = dynamicColors

# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
fp3 = sprintf("%s/05.modules.%s.pdf", diro, gt)
pdf(fp3, width = 7, height = 6)
#sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
xlab = "", sub = "")

MEDissThres = 0.25
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;
dev.off()

#sizeGrWindow(12, 9)
fp4 = sprintf("%s/06.%s.pdf", diro, gt)
pdf(fp4, width = 12, height = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
c("Dynamic Tree Cut", "Merged dynamic"),
dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()

# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs
wg_clus2 = moduleLabels; wg_cols2 = mergedColors
wg_tree = geneTree

fo = sprintf("%s/09.%s.rda", diro, gt)
tm = data.frame(gid = tiw[,1], clu1 = wg_clus1, col1 = wg_cols1, clu2 = wg_clus2, col2 = wg_cols2, stringsAsFactors = F)

save(TOM, wg_tree, tm, file = fo)
}