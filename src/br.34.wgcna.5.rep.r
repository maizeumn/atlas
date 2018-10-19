require(plyr)
require(ape)
require(ggplot2)
require(WGCNA)
require(tidyr)
require(RColorBrewer)
require(pheatmap)

options(stringsAsFactors = FALSE)

dirw = '/home/springer/zhoux379/scratch/briggs2'
diro = file.path(dirw, '52.wgcna', 'control')

fl = file.path(dirw, '00.1.read.correct.tsv')
tl = read.table(fl, header = T, sep = "\t", as.is = T)[,1:5]
tissues = unique(tl$Tissue)

fi = file.path(dirw, "36.long.filtered.tsv")
ti = read.table(fi, sep = "\t", header = T, as.is = T)
gids = unique(ti$gid)

#allowWGCNAThreads()
enableWGCNAThreads()
gts = c("B73", "Mo17", "B73xMo17")

#####
fi = file.path(dirw, '34.fpkm.tsv')
ti = read.table(fi, header = T, sep = "\t", as.is = T)

gt = gts[1]
if (FALSE) {
for (treat in 1:3) {
#treat = 1
sids = tl$SampleID[tl$Genotype == gt & tl$Treatment == treat]

ti2 = ti[ti$gid %in% gids, c('gid',sids)]

datExpr = t(as.matrix(ti2[,-1]))
colnames(datExpr) = ti2[,1]

powers = c(c(1:10), seq(from = 12, to=20, by=2))
powers = 1:20
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
fo = sprintf("%s/01.power.%d.pdf", diro, treat)
pdf(fo, width = 9, height = 5)
#sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signed R^2",type="n", main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
    labels=powers,cex=cex1,col="red")
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
    xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
    main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
}
}


###
softPowers = c(10,9,8)
res = list()

for (treat in 2:3) {
sids = tl$SampleID[tl$Genotype == gt & tl$Treatment == treat]

ti2 = ti[ti$gid %in% gids, c('gid',sids)]

datExpr = t(as.matrix(ti2[,-1]))
colnames(datExpr) = ti2[,1]

softPower = as.numeric(softPowers[treat])
adjacency = adjacency(datExpr, power = softPower, type = "unsigned", corFnc = "cor")
dim(adjacency)

TOM = TOMsimilarity(adjacency, TOMType = "unsigned")
dissTOM = 1-TOM
geneTree = hclust(as.dist(dissTOM), method = "average")

minModuleSize = 30;
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
              deepSplit = 2, pamRespectsDendro = FALSE,
              minClusterSize = minModuleSize);
table(dynamicMods)

dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

fp2 = sprintf("%s/04.tree.%d.pdf", diro, treat)
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

MEDissThres = 0.25
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors
mergedMEs = merge$newMEs

fp4 = sprintf("%s/06.%d.pdf", diro, treat)
pdf(fp4, width = 12, height = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
c("Dynamic Tree Cut", "Merged dynamic"),
dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()

moduleColors = mergedColors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs
wg_clus2 = moduleLabels; wg_cols2 = mergedColors
wg_tree = geneTree

tm = data.frame(gid = colnames(datExpr), clu1 = wg_clus1, col1 = wg_cols1, clu2 = wg_clus2, col2 = wg_cols2, stringsAsFactors = F)
res[[treat]] = tm
}

fo = file.path(diro, "01.rda")
save(res, file = fo)