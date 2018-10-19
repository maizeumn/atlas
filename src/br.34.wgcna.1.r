source("br.fun.R")

dirw = '/home/springer/zhoux379/scratch/briggs2/52.wgcna'

fi = file.path(dirw, "36.long.filtered.tsv")
ti = read.table(fi, sep = "\t", header = T, as.is = T)

for (gt in gts) {

tiw = spread(ti[ti$Genotype == gt, -c(3,4)], Tissue, fpkm)

datExpr = t(as.matrix(tiw[,-1]))
colnames(datExpr) = tiw[,1]

### choose power
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
powers = 1:20
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:
fo = sprintf("%s/01.power.%s.pdf", dirw, gt)
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

##### run br.34.wgcna.hpc.R using soft threshold on MSI

#####
fx = file.path(dirw, "09.B73.rda")
x = load(fx)
x

MEList = moduleEigengenes(datExpr, colors = wg_cols1, excludeGrey = T)
MEs = MEList$eigengenes
MEDiss = 1-abs(cor(MEs))
pdf(file.path(diro, "11.1.pdf"), width = 6, height = 6)
heatmap(MEDiss, Rowv=NA, Colv=NA, symm=TRUE)
#heatmap(as.matrix(dist(t(MEs),diag=TRUE)), Rowv=NA, Colv=NA, symm=TRUE, scale = "column")
dev.off()

MEList = moduleEigengenes(datExpr, colors = wg_cols2, excludeGrey = T)
MEs = MEList$eigengenes
MEDiss = 1-abs(cor(MEs))
pdf(file.path(dirw, "11.2.pdf"), width = 5, height = 5)
heatmap(MEDiss, Rowv=NA, Colv=NA, symm=TRUE)
#heatmap(as.matrix(dist(t(MEs),diag=TRUE)), Rowv=NA, Colv=NA, symm=TRUE, scale = "column")
dev.off()

### camoco clusters
fc = file.path(dirw, "51.camoco/clusters.csv")
tc = read.table(fc, sep = ",", header = T, as.is = T)
tc$cluster = tc$cluster + 1

tc2 = ddply(tc, .(cluster), summarise, ngene = length(X))
clus = tc2$cluster[tc2$ngene >= 30]

tgc = data.frame(gid = toupper(colnames(datExpr)), clu = 0, stringsAsFactors = F)
tgc2 = merge(tgc, tc, by.x = 'gid', by.y = 'X', all.x = T)
stopifnot(identical(tgc$gid, tgc2$gid))
idxs = which(!is.na(tgc2$cluster))
tgc2$clu[idxs] = tgc2$cluster[idxs]
tgc2$clu[!tgc2$clu %in% clus] = NA

identical(tgc2$gid, toupper(colnames(datExpr)))

pdf(file.path(dirw, "15.pdf"), width = 12, height = 7)
plotDendroAndColors(wg_tree, cbind(wg_cols1, wg_cols2, tgc2$clu),
c("WGCNA Raw", "WGCNA Merged", "Camoco"),
dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()

MEList = moduleEigengenes(datExpr, colors = tgc2$clu, excludeGrey = T)
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs)
METree = hclust(as.dist(MEDiss), method = "ward.D")
pdf(file.path(dirw, "31.camoco.tree.pdf"), width = 8, height = 5)
plot(METree, main = "Clustering of Camoco module eigengenes", cex = 0.7,
xlab = "", sub = "")
dev.off()

### 
MEDiss3 = 1-abs(cor(MEs))
pdf(file.path(dirw, "32.pdf"), width = 8, height = 8)
heatmap(MEDiss, Rowv=NA, Colv=NA, symm=TRUE)
dev.off()

