source("functions.R")
dirw = '/home/springer/zhoux379/scratch/briggs2/52.wgcna'
fl = file.path(dirw, '00.1.read.correct.tsv')
tl = read.table(fl, header = T, sep = "\t", as.is = T)[,1:5]
tissues = unique(tl$Tissue)
fi = file.path(dirw, "36.long.filtered.tsv")
ti = read.table(fi, sep = "\t", header = T, as.is = T)
enableWGCNAThreads()
gts = c("B73", "Mo17", "B73xMo17")

#{{{ soft threshold
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
}}}

##### run br.34.wgcna.hpc.R using soft threshold on MSI

#{{{ wgcna plot
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
}}}

#{{{ camoco clusters
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
}}}



#{{{ compare B and M
fx = file.path(dirw, "09.B73.rda")
x = load(fx)
x

tom1 = TOM
tm1 = tm
wg_tree1 = wg_tree

fx = file.path(dirw, "09.Mo17.rda")
x = load(fx)
x

tom2 = TOM
tm2 = tm
wg_tree2 = wg_tree


## module eigengenes
gt = gts[1]
tiw = spread(ti[ti$Genotype == gt, -c(3,4)], Tissue, fpkm)
datExpr = t(as.matrix(tiw[,-1]))
colnames(datExpr) = tiw[,1]
datExprb = datExpr

gt = gts[2]
tiw = spread(ti[ti$Genotype == gt, -c(3,4)], Tissue, fpkm)
datExpr = t(as.matrix(tiw[,-1]))
colnames(datExpr) = tiw[,1]
datExprm = datExpr

gt = gts[3]
tiw = spread(ti[ti$Genotype == gt, -c(3,4)], Tissue, fpkm)
datExpr = t(as.matrix(tiw[,-1]))
colnames(datExpr) = tiw[,1]
datExprh = datExpr

MEList = moduleEigengenes(datExprb, colors = tm1$col1)
MEsb = cbind(tissue = rownames(datExprb), MEList$eigengenes)
teb = gather(MEsb, module, ME, -tissue)
tmp = 1:length(tissues)
names(tmp) = tissues
teb = cbind(x = tmp[teb$tissue], teb)

pb <- ggplot(teb) +
  geom_tile(aes(x = x, y = module, fill = ME), height = 1) +
  theme_bw() + 
  scale_x_continuous(name = 'Tissue', breaks = teb$x, labels = teb$tissue, expand = c(0,0)) +
  scale_y_discrete(name = 'B73 Module', expand = c(0,0)) +
  scale_fill_gradient(name = '', space = "Lab", low = "white", high = "brown1") +
#  theme(legend.position = 'none') +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines")) +
  theme(axis.title.x = element_text(colour = "black", size = 9), axis.ticks.length = unit(0, 'lines')) +
  theme(axis.text.x = element_text(colour = "black", size = 9, angle = 30, hjust = 1)) +
  theme(axis.title.y = element_text(colour = "black", size = 9), axis.ticks.length = unit(0, 'lines')) +
  theme(axis.text.y = element_text(colour = "black", size = 8)) +
  theme(axis.line = element_line(size = 0.3, colour = "grey", linetype = "solid"))
fo = sprintf("%s/24.modules.eigengene.pdf", dirw)
ggsave(pb, filename = fo, width = 7, height = 8)


MEList2 = moduleEigengenes(datExprm, colors = tm1$col1)
MEsb = cbind(tissue = rownames(datExprm), MEList2$eigengenes)
teb = gather(MEsb, module, ME, -tissue)
tmp = 1:length(tissues)
names(tmp) = tissues
teb = cbind(x = tmp[teb$tissue], teb)

pb <- ggplot(teb) +
  geom_tile(aes(x = x, y = module, fill = ME), height = 1) +
  theme_bw() + 
  scale_x_continuous(name = 'Tissue', breaks = teb$x, labels = teb$tissue, expand = c(0,0)) +
  scale_y_discrete(name = 'B73 Module', expand = c(0,0)) +
  scale_fill_gradient(name = '', space = "Lab", low = "white", high = "brown1") +
#  theme(legend.position = 'none') +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines")) +
  theme(axis.title.x = element_text(colour = "black", size = 9), axis.ticks.length = unit(0, 'lines')) +
  theme(axis.text.x = element_text(colour = "black", size = 9, angle = 30, hjust = 1)) +
  theme(axis.title.y = element_text(colour = "black", size = 9), axis.ticks.length = unit(0, 'lines')) +
  theme(axis.text.y = element_text(colour = "black", size = 8)) +
  theme(axis.line = element_line(size = 0.3, colour = "grey", linetype = "solid"))
fo = sprintf("%s/24.modules.eigengene.m.pdf", dirw)
ggsave(pb, filename = fo, width = 7, height = 8)

cor.test(as.matrix(MEList$eigengenes), as.matrix(MEList2$eigengenes))

}}}

#{{{ module preservation btw B and M - cross tabulation
fp4 = sprintf("%s/22.modules_b+m.pdf", dirw)
pdf(fp4, width = 12, height = 6)
plotDendroAndColors(wg_tree1, cbind(tm1$col1, tm2$col1), c("B73", "Mo17"),
dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05,
main = 'B73 gene dendrogram and module colors', ylab = NULL)
dev.off()

cols1 = unique(tm1$col1)
cols2 = unique(tm2$col1)
tc = data.frame(x = rep(1:length(cols1), each = length(cols2)),
	y = rep(1:length(cols2), length(cols1)),
	B73 = rep(cols1, each = length(cols2)), Mo17 = rep(cols2, length(cols1)), 
	nb = 0, nm = 0, no = 0, pval = 0, stringsAsFactors = F)

for (i in 1:nrow(tc)) {
	col1 = tc$B73[i]; col2 = tc$Mo17[i]
	tc$nb[i] = sum(tm1$col1 == col1)
	tc$nm[i] = sum(tm2$col1 == col2)
	tc$no[i] = sum(tm1$col1 == col1 & tm2$col1 == col2)
	nfn = nrow(tm1) - tc$nb[i] - tc$nm[i] - tc$no[i]
	ftest = fisher.test(matrix(c(tc$no[i], tc$nb[i], tc$nm[i], nfn), nrow = 2))
	tc$pval[i] = ftest$p.value
}

tc$pval[tc$pval==0] = min(tc$pval[tc$pval > 0])
tc = cbind(tc, logp = -log(tc$pval), label = sprintf("%d\n%s", tc$no, prettyNum(tc$pval, digits = 1)))

tcx = unique(tc[,c('x','B73','nb')])
tcx = cbind(tcx, label=sprintf("%s\n%d", tcx$B73, tcx$nb))
tcy = unique(tc[,c('y','Mo17','nm')])
tcy = cbind(tcy, label=sprintf("%s\n%d", tcy$Mo17, tcy$nm))

pb <- ggplot(tc) +
  geom_tile(aes(x = x, y = y, fill = logp), height = 1) +
  geom_text(aes(x = x, y = y, label = label), size = 2.5) +
  theme_bw() + 
  scale_x_continuous(name = 'B73', breaks = tcx$x, labels = tcx$label, expand = c(0,0)) +
  scale_y_continuous(name = 'Mo17', breaks = tcy$y, labels = tcy$label, expand = c(0,0)) +
  scale_fill_gradient(name = '', space = "Lab", low = "white", high = "brown1", na.value = 'grey50') +
  theme(legend.position = 'none') +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines")) +
  theme(axis.title.x = element_text(colour = "black", size = 9), axis.ticks.length = unit(0, 'lines')) +
  theme(axis.text.x = element_text(colour = "black", size = 8, angle = 45, hjust = 1)) +
  theme(axis.title.y = element_text(colour = "black", size = 9), axis.ticks.length = unit(0, 'lines')) +
  theme(axis.text.y = element_text(colour = "black", size = 8)) +
  theme(axis.line = element_line(size = 0.3, colour = "grey", linetype = "solid"))
fo = sprintf("%s/23.modules.crosstab.pdf", dirw)
ggsave(pb, filename = fo, width = 12, height = 12)

}}} 

#{{{ module preservation - B & M
nSets = 2

multiExpr = list()
multiExpr[[1]] = list(data = datExprb)
multiExpr[[2]] = list(data = datExprm)
setLabels = c("B73", "Mo17")
names(multiExpr) = setLabels
lapply(multiExpr, lapply, dim)

colorList = list(tm1$col1, tm2$col1)
names(colorList) = setLabels

mp = modulePreservation(multiExpr, colorList, networkType = 'unsigned', corFnc = "cor",
	referenceNetworks = 1, 
	loadPermutedStatistics = FALSE, nPermutations = 200, verbose = 3)
fo = sprintf("%s/30.mp.rda", dirw)
save(mp, file = fo)

## degree correlation
sum(tom1[upper.tri(tom1)] >= 0.75)
sum(tom2[upper.tri(tom2)] >= 0.75)
sum(tom1[upper.tri(tom1)] >= 0.75 & tom2[upper.tri(tom2)] >= 0.75)

degs1 = sapply(1:nrow(tom1), myfunc <- function(i) sum(tom1[i,-i] >= 0.75))
degs2 = sapply(1:nrow(tom2), myfunc <- function(i) sum(tom2[i,-i] >= 0.75))
cor.test(degs1, degs2)

idxs1 = which(degs1 >= quantile(degs1[degs1>0], 0.9))
idxs2 = which(degs2 >= quantile(degs2[degs2>0], 0.9))
sum(idxs1 %in% idxs2)
sum(idxs1 %in% idxs2) / length(unique(c(idxs1, idxs2)))

ptitle = sprintf("PCC = %.03f", cor.test(degs1, degs2)$estimate)
tt = data.frame(b = degs1, m = degs2)
#tt = tt[tt$b > 0 & tt$m > 0,]
p2 = ggplot(tt, aes(x = b, y = m)) +
  geom_point(size = 0.5) +
  geom_smooth(method="lm") +
  scale_x_continuous(name = 'Number of Edges Per Node (B73)', limits = c(0, 1500)) +
  scale_y_continuous(name = 'Number of Edges Per Node (Mo17)', limits = c(0, 1500)) +
  theme_bw() +
  ggtitle(ptitle) + 
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(plot.margin = unit(c(0.3,0.1,0.1,0.1), "lines")) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.title.y = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "black", angle = 0)) + 
  theme(plot.title = element_text(hjust = 0.5))
fp = sprintf("%s/21.degree.corr.pdf", dirw)
ggsave(p2, filename = fp, width = 6, height = 6)

}}}

#{{{ heatmap for selected modules

gids = tm1$gid[degs1 < 1500 & degs1 > 1400]
gids = tm1$gid[tm1$col1 == 'greenyellow']
fo = sprintf("%s/z.pdf", diro)
fpkm_heatmap(gids, ti, tissues, gts, fo)

fo = sprintf("%s/module.gids.txt", dirw)
write(tm$gid[idxs1], fo)

find_enrichment.py --obo $data/db/go/go-basic.obo --outfile module.enrich.tsv module.gids.txt $genome/Zmays_v4/61.interpro/gids.txt $genome/Zmays_v4/61.interpro/11_gene.tsv

j = c("Zm00001d033859", "Zm00001d018742", "Zm00001d024523")


cols = c('blue', 'midnightblue', 'greenyellow', 'red', 'purple')
cols = c('darkmagenta', 'cyan')
for (col1 in cols) {
gids = tm1$gid[tm1$col1 == col1]
fo = sprintf("%s/25.module.heatmap.%s.pdf", dirw, col1)
fpkm_heatmap(gids, ti, tissues, gts, fo)
}

}}} 


#{{{ module preservation B, M and F1
nSets = 3

multiExpr = list()
multiExpr[[1]] = list(data = datExprb)
multiExpr[[2]] = list(data = datExprm)
multiExpr[[3]] = list(data = datExprh)
setLabels = gts
names(multiExpr) = setLabels
lapply(multiExpr, lapply, dim)

colorList = list(tm1$col1, tm2$col1, tm1$col1)
names(colorList) = setLabels

mp = modulePreservation(multiExpr, colorList, networkType = 'unsigned', corFnc = "cor",
	referenceNetworks = 1, parallelCalculation = TRUE, 
	loadPermutedStatistics = FALSE, nPermutations = 200, verbose = 3)
fo = sprintf("%s/30.mp.rda", dirw)
save(mp, file = fo)
}}} 

#{{{ process outputs
fo = sprintf("%s/30.mp.rda", dirw)
x = load(fo)
x

ref = 1
test = 2
ref = 1
test = 2
modNames = rownames(mp$preservation$observed[[ref]][[test]])
modSizes = mp$preservation$Z[[ref]][[test]][, 1];

stats = cbind(modName = modNames, moduleSize = modSizes, 
	mp$quality$observed[[ref]][[test]][,-1], 
	mp$preservation$observed[[ref]][[test]][, -1],
	mp$quality$Z[[ref]][[test]][,-1], 
	mp$preservation$Z[[ref]][[test]][,-1]
)

t1 = stats[!stats$modName %in% c('gold','grey'), c('modName','moduleSize','medianRank.pres','Zsummary.pres')]

p1 = ggplot(t1, aes(x = medianRank.pres, y = Zsummary.pres, size = moduleSize, label = modName)) +
  geom_point(aes(color = modName)) +
  geom_text(vjust = 0, nudge_y = 1) +
  #scale_x_continuous(name = '') +
  #scale_y_continuous(name = '') +
  scale_color_manual(values = t1$modName, guide = F) +
  geom_hline(yintercept = 2, color = 'blue', linetype = 2) + 
  geom_hline(yintercept = 10, color = 'darkgreen', linetype = 2) + 
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "lines")) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.title.y = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "black", angle = 0))
fp = sprintf("%s/31.mp.pdf", dirw)
ggsave(p1, filename = fp, width = 7, height = 6)

}}} 

#{{{ circle plot
pid = 'orangered4'
idxGenes = which(tm1$col1 == pid)
nPathGenes = length(idxGenes)
pathwayAdjs = list()
for (set in 1:nSets)
{
  bc = bicor(multiExpr[[set]]$data[, idxGenes], use = "p")
  pathwayAdjs[[set]] = abs(bc)^4 * sign(bc)
}

conn = matrix(0, nPathGenes, nSets);
for (set in 1:nSets)
  conn[, set] = apply(abs(pathwayAdjs[[set]]), 2, sum)-1
weights = c(3,1,5,1, 3,1,5,1);
wMat = matrix(weights, nPathGenes, nSets, byrow = TRUE);
wconn = apply(conn * wMat, 1, sum);
order = order(-wconn);
# use the gene names as lables
labels = colnames(pathwayAdjs[[ref]]);
labels = substr(labels, 7, 14)

sizeGrWindow(8,4)
fo = file.path(dirw, '35.circleplot.pdf')
pdf(file = fo, wi=12, h=4)
par(mfrow =c(1,3));
par(mar = c(0.3, 0.2, 1.5, 0.2))
for (set in 1:nSets)
{
  circlePlot(pathwayAdjs[[set]], labels, order, main = setLabels[set],
            variable.cex.labels = TRUE,
            radii = c(0.56, 0.62), center = c(0.1, 0.04),
            min.cex.labels = 1.2, max.cex.labels = 1.4, cex.main = 1.4);
}
dev.off();
}}} 

#{{{ clusterRepro
multiMEs = multiSetMEs(multiExpr, universalColors = tm1$col1)

cr = list()
set.seed(10)
fo = file.path(dirw, "30.cr.rda")
for (set in 1:nSets)
{
  printFlush(paste("Working on ", setLabels[set]))
  centr = as.matrix(multiMEs[[set]]$data[, c(2:3)])
  rownames(centr) = rownames(multiExpr[[set]]$data)
  colnames(centr) = c("Random", "CholSynth")
  print(system.time({
    cr[[set]] = clusterRepro(Centroids = centr, New.data = multiExpr[[set]]$data,
                    Number.of.permutations = 200);
     } ))
  save(cr, file = fo)
}
}}}


#{{{ rep
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
}}}



