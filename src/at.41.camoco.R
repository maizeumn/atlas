source("functions.R")
dirw = file.path(dird, "51.camoco")
fi = file.path(dirw, "../36.long.filtered.tsv")
ti = read.table(fi, sep = "\t", header = T, as.is = T)

### for how to run camoco on a workstation - see Note
#{{{ obtain camoco results from scratch
gts = c("B73", "Mo17", "B73xMo17")
for (gt in gts) {
    tiw = spread(ti[ti$Genotype == gt, -c(3,4)], Tissue, fpkm)
    expr = t(as.matrix(tiw[,-1]))
    colnames(expr) = tiw[,1]
    gids = tiw[,1]
    ng = length(gids)

    pcc.matrix = cor(asinh(expr), method = 'pearson')
    pcc = pcc.matrix[lower.tri(pcc.matrix)]

    ii = which(lower.tri(pcc.matrix))
    colidx = as.integer((ii - 1) / nrow(pcc.matrix)) + 1
    rowidx = ii - (colidx - 1) * nrow(pcc.matrix)

    pcc[pcc == 1] = 0.999999
    pcc[pcc == -1] = -0.999999
    #pcc2 = log((1+pcc) / (1-pcc)) / 2
    pcc2 = atanh(pcc)
    coexv = (pcc2 - mean(pcc2)) / sd(pcc2)

    ord = order(-coexv)
    idxs = ord[1:(1*1000000)]
    dw = data.frame(g1 = gids[rowidx[idxs]], g2 = gids[colidx[idxs]], coex = coexv[idxs])

    fo1 = sprintf("%s/01.edges.%s.tsv", dirw, gt)
    write.table(dw, fo1, sep = "\t", row.names = F, col.names = F, quote = F)

    fo2 = sprintf("%s/02.%s.mcl", dirw, gt)
    cmd = sprintf("mcl %s --abc -scheme 7 -I %g -o %s", fo1, 2, fo2)
    system(cmd)

    fo3 = sprintf("%s/03.modules.%s.tsv", dirw, gt)
    cmd = sprintf("mcl2tsv.py %s %s", fo2, fo3)
    system(cmd)

    fm = sprintf("%s/03.modules.%s.tsv", dirw, gt)
    tm = read.table(fm, header = T, as.is = T, sep = "\t")
    cat(sprintf("%s: %d modules with %d genes\n", gt, length(unique(tm$V1)), nrow(tm)))
}
#}}}

#{{{ make B,M,F1 similarity diagram
require(diagram)

fp = sprintf("%s/../61.perm/09.rda", dirw)
x = load(fp)

pme = pc[pc$perm==0, c('gt1','gt2','pcc.expr')]
pmc = pc[pc$perm==0, c('gt1','gt2','pcc.coex')]
pme$pcc.expr = sprintf("%.03f", pme$pcc.expr)
pmc$pcc.coex = sprintf("%.03f", pmc$pcc.coex)
pme = spread(pme, gt2, pcc.expr)
pmc = spread(pmc, gt2, pcc.coex)

sime = as.matrix(pme[,-1])
simc = as.matrix(pmc[,-1])
rownames(sime) = pme[,1]
rownames(simc) = pmc[,1]
sime = sime[gts,gts]
simc = simc[gts,gts]
sime[lower.tri(sime)] = NA
simc[lower.tri(simc)] = NA

fo = file.path(dirw, "81.similarity.pdf")
pdf(fo, width = 9, height = 5)
par(mfrow=c(1,2))
par(mar=c(0.1,0.1,3.1,0.1))
pp1 = plotmat(sime, pos = c(2,1), curve = 0, name = gts,
	lwd = 1, box.lwd = 2, cex.txt = 1, box.cex = 1,
	box.type = "square", box.prop = 0.6, arr.type = "triangle",
	arr.pos = 0.5, arr.width = 0, shadow.size = 0.01, prefix = "",
	main = "Expression Similarity")
pp2 = plotmat(simc, pos = c(2, 1), curve = 0, name = gts,
	lwd = 1, box.lwd = 2, cex.txt = 1, box.cex = 1,
	box.type = "square", box.prop = 0.6, arr.type = "triangle",
	arr.pos = 0.5, arr.width = 0, shadow.size = 0.01, prefix = "",
	main = "Co-expression Similarity")
dev.off()
#}}}

#{{{ construct expr matrix
gt = gts[1]
tiw = spread(ti[ti$Genotype == gt, -c(3,4)], Tissue, fpkm)
datExpr = t(as.matrix(tiw[,-1]))
colnames(datExpr) = tiw[,1]
datExprb = asinh(datExpr)
tib = tiw

gt = gts[2]
tiw = spread(ti[ti$Genotype == gt, -c(3,4)], Tissue, fpkm)
datExpr = t(as.matrix(tiw[,-1]))
colnames(datExpr) = tiw[,1]
datExprm = asinh(datExpr)
tim = tiw

gt = gts[3]
tiw = spread(ti[ti$Genotype == gt, -c(3,4)], Tissue, fpkm)
datExpr = t(as.matrix(tiw[,-1]))
colnames(datExpr) = tiw[,1]
datExprh = asinh(datExpr)
tih = tiw

# read module info
gt = gts[1]
fg = sprintf("%s/13.%s.mcl.tsv", dirw, gt)
tg = read.table(fg, header = T, sep = "\t", as.is = T)
tgs = ddply(tg, .(grp), summarise, size = length(id))
grps = tgs$grp[tgs$size >= 10]
length(grps)
tg = tg[tg$grp %in% grps,]
tg = cbind(tg, gid = tiw[tg$id,1])
nrow(tg)

mids = rep(0, nrow(tiw))
mids[tg$id] = tg$grp
#}}}

#{{{ look at overlap with corncyc pathways
tp2 = tp[tp$gid %in% tiw[,1], c('gid','p2')]
tp3 = unique(tp2)
tz = merge(tp3, tg, by = 'gid', all.x = T)
tzs = ddply(tz, .(p2), summarise, psize = length(gid), gim = sum(!is.na(grp)), ngrp = length(unique(grp[!is.na(grp)])), gnim = sum(is.na(grp)))

total_gene_in_module = nrow(tg)
total_gene_not_in_module = nrow(tib) - total_gene_in_module
pvals = apply(tzs, 1, myfunc <- function(x) dhyper(as.numeric(x[3]), total_gene_in_module - as.numeric(x[3]), total_gene_not_in_module - as.numeric(x[5]), as.numeric(x[2])))
tzs = cbind(tzs, pval = p.adjust(pvals, 'BH'))
#}}}

#{{{ WGCNA datasets
gt = gts[1]
nSets = 3

multiExpr = list()
multiExpr[[1]] = list(data = datExprb[,tg$id])
multiExpr[[2]] = list(data = datExprm[,tg$id])
multiExpr[[3]] = list(data = datExprh[,tg$id])
setLabels = gts
names(multiExpr) = setLabels
lapply(multiExpr, lapply, dim)

## MEs
multiMEs = multiSetMEs(multiExpr, universalColors = tg$grp)
MEs = multiMEs[[1]]$data

# plot MEs
td = t(MEs)
colnames(td) = rownames(datExprb)
td = td[,tissues]

drows1 <- "correlation"
dcols1 <- "correlation"
col.pal <- brewer.pal(9, "Blues")
col.geno = brewer.pal(8, "Paired")[6:4]
col.tissue = c(brewer.pal(8, 'Dark2'), brewer.pal(9, 'Set1'))
ann_colors = list(
    Genotype = col.geno,
    Tissue = col.tissue
)

hm.parameters <- list(td,
  color = col.pal,
  cellwidth = 5, cellheight = 5, scale = "none",
  treeheight_row = 150,
  kmeans_k = NA,
  show_rownames = T, show_colnames = F,
#  main = "Heatmap of asinh(FPKM)",
  clustering_method = "complete",
  cluster_rows = T, cluster_cols = F,
  clustering_distance_rows = drows1,
  clustering_distance_cols = dcols1,
  #annotation_col = ta,
  annotation_colors = ann_colors,
  #gaps_col = c(17,34),
  fontsize_row = 6
)
fo = sprintf("%s/17.eigengene.pdf", dirw)
do.call("pheatmap", c(hm.parameters, filename=fo))

# plot single-module expr heatmap
mid = 20
fp = sprintf("%s/routs/%d.pdf", dirw, mid)
fpkm_heatmap(tim$gid[which(mids == mid)], ti, tissues, gts, fp)

# merge close MEs
MEDiss = 1 - cor(MEs)
METree = hclust(as.dist(MEDiss), method = "average")

fp3 = sprintf("%s/18.mergemodule.pdf", dirw)
pdf(fp3, width = 12, height = 6)
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")

MEDissThres = 0.1
abline(h=MEDissThres, col = "red")
merge = mergeCloseModules(datExprb[,tg$id], tg$grp, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors;
mergedMEs = merge$newMEs;
dev.off()

h = moduleMergeUsingKME(datExprb[,tg$id], tg$grp)
#}}}

#{{{ construct WGCNA datasets
gt = gts[1]
nSets = 3

multiExpr = list()
multiExpr[[1]] = list(data = datExprb[,tg$id])
multiExpr[[2]] = list(data = datExprm[,tg$id])
multiExpr[[3]] = list(data = datExprh[,tg$id])
setLabels = gts
names(multiExpr) = setLabels
lapply(multiExpr, lapply, dim)
multiMEs = multiSetMEs(multiExpr, universalColors = tg$grp)

## module preservation - run on HPC
colorList = list(tg$grp, tg$grp, tg$grp)
names(colorList) = setLabels

mp = modulePreservation(multiExpr, colorList, networkType = 'unsigned', corFnc = "cor",
	referenceNetworks = 1, parallelCalculation = TRUE, 
	loadPermutedStatistics = FALSE, nPermutations = 200, verbose = 3)
fo = sprintf("%s/30.mp.rda", dirw)
save(mp, file = fo)
#}}}

fo = sprintf("%s/30.mp.rda", dirw)
x = load(fo)

#{{{ MP plot
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
#}}}

#{{{ compare Mo17 and BxM
tz = data.frame(modName = modNames, moduleSize = modSizes, 
	Mo17 = as.numeric(mp$preservation$Z[[1]][[2]]$Zsummary.pres), 
	B73xMo17 = as.numeric(mp$preservation$Z[[1]][[3]]$Zsummary.pres)
)
tz = tz[! tz$modName %in% c('0', '0.1'),]

p1 = ggplot(tz, aes(x = Mo17, y = B73xMo17, size = moduleSize, label = modName)) +
  geom_point() +
  geom_text(vjust = 0, nudge_y = 0.2) +
  scale_x_continuous(name = "Mo17 preservation Z-summary", limits = c(0, 20)) +
  scale_y_continuous(name = "B73xMo17 preservation Z-summary", limits = c(0, 20)) +
  geom_hline(yintercept = 2, color = 'blue', linetype = 2) + 
  geom_hline(yintercept = 10, color = 'darkgreen', linetype = 2) + 
  geom_vline(xintercept = 2, color = 'blue', linetype = 2) + 
  geom_vline(xintercept = 10, color = 'darkgreen', linetype = 2) + 
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "lines")) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.title.y = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "black", angle = 0))
fp = sprintf("%s/36.diffpres.z.pdf", dirw)
ggsave(p1, filename = fp, width = 9, height = 8)

tm = data.frame(modName = modNames, moduleSize = modSizes, 
	Mo17 = as.numeric(mp$preservation$observed[[1]][[2]]$medianRank.pres), 
	B73xMo17 = as.numeric(mp$preservation$observed[[1]][[3]]$medianRank.pres)
)
tm = tm[! tm$modName %in% c('0', '0.1'),]

p1 = ggplot(tm, aes(x = Mo17, y = B73xMo17, size = moduleSize, label = modName)) +
  geom_point() +
  scale_x_continuous(name = 'Mo17 preservation MedianRank') +
  scale_y_continuous(name = 'B73xMO17 preservation MedianRank') +
  geom_text(vjust = 0, nudge_y = 3) +
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "lines")) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.title.y = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "black", angle = 0))
fp = sprintf("%s/36.diffpres.medianrank.pdf", dirw)
ggsave(p1, filename = fp, width = 7, height = 6)
#}}}

#{{{ circle plots
pids = c(36,44,207, 90,105,161)
pids = c(38,29)
pids = c(9)

for (pid in pids) {
mid = sprintf("ME%d", pid)
pathlabel1 = sprintf("%d: Zsummary = %.02f; medianRank = %d", pid,
	tz$Mo17[tz$modName == as.character(pid)], 
	tm$Mo17[tm$modName == as.character(pid)])
pathlabel2 = sprintf("%d: Zsummary = %.02f; medianRank = %d", pid,
	tz$B73xMo17[tz$modName == as.character(pid)], 
	tm$B73xMo17[tm$modName == as.character(pid)])

idxGenes = which(tg$grp == pid)
nPathGenes = length(idxGenes)
pathwayAdjs = list()
KMEpathway = matrix(0, nPathGenes, nSets)
for (set in 1:nSets)
{
  bc = bicor(multiExpr[[set]]$data[, idxGenes], use = "p")
  pathwayAdjs[[set]] = abs(bc)^4 * sign(bc)
  KMEpathway[, set] = bicor(multiExpr[[set]]$data[, idxGenes], multiMEs[[set]]$data[, mid], use = "p")
}

conn = matrix(0, nPathGenes, nSets)
for (set in 1:nSets)
  conn[, set] = apply(abs(pathwayAdjs[[set]]), 2, sum)-1
weights = c(3,1,5,1, 3,1,5,1);
wMat = matrix(weights, nPathGenes, nSets, byrow = TRUE)
wconn = apply(conn * wMat, 1, sum)
order = order(-wconn)
# use the gene names as lables
labels = colnames(pathwayAdjs[[ref]])
labels = substr(labels, 7, 14)

myfunc <- function(x,y) {
	echo 'hello world'
}
fo = sprintf("%s/38.circleplot/%d.pdf", dirw, pid)
pdf(file = fo, wi=16, h=4)
par(mfrow =c(1,4));
par(mar = c(0.3, 0.2, 1.5, 0.2))
for (set in 1:nSets)
{
  circlePlot(pathwayAdjs[[set]], labels, order, main = setLabels[set],
            variable.cex.labels = TRUE,
            radii = c(0.56, 0.62), center = c(0.1, 0.04),
            min.cex.labels = 1.2, max.cex.labels = 1.4, cex.main = 1.4)
  if(set == 2) {
    text(0, -1, pathlabel1)
  } else if(set == 3) {
  	text(0, -1, pathlabel2)
  }
}

par(mar = c(3.3, 3.3, 4, 0.5));
par(mgp = c(1.9, 0.6, 0))
verboseScatterplot(KMEpathway[, ref], KMEpathway[, test],
                   xlab = sprintf("KME in %s", setLabels[ref]),
                   ylab = sprintf("KME in %s", setLabels[test]),
                   main = sprintf("%s. KME in %s vs. %s", 
                       LETTERS[1], setLabels[test], setLabels[ref]),
                   cex.main = 1.4, cex.lab = 1.2, cex.axis = 1.2, abline = TRUE
)
dev.off()
}
#}}}

#{{{# convert camoco coex table to R symmetric matrix
fc = file.path(dirw, "coex.csv")
vec = scan(fc, what = numeric())

roots = polyroot(c(-length(vec), -0.5, 0.5))
ng = Re(roots)[Re(roots) > 0]
stopifnot(ng*(ng-1)/2 == length(vec))

coex <- matrix(rep(0, ng*ng), nrow=ng)
coex[lower.tri(coex)] = vec
coex = t(coex)
coex[lower.tri(coex)] = vec
save(coex, file = file.path(dirw, "camoco.rda"))
#}}}


