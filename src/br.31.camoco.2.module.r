source("br.fun.R")

dirw = file.path(Sys.getenv("misc2"), "briggs", "51.camoco")
#dirw = "/home/springer/zhoux379/scratch/briggs/51.camoco"

fi = file.path(dirw, "../36.long.filtered.tsv")
ti = read.table(fi, sep = "\t", header = T, as.is = T)

## construct expr matrix
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

## look at overlap with corncyc pathways
tp2 = tp[tp$gid %in% tiw[,1], c('gid','p2')]
tp3 = unique(tp2)
tz = merge(tp3, tg, by = 'gid', all.x = T)
tzs = ddply(tz, .(p2), summarise, psize = length(gid), gim = sum(!is.na(grp)), ngrp = length(unique(grp[!is.na(grp)])), gnim = sum(is.na(grp)))

total_gene_in_module = nrow(tg)
total_gene_not_in_module = nrow(tib) - total_gene_in_module
pvals = apply(tzs, 1, myfunc <- function(x) dhyper(as.numeric(x[3]), total_gene_in_module - as.numeric(x[3]), total_gene_not_in_module - as.numeric(x[5]), as.numeric(x[2])))
tzs = cbind(tzs, pval = p.adjust(pvals, 'BH'))

## construct WGCNA datasets
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


