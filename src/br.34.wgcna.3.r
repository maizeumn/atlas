source("br.fun.R")

dirw = '/home/springer/zhoux379/scratch/briggs2/52.wgcna'

fi = file.path(dirw, "../36.long.filtered.tsv")
ti = read.table(fi, sep = "\t", header = T, as.is = T)

#####
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


## module preservation btw B and M - cross tabulation
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


## module preservation
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

## make some heatmaps

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