source("br.fun.R")

dirw = '/home/springer/zhoux379/scratch/briggs2/52.wgcna'

fi = file.path(dirw, "../36.long.filtered.tsv")
ti = read.table(fi, sep = "\t", header = T, as.is = T)

## read data ----
fx = file.path(diro, "09.B73.rda")
x = load(fx)
x

tom1 = TOM
tm1 = tm
wg_tree1 = wg_tree

fx = file.path(diro, "09.Mo17.rda")
x = load(fx)
x

tom2 = TOM
tm2 = tm
wg_tree2 = wg_tree

# module eigengenes ----
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

## module preservation
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


### process outputs
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


## circle plot
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


### clusterRepro
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

