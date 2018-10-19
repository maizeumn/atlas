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

mids = rep(0, nrow(tiw))
mids[tg$id] = tg$grp

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
multiMEs = multiSetMEs(multiExpr, universalColors = tg$grp)

## module preservation - run on HPC
colorList = list(tg$grp, tg$grp, tg$grp)
names(colorList) = setLabels

mp = modulePreservation(multiExpr, colorList, networkType = 'unsigned', corFnc = "cor",
	referenceNetworks = 1, parallelCalculation = TRUE, 
	loadPermutedStatistics = FALSE, nPermutations = 200, verbose = 3)
fo = sprintf("%s/30.mp.rda", dirw)
save(mp, file = fo)

### read outputs
fo = sprintf("%s/30.mp.rda", dirw)
x = load(fo)
x

## MP plot
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

## compare Mo17 and BxM
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


## circle plots
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
