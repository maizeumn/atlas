source("br.fun.R")

dirw = file.path(Sys.getenv("misc2"), "briggs", "61.corncyc")
#dirw = "/home/springer/zhoux379/scratch/briggs/61.corncyc"

fi = file.path(dirw, "../36.long.filtered.tsv")
ti = read.table(fi, sep = "\t", header = T, as.is = T)

### read pathway information
fc = '/home/springer/zhoux379/data/genome/Zmays_v4/corncyc/01.tsv'
tc = read.table(fc, header = F, sep = "\t", quote = "")
gids = strsplit(tc$V7, split = ",")
lens = sapply(gids, length)
tc = data.frame(pl = rep(tc$V1, lens), 
	p1 = rep(tc$V2, lens), p2 = rep(tc$V3, lens), 
	p3 = rep(tc$V4, lens), p4 = rep(tc$V5, lens), 
	ps = rep(tc$V6, lens), gid = unlist(gids), stringsAsFactors = F)
table(tc$p1)

tp = tc[tc$p1 %in% c('Biosynthesis', 'Degradation/Utilization/Assimilation'),]
table(tp$p2)

#tp2 = merge(tp[,c('gid','pl')], tm, by = 'gid', all.x = T)
#tp3 = ddply(tp2, .(p2), summarise, psize = length(gid), gim = sum(!is.na(mid)), nmid = length(unique(mid[!is.na(mid)])))

tc = tc[,c('pl','gid')]
colnames(tc) = c("pathway", "gid")
tc = tc[tc$gid %in% unique(ti$gid),]
tcc = ddply(tc, .(pathway), summarize, ngene = length(gid))
tc2 = merge(tc, tcc[,c('pathway','ngene')], by = 'pathway')
tc3 = ddply(tc2, .(gid), summarise, pathway = pathway[which(ngene==max(ngene))[1]])

tcc = ddply(tc3, .(pathway), summarize, ngene = length(gid))
tcc = tcc[tcc$ngene >= 3,]
tcc = tcc[order(tcc$ngene, decreasing = T),]
tcc = cbind(tcc, pid = 1:nrow(tcc))
tc = merge(tc3, tcc, by = 'pathway')

#tp = data.frame(gid = unique(ti$gid), stringsAsFactors = F)
#tp = merge(tp, tc[,c('gid','pid')], all.x = T)
#tp$pid[is.na(tp$pid)] = 0

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

## construct WGCNA datasets
nSets = 3

multiExpr = list()
multiExpr[[1]] = list(data = datExprb[,tc$gid])
multiExpr[[2]] = list(data = datExprm[,tc$gid])
multiExpr[[3]] = list(data = datExprh[,tc$gid])
setLabels = gts
names(multiExpr) = setLabels
lapply(multiExpr, lapply, dim)

colorList = list('1' = tc$pid, "2" = tc$pid, "3" = tc$pid)

mp = modulePreservation(multiExpr, colorList, networkType = 'signed', corFnc = "cor",
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
modNames = rownames(mp$preservation$observed[[ref]][[test]])
modSizes = mp$preservation$Z[[ref]][[test]][, 1];

stats = cbind(modName = modNames, moduleSize = modSizes, 
	mp$quality$observed[[ref]][[test]][,-1], 
	mp$preservation$observed[[ref]][[test]][, -1],
	mp$quality$Z[[ref]][[test]][,-1], 
	mp$preservation$Z[[ref]][[test]][,-1]
)

t1 = stats[!stats$modName %in% c('0','0.1'), c('modName','moduleSize','medianRank.pres','Zsummary.pres')]
p1 = ggplot(t1, aes(x = medianRank.pres, y = Zsummary.pres, size = moduleSize, label = modName)) +
  geom_point() +
  geom_text(vjust = 0, nudge_y = 0.2) +
  #scale_x_continuous(name = '') +
  #scale_y_continuous(name = '') +
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

t1[order(t1$Zsummary.pres),][1:20,]

## compare Mo17 and BxM
tz = data.frame(modName = modNames, moduleSize = modSizes, 
	B73 = as.numeric(mp$preservation$Z[[1]][[2]]$Zsummary.pres), 
	B73xMo17 = as.numeric(mp$preservation$Z[[1]][[3]]$Zsummary.pres)
)
tz = tz[! tz$modName %in% c('0', '0.1'),]

p1 = ggplot(tz, aes(x = B73, y = B73xMo17, size = moduleSize, label = modName)) +
  geom_point() +
  geom_text(vjust = 0, nudge_y = 0.2) +
  #scale_x_continuous(name = '') +
  #scale_y_continuous(name = '') +
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
ggsave(p1, filename = fp, width = 7, height = 6)

tz = data.frame(modName = modNames, moduleSize = modSizes, 
	B73 = as.numeric(mp$preservation$observed[[1]][[2]]$medianRank.pres), 
	B73xMo17 = as.numeric(mp$preservation$observed[[1]][[3]]$medianRank.pres)
)
tz = tz[! tz$modName %in% c('0', '0.1'),]

p1 = ggplot(tz, aes(x = B73, y = B73xMo17, size = moduleSize, label = modName)) +
  geom_point() +
  geom_text(vjust = 0, nudge_y = 3) +
  #scale_x_continuous(name = '') +
  #scale_y_continuous(name = '') +
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "lines")) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.title.y = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "black", angle = 0))
fp = sprintf("%s/36.diffpres.medianrank.pdf", dirw)
ggsave(p1, filename = fp, width = 7, height = 6)


## circle plot + KME plot
if(FALSE) {
multiMEs = multiSetMEs(multiExpr, universalColors = tc$pid)

pids = c(5,10,12,41,45,63,73,89,93,96,102,126)

for (pid in pids) {
#pid = 73
mid = sprintf("ME%d", pid)
pathlabel = sprintf("Zsummary = %.02f; medianRank = %d\n%d. %s", 
	t1$Zsummary.pres[t1$modName == as.character(pid)], 
	t1$medianRank.pres[t1$modName == as.character(pid)], pid, tc$pathway[tc$pid==pid])

idxGenes = which(tp$pid == pid)
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

fo = sprintf("%s/88.circleplot.%d.pdf", dirw, pid)
pdf(file = fo, wi=16, h=4)
par(mfrow =c(1,4));
par(mar = c(0.3, 0.2, 1.5, 0.2))
for (set in 1:nSets)
{
  circlePlot(pathwayAdjs[[set]], labels, order, main = setLabels[set],
            variable.cex.labels = TRUE,
            radii = c(0.56, 0.62), center = c(0.1, 0.04),
            min.cex.labels = 1.2, max.cex.labels = 1.4, cex.main = 1.4)
  if(set == 1)
	text(0,-1, pathlabel)
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

}