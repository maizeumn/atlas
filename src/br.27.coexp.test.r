source("br.fun.R")

dirw = '/home/springer/zhoux379/scratch/briggs2/47.coexp.test'

fi = file.path(dirw, "../36.long.filtered.tsv")
ti = read.table(fi, sep = "\t", header = T, as.is = T)

gt = gts[1]
tiw = spread(ti[ti$Genotype == gt, -c(3,4)], Tissue, fpkm)

expr = t(as.matrix(tiw[,-1]))
colnames(expr) = tiw[,1]
gids = tiw[,1]


ng = length(gids)
dummy <- matrix(rep(NA, ng*ng), nrow=ng)
ii = which(lower.tri(dummy))
colidx = as.integer((ii - 1) / ng) + 1
rowidx = ii - (colidx - 1) * ng


opt = 1
## camoco
if (opt == 1) {

pcc.matrix = cor(asinh(expr), method = 'pearson')
pcc = pcc.matrix[lower.tri(pcc.matrix)]

pcc[pcc == 1] = 0.9999
pcc[pcc == -1] = -0.9999
#pcc2 = log((1+pcc) / (1-pcc)) / 2
pcc2 = atanh(pcc)
pccz = (pcc2 - mean(pcc2)) / sd(pcc2)

fp = sprintf("%s/01.edgeweight/pccz.rda", dirw)
save(pccz, file = fp)

fp = sprintf("%s/01.edgeweight/pccz.pdf", dirw)
pdf(fp, width = 8, height = 8)
hist(pccz, 100)
dev.off()

sigidx = which(pccz >= 3)
dw = data.frame(g1 = gids[rowidx[sigidx]], g2 = gids[colidx[sigidx]], pccz = pccz[sigidx])
fw = sprintf("%s/02.edges/C.tsv", dirw)
write.table(dw, fw, sep = "\t", row.names = F, col.names = F, quote = F)

nets = sprintf("C%d", 1:5)
params = c(1.4, 2, 4, 5, 6)
for (i in 1:length(nets)) {
	fm = sprintf("%s/04.tmp/%s.txt", dirw, nets[i])
	cmd = sprintf("mcl %s --abc -scheme 7 -I %g -o %s", fw, params[i], fm)
	system(cmd)
	
	fo = sprintf("%s/05.modules/%s.tsv", dirw, nets[i])
	cmd = sprintf("mcl2tsv.py %s %s", fm, fo)
	system(cmd)
}

}

## clusterOne
opt = 2
if (opt == 2) {

pcc.matrix = cor(asinh(expr), method = 'pearson')
pcc = pcc.matrix[lower.tri(pcc.matrix)]

ii = which(lower.tri(pcc.matrix))
colidx = as.integer((ii - 1) / nrow(pcc.matrix)) + 1
rowidx = ii - (colidx - 1) * nrow(pcc.matrix)

ng = nrow(tiw)
coex <- matrix(rep(-1, ng*ng), nrow=ng)
coex[lower.tri(coex)] = pcc
coex = t(coex)
coex[lower.tri(coex)] = pcc

rankmatrix = t(apply(-coex, 1, rank))
mr = sapply(1:length(pcc), myfunc <- function(i) 
	sqrt(rankmatrix[rowidx[i], colidx[i]] * rankmatrix[colidx[i], rowidx[i]])
)

nets = sprintf("N%d", 1:5)
params = c(5, 10, 25, 50, 100)
for (i in 1:length(nets)) {
	net = nets[i]
	ew = exp((-(mr-1)/params[i]))
	fe = sprintf("%s/01.edgeweight/%s.rda", dirw, net)
	save(ew, file = fe)

	fp = sprintf("%s/01.edgeweight/%s.pdf", dirw, net)
	pdf(fp, width = 8, height = 8)
	hist(ew[ew>=0.01], 100)
	dev.off()

	idxs = which(ew >= 0.01)
	to = data.frame(g1 = gids[rowidx[idxs]], g2 = gids[colidx[idxs]], ew = ew[idxs])
	fo = sprintf("%s/02.edges/%s.txt", dirw, net)
	write.table(to, fo, sep = " ", row.names = F, col.names = F, quote = F)
	
	fc = sprintf("%s/04.tmp/%s.csv", dirw, net)
	cmd = sprintf("java -jar $src/cluster_one-1.0.jar -f edge_list -F csv %s > %s", fo, fc)
	system(cmd)

	fm = sprintf("%s/05.modules/%s.tsv", dirw, net)
	cmd = sprintf("one2tsv.py %s %s", fc, fm)
	system(cmd)
}

}

## WGCNA
if (opt == 3) {

}

## GO enrichment
fg = '/home/springer/zhoux379/data/genome/Zmays_v4/61.interpro/15_lev2.tsv'
tg = read.table(fg, header = F, as.is = T, sep = "\t")
colnames(tg) = c("gid", "goid", "goterm")
tg = tg[tg$gid %in% gids,]

	net = "C5"
	fm = sprintf("%s/05.modules/%s.tsv", dirw, net)
	tm = read.table(fm, header = T, as.is = T, sep = "\t")

tz = merge(tg, tm, by = 'gid', all.x = T)
tzs = ddply(tz, .(goid), summarise, goterm = goterm[1], msize = length(gid), gim = sum(!is.na(grp)), nmodule = length(unique(grp[!is.na(grp)])), gnim = sum(is.na(grp)))

total_gene_in_module = nrow(tm)
total_gene_not_in_module = length(gids) - total_gene_in_module
pvals = apply(tzs, 1, myfunc <- function(x) dhyper(as.numeric(x[4]), total_gene_in_module - as.numeric(x[6]), total_gene_not_in_module - as.numeric(x[6]), as.numeric(x[3])))
tzs = cbind(tzs, pval = p.adjust(pvals, 'BH'))
tzs = tzs[order(tzs$pval), ]