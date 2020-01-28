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
if('0.1' %in% stats$modName) {
	stats = stats[-which(stats$modName == '0.1'),]
}
dim(stats)
stats$modName = as.integer(stats$modName)
stats = stats[order(stats$modName),]
stats[1:5,1:5]

colors = tc$pid
datExpr1 = datExprb[,tc$gid]
datExpr2 = datExprm[,tc$gid]
ME1 = moduleEigengenes(datExpr1, colors)$eigengenes
ME2 = moduleEigengenes(datExpr2, colors)$eigengenes
k1 = signedKME(datExpr1, ME1, corFnc = 'cor')
k2 = signedKME(datExpr2, ME2, corFnc = 'cor')

p1 = propVarExplained(datExpr1, colors, ME1, corFnc = 'cor')
p2 = propVarExplained(datExpr2, colors, ME2, corFnc = 'cor')

k11 = gather(cbind(gid = rownames(k1), k1), pid, kME, -gid)
k11$pid = as.integer(substr(k11$pid, 4, nchar(k11$pid)))
k12 = merge(tc[,c('gid','pid')], k11, by = c('gid', 'pid'))
stopifnot(nrow(k12) == nrow(tc))

k21 = gather(cbind(gid = rownames(k2), k2), pid, kME, -gid)
k21$pid = as.integer(substr(k21$pid, 4, nchar(k21$pid)))
k22 = merge(tc[,c('gid','pid')], k21, by = c('gid', 'pid'))
stopifnot(nrow(k22) == nrow(tc))

k3 = cbind(k12, kME2 = k22$kME)
k31 = ddply(k3, .(pid), summarise,
	meanSignAwareKME.qual = mean(sign(kME)*kME),
	meanSignAwareKME.pres = abs(mean(sign(kME)*kME2)),
	cor.kME = abs(cor(kME, kME2)))

meanCor.qual = c(); meanCor.pres = c()
meanAdj.qual = c(); meanAdj.pres = c()
cor.kIM = c(); cor.kME = c(); cor.cor = c()
for (mid in sort(unique(colors))) {
	dat1 = datExpr1[,tc$gid[tc$pid == mid]]
	dat2 = datExpr2[,tc$gid[tc$pid == mid]]
	cormat1 = cor(dat1, use = 'p')
	cormat2 = cor(dat2, use = 'p')
	corvec1 = cormat1[lower.tri(cormat1)]
	corvec2 = cormat2[lower.tri(cormat2)]
	
	adjmat1 = ((1+cormat1)/2)^12
	adjmat2 = ((1+cormat2)/2)^12
	adjvec1 = adjmat1[lower.tri(adjmat1)]
	adjvec2 = adjmat2[lower.tri(adjmat2)]
	
	meanCor.qual = c(meanCor.qual, mean(corvec1))
	meanCor.pres = c(meanCor.pres, mean(sign(corvec1) * corvec2))
	meanAdj.qual = c(meanAdj.qual, mean(adjvec1))
	meanAdj.pres = c(meanAdj.pres, mean(adjvec2))
	
	den1 = apply(as.matrix(adjmat1), 1, sum) - 1 
	den2 = apply(as.matrix(adjmat2), 1, sum) - 1 
	cor.kIM = c(cor.kIM, cor(den1, den2))
	cor.cor = c(cor.cor, cor(corvec1, corvec2))
}


ty = data.frame(pid = sort(unique(colors)), 
	propVarExplained.qual = p1,
	propVarExplained = p2, 
	meanSignAwareKME.qual = k31[,2],
	meanSignAwareKME = k31[,3],
	meanCor.qual = meanCor.qual,
	meanCor = meanCor.pres,
	meanAdj.qual = meanAdj.qual,
	meanAdj = meanAdj.pres,
	cor.kIM = cor.kIM,
	cor.kME = k31[,4],
	cor.cor = cor.cor,
	stringsAsFactors = F)
sum(ty$propVarExplained.qual - stats$propVarExplained.qual)
sum(ty$propVarExplained - stats$propVarExplained.pres)
sum(ty$meanSignAwareKME.qual - stats$meanSignAwareKME.qual)
sum(ty$meanSignAwareKME - stats$meanSignAwareKME.pres)
sum(ty$meanCor.qual - stats$meanSignAwareCorDat.qual)
sum(ty$meanCor - stats$meanSignAwareCorDat.pres)
sum(ty$meanAdj.qual - stats$meanAdj.qual)
sum(ty$meanAdj - stats$meanAdj.pres)
sum(ty$cor.kIM - stats$cor.kIM)
sum(ty$cor.kME - stats$cor.kME)
sum(ty$cor.cor - stats$cor.cor)

vnames = c('propVarExplained.pres', 'meanSignAwareKME.pres', 
	'meanSignAwareCorDat.pres', 'meanAdj.pres', 'cor.kIM', 'cor.kME', 'cor.cor')
mnames = c('propVarExplained', 'meanSignAwareKME', 
	'meanCor', 'meanAdj', 'cor.kIM', 'cor.kME', 'cor.cor')
for (mname in mnames) {
	ty = cbind(ty, y = NA)
	colnames(ty)[ncol(ty)] = sprintf("Z.%s", mname)
}
for (mid in sort(unique(colors))) {
	mid2 = as.character(mid)
	for (i in 1:length(vnames)) {
		vname = vnames[i]; mname = mnames[i]
		vnamez = sprintf("Z.%s", vname)
		mnamez = sprintf("Z.%s", mname)
		perm.values = mp$permutationDetails$permutedStatistics[[1]][[2]]$regStats[mid2,vname,]
		perm.mean = mean(perm.values)
		perm.sd = sd(perm.values)
		zvalue = (ty[ty$pid == mid, mname] - perm.mean) / perm.sd
		ty[ty$pid == mid, mnamez] = zvalue
	}
}
for (i in 1:length(vnames)) {
	vname = vnames[i]; mname = mnames[i]
	vnamez = sprintf("Z.%s", vname)
	mnamez = sprintf("Z.%s", mname)
	cat(sum(ty[,mnamez] - stats[,vnamez]), "\n")
}
Zsummary = apply(ty, 1, myfunc <- function(rw) (median(rw[c(13,14,15,16)]) + median(rw[c(17:19)])) / 2)
sum(Zsummary - stats$Zsummary.pres)