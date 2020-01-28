args <- commandArgs(TRUE)
i = as.numeric(args[1])

source("br.fun.R")

dirw = file.path(Sys.getenv("misc2"), "briggs", "61.perm")

fi = file.path(dirw, "../35.long.tsv")
ti = read.table(fi, sep = "\t", header = T, as.is = T)

gts = c("B73", "Mo17", "B73xMo17")
ti2 = ti[ti$Genotype %in% gts,]
ti3 = cbind(ti2[,-c(2:4)], treat = sprintf("%s.%s", ti2$Genotype, ti2$Tissue))
ti4 = spread(ti3, treat, fpkm)

set.seed(i)
ti41 = ti4
if (i > 0) {
	treats = colnames(ti4)[-1]
	tmp = strsplit(treats, , split = "[.]")
	geno = sapply(tmp, "[", 1)
	tiss = sapply(tmp, "[", 2)
	for (j in 1:length(unique(tiss))) {
		idxs = which(tiss == tiss[j])
		geno[idxs] = sample(geno[idxs])
	}
	colnames(ti41)[-1] = sprintf("%s.%s", geno, tiss)
}

ti5 = gather(ti41, treat, fpkm, -gid)
tmp = strsplit(ti5$treat, , split = "[.]")
geno = sapply(tmp, "[", 1)
tiss = sapply(tmp, "[", 2)
ti6 = cbind(ti5[,-2], Genotype = geno, Tissue = tiss)

## QC
grp = dplyr::group_by(ti6, gid, Genotype)
ti7 = dplyr::summarise(grp, nexp = sum(fpkm>=1))
ti8 = spread(ti7, Genotype, nexp)
gids = ti8$gid[ti8$B73 >= 1 & ti8$Mo17 >= 1 & ti8$B73xMo17 >= 1]
cat(sprintf("rep %4d: %d\n", i, length(gids)))
tqc = ti6[ti6$gid %in% gids,]

## construct data matrix
multiExpr = list()
for (gt in gts) {
	tiw = spread(tqc[tqc$Genotype == gt, c('gid','Tissue','fpkm')], Tissue, fpkm)
	datExpr = t(as.matrix(tiw[,-1]))
	colnames(datExpr) = tiw[,1]
	multiExpr[[gt]] = asinh(datExpr)
}

## read module info
fm = file.path(dirw, "../59.allmodules.tsv")
tm = read.table(fm, header = T, sep = "\t", as.is = T, quote = '')
tm = tm[tm$gid %in% gids,]

grp = dplyr::group_by(tm, mid)
tms = dplyr::summarise(grp, size = length(gid))
tms = tms[tms$size >= 5,]
tm = tm[tm$mid %in% tms$mid,]

## compute stats
comp = list()
for (gt1 in gts) {
	datExpr1 = multiExpr[[gt1]]
	coexv1 = get_cormatrix(datExpr1)
	for (gt2 in gts) {
		if (gt1 == gt2) next
		datExpr2 = multiExpr[[gt2]]
		coexv2 = get_cormatrix(datExpr2)
		
		pcc.expr = cor(c(datExpr1), c(datExpr2))
		pcc.coex = cor(coexv1, coexv2)
		comp[[gt1]][[gt2]] = list(
			pcc.expr = pcc.expr,
			pcc.coex = pcc.coex,
			stat = module_preservation_stats(datExpr1, datExpr2, tm)
		)
		cat(gt1, gt2, 'done\n', sep = ' ')
	}
}


fo = sprintf("%s/01.perm/%04d.rda", dirw, i)
save(gids, comp, file = fo)


