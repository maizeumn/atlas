source("br.fun.R")
dirw = file.path(dirp, "51.camoco")

fi = file.path(dirw, "../36.long.filtered.tsv")
ti = read.table(fi, sep = "\t", header = T, as.is = T)

### for how to run camoco on a workstation - see Note

#{{{# obtain camoco results from scratch
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
