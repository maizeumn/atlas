#{{{
require(devtools)
load_all('~/git/rmaize')
require(ggforce)
require(ggpubr)
#require(ggExtra)
dira = '~/projects/assets'
dirp = '~/projects/atlas'
dird = file.path(dirp, 'data')
dirf = file.path(dird, '95_figures')
cols17 = c(brewer.pal(8, 'Dark2'), brewer.pal(9, 'Set1'))
cols23 = c(pal_npg()(10), pal_simpsons()(13))
tissues23 = c("seedlingleaf_11DAS", "blade_v12", "flagleaf_0DAP",
            "auricle_v12", "sheath_v12", "husk_0DAP", "tasselstem_0DAP",
            "internode_v12", "root_0DAP", "silk_0DAP", "floret_0DAP",
            "radicle_root", "seedlingroot_11DAS", "tassel_v12", "spikelets_0DAP",
            "coleoptile_tip", "seedlingmeristem_11DAS", "ear_v14",
            "embryo_27DAP", "embryo_imbibedseed", "endosperm_27DAP",
            "kernel_14DAP", "endosperm_14DAP")
tissues23 = c("seedlingleaf_11DAS", "blade_v12", "flagleaf_0DAP",
            "husk_0DAP", "sheath_v12", "auricle_v12",
            "floret_0DAP", "tasselstem_0DAP",
            "internode_v12",
            "root_0DAP", "seedlingroot_11DAS", "radicle_root",
            "coleoptile_tip",
            "silk_0DAP",
            "tassel_v12", "spikelets_0DAP", "ear_v14",
            "seedlingmeristem_11DAS",
            "embryo_27DAP", "embryo_imbibedseed",
            "endosperm_27DAP", "kernel_14DAP", "endosperm_14DAP")
tiss23 = c("radicle_root", "seedling_root", "seedling_leaf", "seedling_meristem", "coleoptile", "auricle", "blade_leaf", "internode", "tassel", "sheath", "ear", "flag_leaf", "floret", "husk", "root", "silk", "spikelet", "tassel_stem", "endosperm14D", "kernel", "embryo", "seed_imbibed", "endosperm27D")
tiss20 = tissues23[1:20]
tiss2 = c("tasselstem_0DAP", "internode_v12")
gts3 = c("B73", "Mo17", 'BxM')
tiss3 = c("seedling_root", "coleoptile", "leaf_V2")
gts10 = c("B73",'B84',"Mo17",'A682','B73xMo17','Mo17xB73',
          'B84xB73','B84xMo17','A682xB73','A682xMo17')
#require(BiocParallel)
#bpparam <- MulticoreParam()
#bplog(bpparam) <- T

#require(WGCNA)
##allowWGCNAThreads()
#enableWGCNAThreads()
#}}}

cap_bigint <- function(x, int.max=2147483647) ifelse(x>int.max, int.max, x)
trim_range <- function(x, range.min=-2147483647, range.max=2147483647) ifelse(x>range.max, range.max, ifelse(x<range.min, range.min, x))
plot_deseq2 <- function(dds, dirw, tissue) {
    #{{{
    require(DESeq2)
    rld <- DESeq2::rlog(dds, blind = F)
    fp = sprintf("%s/02_pca/%s.pdf", dirw, tissue)
    pdf(fp, width=5, height=5)
    plotPCA(rld, intgroup = "Genotype", ntop = 10000)
    dev.off()
    #
    fp = sprintf("%s/03_disp/%s.pdf", dirw, tissue)
    pdf(fp, width=5, height=5)
    plotDispEsts(dds)
    dev.off()
    #
    resLFC = lfcShrink(dds, coef = 2)
    fp = sprintf("%s/05_ma/%s.pdf", dirw, tissue)
    pdf(fp, width=5, height=5)
    plotMA(resLFC, ylim = c(-3,3))
    dev.off()
    #}}}
}
get_cormatrix <- function(expr) {
    #{{{
    pcc.matrix = cor(expr, method = 'pearson')
    pcc = pcc.matrix[lower.tri(pcc.matrix)]
    pcc[pcc == 1] = 0.999999
    pcc[pcc == -1] = -0.999999
    #pcc2 = log((1+pcc) / (1-pcc)) / 2
    pcc2 = atanh(pcc)
    pcc3 = (pcc2 - mean(pcc2, na.rm = T)) / sd(pcc2, na.rm = T)
    head(pcc3)

    if(FALSE) {
        ng = ncol(expr)
        coex <- matrix(rep(0, ng*ng), nrow=ng)
        coex[lower.tri(coex)] = pcc3
        coex = t(coex)
        coex[lower.tri(coex)] = pcc3
        list(coexv = pcc3, coexm = coex)
    }
    pcc3
    #}}}
}
fpkm_heatmap <- function(gids, ti, tissues, gts, fo) {
    #{{{ make expression heatmap
    ti1 = ti[ti$gid %in% gids & ti$Genotype %in% gts,-4]
    ti1$Tissue = factor(ti1$Tissue, levels = tissues)
    ti1$Genotype = factor(ti1$Genotype, levels = gts)
    ti2 = ti1[order(ti1$Genotype, ti1$Tissue),]
    ti2 = cbind(ti2, cond = sprintf("%s.%s", ti2$Genotype, ti2$Tissue))
    ti3 = ti2[,c(1,4,5)]
    ti3$cond = factor(ti3$cond, levels = unique(ti3$cond))
    ti4 = spread(ti3, cond, fpkm)
    td = ti4
    rownames(td) = td$gid
    td = td[,-1]
    td = asinh(td)

    ta = unique(ti2[,c(2:3,5)])
    rownames(ta) = ta$cond
    ta = ta[,-3]

    drows1 <- "correlation"
    dcols1 <- "correlation"
    col.pal <- brewer.pal(9, "Blues")
    col.geno = brewer.pal(8, "Paired")[6:4]
    col.tissue = c(brewer.pal(8, 'Dark2'), brewer.pal(9, 'Set1'))
    names(col.tissue) = unique(ta$Tissue)
    names(col.geno) = unique(ta$Genotype)
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
        main = "Heatmap of asinh(FPKM)",
        clustering_method = "complete",
        cluster_rows = T, cluster_cols = F,
        clustering_distance_rows = drows1, 
        clustering_distance_cols = dcols1,
        annotation_col = ta,
        annotation_colors = ann_colors,
        gaps_col = c(17,34),
        fontsize_row = 6
    )
    do.call("pheatmap", c(hm.parameters, filename=fo))
    #}}}
}
make_circleplot <- function() {
    #{{{ single-module circle plot
	pathlabel1 = sprintf("%d: Zsummary = %.02f; medianRank = %d", pid,
		tz$Mo17[tz$modName == as.character(pid)], 
		tm$Mo17[tm$modName == as.character(pid)])
	pathlabel2 = sprintf("%d: Zsummary = %.02f; medianRank = %d", pid,
		tz$B73xMo17[tz$modName == as.character(pid)], 
		tm$B73xMo17[tm$modName == as.character(pid)])

	idxGenes = which(mids == pid)
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
    #}}}
}
make_hist <- function(vec, xmin, xmax, xitv) {
    #{{{
    xb = seq(xmin, xmax, by = xitv)
    xm <- xb[-length(xb)] + 0.5 * diff(xb)
    td = data.frame(x = xm)

    x = findInterval(vec, xb)
    d = as.data.frame.table(table(x))
    d$x = as.numeric(levels(d$x))[d$x]
    d = d[d$x > 0 & d$x <= length(xm),]
    d$x = xm[d$x]
    d[!is.na(d$x),]
    td = merge(td, d, by = 'x', all.x = T)
    td[is.na(td[,2]), 2] = 0
    colnames(td)[2] = "freq"
    td
    #}}}
}
module_stats <- function(datExpr, tm) {
    #{{{ compute module statistics
    mids = sort(unique(tm$mid))
    expr = datExpr[,tm$gid]
    ME = moduleEigengenes(expr, tm$mid)$eigengenes
    p1 = propVarExplained(expr, tm$mid, ME, corFnc = 'cor')
    meandens = c()
    for (mid in mids) {
        dat1 = datExpr[,tm$gid[tm$mid == mid]]
        cormat1 = cor(dat1, use = 'p')
        corvec1 = cormat1[lower.tri(cormat1)]
        meandens = c(meandens, mean(corvec1))
    }
    list(propVarExplained = p1, meanCorDensity = meandens)
    #}}}
}
module_preservation_stats <- function(datExpr1, datExpr2, tm) {
    #{{{ calculate module preservation statistics
    mids = sort(unique(tm$mid))
    expr1 = datExpr1[,tm$gid]
    expr2 = datExpr2[,tm$gid]
    ME1 = moduleEigengenes(expr1, tm$mid)$eigengenes
    ME2 = moduleEigengenes(expr2, tm$mid)$eigengenes
    k1 = signedKME(expr1, ME1, corFnc = 'cor')
    k2 = signedKME(expr2, ME2, corFnc = 'cor')
    colnames(k1) = substr(colnames(k1), 4, nchar(colnames(k1)))
    colnames(k2) = substr(colnames(k2), 4, nchar(colnames(k2)))
    kME1 = k1[as.matrix(tm[,c('gid','mid')])]
    kME2 = k2[as.matrix(tm[,c('gid','mid')])]
    stopifnot(length(kME1) == nrow(tm))
    stopifnot(length(kME1) == nrow(tm))

    p1 = propVarExplained(expr1, tm$mid, ME1, corFnc = 'cor')
    p2 = propVarExplained(expr2, tm$mid, ME2, corFnc = 'cor')

    k3 = tibble(mid = tm$mid, kME1 = kME1, kME2 = kME2)
    k31 = k3 %>% group_by(mid) %>%
        summarise(#meanSignAwareKME.qual = mean(sign(kME1)*kME1),
        meanSignAwareKME = abs(mean(sign(kME1)*kME2)),
        cor.kME = abs(cor(kME1, kME2)))
    stopifnot(sum(k31$mid != mids) == 0)

    meanCor.qual = c(); meanCor.pres = c()
    meanAdj.qual = c(); meanAdj.pres = c()
    cor.kIM = c(); cor.cor = c()
    for (mid in mids) {
        dat1 = datExpr1[,tm$gid[tm$mid == mid]]
        dat2 = datExpr2[,tm$gid[tm$mid == mid]]
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

    tibble(mid = mids, 
        #propVarExplained.qual = p1,
        propVarExplained = p2, 
        meanSignAwareKME = k31$meanSignAwareKME,
        #meanCor.qual = meanCor.qual,
        meanCor = meanCor.pres,
        #meanAdj.qual = meanAdj.qual,
        meanAdj = meanAdj.pres,
        cor.kIM = cor.kIM,
        cor.kME = k31$cor.kME,
        cor.cor = cor.cor)
    #}}}
}


