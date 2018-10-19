source("br.fun.R")

dirw = file.path(Sys.getenv("misc2"), "briggs", "61.perm")
fm = file.path(dirw, "../59.allmodules.tsv")
tm = read.table(fm, header = T, sep = "\t", as.is = T, quote = '')
tms = unique(tm[,c(1,3)])
require(stringr)
m2f = word(tms$funcat, 1, 2)
m2f[is.na(m2f)] = word(tms$funcat[is.na(m2f)], 1)
names(m2f) = tms$mid

## pre-process of permutation output
pc = data.frame()
for (i in 0:1000) {
	fi = sprintf("%s/01.perm/%04d.rda", dirw, i)
	if(!file.exists(fi)) next
	x = load(fi)
	
	pp = data.frame()
	for (gt1 in gts) {
	for (gt2 in gts) {
		if(gt1 == gt2) next
		pcc.expr1 = as.numeric(comp[[gt1]][[gt2]]$pcc.expr)
		pcc.coex1 = as.numeric(comp[[gt1]][[gt2]]$pcc.coex)
		pc1 = data.frame(perm = i, gt1 = gt1, gt2 = gt2,
            pcc.expr = pcc.expr1, pcc.coex = pcc.coex1, stringsAsFactors = F)
        pc = rbind(pc, pc1)

        pp1 = comp[[gt1]][[gt2]]$stat
        #pp1 = pp1[substr(pp1$mid, 0, 4) == 'Corn',]
        pp1 = cbind(perm = i, gt1 = gt1, gt2 = gt2, pp1[,c('mid', 'propVarExplained', 'meanSignAwareKME.pres', 'meanCor', 'meanAdj', 'cor.kIM', 'cor.kME', 'cor.cor')])
        pp = rbind(pp, pp1)
	}
	}
	colheader = ifelse(i == 0, TRUE, FALSE)
	fo = sprintf("%s/03.permout/%04d.tsv", dirw, i)
	write.table(pp, fo, sep = "\t", row.names = F, col.names = colheader, quote = F)
	cat(i, "\n")
}
fp = sprintf("%s/05.pcc.rda", dirw)
save(pc, file = fp)
# cat 03.permout/*.tsv > 04.mp.tsv


## calculate z-score for MP statistics
fp = sprintf("%s/04.mp.tsv", dirw)
tp = read.table(fp, header = T, sep = "\t", as.is = T)

grp = dplyr::group_by(tp[tp$perm > 0,], mid)
tp1 = dplyr::summarise(grp, 
	propVarExplained.mean = mean(propVarExplained), 
	propVarExplained.sd = sd(propVarExplained),
	meanSignAwareKME.mean = mean(meanSignAwareKME.pres), 
	meanSignAwareKME.sd = sd(meanSignAwareKME.pres),
	meanCor.mean = mean(meanCor), meanCor.sd = sd(meanCor),
	meanAdj.mean = mean(meanAdj), meanAdj.sd = sd(meanAdj),
	cor.kIM.mean = mean(cor.kIM), cor.kIM.sd = sd(cor.kIM),
	cor.kME.mean = mean(cor.kME), cor.kME.sd = sd(cor.kME),
	cor.cor.mean = mean(cor.cor), cor.cor.sd = sd(cor.cor)
)
tp2 = merge(tp[tp$perm==0, -1], tp1, by = 'mid')
tp3 = within(tp2, {
	Z.propVarExplained = (propVarExplained - propVarExplained.mean) / propVarExplained.sd
	Z.meanSignAwareKME = (meanSignAwareKME.pres - meanSignAwareKME.mean)/meanSignAwareKME.sd
	Z.meanCor = (meanCor - meanCor.mean) / meanCor.sd
	Z.meanAdj = (meanAdj - meanAdj.mean) / meanAdj.sd
	Z.cor.kIM = (cor.kIM - cor.kIM.mean) / cor.kIM.sd
	Z.cor.kME = (cor.kME - cor.kME.mean) / cor.kME.sd
	Z.cor.cor = (cor.cor - cor.cor.mean) / cor.cor.sd
	}
)

fp = sprintf("%s/05.pcc.rda", dirw)
x = load(fp)

pp = tp3
fo = sprintf("%s/09.rda", dirw)
save(pp, pc, file = fo)


### look at pcc.expr and pcc.coex distribution
fp = sprintf("%s/09.rda", dirw)
x = load(fp)

tp = pc
tp$gt1[tp$gt1 == 'B73xMo17'] = 'F1'
tp$gt2[tp$gt2 == 'B73xMo17'] = 'F1'
tp = cbind(tp, comp = sprintf("%s:%s", tp$gt1, tp$gt2))
tp = tp[tp$comp %in% c("B73:Mo17", "F1:B73", "F1:Mo17"),]
p1 = ggplot(tp[tp$perm > 0,]) +
  geom_histogram(aes(x = pcc.expr), bins = 60) +
  geom_point(data = tp[tp$perm == 0,], aes(x = pcc.expr, y = 1, shape=comp, color=comp), size = 3) +
  scale_x_continuous(name = 'Correlation of asinh(FPKM)') +
  scale_y_continuous(name = 'Frequency') +
  scale_shape_discrete(name = "") +
  scale_color_brewer(name = "", palette = "Set1") +
  #facet_grid(gt1 ~ gt2) +
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines")) +
  theme(legend.position = c(0.7, 0.8), legend.direction = "vertical", legend.justification = c(0,0.5), legend.title = element_blank(), legend.key.size = unit(1, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size = 10), legend.background = element_rect(fill=NA, size=0)) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0, hjust = 1)) +
  theme(axis.text.y = element_text(size = 8, colour = "black", angle = 0))
fp = sprintf("%s/11.pcc.expr.pdf", dirw)
ggsave(p1, filename = fp, width = 6, height = 6)

p1 = ggplot(tp[tp$perm > 0,]) +
  geom_histogram(aes(x = pcc.coex), bins = 60) +
  geom_point(data = tp[tp$perm == 0,], aes(x = pcc.coex, y = 1, shape=comp, color=comp), size = 3) +
  scale_x_continuous(name = 'Correlation of co-expression matrix') +
  scale_y_continuous(name = 'Frequency') +
  scale_shape_discrete(name = "") +
  scale_color_brewer(name = "", palette = "Set1") +
  #facet_grid(gt1 ~ gt2) +
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines")) +
  theme(legend.position = c(0.1, 0.8), legend.direction = "vertical", legend.justification = c(0,0.5), legend.title = element_blank(), legend.key.size = unit(1, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size = 10), legend.background = element_rect(fill=NA, size=0)) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0, hjust = 1)) +
  theme(axis.text.y = element_text(size = 8, colour = "black", angle = 0))
fp = sprintf("%s/11.pcc.coex.pdf", dirw)
ggsave(p1, filename = fp, width = 6, height = 6)


### look at pathways
fp = sprintf("%s/09.rda", dirw)
x = load(fp)
pp = pp[,c(1:3,25:ncol(pp))]

require(pheatmap)
td = cor(pp[,-c(1:3)])
pheatmap(td, 
  #color = col.pal,
  cellwidth = 30, cellheight = 30, scale = "none",
  treeheight_row = 200,
  kmeans_k = NA,
  show_rownames = T, show_colnames = F,
#  main = "Heatmap of asinh(FPKM)",
  clustering_method = "complete",
  cluster_rows = T, cluster_cols = T,
  #clustering_distance_rows = drows1, 
  #clustering_distance_cols = dcols1,
  #annotation_col = ta,
  #annotation_colors = ann_colors,
  fontsize_row = 12,
  #width = 6, height = 6,
  filename = file.path(dirw, "10.cor.mpstats.pdf")
)

pp1 = pp
pp1$gt1[pp1$gt1 == 'B73xMo17'] = 'F1'
pp1$gt2[pp1$gt2 == 'B73xMo17'] = 'F1'
pp1 = cbind(pp1, comp = sprintf("%s..%s", pp1$gt1, pp1$gt2))
pp2 = pp1[pp1$comp %in% c("B73..Mo17", "F1..B73", "F1..Mo17"),]

pp3 = pp2#spread(pp2[,c('mid','comp','Z.cor.cor')], comp, Z.cor.cor)
res = strsplit(pp3$mid, split = "[.]")
res1 = sapply(res, "[", 1)
res2 = sapply(res, "[", 2)
idxs = which(res1 == "MCL")
res1[idxs] = sprintf("%s.%s", res1[idxs], res2[idxs])
pp3 = cbind(pp3, opt = res1)
pp4 = pp3[!pp3$opt %in% c("CornCyc", "GO"),]

z_cutoff = -3

p1 = ggplot(pp4) +
  geom_point(aes(x = Z.cor.cor, y = mid, color = comp, shape = comp)) +
  geom_vline(xintercept = z_cutoff, color = 'mediumorchid') + 
  scale_x_continuous() +
  #scale_y_continuous() +
  scale_shape_discrete(name = "") +
  scale_color_brewer(name = "", palette = "Set1") +
  facet_wrap(~opt, scales = 'free') +
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines")) +
  theme(legend.position = 'top', legend.direction = "horizontal", legend.justification = c(0.5,0.5), legend.title = element_blank(), legend.key.size = unit(1, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size = 10), legend.background = element_rect(fill=NA, size=0)) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.title.y = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0, hjust = 1)) +
  theme(axis.text.y = element_blank())
fp = sprintf("%s/22.z1.pdf", dirw)
ggsave(p1, filename = fp, width = 10, height = 12)

mids = pp4$mid[(pp4$comp == 'F1..B73' & pp4$Z.cor.cor <= z_cutoff) | (pp4$comp == 'F1..Mo17' & pp4$Z.cor.cor <= z_cutoff)]
pp5 = pp4[pp4$mid %in% mids,]
p1 = ggplot(pp5) +
  geom_point(aes(x = Z.cor.cor, y = mid, color = comp, shape = comp)) +
  geom_vline(xintercept = z_cutoff, color = 'mediumorchid') + 
  scale_x_continuous() +
  scale_y_discrete(breaks=unique(pp5$mid), labels=m2f[unique(pp5$mid)]) +
  scale_shape_discrete(name = "") +
  scale_color_brewer(name = "", palette = "Set1") +
  facet_wrap(~opt, scales = 'free') +
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines")) +
  theme(legend.position = 'top', legend.direction = "horizontal", legend.justification = c(0.5,0.5), legend.title = element_blank(), legend.key.size = unit(1, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size = 10), legend.background = element_rect(fill=NA, size=0)) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0, hjust = 1)) +
  theme(axis.text.y = element_text(size = 8))
fp = sprintf("%s/22.z2.pdf", dirw)
ggsave(p1, filename = fp, width = 10, height = 12)

unique(tm[tm$mid %in% mids,c('mid','funcat')])