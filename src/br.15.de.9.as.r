require(GenomicRanges)
require(tidyverse)
require(RColorBrewer)
require(pheatmap)
require(DESeq2)

dirw = file.path(Sys.getenv("misc2"), "briggs2")
dirw = '/home/springer/zhoux379/scratch/briggs2/43.deseq'

dirg = '/home/springer/zhoux379/data/genome/Zmays_v4'
f_gtb = file.path(dirg, "51.gtb")
tg = read.table(f_gtb, sep = "\t", header = T, as.is = T)[,c('id','par','cat2','cat3')]

fi = file.path(dirw, '../00.1.read.correct.tsv')
ti = read.table(fi, header = T, sep = "\t", as.is = T)[,1:5]

fc1 = file.path(dirw, "../32.rc.tsv")
tc1 = read.table(fc1, header = T, sep = "\t", as.is = T)
fc2 = file.path(dirw, "../32.rc.as.tsv")
tc2 = read.table(fc2, header = T, sep = "\t", as.is = T)
rownames(tc1) = tc1$gid
rownames(tc2) = tc2$gid

fl = file.path(dirw, "../35.long.tsv")
tl = read.table(fl, sep = "\t", header = T, as.is = T) 

cols = c(brewer.pal(8, 'Dark2'), brewer.pal(9, 'Set1'))

### find DE genes
to = data.frame()
comps = c("B73 vs Mo17", "B73 vs B73xMo17", "Mo17 vs B73xMo17")

for (tissue in unique(ti$Tissue)) {
	#tissue = 'root_0DAP'
	for (comp in comps) {
		gts_comp = unlist(strsplit(comp, split = " vs "))
		tis = ti[ti$Tissue == tissue & ti$Genotype %in% gts_comp,]
		rownames(tis) <- tis$SampleID
		cnts1 = tc1[,tis$SampleID]
		cnts2 = tc2[,tis$SampleID]

## DE sense
dds1 = DESeqDataSetFromMatrix(countData = cnts1, colData = tis, design = ~ Genotype)
#dds1$Genotype <- factor(dds1$Genotype, levels = levels(tis$Genotype))
#dds1 = estimateSizeFactors(dds1)
deseq1 = DESeq(dds1)

#rld <- rlog(deseq1, blind=FALSE)
#fp = sprintf("%s/stats/61.pca.%s.pdf", dirw, tissue)
#pdf(fp, width=6, height=6)
#plotPCA(rld, intgroup = "Genotype", ntop = 10000)
#dev.off()

res1 = as.data.frame(results(deseq1, pAdjustMethod = "fdr"))
is.de1 = ifelse(res1$padj < .05 & res1$log2FoldChange > 1 , "up", ifelse(res1$padj < .05 & res1$log2FoldChange < -1 , "down", 0))

cnts1[rownames(cnts1) %in% rownames(res1)[is.de1 == 'up'],][1:20,]
cnts1[rownames(cnts1) %in% rownames(res1)[is.de1 == 'down'],][1:20,]

		idxs = which(is.de1 %in% c('up','down'))
		tos = data.frame(tissue = tissue, comp = comp, gid = rownames(res1)[idxs], 
			is.de = is.de1[idxs], stringsAsFactors = F)
		to = rbind(to, tos)
	}
}

fo = file.path(dirw, '01.de.tsv')
write.table(to, fo, sep = "\t", row.names = F, col.names = T, quote = F)

stopifnot(1>2)
## Dominance-over-additive value
tis = ti[ti$Genotype %in% c("B73", "Mo17"),]
grp = dplyr::group_by(tis, gid)
tis2 = dplyr::summarise(grp, nts = sum(s.fpkm>=1), nas = sum(a.fpkm>=1))
gids = tis2$gid[tis2$nts>=1 & tis2$nas>=1]

td = ddply(ti[ti$gid %in% gids,], .(gid, Tissue), myfunc <- function(x) {
	s.mid = (x[x[,'Genotype']=="B73",'s.fpkm'] + x[x[,'Genotype']=='Mo17','s.fpkm']) / 2
	s.dif = abs(x[x[,'Genotype']=="B73",'s.fpkm'] - x[x[,'Genotype']=='Mo17','s.fpkm']) / 2
	s.doa = (x[x[,'Genotype']=='B73xMo17','s.fpkm'] - s.mid) / s.dif
	a.mid = (x[x[,'Genotype']=="B73",'a.fpkm'] + x[x[,'Genotype']=='Mo17','a.fpkm']) / 2
	a.dif = abs(x[x[,'Genotype']=="B73",'a.fpkm'] - x[x[,'Genotype']=='Mo17','a.fpkm']) / 2
	a.doa = (x[x[,'Genotype']=='B73xMo17','a.fpkm'] - a.mid) / a.dif
	c('s.doa' = s.doa, 'a.doa' = a.doa)
})
td2 = td[is.finite(td$s.doa) & is.finite(td$a.doa),]
td2$Tissue = factor(td2$Tissue, levels = unique(ti$Tissue))

p1 = ggplot(td2) +
  geom_density_2d(aes(x=s.doa, y=a.doa)) +
  scale_x_continuous(name='sense.DOA', limits = c(-3, 3)) +
  scale_y_continuous(name='antisense.DOA', limits = c(-3, 3)) +
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "lines")) +
  facet_wrap( ~ Tissue, ncol = 5) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.title.y = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "black", angle = 0, hjust = 1)) + 
  geom_rect(xmin=64,xmax=65,ymin=40,ymax=41,fill='black')

fp = file.path(dirw, "stats/51.doa.pdf")
ggsave(p1, filename=fp, width=8, height=7)

## DOA & gene sequence difference
fh = file.path(dirw, 'stats/99.tsv')
th = read.table(fh, sep = "\t", header = T, as.is = T)

th2 = merge(th, ti, by = c("gid", "Tissue"))
th3 = ddply(th2, .(gid, Tissue), myfunc <- function(x) {
	s.mid = (x[x[,'Genotype']=="B73",'s.fpkm'] + x[x[,'Genotype']=='Mo17','s.fpkm']) / 2
	s.dif = abs(x[x[,'Genotype']=="B73",'s.fpkm'] - x[x[,'Genotype']=='Mo17','s.fpkm']) / 2
	s.doa = (x[x[,'Genotype']=='B73xMo17','s.fpkm'] - s.mid) / s.dif
	a.mid = (x[x[,'Genotype']=="B73",'a.fpkm'] + x[x[,'Genotype']=='Mo17','a.fpkm']) / 2
	a.dif = abs(x[x[,'Genotype']=="B73",'a.fpkm'] - x[x[,'Genotype']=='Mo17','a.fpkm']) / 2
	a.doa = (x[x[,'Genotype']=='B73xMo17','a.fpkm'] - a.mid) / a.dif
	c('s.doa' = s.doa, 'a.doa' = a.doa)
})

p1 = ggplot(th3) +
  geom_point(aes(x=s.doa, y=a.doa, color=Tissue)) +
  scale_x_continuous(name='sense.DOA') +
  scale_y_continuous(name='antisense.DOA') +
  scale_color_manual(values = cols) + 
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "lines")) +
#  theme(legend.position = 'right', legend.direction = "horizontal", legend.justification = c(0.3,1), legend.title = element_blank(), legend.key.size = unit(1, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size = 8), legend.background = element_rect(fill=NA, size=0)) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.title.y = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "black", angle = 0, hjust = 1)) + 
  geom_rect(xmin=64,xmax=65,ymin=40,ymax=41,fill='black')

fp = file.path(dirw, "stats/53.de.doa.pdf")
ggsave(p1, filename=fp, width=8, height=6)

dirgm = '/home/springer/zhoux379/data/genome/Mo17'
fd = file.path(dirgm, "62.gene.vnt.tsv")
td = read.table(fd, sep = "\t", header = T, as.is = T)
fm = file.path(dirgm, "21.seq.tsv")
tm = read.table(fm, sep = "\t", header = F, as.is = T)[,1:3]
tm = cbind(tm, gid = sapply(strsplit(tm$V1, split = "_"), "[", 1))
th4 = merge(th3, td, by = 'gid')
t.test(th4$s.doa[th4$nvnt==0], th4$s.doa[th4$nvnt > 0])

## heatmap
fh = file.path(dirw, 'stats/99.tsv')
th = read.table(fh, sep = "\t", header = T, as.is = T)

ths = ddply(th, .(gid), summarise, ntiss = length(Tissue))
table(ths$ntiss)
ths[ths$ntiss > 15,]
gids = ths$gid[ths$ntiss > 1]

fi = file.path(dirw, "39.fpkm.long.tsv")
ti = read.table(fi, sep = "\t", header = T, as.is = T) 

ti1 = ti[ti$gid %in% gids & ti$Genotype %in% c("B73", "Mo17", "B73xMo17"),]
colnames(ti1)[4:5] = c("sense", "antisense")
ti2 = gather(ti1, sense, fpkm, sense, antisense)
ti2$Tissue = factor(ti2$Tissue, levels = unique(ti1$Tissue))
ti2$sense = factor(ti2$sense, levels = c("sense", "antisense"))
ti2 = ti2[order(ti2$sense, ti2$Genotype, ti2$Tissue),]
ti2 = cbind(ti2, cond = sprintf("%s.%s.%s", ti2$sense, ti2$Genotype, ti2$Tissue))
ti3 = ti2[,c(1,5,6)]
ti3$cond = factor(ti3$cond, levels = unique(ti3$cond))
ti4 = spread(ti3, cond, fpkm)
td = ti4
rownames(td) = td$gid
td = td[,-1]
td = asinh(td)

ta = unique(ti2[,c(2:4,6)])
rownames(ta) = ta$cond
ta = ta[,-4]

drows1 <- "correlation"
dcols1 <- "correlation"
col.pal <- brewer.pal(9, "Blues")
col.sense = brewer.pal(8, "Paired")[2:1]
col.geno = brewer.pal(8, "Paired")[6:4]
col.tissue = c(brewer.pal(8, 'Dark2'), brewer.pal(9, 'Set1'))
names(col.tissue) = unique(ta$Tissue)
names(col.geno) = unique(ta$Genotype)
names(col.sense) = unique(ta$sense)
ann_colors = list(
    sense = col.sense,
    Genotype = col.geno,
    Tissue = col.tissue
)

hm.parameters <- list(td, 
  color = col.pal,
  cellwidth = 8, cellheight = 8, scale = "none",
  treeheight_row = 200,
  kmeans_k = NA,
  show_rownames = T, show_colnames = F,
#  main = "Heatmap of asinh(FPKM)",
  clustering_method = "complete",
  cluster_rows = T, cluster_cols = F,
  clustering_distance_rows = drows1, 
  clustering_distance_cols = dcols1,
  annotation_col = ta,
  annotation_colors = ann_colors,
  gaps_col = c(17,34,51,68,85),
  fontsize_row = 8
)
 fo = file.path(dirw, "stats/70.as.heatmap.s.pdf")
do.call("pheatmap", c(hm.parameters, filename=fo))

## single gene plot
plotGeneFPKM <- function(gid, ti, dirw) {
ti1 = ti[ti$gid == gid & ti$Genotype %in% c("B73", "Mo17", "B73xMo17"),]
colnames(ti1)[4:5] = c("sense", "antisense")
ti2 = gather(ti1, sense, fpkm, sense, antisense)
ti2$Tissue = factor(ti2$Tissue, levels = unique(ti1$Tissue))
ti2$sense = factor(ti2$sense, levels = c("sense", "antisense"))
tissues = unique(ti1$Tissue)
tx = data.frame(x = 1:length(tissues), xlab = tissues, stringsAsFactors = F)
ti3 = merge(ti2, tx, by.x = 'Tissue', by.y = 'xlab')

p1 = ggplot(ti3, aes(x = x, y = fpkm, color = Genotype)) +
  geom_point() +
  geom_line() +
  scale_x_continuous(name='Tissue', breaks = tx$x, labels = tx$xlab) +
  scale_y_continuous(name='FPKM') +
  scale_color_brewer(palette = 'Set1') + 
  facet_grid(sense ~ ., scale = 'free') +
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(plot.margin = unit(c(0.5,0.1,0.1,0.1), "lines")) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.title.y = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 30, hjust = 1)) +
  theme(axis.text.y = element_text(size = 8, colour = "black", angle = 0)) + 
  ggtitle(gid)

fp = sprintf("%s/stats/zz.%s.pdf", dirw, gid)
ggsave(p1, filename=fp, width=6, height=5)
}

gids = c("Zm00001d045320",
"Zm00001d013992",
"Zm00001d013639",
"Zm00001d007498",
"Zm00001d025624",
"Zm00001d003713",
"Zm00001d029771",
"Zm00001d043137",
"Zm00001d020761",
"Zm00001d013012"
)
gids = c(
"Zm00001d036755",
"Zm00001d008845")

for (gid in gids) {
	plotGeneFPKM(gid, ti, dirw)
}
