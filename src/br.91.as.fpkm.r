require(dplyr)
require(GenomicRanges)
require(tidyr)
require(ape)
require(ggplot2)
require(plyr)

dirw = file.path(Sys.getenv("misc2"), "briggs2")
dirw = '/home/springer/zhoux379/scratch/briggs2'
diro = '/home/springer/zhoux379/scratch/briggs2/45.antisense'

dirg = '/home/springer/zhoux379/data/genome/Zmays_v4'
f_gtb = file.path(dirg, "51.gtb")
tg = read.table(f_gtb, sep = "\t", header = T, as.is = T)[,c('id','par','cat2','cat3')]

fi = file.path(dirw, '00.1.read.correct.tsv')
ti = read.table(fi, header = T, sep = "\t", as.is = T)[,1:5]

### compute FPM
fi1 = file.path(dirw, '32.rc.tsv')
fi2 = file.path(dirw, '32.rc.as.tsv')

ti1 = read.table(fi1, header = T, sep = "\t", as.is = T)
ti2 = read.table(fi2, header = T, sep = "\t", as.is = T)

total_reads1 = apply(ti1[,-1], 2, sum)
total_reads2 = apply(ti2[,-1], 2, sum)

to1 = ti1
to2 = ti2
for (i in 1:length(total_reads1)) {
	to1[,i+1] = to1[,i+1] / total_reads1[i] * 1000000
	to2[,i+1] = to2[,i+1] / total_reads2[i] * 1000000
}
fo1 = file.path(dirw, '33.fpm.tsv')
write.table(to1, fo1, sep = "\t", row.names = F, col.names = T, quote = F)
fo2 = file.path(dirw, '33.fpm.as.tsv')
write.table(to2, fo2, sep = "\t", row.names = F, col.names = T, quote = F)

### compute FPKM
dirg = '/home/springer/zhoux379/data/genome/Zmays_v4'
f_gtb = file.path(dirg, "51.gtb")
f_tbl = file.path(dirg, "51.tbl")
tg = read.table(f_gtb, sep = "\t", header = T, as.is = T)[,1:2]
tt = read.table(f_tbl, sep = "\t", header = F, as.is = T)
colnames(tg) = c("tid", "gid")
colnames(tt) = c("chr", "beg", "end", "srd", "tid", "type", "fam")
tt2 = tt[tt$type %in% c('cds', 'utr5', 'utr3'),]
tg2 = merge(tg, tt2, by = 'tid')

gr = with(tg2, GRanges(seqnames = chr, ranges = IRanges(beg, end = end), gid = gid))
x = unlist(reduce(split(gr, elementMetadata(gr)$gid)))
tr = data.frame(gid = names(x), chr = seqnames(x), beg = start(x), end = end(x), stringsAsFactors = F)
grp = dplyr::group_by(tr, gid)
tr2 = dplyr::summarise(grp, len = sum(end - beg + 1))

to1 = merge(ti1, tr2, by = 'gid')
stopifnot(nrow(to1) == nrow(ti1))
for (i in 2:(ncol(to1)-1)) {
	to1[,i] = to1[,i] / (to1$len/1000)
}
dim(to1)
fo = file.path(dirw, '34.fpkm.tsv')
write.table(to1[,-ncol(to1)], fo, sep = "\t", row.names = F, col.names = T, quote = F)

to2 = merge(ti2, tr2, by = 'gid')
stopifnot(nrow(to2) == nrow(ti2))
for (i in 2:(ncol(to2)-1)) {
	to2[,i] = to2[,i] / (to2$len/1000)
}
dim(to2)
fo = file.path(dirw, '34.fpkm.as.tsv')
write.table(to2[,-ncol(to2)], fo, sep = "\t", row.names = F, col.names = T, quote = F)

### get overlapping gene models
source("Location.R")
dirg = '/home/springer/zhoux379/data/genome/Zmays_v4'
f_gtb = file.path(dirg, "51.gtb")
tg = read.table(f_gtb, sep = "\t", header = T, as.is = T)[,c('id','par','chr','beg','end','srd','cat2','cat3')]

grp = dplyr::group_by(tg, par)
tg2 = dplyr::summarise(grp, chr = chr[1], beg = min(beg), end = max(end), srd = srd[1])
grg = with(tg2, GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))
tr = intersect_idx(grg, grg)
tr = tr[tr$idx != tr$qidx & tr$idx < tr$qidx,]

tg3 = cbind(tg2[,c(1,5)], idx = 1:nrow(tg2), len = tg2$end - tg2$beg + 1)
colnames(tg3) = c("gid1", "srd1", "idx", "len1")
tr2 = merge(tr, tg3, by = 'idx')
colnames(tg3) = c("gid2", "srd2", "qidx", "len2")
tr3 = merge(tr2, tg3, by = 'qidx')
tr3 = tr3[,-c(1:2)]

tr4 = tr3[tr3$olen / tr3$len1 >= 0.1 | tr3$olen / tr3$len2 >= 0.1,]


fi = file.path(dirw, "35.long.tsv")
ti = read.table(fi, sep = "\t", header = T, as.is = T)
tib = spread(ti[ti$Genotype == 'B73',c(1,2,4)], Tissue, fpm)
tim = spread(ti[ti$Genotype == 'Mo17',c(1,2,4)], Tissue, fpm)
nb = apply(tib[,-1], 1, myfunc <- function(x) sum(x >=1 ))
nm = apply(tim[,-1], 1, myfunc <- function(x) sum(x >=1 ))
gids = tib$gid[ nb >= 1 | nm >= 1]

tr5 = tr4[tr4$gid1 %in% gids & tr4$gid2 %in% gids,]
tt = apply(tr5, 1, pcc <- function(x) {
	e1 = as.numeric(c(tib[tib$gid == x['gid1'], -1], tim[tim$gid == x['gid1'], -1]))
	e2 = as.numeric(c(tib[tib$gid == x['gid2'], -1], tim[tim$gid == x['gid2'], -1]))
	y = cor.test(e1, e2)
	c(pval = y$p.value, pcc = cor(e1, e2))
})
tr5 = cbind(tr5, t(tt))
tr6 = tr5[which(tr5$pval < 0.05), -c(1:3,11)]
tr6 = cbind(i = 1:nrow(tr6), tr6)
tr61 = tr6[tr6$srd1 != tr6$srd2,]
tr62 = tr6[tr6$srd1 == tr6$srd2,]

## gene pair plot
plotGenePairFPKM <- function(gid1, gid2, ti, diro) {
ti1 = ti[ti$gid %in% c(gid1, gid2) & ti$Genotype %in% c("B73", "Mo17", "B73xMo17"),]
ti2 = ti1
ti2$Tissue = factor(ti2$Tissue, levels = unique(ti1$Tissue))
tissues = unique(ti1$Tissue)
tx = data.frame(x = 1:length(tissues), xlab = tissues, stringsAsFactors = F)
ti3 = merge(ti2, tx, by.x = 'Tissue', by.y = 'xlab')

p1 = ggplot(ti3, aes(x = x, y = fpkm, color = Genotype)) +
  geom_point() +
  geom_line() +
  scale_x_continuous(name='Tissue', breaks = tx$x, labels = tx$xlab) +
  scale_y_continuous(name='FPKM') +
  scale_color_brewer(palette = 'Set1') + 
  facet_grid(gid ~ ., scale = 'free') +
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(plot.margin = unit(c(0.5,0.1,0.1,0.1), "lines")) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.title.y = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 30, hjust = 1)) +
  theme(axis.text.y = element_text(size = 8, colour = "black", angle = 0))

fp = sprintf("%s/zz.%s.pdf", diro, gid1)
ggsave(p1, filename=fp, width=6, height=5)
}

for (i in c(70,88,99)) {
	gid1 = tr6$gid1[i]; gid2 = tr6$gid2[i]
	plotGenePairFPKM(gid1, gid2, ti, diro)
}
for (i in c(3,34,36,92)) {
	gid1 = tr6$gid1[i]; gid2 = tr6$gid2[i]
	plotGenePairFPKM(gid1, gid2, ti, diro)
}
