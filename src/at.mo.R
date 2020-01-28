require(tidyverse)
require(GenomicRanges)
require(igraph)

dirg = '/home/springer/zhoux379/data/genome/Zmays_v4'
dirw = '/home/springer/zhoux379/scratch/mo17vnt'
dirw = '/home/springer/zhoux379/scratch/mo17vnt/62.rnaseq.filter'
dirb = '/home/springer/zhoux379/scratch/briggs'

fm = file.path('/home/springer/zhoux379/scratch/briggs', '00.1.read.correct.tsv')
tm = read.table(fm, header = T, sep = "\t", stringsAsFactors = F)[,1:5]

venn.txt <- function(a, b) {
	ovlp = sum(a %in% b)
	cat(sprintf("%10d: Unique to A\n", length(a)-ovlp))
	cat(sprintf("%10d: Overlap\n", ovlp))
	cat(sprintf("%10d: Unique to B\n", length(b)-ovlp))
}

#{{{ filter
#fi = file.path(dirw, "52.vnt/Mo17.filtered.tsv")
fi = file.path(dirw, "52.vnt/Mo17.tsv")
ti = read.table(fi, header = T, sep = "\t", as.is = T)
tf = ti[ti$PASS==1,]
t0 = ti

t1 = t0[(is.na(t0$GT) | t0$GT==2) & (
	(t0$IS_SNP==1 &
		(t0$QD>=2 | is.na(t0$QD)) & 
		(t0$FS<=60 | is.na(t0$FS)) & 
		(t0$MQ>=40 | is.na(t0$MQ)) & 
		(t0$MQRankSum>=-12.5 | is.na(t0$MQRankSum)) & 
		(t0$ReadPosRankSum>=-8 | is.na(t0$ReadPosRankSum)) & 
		(t0$SOR<=4 | is.na(t0$SOR))
	) |
	(t0$IS_SNP==0 & 
		(t0$QD>=2 | is.na(t0$QD)) &
		(t0$FS<=200 | is.na(t0$FS)) & 
		(t0$ReadPosRankSum>=-20 | is.na(t0$ReadPosRankSum)) & 
		(t0$SOR<=10 | is.na(t0$SOR))
	)),]

grv = with(t1, GRanges(seqnames = chr, ranges = IRanges(pos, end = pos+nchar(t1$ref)-1)))
tov = intersect_idx(grv, grv)
tov = tov[tov$idx != tov$qidx,]

edges = as.matrix(tov[,c('idx','qidx')])
edges[,1] = as.character(edges[,1])
edges[,2] = as.character(edges[,2])
graph = graph_from_edgelist(edges)
clus = clusters(graph)
mb = clus$membership
tm = data.frame(idx = as.integer(names(mb)), clu = as.integer(mb))
tms = ddply(tm, .(clu), summarise, clusize = length(clu))
clus_large = tms$clu[tms$clusize>2]
#}}}

#{{{ apply filter 2
dp_mean = mean(ti$DP)
dp_sd = sd(ti$DP)
tf = ti[ti$PASS == 1 & tf$DP < dp_mean+2*dp_sd & tf$GT == 2,]

tf2 = data.frame(chr = tf$chr, beg = tf$pos-1, end = tf$pos+nchar(tf$ref)-1, stringsAsFactors = F)
fo = file.path(dirw, '55.vnt.bed')
write.table(tf2, fo, sep = "\t", row.names = F, col.names = F, quote = F)

fg = file.path(dirg, "51.tbl")
tg = read.table(fg, sep = "\t", header = F, as.is = T)

fg2 = file.path(dirg, "51.gtb")
tg2 = read.table(fg2, sep = "\t", header = T, as.is = T)[,1:2]
colnames(tg2) = c("tid", "gid")

tg3 = tg[tg$V6 == 'cds',][,c(1:3,5)]
colnames(tg3) = c("chr", "beg", "end", "tid")
tg4 = merge(tg3, tg2, by = 'tid')
stopifnot(nrow(tg3) == nrow(tg4))
tg5 = tg4[order(tg4$chr, tg4$beg), c(2:4,1,5)]
tg5$beg = tg5$beg - 1
fo = file.path(dirw, '56.gene.bed')
write.table(tg5, fo, sep = "\t", row.names = F, col.names = F, quote = F)
#intersectBed -wao -a 56.gene.bed -b 55.vnt.bed > 58.gene.vnt.bed

ft = file.path(dirw, "58.gene.vnt.bed")
tt = read.table(ft, sep = "\t", header = F, as.is = T)
colnames(tt) = c('chr', 'beg', 'end', 'tid', 'gid', 'qchr', 'qbeg', 'qend', 'olen')

gp = group_by(tt, tid)
t2 = dplyr::summarise(gp, cnt = sum(qchr != '.'))
t2$cnt[t2$cnt > 5] = 5
table(t2$cnt)

gp = group_by(tt, gid)
t3 = dplyr::summarise(gp, cnt = sum(qchr != '.'))
t3$cnt[t3$cnt > 5] = 5
table(t3$cnt)
#}}}

#{{{ read in genomic intervals
flen = file.path(dirg, "15.sizes")
tlen = read.table(flen, sep = "\t", header = F, as.is = T)
grt = with(tlen, GRanges(seqnames = V1, ranges = IRanges(1, end = V2)))
tt = data.frame(chr = tlen$V1, beg = 1, end = tlen$V2)

fgap = file.path(dirg, "16.gap.bed")
tgap = read.table(fgap, sep = "\t", header = F, as.is = T)
grp = with(tgap, GRanges(seqnames = V1, ranges = IRanges(V2, end = V3)))
tp = data.frame(chr = tgap$V1, beg = tgap$V2, end = tgap$V3)
#}}}

#{{{ create sliding window table using 100kb sliding windows
tt = tt[tt$chr %in% 1:10,]
x = tt$end
names(x) = tt$chr
gr = tileGenome(x, tilewidth = 100000, cut.last.tile.in.chrom = T)

tw = data.frame(chr = as.integer(seqnames(gr)), beg = start(gr), end = end(gr), 
  len = width(gr), stringsAsFactors = F)
tw = cbind(tw, gbeg = tw$chr*1000000000+tw$beg, gend=tw$chr*1000000000+tw$end)
brks = c(tw$gbeg[1], tw$gend)

bp_gap = intersect_basepair(gr, grp)
bp_nogap = tw$len - bp_gap
#}}}

#{{{ plot variant density along chr
tv = ti[ti$chr %in% 1:10,]
tv = tf[tf$chr %in% 1:10,]
tve = tv[tv$GT==1,]
tvo = tv[tv$GT==2,]

gpose = as.integer(tve$chr) * 1000000000 + tve$pos
xe = table(cut(gpose, breaks = brks, right = F, labels = F))
nve = rep(0, nrow(tw))
names(nve) = 1:nrow(tw)
nve[names(xe)] = as.integer(xe)

gposo = as.integer(tvo$chr) * 1000000000 + tvo$pos
xo = table(cut(gposo, breaks = brks, right = F, labels = F))
nvo = rep(0, nrow(tw))
names(nvo) = 1:nrow(tw)
nvo[names(xo)] = as.integer(xo)

nv = nve + nvo
tw2 = cbind(tw, dve = nve/bp_nogap, dvo = nvo/bp_nogap, dv = nv/bp_nogap)

p1 <- ggplot(tw2) +
#  geom_line(aes(x = beg/1000000, y = dv, color = 'All')) +
  geom_line(aes(x = beg/1000000, y = -dve, color = 'Het')) +
  geom_line(aes(x = beg/1000000, y = dvo, color = 'Hom')) +
  facet_grid(chr ~ .) + 
  theme_bw() + 
  scale_x_continuous(name = 'Chr Position (Mbp)', expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0), name = 'Variant Density') +
  scale_color_brewer(palette='Dark2') +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(plot.margin = unit(c(0.5,0.1,0.1,0.1), "lines")) +
  theme(legend.position = c(0.9,0.5), legend.direction = "vertical", legend.justification = c(0.3,1), legend.title = element_blank(), legend.key.size = unit(1, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size = 8), legend.background = element_rect(fill=NA, size=0)) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.title.y = element_text(size = 9, angle = 90)) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 7, colour = "brown", angle = 0, hjust = 1))

fo = file.path(dirw, "stats/60.filtered.pdf")
ggsave(p1, filename=fo, width=8, height=8)
#}}}

#{{{ check for genic SNPs
tv = ti[ti$chr %in% 1:10,]
tv = tf[tf$chr %in% 1:10,]

fg = file.path(dirg, "52.itv.tsv")
tg = read.table(fg, sep = "\t", header = T, as.is = T)
tg = cbind(tg, gbeg = tg$chr*1000000000+tg$beg, gend=tg$chr*1000000000+tg$end)
brks = c(tg$gbeg[1], tg$gend)

gpos = as.integer(tv$chr) * 1000000000 + tv$pos
x = table(cut(gpos, breaks = brks, right = F, labels = F))
nv = rep(0, nrow(tg))
names(nv) = 1:nrow(tg)
nv[names(x)] = as.integer(x)
tg2 = cbind(tg, nv = nv)
to = ddply(tg2, .(type), summarise, len=sum(end-beg+1), nv=sum(nv), vd = nv/len)
to
#}}}

#{{{ check for homo/hetero SNPs in CDSs
tv1 = ti[ti$chr %in% 1:10 & ti$GT == 1,]
tv2 = ti[ti$chr %in% 1:10 & ti$GT == 2,]
tv3 = ti[ti$chr %in% 1:10 & ti$GT == 1 & ti$PASS==1,]
tv4 = ti[ti$chr %in% 1:10 & ti$GT == 2 & ti$PASS==1,]

fg = file.path(dirg, "52.itv.tsv")
tg = read.table(fg, sep = "\t", header = T, as.is = T)
tg = cbind(tg, gbeg = tg$chr*1000000000+tg$beg, gend=tg$chr*1000000000+tg$end)
brks = c(tg$gbeg[1], tg$gend)

tv = tv4
gpos = as.integer(tv$chr) * 1000000000 + tv$pos
x = table(cut(gpos, breaks = brks, right = F, labels = F))
nv = rep(0, nrow(tg))
names(nv) = 1:nrow(tg)
nv[names(x)] = as.integer(x)
tg2 = cbind(tg, nv = nv)
to = ddply(tg2, .(type), summarise, nv=sum(nv))
to
#}}}

#{{{ check TEs - useless
fg2 = file.path(dirg, "51.gtb")
tg2 = read.table(fg2, sep = "\t", header = T, as.is = T)[,c(1:5,16)]
colnames(tg2) = c("id", 'pa', "chr", 'beg', 'end', 'cat')

grc = with(tg[tg$type=='CDS',], GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))
fte = file.path(dirg, "TEs.gff")
tte = read.table(fte, sep = "\t", header = F, as.is = T)[,c(1,4,5)]
gte = with(tte, GRanges(seqnames = V1, ranges = IRanges(V4, end = V5)))
#}}}


#{{{ read RNA-Seq
fi = file.path(dirw, "21.rnaseq.tsv")
ti = read.table(fi, header = T, sep = "\t", stringsAsFactors = F)

sms = substr(colnames(ti)[15:ncol(ti)], 1, 5)
sts = substr(colnames(ti)[15:ncol(ti)], 6, 7)

#dps = apply(ti[,idxs_dp], 1, sum, na.rm = T)

sms_b = tm$SampleID[tm$Genotype == 'B73']
idxs_gt_b = 14 + which(sms %in% sms_b & sts == 'GT')
sms_m = tm$SampleID[tm$Genotype == 'Mo17']
idxs_gt_m = 14 + which(sms %in% sms_m & sts == 'GT')

txb = apply(ti[,idxs_gt_b], 1, myfunc <- function(dx) {
	c('btot'=sum(!is.na(dx)), 'bref'=sum(dx==0, na.rm=T))
})
txm = apply(ti[,idxs_gt_m], 1, myfunc <- function(dx) {
	c('mtot'=sum(!is.na(dx)), 'malt'=sum(dx==2, na.rm=T), 'mhet'=sum(dx==1, na.rm=T))
})
tx = data.frame(t(rbind(txb, txm)))
tx = cbind(ti[,1:10], tx, pctb = tx$bref/tx$btot, pctm = tx$malt/tx$mtot, pcth = tx$mhet/tx$mtot)


tv1 = tv
tx1 = tx#[tx$btot > 1 & tx$mtot > 1 & tx$pctb > 0.8 & (tx$pctm > 0.8),]
posa = sprintf("%s_%d", tv1$chr, tv1$pos)
posb = sprintf("%s_%d", tx1$chr, tx1$pos)
venn.txt(posa, posb)

txc = tx1[posb %in% posa,]
txe = tx1[! posb %in% posa,]
tvc = tv1[posa %in% posb,]
tve = tv1[! posa %in% posb,]

tvd = tv[tv$IS_SNP==0 & nchar(tv$ref) > 1,]
tvd = cbind(tvd, end = tvd$pos + nchar(tvd$ref) - 1)
tvi = tv[tv$IS_SNP==0 & nchar(tv$alt) > 1,]
tvi = cbind(tvi, end = tvi$pos + 1)
grd = with(rbind(tvd, tvi), GRanges(seqnames = chr, ranges = IRanges(pos, end = end)))
#}}}

#{{{ find conflict txc
ib = (txc$btot == 0 | txc$mtot == 0)
sum(ib)
txc = txc[!ib,]

sum(txc$pctb >= 0.8 & txc$pctm >= 0.8)
sum(txc$pctb >= 0.8 & txc$pctm >= 0.8)/nrow(txc)
sum(txc$pctb >= 0.5 & txc$pctm >= 0.5)
sum(txc$pctb >= 0.5 & txc$pctm >= 0.5)/nrow(txc)

idxfs = which(txc$pctb < 0.5 | txc$pctm < 0.5)
tc = txc[idxfs,]
nrow(tc)
grc = with(tc, GRanges(seqnames = chr, ranges = IRanges(pos, end = pos)))
x = intersect_count(grc, grd)
sum(x>0)
tc3 = tc[x == 0,]
nrow(tc3)
tc4 = tc3[tc3$pctb < 0.5 & tc3$pctm > 0.8,]

tc5 = tc3[tc3$pctb > 0.5 & tc3$pctm+tc3$pcth < 0.5,]
#}}}

#{{{ look at txe (variants called by rnaseq but not by genomic-sequencing
sum(txe$mtot == 0)
tx1 = txe[txe$mtot > 0,]
sum(tx1$pctm < 0.5, na.rm = T)
tx2 = tx1[tx1$pctm > 0.5,]
grc = with(tx2, GRanges(seqnames = chr, ranges = IRanges(pos, end = pos)))
x = intersect_count(grc, grd)
sum(x>0)
tx3 = tx2[x == 0,]
nrow(tx3)
#}}}

#{{{ final fitering
nrow(tc3)
fx = file.path(dirw, "30.filter.tsv")
write.table(tc3[,1:2], fx, sep = "\t", row.names = F, col.names = F, quote = F)

nrow(tv)
posf = sprintf("%s_%s", tc3$chr, tc3$pos)
tvf = tv[!posa %in% posf,]
nrow(tvf)
ff = file.path(dirw, '31.mo17.final.tsv')
write.table(tvf, ff, sep = "\t", row.names = F, col.names = T, quote = F)

tob = data.frame(chr=tvf$chr, beg=tvf$pos-1, end=tvf$pos+nchar(tvf$ref)-1, ref=tvf$ref, alt=tvf$alt, qd = tvf$QD, stringsAsFactors = F)
fob = file.path(dirw, "31.mo17.final.bed")
write.table(tob, fob, sep = "\t", row.names = F, col.names = F, quote = F)
# cd ../53.vnt.final
# vcf.filter.py 01.vcf 02.vcf ../62.rnaseq.filter/30.filter.tsv
#}}}

# validation

#{{{ obtain gene IDs for validation
fc = file.path(dirb, '32.cov.tsv')
tc = read.table(fc, header = T, sep = "\t", as.is = T)

npass = apply(tc[,-1], 1, myfunc <- function(x) sum(x>=2))
table(npass)
length(npass >= 10)

gids = tc$gid[npass>=10]

fg = file.path(dirg, "51.tbl")
tg = read.table(fg, sep = "\t", header = F, as.is = T)

fg2 = file.path(dirg, "51.gtb")
tg2 = read.table(fg2, sep = "\t", header = T, as.is = T)[,1:2]
colnames(tg2) = c("tid", "gid")
tids = tg2$tid[tg2$gid %in% gids]

tgs = tg[tg$V5 %in% tids,]
tgs = tgs[tgs$V6 %in% c("cds", "utr5", "utr3"),]

grt = with(tgs, GRanges(seqnames = V1, ranges = IRanges(V2, end = V3)))
grt = reduce(grt)
tl = data.frame(chr = seqnames(grt), beg = start(grt), end = end(grt), stringsAsFactors = F)
#tl$beg = tl$beg - 1
fl = file.path(dirw, '01.loc.tsv')
write.table(tl, fl, sep = "\t", row.names = F, col.names = F, quote = F)

# awk 'BEGIN {OFS="\t"} {$2=$2-1; print}' 01.loc.tsv > 01.loc.bed
# awk 'BEGIN {OFS="\t"} {if ($4=="Intergenic") {$2=$2-1; print}}' 05.itv.tsv > 05.itv.bed
# bcftools view -T 01.loc.tsv 11.mo17.vcf -O v -o 12.vcf
# vcf2tsv.py 12.vcf 12.tsv
# bcftools view -T 01.loc.tsv 21.rnaseq.vcf -O v -o 22.vcf
# vcf2tsv.py 22.vcf 22.tsv
#}}}

#{{{ read in
fi = file.path(dirw, "22.tsv")
ti = read.table(fi, header = T, sep = "\t", stringsAsFactors = F)

sms = substr(colnames(ti)[15:ncol(ti)], 1, 5)
sts = substr(colnames(ti)[15:ncol(ti)], 6, 7)

#dps = apply(ti[,idxs_dp], 1, sum, na.rm = T)

sms_b = tm$SampleID[tm$Genotype == 'B73']
idxs_gt_b = 14 + which(sms %in% sms_b & sts == 'GT')
sms_m = tm$SampleID[tm$Genotype == 'Mo17']
idxs_gt_m = 14 + which(sms %in% sms_m & sts == 'GT')

txb = apply(ti[,idxs_gt_b], 1, myfunc <- function(dx) {
	c('btot'=sum(!is.na(dx)), 'bref'=sum(dx==0, na.rm=T))
})
txm = apply(ti[,idxs_gt_m], 1, myfunc <- function(dx) {
	c('mtot'=sum(!is.na(dx)), 'malt'=sum(dx==2, na.rm=T), 'mhet'=sum(dx==1, na.rm=T))
})
tx = data.frame(t(rbind(txb, txm)))
tx = cbind(ti[,1:10], tx, pctb = tx$bref/tx$btot, pctm = tx$malt/tx$mtot, pcth = tx$mhet/tx$mtot)
#}}}

#{{{ validate sample label mix-up
tis = ti[tx$btot > 1 & tx$mtot > 1 & tx$pctb >= 0.8 & tx$pctm >= 0.8, 14 + which(sts == 'GT')]
colnames(tis) = tm$SampleID
txx = apply(tis, 2, myfunc <- function(dx) {
	c('nsam'=sum(!is.na(dx)), 'nb'=sum(dx==0, na.rm=T), 'nm'=sum(dx==2, na.rm=T))
})
tp = cbind(tm, t(txx))
tp = cbind(tp, pctb = tp$nb/tp$nsam, pctm = tp$nm/tp$nsam)
tp1 = tp[,c(1,4,9,10)]
#tp2 = reshape(tp1, direction='long', varying=list(2:ncol(tp1)), idvar=c('SampleID', "Genotype"), timevar="Location", v.names='Value', times=colnames(tp1)[2:ncol(tp1)])

cols = c(brewer.pal(8, 'Dark2'), brewer.pal(9, 'Set1'))

p1 = ggplot(tp) +
  geom_point(aes(x = pctb, y = pctm, shape = Genotype, color = Tissue)) +
  scale_x_continuous(name = 'Proportion Variants with Homozygous B73 Allele', limits = c(0,1)) +
  scale_y_continuous(name = 'Proportion Variants with Homozygous Mo17 Allele', limits = c(0,1)) +
  scale_color_manual(name = "", values = cols) +
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(plot.margin = unit(c(0.3,0.1,0.1,0.1), "lines")) +
  theme(legend.position = 'right', legend.direction = "vertical", legend.justification = c(0,0.5), legend.title = element_blank(), legend.key.size = unit(1, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size = 8), legend.background = element_rect(fill=NA, size=0)) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.title.y = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "black", angle = 0, hjust = 1))
fp = sprintf("%s/genotype-valid.pdf", dirw)
ggsave(p1, filename = fp, width = 8, height = 7)
#}}}

#{{{ check overlap with CDS
fg = file.path(dirg, "52.itv.tsv")
tg = read.table(fg, sep = "\t", header = T, as.is = T)
tg = cbind(tg, gbeg = tg$chr*1000000000+tg$beg, gend=tg$chr*1000000000+tg$end)
tg = tg[order(tg$gbeg),]
brks = c(tg$gbeg[1], tg$gend)


tt = tx
gpos = as.integer(tt$chr) * 1000000000 + tt$pos
x = table(cut(gpos, breaks = brks, right = T, labels = F))
nv = rep(0, nrow(tg))
names(nv) = 1:nrow(tg)
nv[names(x)] = as.integer(x)
tg2 = cbind(tg, nv = nv)
to = ddply(tg2, .(type), summarise, nv=sum(nv))
to
#tg2[tg2$type=="Intergenic" & tg2$nv > 0,][1:10,]
#

tv1 = tv
tx1 = tx#[tx$btot > 1 & tx$mtot > 1 & tx$pctb > 0.8 & (tx$pctm > 0.8),]
posa = sprintf("%s_%d", tv1$chr, tv1$pos)
posb = sprintf("%s_%d", tx1$chr, tx1$pos)
venn.txt(posa, posb)

txc = tx1[posb %in% posa,]
txe = tx1[! posb %in% posa,]
tvc = tv1[posa %in% posb,]
tve = tv1[! posa %in% posb,]

tvd = tv[tv$IS_SNP==0 & nchar(tv$ref) > 1,]
tvd = cbind(tvd, end = tvd$pos + nchar(tvd$ref) - 1)
tvi = tv[tv$IS_SNP==0 & nchar(tv$alt) > 1,]
tvi = cbind(tvi, end = tvi$pos + 1)
grd = with(rbind(tvd, tvi), GRanges(seqnames = chr, ranges = IRanges(pos, end = end)))
#}}}

#{{{ find conflict txc
ib = (txc$btot == 0 | txc$mtot == 0)
sum(ib)
txc = txc[!ib,]

sum(txc$pctb >= 0.8 & txc$pctm >= 0.8)
sum(txc$pctb >= 0.8 & txc$pctm >= 0.8)/nrow(txc)
sum(txc$pctb >= 0.5 & txc$pctm >= 0.5)
sum(txc$pctb >= 0.5 & txc$pctm >= 0.5)/nrow(txc)

idxfs = which(!(txc$pctb >= 0.5 & txc$pctm >= 0.5))
tc = txc[idxfs,]
nrow(tc)
sum(tc$btot<1)
sum(tc$mtot<1)
tc2 = tc[tc$btot > 0 & tc$mtot > 0,]
grc = with(tc2, GRanges(seqnames = chr, ranges = IRanges(pos, end = pos)))
x = intersect_count(grc, grd)
tc3 = tc2[x == 0,]
nrow(tc3)
tc4 = tc3[tc3$btot > 1 & tc3$pctb < 0.5 & tc3$mtot > 1 & tc3$pctm > 0.8,]

tc5 = tc3[tc3$btot > 1 & tc3$pctb > 0.5 & tc3$mtot > 1 & tc3$pctm+tc3$pcth < 0.5,]
#}}}

#{{{ split txe
#tve[tve$IS_SNP==1,][100:130,]
sum(txe$btot > 1 & txe$mtot > 1 & txe$pctb >= 0.8 & txe$pctm <= 0.2)
txe1 = txe[!(txe$btot > 1 & txe$mtot > 1 & txe$pctb >= 0.8 & txe$pctm <= 0.2),]
tvd = tv[nchar(tv$ref) > 1,]
tvd = cbind(tvd, end = tvd$pos + nchar(tvd$ref) - 1)
gr1 = with(txe1, GRanges(seqnames = chr, ranges = IRanges(pos, end = pos)))
gr2 = with(tvd, GRanges(seqnames = chr, ranges = IRanges(pos, end = end)))
x = intersect_count(gr1, gr2)

plot_hist(txc$QD, file.path(dirw, "solid.rna.ovlp.pdf"), "QD", "count")
plot_hist(tvc$QD, file.path(dirw, "solid.ovlp.pdf"), "QD", "count")
plot_hist(tve$QD, file.path(dirw, "solid.dnaonly.pdf"), "QD", "count")

plot_hist(tve$QD[tve$IS_SNP==0], file.path(dirw, "solid.2.indel.pdf"), "QD", "count")
plot_hist(tve$QD[tve$IS_SNP==1], file.path(dirw, "solid.2.snp.pdf"), "QD", "count")

plot_hist(tve$QD[tve$GT==1], file.path(dirw, "solid.3.het.pdf"), "QD", "count")
plot_hist(tve$QD[tve$GT!=1], file.path(dirw, "solid.3.nhet.pdf"), "QD", "count")


plot_hist(txc$pctb, file.path(dirw, "solid.6.pctb.pdf"), "percent_b73", "count")
plot_hist(txc$pctm, file.path(dirw, "solid.6.pctm.pdf"), "percent_mo17", "count")
plot_hist(txe$pctb, file.path(dirw, "solid.6.e.pctb.pdf"), "percent_b73", "count")
plot_hist(txe$pctm, file.path(dirw, "solid.6.e.pctm.pdf"), "percent_mo17", "count")

x = sum(txc$btot > 1 & txc$mtot > 1 & txc$pctb >= 0.8 & txc$pctm >= 0.8)
cat(sprintf("%g (%.02f%%) out of %g\n", x, x/nrow(txc)*100, nrow(txc)))
x = sum(txc$btot > 1 & txc$mtot > 1 & txc$pctb >= 0.5 & txc$pctm >= 0.5)
cat(sprintf("%g (%.02f%%) out of %g\n", x, x/nrow(txc)*100, nrow(txc)))
x = sum(txc$btot > 0 & txc$mtot > 0 & txc$pctb >= 0.5 & txc$pctm >= 0.5)
cat(sprintf("%g (%.02f%%) out of %g\n", x, x/nrow(txc)*100, nrow(txc)))
#}}}


