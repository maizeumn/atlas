require(plyr)
require(dplyr)
require(GenomicRanges)
require(ggplot2)
require(RColorBrewer)
source('Location.R')

dirg = '/home/springer/zhoux379/data/genome/Zmays_v4'
dirw = '/home/springer/zhoux379/scratch/mo17vnt/61.rnaseq.solid'
dirb = '/home/springer/zhoux379/scratch/briggs'

###
fm = file.path('/home/springer/zhoux379/scratch/briggs', '00.1.read.correct.tsv')
tm = read.table(fm, header = T, sep = "\t", stringsAsFactors = F)[,1:5]

venn.txt <- function(a, b) {
	ovlp = sum(a %in% b)
	cat(sprintf("%10d: Unique to A\n", length(a)-ovlp))
	cat(sprintf("%10d: Overlap\n", ovlp))
	cat(sprintf("%10d: Unique to B\n", length(b)-ovlp))
}

### obtain gene IDs for validation
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


	
### 
fv = file.path(dirw, "12.tsv")
t0 = read.table(fv, header = T, sep = "\t", stringsAsFactors = F)
tv = t0[(is.na(t0$GT) | t0$GT==2) & (
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
nrow(t0)
nrow(tv)

### 
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

### validate sample label mix-up
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



# check overlap with CDS
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

### find conflict txc
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

### split txe
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



