require(plyr)
require(dplyr)
require(GenomicRanges)
require(ggplot2)
require(RColorBrewer)
source('Location.R')

dirg = '/home/springer/zhoux379/data/genome/Zmays_v4'
dirw = '/home/springer/zhoux379/scratch/mo17vnt/62.rnaseq.filter'
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
	
### 
fv = file.path(dirw, "11.mo17.tsv")
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

### find conflict txc
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

### look at txe (variants called by rnaseq but not by genomic-sequencing
sum(txe$mtot == 0)
tx1 = txe[txe$mtot > 0,]
sum(tx1$pctm < 0.5, na.rm = T)
tx2 = tx1[tx1$pctm > 0.5,]
grc = with(tx2, GRanges(seqnames = chr, ranges = IRanges(pos, end = pos)))
x = intersect_count(grc, grd)
sum(x>0)
tx3 = tx2[x == 0,]
nrow(tx3)


### final fitering
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
