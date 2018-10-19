require(tidyverse)

dirg = '/home/springer/zhoux379/data/genome/Zmays_v4'
dirw = '/home/springer/zhoux379/scratch/briggs'

fm = file.path(dirw, '00.1.read.correct.tsv')
tm = read_tsv(fm)[,1:5]

venn.txt <- function(a, b) {
    ovlp = sum(a %in% b)
    cat(sprintf("%10d: Unique to A\n", length(a)-ovlp))
    cat(sprintf("%10d: Overlap\n", ovlp))
    cat(sprintf("%10d: Unique to B\n", length(b)-ovlp))
}
	
# 
fv = file.path(dirw, "../mo17vnt/52.vnt/Mo17.tsv")
tv = read.table(fv, header = T, sep = "\t", stringsAsFactors = F)

fi = file.path(dirw, "42.ase/21.raw.tsv")
ti = read.table(fi, header = T, sep = "\t", stringsAsFactors = F)

sms = substr(colnames(ti)[15:ncol(ti)], 1, 5)
sts = substr(colnames(ti)[15:ncol(ti)], 6, 7)

sms_b = tm$SampleID[tm$Genotype == 'B73']
idxs_gt_b = 14 + which(sms %in% sms_b & sts == 'GT')
sms_m = tm$SampleID[tm$Genotype == 'Mo17']
idxs_gt_m = 14 + which(sms %in% sms_m & sts == 'GT')

txb = apply(ti[,idxs_gt_b], 1, myfunc <- function(dx) {
	c('btot'=sum(!is.na(dx)), 'bref'=sum(dx==0, na.rm=T))
})
txm = apply(ti[,idxs_gt_m], 1, myfunc <- function(dx) {
	c('mtot'=sum(!is.na(dx)), 'malt'=sum(dx==2, na.rm=T))
})
tx = data.frame(t(rbind(txb, txm)))
tx = cbind(ti[,1:10], tx, pctb = tx$bref/tx$btot, pctm = tx$malt/tx$mtot)

# check overlap with CDS
fg = file.path(dirg, "52.itv.tsv")
tg = read.table(fg, sep = "\t", header = T, as.is = T)
tg = cbind(tg, gbeg = tg$chr*1000000000+tg$beg, gend=tg$chr*1000000000+tg$end)
brks = c(tg$gbeg[1], tg$gend)

tt = tx
gpos = as.integer(tt$chr) * 1000000000 + tt$pos
x = table(cut(gpos, breaks = brks, right = F, labels = F))
nv = rep(0, nrow(tg))
names(nv) = 1:nrow(tg)
nv[names(x)] = as.integer(x)
tg2 = cbind(tg, nv = nv)
to = ddply(tg2, .(type), summarise, nv=sum(nv))
to
#

pos1 = sprintf("%s_%d", tv$chr, tv$pos)
pos2 = sprintf("%s_%d", tx$chr, tx$pos)

tv1 = tv#[tv$PASS==1,]
tx1 = tx[tx$btot > 2 & tx$mtot > 2 & tx$pctb > 0.9 & tx$pctm > 0.9,]
posa = sprintf("%s_%d", tv1$chr, tv1$pos)
posb = sprintf("%s_%d", tx1$chr, tx1$pos)
venn.txt(posa, posb)

tx1c = tx1[posb %in% posa,]
tx1e = tx1[! posb %in% posa,]

p1 = ggplot(tx1e) +
    geom_histogram(aes(x = QD)) +
    theme_bw() +
    theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(0.3,0.1,0.1,0.1), "lines")) +
    theme(legend.position = 'right', legend.direction = "vertical", legend.justification = c(0,0.5), legend.title = element_blank(), legend.key.size = unit(1, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size = 8), legend.background = element_rect(fill=NA, size=0)) +
    theme(axis.title.x = element_text(size = 9)) +
    theme(axis.title.y = element_text(size = 9)) +
    theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
    theme(axis.text.y = element_text(size = 8, colour = "black", angle = 0, hjust = 1))
fp = sprintf("%s/stats/tmp.pdf", dirw)
ggsave(p1, filename = fp, width = 6, height = 6)


