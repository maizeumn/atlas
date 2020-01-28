require(plyr)
require(dplyr)
require(GenomicRanges)
require(ggplot2)
require(seqinr)
require(stringr)
source('Location.R')

dirg = '/home/springer/zhoux379/data/genome/Zmays_v4'
dirw = '/home/springer/zhoux379/scratch/mo17vnt/53.vnt.final'

### create final variant BED: mo17vnt/53.vnt.final/01.bed

### make cds interval BED
fg = file.path(dirg, "51.tbl")
tg = read.table(fg, sep = "\t", header = F, as.is = T)
colnames(tg) = c("chr", 'beg', 'end', 'srd', 'id', 'type', 'cat')
unique(tg$chr)

tgc = tg[tg$type == 'cds',]
to = data.frame(chr=as.integer(tgc$chr), beg=tgc$beg-1, end=tgc$end, srd=tgc$srd, tid=tgc$id, stringsAsFactors = F)
to = to[order(to$chr, to$beg),]
to = cbind(to, idx = 1:nrow(to))
fo = file.path(dirw, "05.cds.bed")
#write.table(to, fo, sep = "\t", row.names = F, col.names = F, quote = F)

# bedseq.pl -db $genome/Zmays_v4/11_genome.fas -i 05.cds.bed -o 06.cds.seq.bed
# intersectBed -wao -a 05.cds.bed -b 01.bed > 09.cds.vnt.bed

fi = file.path(dirw, "09.cds.vnt.bed")
ti = read.table(fi, sep = "\t", header = F, as.is = T)
colnames(ti) = c("chr", "beg", "end", "srd", "tid", "idx", "vchr", "vbeg", "vend", "ref", "alt", 'qd', "olen")
ti$beg = ti$beg + 1
ti$vbeg = ti$vbeg + 1
ti2 = ti[ti$vchr != '.',]

grp = dplyr::group_by(ti2, idx)
tv = dplyr::summarise(grp, chr=chr[1], beg=beg[1], end=end[1], srd=srd[1], tid=tid[1], 
	vnt = paste(sprintf("%s:%d:%d:%s:%s:%s", vchr, vbeg, vend, ref, alt, qd), collapse=" "))

fc = file.path(dirw, "06.cds.seq.bed")
tc = read.table(fc, sep = "\t", header = F, as.is = T)
colnames(tc) = c("chr", "beg", "end", "srd", "tid", "idx", 'seq')

tc2 = merge(tc, tv[,c('idx','vnt')], by = 'idx', all.x = T)
stopifnot(tc$idx == tc2$idx)

fo = file.path(dirw, "11.tsv")
tc2$beg = tc2$beg + 1
#write.table(tc2, fo, sep = "\t", row.names = F, col.names = F, quote = F, na = '')

# run python script to recover mo17 cds sequences
# vnt.recover.py 11.tsv 13.new.tsv

### stats
nvnt = sapply(strsplit(tv$vnt, split=" "), myfunc <- function(y) sum(!is.na(y)))
tv2 = cbind(tv, nvnt=nvnt)
tv3 = merge(ti, tv2[,c('idx','nvnt')], by = 'idx', all.x = T)
stopifnot(nrow(ti) == nrow(tv3))

grp = dplyr::group_by(tv3, tid)
tv4 = dplyr::summarise(grp, nvnt = sum(nvnt, na.rm = T))
cat(sprintf("nvnt = 0: %d\n", sum(tv4$nvnt == 0)))
cat(sprintf("nvnt = 1: %d\n", sum(tv4$nvnt == 1)))
cat(sprintf("nvnt>= 2: %d\n", sum(tv4$nvnt >= 2)))

fl = file.path(dirg, "57.longest.tsv")
tl = read.table(fl, sep = "\t", header = F, as.is = T)

tv5 = merge(tv4, tl, by.x = 'tid', by.y = 'V2')
stopifnot(nrow(tv5) == nrow(tl))
tv5 = tv5[,3:2]
colnames(tv5)[1] = 'gid'
cat(sprintf("nvnt = 0: %d\n", sum(tv5$nvnt == 0)))
cat(sprintf("nvnt = 1: %d\n", sum(tv5$nvnt == 1)))
cat(sprintf("nvnt>= 2: %d\n", sum(tv5$nvnt >= 2)))

dirgm = '/home/springer/zhoux379/data/genome/Mo17'
fo1 = file.path(dirgm, "61.transcript.vnt.tsv")
#write.table(tv4, fo1, sep = "\t", row.names = F, col.names = T, quote = F, na = '')
fo2 = file.path(dirgm, "62.gene.vnt.tsv")
#write.table(tv5, fo2, sep = "\t", row.names = F, col.names = T, quote = F, na = '')

### construct transcript seqs from cds seqs
translate_srd <- function(dx) {
	seqstr = s2c(dx['seq']); srd = dx['srd']
	sense = ifelse(srd == '-', 'R', 'F')
	c2s(translate(seqstr, sens = sense))
}

fi = file.path(dirw, "11.tsv")
ti = read.table(fi, sep = "\t", header = F, as.is = T)
colnames(ti) = c("idx", "chr", "beg", "end", "srd", "tid", 'seq', 'vntstr')
grp = dplyr::group_by(ti, tid)
ti2 = dplyr::summarise(grp, srd = srd[1], seq = paste(seq[order(beg)], collapse=''))

#prostrs = apply(ti2, 1, translate_srd)
#stop_count = str_counts("*", prostrs)

fm = file.path(dirw, "13.new.tsv")
tm = read.table(fm, sep = "\t", header = F, as.is = T)
colnames(tm) = c("idx", "mseq")
tm = merge(tm, ti[,c('idx','tid','beg','srd')], by = 'idx')
grp = dplyr::group_by(tm, tid)
tm2 = dplyr::summarise(grp, srd = srd[1], seq = paste(mseq[order(beg)], collapse=''))

stopifnot(tm2$tid == ti2$tid)
to = cbind(ti2, mseq = tm2$seq)
fo = file.path(dirw, "15.seq.tsv")
write.table(to, fo, sep = "\t", row.names = F, col.names = F, quote = F)

### run mo17.translate.py 15.seq.tsv 21.seq.tsv 22.mo17.fasta
fi = file.path(dirw, "21.seq.tsv")
ti = read.table(fi, sep = "\t", header = F, as.is = T)
colnames(ti) = c("tid", 'bstop', 'mstop', 'bseq', 'mseq')

table(ti$bstop)
n1 = sum(ti$bseq == ti$mseq)
n2 = sum(ti$bseq != ti$mseq)
n3 = sum(ti$mstop > 1)
cat(sprintf("no AA change: %d\nat least 1 AA change: %d\n  Premature stop: %d\n", n1, n2, n3))

fl = file.path(dirg, "57.longest.tsv")
tl = read.table(fl, sep = "\t", header = F, as.is = T)

ti = ti[ti$tid %in% tl$V2,]
table(ti$bstop)
n1 = sum(ti$bseq == ti$mseq)
n2 = sum(ti$bseq != ti$mseq)
n3 = sum(ti$mstop > 1)
cat(sprintf("no AA change: %d\nat least 1 AA change: %d\n  Premature stop: %d\n", n1, n2, n3))
