require(plyr)
require(dplyr)
require(GenomicRanges)
require(ggplot2)
source('Location.R')
require(igraph)

dirg = '/home/springer/zhoux379/data/genome/Zmays_v4'
dirw = '/home/springer/zhoux379/scratch/mo17vnt'
#fi = file.path(dirw, "52.vnt/Mo17.filtered.tsv")
fi = file.path(dirw, "52.vnt/Mo17.tsv")
ti = read.table(fi, header = T, sep = "\t", as.is = T)
tf = ti[ti$PASS==1,]

### filter ###
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

### apply filter 2
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

### read in genomic intervals
flen = file.path(dirg, "15.sizes")
tlen = read.table(flen, sep = "\t", header = F, as.is = T)
grt = with(tlen, GRanges(seqnames = V1, ranges = IRanges(1, end = V2)))
tt = data.frame(chr = tlen$V1, beg = 1, end = tlen$V2)

fgap = file.path(dirg, "16.gap.bed")
tgap = read.table(fgap, sep = "\t", header = F, as.is = T)
grp = with(tgap, GRanges(seqnames = V1, ranges = IRanges(V2, end = V3)))
tp = data.frame(chr = tgap$V1, beg = tgap$V2, end = tgap$V3)

##### create sliding window table using 100kb sliding windows
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

### plot variant density along chr
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

### check for genic SNPs
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

### check for homo/hetero SNPs in CDSs
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

### check TEs - useless
fg2 = file.path(dirg, "51.gtb")
tg2 = read.table(fg2, sep = "\t", header = T, as.is = T)[,c(1:5,16)]
colnames(tg2) = c("id", 'pa', "chr", 'beg', 'end', 'cat')

grc = with(tg[tg$type=='CDS',], GRanges(seqnames = chr, ranges = IRanges(beg, end = end)))
fte = file.path(dirg, "TEs.gff")
tte = read.table(fte, sep = "\t", header = F, as.is = T)[,c(1,4,5)]
gte = with(tte, GRanges(seqnames = V1, ranges = IRanges(V4, end = V5)))


