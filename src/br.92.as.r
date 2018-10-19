require(dplyr)
require(GenomicRanges)
require(tidyr)
require(ape)
require(ggplot2)
require(plyr)

dirw = file.path(Sys.getenv("misc2"), "briggs2")
dirw = '/home/springer/zhoux379/scratch/briggs2'

dirg = '/home/springer/zhoux379/data/genome/Zmays_v4'
f_gtb = file.path(dirg, "51.gtb")
tg = read.table(f_gtb, sep = "\t", header = T, as.is = T)[,c('id','par','cat2','cat3')]

ft = file.path(dirg, "57.longest.tsv")
tt = read.table(ft, sep = "\t", header = F, as.is = T)
tt = tg[tg$id %in% tt$V2,]

fi = file.path(dirw, '00.1.read.correct.tsv')
ti = read.table(fi, header = T, sep = "\t", as.is = T)[,1:5]

### stats
fi1 = file.path(dirw, '34.fpkm.tsv')
fi2 = file.path(dirw, '34.fpkm.as.tsv')

ti1 = read.table(fi1, header = T, sep = "\t", as.is = T)
ti2 = read.table(fi2, header = T, sep = "\t", as.is = T)
e1 = ti1[,-1]
e2 = ti2[,-1]

tl = data.frame()
tw1 = data.frame(gid = ti1$gid, stringsAsFactors = F)
tw2 = data.frame(gid = ti1$gid, stringsAsFactors = F)
tiss = unique(ti$Tissue); gts = unique(ti$Genotype)
for (tis in tiss) {
	for (gt in gts) {
		sms = ti$SampleID[ti$Tissue == tis & ti$Genotype == gt]
		if(length(sms) == 0) next
		x1 = ti1[,sms]
		x2 = ti2[,sms]
		s.fpkm = apply(x1, 1, mean)
		a.fpkm = apply(x2, 1, mean)
		t1 = data.frame(gid = ti1$gid, Tissue = tis, Genotype = gt, s.fpkm = s.fpkm, a.fpkm = a.fpkm, stringsAsFactors = F)
		tl = rbind(tl, t1)
		tw1 = cbind(tw1, fpkm = s.fpkm)
		colnames(tw1)[ncol(tw1)] = sprintf("%s.%s", tis, gt)
		tw2 = cbind(tw2, fpkm = a.fpkm)
		colnames(tw2)[ncol(tw2)] = sprintf("%s.%s", tis, gt)
	}
}

#tl1 = gather(ti1, SampleID, s.fpkm, -gid)
#tl2 = gather(ti2, SampleID, a.fpkm, -gid)
#tl3 = rbind(cbind(tl1, type='sense'), cbind(tl2, type='antisense'))
#tl4 = merge(tl3, ti[,c(1,3:5)], by='SampleID')
#grp = dplyr::group_by(tl4, SampleID, gid, type, Tissue, Genotype)
#tl = dplyr::summarise(grp, fpkm = median(fpkm))

fl = file.path(dirw, '39.fpkm.long.tsv')
write.table(tl, fl, sep = "\t", row.names = F, col.names = T, quote = F)
fw1 = file.path(dirw, '38.fpkm.s.wide.tsv')
write.table(tw1, fw1, sep = "\t", row.names = F, col.names = T, quote = F)
fw2 = file.path(dirw, '38.fpkm.a.wide.tsv')
write.table(tw2, fw2, sep = "\t", row.names = F, col.names = T, quote = F)


## prop of antisense transcription
sum1 = apply(ti1[,-1], 2, sum)
sum2 = apply(ti2[,-1], 2, sum)
stopifnot(names(sum1) == ti$SampleID)
tp = cbind(ti[,1:5], x = 1:nrow(ti), pct_as = sum2/(sum1+sum2))
tpx = ddply(tp, .(Tissue), summarise, x = mean(x))
p1 = ggplot(tp) +
  geom_bar(aes(x=x, y=pct_as, fill=Genotype), stat='identity') +
  scale_x_continuous(name='Tissue', breaks = tpx$x, labels = tpx$Tissue, expand = c(0,0), limits = c(0,max(tp$x)+1)) +
  scale_y_continuous(name='Proportion of AntiSense Transcription', expand = c(0.01,0)) +
  coord_flip() +
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "lines")) +
  theme(legend.position = 'top', legend.direction = "horizontal", legend.justification = c(0.3,1), legend.title = element_blank(), legend.key.size = unit(1, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size = 8), legend.background = element_rect(fill=NA, size=0)) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "black", angle = 0, hjust = 1)) + 
  geom_rect(xmin=64,xmax=65,ymin=40,ymax=41,fill='black')

fp = file.path(dirw, "stats/41.as.prop.pdf")
ggsave(p1, filename=fp, width=5, height=8)


##### read in long expression matrix
fl = file.path(dirw, '39.fpkm.long.tsv')
tl = read.table(fl, header = T, sep = "\t", as.is = T)
#tl$s.fpkm = asinh(tl$s.fpkm)
#tl$a.fpkm = asinh(tl$a.fpkm)

## look at S/AS correlation
tl0 = tl[tl$Genotype == 'B73' & tl$a.fpkm > 0,]
tl0 = cbind(tl0, s.no = tl0$s.fpkm == 0)
p1 = ggplot(tl0) +
  geom_boxplot(aes(x=s.no, y=a.fpkm), notch = T, outlier.shape = NA) +
  scale_x_discrete(name='No expression (sense transcript)') +
  scale_y_continuous(name='asinh(antisense FPKM)', expand = c(0.01,0), limits = c(0, 3)) +
  facet_wrap( ~ Tissue, ncol=3) +
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "lines")) +
  theme(legend.position = 'top', legend.direction = "horizontal", legend.justification = c(0.3,1), legend.title = element_blank(), legend.key.size = unit(1, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size = 8), legend.background = element_rect(fill=NA, size=0)) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.title.y = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "black", angle = 0, hjust = 1)) + 
  geom_rect(xmin=64,xmax=65,ymin=40,ymax=41,fill='black')

fp = file.path(dirw, "stats/42.as.ttest.pdf")
ggsave(p1, filename=fp, width=5, height=8)


t.test(tl0$a.fpkm[tl0$s.fpkm==0], tl0$a.fpkm[tl0$s.fpkm>0])
cor.test(tl0$s.fpkm, tl0$a.fpkm) ## -0.08 pval < 2.2E-16

### look at S/AS ratio distribution
p1 = ggplot(tl0) +
  geom_histogram(aes(x=asinh(s.fpkm/a.fpkm))) +
  scale_x_discrete(name='asinh(sense FPKM / antisense FPKM)') +
  scale_y_continuous(name='Frequency') +
  facet_wrap( ~ Tissue, ncol=3) +
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "lines")) +
  theme(legend.position = 'top', legend.direction = "horizontal", legend.justification = c(0.3,1), legend.title = element_blank(), legend.key.size = unit(1, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size = 8), legend.background = element_rect(fill=NA, size=0)) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.title.y = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "black", angle = 0, hjust = 1)) + 
  geom_rect(xmin=64,xmax=65,ymin=40,ymax=41,fill='black')

fp = file.path(dirw, "stats/43.as.ratio.pdf")
ggsave(p1, filename=fp, width=5, height=8)



## hclust + heatmap
tl1 = tl[tl$Genotype %in% c('B73', "Mo17"),]

tl21 = tl1[,c('gid','Tissue','Genotype','rpm1')] %>% spread(Genotype, rpm1)
tl22 = tl1[,c('gid','Tissue','Genotype','rpm2')] %>% spread(Genotype, rpm2)
colnames(tl22)[3:ncol(tl22)] = paste(colnames(tl22)[3:ncol(tl22)], '_as', sep = '')
tl2 = cbind(tl21, tl22[,-c(1:2)])

tl3 = tl2[(tl2$B73_as > 0 | tl2$Mo17_as > 0) & tl2$Tissue == 'sheath_v12',]
tl4 = tl3[(tl3$B73_as == 0 & tl3$Mo17_as > 1) | (tl3$B73_as > 1 & tl3$Mo17_as == 0),]
#tl5 = ddply(tl4, .(gid), summarise, ntis = length(Tissue))
#tl6 = merge(tl5, tt[,-1], by.x = 'gid', by.y = 'par')
#tl6[tl6$ntis >=2,][1:100,]

e = tl4[,-c(1:2)]
cor_opt = "pearson"
hc_opt = "ward.D"
e.r.dist <- as.dist(1-cor(t(e), method = cor_opt))
e.r.hc <- hclust(e.r.dist, method = hc_opt)
hc = e.r.hc
idxs = hc$order

tl5 = gather(tl4[,-2], type, rpm, -gid)
tl5$gid = factor(tl5$gid, levels = tl5$gid[idxs])

pb <- ggplot(tl5) +
  geom_tile(aes(x = type, y = gid, fill = rpm), height = 1) +
  theme_bw() + 
  scale_x_discrete(name = '') +
  scale_y_discrete(expand = c(0, 0), name = '') +
#  scale_fill_gradient(name = 'RPM', space = "Lab", low = 'firebrick1', high = 'dodgerblue', na.value = 'grey50') +
#  theme(legend.position = 'top', legend.direction = "horizontal", legend.justification = c(0, 1), legend.title = element_text(size = 8), legend.key.size = unit(0.5, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size = 7), legend.background = element_rect(fill=NA, size=0)) +
  theme(plot.margin = unit(c(0.5,0.5,0,1), "lines")) +
  theme(axis.title.x = element_blank(), axis.ticks.length = unit(0, 'lines')) +
  theme(axis.text.x = element_text(colour = "black", size = 8, angle = 60, hjust = 1)) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.y = element_blank()) +
  theme(axis.line = element_line(size = 0.3, colour = "grey", linetype = "solid"))
fo = sprintf("%s/stats/42.heatmap.pdf", dirw)
ggsave(pb, filename = fo, width = 8, height = 10)
