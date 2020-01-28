source("br.fun.r")
sid = 'me99b'
sid = 'me99b.m'
dirw = file.path(dirp, ifelse(sid == 'me99b.m', '41_qc_m', "41_qc"))
genome = ifelse(sid == 'me99b.m', 'Mo17', 'B73')
x = load(file.path(dirg, genome, '55.rda'))
fm = file.path(dirw, '10.rda')
y = load(fm)

#{{{ CPM distri.
to1 = tmm %>% select(Tissue, Genotype, gid, CPM) %>%
    filter(CPM > 0) %>%
    mutate(lCPM = log2(CPM+1)) %>%
    mutate(lCPM = ifelse(lCPM > 10, 10, lCPM))
to2 = to1 %>% 
    mutate(x_bin = cut(lCPM, breaks=seq(0,10,0.2), include.lowest = T, 
        labels = seq(0.1, 9.9, by = 0.2)))
sum(is.na(to2$x_bin))
to2$x_bin = as.numeric(as.character(to2$x_bin))
to3 = to2 %>%
    group_by(Tissue, Genotype, x_bin) %>%
    summarise(num_genes = n(), total_lcpm = sum(lCPM))
to4 = gather(to3, type, value, -Tissue, -Genotype, -x_bin)

labs = c(0, 1, 10, 100, 1000, 10000)
brks = log2(labs+1)
p1 = ggplot(to3) +
    geom_bar(aes(x = x_bin, y = num_genes, fill = Genotype), stat = 'identity', position = 'dodge', alpha=.9) +
    #geom_line(aes(color = Genotype, linetype = Genotype)) +
    #geom_point(aes(color = Genotype, shape = Genotype), size = 0.5) +
    scale_x_continuous(name = 'CPM', limits = c(0, 10), breaks = brks, labels = labs) +
    scale_y_continuous(name = 'Num. Genes') +
    scale_fill_aaas() +
    facet_wrap(~Tissue, ncol = 3, scale = 'free') +
    theme_bw() +
    theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(.3,.3,.3,.3), "lines")) +
    theme(legend.position = c(.85,.03), legend.direction = "vertical", legend.justification = c(.5,0)) +
    theme(legend.title = element_blank(), legend.key.size = unit(1, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size = 8)) 
fp = sprintf("%s/12.cpm.den.pdf", dirw)
ggsave(p1, filename = fp, width = 8, height = 12)
#
p2 = ggplot(to3) +
    geom_bar(aes(x = x_bin, y = total_lcpm, fill = Genotype), stat = 'identity', position = 'dodge', alpha=.7) +
    scale_x_continuous(name = 'CPM', limits = c(0, 10), breaks = brks, labels = labs) +
    scale_y_continuous(name = 'Num. Genes') +
    scale_fill_aaas() +
    facet_wrap( ~Tissue, ncol = 3, scale = 'free') +
    theme_bw() +
    theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(.3,.3,.3,.3), "lines")) +
    theme(legend.position = c(.85,0.03), legend.direction = "vertical", legend.justification = c(.5,0)) +
    theme(legend.title = element_blank(), legend.key.size = unit(1, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size = 8))
fp = sprintf("%s/12.cpm.cumsum.pdf", dirw)
ggsave(p2, filename = fp, width = 8, height = 12)

g1 <- ggplotGrob(p1)
g2 <- ggplotGrob(p2)
gs = list(g1, g2)
wids = c(4,4)
g <- gtable_matrix(name = 'demo', grobs = matrix(gs, ncol = length(gs)), widths = wids, heights = 1)

fp = sprintf("%s/12.cpm.den.pdf", dirw)
pdf(file = fp, width = 10, height = 30, bg = 'transparent')
grid.newpage()
grid.draw(g)
dev.off()
#}}}

#{{{ num expressed genes
etags = c('silent','lowly-expressed','expressed')
tp = tmm %>%
    mutate(etag = ifelse(CPM == 0, 'silent', 
                  ifelse(CPM < 1, 'lowly-expressed', 'expressed'))) %>%
    group_by(Tissue, Genotype, etag) %>%
    summarise(ngene = n()) %>%
    ungroup() %>%
    mutate(x = sprintf("%s:%s", Tissue, Genotype)) %>%
    mutate(Tissue = factor(Tissue, levels = tissues23)) %>%
    mutate(etag = factor(etag, levels = rev(etags))) %>%
    mutate(Genotype = factor(Genotype, levels=gts))

p1 = ggplot(tp) +
    geom_bar(aes(x = Genotype, y = ngene/1000, fill = etag), stat = 'identity', width = 1, alpha = 1) +
    coord_polar('y') +
    scale_fill_d3() +
    facet_wrap(Tissue ~ ., ncol = 4) +
    otheme(xticks = T, yticks = T, xgrid = T, ygrid = T,
           xtext = T, ytext = T) +
    theme(strip.background = element_blank(), strip.placement = "outside") +
    theme(strip.text = element_text(size = 8)) +
    theme(legend.position = c(.5,1), legend.justification = c(.5,-.5)) +
    guides(direction = 'horizontal', fill = guide_legend(nrow = 1)) +
    theme(legend.title = element_blank(), legend.key.size = unit(.8, 'lines'), legend.text = element_text(size = 8)) +
    theme(plot.margin = unit(c(1.5,.5,.5,0), "lines")) 
fo = file.path(dirw, "15.expressed.pdf")
ggsave(p1, filename = fo, width = 6, height = 8)

#{{{ #old bar plot
p1 = ggplot(tp) +
    geom_bar(aes(x = Genotype, y = ngene, fill = etag), position = position_stack(reverse = T), stat = 'identity', width = .8, alpha = .7) +
    scale_x_discrete(expand = c(.01,0)) +
    scale_y_continuous(name = 'Num. Genes', expand = c(0,0)) +
    coord_flip() +
    scale_fill_d3() +
    facet_grid(Tissue ~ ., switch = 'y') +
    otheme(xtitle = T, ytitle = F, yticks = T, xgrid = F, ygrid = F,
           xtext = T, ytext = T) +
    theme(strip.background = element_blank(), strip.placement = "outside") +
    theme(strip.text.y = element_text(size = 8, angle = 180, hjust=1)) +
    theme(legend.position = c(.5,1), legend.justification = c(.5,0)) +
    guides(direction = 'horizontal', fill = guide_legend(nrow = 1)) +
    theme(legend.title = element_blank(), legend.key.size = unit(0.8, 'lines'), legend.text = element_text(size = 8)) +
    theme(panel.border = element_blank(), panel.grid = element_blank()) +
    theme(plot.margin = unit(c(2,.5,.5,0), "lines"))
#}}}
#}}}

#{{{ tissue-specificity of top100 expresed genes
tz = t_rep_merged %>% group_by(Tissue, Genotype) %>%
    top_n(100, CPM)

ty = tz %>%
    group_by(Genotype, gid) %>% summarise(ntis = n()) %>%
    count(Genotype, ntis) %>% group_by(Genotype) %>%
    summarise(prop = sum(n[ntis<5])/sum(n))
ty = tz %>%
    group_by(Tissue, gid) %>% summarise(ngt = n()) %>%
    count(Tissue, ngt) %>% group_by(Tissue) %>% 
    summarise(prop = sum(n[ngt>2])/sum(n)) %>% print(n=23)

#}}}

#{{{# look at top 50 gene proportion
summary(trl$fpm)
p1 = ggplot(trl) +
    geom_boxplot(aes(x = sid, y = fpm), outlier.shape = NA) + #, draw_quantiles = c(0.25, 0.5, 0.75)) + 
    coord_flip() +
    scale_x_discrete(name = '', breaks = tm$SampleID, labels = sprintf("%s|%s", tm$Tissue, tm$Genotype)) +
    scale_y_continuous(name = 'FPM', limits = c(0, 30)) +
    theme_bw() +
    theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "lines")) +
    theme(axis.title.x = element_text(size = 9)) +
    theme(axis.title.y = element_text(size = 9)) +
    theme(axis.text.x = element_text(size = 8, color = "black", angle = 0)) +
    theme(axis.text.y = element_text(size = 8, color = "black", angle = 0, hjust = 1))
fp = sprintf("%s/11.fpm.pdf", diro)
ggsave(p1, filename = fp, width = 6, height = 12)

y = apply(tr[,2:ncol(tr)], 2, myfun <- function(x) {
x1 = sort(x, decreasing = T)
c('top[01-05]' = sum(x1[1:5]),
    'top[06-10]' = sum(x1[6:10]),
    'top[11-15]' = sum(x1[11:15]),
    'top[16-20]' = sum(x1[16:20]),
    'top[21-30]' = sum(x1[21:30]),
    'top[31-40]' = sum(x1[31:40]),
    'top[41-50]' = sum(x1[41:50])
)
})
y = cbind(tag = rownames(y), as.data.frame(y))
yl = reshape(y, direction = 'long', varying = list(2:ncol(y)), idvar = c("tag"), timevar = "sid", v.names = 'fpm', times = colnames(y)[2:ncol(y)])

yl$tag = factor(yl$tag, levels = rev(sort(unique(yl$tag))))
p1 = ggplot(yl) +
    geom_bar(aes(x = sid, y = fpm/1000000, fill = tag), position = 'stack', stat = 'identity', width = 0.7) +
    coord_flip() +
    scale_x_discrete(name = '', breaks = tm$SampleID, labels = sprintf("%s|%s", tm$Tissue, tm$Genotype)) +
    scale_y_continuous(name = 'Proportion Reads', expand = c(0,0), limits = c(0, 1)) +
    scale_fill_brewer(palette = "Accent") +
    theme_bw() +
    theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(1,1,0.1,0.1), "lines")) +
    theme(legend.position = c(0.7, 0.7), legend.direction = "vertical", legend.justification = c(0,0), legend.title = element_blank(), legend.key.size = unit(1, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size = 8), legend.background = element_rect(fill=NA, size=0)) +
    theme(axis.title.x = element_text(size = 9)) +
    theme(axis.title.y = element_text(size = 9)) +
    theme(axis.text.x = element_text(size = 8, color = "black", angle = 0)) +
    theme(axis.text.y = element_text(size = 8, color = "black", angle = 0, hjust = 1))
fp = sprintf("%s/11.fpm.top50.pdf", diro)
ggsave(p1, filename = fp, width = 6, height = 12)
#}}}
#{{{# get top10 expressed genes w. functions
fg = file.path(dirg, "51.gtb")
tg = read.table(fg, sep = "\t", header = T, as.is = T)[,c(1:6,16:18)]
gb = group_by(tg, par)
tg2 = summarise(gb, fam = names(sort(table(cat3), decreasing = T))[1])

z = ddply(trl, .(sid), myfun <- function(x) {
	x1 = x[order(x[,'fpm'], decreasing = T),][1:10,]
	x2 = cbind.data.frame(x1, crpm = as.numeric(cumsum(x1[,'fpm'])))
	x2
})
z2 = merge(z, tg2, by.x = 'gid', by.y = 'par')
z3 = merge(z2, tm[,c(1,3)], by.x = 'sid', by.y = 'SampleID')
z4 = z3[order(z3$sid, -z3$fpm),]

fo = file.path(diro, '11.fpm.top10.tsv')
write.table(z4, fo, sep = "\t", row.names = F, col.names = T, quote = F)

## FPM correction
trw = spread(trl4, sid, fpm)
yt = ddply(trl4, .(sid), summarise, fpm = sum(fpm))
y = apply(trw[,2:ncol(trw)], 2, myfun <- function(x) {
	x1 = sort(x, decreasing = T)
	c(
		'top0' = 0,
		'top5' = sum(x1[1:5]),
		'top10' = sum(x1[1:10]),
		'top20' = sum(x1[1:20]),
		'top30' = sum(x1[1:30]),
		'top50' = sum(x1[1:50]),
		'top100' = sum(x1[1:100])
	)
})
y = cbind.data.frame(opt = rownames(y), y, stringsAsFactors = F)
#yl2 = reshape(y, direction = 'long', varying = list(2:ncol(y)), idvar = c("opt"), timevar = "sid", v.names = 'fpm', times = colnames(y)[2:ncol(y)])
yl = gather(y, sid, fpm, -opt)
#}}}
#{{{# check housekeeping genes
fh = '/home/springer/zhoux379/data/genome/Zmays_v4/housekeeping/11.tsv'
thk = read.table(fh, sep = "\t", header = T, as.is = T)
hgids = thk$ngid
th = trl4[trl4$gid %in% hgids,]
th$gid = factor(th$gid, levels = thk$ngid)

opts = c("top0", "top5", "top10", "top20", "top30", "top50", "top100")
opt = opts[6]
z = yl[yl$opt == opt,]
z2 = merge(z, yt, by = 'sid')
z3 = cbind(z2[,c('sid','opt')], sf = z2$fpm.y / (z2$fpm.y - z2$fpm.x))
th2 = merge(th, z3, by = 'sid')
th3 = cbind(th2, fpmc = th2$fpm * th2$sf)

cols = brewer.pal(nrow(thk), 'Set1')
p1 = ggplot(th3) +
    geom_bar(aes(x = sid, y = fpmc, fill = gid), position = 'stack', stat = 'identity', width = 0.7) +
    coord_flip() +
    scale_x_discrete(name = '') +
    scale_y_continuous(name = 'FPM', expand = c(0,0)) +
    scale_fill_manual(values = cols, breaks = thk$ngid, labels = thk$name) +
    theme_bw() +
    theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(0.1,2,0.1,0.1), "lines")) +
    theme(legend.position = 'top', legend.direction = "horizontal", legend.justification = c(0,0), legend.title = element_blank(), legend.key.size = unit(1, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size = 8), legend.background = element_rect(fill=NA, size=0)) +
    guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
    theme(axis.title.x = element_text(size = 9)) +
    theme(axis.title.y = element_text(size = 9)) +
    theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
    theme(axis.text.y = element_text(size = 8, colour = "brown", angle = 0, hjust = 1))
fp = sprintf("%s/13.hk.%s.pdf", diro, opt)
ggsave(p1, filename = fp, width = 6, height = 8)

tb = data.frame()
for (opt in opts) {
	z = yl[yl$opt == opt,]
	z2 = merge(z, yt, by = 'sid')
	z3 = cbind(z2[,c('sid','opt')], sf = z2$fpm.y / (z2$fpm.y - z2$fpm.x))
	th2 = merge(th, z3, by = 'sid')
	th3 = cbind(th2, fpmc = th2$fpm * th2$sf)
	tb = rbind(tb, th3[,c('sid','gid','opt','fpmc')])
}
tb2 = ddply(tb, .(opt, gid), summarize, cv = (sd(fpmc)/mean(fpmc))*100)
tb2$opt = factor(tb2$opt, levels = opts)

cols = brewer.pal(length(opts), 'Paired')
p1 = ggplot(tb2) +
    geom_bar(aes(x = gid, y = cv, fill = opt), position = 'dodge', stat = 'identity', width = 0.8) +
    coord_flip() +
    scale_x_discrete(name = '', breaks = thk$ngid, labels = thk$name) +
    scale_y_continuous(name = 'C.V. of FPM across tissues') +
    scale_fill_manual(values = cols) +
    theme_bw() +
    theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(0.1,2,0.1,0.1), "lines")) +
    theme(legend.position = 'top', legend.direction = "horizontal", legend.justification = c(0,0), legend.title = element_blank(), legend.key.size = unit(1, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size = 8), legend.background = element_rect(fill=NA, size=0)) +
    guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
    theme(axis.title.x = element_text(size = 9)) +
    theme(axis.title.y = element_text(size = 9)) +
    theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
    theme(axis.text.y = element_text(size = 8, colour = "brown", angle = 0, hjust = 1))
fp = sprintf("%s/14.hk.pdf", diro)
ggsave(p1, filename = fp, width = 7, height = 7)
#}}}

