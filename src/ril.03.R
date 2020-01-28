source("functions.R")

#{{{ plot RIL heights
dirw = file.path(dird, "81_ril_heights")
fi = file.path(dirw, 'RIL_heights.tsv')
ti = read_tsv(fi) %>%
    rename(gt=1, ht_B=2, ht_M=3, ht=4) %>%
    filter(!is.na(ht) & !is.na(ht_B) & !is.na(ht_M))

tp = ti %>% mutate(dht_B = ht_B-ht, dht_M = ht_M-ht) %>%
    mutate(dht_B = pmax(dht_B, 0), dht_M = pmax(dht_M, 0)) %>%
    mutate(dht = dht_B + dht_M, prop_dht_B = dht_B / dht) %>%
    replace_na(list(prop_dht_B = .5))

fh = '~/projects/rnaseq/data/11_qc/rn13a/meta.tsv'
th = read_tsv(fh)
gts = unique(th$Genotype)

tp %>% filter(gt %in% gts, ht >= 192, ht <= 208) %>%
    arrange(prop_dht_B, dht) %>%
    print(n=50)

tps = tp %>% filter(gt %in% c("M0317",'M0021','M0016','M0014'))

jit = position_jitter(width=10, height=10)
p = ggplot(tp) +
    geom_point(aes(ht, dht, color=prop_dht_B), size=2) +
    geom_text_repel(data=tps,aes(ht, dht, label=gt), nudge_x=10) +
    scale_x_continuous(name='self_height') +
    scale_y_continuous(name='diff(B_OC, self) + diff(M_OC, self)') +
    #scale_color_viridis(limits=c(0,1), option='magma') +
    scale_color_gradientn(name='prop. diff(B_OC, self)', colors = brewer.pal(12, "RdYlBu"), limits=c(0,1)) +
    otheme(xtitle=T,ytitle=T,xtext=T,ytext=T,xtick=T, ytick=T,
        legend.pos='top.right', legend.title=T)
fo = file.path(dirw, '10.pdf')
ggsave(fo, p, width=8, height=8)

p = ggplot(tp) +
    geom_point(aes(dht_B, dht_M, color=ht), position=jit, size=3) +
    scale_x_continuous(name='B_OC - self') +
    scale_y_continuous(name='M_OC - self') +
    scale_color_viridis(name='self height', option='viridis') +
    otheme(xtitle=T,ytitle=T,xtext=T,ytext=T,xtick=T, ytick=T,
        legend.pos='top.left', legend.title=T)
fo = file.path(dirw, '11.pdf')
ggsave(fo, p, width=8, height=8)
#}}}

#{{{ plot and output genotype
dirw = file.path(dird, '82_ril_genotype')
fi = '/home/springer/zhoux379/projects/rnaseq/data/11_qc/rn13a/meta.tsv'
th = read_tsv(fi) %>% mutate(lab = sprintf("%s_%s", Genotype, Replicate))

#{{{ plot genotype
fw = '~/projects/genome/data/Zmays_B73/15_intervals/20.win11.tsv'
tw = read_tsv(fw)
offs = c(0, cumsum(tw$size)[-nrow(tw)]) + 0:10 * 10e6
tx = tw %>% mutate(off = offs) %>%
    mutate(gstart=start+off, gend=end+off, gpos=(gstart+gend)/2) %>%
    select(rid,chrom,gstart,gend,gpos,off) %>% filter(chrom!='B99')
fi = '/home/springer/zhoux379/projects/rnaseq/data/raw/rn13a/ril.rds'
res = readRDS(fi)

ti = res$cp
tps = th %>% select(sid=SampleID,Genotype,Replicate) %>% arrange(Genotype) %>%
    mutate(y = 1:n())
tp = ti %>% inner_join(tx, by='rid') %>%
    mutate(start=start+off, end=end+off) %>%
    inner_join(tps, by='sid')
tz = ti %>% mutate(size=end-start) %>% group_by(sid, gt) %>%
    summarise(size = sum(size)) %>%
    mutate(total_size = sum(size)) %>%
    mutate(prop = size / total_size) %>% ungroup() %>%
    select(sid, gt, prop) %>% spread(gt, prop) %>%
    replace_na(list(a=0,b=0,h=0)) %>%
    inner_join(tps, by='sid') %>% arrange(Genotype) %>%
    filter(Replicate==1) %>% select(-y)

xmax = max(tp$end)
gts = c("B73","Mo17","M0021","M0317","M0014","M0016")
tcol = pal_startrek()(1)
ty = tps %>% group_by(Genotype) %>% summarise(y=mean(y)) %>% ungroup() %>%
    mutate(col=ifelse(Genotype %in% gts, tcol, 'black')) %>%
    inner_join(tz, by='Genotype') %>% mutate(lab=sprintf("%.1f", a*100))
p = ggplot(tp) +
    geom_rect(aes(xmin=start,xmax=end,ymin=y-.3,ymax=y+.4, fill=gt)) +
    geom_text(data=ty, aes(x=xmax+5e6,y=y, label=lab), hjust=0, size=2.5, color=ty$col) +
    scale_x_continuous(breaks=tx$gpos, labels=tx$chrom, expand=expand_scale(mult=c(.001,.03))) +
    scale_y_continuous(breaks=ty$y, labels=ty$Genotype, expand=expand_scale(mult=c(.001,.001))) +
    scale_fill_manual(values=pal_simpsons()(8)[c(1,2,5)], labels=c('B73','Mo17','het')) +
    otheme(legend.pos='top.center.out', legend.dir='h',
           xtext=T, xtick=T, ytext=T, ytick=T) +
    theme(axis.text.y = element_text(color=ty$col))
fo = file.path(dirw, '01.ril.gt.pdf')
ggsave(fo, p, width=10, height=10)
#}}}

gts = c("M0317",'M0021','M0016','M0014')
to = tp %>% filter(Genotype %in% gts) %>%
    mutate(sid=str_c(Genotype, Replicate, sep='_')) %>%
    select(sid, chrom,start,end, gt)

fo = file.path(dirw, '04.selected.gt.tsv')
write_tsv(to, fo)
#}}}


