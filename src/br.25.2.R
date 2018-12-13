#{{{ head
source("functions.R")
source(file.path(dirr, "enrich.R")
sid = 'me99b'
genome = 'B73'
x = load(file.path(dirg, genome, '55.rda'))
fi = file.path(dird, '49_coop', "01.master.rda")
x = load(fi)
# color palettes
cols.ts = pal_npg()(3)
cols.expr = pal_npg()(4)
cols.mix3 = pal_d3()(10)[c(1,2,7)]
cols.mix4 = pal_d3()(10)[c(8,1,2,7)]
cols.dom = pal_startrek()(7)
cols.reg1 = pal_d3()(5)
cols.reg2 = pal_d3()(10)[c(8,6,9)]
#}}}

#{{{ remake tm w.o. using DE 2-FC cutoff
fi = file.path(dirp, "41.qc/10.rep.merged.RData")
x = load(fi)
nrow(t_rep_merged)/69
fi = file.path(dirp, "42.de/11.de.RData")
x = load(fi)
fi = file.path(dirp, "44.ase/10.ase.2.RData")
x = load(fi)
#fi = file.path(dirp, "45.doa/10.dom.RData")
#x = load(fi)
t_de = t_de %>% mutate(Tissue = factor(Tissue, levels = tissues23))
t_ase = t_ase %>% mutate(Tissue = factor(Tissue, levels = tissues23))
#t_dom = t_dom %>% mutate(Tissue = factor(Tissue, levels = tissues23)) %>% select(Tissue, gid, Dom, DoA)
taglst = list(
    pDE = c("DE_B", "DE_M", "non_DE"),
    hDE = levels(t_de$hDE),
    Dom = levels(t_de$Dom),
    Reg1 = levels(t_ase$Reg1),
    Reg2 = levels(t_ase$Reg2)
)

tm = t_rep_merged %>% filter(Genotype != 'MxB') %>%
    select(Tissue, gid, Genotype, CPM) %>%
    spread(Genotype, CPM)
nrow(tm)/23
tm2 = tm %>%
    left_join(t_de, by = c("Tissue", 'gid')) %>%
    #mutate(pDE = pDE2) %>%  #!!! use 2-FC DE cutoff
    mutate(silent = ifelse(is.na(pDE), 1, 0)) %>%
    select(Tissue, gid, B73, Mo17, BxM, silent, log2MB, pDE, log2FM, hDE, Dom, DoA)
tm2 %>% filter(is.nan(DoA) | is.infinite(DoA))
tm2 %>% count(pDE, Dom) %>% spread(Dom, n)
tm2 %>% count(pDE, hDE) %>% spread(hDE, n)

tm3 = tm2 %>% left_join(t_ase, by = c("Tissue", 'gid')) 
nrow(tm3)/23
tm3 %>% count(pDE, Reg1) %>% spread(Reg1, n)
tm3 %>% count(hDE, Reg2) %>% spread(Reg2, n)
tm3 = tm3 %>% 
    mutate(Reg1 = as.character(Reg1),
           Reg2 = as.character(Reg2)) %>%
    mutate(Reg1 = ifelse(!is.na(pDE) & pDE != 'non_DE', Reg1, NA),
           Reg2 = ifelse(!is.na(pDE) & pDE == 'non_DE', Reg2, NA)) %>%
    mutate(Reg1 = factor(Reg1, levels = taglst$Reg1),
           Reg2 = factor(Reg2, levels = taglst$Reg2))
tm3 %>% count(pDE, Reg1) %>% spread(pDE, n)
tm3 %>% count(pDE, Reg2) %>% spread(pDE, n)

tm4 = tm3 %>% 
    mutate(MP = (B73 + Mo17) / 2,
           SPE = ifelse(silent == 1, NA,
                 ifelse(B73>=1 & Mo17<0.1 & pDE=='DE_B', 'SPE_B',
                 ifelse(B73<0.1 & Mo17>=1 & pDE=='DE_M', 'SPE_M', 'non_SPE'))),
           HC = ifelse(is.na(SPE), NA,
                ifelse(SPE=='SPE_B' & (BxM>=1 | BxM >= MP), 'HC_B',
                ifelse(SPE=='SPE_M' & (BxM>=1 | BxM >= MP), 'HC_M', 'non_HC')))) %>%
    select(-MP)
tm4 %>% count(Tissue, SPE) %>% spread(SPE, n) %>% print(n=23)
tm = tm4
tm %>% count(Tissue, pDE) %>% spread(pDE, n) %>% print(n=23)
#}}}

#{{{ eQTL
fq = '~/data/misc1/li2013/10.eQTL.v4.tsv'
tq = read_tsv(fq)

types = c("cis", "cis+trans", "trans")
tq2 = tq %>% group_by(gid) %>%
    summarise(ntype = length(unique(type)), type = type[1]) %>%
    mutate(type = ifelse(ntype == 1, type, 'cis+trans')) %>%
    mutate(type = factor(type, levels = rev(types))) %>%
    select(-ntype)
tq2 %>% count(type)

#{{{ baseline piechart - p1
tx1 = tm %>% filter(!is.na(pDE) & pDE != 'non_DE', !is.na(Reg1)) %>%
    mutate(reg = Reg1) %>% select(Tissue, gid, pDE, reg) %>%
    filter(Tissue %in% tissues20) %>%
    mutate(Tissue = factor(Tissue, levels = tissues20)) 
tph = tx1 %>% mutate(ctag = 1, tag = 'All (DE2+)') %>% 
    select(ctag, tag, Tissue, gid, reg) %>% 
    count(reg) %>%
    mutate(prop.baseline = n/sum(n)) %>%
    select(reg, prop.baseline) %>% print(n=10)
p1 = ggplot(tph, aes(x = '', y = prop.baseline, fill = reg)) + 
    geom_bar(width = 1, alpha = 0.8, stat = 'identity') + 
    coord_polar("y", direction = -1) +
    geom_text(aes(label=scales::percent(prop.baseline)), position = position_stack(vjust = .5), size = 3) +
    scale_fill_manual(values = cols.reg1) +
    ggtitle("All (DE2+)") +
    theme_bw() +
    theme(plot.title = element_text(hjust = .5, size = 10)) +
    theme(legend.position = 'none') +
    theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(0,0,1,0), "lines")) +
    theme(panel.border = element_blank(), panel.grid = element_blank()) +
    theme(axis.title = element_blank()) +
    theme(axis.text = element_blank())
#}}}

#{{{ p2: fold change
tissues23
tissues = tissues23[c(17,3,9,18)]
tp = tibble()
for (i in 1:length(tissues)) {
    tissue = tissues[i]
    tp1 = tm %>% filter(Tissue == tissue) %>%
        filter(!is.na(Reg1)) %>%
        select(gid, Reg1) %>% 
        inner_join(tq2, by = 'gid') %>%
        group_by(type, Reg1) %>%
        summarise(n = n()) %>%
        mutate(prop = n/sum(n)) %>% ungroup() %>%
        mutate(tissue = tissue, tag = type, reg = Reg1) %>%
        select(tissue, tag, reg, n, prop) %>%
        inner_join(tph, by = 'reg') %>%
        mutate(fc = log2(prop/prop.baseline))
        #mutate(fc = ifelse(fc < -2, -2, fc))
    tp = rbind(tp, tp1)
}
tp = tp %>% mutate(tissue = factor(tissue, levels = tissues))
tp %>% select(tissue, tag, reg, fc) %>% spread(reg, fc)
tps = tp %>% group_by(tissue, tag, reg) %>% summarise(n = sum(n)) %>%
    mutate(lab = sprintf("%s (%d)", reg, n))
p2 = ggplot(tp) +
    geom_bar(aes(x = tag, y = fc, fill = (as.numeric(tag) %% 2 == 0)), stat = 'identity', alpha = .8, width = .7) +
    geom_hline(yintercept = 0, color = 'gray75') +
    scale_x_discrete(name = 'Li et al eQTL', expand = c(.02,0)) + 
    scale_y_continuous(name = 'Observed proportion / Expected proportion based on DE2+', expand = c(.03,0), breaks = -2:3, labels = 2^(-2:3), limits = c(-2.2, 1)) +
    coord_flip() +
    scale_fill_manual(values = rep(c('gray50','gray50'),2)) +
    facet_grid(tissue ~ reg) +
    theme_bw() +
    theme(strip.placement = "outside") +#, strip.background.y = element_blank()) +
    theme(strip.text.x = element_text(size = 8)) +#, margin=margin(0,0,.1,.5,'lines'))) +
    theme(strip.text.y = element_text(size = 8, angle = 0)) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank()) +
    theme(legend.position = 'none') +
    theme(plot.margin = unit(c(.3,.3,.3,0), "lines")) +
    theme(axis.title = element_text(size = 9)) +
    theme(axis.text = element_text(size = 8))
pg = ggplot_gtable(ggplot_build(p2))
idxs = grep(pattern="strip", pg$layout$name)
for (i in 1:5) {
    idx = idxs[i]
    pg$grobs[[idx]]$grobs[[1]]$children[[1]]
    pg$grobs[[idx]]$grobs[[1]]$children[[1]]$gp$fill = cols.reg1[i]
    pg$grobs[[idx]]$grobs[[1]]$children[[1]]$gp$alpha = .9
    pg$grobs[[idx]]$grobs[[1]]$children[[2]]$children[[1]]$gp$col = 'white'
}
#}}}

fo = file.path(dirw, '85.eQTL.pdf')
ggarrange(p1, pg,
    nrow = 1, ncol = 2, widths = c(1,5))  %>% 
    ggexport(filename = fo, width = 10, height = 3)
#}}}

#{{{ eQTL & GRN
fg = '~/data/genome/B73/v32/t5.gtb'
tg = read_tsv(fg) %>% transmute(gid = par, chrom = chr, start = beg, end = end, srd = srd)

fq = '~/data/misc1/li2013/10.eQTL.v4.tsv'
tq = read_tsv(fq)
tql = tq %>%
    transmute(t.gid = gid, r.chr = qchr, r.pos = qpos, reg = type, addi = ADD) %>%
    filter(r.chr %in% 1:10) %>%
    mutate(r.chr = as.integer(r.chr)) %>%
    arrange(r.chr, r.pos)
tql %>% count(reg)
tqlo = tql %>% transmute(chrom = r.chr, 
                         start = r.pos - 1, end = r.pos, 
                         t.gid = t.gid, reg = reg)

#
grn = 'huang_sam'
tn = t_grn %>% filter(ctag == grn) %>% select(-ctag) %>% 
    top_n(100000, score) %>%
    transmute(r.gid = regulator, t.gid = target, r.fam = fam)
tno = tn %>% distinct(r.gid) %>% inner_join(tg, by = c('r.gid'='gid')) %>%
    arrange(chrom, start) %>%
    transmute(chrom = chrom, 
              start = pmax(0, start - 10000001), 
              end = end + 10000000, gid = r.gid)

fo = file.path(dirw, 't1.tf.bed')
write_tsv(tno, fo, col_names = F)
fo = file.path(dirw, 't2.eqtl.bed')
write_tsv(tqlo, fo, col_names = F)
system("intersectBed -a t1.tf.bed -b t2.eqtl.bed -wo > t3.bed")

fx = file.path(dirw, 't3.bed')
tx = read_tsv(fx, col_names = F)[,c(4,5,7,8,9)]
colnames(tx) = c('r.gid', 'r.chr', 'r.pos', 't.gid', 'reg')
tx %>% distinct(r.gid)
tx %>% distinct(t.gid)
tn2 = tn %>% left_join(tx, by = c('r.gid','t.gid'))
tn2 %>% count(reg)
#}}}

#{{{ marcon2016 DE comparison
dirw = file.path(dird, 'marcon2016')
fi = file.path(dirw, 'marcon2016.tsv')
ti = read_tsv(fi)

fmap = '~/data/genome/B73/gene_mapping/maize.v3TOv4.geneIDhistory.txt'
tmap = read_tsv(fmap, col_names = c("gid", "ngid", "note", "method", "type"))
tmap %>% count(type)

ti2 = ti %>% inner_join(tmap, by = 'gid')
ti2 %>% count(type)

t_mc = ti2 %>% filter(type == '1-to-1') %>%
    transmute(gid = ngid, e.b = B73_wd, e.m = Mo17_wd, e.h = B73xMo17_wd) %>%
    mutate(SPE = ifelse(e.b==1 & e.m==0, 'SPE_B',
                 ifelse(e.b==0 & e.m==1, 'SPE_M', 'non_SPE'))) %>%
    mutate(HC = ifelse(SPE == 'SPE_B' & e.h == 1, 'HC_B',
                ifelse(SPE == 'SPE_M' & e.h == 1, 'HC_M', 'non_HC'))) %>%
    transmute(gid = gid, mSPE = SPE, mHC = HC)
t_mc %>% count(mSPE, mHC) %>% mutate(prop = n/sum(n))

#t_br = tm %>% filter(Tissue == 'root_0DAP', gid %in% t_mc$gid) %>% 
tissue1 = 'seedlingroot_11DAS'
t_br = tm %>% filter(Tissue == tissue1, gid %in% t_mc$gid) %>% 
    select(gid, B73, Mo17, silent, SPE, HC) %>%
    inner_join(t_mc, by = 'gid')
t_br %>% count(silent, SPE)
t_br %>% count(mSPE)
t_br %>% count(SPE, mSPE)
t_br %>% count(SPE, mSPE) %>% 
    filter(!is.na(SPE)) %>% mutate(con = SPE==mSPE) %>%
    group_by(SPE) %>% 
    summarise(ntot = sum(n), pcon = sum(n[which(con)])/ntot)
t_br %>% count(SPE, mSPE) %>% 
    filter(!is.na(SPE)) %>% mutate(con = SPE==mSPE) %>%
    group_by(mSPE) %>% 
    summarise(ntot = sum(n), pcon = sum(n[which(con)])/ntot)

gids_b = t_br %>% filter(mSPE=='SPE_B') %>% pull(gid)
gids_m = t_br %>% filter(mSPE=='SPE_M') %>% pull(gid)
describe(tm %>% filter(Tissue==tissue1, gid %in% gids_b) %>% pull(log2mb))
describe(tm %>% filter(Tissue==tissue1, gid %in% gids_m) %>% pull(log2mb))

fo = file.path(dirw, 'marcon.briggs.tsv')
write_tsv(t_br, fo)
#}}}

#{{{ gene effect
impacts = c('no_change','modifier','low','moderate','high', 'non-syntenic')
fv = '~/projects/wgc/data/05_stats/10.gene.eff.tsv'
tv = read_tsv(fv) %>% filter(qry == 'B73', tgt == 'Mo17') %>%
    select(gid, impact, eff) %>% 
    mutate(impact = factor(impact, levels = impacts))
tv %>% count(impact)

# vnteff on expression broadth
tsh_e %>% inner_join(tv, by = 'gid') %>%
    group_by(impact, etag) %>%
    summarise(ng = n()) %>%
    mutate(pg = ng/sum(ng)) %>%
    select(-ng) %>% spread(etag, pg)

# vnteff on DE
tsh_d %>% filter(ctag == 'pDE') %>% inner_join(tv, by = 'gid') %>%
    group_by(impact, tsTag) %>%
    summarise(ng = n()) %>%
    mutate(pg = ng/sum(ng)) %>%
    select(-ng) %>% spread(tsTag, pg)
tsh_d %>% filter(ctag == 'SPE') %>% inner_join(tv, by = 'gid') %>%
    group_by(impact, tsTag) %>%
    summarise(ng = n()) %>%
    mutate(pg = ng/sum(ng)) %>%
    select(-ng) %>% spread(tsTag, pg)
tsh_d %>% filter(ctag == 'hDE') %>% inner_join(tv, by = 'gid') %>%
    group_by(impact, tsTag) %>%
    summarise(ng = n()) %>%
    mutate(pg = ng/sum(ng)) %>%
    select(-ng) %>% spread(tsTag, pg)
tsh_d %>% filter(ctag == 'pDE') %>% inner_join(tv, by = 'gid') %>%
    group_by(impact, tag) %>%
    summarise(ng = n()) %>%
    mutate(pg = ng/sum(ng)) %>%
    select(-ng) %>% spread(tag, pg)

# vnteff on additivity
tsh_r %>% filter(ctag == 'Dom') %>% inner_join(tv, by = 'gid') %>%
    group_by(impact, tag) %>%
    summarise(ng = n()) %>%
    mutate(pg = ng/sum(ng)) %>%
    select(-ng) %>% spread(tag, pg)
tsh_r %>% filter(ctag == 'Reg1') %>% inner_join(tv, by = 'gid') %>%
    group_by(impact, tag) %>%
    summarise(ng = n()) %>%
    mutate(pg = ng/sum(ng)) %>%
    select(-ng) %>% spread(tag, pg)

#}}}



