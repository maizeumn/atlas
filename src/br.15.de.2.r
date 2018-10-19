source("br.fun.r")
sid = 'me99b'
#sid = 'me99b.m'
dirw = file.path(dirp, ifelse(sid == 'me99b.m', '42_de_m', "42_de"))
genome = ifelse(sid == 'me99b.m', 'Mo17', 'B73')
x = load(file.path(dirg, genome, '55.rda'))
diri = file.path(dirp, ifelse(sid == 'me99b.m', '41_qc_m', "41_qc"))
fm = file.path(diri, '10.rda')
y = load(fm)
fi = file.path(dirw, '10.rda')
x = load(fi)

#{{{ call pDE, hDE, D/A, Dom
tm1 = tmm %>% select(Tissue, Genotype, gid, CPM) %>%
    spread(Genotype, CPM) %>%
    mutate(LP = pmin(B73, Mo17),
           HP = pmax(B73, Mo17),
           MP = (B73 + Mo17) / 2,
           DoA = (BxM - MP) / (HP - MP)) %>%
    select(Tissue, gid, B73, Mo17, BxM, DoA) %>%
    mutate(DoA = ifelse(is.nan(DoA)|is.infinite(DoA), NA, DoA))
describe(tm1$DoA)

pDEs = c("DE_B", "DE_M", "non_DE")
tm2 = t_de %>% select(Tissue, deseq) %>% unnest() %>%
    #replace_na(list(log2MB = 0, log2HB = 0, log2HM = 0, log2FM = 0)) %>%
    mutate(tag.mb = ifelse(padj.mb < .01, ifelse(log2mb < 0, -1, 1), 0), 
           tag.hb = ifelse(padj.hb < .01, ifelse(log2hb < 0, -1, 1), 0), 
           tag.hm = ifelse(padj.hm < .01, ifelse(log2hm < 0, -1, 1), 0), 
           tag.fm = ifelse(padj.fm < .01, ifelse(log2fm < 0, -1, 1), 0)) %>%
    mutate(pDE = ifelse(is.na(tag.mb), NA,
                 ifelse(tag.mb == -1, 'DE_B',
                 ifelse(tag.mb == 1, 'DE_M', 'non_DE')))) %>%
    mutate(pDE = factor(pDE, levels = pDEs))
tm2 %>% group_by(Tissue) %>%
    summarise(ng.tot = n(), ng.b = sum(tag.mb==-1), ng.m = sum(tag.mb==1),
              pg.b = ng.b/ng.tot, pg.m = ng.m/ng.tot) %>%
    print(n=23)
#describe(tm2$log2mb)
#describe(tm2$log2hm)

doms = c("BLP", "LP", "PD_L", "MP", "PD_H", "HP", "AHP")
tm3 = tm2 %>%
    filter(tag.mb != 0) %>%
    mutate(tag.lp = ifelse(log2mb > 0, tag.hb, tag.hm),
           tag.hp = ifelse(log2mb > 0, tag.hm, tag.hb)) %>%
    mutate(Dom = "MP") %>%
    mutate(Dom = ifelse(tag.fm == -1 & tag.lp == -1, 'BLP', Dom)) %>%
    mutate(Dom = ifelse(tag.fm == -1 & tag.lp == 0, 'LP', Dom)) %>%
    mutate(Dom = ifelse(tag.fm == -1 & tag.lp == 1, 'PD_L', Dom)) %>%
    mutate(Dom = ifelse(tag.fm == 1 & tag.hp == -1, 'PD_H', Dom)) %>%
    mutate(Dom = ifelse(tag.fm == 1 & tag.hp == 0, 'HP', Dom)) %>%
    mutate(Dom = ifelse(tag.fm == 1 & tag.hp == 1, 'AHP', Dom)) %>%
    mutate(Dom = factor(Dom, levels = doms)) %>%
    select(Tissue, gid, Dom)
tm3 %>% dplyr::count(Tissue, Dom) %>% spread(Dom, n) %>% print(n=23)

doms2 = c("BP", "PL", "AP")
tm4 = tm2 %>%
    filter(tag.mb == 0) %>%
    mutate(hDE = "PL") %>%
    mutate(hDE = ifelse(tag.fm+tag.hb+tag.hm == -3, "BP", hDE)) %>%
    mutate(hDE = ifelse(tag.fm+tag.hb+tag.hm == 3, "AP", hDE)) %>%
    mutate(hDE = factor(hDE, levels = doms2)) %>%
    select(Tissue, gid, hDE)
tm4 %>% dplyr::count(Tissue, hDE) %>% spread(hDE, n) %>% print(n=23)

tx = tm1 %>% left_join(tm2, by = c('Tissue', 'gid')) %>%
    left_join(tm3, by = c('Tissue', 'gid')) %>%
    left_join(tm4, by = c('Tissue', 'gid'))
tx %>% count(pDE, Dom, hDE)

dd = tx
fo = file.path(dirw, '11.de.dom.rda')
save(dd, file = fo)
#}}}

#{{{ log2mb distri.
x_min = -5; x_max = 5; x_itv = .2
tp = dd %>%
    filter(!is.na(pDE), Tissue %in% tissues20) %>%
    mutate(Tissue = factor(Tissue, levels = tissues20)) %>%
    mutate(tag = factor(pDE, levels = levels(dd$pDE)[c(1,3,2)])) %>%
    mutate(x_bin = cut(log2mb, breaks=seq(x_min,x_max,by=x_itv), include.lowest = T, 
        labels = seq(x_min+x_itv/2,x_max-x_itv/2,by=x_itv))) %>%
    mutate(x_bin = as.numeric(as.character(x_bin))) %>%
    count(Tissue, tag, x_bin)
sum(is.na(tp$x_bin))
p1 = ggplot(tp) +
    geom_area(aes(x = x_bin, y = n, fill = tag), stat = 'identity', position = 'identity', alpha = .5) +
    scale_x_continuous(name = 'log2(Mo17/B73)', limits=c(x_min,x_max), expand = c(.05,0)) +
    scale_y_continuous(expand = expand_scale(mult = c(0, .05))) +
    scale_fill_uchicago() +
    facet_wrap(~Tissue, ncol = 4, scale = 'free_y', strip.position = 'top') +
    otheme(legend.pos = 'top.center.out', legend.dir = 'h',
           xtext = T, xtitle = T, ytext = T, xgrid = T) +
    theme(strip.placement = "outside") +
    theme(legend.justification = c(.5, -.3))
fo = file.path(dirw, "15.density.log2mb.pdf")
ggsave(p1, filename = fo, width = 8, height = 9)
#}}}

#{{{ log2hm distri.
x_min = -5; x_max = 5; x_itv = .4
tp = dd %>%
    filter(!is.na(hDE), Tissue %in% tissues20) %>%
    mutate(Tissue = factor(Tissue, levels = tissues20)) %>%
    mutate(tag = hDE) %>%
    mutate(x_bin = cut(log2fm, breaks=seq(x_min,x_max,by=x_itv), include.lowest = T, 
        labels = seq(x_min+x_itv/2,x_max-x_itv/2,by=x_itv))) %>%
    mutate(x_bin = as.numeric(as.character(x_bin))) %>%
    count(Tissue, tag, x_bin)
sum(is.na(tp$x_bin))
p1 = ggplot(tp) +
    geom_area(aes(x = x_bin, y = n, fill = tag), stat = 'identity', position = 'identity', alpha = .5) +
    scale_x_continuous(name = 'log2(Hybrid/Parent)', limits=c(x_min,x_max), expand = c(.05,0)) +
    scale_y_continuous(expand = expand_scale(mult = c(0, .05))) +
    scale_fill_uchicago() +
    facet_wrap(~Tissue, ncol = 4, scale = 'free_y', strip.position = 'top') +
    otheme(legend.pos = 'top.center.out', legend.dir = 'h',
           xtext = T, xtitle = T, ytext = T, xgrid = T) +
    theme(strip.placement = "outside") +
    theme(legend.justification = c(.5, -.3))
fo = file.path(dirw, "15.density.log2hm.pdf")
ggsave(p1, filename = fo, width = 8, height = 9)
#}}}

#{{{ D/A density for different Dom categories
tp0 = dd %>%
    filter(Tissue %in% tissues20) %>%
    mutate(Tissue = factor(Tissue, levels = tissues20)) %>%
    filter(!is.na(Dom))
tk = tp0 %>% count(Tissue, Dom) %>% replace_na(list(n=0)) %>% spread(Dom, n) %>% print(n=23)

x_min = -3; x_max = 3; x_itv = .1
itv_labs = seq(x_min+x_itv/2,x_max-x_itv/2,by=x_itv)
tp = tp0 %>%
    mutate(x_bin = cut(DoA, breaks=seq(x_min,x_max,by=x_itv), include.lowest = T, 
        labels = itv_labs)) %>%
    mutate(x_bin = as.numeric(as.character(x_bin))) %>%
    count(Tissue, Dom, x_bin) #%>% 
    #right_join(t_itv, by = c("Tissue","Dom",'x_bin')) %>% replace_na(list(n=0))
sum(is.na(tp$x_bin))
p1 = ggplot(tp, aes(x = x_bin, y = n, color = Dom, fill = Dom)) +
    geom_area(aes(x = x_bin, y = n), fill = NA, stat = 'identity', position = 'identity', alpha = .5) +
    scale_x_continuous(name = 'Scaled difference: (hybrid-midparent)/(highparent-midparent)', breaks = seq(-2,2,by=1), limits=c(x_min,x_max)) +
    scale_y_continuous(expand = expand_scale(mult = c(0, .05))) +
    scale_color_nejm() +
    facet_wrap(~Tissue, ncol = 4, scale = 'free', strip.position = 'top') +
    otheme(legend.pos = 'top.center.out', legend.dir = 'h',
           yticks = T, xtext = T, xtitle = T, ytext = T, xgrid = T) +
    theme(strip.placement = "outside") +
    guides(color = guide_legend(nrow = 1, byrow = F)) +
    theme(legend.justification = c(.5, -.3))
fo = file.path(dirw, "15.density.doa.pdf")
ggsave(p1, filename = fo, width = 8, height = 9)
#}}}


#{{{ selected examples of mixed-DE patterns
gts2 = c("B73", "Mo17")
tpx = th %>% distinct(Tissue) %>%
    mutate(Tissue = factor(Tissue, levels = tissues23)) %>%
    arrange(Tissue) %>%
    mutate(x = 1:length(Tissue))

ta = dd %>% group_by(gid) %>%
    summarise(n.exp = n(), 
              n.de.b = sum(pDE == 'DE_B'),
              n.de.m = sum(pDE == 'DE_M')) %>% ungroup()
td = ta %>% filter(n.de.b > 4, n.de.m > 4, n.exp == 23) 

to = tm %>% filter(gid %in% td$gid) %>%
    inner_join(th[,c('SampleID','Tissue','Genotype')], by = 'SampleID') %>%
    select(gid, Tissue, Genotype, CPM) %>%
    #mutate(CPM = asinh(CPM)) %>%
    group_by(gid, Tissue, Genotype) %>%
    summarise(cpm.me = mean(CPM), cpm.sd = sd(CPM)) %>% ungroup() #%>%
    #mutate(cond = sprintf("%s %s", Genotype, Tissue)) %>%
    #select(-Genotype, -Tissue) %>%
    #spread(cond, CPM)
fo = file.path(dirw, "30.mixed.de.tsv")
#write_tsv(to, fo)

gids = td$gid[c(2,3,9,10)]
tp = to %>% inner_join(tpx, by = 'Tissue') %>% 
    filter(gid %in% gids, Genotype %in% gts2)
tps = dd %>% filter(gid == !!gid, pDE != 'non_DE') %>%
    select(Tissue, pDE) %>% mutate(y = ymax, lab = '*')
p1 = ggplot(tp) +
    geom_line(aes(x = x, y = cpm.me, color = Genotype, group = Genotype), width = .3, linetype = 'solid') +  
    geom_pointrange(aes(x = x, y = cpm.me, ymin = cpm.me-cpm.sd, ymax = cpm.me+cpm.sd, color = Genotype, group = Genotype)) +
    #geom_vline(xintercept = seq(1.5,22.5), alpha= .5, linetype='dotted') +
    scale_x_continuous(name = 'Tissue', expand = expand_scale(mult=c(.03,.03))) +
    scale_y_continuous(expand = expand_scale(mult=c(.05,.05))) +
    scale_fill_npg() +
    scale_color_npg() +
    facet_wrap(~gid, scale = 'free', ncol = 2) +
    otheme(xtitle = T, ytitle = T, ytext = T, yticks = T, ygrid = T, 
           legend.pos = 'top.center.out') +
    theme(plot.margin = unit(c(1.5,.5,.5,.5), "lines")) +
    guides(color = guide_legend(nrow = 1, byrow = T))
fo = sprintf("%s/31_mixed_de.pdf", dirw)
ggsave(p1, filename = fo, width = 6, height = 5)
#}}}

