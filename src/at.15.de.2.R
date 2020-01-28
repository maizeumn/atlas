source("functions.R")
sid = 'me99b'
#sid = 'me99b.m'
genome = ifelse(sid == 'me99b.m', 'Mo17', 'B73')
dirw = file.path(dird, ifelse(sid == 'me99b.m', '42_de_m', "42_de"))
gcfg = read_genome_conf(genome)
diri = file.path(dird, ifelse(sid == 'me99b.m', '41_qc_m', "41_qc"))
fm = file.path(diri, '10.rda')
y = load(fm)
fi = file.path(dirw, '10.rda')
x = load(fi)

dd = call_de_dom(t_de, tmm)
fo = file.path(dirw, '11.de.dom.rda')
save(dd, file = fo)

fi = file.path(dirw, '11.de.dom.rda')
x = load(fi)

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
tpx = th %>% distinct(Tissue) %>%
    mutate(Tissue = factor(Tissue, levels = tissues23)) %>%
    arrange(Tissue) %>%
    mutate(x = 1:length(Tissue))

ta = dd %>% group_by(gid) %>%
    summarise(n.exp = n(),
              n.de.b = sum(pDE == 'DE_B'),
              n.de.m = sum(pDE == 'DE_M')) %>% ungroup()
#
gids1 = ta %>% filter(n.de.b == 0, n.de.m == 0, n.exp == 23) %>% pull(gid)
gids2 = ta %>% filter(n.de.b == 0, n.de.m > 10, n.exp == 23) %>% pull(gid)
gids3 = ta %>% filter(n.de.b > 10, n.de.m == 0, n.exp == 23) %>% pull(gid)
gids4 = ta %>% filter(n.de.b > 4, n.de.m > 4, n.exp == 23) %>% pull(gid)
gids = c(gids1, gids2, gids3, gids4)
#
to = tm %>% filter(gid %in% gids) %>%
    inner_join(th[,c('SampleID','Tissue','Genotype')], by = 'SampleID') %>%
    select(gid, Tissue, Genotype, CPM) %>%
    group_by(gid, Tissue, Genotype) %>%
    summarise(cpm.me=mean(CPM), cpm.sd = sd(CPM)) %>% ungroup()

gts = c("B73")
gts = c("B73","Mo17",'BxM')
gts = c("B73","Mo17")
gids = c(gids1[c(21,252,290)], gids3[c(3,284,292)],
         gids2[c(4,23,25)],gids4[c(159,94,110)])
fo = sprintf("%s/31_mixed_de.%d.pdf", dirw, length(gts))
tp = to %>%
    mutate(Tissue = factor(Tissue, levels=tissues23)) %>%
    inner_join(tpx, by = 'Tissue') %>%
    filter(gid %in% gids, Genotype %in% gts) %>%
    mutate(Genotype = factor(Genotype, levels=gts)) %>%
    mutate(gid = factor(gid, levels=gids)) %>%
    replace_na(list(cpm.sd=.1))
#tps = dd %>% filter(gid %in% gids, pDE != 'non_DE') %>%
#    select(Tissue, pDE) %>% mutate(y = ymax, lab = '*')
p1 = ggplot(tp) +
    #geom_line(aes(x=x, y=cpm.me, color=Genotype, group=Genotype), width=.3, linetype='solid') +
    geom_point(aes(x=Tissue, y=cpm.me, color=Genotype, group=Genotype), size=.8) +
    geom_errorbar(aes(x=Tissue, ymin=cpm.me-cpm.sd, ymax=cpm.me+cpm.sd, color = Genotype, group = Genotype), width=.1, alpha=.8) +
    geom_ribbon(aes(x=Tissue, ymin=cpm.me-cpm.sd, ymax=cpm.me+cpm.sd, fill=Genotype, group = Genotype), alpha=.3) +
    #geom_vline(xintercept = seq(1.5,22.5), alpha= .5, linetype='dotted') +
    scale_x_discrete(name = 'Tissue', expand = expand_scale(mult=c(.02,.02))) +
    scale_y_continuous(name = 'CPM (Counts Per Million)', expand = expand_scale(mult=c(.05,.05))) +
    scale_fill_npg() +
    scale_color_npg() +
    facet_wrap(~gid, scale='free_y', ncol=2, dir='v') +
    otheme(xtext=T, xticks=T, ytitle = T, ytext = T, yticks = T, ygrid = T,
           legend.pos = 'top.left', legend.dir='v') +
    theme(axis.text.x=element_text(angle = 30, hjust=1, vjust=1)) +
    guides(color = guide_legend(nrow = 1, byrow = T))
ggsave(p1, filename = fo, width = 12, height = 10)
#}}}


