source("br.fun.r")
dirw = file.path(dirp, "45.doa")

#{{{ process dom results
fd = file.path(dirw, "09.dom.RData")
load(fd)
td1 = t_dom %>% select(Tissue, gid, bica, bicd, bicr, bico) %>%
    gather(dmode, bic, -Tissue, -gid) %>%
    arrange(Tissue, gid, bic) %>%
    group_by(Tissue, gid) %>%
    summarise(dmode1 = dmode[1], dmode2 = dmode[2]) %>%
    ungroup()

doms = c("BLP", "LP", "PD_L", "MP", "PD_H", "HP", "AHP")
td2 = t_dom %>% select(Tissue, gid, B73, Mo17, BxM, log2MB, pDE, pDE2) %>%
    inner_join(td1, by = c("Tissue", "gid")) %>%
    mutate(LP = pmin(B73, Mo17),
           HP = pmax(B73, Mo17),
           MP = (B73 + Mo17) / 2,
           DoA = (BxM - MP) / (HP - MP)) %>%
    mutate(DoA = ifelse(DoA < -3, -3, DoA)) %>%
    mutate(DoA = ifelse(DoA > 3, 3, DoA)) %>%
    mutate(Dom = ifelse(dmode1 == 'bico',
                 ifelse(BxM > MP & BxM < HP, 'PD_H',
                 ifelse(BxM < MP & BxM > LP, 'PD_L',
                 ifelse(BxM < LP, 'BLP',
                 ifelse(BxM > HP, 'AHP', 'unknown')))), dmode1)) %>%
    mutate(Dom = ifelse(Dom == 'bica', 'MP', Dom)) %>%
    mutate(Dom = ifelse(Dom == 'bicr', 'LP', Dom)) %>%
    mutate(Dom = ifelse(Dom == 'bicd', 'HP', Dom)) %>%
    select(-HP, -LP, -MP) %>%
    mutate(Dom = factor(Dom, levels = doms))
td2 %>% dplyr::count(Tissue, Dom) %>% spread(Dom, n) %>% print(n=23)

td3 = td2 %>%
    mutate(Dom = ifelse(Dom == 'PD_H', ifelse(dmode2 == 'bica', 'MP', 'HP'), Dom)) %>%
    mutate(Dom = ifelse(Dom == 'PD_L', ifelse(dmode2 == 'bica', 'MP', 'LP'), Dom)) %>%
    select(-dmode1, -dmode2)
table(td3$Dom)

doms = c("BLP", "LP", "MP", "HP", "AHP")
t_dom = td3 %>%
    mutate(Tissue = factor(Tissue, levels = tissues23),
           Dom = factor(Dom, levels = doms))
table(t_dom$Dom)

fo = file.path(dirw, "10.dom.RData")
save(t_dom, file = fo)
#}}}

fi = file.path(dirp, "42.de/11.de.RData")
x = load(fi)
td = t_dom %>%
    filter(Tissue %in% tissues20) %>%
    mutate(Tissue = factor(Tissue, levels = tissues20))


x_min = -3; x_max = 3; x_itv = .1
tp = td2 %>%
    filter(Tissue %in% tissues20) %>%
    mutate(Tissue = factor(Tissue, levels = tissues20)) %>%
    filter(pDE != 'non_DE', abs(log2MB) >= 0) %>%
    mutate(x_bin = cut(DoA, breaks=seq(x_min,x_max,by=x_itv), include.lowest=T,
        labels = seq(x_min+x_itv/2,x_max-x_itv/2,by=x_itv))) %>%
    mutate(x_bin = as.numeric(as.character(x_bin))) %>%
    dplyr::count(Tissue, Dom, x_bin)
sum(is.na(tp$x_bin))
p1 = ggplot(tp) +
    geom_area(aes(x = x_bin, y = n, fill = Dom), stat = 'identity', position = 'identity', alpha = .5) +
    geom_vline(xintercept = 0, color = 'gray50') +
    geom_vline(xintercept = c(-1,1), color = 'gray75') +
    scale_x_continuous(name = 'DOA: (F1-MP)/(HP-MP)', limits=c(x_min,x_max)) +
    scale_y_continuous(expand = expand_scale(mult = c(0, .05))) +
    scale_fill_nejm() +
    facet_wrap(~Tissue, ncol = 4, scale = 'free') +
    theme_bw() +
    #theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(.5,.5,.5,.5), "lines")) +
    theme(legend.position = 'top') + 
    guides(fill = guide_legend(nrow = 1, byrow = F)) +
    theme(legend.title = element_blank(), legend.key.size = unit(1, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size = 8)) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank()) +
    theme(axis.title.x = element_text(size = 9)) +
    theme(axis.title.y = element_blank()) +
    theme(axis.text = element_text(size = 8))
fo = file.path(dirw, "15.density.doa.pdf")
ggsave(p1, filename = fo, width = 8, height = 10)
#


x_min = -3; x_max = 3; x_itv = .1
tp = td2 %>%
    filter(Tissue %in% tissues20) %>%
    mutate(Tissue = factor(Tissue, levels = tissues20)) %>%
    filter(pDE != 'non_DE', abs(log2MB) >= 0) %>%
    mutate(x_bin = cut(DoA, breaks=seq(x_min,x_max,by=x_itv), include.lowest = T, 
        labels = seq(x_min+x_itv/2,x_max-x_itv/2,by=x_itv))) %>%
    mutate(x_bin = as.numeric(as.character(x_bin))) %>%
    dplyr::count(Tissue, x_bin)
sum(is.na(tp$x_bin))
p1 = ggplot(tp) +
    geom_area(aes(x = x_bin, y = n), stat = 'identity', position = 'identity', alpha = .5) +
    geom_vline(xintercept = 0, color = 'gray50', size = .2) +
    geom_vline(xintercept = c(-1,1), color = 'gray75', size = .2) +
    scale_x_continuous(name = 'DOA: (F1-MP)/(HP-MP)', limits=c(x_min,x_max)) +
    scale_y_continuous(expand = expand_scale(mult = c(0, .05))) +
    facet_wrap(~Tissue, ncol = 4, scale = 'free') +
    theme_bw() +
    #theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(.5,.5,.5,.5), "lines")) +
    theme(legend.position = 'top') + 
    guides(fill = guide_legend(nrow = 1, byrow = F)) +
    theme(legend.title = element_blank(), legend.key.size = unit(1, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size = 8)) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank()) +
    theme(axis.title.x = element_text(size = 9)) +
    theme(axis.title.y = element_blank()) +
    theme(axis.text = element_text(size = 8))
fo = file.path(dirw, "12.doa.pdf")
ggsave(p1, filename = fo, width = 8, height = 10)
#
