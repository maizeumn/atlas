#{{{ head
source("functions.R")
library(ggpubr)
#library(ggridges)
sid = 'me99b'
#sid = 'me99b.m'
dirw = file.path(dird, ifelse(sid == 'me99b.m', '49_coop_m', "49_coop"))
genome = ifelse(sid == 'me99b.m', 'Mo17', 'B73')
genome_cfg = readRDS(file.path(genome_dir(genome), '55.rds'))
diri = file.path(dirp, ifelse(sid == 'me99b.m', '42_de_m', "42_de"))
fi = file.path(dirw, "01.master.rda")
x = load(fi)
tis.prop = 0.5
# tissue distance
#ftd = sprintf("%s/41_qc/50.tissue.dist.RData", dirp)
#load(ftd)
# syntenic proportion
fsyn = '~/data/genome/B73/gene_mapping/syn.gid.tsv'
tsyn = read_tsv(fsyn) %>% distinct(gid) %>% add_column(atag = 'syntenic')
# TF
ftf = '~/data/genome/B73/61_functional/06.tf.tsv'
ttf = read_tsv(ftf) %>% transmute(gid, atag = 'TF')
# GO
fgo = '~/data/genome/B73/61_functional/01.go.tsv'
tgo = read_tsv(fgo)
tdom = tgo %>% filter(ctag %in% c('Interproscan5', 'pannzer')) %>%
    distinct(gid) %>% add_column(atag = 'domain')
thom = tgo %>% filter(ctag %in% c('arabidopsis', 'uniprot.plants')) %>%
    distinct(gid) %>% add_column(atag = 'homology')
tann = rbind(tsyn, tdom, thom)
# color palettes
cols.ts = pal_npg()(3)
cols.expr = pal_npg()(4)
cols.mix3 = pal_d3()(10)[c(1,2,7)]
cols.mix4 = pal_d3()(10)[c(8,1,2,7)]
cols.dom = pal_startrek()(7)
cols.reg1 = pal_d3()(5)
cols.reg2 = pal_d3()(10)[c(8,6,9)]
#}}}

#{{{ num genes w. expression
#{{{ p1 - blank
#require(png)
#require(magick)
#require(cowplot)
#ff = file.path(dirw, "../../Rmd/flowchart.png")
#img = readPNG(ff)
#p0 = rasterGrob(img, interpolate=T)
ff = "http://jeroen.github.io/images/tiger.svg"
#ff = file.path(dird, "../Rmd/flowchart.svg")
#p1 = ggdraw() + draw_image(image_read_svg(ff))
#ff = file.path(dird, "../Rmd/flowchart.pdf")
#img = image_read_pdf(ff, density = 300)
#p1 = ggdraw() + draw_image(image_trim(img))
p1 = ggdraw()
#}}}

#{{{ p2 - #genes expressed in 0-23 tissues
tp = tsh_e %>% group_by(n.tis, etag) %>%
    summarise(num_genes = n())
cat("genes expressed in >=1 tissues:\n")
sum(tp %>% filter(n.tis > 0) %>% pull(num_genes))
cat("prop. genes silent, constitutive. etc:\n")
tp %>% group_by(etag) %>% summarise(n = sum(num_genes)) %>% mutate(p=n/sum(n))
p2 = ggplot(tp) +
    geom_bar(aes(x = n.tis, y = num_genes, fill = etag), stat = 'identity', width = .8) +
    scale_x_continuous(name = 'Number Tissues with Expression', expand=c(0,0)) +
    scale_y_continuous(name = "Number Genes", expand=expand_scale(mult=c(0,.03))) +
    scale_fill_npg() +
    theme_bw() +
    theme(legend.position = c(0,1), legend.justification = c(0,1)) +
    theme(legend.background = element_blank(), legend.direction = "vertical") +
    theme(legend.title = element_blank(), legend.key.size = unit(.8, 'lines'), legend.text = element_text(size = 8)) +
    #theme(axis.ticks.length = unit(0, 'lines')) +
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
    theme(plot.margin = unit(c(.3,1.3,.3,.3), "lines")) +
    theme(axis.title = element_text(size = 9))
#}}}

#{{{ p3: pie chart
tp = tsh_e %>% count(etag)
p3 = ggplot(tp, aes(x = '', y = n, fill = etag)) + 
    geom_bar(width = 1, stat = 'identity') + 
    coord_polar("y", direction = 1) +
    geom_text(aes(label=n), position = position_stack(vjust = 0.5)) +
    scale_fill_npg() +
    theme_bw() +
    theme(legend.position = 'none')+
    theme(panel.border = element_blank(), panel.grid = element_blank()) +
    theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(.3,.3,.3,.3), "lines")) +
    theme(axis.title = element_blank()) +
    theme(axis.text = element_blank())
#}}}

#{{{ p4 a-c: GO/nonsyn
tc = tsh_e %>% select(gid) %>% mutate(tag = 'All genes')
tags = c("All genes", levels(tsh_e$etag))
td = tsh_e %>% 
    transmute(gid = gid, tag = as.character(etag)) %>%
    bind_rows(tc) %>%
    mutate(tag = factor(tag, levels = c(tags)))
tds = td %>% count(tag) %>%
    mutate(lab = sprintf("%s (%5d)", tag, n))
tjs = list(tsyn, tdom, thom)
legend.titles = c("Non-syntenic", "w.o. Known Domain", "w.o. Homologs")
for (i in 1:2) {
    tj = tjs[[i]]
    legend.title = sprintf("Genes %s (%%)", legend.titles[i])
    tp = td %>%
        left_join(tj, by = 'gid') %>% group_by(tag) %>%
        summarise(ngene = n(),
                  n.na = sum(is.na(atag)),
                  p.na = n.na/ngene * 100) %>% ungroup() %>%
        mutate(p.lab = sprintf("%.0f%%", p.na))
    ymax = max(tp$p.na)
    p = ggplot(tp) +
        geom_bar(aes(x = tag, y = p.na, fill = tag), stat = 'identity', width = .8) +
        geom_text(aes(x = tag, y = p.na-1, label = p.lab), color = 'white', size = 2, vjust = 1) +
        scale_x_discrete(breaks = tds$tag, labels = tds$lab, expand = c(.02,0)) +
        scale_y_continuous(name = legend.title, expand = expand_scale(mult=c(0,.03))) +
        scale_fill_manual(values = c(c('gray20', pal_npg()(4)))) +
        #coord_flip() +
        theme_bw() +
        theme(plot.margin = unit(c(.5,.5,.5,.5), "lines")) +
        theme(legend.position = 'none') +
        theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank()) +
        theme(axis.ticks.y = element_blank()) +
        theme(axis.title.x = element_blank()) +
        theme(axis.title.y = element_text(size = 8)) +
        theme(axis.text.y = element_text(size = 8)) +
        theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1))
    if(i == 1) pa = p
    if(i == 2) pb = p
    if(i == 3) pc = p
}
#p4 =  ggarrange(pa, pb, pc, nrow = 1, ncol = 3, widths = c(4.5,2,2), labels = 'D')
p4 =  ggarrange(pa, pb, nrow = 1, ncol = 2, widths = c(2,2))
#}}}

#{{{ enrichment test
tge = tsh_e %>% transmute(tag = etag, gid = gid) %>% distinct(tag, gid)
te = go_enrich_genesets(tge, 1e-4)
fo = file.path(dirw, "51.go.enrich/01.exp.tsv")
write_tsv(te, fo)
#}}}

#{{{ create plot
p_tsne = ggdraw()
fo = file.path(dirw, "11.expr.pdf")
ggarrange(
    p1, p_tsne,
    ggarrange(p2, p4, nrow=1, ncol=2, widths = c(3,2), labels = c('C', 'D')),
    nrow = 3, ncol = 1, labels = c("A",'B'), heights = c(4,5,4)) %>%
    ggexport(filename = fo, width = 8, height = 10)
#}}}
#}}}

#{{{ sum of pDE SPE HC
ctags = c("pDE", "SPE")
#{{{ stats
cat("-----pct SPE that are HC-----\n")
tm %>% count(Tissue, SPE, HC) %>%
    filter(!is.na(SPE), SPE!='non_SPE') %>%
    group_by(Tissue,SPE) %>%
    summarise(prop = n[HC!='non_HC']/sum(n)) %>%
    spread(SPE, prop) %>% print(n=23)
tm %>% #filter(Tissue %in% tissues20) %>%
    count(SPE, HC) %>%
    filter(!is.na(SPE), SPE!='non_SPE') %>%
    group_by(1) %>%
    summarise(prop = sum(n[HC!='non_HC'])/sum(n))
#}}}

#{{{ p1: #genes DE/SPE/HC in each tissue
t_num %>% filter(ctag %in% c("pDE", "SPE", "HC")) %>%
    select(-ctag) %>% spread(tag, ngene) %>%
    mutate(props = (SPE_B+SPE_M)/(DE_B+DE_M), proph = (HC_B+HC_M)/(SPE_B+SPE_M)) %>%
    group_by(1) %>%
    summarise(minDE = min(DE_B+DE_M), maxDE = max(DE_B+DE_M),
              minSPE = min(SPE_B+SPE_M), maxSPE = max(SPE_B+SPE_M),
              minHC = min(HC_B+HC_M), maxHC = max(HC_B+HC_M),
              lprops = min(props), mprops = max(props),
              lproph = min(proph), mproph = max(proph)) %>%
    ungroup()
tz = tsh_d %>% filter(n.tis >= 1, ctag != 'hDE') %>% distinct(ctag, gid) %>% count(ctag)
print(tz)
print(tz$n[tz$ctag=='SPE']/tz$n[tz$ctag=='pDE'])
print(tz$n[tz$ctag=='HC']/tz$n[tz$ctag=='SPE'])

t_num_de = tm %>% filter(!is.na(pDE), pDE != 'non_DE') %>%
    mutate(afc = abs(log2mb)) %>%
    mutate(fc = ifelse(afc < 1, 'd12', ifelse(afc < 2, 'd24',
                ifelse(afc < 3, 'd48', 'd8')))) %>%
    count(Tissue, pDE, fc)
t_num_de %>% group_by(Tissue) %>%
    summarise(prop12 = sum(n[fc=='d12'])/sum(n)) %>% print(n=23)
t_num_spe = t_num %>% filter(ctag == 'SPE') %>%
    mutate(pDE = ifelse(tag == 'SPE_B', 'DE_B', 'DE_M'),
           fc = 'SPE') %>%
    #mutate(pDE = factor(pDE)) %>%
    select(Tissue, pDE, fc, n = ngene)
tp = t_num_de %>% mutate(pDE = as.character(pDE)) %>%
    bind_rows(t_num_spe) %>%
    spread(fc, n) %>%
    mutate(d4 = d48 +d8, d2 = d24+d4, d1 = d12+d2) %>%
    select(Tissue, pDE, d1,d2,d4,d8,SPE) %>%
    gather(fc, n, -Tissue, -pDE) %>%
    mutate(Tissue=factor(Tissue,levels=tissues23))
#tps = tp0 %>% distinct(pDE, fc, tag) %>% arrange(pDE,fc)
#tp = tp %>% mutate(tag = factor(tag, levels = tps$tag))
#
tags = c('DE (All)', 'DE (2+)', 'DE (4+)', 'DE (8+)', 'SPE')
p1 = ggplot(tp) +
    geom_point(aes(x = Tissue, y = n, shape = fc, group = pDE, color = pDE), size = 2, position = position_dodge(width=1)) +
    geom_vline(xintercept = seq(1.5, 22.5, by = 1), linetype=2, color='grey60') +
    scale_x_discrete(name = '') +
    scale_y_continuous(name = 'Number Genes', limits = c(0,max(tp$n)), expand = expand_scale(mult=c(0,.05))) +
    scale_shape_manual(values = c(15,24,23,25,19), labels = tags) +
    scale_color_d3(labels = c("B73 higher", "Mo17 higher")) +
    otheme(legend.pos = 'top.center.out', legend.dir = 'h',
           xtext = T, ytitle = T, ytext = T, xticks = T, yticks = T,
           ygrid = T) +
    guides(direction = 'horizontal', color = guide_legend(nrow = 1), shape = guide_legend(nrow = 1)) +
    theme(legend.box = 'horizontal') +
    theme(plot.margin = unit(c(1.5,.5,.5,1.5), "lines")) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
#ggsave(p1, filename = file.path(dirw, 'test.pdf'), width = 9, height = 4)
#}}}

#{{{ p2: #genes DE/SPE/HC shared in 1-23 tissues
tsh_d %>% filter(ctag == 'SPE', tag != 'non-SPE', n.tis >= 10)
tsh_d %>% filter(ctag == 'SPE', tag != 'non-SPE', n.tis == 23)
for (ctag1 in ctags) {
    tp = tsh_d %>% filter(ctag == ctag1, n.tis >= 1) %>% 
        count(n.tis, tag) %>%
        mutate(fac = ifelse(n.tis <= 2, 1, 2))
    print(sum(tp$n))
    print(sum(tp$n[tp$n.tis >= 10]))
    xtag = ifelse(ctag1 == 'pDE', 'DE', ctag1)
    xlabel = sprintf("Number %s Tissues", xtag)
    tp2 = tp %>% filter(n.tis >= 3) %>%
        group_by(tag) %>% summarise(n = sum(n)) %>% mutate(p = n/sum(n))
    print(tp2)
    print(sum(tp2$n))
    p = ggplot(tp) +
        geom_bar(aes(x = n.tis, y = n, fill = tag), position = 'stack', stat='identity', width=0.7) +
        scale_x_continuous(name = xlabel, breaks = c(1,2,5,10,15,20), expand = c(0.01,0)) +
        scale_y_continuous(name = "Number Genes", expand = c(0.02, 0)) +
        scale_fill_manual(values = cols.mix3) +
        facet_wrap(~fac, scales='free') +
        theme_bw() +
        theme(legend.position = c(1,1), legend.justification=c(1,1)) +
        guides(direction = 'vertical', fill = guide_legend(ncol = 1, byrow = F)) +
        theme(legend.title = element_blank(), legend.key.size = unit(0.8, 'lines'), legend.text = element_text(size = 8)) +
        theme(panel.border = element_blank(), panel.grid = element_blank()) +
        theme(strip.background = element_blank(), strip.text.x = element_blank()) +
        theme(plot.margin = unit(c(.5,.5,.5,.5), "lines")) +
        theme(axis.title.x = element_text(size = 9)) +
        theme(axis.title.y = element_text(size = 9)) +
        theme(axis.text.x = element_text(size = 8, color = "black", angle = 0)) +
        theme(axis.text.y = element_text(size = 8, color = "black", angle = 0))
    gt = ggplotGrob(p)
    N <- tp %>% group_by(fac) %>% summarise(count = length(unique(n.tis))) %>% `[[`(2)
    panelI <- gt$layout$l[grepl("panel", gt$layout$name)]
    gt$widths[panelI] <- unit(N, "null")
    gt$widths[panelI[1] + 1] = unit(0.05, "cm")
    p = gt
    if(ctag1 == 'pDE') p2a = p
    if(ctag1 == 'SPE') p2b = p
    if(ctag1 == 'HC') p2c = p
}
p2 = ggarrange(p2a, p2b, nrow = 1, ncol = 2, labels = c("B","C"))
#}}}

#{{{ check 4-FC DE
gids.de2 = tm %>% filter(!is.na(pDE) & pDE != 'non_DE') %>% distinct(gid) %>% pull(gid)
tm4 = tm %>% mutate(pDE = ifelse(abs(log2mb) < 2, 'non_DE', pDE))
    tr1 = tm4 %>% mutate(tag = pDE) %>%
        filter(!is.na(tag)) %>%
        group_by(gid, tag) %>%
        summarise(n.tis = n()) %>%
        ungroup()
    trs = tr1 %>%
        group_by(gid) %>%
        summarise(n.tis.tot = sum(n.tis)) %>% ungroup()
    tr1 = tr1 %>% filter(! tag %in% c("non_DE")) %>%
        arrange(gid, n.tis) %>%
        group_by(gid) %>%
        summarise(n.tis = sum(n.tis),
                  tag = ifelse(length(unique(tag)) == 1, 
                               sprintf("consis. %s", tag[1]), 
                               sprintf("mix of %s", "DE_B/DE_M"))) %>%
        ungroup() %>%
        right_join(trs, by = 'gid') %>%
        replace_na(list(n.tis = 0, tag = 'nothing')) %>%
        mutate(prop.tis = n.tis / n.tis.tot,
               tsTag = ifelse(prop.tis == 0, 'silent',
                       ifelse(prop.tis <= 0.2, 'tis-specific',
                       ifelse(prop.tis < 0.8, 'intermediate', 'constit.')))) %>%
        select(gid, tag, tsTag, n.tis, n.tis.tot, prop.tis)
    tr = tr1 %>% filter(n.tis >= 1) %>% 
        count(n.tis, tag) %>%
        filter(n.tis >= 2) %>%
        group_by(tag) %>% summarise(n = sum(n)) %>% mutate(p = n/sum(n)) %>%
        print(n=20)
    sum(tr$n)
    length(gids.de2)
    sum(tr$n)/length(gids.de2)
#}}}

#{{{ p3a: per tissue stats
ntissue = max(tsh_d$n.tis.tot)
tis.prop = .5
tp = tsh_d %>% filter(n.tis.tot >= tis.prop * ntissue, n.tis >= 1) %>%
    filter(ctag %in% c("pDE", "SPE", "HC"))  %>%
    mutate(ctag = as.character(ctag)) %>%
    mutate(ctag = ifelse(ctag=='pDE', 'DE', ctag)) %>%
    mutate(ctag = factor(ctag, levels = rev(c("DE", "SPE", "HC")))) %>%
    count(ctag, tsTag)
tps = tp %>% group_by(ctag) %>% summarise(n = sum(n)) %>%
    mutate(lab = sprintf("%s (%d)", ctag, n))
p3a = ggplot(tp) +
    geom_bar(aes(x = ctag, y = n, fill = tsTag), position = position_fill(reverse = T), stat='identity', color='white', size=0.1, width = 0.7) +
    scale_x_discrete(breaks = tps$ctag, labels = tps$lab, expand = c(0.01,0)) +
    scale_y_continuous(name = sprintf("Proportion Genes"), breaks = seq(0,1,by=0.25), expand = c(0.01, 0)) +
    coord_flip() +
    scale_fill_npg() +
    theme_bw() +
    theme(legend.position = 'top') +
    theme(legend.position = c(1,1),  legend.justification=c(1,0)) +
    guides(direction = 'horizontal', fill = guide_legend(nrow = 1, byrow = F)) +
    theme(legend.title = element_blank(), legend.key.size = unit(0.8, 'lines'), legend.text = element_text(size = 8)) +
    theme(panel.border = element_blank(), panel.grid = element_blank()) +
    theme(plot.margin = unit(c(3.5,.5,.5,1.5), "lines")) +
    theme(axis.title.x = element_text(size = 9)) +
    theme(axis.title.y = element_blank()) +
    theme(axis.text = element_text(size = 8))
#}}}

#{{{ p3: tissue specificity
tsmap = c("No data" = "not DE/SPE in any tissue",
          "Tissue specific" = "DE/SPE in <20% expressed tissues",
          "Intermediate frequency" = "DE/SPE in 20-80% expressed tissues",
          "Constitutive" = "DE/SPE in >80% expressed tissues")
tp = tsh_d %>% filter(ctag %in% !!ctags, n.tis >= 1, n.tis.tot >= 10) %>% 
    group_by(ctag, tsTag) %>%
    summarise(n = n()) %>% mutate(prop = n/sum(n)) %>% ungroup() %>%
    mutate(ctag = as.character(ctag)) %>%
    mutate(ctag = ifelse(ctag == 'pDE', 'DE', ctag)) %>%
    mutate(tsTag = tsmap[tsTag]) %>%
    mutate(tsTag = factor(tsTag, levels = as.character(tsmap)))
p3 = ggplot(tp, aes(x = '', y = prop, fill = tsTag)) + 
    geom_bar(width = 1, stat = 'identity', alpha = .8) + 
    coord_polar("y", direction = 1) +
    geom_text(aes(label=n), position = position_stack(vjust = .5)) +
    facet_wrap(~ctag, nrow = 1) +
    scale_fill_startrek() +
    theme_bw() +
    theme(legend.position = c(.5,1),  legend.justification=c(.5,-.5)) +
    theme(legend.background = element_blank()) +
    theme(legend.title = element_blank(), legend.key.size = unit(.8, 'lines'), legend.text = element_text(size = 8)) +
    theme(panel.border = element_blank(), panel.grid = element_blank()) +
    theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(1.5,.5,.5,.5), "lines")) +
    theme(axis.title = element_blank()) +
    theme(axis.text = element_blank())
#}}}

#{{{ p4: expression v.s. DE/SPE specificity
clist = list(
    'pDE' = c("consis. non-DE", "consis. DE_B", "consis. DE_M", "mix"),
    'SPE' = c("consis. non-SPE", "consis. SPE_B", "consis. SPE_M", "mix"),
    'HC' = c("consis. non-HC", "consis. HC_B", "consis. HC_M", "mix")
)
itv.brks = seq(0, 1, by = .05)
itv.labs = sprintf("%.02f-%.02f", itv.brks[-length(itv.brks)], itv.brks[-1])
itv.labs.all = c('0', itv.labs)
for (ctag1 in names(clist)) {
    mtags = clist[[ctag1]]
    otag = ifelse(ctag1 == 'pDE', 'DE', ctag1)
    tp0 = tsh_d %>% filter(ctag == ctag1) %>%
        mutate(tag = ifelse(str_detect(tag, 'non-'), mtags[1], tag)) %>%
        mutate(tag = ifelse(str_detect(tag, 'mix'), 'mix', tag)) %>%
        mutate(tag = factor(tag, levels = mtags)) %>%
        mutate(ptis.tag = cut(prop.tis, breaks = itv.brks, labels = itv.labs,
                              include.lowest = F)) %>%
        mutate(ptis.tag = as.character(ptis.tag)) %>%
        replace_na(list(ptis.tag = '0')) %>%
        mutate(ptis.tag = factor(ptis.tag, levels = itv.labs.all)) %>%
        transmute(gid=gid, ptis.tag = ptis.tag, ntis.exp = n.tis.tot, mtag = tag)
    tp1 = tp0 %>% count(ntis.exp, ptis.tag)
    tpm = tp0 %>% count(ntis.exp, mtag)
    tp1 %>% filter(ntis.exp == 3) %>% mutate(prop = n/sum(n)) %>% print(n = 20)
    tps1 = tp1 %>% group_by(ntis.exp) %>%
        summarise(ng.tot = sum(n)) %>% ungroup() %>%
        mutate(lab = sprintf("%2d (%4d)", ntis.exp, ng.tot))
    tpsm = tpm %>% group_by(ntis.exp) %>%
        summarise(ng.tot = sum(n)) %>% ungroup() %>%
        mutate(lab = sprintf("%2d (%4d)", ntis.exp, ng.tot))
    pm = ggplot(tpm) +
        geom_bar(aes(x = ntis.exp, y = n, fill = mtag), position = position_fill(reverse = T), stat = 'identity', width = 0.8) +
        scale_x_continuous(name = 'Number Expressed Tissues', breaks = tpsm$ntis.exp, labels = tpsm$lab, expand = c(0,0)) +
        scale_y_continuous(name = 'Proportion Genes', expand = c(0,0)) +
        scale_fill_npg() +
        coord_flip() +
        theme_bw() +
        theme(legend.position = c(1,1), legend.justification=c(1,0), legend.background = element_blank()) +
        guides(direction = 'horizontal', fill = guide_legend(nrow = 1, byrow = F)) +
        theme(legend.title = element_blank(), legend.key.size = unit(0.7, 'lines'), legend.text = element_text(size = 8)) +
        #theme(legend.margin = margin(0,0,0,0, unit='pt')) +
        theme(panel.grid = element_blank()) +
        theme(plot.margin = unit(c(1.5,.5,.5,.5), "lines")) +
        theme(axis.title.x = element_text(size=9)) +
        theme(axis.title.y = element_text(size=9)) +
        theme(axis.text.x = element_text(size=8)) +
        theme(axis.text.y = element_text(size=8))
    ncolor = length(unique(tp1$ptis.tag))
    cols21 = c('gray80', viridis_pal(direction = -1)(20))
    ntiss = as.character(sort(unique(tp1$ptis.tag)))
    tp1 = tp1 %>% mutate(ptis.tag = factor(as.character(ptis.tag), levels = rev(ntiss)))
    title.1 = sprintf("Proportion Tissues %s", otag)
    pa = ggplot(tp1) +
        geom_bar(aes(x = ntis.exp, y = n, fill = ptis.tag), position = position_fill(), stat = 'identity', width = .8) + 
        scale_x_continuous(name = 'Number Expressed Tissues', breaks = tps1$ntis.exp, labels = tps1$lab, expand = c(0,0)) +
        scale_y_continuous(name = 'Proportion Genes', expand = c(0,0)) +
        scale_fill_manual(name = title.1, values = rev(cols21), breaks = itv.labs.all[seq(1,21,by=5)]) +
        coord_flip() +
        theme_bw() +
        theme(legend.position = c(1,1), legend.justification=c(1,0)) +
        theme(legend.direction = 'vertical', legend.background = element_blank()) +
        guides(direction = 'horizontal', fill = guide_legend(title.hjust = .5, nrow = 1, byrow = F)) +
        theme(legend.title = element_text(size = 9), legend.key.size = unit(0.7, 'lines'), legend.text = element_text(size = 8)) +
        #theme(legend.key = element_rect(color = 'black', size = .5)) +
        #theme(legend.margin = margin(0,0,0,0, unit='pt')) +
        theme(panel.grid = element_blank()) +
        theme(plot.margin = unit(c(2.3,.5,.5,.5), "lines")) +
        theme(axis.title.x = element_text(size=9)) +
        theme(axis.title.y = element_text(size=9)) +
        theme(axis.text.x = element_text(size=8)) +
        theme(axis.text.y = element_text(size=8))
    if(ctag1 == 'pDE') {ppa = pa; ppm = pm}
    if(ctag1 == 'SPE') {psa = pa; psm = pm}
    if(ctag1 == 'HC') {pha = pa; phm = pm}
}
#}}}

#{{{ p5: non-syn
tc1 = tsh_e %>% select(gid) %>% 
    mutate(ctag = 'control', tag = 'All genes')
tc2 = tsh_e %>% filter(n.tis >= 10) %>% select(gid) %>% 
    mutate(ctag = 'control', tag = 'All expr. genes')
tc3 = tsh_d %>% filter(tag == 'non-pDE', n.tis.tot >= 10) %>% select(gid) %>%
    mutate(ctag = 'control', tag = "non-DE genes")
tc = rbind(tc1, tc2, tc3) %>% select(ctag, tag, gid)
tr = tsh_d %>% filter(n.tis.tot >= 10, tsTag != 'No data', ctag %in% ctags) %>%
    arrange(ctag, desc(tsTag)) %>%
    transmute(ctag = as.character(ctag), tag = as.character(tsTag), gid = gid) %>%
    mutate(ctag = ifelse(ctag == 'pDE', 'DE', ctag)) %>%
    mutate(tag = ifelse(tag == 'Intermediate frequency', 'Intermediate', tag)) %>%
    mutate(tag = ifelse(tag == 'Tissue specific', 'Tissue-specific', tag))
ctags1 = unique(tr$ctag)
tra = tr %>% mutate(tag = 'All')
trb = tr %>% filter(tag == 'Constitutive')
tr = rbind(tra, trb) %>%
    mutate(ctag = factor(ctag, levels = ctags1)) %>%
    arrange(ctag, tag) %>%
    mutate(tag = sprintf("%s|%12s", ctag, tag))
tags = c(unique(tc$tag), unique(tr$tag))
td = rbind(tc, tr)
ctags = unique(td$ctag)
td = td %>%
    mutate(tag = factor(tag, levels = rev(tags))) %>%
    mutate(ctag = factor(ctag, levels = ctags))
tds = td %>% count(tag) %>%
    mutate(lab = sprintf("%s (%5d)", tag, n))
tjs = list(tsyn, tdom, thom)
legend.titles = c("Non-syntenic", "w.o. Known Domain", "w.o. Homologs")
for (i in 1:3) {
    tj = tjs[[i]]
    legend.title = sprintf("Genes %s (%%)", legend.titles[i])
    tp = td %>%
        left_join(tj, by = 'gid') %>% group_by(ctag, tag) %>%
        summarise(ngene = n(),
                  n.na = sum(is.na(atag)),
                  p.na = n.na/ngene * 100) %>% ungroup() %>%
        mutate(p.lab = sprintf("%.0f%%", p.na))
    ymax = max(tp$p.na)
    p = ggplot(tp) +
        geom_bar(aes(x = tag, y = p.na, fill = ctag), stat = 'identity', alpha = .8, width = .8) +
        geom_text(aes(x = tag, y = p.na-1, label = p.lab), color = 'white', size = 2.5, hjust = 1) +
        scale_x_discrete(breaks = tds$tag, labels = tds$lab, expand = c(0.02,0)) +
        scale_y_continuous(name = legend.title, expand = expand_scale(mult=c(0,.03))) +
        scale_fill_manual(values = c("gray20", pal_nejm()(3))) +
        coord_flip() +
        theme_bw() +
        theme(plot.margin = unit(c(.5,.5,.5,.5), "lines")) +
        theme(legend.position = 'none') +
        theme(panel.grid = element_blank()) +
        theme(axis.ticks.y = element_blank()) +
        theme(axis.title.x = element_text(size = 9, hjust = .5)) +
        theme(axis.title.y = element_blank()) +
        theme(axis.text.x = element_text(size=8)) +
        theme(axis.text.y = element_text(size=8, family='mono', face='bold'))
    if(i != 1) p = p + theme(axis.text.y = element_blank())
    if(i == 1) p5a = p
    if(i == 2) p5b = p + theme(plot.margin = unit(c(.5,.1,.5,.1), 'lines'))
    if(i == 3) p5c = p
}
#}}}

#{{{ enrichment
tge = tsh_d %>% filter(ctag == 'pDE', n.tis.tot >= 10, n.tis >= 1) %>% 
    transmute(tag = tsTag, gid = gid) %>%
    distinct(gid, tag)
te = go_enrich_genesets(tge, 1e-4)
fo = file.path(dirw, "51.go.enrich/10.pDE.tsv")
write_tsv(te, fo)

tge = tsh_d %>% filter(ctag == 'pDE', n.tis.tot >= 10, n.tis >= 1) %>% 
    transmute(tag = tag, gid = gid) %>%
    distinct(gid, tag)
te = go_enrich_genesets(tge, 1e-4)
fo = file.path(dirw, "51.go.enrich/11.pDE.mix.tsv")
write_tsv(te, fo)

tge = tsh_d %>% filter(ctag == 'pDE', n.tis.tot == 23) %>%
    mutate(tag = ifelse(n.tis == 0, 'a: DE0', ifelse(n.tis <=3, 'b: DE1-3',
                 ifelse(n.tis <= 22, 'c: DE4-22', 'd: DE23')))) %>%
    select(tag, gid) %>% arrange(tag, gid)
te = go_enrich_genesets(tge, 1e-4)
fo = file.path(dirw, "51.go.enrich/12.pDE.exp23.tsv")
write_tsv(te, fo)
#}}}

#{{{ plot
fo = file.path(dirw, "25.DE.pdf")
p2 = ggarrange(p2a, p2b, p3, nrow = 1, ncol = 3, labels = c("B","C","D"))
p4 = ggarrange(ppa, psa, nrow = 1, ncol = 2, labels = c("E","F"))
p5 =  ggarrange(p5a, p5b, p5c, nrow = 1, ncol = 3, widths = c(3.5,2,2), labels = 'G')
ggarrange(p1, p2, p4, p5,
    nrow = 4, ncol = 1, labels = "A", heights = c(4,3,4,2)) %>%
    ggexport(filename = fo, width = 9, height = 11)
#
fo = file.path(dirw, "25.DE.sup.pdf")
ggarrange(ppa, psa, ppm, psm, nrow = 2, ncol = 2, labels = LETTERS[1:4]) %>%
    ggexport(filename = fo, width = 8, height = 7)

if(genome == 'Mo17') {
fo = file.path(dirw, "25.DE.pdf")
p2 = ggarrange(p2a, p2b, p3, nrow = 1, ncol = 3, labels = c("B","C","D"))
ggarrange(p1, p2,
    nrow = 2, ncol = 1, labels = "A", heights = c(4,3)) %>%
    ggexport(filename = fo, width = 9, height = 6)
}
#}}}
#}}}

#{{{ sum of hDE
#{{{ p1: per tissue stats
tm %>% count(Tissue, pDE) %>% filter(!is.na(pDE) & pDE == 'non_DE') %>%
    group_by(1) %>%
    summarise(nde.min = min(n), nde.max = max(n))
tp = t_num %>% filter(ctag %in% c("hDE")) %>% select(-ctag)
tps = tp %>% group_by(Tissue) %>% summarise(ngene = sum(ngene)) %>%
    mutate(lab = sprintf("%s (%d)", Tissue, ngene))
t_num %>% filter(ctag %in% c("hDE")) %>%
    select(-ctag) %>% spread(tag, ngene) %>%
    replace_na(list(AP=0, BP=0, PL=0)) %>%
    group_by(1) %>%
    summarise(minDE = min(BP+AP), maxDE = max(BP+AP)) %>%
    ungroup()
tsh_d %>% filter(ctag == 'hDE', n.tis >= 1) %>% distinct(gid) %>% count(1)
p1 = ggplot(tp) +
    geom_bar(aes(x = Tissue, y = ngene, fill = tag), position = position_stack(reverse = T), stat = 'identity', width = 0.8) +
    scale_x_discrete(breaks = tps$Tissue, labels = tps$lab, expand = c(0,0)) +
    scale_y_continuous(name = sprintf("Number Genes"), expand = expand_scale(mult=c(0,.03))) +
    scale_fill_d3()+
    coord_flip() +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    theme(plot.margin = unit(c(1.5,.5,.5,.5), "lines")) +
    theme(legend.position = c(1,1), legend.justification=c(1,0), legend.background = element_blank()) +
    guides(direction = 'vertical', fill = guide_legend(nrow = 1, byrow = F)) +
    theme(legend.title = element_blank(), legend.key.size = unit(0.8, 'lines'), legend.text = element_text(size = 8)) +
    theme(axis.title.x = element_text(size=9)) +
    theme(axis.title.y = element_blank())

p1b = ggplot(tp) +
    geom_point(aes(x = Tissue, y = ngene, group = tag, color = tag, shape = tag), size = 2) +
    scale_x_discrete(breaks = tps$Tissue, labels = tps$lab) +
    scale_y_continuous(name = 'Number Genes', expand = expand_scale(mult=c(0,.05))) +
    scale_shape_manual(values = c(16, 4)) +
    scale_color_d3() +
    theme_bw() +
    theme(plot.margin = unit(c(.5,.5,.5,1.5), "lines")) +
    theme(legend.position = c(1,1), legend.direction = 'vertical', legend.justification = c(1,1)) +
    theme(legend.background = element_blank()) +
    theme(legend.title = element_blank(), legend.key.size = unit(0.8, 'lines'), legend.text = element_text(size = 8)) +
    theme(panel.grid.minor = element_blank()) +
    theme(panel.grid.major.x = element_line(linetype='dashed')) +
    theme(axis.title.x = element_blank()) +
    theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1))
#}}}

#{{{ p2
tsh_d %>% filter(ctag == 'hDE', n.tis >= 1) %>%
    distinct(gid, n.tis) %>% count(n.tis) %>% mutate(p.tis = n/sum(n))
tp = tsh_d %>% filter(ctag == 'hDE', n.tis >= 1) %>% count(n.tis, tag)
sum(tp %>% filter(n.tis >= 3) %>% pull(n))
sum(tp %>% filter(n.tis >= 3) %>% pull(n)) / sum(tp$n)
tp2 = tp %>% filter(n.tis >= 2) %>%
    group_by(tag) %>% summarise(n = sum(n)) %>% mutate(p = n/sum(n))
sum(tp2$n)
print(tp2)
ymax = max(tp %>% group_by(n.tis) %>% summarise(n=sum(n)) %>% pull(n)) * 1.05
p2 = ggplot(tp) +
    geom_bar(aes(x = n.tis, y = n, fill = tag), position = 'stack', stat='identity', width=0.7) +
    scale_x_continuous(name = 'Number Tissues', breaks = seq(1,10,by=1), expand = c(0.01,0)) +
    scale_y_continuous(name = "Number Genes", limits = c(0,ymax), expand = c(0,0)) +
    scale_fill_manual(values = cols.mix3) +
    theme_bw() +
    theme(legend.position = c(1,1), legend.justification=c(1,1), legend.background = element_blank()) +
    theme(legend.title = element_blank(), legend.key.size = unit(0.8, 'lines'), legend.text = element_text(size = 8)) +
    theme(panel.grid = element_blank()) +
    #theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(.3,.3,.3,.3), "lines")) +
    theme(axis.title = element_text(size = 9)) +
    theme(axis.text = element_text(size = 8))
#}}}

#{{{ p3: pie-chart of tissue specificity
tp = tsh_d %>% filter(ctag == 'hDE', n.tis >= 1, n.tis.tot >= 10) %>% count(tsTag)
p3 = ggplot(tp, aes(x = '', y = n, fill = tsTag)) + 
    geom_bar(width = 1, stat = 'identity', alpha = .8) + 
    coord_polar("y", direction = 1) +
    geom_label_repel(aes(label=n), nudge_y = -1) +# position = position_stack(vjust = 0.5)) +
    scale_fill_npg(labels = c("DE in <20% tissues", "DE in 20-80% tissues")) +
    theme_bw() +
    theme(legend.position = c(0.5,1),  legend.justification=c(0.5,0.3)) +
    theme(legend.title = element_blank(), legend.key.size = unit(0.8, 'lines'), legend.text = element_text(size = 8)) +
    theme(panel.border = element_blank(), panel.grid = element_blank()) +
    theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(1.5,.5,.5,.5), "lines")) +
    theme(axis.title = element_blank()) +
    theme(axis.text = element_blank())
#}}}

#{{{ p4: pie-chart of fold change
tp = tsh_d %>% filter(ctag == 'hDE', n.tis >= 1, n.tis.tot >= 10)
tp = tm %>% filter(gid %in% tp$gid) %>%
    select(Tissue, gid, hDE, log2fm) %>%
    filter(!is.na(hDE) & hDE != 'PL') %>%
    group_by(gid) %>% 
    summarise(log2FMa=min(abs(log2fm)), log2FMb=max(abs(log2fm))) %>%
    mutate(itv = cut(log2FMa, breaks = c(0,1,log2(3),2,log2(5),Inf), 
                     include.lowest = T, 
                     labels = c("1-2", "2-3", "3-4", "4-5", "5+"))) %>%
    count(itv)
tp$n[tp$itv == '5+'] 
tp$n[tp$itv == '5+']/sum(tp$n) 
sum(is.na(tp$itv))
tp %>% mutate(prop = n/sum(n))
p4 = ggplot(tp, aes(x = '', y = n, fill = itv)) + 
    geom_bar(width = 1, stat = 'identity', alpha = .8) + 
    coord_polar("y", direction = 1) +
    geom_label_repel(aes(label=n), nudge_y = -2) +
    scale_fill_npg(name = 'Fold change') +
    theme_bw() +
    theme(legend.position = c(.5,1), legend.justification=c(.5,0), legend.direction = 'horizontal') +
    theme(legend.title = element_text(size=9), legend.key.size = unit(0.8, 'lines'), legend.text = element_text(size = 8)) +
    theme(panel.border = element_blank(), panel.grid = element_blank()) +
    theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(1.5,.5,.5,.5), "lines")) +
    theme(axis.title = element_blank()) +
    theme(axis.text = element_blank())
#}}}

#{{{ enrichment
tge = tm %>% filter(!is.na(hDE) & hDE != 'PL') %>%
    select(Tissue, gid, hDE) %>%
    transmute(tag = sprintf("%s %s", Tissue, hDE), gid = gid)
tge2 = tm %>% filter(!is.na(hDE) & hDE != 'PL') %>%
    distinct(gid, hDE) %>%
    transmute(tag = hDE, gid = gid)
tge = rbind(tge, tge2)
te = go_enrich_genesets(tge, pval.cutoff = .005)
fo = file.path(dirw, "51.go.enrich/16.hDE.tsv")
write_tsv(te, fo)
#}}}

#{{{ plot
fo = file.path(dirw, "25.hDE.pdf")
ggarrange(p1, p2, p3, p4,
    nrow = 2, ncol = 2, labels = LETTERS[1:4], heights = c(3,2)) %>%
    ggexport(filename = fo, width = 8, height = 6)
#}}}
#}}}

#{{{ sum of Dom
#{{{
ntissue = max(tsh_r$n.tis.tot)
de.prop = .25
doms = taglst$Dom
listd = list("LP/BLP" = c("LP", "BLP"), "HP/AHP" = c("HP", "AHP"))
tnd = t_num %>% filter(ctag %in% c("Dom")) %>%
    select(-ctag) %>%
    filter(Tissue %in% tissues20) %>%
    mutate(Tissue = factor(Tissue, levels = tissues20))
tsd = tsh_r %>% filter(ctag %in% c("Dom")) %>% select(-ctag) %>%
    mutate(tag = factor(tag, levels = taglst$Dom))
tsd %>% distinct(gid, n.tis.tot) %>% count(n.tis.tot) %>% summarise(ntot = sum(n))
tz = tsd %>% filter(tag != 'MP', n.tis >= 1) %>% group_by(gid) %>%
    summarise(n.tis.dom = sum(n.tis)) %>%
    ungroup() %>% count(n.tis.dom)
tz %>% summarise(ntot = sum(n))
tags.mix.d = c("consis. MP", "consis. LP/BLP", "consis. HP/AHP", "mix of HP/LP")
txd = tsd %>% filter(n.tis.tot >= 0 * ntissue) %>%
    select(-prop.tis) %>%
    mutate(tag = as.character(tag)) %>%
    mutate(tag = ifelse(tag %in% c("LP", "BLP"), "LP", tag)) %>%
    mutate(tag = ifelse(tag %in% c("HP", "AHP"), "HP", tag)) %>%
    mutate(tag = ifelse(tag %in% c("PD_L", "PD_H"), "MP", tag)) %>%
    group_by(gid, n.tis.tot, tag) %>%
    summarise(n.tis = sum(n.tis)) %>% ungroup() %>%
    spread(tag, n.tis) %>%
    replace_na(list(LP=0,HP=0,MP=0)) %>%
    mutate(tag = ifelse(LP > 0 & HP > 0, 'mix of HP/LP',
                 ifelse(LP > 0, 'consis. LP/BLP',
                 ifelse(HP > 0, 'consis. HP/AHP', 'consis. MP')))) %>%
    mutate(tag = factor(tag, levels = tags.mix.d))
txd %>% filter(n.tis.tot >= 5) %>% filter(HP >= 2)
txd %>% filter(n.tis.tot >= 5) %>% filter(LP >= 2)
#}}}

#{{{ p1: per-tissue stats
tz = tnd %>% group_by(Tissue) %>%
    summarise(ng.tot = sum(ngene), ng.add = ngene[tag=='MP'],
              ng.dom = sum(ngene[tag %in% c("LP","HP")]),
              ng.odom = sum(ngene[tag %in% c("BLP","AHP")])) %>%
    mutate(p.add = ng.add/ng.tot, p.dom = ng.dom/ng.tot, p.odom=ng.odom/ng.tot) %>%
    print(n = 20)
tnd %>% group_by(tag) %>% summarise(n=sum(ngene)) %>% mutate(p=n/sum(n))
#
tz %>% group_by(1) %>%
    summarise(p.add.l=min(p.add), p.add.m=max(p.add),
              p.dom.l=min(p.dom), p.dom.m=max(p.dom),
              p.odom.l=min(p.odom), p.odom.m=max(p.odom))
tp = tnd %>% 
    mutate(tag = factor(tag, levels = taglst$Dom))
tps = tp %>% group_by(Tissue) %>% summarise(ng.tot = sum(ngene)) %>%
    mutate(lab = sprintf("%s (%d)", Tissue, ng.tot)) %>%
    arrange(desc(Tissue))
p1 = ggplot(tp) +
    geom_bar(aes(x = Tissue, y = ngene, fill = tag), position = position_stack(reverse = T), stat = 'identity', width = 0.8) +
    scale_x_discrete(breaks = tps$Tissue, labels = tps$lab, expand = c(0,0)) +
    scale_y_continuous(name = sprintf("Number Genes"), expand = expand_scale(mult=c(0,.03))) +
    scale_fill_manual(values = cols.dom)+
    coord_flip() +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    theme(plot.margin = unit(c(1.5,.5,.5,.5), "lines")) +
    theme(legend.position = c(1,1), legend.justification=c(1,0), legend.background = element_blank()) +
    guides(direction = 'vertical', fill = guide_legend(nrow = 1, byrow = F)) +
    theme(legend.title = element_blank(), legend.key.size = unit(0.8, 'lines'), legend.text = element_text(size = 8)) +
    theme(axis.title.x = element_text(size=8)) +
    theme(axis.title.y = element_blank())
#}}}

#{{{ p2 a-c: composition of MP/LP/HP
tx1 = tm %>% filter(!is.na(pDE) & pDE != 'non-DE', abs(log2mb) >= 1,
                    !is.na(SPE), !is.na(Dom)) %>%
    mutate(dtag = Dom) %>%
    select(Tissue, gid, pDE, log2mb, DoA, SPE, dtag)
sum(is.na(tx1$dtag))
tx2 = tx1 %>%
    filter(Tissue %in% tissues20) %>%
    mutate(Tissue = factor(Tissue, levels = tissues20))
tx = tx2 %>% mutate(SPE = ifelse(SPE == 'non_SPE', SPE, 'SPE')) %>%
    mutate(DE = ifelse(abs(log2mb) < 2, 'DE2-4',
                ifelse(abs(log2mb) < 3, 'DE4-8', 'DE8+')))
tp0 = tx %>% mutate(tag = 'All (DE 2+)') %>% select(tag, Tissue, gid, dtag, DoA)
tp1 = tx %>% mutate(tag = DE) %>% select(tag, Tissue, gid, dtag, DoA)
tp2 = tx %>% mutate(tag = SPE) %>% select(tag, Tissue, gid, dtag, DoA)
tp3 = tx %>% left_join(tsyn, by = 'gid') %>%
    replace_na(list(atag='non-syntenic')) %>%
    mutate(tag = atag) %>% select(tag, Tissue, gid, dtag, DoA)
tp4 = tx %>% left_join(ttf, by = 'gid') %>%
    replace_na(list(atag='non-TF')) %>%
    mutate(tag = atag) %>% select(tag, Tissue, gid, dtag, DoA)
tags0 = c("All (DE 2+)")
tags1 = c("DE1-2", "DE2-4", "DE4-8", "DE8+")
tags2 = c("SPE", "non_SPE")
tags3 = c("syntenic", "non-syntenic")
tags4 = c("TF", "non-TF")
tags = c(tags0, tags1, tags2)
#tags = c(tags0, tags1, tags2, tags3, tags4)
tp11 = rbind(tp0, tp1) %>% mutate(tag = factor(tag, levels = c(tags0, tags1)))
tp21 = rbind(tp0, tp2) %>% mutate(tag = factor(tag, levels = c(tags0, tags2)))
tp31 = rbind(tp0, tp3) %>% mutate(tag = factor(tag, levels = c(tags0, tags3)))
tp41 = rbind(tp0, tp4) %>% mutate(tag = factor(tag, levels = c(tags0, tags4)))
#
lstd = list(tp11, tp21)
#lstd = list(tp31, tp41)
doa_min = -3
doa_max = 3
for (i in 1:length(lstd)) {
    tp = lstd[[i]] %>%
        mutate(DoA = ifelse(DoA < doa_min, doa_min, DoA)) %>%
        mutate(DoA = ifelse(DoA > doa_max, doa_max, DoA))
    pr = ggplot(tp) +
        geom_density_ridges_gradient(aes(x = DoA, y = tag, fill = ..x..)) +
        scale_x_continuous(limits = c(doa_min,doa_max), expand = c(.05,0)) +
        scale_y_discrete(expand = c(.01,0)) +
        scale_fill_viridis() +
        theme_bw() +
        theme(legend.position = 'none') +
        theme(axis.ticks.length = unit(0, 'lines')) +
        theme(plot.margin = unit(c(.3,.1,.3,.1), "lines")) +
        theme(axis.title.x = element_text(size = 9)) +
        theme(axis.title.y = element_blank()) +
        theme(axis.text.x = element_text(size = 8, color = "black")) +
        theme(axis.text.y = element_text(size = 8, color = 'black'))
    p = ggplot(tp) +
        geom_density(aes(x = DoA, color = tag), alpha = .3) +
        scale_x_continuous(name = 'Scaled difference (btw. hybrid and mid-parent)', limits = c(doa_min, doa_max), expand = c(0.01,0)) +
        scale_y_continuous(name = 'Density', expand = c(.01,0)) +
        #scale_fill_manual(values = cols.dom) +
        scale_fill_locuszoom() +
        theme_bw() +
        theme(legend.position = c(1,1), legend.justification=c(1,1)) +
        guides(direction = 'horizontal', fill = guide_legend(ncol = 1, byrow = F)) +
        theme(legend.title = element_blank(), legend.key.size = unit(0.8, 'lines'), legend.text = element_text(size = 8)) +
        theme(legend.margin = margin(t = .1, r = .5, unit = 'lines')) +
        theme(panel.grid.minor = element_blank()) +
        theme(legend.background = element_blank()) +
        theme(plot.margin = unit(c(.5,.5,.5,.5), "lines")) +
        theme(axis.title = element_text(size = 9))
    if(i == 1) p2a = p
    if(i == 2) p2b = p
}
#
tp = rbind(tp0, tp1, tp2) %>%
    count(tag, dtag) %>%
    mutate(tag = factor(tag, levels = rev(tags))) %>%
    mutate(dtag = factor(dtag, levels = taglst$Dom))
tpx = tp %>% group_by(tag) %>% summarise(n.tot = sum(n)) %>%
    mutate(lab = sprintf("%s (%d)", tag, n.tot)) %>% arrange(tag)
tp = tp %>% left_join(tpx, by = 'tag') %>% mutate(prop = n/n.tot) 
tp %>% select(-n, -n.tot) %>% spread(dtag, prop)
tp = tp %>% filter(tag != 'non_SPE')
p2c = ggplot(tp) +
    geom_bar(aes(x = tag, y = prop, fill = dtag), position = position_stack(reverse = T), stat = 'identity', width = 0.8) +
    scale_x_discrete(breaks = tpx$tag, labels = tpx$lab, expand = c(0,0)) +
    scale_y_continuous(name = sprintf("Proportion Genes"), breaks = seq(0,1,by=0.25), expand=c(0,0)) +
    coord_flip() +
    scale_fill_manual(values = cols.dom) +
    theme_bw() +
    theme(legend.position = c(1,1), legend.justification=c(1,0), legend.background = element_blank()) +
    guides(direction = 'horizontal', fill = guide_legend(nrow = 1, byrow = F)) +
    theme(legend.title = element_blank(), legend.key.size = unit(0.8, 'lines'), legend.text = element_text(size = 8)) +
    theme(legend.margin = margin(b = .2, unit = 'lines')) +
    theme(panel.grid = element_blank()) +
    theme(plot.margin = unit(c(1.5,.5,.5,.5), "lines")) +
    theme(axis.title.x = element_text(size = 9)) +
    theme(axis.title.y = element_blank())
#}}}

#{{{ p3: #non-MP genes per DE category
dtags = c("HP","AHP",'LP','BLP')
dtag = "HP/AHP/LP/BLP"
tp = tsd %>% mutate(tag = as.character(tag)) %>%
    group_by(gid, n.tis.tot) %>%
    summarise(tag = ifelse(sum(dtags %in% tag), dtag, 'MP')) %>% ungroup() %>%
    count(n.tis.tot, tag) %>%
    mutate(fac = ifelse(n.tis.tot <= 2, 1, 2))
tps = tp %>% group_by(n.tis.tot) %>% summarise(ng.tot = sum(n)) %>%
    mutate(lab = sprintf("%2d (%4d)", n.tis.tot, ng.tot))
p = ggplot(tp) +
    geom_bar(aes(x = n.tis.tot, y = n, fill = tag), position = 'stack', stat='identity', width=0.7) +
    scale_x_continuous(name = 'Number DE Tissues', breaks = c(1,2,5,10,15,20), expand = c(0.01,0)) +
    scale_y_continuous(name = "Number Genes", expand = c(0.02, 0)) +
    scale_fill_d3() +
    facet_wrap(~fac, scales='free') +
    theme_bw() +
    theme(legend.position = c(1,1), legend.justification=c(1,1), legend.background = element_blank()) +
    guides(direction = 'vertical', fill = guide_legend(ncol = 1, byrow = F)) +
    theme(legend.title = element_blank(), legend.key.size = unit(.8, 'lines'), legend.text = element_text(size = 8)) +
    theme(panel.border = element_blank(), panel.grid = element_blank()) +
    theme(strip.background = element_blank(), strip.text.x = element_blank()) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
    theme(plot.margin = unit(c(.5,.5,.5,1.5), "lines")) +
    theme(axis.title.x = element_text(size = 9)) +
    theme(axis.title.y = element_text(size = 9)) +
    theme(axis.text.x = element_text(size = 8, color = "black")) +
    theme(axis.text.y = element_text(size = 8, color = 'black'))
gt = ggplotGrob(p)
N <- tp %>% group_by(fac) %>% summarise(count = length(unique(n.tis.tot))) %>% `[[`(2)
panelI <- gt$layout$l[grepl("panel", gt$layout$name)]
gt$widths[panelI] <- unit(N, "null")
gt$widths[panelI[1] + 1] = unit(0.05, "cm")
p3 = gt
#}}}

#{{{ p4: #non-MP genes broken by #DE tissues
dtags = c("HP","AHP",'LP','BLP')
tp = tsd %>% filter(tag %in% dtags, n.tis >= 1) %>% group_by(gid, n.tis.tot) %>%
    summarise(n.tis.dom = sum(n.tis)) %>%
    ungroup() %>% 
    group_by(n.tis.dom, n.tis.tot) %>%
    summarise(ng = n()) %>%
    mutate(fac = ifelse(n.tis.dom <= 3, 1, 2)) %>%
    mutate(tag = sprintf("%2d", n.tis.tot))
tps = tp %>% group_by(n.tis.dom) %>% summarise(ng.tot = sum(ng)) %>%
    mutate(lab = sprintf("%2d (%4d)", n.tis.dom, ng.tot))
sum(tp %>% filter(n.tis.dom >= 5) %>% pull(ng))
sum(tp %>% filter(n.tis.dom > 5) %>% pull(ng))/sum(tp %>% filter(n.tis.dom > 0) %>% pull(ng))
tsd %>% filter(tag %in% c("LP",'BLP'), n.tis >= 2) %>% count(1)
tsd %>% filter(tag %in% c("HP",'AHP'), n.tis >= 2) %>% count(1)
tsd %>% filter(tag %in% dtags) %>% group_by(gid) %>% summarise(n.tis = sum(n.tis)) %>% filter(n.tis >= 2) %>% count(1)
p = ggplot(tp) +
    geom_bar(aes(x = n.tis.dom, y = ng, fill = n.tis.tot), position = 'stack', stat='identity', width=0.7) +
    scale_x_continuous(name = 'Number non-MP Tissues', breaks = 1:14, expand = c(0.01,0)) +
    scale_y_continuous(name = "Number Genes", expand = c(0.02, 0)) +
    scale_fill_viridis(name = 'Number DE Tissues', breaks = c(1,2,5,10,15,20), direction = -1) +
    facet_wrap(~fac, scales='free') +
    theme_bw() +
    theme(legend.position = c(1,1), legend.justification=c(1,1), legend.background = element_blank()) +
    guides(direction = 'vertical', fill = guide_legend(ncol = 1, byrow = F)) +
    theme(legend.title = element_text(size=8), legend.key.size = unit(.8, 'lines'), legend.text = element_text(size = 8)) +
    theme(panel.border = element_blank(), panel.grid = element_blank()) +
    theme(strip.background = element_blank(), strip.text.x = element_blank()) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
    theme(plot.margin = unit(c(.5,.5,.5,1.5), "lines")) +
    theme(axis.title = element_text(size = 9)) +
    theme(axis.text = element_text(size = 8)) 
gt = ggplotGrob(p)
N <- tp %>% group_by(fac) %>% summarise(count = length(unique(n.tis.dom))) %>% `[[`(2)
panelI <- gt$layout$l[grepl("panel", gt$layout$name)]
gt$widths[panelI] <- unit(N, "null")
gt$widths[panelI[1] + 1] = unit(0.05, "cm")
p4 = gt
#}}}

#{{{ p5 a-b: tissue specificity
de.prop = .25
tsd %>% filter(n.tis.tot >= de.prop * ntissue) %>% distinct(gid) %>% count(1)
tp = tsd %>% filter(n.tis.tot >= de.prop * ntissue) %>%
    filter(! tag %in% c("PD_H", "PD_L"))
tp %>% group_by(tag) %>% summarise(prop.tis.mean = mean(prop.tis))
tps = tp %>% count(tag) %>%
    mutate(lab = sprintf("%s (%d)", tag, n)) %>%
    arrange(desc(tag))
p5a = ggplot(tp) +
    geom_boxplot(aes(x = tag, y = prop.tis, fill = tag), notch = T, outlier.size = 1, width = .8) +
    scale_x_discrete(breaks = tps$tag, labels = tps$lab, expand = c(0,.5)) +
    scale_y_continuous(name = 'Proportion (DE) Tissues', limits = c(0,1), expand = c(0,0)) +
    scale_fill_startrek() +
    coord_flip() +
    theme_bw() +
    theme(legend.position = 'none') +
    theme(plot.margin = unit(c(.5,.5,.5,1), "lines")) +
    theme(axis.title.x = element_text(size = 9)) +
    theme(axis.title.y = element_blank()) +
    theme(axis.text = element_text(size = 8))
#
tp = tsd %>% filter(n.tis.tot >= de.prop * ntissue)  %>%
    count(tag, tsTag)
tps = tp %>% group_by(tag) %>% summarise(n = sum(n)) %>% 
    mutate(lab = sprintf("%s (%d)", tag, n))
tp %>% group_by(tag) %>%
    summarise(p.rare = n[tsTag == 'Tissue specific']/sum(n)) %>%
    ungroup()
sum(tps %>% filter(tag %in% c("BLP", "LP")) %>% select(n))
sum(tps %>% filter(tag %in% c("AHP", "HP")) %>% select(n))
p5b = ggplot(tp) +
    geom_bar(aes(x=tag, y=n, fill=tsTag), position=position_fill(reverse=T), width=0.8) +
    scale_x_discrete(breaks = tps$tag, labels = tps$lab, expand = c(0,0)) +
    scale_y_continuous(name = sprintf("Prop. Genes"), breaks = seq(0,1,by=0.25), expand = c(0, 0)) +
    coord_flip() +
    scale_fill_npg() +
    theme_bw() +
    theme(element_blank(), panel.grid = element_blank()) +
    theme(legend.position = c(1,1),  legend.justification=c(1,0), legend.background = element_blank()) +
    guides(direction = 'horizontal', fill = guide_legend(nrow = 1, byrow = F)) +
    theme(legend.title = element_blank(), legend.key.size = unit(0.8, 'lines'), legend.text = element_text(size = 8)) +
    theme(legend.margin = margin(b = .2, unit = 'lines')) +
    theme(panel.grid = element_blank()) +
    theme(plot.margin = unit(c(1.5,.5,.5,.5), "lines")) +
    theme(axis.title.x = element_text(size = 9)) +
    theme(axis.title.y = element_blank()) +
    theme(axis.text = element_text(size = 8))
#}}}

#{{{ p6 a-b: #DE Tissues ~ #LP/BLP | HP/AHP Tissues
for (dtag in names(listd)) {
    title.p = sprintf("Num. %s Tissues", dtag)
    tp = tsd %>% group_by(gid) %>%
        summarise(n.tis.tot = n.tis.tot[1],
                  n.tis = sum(n.tis[tag %in% listd[[dtag]]])) %>%
        ungroup() %>%
        mutate(tag = sprintf("%2d", n.tis)) %>%
        count(n.tis.tot, tag)
    print(tp %>% filter(as.numeric(tag) >= 1) %>% summarise(ng = sum(n)))
    tps = tp %>% group_by(n.tis.tot) %>%
        summarise(ng.tot = sum(n)) %>%
        ungroup() %>%
        mutate(lab = sprintf("%2d (%4d)", n.tis.tot, ng.tot))
    p = ggplot(tp) +
        geom_bar(aes(x = n.tis.tot, y = n, fill = tag), position = position_fill(reverse = T), stat='identity', width=0.8) +
        scale_x_continuous(name = "Number DE Tissue", breaks = tps$n.tis.tot, labels = tps$lab, expand = c(0,0)) +
        scale_y_continuous(name = "Proportion Genes", breaks = c(.25,.5,.75), expand = c(0,0)) +
        coord_flip() +
        scale_fill_simpsons(name = title.p) +
        theme_bw() +
        theme(legend.position = c(1,1), legend.justification=c(1,0), legend.background = element_blank()) +
        guides(direction = 'horizontal', fill = guide_legend(title.hjust = 0, nrow = 2, byrow = T)) +
        theme(legend.title = element_text(size=9), legend.key.size = unit(0.7, 'lines'), legend.text = element_text(size = 8)) +
        theme(legend.margin = margin(b = .2, unit='lines')) +
        theme(panel.grid = element_blank()) +
        #theme(axis.ticks.length = unit(0, 'lines')) +
        theme(plot.margin = unit(c(2.5,.5,.5,.5), "lines")) +
        theme(axis.title = element_text(size = 9)) +
        theme(axis.text = element_text(size = 8))
    if(dtag == 'LP/BLP') p6a = p
    if(dtag == 'HP/AHP') p6b = p
}
#}}}

#{{{ p7: mix LP/HP
de.prop = 0
txd %>% filter(HP+LP>=2) %>% count(tag) %>% mutate(p = n/sum(n))
txd %>% filter(HP+LP>=2) %>% summarise(ntot = n())
tp = txd %>% mutate(n.tis = LP+HP) %>% filter(n.tis >= 1) %>%
    count(n.tis, tag) %>%
    mutate(fac = ifelse(n.tis <= 2, 1, 2))
xlabel = sprintf("Number non-MP Tissues")
tp2 = tp %>% filter(n.tis >= 2) %>%
    group_by(tag) %>% summarise(n = sum(n)) %>% mutate(p = n/sum(n)) %>% print(n=5)
p = ggplot(tp) +
    geom_bar(aes(x = n.tis, y = n, fill = tag), position = 'stack', stat='identity', width=0.7) +
    scale_x_continuous(name = xlabel, breaks = 1:10, expand = c(.01,0)) +
    scale_y_continuous(name = "Number Genes", expand = c(.02, 0)) +
    scale_fill_manual(values = cols.mix3) +
    facet_wrap(~fac, scales='free') +
    theme_bw() +
    theme(legend.position = c(1,1), legend.justification=c(1,1)) +
    guides(direction = 'vertical', fill = guide_legend(ncol = 1, byrow = F)) +
    theme(legend.title = element_blank(), legend.key.size = unit(.8, 'lines'), legend.text = element_text(size = 8)) +
    theme(panel.border = element_blank(), panel.grid = element_blank()) +
    theme(strip.background = element_blank(), strip.text.x = element_blank()) +
    theme(plot.margin = unit(c(.5,.5,.5,.5), "lines")) +
    theme(axis.title = element_text(size = 9)) +
    theme(axis.text = element_text(size = 8))
gt = ggplotGrob(p)
N <- tp %>% group_by(fac) %>% summarise(count = length(unique(n.tis))) %>% `[[`(2)
panelI <- gt$layout$l[grepl("panel", gt$layout$name)]
gt$widths[panelI] <- unit(N, "null")
gt$widths[panelI[1] + 1] = unit(0.05, "cm")
p7a = gt
#
tp = txd %>% filter(n.tis.tot >= de.prop * ntissue) %>%
    count(n.tis.tot, tag) 
tps = tp %>% group_by(n.tis.tot) %>%
    summarise(ng.tot = sum(n)) %>%
    ungroup() %>%
    mutate(lab = sprintf("%2d (%4d)", n.tis.tot, ng.tot))
p7 = ggplot(tp) +
    geom_bar(aes(x = n.tis.tot, y = n, fill = tag), position = position_fill(reverse = T), stat='identity', width=0.8) +
    scale_x_continuous(name = "Number DE Tissue", breaks = tps$n.tis.tot, labels = tps$lab, expand = c(0,0)) +
    scale_y_continuous(name = "Proportion Genes", expand = c(0, 0)) +
    scale_fill_manual(values = cols.mix4) +
    coord_flip() +
    theme_bw() +
    theme(legend.position = c(1,1), legend.justification=c(1,0), legend.background = element_blank()) +
    guides(direction = 'horizontal', fill = guide_legend(nrow = 1, byrow = F)) +
    theme(legend.title = element_blank(), legend.key.size = unit(.8, 'lines'), legend.text = element_text(size = 8)) +
    theme(legend.margin = margin(b=0.2, unit = 'lines')) +
    theme(panel.grid = element_blank()) +
    #theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(1,.5,.5,.5), "lines")) +
    theme(axis.title = element_text(size = 9)) +
    theme(axis.text = element_text(size = 8))
#}}}

#{{{ enrichment
tge = tm %>% filter(Tissue %in% tissues20, !is.na(pDE) & pDE != 'non-DE') %>%
    filter(!is.na(Dom) & Dom %in% c("LP", "BLP", "HP", "AHP")) %>%
    select(Tissue, gid, Dom) %>%
    mutate(Dom = as.character(Dom)) %>%
    mutate(Dom = ifelse(Dom == 'AHP', 'HP', Dom)) %>%
    mutate(Dom = ifelse(Dom == 'BLP', 'LP', Dom)) %>%
    transmute(tag = sprintf("%s %s", Tissue, Dom), gid = gid)
te = go_enrich_genesets(tge, pval.cutoff = 1e-4)
fo = file.path(dirw, "51.go.enrich/21.dom.tsv")
write_tsv(te, fo)

tge = txd %>% filter(n.tis.tot >= 2) %>% select(tag, gid)
te = go_enrich_genesets(tge, pval.cutoff = 1e-4)
fo = file.path(dirw, "51.go.enrich/22.dom.mix.tsv")
write_tsv(te, fo)
#}}}

#{{{ plot
fo = file.path(dirw, "25.Dom.pdf")
ggarrange(p2a, p2b, p2c, p5a,
    nrow = 2, ncol = 2, labels = LETTERS[1:4], heights = c(2,2)) %>%
    ggexport(filename = fo, width = 8, height = 4)
#
fo = file.path(dirw, "25.Dom.sup.pdf")
ggarrange(p1, p4, p6a, p6b, p7, p7a,
    nrow = 3, ncol = 2, labels = LETTERS[1:6], heights = c(2.5,2,2)) %>%
    ggexport(filename = fo, width = 8, height = 9)
#}}}

#{{{# DoA heatmap
de.prop = 0.5
gids = tsd %>% filter(tag != 'MP', n.tis >= 1) %>% group_by(gid) %>%
    summarise(n.tis.dom = sum(n.tis), n.tis.tot = n.tis.tot[1]) %>%
    ungroup() %>%
    filter(n.tis.tot >= de.prop * ntissue, n.tis.dom >= 2) %>%
    pull(gid)
length(gids)
tkl = tm %>% filter(gid %in% gids, Tissue %in% tissues20) %>%
    select(gid, Tissue, DoA)
#tk[is.na(tk)] = 0
describe(tkl$DoA)
tkl = tkl %>%
    mutate(DoA = ifelse(!is.na(DoA) & DoA < -2, -2, DoA)) %>%
    mutate(DoA = ifelse(!is.na(DoA) & DoA > 2, 2, DoA))
describe(tkl$DoA)
tkw = tkl %>% spread(Tissue, DoA)
dis = daisy(tkw[,-1], metric = "gower")
dis[is.na(dis)] = median(dis, na.rm = T)
hc = hclust(dis, method = "ward.D")
gidsO = tkw$gid[hc$order]
#
tkl = tkl %>% mutate(gid = factor(gid, levels = gidsO),
                     Tissue = factor(Tissue, levels = tissues20))
p1 = ggplot(tkl) +
    geom_tile(aes(x = Tissue, y = gid, fill = DoA)) + 
    #scale_x_discrete(name = '') +
    scale_y_discrete(name = 'Genes') +
    scale_fill_gradient2(breaks=c(-1,0,1), labels=c("F1=LP", "F1=MP", "F1=HP"), na.value = "gray50") + 
    theme_bw() +
    theme(plot.margin = unit(c(.3,.3,.3,.3), "lines")) +
    theme(panel.border = element_blank(), panel.grid = element_blank()) +
    theme(legend.position = 'right', legend.direction = "vertical", legend.justification = c(.5,.5)) +
    theme(legend.title = element_blank(), legend.key.size = unit(0.5, 'lines')) +
    theme(axis.ticks.length = unit(0, 'lines')) +
    theme(axis.title.x = element_blank()) +
    theme(axis.title.y = element_blank()) +
    theme(axis.text.x = element_text(size = 8, colour = "black", angle = 70, hjust = 1)) +
    theme(axis.text.y = element_blank())
fo = sprintf("%s/25.Dom.heatmap.pdf", dirw)
ggsave(p1, filename = fo, width = 4, height = 9)
#}}}

#{{{# Dom reproducibility, tissue distance
tp = tdt %>% filter(t1 %in% tissues20, t2 %in% tissues20) %>%
    filter(t1 != t2) %>%
    mutate(lab = 'bg_random') %>% select(lab, distance)
summary(tp$distance)
for (dtag in names(listd)) {
    gids = tsd %>% filter(tag %in% listd[[dtag]]) %>%
        group_by(gid) %>% summarise(n.tis = sum(n.tis)) %>% ungroup() %>%
        filter(n.tis == 2) %>%
        distinct(gid) %>% pull(gid)
    tx = tm %>% filter(gid %in% gids, Dom %in% listd[[dtag]], Tissue %in% tissues20) %>%
        select(Tissue, gid, Dom) %>%
        distinct(Tissue, gid)
    table(tx %>% count(gid) %>% pull(n))
    ty = tx %>% mutate(Tissue = as.character(Tissue)) %>%
        group_by(gid) %>%
        summarise(t1 = Tissue[1], t2 = Tissue[2]) %>% ungroup() %>%
        left_join(tdt, by = c('t1', 't2')) %>%
        mutate(lab = sprintf("%s (%d)", dtag, length(gids))) %>%
        select(lab, distance)
    tp = rbind(tp, ty)
}
for (dtag in c("DE_B", "DE_M")) {
    gids = tsh_d %>% filter(ctag == 'pDE') %>%
        filter(tag == dtag, n.tis == 2) %>% distinct(gid) %>% pull(gid)
    tx = tm %>% filter(gid %in% gids, pDE2 %in% dtag) %>%
        distinct(Tissue, gid)
    table(tx %>% count(gid) %>% pull(n))
    ty = tx %>% mutate(Tissue = as.character(Tissue)) %>%
        group_by(gid) %>%
        summarise(t1 = Tissue[1], t2 = Tissue[2]) %>% ungroup() %>%
        left_join(tdt, by = c('t1', 't2')) %>%
        mutate(lab = sprintf("%s (%d)", dtag, length(gids))) %>%
        select(lab, distance)
    tp = rbind(tp, ty)
}

p1 = ggplot(tp) +
    geom_density(aes(x = distance, color = lab)) +
    scale_x_continuous(name = 'Mean Pairwise Tissue Distance', expand = c(0.01,0)) +
    scale_y_continuous(expand = c(0.01, 0)) +
    scale_color_npg() +
    theme_bw() +
    theme(panel.border = element_blank()) +
    #theme(panel.grid = element_blank()) +
    theme(legend.position = c(1,1),  legend.justification=c(1,1)) +
    guides(direction = 'vertical', fill = guide_legend(nrow = 3, byrow = F)) +
    theme(legend.title = element_blank(), legend.key.size = unit(0.8, 'lines'), legend.text = element_text(size = 8)) +
    #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
    theme(plot.margin = unit(c(.5,.5,.5,.5), "lines")) +
    theme(axis.title.x = element_text(size = 9)) +
    theme(axis.title.y = element_blank()) +
    theme(axis.text.x = element_text(size = 8, color = "black")) +
    theme(axis.text.y = element_text(size = 8, color = 'black')) 
fo = file.path(dirw, "test.dom.pdf")
ggsave(p1, file=fo, width=6, height=6)
#}}}
#}}}

#{{{ sum of Reg1
#{{{
ng.all = length(unique(tm$gid))
tm %>% group_by(Tissue) %>%
    summarise(n = sum(!is.na(prop.h)), n.exp = sum(silent == 1)) %>%
    mutate(p = n/n.exp) %>%
    group_by(1) %>% 
    summarise(ng.min = min(n), ng.max = max(n),
              p.min = min(p), p.max = max(p))
ntissue = max(tsh_r$n.tis.tot)
regs = taglst$Reg1
tnr = t_num %>% filter(ctag %in% c("Reg1")) %>%
    select(-ctag) %>%
    filter(Tissue %in% tissues20) %>%
    mutate(Tissue = factor(Tissue, levels = tissues20)) %>%
    mutate(tag = factor(tag, levels = taglst$Reg1))
tsr = tsh_r %>% filter(ctag %in% c("Reg1")) %>% select(-ctag)
#}}}

#{{{ p1 a-b: root case study
require(ggExtra)
tissue = 'root_0DAP'
tp1 = tm %>% filter(Tissue == tissue, !is.na(Reg1)) %>%
    mutate(tag = Reg1) 
tp2 = tm %>% filter(Tissue == tissue, !is.na(Reg2)) %>%
    mutate(tag = Reg2)
pa = ggplot(tp1) +
    geom_point(aes(x = prop.p, y = prop.h, color = tag), size = .5) +
    geom_vline(xintercept = .5) +
    geom_hline(yintercept = .5) +
    geom_abline(intercept = 0, slope = 1) +
    scale_x_continuous(name = 'Parental B73 Proportion', expand = c(.01, 0), limits = c(0,1)) +
    scale_y_continuous(name = 'F1 B73 Proportion', expand = c(.01, 0)) +
    scale_color_manual(values = cols.reg1) +
    theme_bw() +
    theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(4,3,.5,.5), "lines")) +
    theme(legend.position = c(1,1),  legend.justification=c(.1,0)) +
    theme(legend.key.size = unit(1, 'lines'), legend.key.height = unit(.5, 'lines')) +
    theme(legend.title = element_blank(), legend.text = element_text(size = 8)) +
    guides(color = guide_legend(override.aes = list(size=2.5), ncol = 1, byrow = F)) +
    theme(axis.title = element_text(size = 9)) +
    theme(axis.text = element_text(size = 8))
pb = ggplot(tp2) +
    geom_point(aes(x = prop.p, y = prop.h, color = tag), size = .5) +
    geom_vline(xintercept = .5) +
    geom_hline(yintercept = .5) +
    geom_abline(intercept = 0, slope = 1) +
    scale_x_continuous(name = 'Parental B73 Proportion', expand = c(.01, 0), limits = c(0,1)) +
    scale_y_continuous(name = 'F1 B73 Proportion', expand = c(.01, 0)) +
    scale_color_manual(values = cols.reg2) +
    theme_bw() +
    theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(4,3,.5,.5), "lines")) +
    theme(legend.position = c(1,1),  legend.justification=c(0.1,0)) +
    theme(legend.key.size = unit(1, 'lines'), legend.key.height = unit(1, 'lines')) +
    theme(legend.title = element_blank(), legend.text = element_text(size = 8)) +
    guides(color = guide_legend(override.aes = list(size=2.5), ncol = 1, byrow = F)) +
    theme(axis.title = element_text(size = 9)) +
    theme(axis.text = element_text(size = 8))
p1a = ggMarginal(pa, type = 'histogram', margins = 'both', size = 3,
                groupColour = T, groupFill = T,
                xparams = list(size=0), yparams = list(size=0))
p1b = ggMarginal(pb, type = 'histogram', margins = 'both', size = 3,
                groupColour = T, groupFill = T,
                xparams = list(size=0), yparams = list(size=0))
#}}}

#{{{ p2: per-tissue stats
tnr1 = tnr %>% group_by(Tissue) %>% summarise(ng.tot = sum(ngene))
tnr %>%
    left_join(tnr1, by = 'Tissue') %>%
    mutate(pg = ngene/ng.tot) %>%
    group_by(tag) %>%
    summarise(ng.de.min = min(ng.tot),
              ng.de.max = max(ng.tot),
              ng.min = min(ngene),
              ng.max = max(ngene),
              pg.min = min(pg),
              pg.max = max(pg))
tp = tnr
tps = tp %>% group_by(Tissue) %>% summarise(ng.tot = sum(ngene)) %>%
    mutate(lab = sprintf("%s (%d)", Tissue, ng.tot)) %>%
    arrange(desc(Tissue))
p2 = ggplot(tp) +
    geom_bar(aes(x = Tissue, y = ngene, fill = tag), position = position_fill(reverse = T), stat = 'identity', alpha = 0.8, width = 0.8) +
    scale_x_discrete(breaks = tps$Tissue, labels = tps$lab, expand = c(0,0)) +
    scale_y_continuous(name = sprintf("Proportion Genes"), breaks = c(.25,.5,.75), expand = c(0,0)) +
    scale_fill_manual(values = cols.reg1) +
    coord_flip() +
    theme_bw() +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank()) +
    theme(plot.margin = unit(c(.5,.5,.5,.5), "lines")) +
    theme(legend.position = c(.5,1),  legend.justification=c(.5,0)) +
    guides(direction = 'vertical', fill = guide_legend(nrow = 1, byrow = F)) +
    theme(legend.title = element_blank(), legend.key.size = unit(0.8, 'lines'), legend.text = element_text(size = 8)) +
    theme(axis.title.x = element_text(size = 9)) +
    theme(axis.title.y = element_blank()) +
    theme(axis.text = element_text(size = 8))
#}}}

#{{{ p3 a-e: cis/trans composition
#{{{
tx1 = tm %>%
    filter(!is.na(pDE) & pDE != 'non_DE', abs(log2mb) > 1,
           !is.na(SPE), !is.na(Dom), !is.na(Reg1)) %>%
    filter(Tissue %in% tissues20) %>%
    mutate(Tissue = factor(Tissue, levels = tissues20)) %>%
    mutate(dtag = Dom, rtag = Reg1, bic.diff = bic.c - bic.t) %>%
    select(Tissue, gid, pDE, log2mb, SPE, dtag, rtag, prop.p, prop.h, bic.diff)
sum(is.na(tx1$dtag))
tx = tx1 %>% mutate(SPE = ifelse(SPE == 'non_SPE', SPE, 'SPE')) %>%
    mutate(DE = ifelse(abs(log2mb) < 2, 'DE2-4',
                ifelse(abs(log2mb) < 3, 'DE4-8', 'DE8+')))
#}}}
#
#{{{ read Li et al eQTL
fq = file.path(dird, 'li2013/10.eQTL.v4.tsv')
tq = read_tsv(fq, col_types = 'ccccicicii') %>% rename(gid = gid.v4)
types = c("cis", "cis+trans", "trans")
tq2 = tq %>% group_by(gid) %>%
    summarise(ntype = length(unique(type)), type = type[1]) %>%
    mutate(type = ifelse(ntype == 1, type, 'cis+trans')) %>%
    mutate(type = factor(type, levels = rev(types))) %>%
    select(-ntype)
tq2 %>% count(type)
#tissues23
tissues = tissues23[c(17,3,9,18)]
tq = tibble()
for (i in 1:length(tissues)) {
    tissue = tissues[i]
    tq1 = tx %>% filter(Tissue == tissue) %>%
        inner_join(tq2, by = 'gid') %>%
        mutate(ctag = "Li et al eQTL", tag = type) %>%
        select(ctag, tag, Tissue, gid, rtag, bic.diff)
    tq = rbind(tq, tq1)
}
tags1 = c("cis", "cis+trans", "trans")
tp1 = tq %>% filter(Tissue == tissues[1])
#}}}
#
#{{{
tags0 = c("(DE2+)")
tags2 = c("DE2-4", "DE4-8", "DE8+", "SPE")
tags3 = c("BLP","LP","MP","HP","AHP")
tags4 = c("syntenic","non-syntenic","TF","non-TF")
tp0 = tx %>% mutate(ctag = 'All', tag = '(DE2+)') %>% 
    select(ctag, tag, Tissue, gid, rtag, bic.diff)
tp2a = tx %>% mutate(ctag = 'DE', tag = DE) %>% 
    select(ctag, tag, Tissue, gid, rtag, bic.diff)
tp2b = tx %>% mutate(ctag = 'DE', tag = SPE) %>% filter(tag == 'SPE') %>%
    select(ctag, tag, Tissue, gid, rtag, bic.diff)
tp2 = rbind(tp2a, tp2b)
tp3 = tx %>% mutate(ctag = 'Additivity', tag = dtag) %>% 
    select(ctag, tag, Tissue, gid, rtag, bic.diff) %>% 
    mutate(tag = as.character(tag)) %>% filter(tag %in% tags3)
tp4a = tx %>% left_join(tsyn, by = 'gid') %>% 
    replace_na(list(atag='non-syntenic')) %>%
    mutate(ctag = 'misc', tag = atag) %>%
    select(ctag, tag, Tissue, gid, rtag, bic.diff) 
tp4b = tx %>% left_join(ttf, by = 'gid') %>% 
    replace_na(list(atag='non-TF')) %>%
    mutate(ctag = 'misc', tag = atag) %>%
    select(ctag, tag, Tissue, gid, rtag, bic.diff)
tp4 = rbind(tp4a, tp4b)
#
tags = c(tags0, tags1, tags2, tags3)
#tags = c(tags0, tags1, tags2, tags3, tags4)
tr = rbind(tp0, tp1, tp2, tp3) %>%
    mutate(tag = factor(tag, levels = tags))
ctags = tr %>% distinct(ctag) %>% pull(ctag)
tr = tr %>% mutate(ctag = factor(ctag, levels = ctags))
#describe(tp$bic.diff)
bic_min = -15
bic_max = 20
#}}}

#{{{ p3a: proportion composition
tags.r = c(rev(tags0), rev(tags1), rev(tags2), rev(tags3))
#tags.r = c(rev(tags0), rev(tags1), rev(tags2), rev(tags3), rev(tags4))
tp = tr %>%
    mutate(tag = factor(tag, levels = tags.r)) %>%
    count(ctag, tag, rtag)
tps = tp %>% group_by(ctag, tag) %>% summarise(n.tot = sum(n)) %>%
    mutate(lab = sprintf("%s (%d)", tag, n.tot)) %>% arrange(tag)
tp = tp %>% left_join(tps, by = c('ctag','tag')) %>% mutate(prop = n/n.tot)
tp %>% select(-n, -n.tot, -lab) %>% spread(rtag, prop)
p3a = ggplot(tp) +
    geom_bar(aes(x = tag, y = prop, fill = rtag), position = position_stack(reverse = T), stat = 'identity', alpha = .8, width = .8) +
    scale_x_discrete(breaks = tps$tag, labels = tps$lab, expand = c(0,0)) +
    scale_y_continuous(name = sprintf("Proportion Genes"), breaks = c(.25,.5,.75), expand = c(0,0)) +
    scale_fill_manual(values = cols.reg1) +
    coord_flip() +
    facet_grid(ctag~., scale='free', space = 'free') +
    theme_bw() +
    #theme(strip.placement = "outside") +
    theme(strip.text.y = element_text(size = 8, angle = 0)) +
    theme(legend.position = c(.5,1),  legend.justification=c(.5,0)) +
    guides(direction = 'horizontal', fill = guide_legend(nrow = 1, byrow = F)) +
    theme(legend.title = element_blank(), legend.key.size = unit(0.8, 'lines'), legend.text = element_text(size = 8)) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank()) +
    theme(plot.margin = unit(c(.5,.5,.5,1.5), "lines")) +
    theme(axis.title.x = element_text(size = 9)) +
    theme(axis.title.y = element_blank()) +
    theme(axis.text = element_text(size = 8))
#}}}

#{{{ p3b - piechart
tph = tp0 %>% count(rtag) %>%
    mutate(hitInPop = n, popSize = sum(n), prop.baseline = n/popSize) %>%
    select(rtag, prop.baseline, hitInPop, popSize)
p3b = ggplot(tph, aes(x = '', y = prop.baseline, fill = rtag)) + 
    geom_bar(width = 1, alpha = .8, stat = 'identity') + 
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

#{{{ p3c - fold change
tags.r = c(rev(tags1), rev(tags2), rev(tags3))
#tags.r = c(rev(tags1), rev(tags2), rev(tags3), rev(tags4))
tp = tr %>% filter(ctag != 'All') %>%
    mutate(tag = factor(tag, levels = tags.r)) %>%
    group_by(ctag, tag, rtag) %>%
    summarise(hitInSample = n()) %>%
    mutate(sampleSize = sum(hitInSample), prop = hitInSample/sampleSize) %>% 
    ungroup() %>%
    inner_join(tph, by = 'rtag') %>%
    mutate(fc = log2(prop/prop.baseline)) %>%
    mutate(fc = ifelse(fc < -2, -2, fc)) %>%
    mutate(pval = ifelse(prop < prop.baseline,
        phyper(hitInSample, hitInPop, popSize-hitInPop, sampleSize, lower.tail = T),
        phyper(hitInSample-1, hitInPop, popSize-hitInPop, sampleSize, lower.tail = F))) %>%
    mutate(pval.lab = ifelse(pval < .01, ifelse(pval < .001, "**", "*"), "")) %>%
    print(n=10, width=Inf)
p = ggplot(tp) +
    geom_bar(aes(x = tag, y = fc, fill = (as.numeric(tag) %% 2 == 0)), stat = 'identity', alpha = .8, width = .7) +
    geom_text(aes(x = tag, y = ifelse(fc < 0, fc, fc + .25), label = pval.lab), hjust = .5, vjust = .5, angle = 90) +
    geom_hline(yintercept = 0, color = 'gray25') +
    scale_x_discrete(expand = c(0,.5)) +
    scale_y_continuous(name = 'Observed proportion / Expected proportion based on DE2+', expand = c(.03,0), limits =c(-2.2, 2), breaks = -2:3, labels = 2^(-2:3)) +
    coord_flip() +
    scale_fill_manual(values = rep(c('gray40','gray40'),2)) +
    facet_grid(ctag ~ rtag, scale = 'free_y', space = "free") +
    theme_bw() +
    theme(strip.placement = "outside") +
    theme(strip.text.x = element_text(size = 8)) +#, margin=margin(0,0,.1,.5,'lines'))) +
    theme(strip.text.y = element_text(size = 8, angle = 0)) +
    theme(panel.grid.minor = element_blank()) +#, panel.grid.major.y = element_blank()) +
    theme(legend.position = 'none') +
    theme(plot.margin = unit(c(.3,.3,.3,0), "lines")) +
    theme(axis.title.x = element_text(size = 9)) +
    theme(axis.title.y = element_blank()) +
    theme(axis.text = element_text(size = 8))
pg = ggplot_gtable(ggplot_build(p))
idxs = grep(pattern="strip", pg$layout$name)
for (i in 1:5) {
    idx = idxs[i]
    pg$grobs[[idx]]$grobs[[1]]$children[[1]]
    pg$grobs[[idx]]$grobs[[1]]$children[[1]]$gp$fill = cols.reg1[i]
    pg$grobs[[idx]]$grobs[[1]]$children[[1]]$gp$alpha = .9
    pg$grobs[[idx]]$grobs[[1]]$children[[2]]$children[[1]]$gp$col = 'white'
}
p3c = pg
#}}}

#{{{ p3x: Li et al eQTL fold change - 4 tissues
tp = tq %>% mutate(tissue = factor(Tissue, levels = tissues)) %>%
    group_by(tissue, tag, rtag) %>%
    summarise(hitInSample = n()) %>%
    mutate(sampleSize = sum(hitInSample), prop = hitInSample/sampleSize) %>% 
    ungroup() %>%
    inner_join(tph, by = 'rtag') %>%
    mutate(fc = log2(prop/prop.baseline)) %>%
    #mutate(fc = ifelse(fc < -2, -2, fc)) %>%
    mutate(pval = ifelse(prop < prop.baseline,
        phyper(hitInSample, hitInPop, popSize-hitInPop, sampleSize, lower.tail = T),
        phyper(hitInSample-1, hitInPop, popSize-hitInPop, sampleSize, lower.tail = F))) %>%
    mutate(pval.lab = ifelse(pval < .01, ifelse(pval < .001, "**", "*"), "")) %>%
    print(n=10, width=Inf)
tp %>% select(tissue, tag, rtag, fc) %>% spread(rtag, fc)
p = ggplot(tp) +
    geom_bar(aes(x = tag, y = fc, fill = (as.numeric(tag) %% 2 == 0)), stat = 'identity', alpha = .8, width = .7) +
    geom_text(aes(x = tag, y = ifelse(fc < 0, fc, fc + .25), label = pval.lab), hjust = .5, vjust = .5, angle = 90) +
    geom_hline(yintercept = 0, color = 'gray25') +
    scale_x_discrete(name = 'Li et al eQTL', expand = c(0,.5)) + 
    scale_y_continuous(name = 'Observed proportion / Expected proportion based on DE2+', expand = c(.03,0), breaks = -2:3, labels = 2^(-2:3), limits = c(-2.2, 1)) +
    coord_flip() +
    scale_fill_manual(values = rep(c('gray35','gray35'),2)) +
    facet_grid(tissue ~ rtag) +
    theme_bw() +
    theme(strip.placement = "outside") +#, strip.background.y = element_blank()) +
    theme(strip.text.x = element_text(size = 8)) +#, margin=margin(0,0,.1,.5,'lines'))) +
    theme(strip.text.y = element_text(size = 8, angle = 0)) +
    theme(panel.grid.minor = element_blank()) +
    theme(legend.position = 'none') +
    theme(plot.margin = unit(c(.3,.3,.3,0), "lines")) +
    theme(axis.title = element_text(size = 9)) +
    theme(axis.text = element_text(size = 8))
pg = ggplot_gtable(ggplot_build(p))
idxs = grep(pattern="strip", pg$layout$name)
for (i in 1:5) {
    idx = idxs[i]
    pg$grobs[[idx]]$grobs[[1]]$children[[1]]
    pg$grobs[[idx]]$grobs[[1]]$children[[1]]$gp$fill = cols.reg1[i]
    pg$grobs[[idx]]$grobs[[1]]$children[[1]]$gp$alpha = .9
    pg$grobs[[idx]]$grobs[[1]]$children[[2]]$children[[1]]$gp$col = 'white'
}
p3x = pg
#}}}

#{{{ p3d - obsolete: cis-trans density plot
for (i in 1:length(lstd)) {
    tp = lstd[[i]]
    tp = tp %>%
        mutate(tag = factor(tag, levels = rev(levels(tp$tag)))) %>%
        mutate(bic.diff = ifelse(bic.diff < bic_min, bic_min, bic.diff)) %>%
        mutate(bic.diff = ifelse(bic.diff > bic_max, bic_max, bic.diff))
    p = ggplot(tp) +
        geom_density_ridges_gradient(aes(x = bic.diff, y = tag, fill = ..x..)) +
        scale_x_continuous(name = expression(cis %<->% trans), limits = c(bic_min,bic_max), expand = c(0,0)) +
        scale_y_discrete(expand = c(0,0)) +
        scale_fill_viridis() +
        theme_bw() +
        theme(legend.position = 'none') +
        #theme(axis.ticks.length = unit(0, 'lines')) +
        theme(plot.margin = unit(c(.3,.3,.3,.3), "lines")) +
        theme(axis.title.x = element_text(size = 9)) +
        theme(axis.title.y = element_blank()) +
        theme(axis.text.x = element_text(size = 8)) +
        theme(axis.text.y = element_text(size = 8))
    if(i == 1) pa = p
    if(i == 2) pb = p
    if(i == 3) pc = p
}
p3c = ggarrange(pa, pb, pc, nrow = 1, ncol = 3, labels = c('C'))
#}}}
#}}}

#{{{ p4 a-b: tissue specificity
de.prop = .25
tp = tsr %>% filter(n.tis.tot >= de.prop * ntissue) %>%
    mutate(tag = factor(tag, levels = rev(regs)))
tp %>% group_by(tag) %>% summarise(n.tis.median = mean(prop.tis))
tps = tp %>% count(tag) %>% 
    mutate(lab = sprintf("%s (%d)", tag, n))
p4a = ggplot(tp) +
    geom_boxplot(aes(x = tag, y = prop.tis, fill = tag), notch = T, outlier.size = 0.5, width = .8, alpha = .8) +
    scale_x_discrete(breaks = tps$tag, labels = tps$lab, expand = c(0,.5)) +
    scale_y_continuous(name = 'Proportion (DE) Tissues', limits = c(0,1), expand = c(0,0)) +
    scale_fill_manual(values = rev(cols.reg1)) +
    coord_flip() +
    theme_bw() +
    theme(legend.position = 'none') +
    theme(plot.margin = unit(c(.5,.5,.5,.5), "lines")) +
    theme(axis.title.x = element_text(size = 9)) +
    theme(axis.title.y = element_blank()) +
    theme(axis.text.x = element_text(size = 8)) +
    theme(axis.text.y = element_text(size = 8)) 
tp = tsr %>% 
    filter(n.tis.tot >= de.prop * ntissue, n.tis >= 1) %>%
    count(tag, tsTag) %>%
    mutate(tag = factor(tag, levels = rev(regs)))
tp %>% group_by(tag) %>% summarise(ntot=sum(n)) %>%
    left_join(tp[tp$tsTag == 'Constitutive',], by = 'tag') %>%
    mutate(p.con = n/ntot)
tps = tp %>% group_by(tag) %>% summarise(n = sum(n)) %>% 
    mutate(lab = sprintf("%s (%d)", tag, n))
p4b = ggplot(tp) +
    geom_bar(aes(x = tag, y = n, fill = tsTag), position = position_fill(reverse = T), stat = 'identity', width = 0.8) +
    scale_x_discrete(breaks = tps$tag, labels = tps$lab, expand = c(0,0)) +
    scale_y_continuous(name = sprintf("Proportion Genes"), breaks = seq(.25,.75,by=.25), expand = c(0,0)) +
    coord_flip() +
    scale_fill_manual(values = cols.ts) +
    theme_bw() +
    theme(legend.position = c(1,1),  legend.justification=c(1,0), legend.background = element_blank()) +
    guides(direction = 'horizontal', fill = guide_legend(nrow = 1, byrow = F)) +
    theme(legend.title = element_blank(), legend.key.size = unit(0.7, 'lines'), legend.text = element_text(size = 8)) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank()) +
    theme(plot.margin = unit(c(1.5,.5,.5,.5), "lines")) +
    theme(axis.title.x = element_text(size = 9)) +
    theme(axis.title.y = element_blank()) +
    theme(axis.text.x = element_text(size = 8)) +
    theme(axis.text.y = element_text(size = 8))
#}}}

#{{{ p5: non-syn/GO enrichment
tc1 = tsh_e %>% select(gid) %>% 
    mutate(ctag = 'control', tag = 'All genes')
tc2 = tsh_e %>% filter(n.tis >= 10) %>% select(gid) %>% 
    mutate(ctag = 'control', tag = 'All expr. genes')
tc3 = tsh_d %>% filter(tag == 'non-pDE', n.tis.tot >= 10) %>% select(gid) %>%
    mutate(ctag = 'control', tag = "non-DE genes")
tc = rbind(tc1, tc2, tc3) %>% select(ctag, tag, gid)
de.prop = .25
tr = tsh_r %>% filter(n.tis.tot >= de.prop * ntissue, tsTag != 'No data', ctag == 'Reg1') %>%
    arrange(tag, desc(tsTag)) %>%
    transmute(ctag = as.character(tag), tag = as.character(tsTag), gid = gid) %>%
    mutate(tag = ifelse(tag == 'Intermediate frequency', 'Intermediate', tag)) %>%
    mutate(tag = ifelse(tag == 'Tissue specific', 'Tissue-specific', tag))
ctags1 = unique(tr$ctag)
tra = tr %>% mutate(tag = 'All')
trb = tr %>% filter(tag == 'Constitutive')
tr = rbind(tra, trb) %>% 
    mutate(ctag = factor(ctag, levels = ctags1)) %>%
    arrange(ctag, tag) %>%
    mutate(tag = sprintf("%s|%12s", ctag, tag))
tags = c(unique(tc$tag), unique(tr$tag))
td = rbind(tc, tr)
ctags = unique(td$ctag)
td = td %>%
    mutate(tag = factor(tag, levels = rev(tags))) %>%
    mutate(ctag = factor(ctag, levels = ctags))
tds = td %>% count(tag) %>%
    mutate(lab = sprintf("%s (%4d)", tag, n))
tjs = list(tsyn, tdom, thom)
legend.titles = c("Non-syntenic", "w.o. Known Domain", "w.o. Homologs")
for (i in 1:3) {
    tj = tjs[[i]]
    legend.title = sprintf("Genes %s (%%)", legend.titles[i])
    tp = td %>%
        left_join(tj, by = 'gid') %>% group_by(ctag, tag) %>%
        summarise(ngene = n(),
                  n.na = sum(is.na(atag)),
                  p.na = n.na/ngene * 100) %>% ungroup() %>%
        mutate(p.lab = sprintf("%.0f%%", p.na))
    ymax = max(tp$p.na)
    p = ggplot(tp) +
        geom_bar(aes(x = tag, y = p.na, fill = ctag), stat = 'identity', alpha = .8, width = .8) +
        geom_text(aes(x = tag, y = p.na-1, label = p.lab), color = 'white', size = 2.5, hjust = 1) +
        scale_x_discrete(breaks = tds$tag, labels = tds$lab, expand = c(0.02,0)) +
        scale_y_continuous(name = legend.title, limits = c(0,ymax+1), expand = c(0,0)) +
        scale_fill_manual(values = c("gray20", cols.reg1)) +
        coord_flip() +
        theme_bw() +
        theme(plot.margin = unit(c(.5,.5,.5,.5), "lines")) +
        theme(legend.position = 'none') +
        theme(panel.grid = element_blank()) +
        theme(axis.ticks.y = element_blank()) +
        theme(axis.title.x = element_text(size = 8, hjust=.5)) +
        theme(axis.title.y = element_blank()) +
        theme(axis.text.x = element_text(size = 8)) +
        theme(axis.text.y = element_text(size=8, family='mono',face='bold'))
    if(i != 1) p = p + theme(axis.text.y = element_blank())
    if(i == 1) pa = p
    if(i == 2) pb = p + theme(plot.margin = unit(c(.5,0,.5,0), 'lines'))
    if(i == 3) pc = p
}
p5 =  ggarrange(pa, pb, pc, nrow = 1, ncol = 3, widths = c(6,2,2))
#}}}

#{{{ plot main figure
fo = file.path(dirw, "25.Reg.pdf")
ggarrange(
    ggarrange(p1a, p1b, nrow = 1, ncol = 2, labels = c('A','B')),
    ggarrange(p3b, p3c, nrow = 1, ncol = 2, widths = c(1,5), labels = 'C'),
    ggarrange(p4a, p4b, nrow = 1, ncol = 2, labels=c('D','E')),
    nrow = 3, ncol = 1, heights = c(4,3.5,2)) %>%
    ggexport(filename = fo, width = 8.5, height = 9)
#
fo = file.path(dirw, "25.Reg.sup.pdf")
ggarrange(
    ggarrange(p2, p3a, nrow = 1, ncol = 2, labels = c('A','B'), common.legend=T),
    ggarrange(p5, nrow = 1, ncol = 1, labels = 'C'),
    ggarrange(p_con1, NULL, nrow = 1, ncol = 2, labels = "D"),
    nrow = 3, ncol = 1, heights = c(3,2.5,1.5)) %>%
    ggexport(filename = fo, width = 8, height = 8)
#
fo = file.path(dirw, '85.eQTL.pdf')
ggarrange(p3b, p3x,
    nrow = 1, ncol = 2, widths = c(1,6))  %>%
    ggexport(filename = fo, width = 9, height = 3)
#}}}
#}}}

#{{{ sum of Reg2
#{{{
regs2 = c('conserved',"B73 biased","Mo17 biased")
ntissue = max(tsh_r$n.tis.tot) 
tnb = t_num %>% filter(ctag %in% c("Reg2")) %>%
    select(-ctag) %>%
    filter(Tissue %in% tissues20) %>%
    mutate(Tissue = factor(Tissue, levels = tissues20)) %>%
    mutate(tag = factor(tag, levels = regs2))
tnb %>% group_by(Tissue) %>% 
    summarise(pb = 1-ngene[tag=='conserved']/sum(ngene)) %>% ungroup() %>%
    group_by(1) %>% summarise(pb.min = min(pb), pb.max = max(pb))
tnb %>% group_by(tag) %>%
    summarise(min.ng = min(ngene),
              max.ng = max(ngene))
tnb %>% filter(tag != 'conserved') %>% group_by(Tissue) %>%
    summarise(ng = sum(ngene)) %>%
    group_by(1) %>% summarise(ng.min = min(ng), ng.max = max(ng))
tsb = tsh_r %>% filter(ctag %in% c("Reg2"))  %>%
    select(-ctag) %>%
    mutate(tag = factor(tag, levels = regs2))
tags.mix.r2 = c("consis. conserved", 
                "consis. B73-biased",
                "consis. Mo17-biased",
                "mix of B73-biased/Mo17-biased")
txb = tsh_r %>% filter(ctag == 'Reg2') %>%
    select(-ctag, -prop.tis) %>%
    mutate(tag = as.character(tag)) %>%
    mutate(tag = ifelse(tag == 'B73 biased', 'b', tag)) %>%
    mutate(tag = ifelse(tag == 'Mo17 biased', 'm', tag)) %>%
    group_by(gid, n.tis.tot, tag) %>%
    summarise(n.tis = sum(n.tis)) %>% ungroup() %>%
    spread(tag, n.tis) %>%
    replace_na(list(b=0,m=0,conserved=0)) %>% 
    mutate(tag = ifelse(b > 0 & m > 0, 'mix of B73-biased/Mo17-biased',
                 ifelse(b > 0, 'consis. B73-biased', 
                 ifelse(m > 0, 'consis. Mo17-biased', 'consis. conserved')))) %>%
    mutate(tag = factor(tag, levels = tags.mix.r2))
#}}}

#{{{ p1: per tissue stats
tps = tnb %>% group_by(Tissue) %>% summarise(ng.tot = sum(ngene)) %>%
    mutate(lab = sprintf("%s (%d)", Tissue, ng.tot)) %>%
    arrange(desc(Tissue))
p1 = ggplot(tnb) +
    geom_bar(aes(x = Tissue, y = ngene, fill = tag), position = position_fill(reverse = T), stat = 'identity', width = 0.8, alpha = .8) +
    scale_x_discrete(breaks = tps$Tissue, labels = tps$lab, expand = c(0,0)) +
    scale_y_continuous(name = sprintf("Proportion Genes"), breaks = c(.25,.5,.75), expand = c(0,0)) +
    scale_fill_manual(values = cols.reg2) +
    coord_flip() +
    theme_bw() +
    theme(plot.margin = unit(c(1.5,.5,.5,.5), "lines")) +
    theme(legend.position = c(.5,1),  legend.justification=c(.5,0), legend.background = element_blank()) +
    guides(direction = 'vertical', fill = guide_legend(nrow = 1, byrow = F)) +
    theme(legend.title = element_blank(), legend.key.size = unit(0.7, 'lines'), legend.text = element_text(size = 8)) +
    theme(panel.grid = element_blank()) +
    theme(axis.title.x = element_text(size = 9)) +
    theme(axis.title.y = element_blank()) +
    theme(axis.text.x = element_text(size = 8)) +
    theme(axis.text.y = element_text(size = 8))
#}}}

#{{{ p2: composition of Reg2 in hDE
tx1 = tm %>% filter(!is.na(hDE)) %>%
    select(Tissue, gid, hDE, log2fm)
tx1 %>%
    group_by(hDE) %>%
    summarise(mean = mean(log2fm), median=median(log2fm), sd = sd(log2fm))
tx2 = tm %>% filter(!is.na(Reg2)) %>%
    select(Tissue, gid, Reg2)
tx3 = tx1 %>% inner_join(tx2, by = c('Tissue', 'gid'))
tp = tx3 %>% mutate(htag = hDE, btag = Reg2) %>%
    group_by(htag, btag) %>%
    summarise(n = n()) %>% ungroup()
tps = tp %>% group_by(htag) %>% summarise(ntot = sum(n)) %>% ungroup() %>%
    mutate(lab = sprintf("%s (%d)", htag, ntot))
pb = ggplot(tp) +
    geom_bar(aes(x = htag, y = n, fill = btag), position = position_fill(reverse = T), stat = 'identity', width = 0.7, alpha = .8) +
    scale_x_discrete(breaks = tps$htag, labels = tps$lab, expand = c(0,0)) +
    scale_y_continuous(name = sprintf("Proportion Genes"), breaks = c(.25,.5,.75), expand=c(0,0)) +
    coord_flip() +
    scale_fill_manual(values = cols.reg2) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    theme(legend.position = c(.5,1),  legend.justification=c(.5,0), legend.background = element_blank()) +
    guides(direction = 'horizontal', fill = guide_legend(nrow = 1, byrow = F)) +
    theme(legend.title = element_blank(), legend.key.size = unit(0.8, 'lines'), legend.text = element_text(size = 8)) +
    theme(plot.margin = unit(c(1.5,.5,.5,.5), "lines")) +
    theme(axis.title.x = element_text(size = 9)) +
    theme(axis.title.y = element_blank()) +
    theme(axis.text = element_text(size = 8))
p2 = pb
#}}}

#{{{ p3: mix pattern
txb %>% filter(b+m>=2) %>% count(tag) %>% mutate(p = n/sum(n))
txb %>% filter(b+m>=2) %>% summarise(ntot = n())
nde.prop = 0
tp = txb %>% filter(n.tis.tot >= nde.prop * ntissue) %>%
    count(n.tis.tot, tag)
ymax = tp %>% group_by(n.tis.tot) %>% summarise(n = sum(n)) %>%
    ungroup() %>% group_by(1) %>%
    summarise(nmax = max(n)) %>% pull(nmax) + 50
p3 = ggplot(tp) +
    geom_bar(aes(x = n.tis.tot, y = n, fill = tag), position = position_stack(reverse = F), stat='identity', width = .8) +
    scale_x_continuous(name = 'Number non-DE tissues', expand = c(0,0)) +
    scale_y_continuous(name = "Number Genes", expand = c(0,0)) +
    scale_fill_manual(values = cols.mix4) +
    theme_bw() +
    theme(legend.position = c(1,1), legend.direction = 'vertical', legend.justification=c(1,1), legend.background = element_blank()) +
    theme(legend.title = element_blank(), legend.key.size = unit(0.8, 'lines'), legend.text = element_text(size = 8)) +
    theme(panel.border = element_blank(), panel.grid = element_blank()) +
    theme(plot.margin = unit(c(.5,.5,.5,.5), "lines")) +
    theme(axis.title.x = element_text(size = 9)) +
    theme(axis.title.y = element_text(size = 9)) +
    theme(axis.text.x = element_text(size = 8)) +
    theme(axis.text.y = element_text(size = 8))
#}}}

#{{{ p7: consis. Reg2 composition of cis/trans
gids.r = tsh_r %>% filter(ctag == 'Reg1', n.tis.tot >= 3) %>%
    distinct(gid) %>% pull(gid)
gids.b = txb %>% filter(b+m >= 3) %>%
    distinct(gid) %>% pull(gid)
gids = gids.r[gids.r %in% gids.b]
tr1 = tm %>% filter(gid %in% gids, !is.na(Reg1), Tissue %in% tissues20) %>%
    transmute(gid = gid, rtag = Reg1) %>%
    mutate(rtag = factor(rtag, levels = taglst$Reg1)) %>%
    group_by(gid, rtag) %>%
    summarise(n = n()) %>%
    mutate(nc = n/sum(n))
sum(tr1$nc)
tp = txb %>% filter(gid %in% gids)
tps = tp %>% count(tag) %>% mutate(lab = sprintf("%s (%d)", tag, n))
tp = tp %>%
    transmute(gid = gid, tag = tag) %>%
    inner_join(tr1, by = 'gid') %>%
    group_by(tag, rtag) %>%
    summarise(n = sum(nc)) %>%
    mutate(prop = n/sum(n)) %>%
    print(n=20)
sum(tp$n) 
p7 = ggplot(tp) +
    geom_bar(aes(x = tag, y = prop, fill = rtag), position = position_fill(reverse = T), stat = 'identity', width = .8) +
    scale_x_discrete(breaks = tps$tag, labels = tps$lab, expand = c(0,0)) +
    scale_y_continuous(name = sprintf("Proportion Genes"), breaks = c(.25,.5,.75), expand = c(0,0)) +
    scale_fill_manual(values = cols.reg1) +
    coord_flip() +
    theme_bw() +
    theme(legend.position = c(1,1),  legend.justification=c(1,0), legend.background = element_blank()) +
    guides(direction = 'horizontal', fill = guide_legend(nrow = 2, byrow = T)) +
    theme(legend.title = element_blank(), legend.key.size = unit(.8, 'lines'), legend.text = element_text(size = 8)) +
    theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank()) +
    theme(plot.margin = unit(c(2.5,.5,.5,.5), "lines")) +
    theme(axis.title.x = element_text(size = 9)) +
    theme(axis.title.y = element_blank()) +
    theme(axis.text = element_text(size = 8))
#}}}

#{{{ p8: consis. Reg2 <> consis. Reg1 (density)
gids = tsh_r %>% filter(ctag == 'Reg1', n.tis.tot >= 5) %>%
    distinct(gid) %>% pull(gid)
tm2 = tm %>% filter(gid %in% gids, !is.na(Reg1), Tissue %in% tissues20) %>%
    filter(!is.na(bic.c)) %>%
    group_by(gid) %>%
    summarise(sd.bic = sd(bic.c-bic.t))
tp = txb %>% filter(b+m >= 3) %>%
    transmute(gid = gid, tag = tag) %>%
    inner_join(tm2, by = 'gid')
p8 = ggplot(tp) +
    geom_density(aes(x = sd.bic, color = tag)) +
    scale_x_continuous(name = 'sd(loglik(cis) - loglik(trans))', limits = c(0,15), expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    scale_color_aaas() +
    theme_bw() +
    theme(legend.position = c(1,1),  legend.justification=c(1,1), legend.background = element_blank()) +
    guides(direction = 'vertital', color = guide_legend(ncol = 1, byrow = F)) +
    theme(legend.title = element_blank(), legend.key.size = unit(.7,'lines'), legend.text = element_text(size = 8)) +
    theme(legend.margin = margin(t=.1, r=.1, unit='lines')) +
    theme(plot.margin = unit(c(.5,.5,.5,.5), "lines")) +
    theme(axis.title.x = element_text(size = 9)) +
    theme(axis.title.y = element_blank()) +
    theme(axis.text = element_text(size = 8))
#}}}

#{{{ p4: #B/M-biased genes broken by #non-DE tissues
tp = tsb %>% filter(tag != 'conserved', n.tis >= 1) %>% group_by(gid, n.tis.tot) %>%
    summarise(n.tis = sum(n.tis)) %>%
    ungroup() %>%
    count(n.tis, n.tis.tot)
tps = tp %>% group_by(n.tis) %>% summarise(ng.tot = sum(n)) %>%
    mutate(lab = sprintf("%2d (%4d)", n.tis, ng.tot))
p4 = ggplot(tp) +
    geom_bar(aes(x = n.tis, y = n, fill = n.tis.tot), position = 'stack', stat='identity', width=0.7) +
    scale_x_continuous(name = 'Number B73/M17-biased Tissues', c(5,10,15), expand = c(0,0)) +
    scale_y_continuous(name = "Number Genes", expand = c(0, 0)) +
    scale_fill_viridis(name = 'Number non-DE Tissues', breaks = c(1,2,5,10,15,20), direction = -1) +
#    facet_wrap(~fac, scales='free') +
    theme_bw() +
    theme(legend.position = c(1,1), legend.justification=c(1,1), legend.background = element_blank()) +
    guides(direction = 'vertical', fill = guide_legend(ncol = 1, byrow = F)) +
    theme(legend.title = element_text(size=8), legend.key.size = unit(.8, 'lines'), legend.text = element_text(size = 8)) +
    theme(panel.border = element_blank(), panel.grid = element_blank()) +
    theme(strip.background = element_blank(), strip.text.x = element_blank()) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
    theme(plot.margin = unit(c(.5,.5,.5,.5), "lines")) +
    theme(axis.title = element_text(size = 9)) +
    theme(axis.text = element_text(size = 8))
#}}}

#{{{ p5 a-b:tissue specificity
tp = tsb %>% mutate(tag = factor(tag, levels = rev(regs2)))
tp %>% group_by(tag) %>% summarise(n.tis.median = mean(prop.tis))
tps = tp %>% count(tag) %>% 
    mutate(lab = sprintf("%s (%d)", tag, n))
p5a = ggplot(tp) +
    geom_boxplot(aes(x = tag, y = prop.tis, fill = tag), notch = T, outlier.size = 0.5, width = .8) +
    scale_x_discrete(breaks = tps$tag, labels = tps$lab, expand = c(0,.5)) +
    scale_y_continuous(name = 'Proportion (non-DE) Tissues', limits = c(0,1), breaks = c(.25,.5,.75), expand = c(0,0)) +
    scale_fill_manual(values = cols.reg2) +
    coord_flip() +
    theme_bw() +
    theme(legend.position = 'none') +
    theme(plot.margin = unit(c(2.5,.5,.5,1.5), "lines")) +
    theme(axis.title.x = element_text(size = 9)) +
    theme(axis.title.y = element_blank()) +
    theme(axis.text.x = element_text(size = 8)) +
    theme(axis.text.y = element_text(size = 8))
tp = tsh_r %>% filter(ctag %in% c("Reg2"))  %>%
    count(tag, tsTag) %>%
    mutate(tag = factor(tag, levels = rev(regs2)))
tps = tp %>% group_by(tag) %>% summarise(n = sum(n)) %>%
    mutate(lab = sprintf("%s (%d)", tag, n))
p5b = ggplot(tp) +
    geom_bar(aes(x = tag, y = n, fill = tsTag), position = position_fill(reverse = T), stat = 'identity', width = .8) +
    scale_x_discrete(breaks = tps$tag, labels = tps$lab, expand = c(0,0)) +
    scale_y_continuous(name = sprintf("Proportion Genes"), breaks = seq(.25,.75,by=0.25), expand = c(0,0)) +
    coord_flip() +
    scale_fill_manual(values = cols.ts) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    theme(legend.position = c(1,1),  legend.justification=c(1,0), legend.background = element_blank()) +
    guides(direction = 'horizontal', fill = guide_legend(nrow = 1, byrow = F)) +
    theme(legend.title = element_blank(), legend.key.size = unit(0.8, 'lines'), legend.text = element_text(size = 8)) +
    theme(panel.grid = element_blank()) +
    theme(plot.margin = unit(c(2.5,.5,.5,.5), "lines")) +
    theme(axis.title.x = element_text(size = 9)) +
    theme(axis.title.y = element_blank()) +
    theme(axis.text = element_text(size = 8))
#}}}

#{{{ plot
fo = file.path(dirw, "25.Reg2.pdf")
#ggarrange(p1, p2, p3, p4, p5a, p5b, p7, p8,
#    nrow = 4, ncol = 2, labels=LETTERS[1:8], heights = c(4,4,2,2)) %>%
ggarrange(p1, p2,
    nrow = 2, ncol = 1, labels=LETTERS[1:2], heights = c(4,1.5)) %>%
    ggexport(filename = fo, width = 6, height = 7)
#}}}
#}}}

#{{{ mix-pDE <> mix-Reg1; mix-pDE <> mix-Reg2
#{{{
# pDE
ntissue.pde = max(tsh_d$n.tis.tot)
tags.mix.pde = c("consis. non-DE", "consis. DE_B", "consis. DE_M", "mix of DE_B/DE_M")
tx.pde = tsh_d %>% filter(ctag %in% c("pDE")) %>%
    transmute(gid = gid, tag = tag, n.tis.expr = n.tis.tot, n.tis.pde = n.tis) %>%
    mutate(tag = ifelse(tag == 'non-pDE', 'consis. non-DE',
                 ifelse(tag == 'DE mix', 'mix of DE_B/DE_M', tag))) %>%
    mutate(tag = factor(tag, levels = tags.mix.pde))
# Dom
ntissue.dom = max(tsh_r$n.tis.tot)
tags.mix.dom = c("consis. MP", "consis. LP/BLP", "consis. HP/AHP", "mix of HP/LP")
tx.dom = tsh_r %>% filter(ctag %in% c("Dom")) %>%
    select(-ctag, -prop.tis) %>%
    mutate(tag = as.character(tag)) %>%
    mutate(tag = ifelse(tag %in% c("LP", "BLP"), "LP", tag)) %>%
    mutate(tag = ifelse(tag %in% c("HP", "AHP"), "HP", tag)) %>%
    group_by(gid, n.tis.tot, tag) %>%
    summarise(n.tis = sum(n.tis)) %>% ungroup() %>%
    spread(tag, n.tis) %>%
    replace_na(list(LP=0,HP=0,MP=0)) %>% 
    mutate(tag = ifelse(LP > 0 & HP > 0, 'mix of HP/LP',
                 ifelse(LP > 0, 'consis. LP/BLP', 
                 ifelse(HP > 0, 'consis. HP/AHP', 'consis. MP')))) %>%
    mutate(tag = factor(tag, levels = tags.mix.dom)) %>%
    transmute(gid = gid, tag = tag, n.tis.de = n.tis.tot)
# Reg2
tags.reg2 = c('conserved',"B73 biased","Mo17 biased")
tags.mix.reg2 = c("consis. conserved", "consis. B73-biased", "consis. Mo17-biased", "mix of B73-biased/Mo17-biased")
ntissue.reg2 = max(tsh_r$n.tis.tot)
tx.reg2 = tsh_r %>% filter(ctag == 'Reg2') %>%
    select(-ctag, -prop.tis) %>%
    mutate(tag = as.character(tag)) %>%
    mutate(tag = ifelse(tag == 'B73 biased', 'b', tag)) %>%
    mutate(tag = ifelse(tag == 'Mo17 biased', 'm', tag)) %>%
    group_by(gid, n.tis.tot, tag) %>%
    summarise(n.tis = sum(n.tis)) %>% ungroup() %>%
    spread(tag, n.tis) %>%
    replace_na(list(b=0,m=0,conserved=0)) %>% 
    mutate(tag = ifelse(b > 0 & m > 0, 'mix of B73-biased/Mo17-biased',
                 ifelse(b > 0, 'consis. B73-biased',
                 ifelse(m > 0, 'consis. Mo17-biased', 'consis. conserved')))) %>%
    mutate(tag = factor(tag, levels = tags.mix.reg2)) %>%
    transmute(gid = gid, tag = tag, n.tis.nde = n.tis.tot)
#}}}

#{{{ p1: cis/trans composition of consis. DE
gids.p = tx.pde %>% filter(n.tis.pde >= .25 * ntissue.pde) %>%
    distinct(gid) %>% pull(gid)
gids.r = tsh_r %>% filter(ctag == 'Reg1', n.tis.tot >= 3) %>%
    distinct(gid) %>% pull(gid)
gids = gids.r[gids.r %in% gids.p]
tr1 = tm %>% filter(gid %in% gids, !is.na(Reg1), Tissue %in% tissues20) %>%
    transmute(gid = gid, rtag = Reg1) %>%
    group_by(gid, rtag) %>%
    summarise(n = n()) %>%
    mutate(nc = n/sum(n))
sum(tr1$nc)
tp = tx.pde %>% filter(gid %in% gids)
tps = tp %>% distinct(tag, gid) %>% count(tag) %>%
    mutate(lab = sprintf("%s (%d)", tag, n))
tp = tp %>%
    transmute(gid = gid, tag = tag) %>%
    inner_join(tr1, by = 'gid') %>%
    group_by(tag, rtag) %>%
    summarise(n = sum(nc)) %>%
    mutate(prop = n/sum(n)) %>%
    print(n=20)
sum(tp$n)
p1 = ggplot(tp) +
    geom_bar(aes(x = tag, y = prop, fill = rtag), position = position_fill(reverse = T), stat = 'identity', width = .8, alpha = .8) +
    scale_x_discrete(breaks = tps$tag, labels = tps$lab, expand = c(0,0)) +
    scale_y_continuous(name = sprintf("Proportion Genes"), breaks = c(.25,.5,.75), expand = c(0,0)) +
    scale_fill_manual(values = cols.reg1) +
    coord_flip() +
    theme_bw() +
    theme(legend.position = c(1,1),  legend.justification=c(1,0), legend.background = element_blank()) +
    guides(direction = 'horizontal', fill = guide_legend(nrow = 2, byrow = T)) +
    theme(legend.title = element_blank(), legend.key.size = unit(.8, 'lines'), legend.text = element_text(size = 8)) +
    theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank()) +
    theme(plot.margin = unit(c(2.5,.5,.5,.5), "lines")) +
    theme(axis.title.x = element_text(size = 9)) +
    theme(axis.title.y = element_blank()) +
    theme(axis.text = element_text(size = 8))
p_con1 = p1
#}}}

#{{{ p2 : consis. DE <> consis. cis/trans 
tm2 = tm %>% filter(Tissue %in% tissues20, !is.na(Reg1), !is.na(bic.c), !is.na(bic.t)) %>%
    select(Tissue, gid, bic.c, bic.t)
tp = tx.pde %>% filter(n.tis.pde >= .25 * ntissue.pde) %>%
    inner_join(tm2, by = 'gid') %>%
    group_by(gid) %>%
    summarise(tag = tag[1], sd.bic = sd(bic.c - bic.t)) 
tp %>% group_by(tag) %>% summarise(
    q25 = quantile(sd.bic, .25, na.rm = T),
    q50 = quantile(sd.bic, .50, na.rm = T),
    q75 = quantile(sd.bic, .75, na.rm = T))
p2a = ggplot(tp) +
    geom_density(aes(x = sd.bic, color = tag)) +
    scale_x_continuous(name = 'sd(loglik(cis) - loglik(trans))', limits = c(0,15), expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    scale_color_aaas() +
    theme_bw() +
    theme(legend.position = c(1,1),  legend.justification=c(1,1), legend.background = element_blank()) +
    guides(direction = 'vertital', color = guide_legend(ncol = 1, byrow = F)) +
    theme(legend.title = element_blank(), legend.key.size = unit(0.8, 'lines'), legend.text = element_text(size = 8)) +
    theme(plot.margin = unit(c(.5,.5,.5,.5), "lines")) +
    theme(axis.title.x = element_text(size = 9)) +
    theme(axis.title.y = element_blank()) +
    theme(axis.text.x = element_text(size = 8)) +
    theme(axis.text.y = element_text(size = 8)) 
p2b = ggplot(tp) +
    geom_boxplot(aes(x =tag, y = sd.bic, fill = tag), notch = T, outlier.size = .5, width = .8) +
    scale_x_discrete(expand = c(0,.5)) +
    scale_y_continuous(name = 'sd(loglik(cis) - loglik(trans))', expand = c(0,0), limits = c(0,15)) +
    coord_flip() +
    scale_fill_aaas() +
    theme_bw() +
    theme(legend.position = 'none') +
    theme(plot.margin = unit(c(.5,.5,.5,.5), "lines")) +
    theme(axis.title.x = element_text(size = 9)) +
    theme(axis.title.y = element_blank()) +
    theme(axis.text = element_text(size = 8))
#ggsave(p, filename = file.path(dirw, 'test.pdf'), width = 6, height = 3)
#}}}

#{{{ p3: pDE <> Reg2
tx1 = tx.pde %>% filter(n.tis.pde >= 2) %>%
    transmute(gid = gid, tag.pde = tag)
tx2 = tx.reg2 %>% filter(n.tis.nde >= 2) %>%
    transmute(gid = gid, tag.reg2 = tag)
tp = tx1 %>% inner_join(tx2, by = 'gid') %>%
    group_by(tag.pde, tag.reg2) %>%
    summarise(n = n()) %>%
    mutate(prop = n/sum(n))
tps = tp %>% group_by(tag.pde) %>%
    summarise(ntot = sum(n)) %>% 
    mutate(lab = sprintf("%s (%d)", tag.pde, ntot))
p3 = ggplot(tp) +
    geom_bar(aes(x = tag.pde, y = n, fill = tag.reg2), position = position_fill(reverse = T), stat = 'identity', width = 0.8) +
    scale_x_discrete(breaks = tps$tag.pde, labels = tps$lab, expand = c(0,0)) +
    scale_y_continuous(name = sprintf("Prop. Genes"), breaks = c(.25,.5,.75), expand = c(0, 0)) +
    coord_flip() +
    scale_fill_manual(values = cols.mix4) +
    theme_bw() +
    theme(legend.position = c(1,1),  legend.justification=c(1,0)) +
    guides(direction = 'horizontal', fill = guide_legend(nrow = 2, byrow = F)) +
    theme(legend.title = element_blank(), legend.key.size = unit(0.8, 'lines'), legend.text = element_text(size = 8)) +
    theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank()) +
    theme(plot.margin = unit(c(5.5,.5,.5,.5), "lines")) +
    theme(axis.title.x = element_text(size = 9)) +
    theme(axis.title.y = element_blank()) +
    theme(axis.text = element_text(size = 8)) 
#}}}

#{{{ plot
fo = file.path(dirw, "26.cons1.pdf")
ggarrange(
    ggarrange(p1, p2b, nrow = 1, ncol = 2, labels = c('A', 'C')),
    ggarrange(p2a, p3, nrow = 1, ncol = 2, labels = c('B', 'D')),
    nrow = 2, ncol = 1, heights = c(2,3)) %>% 
    ggexport(filename = fo, width = 8, height = 4)
#}}}
#}}}

#{{{ create supplemental table for GO
fi = file.path(dirw, '51.go.enrich', '01.exp.tsv')
t01 = read_tsv(fi) %>%
    filter(tag %in% c('Constitutive', 'Silent')) %>%
    mutate(ctag = 'All genes') %>%
    filter(source == 'uniprot.plants', gotype == c("P", "C"))
fi = file.path(dirw, '51.go.enrich', '10.pDE.tsv')
t10 = read_tsv(fi)
fi = file.path(dirw, '51.go.enrich', '11.pDE.mix.tsv')
t11 = read_tsv(fi)
fi = file.path(dirw, '51.go.enrich', '16.hDE.tsv')
t16 = read_tsv(fi) %>%
    filter(tag %in% c("AP", "BP")) %>%
    mutate(tag = ifelse(tag == 'AP', 'Above-Parent', tag)) %>%
    mutate(tag = ifelse(tag == 'BP', 'Below-Parent', tag)) %>%
    arrange(tag, pval.adj) %>%
    mutate(ctag = 'non-DEGs btw. B73 and Mo17') %>%
    filter(source == 'uniprot.plants', gotype %in% c("P", "C"))
fi = file.path(dirw, '51.go.enrich', '21.dom.tsv')
t21 = read_tsv(fi)
fi = file.path(dirw, '51.go.enrich', '22.dom.mix.tsv')
t22 = read_tsv(fi) %>%
    filter(tag %in% c("consis. HP/AHP", "consis. LP/BLP")) %>%
    mutate(tag = ifelse(tag == 'consis. HP/AHP', 'HP/AHP', tag)) %>%
    mutate(tag = ifelse(tag == 'consis. LP/BLP', 'LP/BLP', tag)) %>%
    mutate(ctag = 'DEGs btw. B73 and Mo17') %>%
    filter(source == 'uniprot.plants', gotype %in% c('P','C'))

groupnames = c('All genes' = nrow(t01),
               'non-DEGs btw. B73 and Mo17' = nrow(t16),
               'DEGs btw. B73 and Mo17' = nrow(t22))
tx = rbind(t01, t16, t22) %>%
    mutate(pval.adj = sprintf("%.02e",pval.adj)) %>%
    transmute(' ' = ctag, 'Gene set' = tag, 
              'GO' = goid, 'P-value' = pval.adj, 'GO name' = goname)
fx = file.path(dirw, '51.go.enrich.rda')
save(tx, groupnames, file = fx)
#}}}

#{{{# tissue specificity btw genotypes - tau
fg = '/home/springer/zhoux379/data/genome/B73/GO/maize.B73.AGPv4.aggregate.gaf'
tg = read_tsv(fg, skip = 1)
gids = unique(tg$db_object_id[tg$term_accession == 'GO:0051707'])

tts = tm %>% filter(silent==0) %>% count(gid)
tt = tm %>% filter(gid %in% tts$gid[tts$n >= 5]) %>%
#    filter(gid %in% gids) %>%
    group_by(gid) %>%
    summarise(B73 = sum(1-B73/max(B73))/(length(B73)-1),
              Mo17 = sum(1-Mo17/max(Mo17))/(length(Mo17)-1),
              BxM = sum(1-BxM/max(BxM))/(length(BxM)-1)) %>%
    ungroup() %>%
    gather(gt, tau, -gid)

p1 = ggplot(tt) +
    geom_density(aes(x = tau, color = gt), alpha = 0.9) +
    scale_x_continuous(name = 'tau', ) +
    scale_y_continuous(name = '', ) +
    scale_color_npg() +
    #theme_minimal() +
    theme(legend.position = 'top', legend.direction = "horizontal") +
    theme(legend.title = element_blank(), legend.key.size = unit(0.8, 'lines'), legend.text = element_text(size = 8)) +
    #theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(0.3,0.3,0.3,0.3), "lines")) +
    theme(axis.title.x = element_text(size = 9)) +
    theme(axis.title.y = element_text(size = 9)) +
    theme(axis.text.x = element_text(size = 8, color = "black", angle = 0)) +
    theme(axis.text.y = element_text(size = 8, color = "black", angle = 0))
fo = file.path(dirw, '15.tau.pdf') 
ggsave(p1, filename = fo, width = 4, height = 4)
#}}}

#{{{# QC: ePAV gene expression in F1
tq = tm %>% filter(B_only) %>%
    mutate(log2HB = log2(BxM/B73), BxM = asinh(BxM)) %>%
    select(Tissue, gid, B_and_H, log2HB, BxM)
describe(tq$log2HB)
tq$log2HB[tq$log2HB == '-Inf'] = -8

p1 = ggplot(tq) +
    geom_histogram(aes(x = BxM, fill = B_and_H), col = 'white', position = 'dodge') +
    geom_histogram(aes(x = , fill = B_and_H), col = 'white', position = 'dodge') +
    scale_x_continuous(name = 'log2(BxM)') +
    scale_y_continuous(name = 'Num Genes') +
    scale_fill_brewer(palette = "Set2") +
    theme_bw() +
    theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(.5,.5,.5,.5), "lines")) +
    theme(legend.position = c(.7,.8), legend.direction = "vertical", legend.justification = c(.5,.5), legend.background = element_blank()) +
    theme(legend.title = element_text(size=9), legend.key.size = unit(.5, 'lines'), legend.key.width = unit(.5, 'lines'), legend.text = element_text(size = 8)) +
    theme(axis.title.x = element_text(size = 9)) +
    theme(axis.title.y = element_blank()) +
    theme(axis.text.x = element_text(size = 8)) +
    theme(axis.text.y = element_text(size = 8, hjust = 1))
fo = sprintf("%s/91.qc.compb.pdf", dirw)
ggsave(p1, filename = fo, width = 5, height = 5)
p1 = ggplot(tq) +
    geom_histogram(aes(x = log2HB, fill = B_and_H), col = 'white', position = 'dodge') +
    scale_x_continuous(name = 'log2(BxM / B)') +
    scale_y_continuous(name = 'Num Genes') +
    scale_fill_brewer(palette = "Set2") +
    theme_bw() +
    theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines")) +
    theme(legend.position = c(0.7,0.8), legend.direction = "vertical", legend.justification = c(0.5,0.5)) +
    theme(legend.title = element_text(size=9), legend.key.size = unit(0.5, 'lines'), legend.key.width = unit(0.5, 'lines'), legend.text = element_text(size = 8), legend.background = element_rect(fill=NA, size=0)) +
    theme(axis.title.x = element_text(size = 9)) +
    theme(axis.title.y = element_blank()) +
    theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
    theme(axis.text.y = element_text(size = 8, colour = "black", angle = 0, hjust = 1))
fo = sprintf("%s/91.qc.compb1.pdf", dirw)
ggsave(p1, filename = fo, width = 5, height = 5)

tq = tm %>% filter(M_only) %>%
    mutate(log2HM = log2(BxM/Mo17), BxM = asinh(BxM)) %>%
    select(Tissue, gid, M_and_H, log2HM, BxM)
describe(tq$log2HM)
tq$log2HM[tq$log2HM == '-Inf'] = -8

p1 = ggplot(tq) +
    geom_histogram(aes(x = BxM, fill = M_and_H), col = 'white', position = 'dodge') + 
    scale_x_continuous(name = 'log2(BxM)') +
    scale_y_continuous(name = 'Num Genes') +
    scale_fill_brewer(palette = "Set2") +
    theme_bw() +
    theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines")) +
    theme(legend.position = c(0.7,0.8), legend.direction = "vertical", legend.justification = c(0.5,0.5)) +
    theme(legend.title = element_text(size=9), legend.key.size = unit(0.5, 'lines'), legend.key.width = unit(0.5, 'lines'), legend.text = element_text(size = 8), legend.background = element_rect(fill=NA, size=0)) +
    theme(axis.title.x = element_text(size = 9)) +
    theme(axis.title.y = element_blank()) +
    theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
    theme(axis.text.y = element_text(size = 8, colour = "black", angle = 0, hjust = 1))
fo = sprintf("%s/91.qc.compm.pdf", dirw)
ggsave(p1, filename = fo, width = 5, height = 5)
p1 = ggplot(tq) +
    geom_histogram(aes(x = log2HM, fill = M_and_H), col = 'white', position = 'dodge') +
    scale_x_continuous(name = 'log2(BxM / M)') +
    scale_y_continuous(name = 'Num Genes') +
    scale_fill_brewer(palette = "Set2") +
    theme_bw() +
    theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines")) +
    theme(legend.position = c(0.7,0.8), legend.direction = "vertical", legend.justification = c(0.5,0.5)) +
    theme(legend.title = element_text(size=9), legend.key.size = unit(0.5, 'lines'), legend.key.width = unit(0.5, 'lines'), legend.text = element_text(size = 8), legend.background = element_rect(fill=NA, size=0)) +
    theme(axis.title.x = element_text(size = 9)) +
    theme(axis.title.y = element_blank()) +
    theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
    theme(axis.text.y = element_text(size = 8, colour = "black", angle = 0, hjust = 1))
fo = sprintf("%s/91.qc.compm1.pdf", dirw)
ggsave(p1, filename = fo, width = 5, height = 5)
#}}}

#{{{ DE heatmap
tissue1 = tissues23[1]
fo = sprintf("%s/24.DE.heatmap1.pdf", dirw)
gids_e = tm %>% filter(silent == 0) %>%
    count(gid) %>%
    filter(n >= 20) %>% pull(gid)
td = tm %>% filter(Tissue == tissue1) %>%
    filter(pDE != 'non_DE', abs(log2mb) >= 1) %>%
    filter(gid %in% gids_e)

th = tm %>% filter(gid %in% td$gid) %>%
    select(Tissue, gid, log2mb) %>%
    mutate(log2mb = pmin(log2mb, 2)) %>%
    mutate(log2mb = pmax(log2mb, -2)) %>%
    spread(Tissue, log2mb)
hdist = daisy(th[,-1], metric = 'gower')
hdist[is.na(hdist)] = 0
hcl = hclust(hdist, method = "ward.D")
tsh_d0 = tsh_d %>% filter(ctag == 'pDE') %>%
    select(gid, tsTag, mixTag = tag)
ths = th %>% select(gid) %>%
    mutate(clu = hcl$order) %>%
    inner_join(tsh_d0, by = 'gid') %>%
    arrange(tsTag, mixTag, clu) %>%
    mutate(y = 1:length(gid))
ths1 = ths %>% group_by(tsTag) %>%
    summarise(ymin = min(y), ymax = max(y), y = mean(y)) %>%
    ungroup()
ths2 = ths %>% group_by(tsTag, mixTag) %>%
    summarise(ymin = min(y), ymax = max(y), y = mean(y)) %>%
    ungroup()
tt = tibble(tissue = tissues23, x = 1:length(tissues23))

hcl = hclust(daisy(t(th[,-1])), method = "ward.D")
tl = th %>%
    inner_join(ths[,c('gid','y')], by = 'gid') %>% select(-gid) %>%
    gather(tissue, log2mb, -y) %>%
    inner_join(tt, by = 'tissue') %>% select(-tissue)
#
p1 = ggplot(tl) +
    geom_tile(aes(x = x, y = y, fill = log2mb)) +
    geom_segment(aes(x = 0, xend = 0, y = 0, yend = max(tl$y)), size = .1) +
    geom_segment(aes(x = -5, xend = -5, y = 0, yend = max(tl$y)), size = .1) +
    geom_segment(data = ths2, aes(x = -.2, xend = 0, y = ymin, yend = ymin), size = .2) +
    geom_segment(data = ths2, aes(x = -.2, xend = 0, y = ymax, yend = ymax), size = .2) +
    geom_text(data = ths2, aes(x = -.1, y = y, label = mixTag), hjust = 1, size = 2.5) +
    geom_segment(data = ths1, aes(x = -5.2, xend = -5, y = ymin, yend = ymin), size = .2) +
    geom_segment(data = ths1, aes(x = -5.2, xend = -5, y = ymax, yend = ymax), size = .2) +
    geom_text(data = ths1, aes(x = -5.1, y = y, label = tsTag), hjust = 1, size = 2.5) +
    scale_x_continuous(breaks = tt$x, labels = tt$tissue, limits = c(-10, 23), expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    scale_fill_gradient2(breaks = c(-1.5, 1.5), labels = c("B > M", "B < M")) + 
    otheme() +
    theme(panel.border = element_blank()) +
    theme(legend.position = 'right', legend.direction = "vertical", legend.justification = c(.5,.5)) +
    theme(axis.text.x = element_text(size = 8, angle = -45, hjust = 0))
ggsave(p1, filename = fo, width = 8, height = 12)
#}}}

#{{{# hDE heatmap
tf1 = tg %>% filter(comp >= 1)
th = tm %>% filter(Tissue %in% tissues2, gid %in% tf1$gid) %>%
    select(Tissue, gid, log2HM) %>%
    spread(Tissue, log2HM)

hdist = daisy(th[,-1], metric = 'gower')
hdist[is.na(hdist)] = 0
hcl = hclust(hdist, method = "ward.D")
gidsO = th$gid[hcl$order]

hcl = hclust(daisy(t(th[,-1])), method = "ward.D")
tissuesO = tissues2[hcl$order] 

tl = th %>%
    gather(tissue, log2HM, -gid) %>%
    mutate(gid = factor(gid, levels = gidsO),
           tissue = factor(tissue, levels = tissuesO))

p1 = ggplot(tl) +
    geom_tile(aes(x = tissue, y = gid, fill = log2HM)) + 
    #scale_x_discrete(name = '') +
    #scale_y_discrete(name = 'Genes') +
    scale_fill_gradient2(breaks = c(-4, 4), labels = c("F1 < MP", "F1 > MP")) + 
    theme_bw() +
    theme(panel.grid = element_blank(), panel.border = element_rect(fill=NA, linetype=0)) +
    theme(plot.margin = unit(c(0.3,0.3,0.3,0.3), "lines")) +
    theme(legend.position = 'right', legend.direction = "vertical", legend.justification = c(0.5,0.5)) + 
    theme(legend.title = element_blank(), legend.key.size = unit(0.5, 'lines'), legend.key.width = unit(0.5, 'lines'), legend.text = element_text(size = 7)) +
    theme(axis.ticks.length = unit(0, 'lines')) +
    theme(axis.title.x = element_blank()) +
    theme(axis.title.y = element_blank()) +
    theme(axis.text.x = element_text(size = 8, colour = "black", angle = -45, hjust = 0)) +
    theme(axis.text.y = element_blank())
fo = sprintf("%s/29.comp.heatmap.pdf", dirw)
ggsave(p1, filename = fo, width = 6, height = 12)
#}}}

#{{{# FPKM based tissue similarity
grp = group_by(ti, gid)
tog = summarise(grp, ntis_b = sum(B73>=1), ntis_m = sum(Mo17>=1))
gids = tog$gid[tog$ntis_b > 0 & tog$ntis_m > 0]

to = ti
to2 = to[to$gid %in% gids, c("Tissue","gid","B73","Mo17")]
to3 = gather(to2, gt, FPKM, -Tissue, -gid)
to3 = within(to3, {
	tis_gt = sprintf("%s|%s", Tissue, gt)
	rm(Tissue, gt)
})
to4 = spread(to3, tis_gt, FPKM)

e = asinh(to4[,-1])
pca <- prcomp(e, center = F, scale. = F)
x = pca['rotation'][[1]]
y = summary(pca)$importance
y[,1:5]
xlab = sprintf("PC1 (%.01f%%)", y[2,1]*100)
ylab = sprintf("PC2 (%.01f%%)", y[2,2]*100)
tp = cbind.data.frame(tis_gt = rownames(x), x[,1:5], stringsAsFactors = F)
res = strsplit(tp$tis_gt, split = "[|]")
tp = cbind(tp, Tissue = sapply(res, "[", 1), Genotype = sapply(res, "[", 2))
tp$Tissue = factor(tp$Tissue, levels = unique(tl$Tissue))
tp$Genotype = factor(tp$Genotype, levels = unique(tp$Genotype))

cols = c(brewer.pal(8, 'Dark2'), brewer.pal(9, 'Set1'))

p1 = ggplot(tp) +
  geom_point(aes(x = PC1, y = PC2, shape = Genotype, color = Tissue), size = 4) +
  scale_x_continuous(name = xlab) +
  scale_y_continuous(name = ylab) +
  scale_color_manual(name = "", values = cols) +
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(plot.margin = unit(c(0.3,0.1,0.1,0.1), "lines")) +
  theme(legend.position = 'right', legend.direction = "vertical", legend.justification = c(0,0.5), legend.title = element_blank(), legend.key.size = unit(1, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size = 8), legend.background = element_rect(fill=NA, size=0)) +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.title.y = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "black", angle = 0, hjust = 1))
fp = sprintf("%s/10.pca.FPKM.pdf", dirw)
ggsave(p1, filename = fp, width = 6, height = 5.5)
#}}}

#{{{# DoA similarity
ti2 = ti[!ti$Tissue %in% c("endosperm_14DAP", "endosperm_27DAP", "kernel_14DAP"),]

grp = group_by(ti2, gid)
tig = summarise(grp, 
	ntis_de = sum(is.de != 'non-DE')
)
gids = tig$gid[tig$ntis_de >= 10]

ti2 = ti2[ti2$gid %in% gids, c("Tissue","gid","is.de","DoA")]
ti2$DoA[ti2$is.de == 'non-DE'] = NA

describe(ti2$DoA)
ti2$DoA[!is.na(ti2$DoA) & ti2$DoA < -3] = -3
ti2$DoA[!is.na(ti2$DoA) & ti2$DoA > 3] = 3

ti3 = spread(ti2[,-3], Tissue, DoA)

e = ti3[,-1]
e[is.na(e)] = 0

pca <- prcomp(e, center = F, scale. = F)
x = pca['rotation'][[1]]
y = summary(pca)$importance
y[,1:5]
xlab = sprintf("PC1 (%.01f%%)", y[2,1]*100)
ylab = sprintf("PC2 (%.01f%%)", y[2,2]*100)
tp = cbind.data.frame(Tissue = rownames(x), x[,1:5], stringsAsFactors = F)
tp$Tissue = factor(tp$Tissue, levels = unique(tl$Tissue))

cols = c(brewer.pal(8, 'Dark2'), brewer.pal(9, 'Set1'))

p1 = ggplot(tp) +
  geom_point(aes(x = PC1, y = PC2, color = Tissue), size = 4) +
  geom_text(aes(x = PC1, y = PC2, label = Tissue), size = 3, nudge_y = 0.03, check_overlap = T) +
  scale_x_continuous(name = xlab) +
  scale_y_continuous(name = ylab) +
  scale_color_manual(name = "", values = cols) +
  theme_bw() +
  theme(axis.ticks.length = unit(0, 'lines')) +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines")) +
  theme(legend.position = 'none') +
  theme(axis.title.x = element_text(size = 9)) +
  theme(axis.title.y = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
  theme(axis.text.y = element_text(size = 8, colour = "black", angle = 0, hjust = 1))
fp = sprintf("%s/12.pca.DoA.pdf", dirw)
ggsave(p1, filename = fp, width = 5.5, height = 5.5)
#}}}

#{{{# snp density v.s. DE 
fv = '~/data/genome/Mo17/62.gene.vnt.tsv'
tv = read_tsv(fv) %>%
    mutate(vtag = ifelse(nvnt > 0, ifelse(nvnt >= 5, 'vnt', 'svnt'), 'ident'))

tp = tm %>% mutate(ase = ifelse(is.na(Reg1) & is.na(Reg2), 'unk', 'kno')) %>%
    select(Tissue, gid, expr, log2MB, pDE, pDE2, ase) %>%
    count(Tissue, pDE, ase) %>% 
    group_by(Tissue, pDE) %>%
    summarise(prop.unk = n[ase=='unk']/sum(n)) %>%
    replace_na(list(pDE = 'non-expr')) %>%
    spread(pDE, prop.unk) %>%
    print(n = 40)

vmap = c("ident"=0, 'svnt'=1, 'vnt'=2)
tp = tm %>% 
    filter(expr == 'expressed') %>%
    select(Tissue, gid, expr, log2MB, pDE, pDE2) %>%
    mutate(tag = ifelse(is.na(pDE), expr, pDE)) %>%
    inner_join(tv, by = 'gid') %>%
    mutate(vtag2 = vmap[vtag]) %>%
    group_by(Tissue) %>%
    do(tidy(lm(abs(log2MB) ~ vtag2, .)))

tp %>%
    group_by(Tissue, vtag) %>%
    summarise(q25 = quantile(abs(log2MB), .25),
              q50 = quantile(abs(log2MB), .5),
              q75 = quantile(abs(log2MB), .75)) %>%
    print(n=40)
p = ggplot(tp) +
    geom_boxplot(aes(x = vtag, y = abs(log2MB)), notch = T, outlier.shape = NA, width = .8) +
    scale_x_discrete(expand = c(0,.5), labels = c("no variant", "1-5 vaariants", ">5 variant")) +
    scale_y_continuous(name = 'abs(log2(B73/Mo17))', limits = c(0,1.5), breaks = c(.5,1), expand = c(0,0)) +
    coord_flip() +
    facet_wrap(~Tissue, nrow = 5) +
    theme_bw() +
    theme(legend.position = 'none') +
    theme(plot.margin = unit(c(.5,.5,.5,.5), "lines")) +
    theme(axis.title.x = element_text(size = 9)) +
    theme(axis.title.y = element_blank()) +
    theme(axis.text.x = element_text(size = 8)) +
    theme(axis.text.y = element_text(size = 8)) 
ggsave(p, filename=file.path(dirw, 'test.snpden.pdf'), width = 8, height = 6)

tp = tv %>% inner_join(tsh_e, by = 'gid') %>%
    group_by(vtag, etag) %>%
    summarise(n = n()) %>%
    mutate(prop = n/sum(n))

tsh_pde = tsh_d %>% filter(n.tis.tot >= 5, ctag == 'pDE')
tp = tv %>% inner_join(tsh_pde, by = 'gid') %>%
    group_by(vtag, tag) %>%
    summarise(n = n()) %>%
    mutate(prop = n/sum(n))

tm %>% filter(B73 < .1, Mo17 < .1, BxM > 1) %>%
    select(Tissue, gid, B73, Mo17, BxM, hDE, prop.h) %>% print(n=50)
#}}}


#{{{ test ASE with gene sequence changes
tt = tm %>% filter(pDE == 'non_DE', !is.na(Reg2)) %>%
    select(Tissue, gid, hDE, prop.h, Reg2)
fv = '~/projects/wgc/data/05_stats/10.B73_Mo17.tsv'
impacts = c("no_change","low","modifier",'moderate','high')
tv = read_tsv(fv) %>% filter(impact %in% c("low",'modifier','moderate')) %>% select(gid, impact, eff)

tt2 = tt %>% inner_join(tv, by='gid') %>%
    #mutate(Reg2=hDE) %>%
    mutate(impact=eff)
tt3 = tt2 %>% count(Tissue, impact) %>% rename(n_tot = n) %>% filter(n_tot >= 100)
tt4 = tt2 %>% count(Tissue, impact, Reg2) %>%
    inner_join(tt3, by=c("Tissue","impact")) %>%
    mutate(prop = n/n_tot) %>%
    select(Tissue,impact,Reg2,prop) %>%
    #mutate(impact = factor(impact, levels=impacts)) %>%
    spread(Reg2, prop)
tt4 %>% print(n=100)
#}}}

