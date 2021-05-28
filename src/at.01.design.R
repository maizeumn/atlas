source("functions.R")
require(grImport)
require(magick)
require(ggpubr)
load_fonts()
#
f1 = file.path(dira, 'images/lowres-hiseq2500-right.jpg')
f2 = file.path(dira, 'images/hplc.jpg')
f3 = file.path(dira, 'images/LTQ.jpg')
fac0='\ue901'; fac1='\ue902'; fac2='\ue903'
fao0='\ue904'; fao1='\ue905'; fao2='\ue906'
arrow1 = arrow(length = unit(0.02, "npc"), type="closed")

#{{{ atlas design
tp = crossing(tissue=tissues23, genotype=gts3) %>%
    mutate(tissue=factor(tissue, levels=tissues23)) %>%
    mutate(genotype=factor(genotype,levels=gts3))
tpx = tp %>% distinct(genotype) %>% arrange(genotype) %>% mutate(x=1:n())
tpy = tp %>% distinct(tissue) %>% arrange(tissue) %>% mutate(y=1:n())
tp = tp %>% inner_join(tpx, by='genotype') %>%
    inner_join(tpy, by = 'tissue')

colR = pal_aaas()(2)[1]; colP = pal_aaas()(2)[2]
x.off1 = .6; x.off2 = .7; y.off1 = 1; y.off2 = .4
px = ggplot(tp) +
    geom_text(aes(x-.2, y, color='mRNA'), label=fac1, family='icm',size=5)+
    geom_text(aes(x+.2, y, color='protein'), label=fac2, family='icm',size=6)+
    geom_text(aes(x-.2, y, color='mRNA'), label='x3', size=2, hjust=-.6, vjust=0, show.legend=F) +
    geom_text(aes(x+.2, y, color='protein'), label='x3', size=2, hjust=-.6, vjust=0, show.legend=F) +
    scale_y_reverse(breaks=tpy$y-.2, labels=tpy$tissue, expand=c(0,0)) +
    scale_x_continuous(breaks=tpx$x, labels=tpx$genotype, expand=expand_scale(mult=c(0,0)), position="top") +
    scale_color_manual(values=pal_aaas()(5)) +
    expand_limits(x = c(1-x.off1-.05,15)) +
    expand_limits(y = c(1-y.off1-.05,23+y.off2+.05)) +
    annotate('rect',xmin=1-x.off1,xmax=3+x.off2,ymin=1-y.off1,ymax=23+y.off2,color='black',fill=NA, size=.2) +
    annotate('segment',x=4,xend=7,y=7,yend=7,arrow=arrow1,color=colR,size=1) +
    annotate('segment',x=4,xend=7,y=16,yend=16,arrow=arrow1,color=colP,size=1) +
    annotate('text',x=5.5,y=7,color=colR,size=5,vjust=-.5,label='mRNA') +
    annotate('text',x=5.5,y=16,color=colP,size=5,vjust=-.5,label='protein') +
    annotate('text',x=11,y=4,size=3,label=expression("Illumina"^"\U00AE"~HiSeq^"\U2122"*2500), parse=T) +
    annotate('text',x=11,y=12.5,size=3,label=expression('Agilent 1100 HPLC / Thermo Scientific'^"\U2122"~'TMT10plex'), parse=T) +
    #scale_fill_manual(values=c('springgreen','steelblue1','firebrick1')) +
    otheme(legend.pos='none',legend.dir='v',legend.title=F,
           margin = c(.1,.1,.1,.1),
           xtitle=F, xtext=T, ytext=T, xgrid=F, ygrid=F) +
    theme(axis.text.y = element_text(color=cols23)) +
    theme(axis.text.x = element_text(color=pal_futurama()(3))) +
    theme(panel.border = element_blank()) +
    theme(panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA))
p = ggdraw() +
    draw_image(f1, scale=.2, x=.78,y=.68,hjust=.5,vjust=.5) +
    draw_image(f2, scale=.2, x=.65,y=.3,hjust=.5,vjust=.5) +
    draw_image(f3, scale=.2, x=.85,y=.3,hjust=.5,vjust=.5) +
    draw_plot(px)
fo = file.path(dirf, '00.design.atlas.pdf')
ggsave(p, file=fo, width=8, height=6)
#}}}

#{{{ 6-hybrid design
tp = crossing(tissue=tissues3, genotype=gts10) %>%
    mutate(tissue=factor(tissue, levels=tissues3)) %>%
    mutate(genotype=factor(genotype,levels=gts10))
tpx = tp %>% distinct(tissue) %>% arrange(tissue) %>% mutate(x=1:n())
tpy = tp %>% distinct(genotype) %>% arrange(genotype) %>% mutate(y=1:n())
tp = tp %>% inner_join(tpx, by='tissue') %>%
    inner_join(tpy, by = 'genotype')

cols3 = pal_npg()(3)
colR = cols3[1]; colT = cols3[2]; colP = cols3[3]
x.off1 = .5; x.off2 = .5; y.off1 = .7; y.off2 = .3
px = ggplot(tp) +
    geom_text(aes(x-.3, y), color=colR, label=fac1, family='icm',size=5)+
    geom_text(aes(x, y), color=colT, label=fac1, family='icm',size=5)+
    geom_text(aes(x+.3, y), color=colP, label=fac2, family='icm',size=6)+
    geom_text(aes(x-.3, y), color=colR, label='x3', size=2, hjust=-.6, vjust=0, show.legend=F) +
    geom_text(aes(x, y), color=colT, label='x3', size=2, hjust=-.6, vjust=0, show.legend=F) +
    geom_text(aes(x+.3, y), color=colP, label='x3', size=2, hjust=-.6, vjust=0, show.legend=F) +
    geom_segment(data=tpx[-nrow(tpx),], aes(x=x+.55,xend=x+.55,y=1-y.off1,yend=10+y.off2), size=.1) +
    scale_y_reverse(breaks=tpy$y-.1, labels=tpy$genotype, expand=c(0,0)) +
    scale_x_continuous(breaks=tpx$x, labels=tpx$tissue, expand=expand_scale(mult=c(0,0)), position="top") +
    scale_color_manual(values=pal_npg()(5)) +
    expand_limits(x = c(1-x.off1-.05,10)) +
    expand_limits(y = c(1-y.off1-.05,10+y.off2+.05)) +
    annotate('rect',xmin=1-x.off1,xmax=3+x.off2,ymin=1-y.off1,ymax=10+y.off2,color='black',fill=NA, size=.2) +
    annotate('segment',x=4,xend=6,y=3,yend=3,arrow=arrow1,color=colR,size=1) +
    annotate('segment',x=4,xend=6,y=3.2,yend=3.2,arrow=arrow1,color=colT,size=1) +
    annotate('segment',x=4,xend=6,y=8,yend=8,arrow=arrow1,color=colP,size=1) +
    annotate('text',x=5,y=3,color=colR,size=5,vjust=-.5,label='mRNA') +
    annotate('text',x=5,y=3.2,color=colT,size=5,vjust=1.5,label='total RNA') +
    annotate('text',x=5,y=8,color=colP,size=5,vjust=-.5,label='protein') +
    annotate('text',x=7.7,y=2,size=3,label=expression("Illumina"^"\U00AE"~HiSeq^"\U2122"*2500), parse=T) +
    annotate('text',x=7.7,y=6,size=3,label=expression('Agilent 1100 HPLC / Thermo Scientific'^"\U2122"~'TMT10plex'), parse=T) +
    #scale_fill_manual(values=c('springgreen','steelblue1','firebrick1')) +
    otheme(legend.pos='none',legend.dir='v',legend.title=F,
           margin = c(.1,.1,.1,.1),
           xtitle=F, xtext=T, ytext=T, xtick=F, ygrid=F) +
    theme(axis.text.y = element_text(color=pal_aaas()(10))) +
    theme(axis.text.x = element_text(color=pal_futurama()(3))) +
    theme(panel.border = element_blank()) +
    theme(panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA))
p = ggdraw() +
    draw_image(f1, scale=.2, x=.78,y=.65,hjust=.5,vjust=.5) +
    draw_image(f2, scale=.2, x=.68,y=.25,hjust=.5,vjust=.5) +
    draw_image(f3, scale=.2, x=.85,y=.25,hjust=.5,vjust=.5) +
    draw_plot(px)
fo = file.path(dirf, '00.design.hybrid.pdf')
ggsave(p, file=fo, width=8, height=4)
#}}}



#{{{ make B,M,F1 similarity diagram
require(diagram)

gts4 = gts10[1:4]
M = matrix(nrow = 4, ncol = 4, byrow = T, data = c(
    NA, NA, 'B73xMo17', NA,
    'B84xB73', NA, 'B84xMo17', NA,
    'Mo17xB73', NA, NA, NA,
    'A682xB73', NA, 'A682xMo17', NA
))

fo = file.path(dirf, "00.design.cross.pdf")
pdf(fo, width = 5, height = 5)
par(mfrow=c(1,2))
par(mar=c(0.1,0.1,3.1,0.1))
pp1 = plotmat(M, pos = c(2,2), curve = 0, name = gts4,
    lwd = 1, box.lwd = 2, cex.txt = 1, box.cex = 1,
    box.type = "square", box.prop = 0.6, arr.type = "triangle",
    arr.pos = 0.5, arr.width = 0, shadow.size = 0.01, prefix = "",
    main = "Study design")
dev.off()
#}}}

#{{{ match 6-hybrid samplelist to UMGC files
dirw = file.path(dird, '01_exp_design')
fh = file.path(dirw, '61.BR5.meta.tsv')
th = read_tsv(fh)

dir1 = '/home/springer/data_release/umgc/novaseq'
dir2 = c('191010_A00223_0227_AHG3TVDRXX', '191014_A00223_0231_AHG22WDRXX',
    '191023_A00223_0238_BHG57FDRXX', '191023_A00223_0239_AHGLLLDRXX')
dir3 = str_c("Springer_Project_065_Pool", 1:8, sep='')
ti = crossing(dir1=dir1,dir2=dir2,dir3=dir3) %>%
    mutate(directory=file.path(dir1, dir2, dir3)) %>%
    mutate(isfolder = map_lgl(directory, dir.exists)) %>% filter(isfolder) %>%
    select(directory)

to = ti %>%
    mutate(fname = map(directory, list.files, pattern="_R1_")) %>%
    unnest() %>%
    separate(fname, c('SampleID','idx','suf1','suf2'), sep='_') %>%
    mutate(file_prefix=str_c(SampleID, idx, sep='_')) %>%
    mutate(sid = str_sub(SampleID, 2)) %>%
    inner_join(th, by=c('sid'='SampleID')) %>%
    mutate(Treatment = str_sub(SampleID, 1, 1)) %>%
    select(SampleID, Tissue, Genotype, Treatment, Replicate=Rep, directory, file_prefix)

fo = file.path(dirw, '63.BR5.path.tsv')
write_tsv(to, fo)
#}}}

#{{{ RIL design sample list
fi = file.path(dirw, '71.BR7.genotypes.xlsx')
ti = read_xlsx(fi)

to = crossing(Genotype = ti$Genotype, Replicate=1:4) %>%
    mutate(SampleID = sprintf("BR%03d", (700+row_number()))) %>%
    mutate(Tissue = 'leaf_V3') %>%
    select(SampleID, Tissue, Genotype, Replicate)
fo = file.path(dirw, "71.BR7.meta.tsv")
write_tsv(to, fo, na = '')
#}}}


