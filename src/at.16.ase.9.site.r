source("br.fun.r")
dirw = file.path(dirp, "45.coop")

#{{{ read annotation
dirg = '/home/springer/zhoux379/data/genome/Zmays_v4'
f_gtb = file.path(dirg, "v37/t2.gtb")
tg = read.table(f_gtb, sep = "\t", header = T, as.is = T)[,c('id','par','cat1','cat2')]
colnames(tg)[2] = 'gid'
grp = dplyr::group_by(tg, gid)
tg2 = dplyr::summarise(grp, gtype = cat1[1])

fl = file.path(dirw, '00.1.read.correct.tsv')
tl = read.table(fi, header = T, sep = "\t", as.is = T)[,1:5]
tl = within(tl, {label = sprintf("%s: %s %s rep%d", SampleID, Tissue, Genotype, Treatment)})

fv = '/home/springer/zhoux379/data/misc2/mo17vnt/53.vnt.final/65.vnt2gene.bed'
tv = read.table(fv, header = F, sep = "\t", as.is = T)
colnames(tv) = c("chr", "pos", "end", "alleles", "gchr", "gbeg", "gend", "gid", "obp")
tv$pos = tv$pos + 1

ta = data.frame()
sid = 'BR006'
for (i in 1:41) {
    sid = sprintf("BR%03d", i)
    fa = sprintf("%s/25.ase/%s.tsv", dirw, sid)
    ta1 = read.table(fa, sep = "\t", as.is = T, header = T)[,c(1,2,6:8)]
    ta = rbind.data.frame(ta, cbind.data.frame(sid = sid, ta1, stringsAsFactors = F), stringsAsFactors = F)
    cat(sid, "\n")
}

ta2 = ta[ta$depth >= 5,]

to = ta2
to = merge(to, ti[,c(1,6)], by.x = 'sid', by.y = 'SampleID')

p1 = ggplot(to) +
    geom_histogram(aes(x = dpr/depth), bins = 50) +
    scale_x_continuous(name = 'Proportion B73 Allele') +
    scale_y_continuous(name = "# Sites") +
    #scale_color_manual(name = "", values = cols) +
    facet_wrap( ~ label, nrow = 7) +
    theme_bw() +
    theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines")) +
    #theme(legend.position = c(0.4, 0.8), legend.direction = "vertical", legend.justification = c(0,0.5), legend.title = element_blank(), legend.key.size = unit(1, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size = 8), legend.background = element_rect(fill=NA, size=0)) +
    theme(axis.title.x = element_text(size = 9)) +
    theme(axis.title.y = element_blank()) +
    theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0, hjust = 1)) +
    theme(axis.text.y = element_text(size = 8, colour = "black", angle = 0))
fp = file.path(diro, "25.ase.pdf")
ggsave(p1, filename = fp, width = 8, height = 10)


ta3 = merge(ta2, tv[,c(1,2,8)], by = c("chr", "pos"))

grp = dplyr::group_by(ta3, sid, gid)
ta4 = dplyr::summarise(grp, nvnt = n())
ta5 = ta4[ta4$nvnt > 1,]
ta6 = merge(ta5, tg2, by = 'gid')
grp = dplyr::group_by(ta6, sid, gtype)
ta7 = dplyr::summarise(grp, ngene = n())
ta8 = merge(ta7, ti[,c(1,6)], by.x = 'sid', by.y = 'SampleID')

p1 = ggplot(ta8) +
    geom_bar(aes(x = label, y = ngene, fill = gtype), stat = 'identity', position = position_stack(reverse=T)) +
    #scale_x_continuous(name = '') +
    scale_y_continuous(name = "# Genes with >=2 distinguishing sites (depth >= 5)") +
    scale_fill_brewer(palette = "Set1") +
    coord_flip() +
    theme_bw() +
    theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines")) +
    theme(legend.position = 'top', legend.direction = "horizontal", legend.justification = c(0,0.5), legend.title = element_blank(), legend.key.size = unit(1, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size = 8), legend.background = element_rect(fill=NA, size=0)) +
    theme(axis.title.x = element_text(size = 9)) +
    theme(axis.title.y = element_blank()) +
    theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0, hjust = 1)) +
    theme(axis.text.y = element_text(size = 8, colour = "black", angle = 0))
fp = file.path(diro, "26.ase.gene.pdf")
ggsave(p1, filename = fp, width = 8, height = 8)
#}}}

#{{{# look at indel bias
sid = "BR152"
sids = c("BR006", "BR054", "BR151")

to = data.frame()
for (sid in sids) {
fa = sprintf("%s/25.ase.site/%s.tsv", dirw, sid)
ta = read.table(fa, sep = "\t", as.is = T, header = T)[,c(1,2,3,4,6:8)]
ta1 = ta[ta$depth >= 10 & ta$alt != ".",]
ta2 = cbind(ta1, snp = nchar(ta1$ref)==nchar(ta1$alt))
ta2$snp[ta2$snp == T] = 'SNP'
ta2$snp[ta2$snp == F] = 'InDel'

ta3 = cbind(sid = tl$label[tl$SampleID==sid], ta2[,-c(1:4)])
to = rbind(to, ta3)
}
p1 = ggplot(to) +
    geom_histogram(aes(x = dpr/(dpr+dpa)), bins = 30) +
    scale_x_continuous(name = 'B73 Allele Abundance', limits = c(0,1), expand = c(0,0)) +
    scale_y_continuous(name = "# Polymorphic Sites") +
    #scale_color_manual(name = "", values = cols) +
    facet_wrap(snp ~ sid, nrow = 2, scale = 'free') +
    theme_bw() +
    theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines")) +
    #theme(legend.position = c(0.4, 0.8), legend.direction = "vertical", legend.justification = c(0,0.5), legend.title = element_blank(), legend.key.size = unit(1, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size = 8), legend.background = element_rect(fill=NA, size=0)) +
    theme(axis.title.x = element_text(size = 9)) +
    theme(axis.title.y = element_text(size = 9)) +
    theme(axis.text.x = element_text(size = 8, colour = "black", angle = 0)) +
    theme(axis.text.y = element_text(size = 8, colour = "black", angle = 0))
fp = file.path(diro, "25.ase.bias.pdf")
ggsave(p1, filename = fp, width = 12, height = 8)
#}}}

