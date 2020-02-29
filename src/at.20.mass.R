require(DESeq2)
require(edgeR)
source('functions.R')
diri = file.path(dird, '20_mass_raw')
dirw = file.path(dird, '21_mass')
gcfg = read_genome_conf()
t_syn = read_syn(gcfg) %>% arrange(gid, -ftype) %>% group_by(gid) %>% slice(1) %>% ungroup()

#{{{ #[obsolete] read protein table and pre-process
th = read_xlsx(file.path(diri, 'TMT_labels.xlsx')) %>%
    mutate(SampleID = str_c("ittc", 1:10, sep='')) %>%
    select(SampleID, tmt, Genotype=genotype, Replicate=replicate)
fh = file.path(dird, '01_exp_design', '01.BR0.meta.tsv')
gtm = c("B73"="B73","Mo17"="Mo17",'B73xMo17'='BxM','Mo17xB73'='MxB')
th = read_tsv(fh) %>% select(1:5) %>% mutate(Genotype = gtm[Genotype])

#{{{ read Zhouxin's table
fi = file.path(diri, 'spectro_mills.xlsx')
ti1 = read_xlsx(fi, "root") %>%
    select(group=Group,gid=Accession,root, starts_with("BR")) %>%
    gather(SampleID, itt, -group, -gid)
ti2 = read_xlsx(fi, "embryo") %>%
    select(group=Group,gid=Accession,embryo,embryo_2, starts_with("BR")) %>%
    gather(SampleID, itt, -group, -gid)
ti3 = read_xlsx(fi, "coleoptile_tip") %>%
    select(group=Group,gid=Accession,coleoptile_tip, starts_with("BR")) %>%
    gather(SampleID, itt, -group, -gid) %>%
    arrange(SampleID, group, -itt) %>%
    group_by(SampleID, group) %>%
    summarise(gid = gid[1], itt = itt[1]) %>% ungroup()

ti = ti1 %>% bind_rows(ti2) %>% bind_rows(ti3) %>% filter(str_detect(gid, '^Z'))
ti %>% count(SampleID) %>% print(n=50)

tm = ti %>% group_by(SampleID) %>% mutate(total_itt = sum(itt)) %>%
    mutate(CPM = itt/total_itt * 1e6) %>% ungroup() %>% select(-itt, -total_itt)
#}}}

#{{{ read maxquant
fd = '~/Downloads/testMS/combined/txt/proteinGroups.txt'
td = read_tsv(fd)
colnames(td) = c('pids','mpids','npep_all','npep_ru','npep_u',
    'fhead','npro','peps','peps_ru','peps_u','scov','scov_ru','scov_u',
    'molweight','seqlen','seqlens','qval','score',
    str_c('ittc', 1:10, sep=''),
    str_c('ittu', 1:10, sep=''),
    str_c('itt_cnt', 1:10, sep=''),
    'itt', 'msms_cnt', 'site_only','reverse','contam',
    'idx','pep_ids','isRazor','mod_pep_ids','evidence_ids', 'msms_ids',
    'msms_best','oxi_sites', 'oxi_site_pos')
td2 = td %>% select(idx, pids, mpids, npro, peps, peps_ru, scov, scov_ru, seqlen, seqlens)
tm = td %>% select(pids, starts_with("ittc")) %>%
    filter(str_detect(pids, "^Z")) %>%
    rename(gid=pids) %>%
    gather(SampleID, CPM, -gid)
#}}}

#{{{ tSNE
require(Rtsne)
tw = tm %>% select(SampleID, gid, CPM) %>% mutate(CPM=asinh(CPM)) %>% spread(SampleID, CPM)
t_exp = tm %>% group_by(gid) %>% summarise(n.exp = sum(CPM>=1))
gids = t_exp %>% filter(n.exp >= (ncol(tw)-1) * .5) %>% pull(gid)
tt = tw %>% filter(gid %in% gids)
tt[is.na(tt)] = 0
dim(tt)
tsne <- Rtsne(t(as.matrix(tt[-1])), dims=2, verbose=T, perplexity=4,
              pca = T, max_iter = 1000)

tp = as_tibble(tsne$Y) %>%
    add_column(SampleID = colnames(tt)[-1]) %>%
    left_join(th, by = 'SampleID') %>%
    replace_na(list(Tissue='unk', Genotype='unk'))
x.max=max(tp$V1)
p_tsne = ggplot(tp, aes(x=V1,y=V2)) +
    geom_point(aes(color=Tissue,shape=Genotype), size=2) +
    scale_x_continuous(name = 'tSNE-1') +
    scale_y_continuous(name = 'tSNE-2') +
    scale_shape_manual(values = c(0:5)) +
    scale_color_aaas() +
    otheme(legend.pos='top.right', legend.dir='v', legend.title=T,
           xtitle=T, ytitle=T,
           margin = c(.2,.2,.2,.2)) +
    theme(axis.ticks.length = unit(0, 'lines')) +
    guides(fill=F)
fp = file.path(dirw, "test.tsne.pdf")
ggsave(p_tsne, filename = fp, width=8, height=8)
#}}}

#{{{ mess with Zhouxi's protein database
fh = file.path(diri, 'seq.tsv')
th = read_tsv(fh)

tp = th %>% separate(seqid, c('pid','note'), sep='[\\|]', extra='merge') %>%
    select(-desc) %>%
    separate(pid, c('gid', 'suf', 'note1', 'note2'), sep="_", extra='merge', fill='right') %>%
    mutate(note = ifelse(!is.na(note1), note1, note)) %>%
    mutate(note = ifelse(note %in% c("reverse",'pep','B73','Mo17'), note, 'decoy'))
tp %>% count(note)

tp1 = tp %>% filter(note != 'decoy') %>% filter(!(is.na(note1) & note=='pep2')) %>% mutate(suf=str_sub(suf,2))
x = tp1 %>% count(gid, suf)
x1 = x %>% filter(n==1) %>% select(gid, suf)
tp1 %>% inner_join(x1, by=c('gid','suf')) %>% count(note1, note2, note)
x2 = x %>% filter(n==2) %>% select(gid, suf)
tp1 %>% inner_join(x2, by=c('gid','suf')) %>% count(note1, note2, note)
x3 = x %>% filter(n==3) %>% select(gid, suf)
tp1 %>% inner_join(x3, by=c('gid','suf')) %>% count(note1, note2, note)
#}}}
#}}}

#{{{ process the 11 tissue table and store
fi = file.path(diri, 'tissue11.xlsx')
x0 = read_xlsx(fi, n_max=1)
ti = read_xlsx(fi, col_names=as.character(x0[1,]), skip=2) %>%
    rename(grp=1,score=2,n_pep_u=3,n_spec=4,pro_cov=5,tid=6,tids=7,pro=8)

ti1 = ti %>% select(grp,score,n_pep_u,n_spec,pro_cov, tid, tids, pro)
colmap = tibble(sid = colnames(ti)[-c(1:8)]) %>%
    separate(sid,c('sid0','batch'),sep='[.][.][.]',fill='right',remove=F) %>%
    mutate(batch = rep(1:15, each=10))
ti2 = ti %>% select(-grp,-score,-n_pep_u,-n_spec,-pro_cov,-tids,-pro) %>%
    gather(sid, itt, -tid) %>%
    inner_join(colmap, by='sid') %>% select(-sid) %>% rename(sid=sid0)

tn = ti1 %>% select(tid) %>%
    filter(str_detect(tid, "^(Zm)|(Ze)")) %>%
    filter(!str_detect(tid, "reverse$")) %>%
    separate(tid, c('gid', 'iso', 'note1', 'note2', 'note3'), sep="_", extra='merge', fill='right', remove=F) %>%
    replace_na(list(iso='',note1='shared',note2='',note3='')) %>%
    mutate(note1 = ifelse(note1=='shared', note1, str_c(note1,note2,sep="_"))) %>%
    mutate(iso = str_replace(iso, '^P', 'T')) %>%
    select(gid, iso, tid, note=note1) %>%
    arrange(gid, note, iso, tid) %>%
    group_by(gid, note) %>% slice(1) %>% ungroup()
tn %>% count(gid) %>% print(n=50)
tn %>% count(note)

tg = tn %>% inner_join(ti1, by='tid')
itt = ti2 %>% inner_join(tn[,c('tid','gid','note')], by='tid') %>%
    select(sid, batch, gid, note, itt)

fh = file.path(dirw, 'meta.xlsx')
th = read_xlsx(fh)
th = itt %>% distinct(sid, batch) %>% left_join(th, by=c('sid')) %>%
    mutate(tissue = ifelse(is.na(tissue), sid, tissue)) %>%
    mutate(genotype = ifelse(is.na(genotype), 'pool', genotype)) %>%
    rename(sid0 = sid) %>%
    mutate(sid=sprintf("batch%02d_%s", batch, sid0)) %>%
    select(sid, tissue, stage, genotype, rep, batch, sid0)

itt2 = itt %>%
    mutate(sid=sprintf("batch%02d_%s", batch, sid)) %>%
    select(-batch)

x = tg %>% select(gid,note) %>% spread(note, -gid) %>% count(shared, B73_only, Mo17_only, Mo17_unique)
fo = file.path(dirw, '02.cnt.tsv')
write_tsv(x, fo)

res = list(tg=tg, th=th, itt=itt2)
fo = file.path(dirw, '01.rds')
saveRDS(res, fo)
#}}}

#{{{ run DE test, store
fi = file.path(dirw, '01.rds')
res = readRDS(fi)

#{{{ tSNE
t_sne = res$itt %>% mutate(gid = str_c(gid, note, sep='-')) %>%
    select(sid, gid, CPM=itt)
require(Rtsne)
tw = t_sne %>% select(sid, gid, CPM) %>% mutate(CPM=asinh(CPM)) %>% spread(sid, CPM)
t_exp = t_sne %>% group_by(gid) %>% summarise(n.exp = sum(CPM>=1))
gids = t_exp %>% filter(n.exp >= (ncol(tw)-1) * .5) %>% pull(gid)
tt = tw %>% filter(gid %in% gids)
tt[is.na(tt)] = 0
dim(tt)
tsne <- Rtsne(t(as.matrix(tt[-1])), dims=2, verbose=T, perplexity=6,
              pca = T, max_iter = 1200)

tp = as_tibble(tsne$Y) %>%
    add_column(sid = colnames(tt)[-1]) %>%
    left_join(res$th, by = 'sid') %>%
    mutate(clu = sprintf("b%02d_%s", batch, tissue)) %>%
    mutate(txt = ifelse(genotype=='pool', clu, ''))
x.max=max(tp$V1)
p_tsne = ggplot(tp, aes(x=V1,y=V2)) +
    geom_point(aes(shape=genotype, color=genotype), size=2) +
    geom_text_repel(aes(label=txt), size = 3, alpha = .8) +
    stat_ellipse(aes(fill=clu), linetype = 1, alpha = .4) +
    scale_x_continuous(name = 'tSNE-1') +
    scale_y_continuous(name = 'tSNE-2') +
    scale_shape_manual(values = c(0:5)) +
    scale_color_aaas() +
    otheme(legend.pos='top.right', legend.dir='v', legend.title=T,
           xtitle=T, ytitle=T,
           margin = c(.2,.2,.2,.2)) +
    theme(axis.ticks.length = unit(0, 'lines'))
    #guides(color=F)
fp = file.path(dirw, "10.tsne.pdf")
ggsave(p_tsne, filename = fp, width=8, height=8)

pt = ggplot(tp, aes(x=V1,y=V2)) +
    geom_point(aes(shape=genotype, color=genotype), size=2) +
    scale_x_continuous(name = 'tSNE-1') +
    scale_y_continuous(name = 'tSNE-2') +
    scale_shape_manual(values = c(0:5)) +
    scale_color_aaas() +
    facet_wrap(~clu, scale='free', ncol=5) +
    otheme(legend.pos='top.center.out', legend.dir='h', legend.title=T,
           xtitle=T, ytitle=T,
           margin = c(.2,.2,.2,.2)) +
    theme(legend.justification=c(.5,-.3)) +
    theme(axis.ticks.length = unit(0, 'lines'))
    #guides(color=F)
fp = file.path(dirw, "10.tsne2.pdf")
ggsave(pt, filename = fp, width=9, height=6.5)
#}}}

#{{{ read mRNA table
diri = file.path(dird, "19_coop")
fa = file.path(diri, "01.master.rda")
x = load(fa)
#
tismap = c(
    'kernel_14DAP'='kernel',
    'root_0DAP'='root',
    'endosperm_14DAP'='endosperm14D',
    'endosperm_27DAP'='endosperm27D',
    'embryo_27DAP'='embryo',
    'coleoptile_tip'='coleoptile',
    'radicle_root'='radicle_root',
    'embryo_imbibedseed'='seed_imbibed',
    'seedlingroot_11DAS'='seedling_root',
    'seedlingmeristem_11DAS'='seedling_meristem',
    'ear_v14'='ear'
)
tm %>% count(Tissue) %>% print(n=30)
t_rna = tm %>% mutate(Tissue = as.character(Tissue)) %>%
    mutate(Tissue=ifelse(Tissue %in% names(tismap), tismap[Tissue], 'unk')) %>%
    filter(Tissue != 'unk') %>%
    select(tissue=Tissue,gid, B73,Mo17,BxM, silent, pDE, hDE, log2mb, log2fm, Dom, Reg1, Reg2)
t_rna %>% count(tissue)
#}}}

pro_cnt = res$itt %>% filter(note=='shared') %>% select(SampleID=sid, gid, ReadCount=itt)
pro_h = pro$th %>% select(SampleID=sid, Tissue=tissue, Genotype=genotype, Replicate=rep, batch) %>%
    mutate(Genotype=str_replace(Genotype, '^B$', 'B73')) %>%
    mutate(Genotype=str_replace(Genotype, '^M$', 'Mo17')) %>%
    filter(Genotype != 'pool')

#{{{ normalize and DE test
ti = pro_cnt %>% mutate(ReadCount = as.integer(cap_bigint(ReadCount)))
res = readcount_norm(ti)
tmm = res$tm %>% inner_join(pro_h, by='SampleID') %>%
    group_by(batch,Tissue,Genotype,gid) %>%
    summarise(ReadCount=mean(ReadCount), nRC=mean(nRC), rCPM=mean(rCPM), CPM=mean(CPM)) %>% ungroup()
tm1 = res$tm %>% mutate(nRC = cap_bigint(nRC)) %>%
    inner_join(pro_h[,c(1,2,5)], by = 'SampleID') %>%
    group_by(batch,Tissue) %>% nest()
th1 = pro_h %>% inner_join(res$tl, by = 'SampleID') %>%
    group_by(batch,Tissue) %>% nest()
t_de = tm1 %>% inner_join(th1, by = c('batch','Tissue')) %>%
    mutate(res = map2(data.x, data.y, run_de_test)) %>%
    mutate(deseq = map(res, 'deseq'), edger = map(res, 'edger')) %>%
    select(batch, Tissue, deseq, edger)
dd = call_de_dom(t_de, tmm)
#}}}

fo = file.path(dirw, '05.de.dom.rds')
res = list(rna = t_rna, pro = dd)
saveRDS(res, fo)
#}}}


# merge protein and mRNA
fi = file.path(dirw, '05.de.dom.rds')
res = readRDS(fi)

tx1 = res$rna %>% select(tissue,gid,everything())
colnames(tx1)[3:ncol(tx1)] = str_c("m.", colnames(tx1)[3:ncol(tx1)])
tx2 = res$pro %>% rename(tissue=Tissue)
colnames(tx2)[4:ncol(tx2)] = str_c("p.", colnames(tx2)[4:ncol(tx2)])
tx = tx1 %>% inner_join(tx2, by=c('tissue','gid')) %>%
    mutate(pnl = sprintf("b%02d_%s", batch, tissue))

#{{{ pDE
tp = tx %>% count(pnl,m.pDE, p.pDE) %>% filter(!is.na(p.pDE)) %>%
    mutate(m.pDE = as.character(m.pDE), p.pDE = as.character(p.pDE)) %>%
    replace_na(list(m.pDE='silent')) %>%
    mutate(m.pDE = factor(m.pDE), p.pDE = factor(p.pDE)) %>%
    select(pnl, tag1=m.pDE, tag2=p.pDE, n)
#{{{ proportion plot
tp1 = tp %>% filter(pnl == 'b10_coleoptile')
p1 = cmp_proportion1(tp1, xangle=0, oneline=T, xtitle='mRNA DE btw. B & M',
                    legend.title='protein DE btw. B & M:')
fo = file.path(dirw, "12.pde.1.pdf")
ggsave(p1, file=fo, width=5, height=5)

p1 = cmp_proportion(tp, xangle=0, oneline=T, xtitle='mRNA DE btw. B & M',
                    legend.title='protein DE btw. B & M:')
fo = file.path(dirw, "12.pde.pdf")
ggsave(p1, file=fo, width=10, height=8)
#
tp = tp %>% mutate(tmp = tag1) %>% mutate(tag1=tag2, tag2=tmp)
p1 = cmp_proportion(tp, xangle=0, oneline=T, xtitle='protein DE btw. B & M',
                    legend.title='mRNA DE btw. B & M:')
fo = file.path(dirw, "12.pde.alt.pdf")
ggsave(p1, file=fo, width=10, height=8)
#}}}

#{{{ one density plot
fc = 5
tp = tx %>% filter(pnl == 'b10_coleoptile') %>%
    filter(!is.na(m.pDE), !is.na(p.pDE)) %>%
    #filter(m.pDE != 'non_DE') %>%
    mutate(m.log2mb = trim_range(m.log2mb, -fc, fc)) %>%
    mutate(p.log2mb = trim_range(p.log2mb, -fc, fc))
p2a = ggplot(tp, aes(m.log2mb, p.log2mb)) +
    geom_point(aes(color=m.pDE), size=1) +
    geom_vline(xintercept = 0, linetype = 'dashed') +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    scale_x_continuous(name='mRNA log2(M/B)', limits=c(-5,5)) +
    scale_y_continuous(name='protein log2(M/B)', limits=c(-5,5)) +
    scale_color_npg(name = 'mRNA DE btw. B & M') +
    otheme(xtitle=T,ytitle=T,xtext=T,ytext=T,xtick=T,ytick=T,
           legend.pos='top.left', legend.dir='v', legend.vjust=-.3,
           legend.title=T, margin=c(1.5,.2,.2,.2))
p2 = ggMarginal(p2a, type = 'density', margins = 'both', size = 4,
                groupColour = T, groupFill = T,
                xparams = list(size=0), yparams = list(size=0))
fo = file.path(dirw, "13.pde.density.1.pdf")
ggsave(p2, file=fo, width=6, height=6)
#}}}

#{{{ all density plot
fc = 5
tp = tx %>% filter(!is.na(m.pDE), !is.na(p.pDE), m.pDE != 'non_DE') %>%
    mutate(m.log2mb = trim_range(m.log2mb, -fc, fc)) %>%
    mutate(p.log2mb = trim_range(p.log2mb, -fc, fc))
    mutate(pde = sprintf("m=%s | p=%s", m.pDE, p.pDE))
p2a = ggplot(tp, aes(m.log2mb, p.log2mb)) +
    geom_point(aes(color=pde), size=.5) +
    geom_vline(xintercept = 0, linetype = 'dashed') +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    scale_x_continuous(name='mRNA log2(M/B)', limits=c(-5,5)) +
    scale_y_continuous(name='protein log2(M/B)', limits=c(-5,5)) +
    scale_color_npg(name = 'mRNA DE btw. B & M') +
    facet_wrap(~pnl, ncol=5) +
    otheme(xtitle=T,ytitle=T,xtext=T,ytext=T,xtick=T,ytick=T,
           legend.pos='top.center.out', legend.dir='h', legend.vjust=-.3,
           legend.title=T, margin=c(1.5,.2,.2,.2))
fo = file.path(dirw, "13.pde.density.pdf")
ggsave(p2a, file=fo, width=10, height=7)
#}}}

#{{{ kernel density
require(hexbin)
fc=3
tp = tx %>% filter(!is.na(m.pDE), !is.na(p.pDE)) %>%
    filter(m.pDE != 'non_DE') %>%
    mutate(m.log2mb = trim_range(m.log2mb, -fc, fc)) %>%
    mutate(p.log2mb = trim_range(p.log2mb, -fc, fc))
pg = ggplot(tp, aes(m.log2mb, p.log2mb)) +
    geom_hex(bins=80) +
    scale_fill_viridis() +
    #stat_density_2d(aes(fill = ..density..), geom="raster", contour=F) +
    #scale_fill_distiller(palette=4, direction=1) +
    scale_x_continuous(name='mRNA log2(M/B)', limits=c(-fc,fc)) +
    scale_y_continuous(name='protein log2(M/B)', limits=c(-fc,fc)) +
    geom_vline(xintercept = 0, linetype = 'dashed') +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    #scale_color_npg(name = 'mRNA DE btw. B & M') +
    otheme(xtitle=T,ytitle=T,xtext=T,ytext=T,xtick=T,ytick=T,
           legend.pos='top.center.out', legend.dir='h', legend.vjust=-.3,
           legend.title=T, margin=c(.2,.2,.2,.2))
fo = file.path(dirw, "13.pde.kernel.pdf")
ggsave(pg, file=fo, width=6, height=6)
#}}}

#}}}

#{{{ hDE
#{{{ one tissue
fc=2
tp = tx %>% filter(pnl == 'b10_coleoptile') %>%
    filter(!is.na(m.hDE), m.hDE != '1PL') %>%
    mutate(m.log2fm = trim_range(m.log2fm, -fc, fc)) %>%
    mutate(p.log2fm = trim_range(p.log2fm, -fc, fc))
p = ggplot(tp, aes(m.log2fm, p.log2fm)) +
    geom_point(aes(color=m.hDE), size=1) +
    geom_vline(xintercept = 0, linetype = 'dashed') +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    scale_x_continuous(name='mRNA log2(hybrid/midparent)', limits=c(-fc,fc)) +
    scale_y_continuous(name='protein log2(hybrid/midparent)', limits=c(-fc,fc)) +
    scale_color_npg(name = 'F1 level for non_DE genes') +
    otheme(xtitle=T,ytitle=T,xtext=T,ytext=T,xtick=T,ytick=T,
           legend.pos='top.left', legend.dir='v', legend.title=T)
p2a = ggMarginal(p, type = 'density', margins = 'both', size = 4,
                groupColour = T, groupFill = T,
                xparams = list(size=0), yparams = list(size=0))
#
tp = tx %>% filter(pnl == 'b10_coleoptile') %>%
    filter(!is.na(m.Dom), m.Dom != '1MP') %>%
    mutate(m.log2fm = trim_range(m.log2fm, -fc, fc)) %>%
    mutate(p.log2fm = trim_range(p.log2fm, -fc, fc))
p = ggplot(tp, aes(m.log2fm, p.log2fm)) +
    geom_point(aes(color=m.Dom), size=1) +
    geom_vline(xintercept = 0, linetype = 'dashed') +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    scale_x_continuous(name='mRNA log2(hybrid/midparent)', limits=c(-fc,fc)) +
    scale_y_continuous(name='protein log2(hybrid/midparent)', limits=c(-fc,fc)) +
    scale_color_npg(name = 'F1 level for DE genes') +
    otheme(xtitle=T,ytitle=T,xtext=T,ytext=T,xtick=T,ytick=T,
           legend.pos='top.left', legend.dir='v', legend.title=T)
p2b = ggMarginal(p, type = 'density', margins = 'both', size = 4,
                groupColour = T, groupFill = T,
                xparams = list(size=0), yparams = list(size=0))

fo = file.path(dirw, "15.hde.density.1.pdf")
ggarrange(p2a, p2b,
    nrow = 1, ncol = 2, widths=c(1,1)) %>%
    ggexport(filename = fo, width = 11, height = 6)
#}}}

#{{{ all tissues
fc = 3
tp = tx %>%
    filter(!is.na(m.Dom), m.Dom != 'MP') %>%
    mutate(m.log2fm = trim_range(m.log2fm, -fc, fc)) %>%
    mutate(p.log2fm = trim_range(p.log2fm, -fc, fc))
p = ggplot(tp, aes(m.log2fm, p.log2fm)) +
    geom_point(aes(color=m.Dom), size=.5) +
    geom_vline(xintercept = 0, linetype = 'dashed') +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    scale_x_continuous(name='mRNA log2(hybrid/midparent)', limits=c(-fc,fc)) +
    scale_y_continuous(name='protein log2(hybrid/midparent)', limits=c(-fc,fc)) +
    scale_color_npg(name = 'F1 level for DE genes') +
    facet_wrap(~pnl, ncol=5) +
    otheme(xtitle=T,ytitle=T,xtext=T,ytext=T,xtick=T,ytick=T,
           legend.pos='top.center.out', legend.dir='h', legend.vjust=-.3,
           legend.title=T, margin=c(1.5,.2,.2,.2))
fo = file.path(dirw, "15.hde.density.pdf")
ggsave(p, file=fo, width=10, height=7)
#}}}

#{{{ kernel density
require(hexbin)
fc=3
tp = tx %>% filter(!tissue %in% c("kernel",'endosperm14D','endosperm27D')) %>%
    filter(!is.na(m.Dom), m.Dom != 'MP') %>%
    mutate(m.log2fm = trim_range(m.log2fm, -fc, fc)) %>%
    mutate(p.log2fm = trim_range(p.log2fm, -fc, fc))
pg = ggplot(tp, aes(m.log2fm, p.log2fm)) +
    geom_hex(bins=80) +
    scale_fill_viridis() +
    scale_x_continuous(name='mRNA log2(hybrid/midparent)', limits=c(-fc,fc)) +
    scale_y_continuous(name='protein log2(hybrid/midparent)', limits=c(-fc,fc)) +
    geom_vline(xintercept = 0, linetype = 'dashed') +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    otheme(xtitle=T,ytitle=T,xtext=T,ytext=T,xtick=T,ytick=T,
           legend.pos='top.center.out', legend.dir='h', legend.vjust=-.3,
           legend.title=T, margin=c(.2,.2,.2,.2))
fo = file.path(dirw, "15.hde.kernel.pdf")
ggsave(pg, file=fo, width=6, height=6)
#}}}
#}}}

#{{{ enrich
#{{{ collapse mRNA and protein classes
tx1b = tx2 %>% distinct(batch,tissue) %>% inner_join(tx1, by='tissue')
tx = tx1b %>% left_join(tx2, by=c('batch','tissue','gid')) %>%
    mutate(pnl = sprintf("b%02d_%s", batch, tissue))
fc = 5
DEs = c("DE_B",'DE_M')
ctags = c(
    'm=NA p=NA', 'm!=NA p=NA', 'm=NA p!=NA',
    'm=nDE p=nDE', 'm=DE p=nDE', 'm=nDE p=DE',
    'm=DE = p=DE', 'm=DE != p=DE'
)
tp0 = tx %>% mutate(m.pDE=as.character(m.pDE), p.pDE=as.character(p.pDE)) %>%
    replace_na(list(m.pDE='NA',p.pDE='NA')) %>%
    mutate(m.log2mb = trim_range(m.log2mb, -fc, fc)) %>%
    mutate(p.log2mb = trim_range(p.log2mb, -fc, fc)) %>%
    mutate(ctag = ctags[1]) %>%
    mutate(ctag = ifelse(m.pDE!='NA'&p.pDE=='NA', ctags[2], ctag)) %>%
    mutate(ctag = ifelse(m.pDE=='NA'&p.pDE!='NA', ctags[3], ctag)) %>%
    mutate(ctag = ifelse(m.pDE=='non_DE'&p.pDE=='non_DE', ctags[4], ctag)) %>%
    mutate(ctag = ifelse(m.pDE %in% DEs & p.pDE=='non_DE', ctags[5], ctag)) %>%
    mutate(ctag = ifelse(m.pDE=='non_DE' & p.pDE %in% DEs, ctags[6], ctag)) %>%
    mutate(ctag = ifelse(m.pDE %in% DEs & m.pDE==p.pDE, ctags[7], ctag)) %>%
    mutate(ctag = ifelse(m.pDE %in% DEs & p.pDE %in% DEs & m.pDE!=p.pDE, ctags[8], ctag)) %>%
    mutate(ctag=factor(ctag, levels=ctags)) %>%
    select(pnl, gid, ctag)
#}}}

#{{{ show proportion
tp = tp0 %>% filter(ctag != 'm=NA p=NA') %>%
    rename(tag1=pnl, tag2=ctag) %>% count(tag1, tag2)
p1 = cmp_proportion1(tp,xangle=30, oneline=T,legend.title='',
    margin.l=2)
fo = file.path(dirw, "20.prop.pdf")
ggsave(p1, file=fo, width=8, height=7)
#}}}

#{{{ sharing
pnls = unique(tx$pnl)
pnls4 = c('b02_kernel','b05_endosperm14D','b07_embryo','b09_endosperm27D')
pnls11 = pnls[!pnls %in% pnls4]
tp0a = tp0 %>% filter(pnl %in% pnls11)

tp = tp0a %>% filter(ctag != 'm=NA p=NA') %>%
    count(gid, ctag) %>% rename(ntis=n) %>%
    count(ctag, ntis) %>% rename(ngene=n) %>%
    group_by(ctag) %>% mutate(pgene = ngene/sum(ngene)) %>% ungroup()
pg = ggplot(tp, aes(ntis, pgene)) +
    geom_bar(stat='identity', position='stack', width=.8) +
    geom_text(aes(label=number(ngene)),size=2.5,lineheight=.8, vjust=0) +
    facet_wrap(~ctag, ncol=3) +
    scale_x_continuous(name='number tissues observed',expand=expand_scale(mult=c(.05,.05)),breaks=1:11) +
    scale_y_continuous(name='proportion genes', expand=expand_scale(mult=c(0,.05))) +
    otheme(xtitle=T,ytitle=T,xtext=T,ytext=T,xtick=T,ytick=T,
           margin=c(.2,.2,.2,.2))
fo = file.path(dirw, "22.share.pdf")
ggsave(pg, file=fo, width=8, height=6)
#}}}

#{{{ syn
tp = tp0 %>% inner_join(t_syn, by='gid') %>%
    rename(tag1=ctag, tag2=ftype) %>% count(pnl, tag1, tag2)
p1 = cmp_proportion(tp,xangle=30, oneline=T,legend.title='synteny:',nc=3,
    margin.l=1.5)
fo = file.path(dirw, "31.syn.pdf")
ggsave(p1, file=fo, width=8, height=10)

tp = tp0 %>% filter(pnl == 'b10_coleoptile') %>%
    inner_join(t_syn, by='gid') %>%
    rename(tag1=ctag, tag2=ftype) %>% count(pnl, tag1, tag2)
p1 = cmp_proportion1(tp,xangle=30, oneline=T,legend.title='synteny:',
    margin.l=2)
fo = file.path(dirw, "31.syn.1.pdf")
ggsave(p1, file=fo, width=5, height=5)
#}}}

#}}}

opt = 'Mo17_only'
gids = res$tg %>% filter(note==opt) %>% pull(gid)
res$itt %>% filter(gid %in% gids, note==opt) %>%
    inner_join(res$th, by='sid') %>% replace_na(list(rep=1)) %>%
    mutate(gt = str_c(genotype, rep, sep='_')) %>%
    filter(tissue=='coleoptile') %>%
    select(gid,gt,itt) %>% spread(gt, itt) %>% print(n=150)


#{{{ coleoptile
#{{{ read protein table
fi = file.path(diri, 'spectro_mills.xlsx')
ti = read_xlsx(fi, "coleoptile_tip") %>%
    select(group=Group,gid=Accession,coleoptile_tip, starts_with("BR")) %>%
    group_by(group) %>%
    mutate(gidx = 1:n()) %>%
    mutate(leader= gidx == 1) %>% ungroup() %>% select(-gidx)

tn = ti %>% select(tid=gid) %>%
    separate(tid, c('gid', 'iso', 'note1', 'note2', 'note3'), sep="_", extra='merge', fill='right', remove=F) %>%
    replace_na(list(note1='shared',note2='',note3='')) %>%
    filter(str_detect(gid, "^(GRM|Zm)"), note1!='reverse', note3 != 'reverse') %>%
    mutate(iso = str_replace(iso, '^P', 'T')) %>%
    select(gid, iso, tid, note=note1)
tn %>% count(gid) %>% print(n=50)
tn %>% count(note)

t_pro = ti %>% rename(tid=gid) %>% inner_join(tn, by='tid') %>%
    arrange(gid, note, -leader, group) %>%
    group_by(gid, note) %>%
    slice(n=1) %>% ungroup() %>% arrange(group, leader) %>%
    select(group, leader, tid, gid, iso, note, pool='coleoptile_tip', everything())
t_pro %>% count(leader)
t_pro %>% count(leader, note)
#}}}

#{{{ tSNE
t_sne = t_pro %>% filter(leader) %>%
    mutate(tid=str_c(gid, iso, note, sep='_')) %>% mutate(gid=tid) %>%
    select(-group, -leader, -tid, -iso, -note) %>% gather(SampleID, CPM, -gid)
require(Rtsne)
tw = t_sne %>% select(SampleID, gid, CPM) %>% mutate(CPM=asinh(CPM)) %>% spread(SampleID, CPM)
t_exp = t_sne %>% group_by(gid) %>% summarise(n.exp = sum(CPM>=1))
gids = t_exp %>% filter(n.exp >= (ncol(tw)-1) * .5) %>% pull(gid)
tt = tw %>% filter(gid %in% gids)
tt[is.na(tt)] = 0
dim(tt)
tsne <- Rtsne(t(as.matrix(tt[-1])), dims=2, verbose=T, perplexity=3,
              pca = T, max_iter = 1500)

tp = as_tibble(tsne$Y) %>%
    add_column(SampleID = colnames(tt)[-1]) %>%
    left_join(th, by = 'SampleID') %>%
    replace_na(list(Tissue='unk', Genotype='pool'))
x.max=max(tp$V1)
p_tsne = ggplot(tp, aes(x=V1,y=V2)) +
    geom_point(aes(color=Genotype,shape=Genotype), size=2) +
    scale_x_continuous(name = 'tSNE-1') +
    scale_y_continuous(name = 'tSNE-2') +
    scale_shape_manual(values = c(0:5)) +
    scale_color_aaas() +
    otheme(legend.pos='top.right', legend.dir='v', legend.title=T,
           xtitle=T, ytitle=T,
           margin = c(.2,.2,.2,.2)) +
    theme(axis.ticks.length = unit(0, 'lines')) +
    guides(fill=F)
fp = file.path(dirw, "ctip.tsne.pdf")
ggsave(p_tsne, filename = fp, width=4, height=4)
#}}}

#{{{ read mRNA table
diri = file.path(dird, "19_coop")
fa = file.path(diri, "01.master.rda")
x = load(fa)

t_rna = tm %>% filter(Tissue == 'coleoptile_tip') %>%
    select(gid, B73,Mo17,BxM, silent, pDE, hDE, log2mb, log2fm, Dom, Reg1, Reg2)
#}}}

# merge protein and mRNA
tx = t_pro %>% arrange(gid, -leader) %>%
    group_by(gid) %>% slice(1) %>% ungroup() %>%
    mutate(protein=ifelse(leader,'leader','non_leader')) %>%
    right_join(t_rna, by='gid') %>% replace_na(list(protein='none')) %>%
    left_join(t_syn, by='gid')

#{{{ CPM distri
tx %>% count(protein, silent) %>%
    group_by(protein) %>% mutate(n_tot = sum(n)) %>%
    mutate(prop = n/n_tot) %>%
    mutate(lab = str_c(number(n,big.mark=','), ' (', percent(prop), ')')) %>%
    ungroup() %>% select(protein, n_tot, silent, lab) %>% spread(silent, lab)

tp = tx %>% mutate(CPM=asinh(B73))
tpx = tibble(cpm = c(0,1,10,100,1e3,1e4,1e5), brk=asinh(cpm))
p = ggplot(tp) +
    geom_histogram(aes(CPM, fill=protein), position='dodge', alpha=.8, bins=50) +
    geom_vline(xintercept = asinh(1), color='black', size=1, linetype='dashed') +
    scale_x_continuous(name='CPM', breaks=tpx$brk, labels=tpx$cpm, expand=expand_scale(mult=c(.04,.01))) +
    scale_y_continuous(name='number genes', expand=expand_scale(mult=c(.01,.01))) +
    scale_color_aaas() +
    scale_fill_aaas(labels=c('detected: leader', 'detected: non_leader', 'not detected')) +
    otheme(legend.pos='top.right', legend.dir='v', legend.title=T,
           xtitle=T, ytitle=T, xtext=T, xtick=T, ytext=T, ytick=T,
           margin = c(.2,.2,.2,.2), nostrip=F) +
    facet_zoom(y = protein != 'none')
fp = file.path(dirw, "ctip.cpm.hist.pdf")
ggsave(p, filename = fp, width=8, height=6)

rnamap = c('0'='expressed', '1'='not expressed')
rnas = as.character(rnamap)
promap = c('leader'='translated (leader)', 'non_leader'='translated (non_leader)', 'none'='not translated')
pros = as.character(promap)
tp = tx %>%
    mutate(rna = rnamap[as.character(silent)]) %>%
    mutate(pro = promap[protein]) %>%
    mutate(rna = factor(rna, levels=rnas)) %>%
    mutate(pro = factor(pro, levels=pros)) %>%
    count(rna, pro) %>% select(tag1=pro, tag2=rna, n)

p1 = cmp_proportion(tp, xtitle='protein expression', legend.title='mRNA:', xangle=0, barwidth=.7, oneline=T)
fo = file.path(dirw, "ctip.cpm.pdf")
ggsave(p1, file=fo, width=3.5, height=4)
#}}}

#{{{ syn proportion
rnamap = c('0'='expressed', '1'='not expressed')
rnas = as.character(rnamap)
promap = c('leader'='translated (leader)', 'non_leader'='translated (non_leader)', 'none'='not translated')
pros = as.character(promap)
rna_pros = crossing(rna=rnas, pro=pros) %>%
    mutate(rna=factor(rna,levels=rnas), pro=factor(pro,levels=pros)) %>%
    arrange(rna, pro) %>%
    mutate(rna_pro=str_c(rna, pro, sep=' + ')) %>% pull(rna_pro)
tp = tx %>% mutate(rna = rnamap[as.character(silent)]) %>%
    mutate(pro = promap[protein]) %>%
    mutate(rna_pro = str_c(rna, pro, sep=' + ')) %>%
    mutate(rna = factor(rna, levels=rnas)) %>%
    mutate(pro = factor(pro, levels=pros)) %>%
    mutate(rna_pro = factor(rna_pro, levels=rna_pros)) %>%
    count(rna, pro, rna_pro, ftype) %>%
    select(tag1 = rna_pro, tag2 = ftype, n)

p1 = cmp_proportion(tp, xangle=10, legend.pos='left', legend.dir='v')
fo = file.path(dirw, 'ctip.syn.pdf')
ggsave(p1, file=fo, width=6, height=6)
tp = tp %>% mutate(tmp=tag1) %>% mutate(tag1=tag2,tag2=tmp)
p1 = cmp_proportion(tp, xangle=10, legend.pos='left', legend.dir='v')
fo = file.path(dirw, 'ctip.syn2.pdf')
ggsave(p1, file=fo, width=5, height=6)
#}}}

# get protein count/itt table
t_pro_cnt = tx %>% filter(protein == 'leader') %>%
    select(gid, starts_with("BR")) %>%
    gather(SampleID, ReadCount, -gid)

#{{{ normalize and DE test
ti = t_pro_cnt %>% mutate(ReadCount = as.integer(cap_bigint(ReadCount)))
res = readcount_norm(ti)
tmm = res$tm %>% inner_join(th, by='SampleID') %>%
    group_by(Tissue,Genotype,gid) %>%
    summarise(ReadCount=mean(ReadCount), nRC=mean(nRC), rCPM=mean(rCPM), CPM=mean(CPM)) %>% ungroup()
tm1 = res$tm %>% mutate(nRC = cap_bigint(nRC)) %>%
    inner_join(th[,1:2], by = 'SampleID') %>%
    group_by(Tissue) %>% nest()
th1 = th %>% inner_join(res$tl, by = 'SampleID') %>%
    group_by(Tissue) %>% nest()
t_de = tm1 %>% inner_join(th1, by = 'Tissue') %>%
    mutate(res = map2(data.x, data.y, run_de_test)) %>%
    mutate(deseq = map(res, 'deseq'), edger = map(res, 'edger')) %>%
    select(Tissue, deseq, edger)
dd = call_de_dom(t_de, tmm)
#}}}

# merge protein and mRNA
tx1 = t_rna %>% mutate(Tissue='coleoptile') %>% select(Tissue,gid,everything())
colnames(tx1)[3:ncol(tx1)] = str_c("m.", colnames(tx1)[3:ncol(tx1)])
tx2 = dd
colnames(tx2)[3:ncol(tx2)] = str_c("p.", colnames(tx2)[3:ncol(tx2)])
tx = tx1 %>% inner_join(tx2, by=c('Tissue', 'gid'))

#{{{ pDE
tp = tx %>% dplyr::count(m.pDE, p.pDE) %>% filter(!is.na(p.pDE)) %>%
    mutate(m.pDE = as.character(m.pDE), p.pDE = as.character(p.pDE)) %>%
    replace_na(list(m.pDE='silent')) %>%
    mutate(m.pDE = factor(m.pDE), p.pDE = factor(p.pDE)) %>%
    select(tag1=m.pDE, tag2=p.pDE, n)
p1 = cmp_proportion(tp, xangle=0, oneline=T, xtitle='mRNA DE btw. B & M',
                    legend.title='protein DE btw. B & M:')
fo = file.path(dirw, "ctip.pde.pdf")
ggsave(p1, file=fo, width=4, height=4)
tp = tp %>% mutate(tmp = tag1) %>% mutate(tag1=tag2, tag2=tmp)
p1 = cmp_proportion(tp, xangle=0, oneline=T, xtitle='protein DE btw. B & M',
                    legend.title='mRNA DE btw. B & M:')
fo = file.path(dirw, "ctip.pde2.pdf")
ggsave(p1, file=fo, width=4, height=4)

tp = tx %>% filter(!is.na(m.pDE), !is.na(p.pDE), m.pDE != '1non_DE') %>%
    mutate(m.log2mb = trim_range(m.log2mb, -5, 5)) %>%
    mutate(p.log2mb = trim_range(p.log2mb, -5, 5))
p2a = ggplot(tp, aes(m.log2mb, p.log2mb)) +
    geom_point(aes(color=m.pDE), size=1) +
    geom_vline(xintercept = 0, linetype = 'dashed') +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    scale_x_continuous(name='mRNA log2(M/B)', limits=c(-5,5)) +
    scale_y_continuous(name='protein log2(M/B)', limits=c(-5,5)) +
    scale_color_npg(name = 'mRNA DE btw. B & M') +
    otheme(xtitle=T,ytitle=T,xtext=T,ytext=T,xtick=T,ytick=T,
           legend.pos='top.left', legend.dir='v', legend.title=T)
p2 = ggMarginal(p2a, type = 'density', margins = 'both', size = 4,
                groupColour = T, groupFill = T,
                xparams = list(size=0), yparams = list(size=0))
fo = file.path(dirw, "ctip.pde.density.pdf")
ggsave(p2, file=fo, width=6, height=6)

range.min=-2; range.max=2
tp = tx %>% filter(!is.na(m.hDE), m.hDE != '1PL') %>%
    mutate(m.log2fm = trim_range(m.log2fm, range.min, range.max)) %>%
    mutate(p.log2fm = trim_range(p.log2fm, range.min, range.max))
p = ggplot(tp, aes(m.log2fm, p.log2fm)) +
    geom_point(aes(color=m.hDE), size=1) +
    geom_vline(xintercept = 0, linetype = 'dashed') +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    scale_x_continuous(name='mRNA log2(hybrid/midparent)', limits=c(range.min, range.max)) +
    scale_y_continuous(name='protein log2(hybrid/midparent)', limits=c(range.min, range.max)) +
    scale_color_npg(name = 'F1 level for non_DE genes') +
    otheme(xtitle=T,ytitle=T,xtext=T,ytext=T,xtick=T,ytick=T,
           legend.pos='top.left', legend.dir='v', legend.title=T)
p2a = ggMarginal(p, type = 'density', margins = 'both', size = 4,
                groupColour = T, groupFill = T,
                xparams = list(size=0), yparams = list(size=0))
#
tp = tx %>% filter(!is.na(m.Dom), m.Dom != '1PL') %>%
    mutate(m.log2fm = trim_range(m.log2fm, range.min, range.max)) %>%
    mutate(p.log2fm = trim_range(p.log2fm, range.min, range.max))
p = ggplot(tp, aes(m.log2fm, p.log2fm)) +
    geom_point(aes(color=m.Dom), size=1) +
    geom_vline(xintercept = 0, linetype = 'dashed') +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    scale_x_continuous(name='mRNA log2(hybrid/midparent)', limits=c(range.min, range.max)) +
    scale_y_continuous(name='protein log2(hybrid/midparent)', limits=c(range.min, range.max)) +
    scale_color_npg(name = 'F1 level for DE genes') +
    otheme(xtitle=T,ytitle=T,xtext=T,ytext=T,xtick=T,ytick=T,
           legend.pos='top.left', legend.dir='v', legend.title=T)
p2b = ggMarginal(p, type = 'density', margins = 'both', size = 4,
                groupColour = T, groupFill = T,
                xparams = list(size=0), yparams = list(size=0))

fo = file.path(dirw, "ctip.hde.density.pdf")
ggarrange(p2a, p2b,
    nrow = 1, ncol = 2, widths=c(1,1)) %>%
    ggexport(filename = fo, width = 11, height = 6)
#}}}

#{{{ protein and mRNA
tp1 = tx %>% dplyr::count(m.pDE) %>% filter(!is.na(m.pDE)) %>%
    mutate(tag1='mRNA') %>% rename(tag2=m.pDE)
tp2 = tx %>% dplyr::count(p.pDE) %>% filter(!is.na(p.pDE)) %>%
    mutate(tag1='protein') %>% rename(tag2=p.pDE)
tp = rbind(tp1, tp2) %>%
    mutate(tag1 = as.character(tag1), tag2 = as.character(tag2)) %>%
    mutate(tag1 = factor(tag1), tag2 = factor(tag2))
p1 = cmp_proportion(tp, xangle=0, oneline=T,
                    legend.title='DE btw. B & M:')
fo = file.path(dirw, "ctip.mp.pde.pdf")
ggsave(p1, file=fo, width=3, height=4)

tags = c(levels(tx$p.hDE), levels(tx$p.Dom))
tp1 = tx %>% rename(tag1=m.pDE) %>%
    mutate(tag2=ifelse(tag1=='non_DE',as.character(m.hDE),as.character(m.Dom))) %>%
    count(tag1, tag2) %>% filter(!is.na(tag1)) %>%
    mutate(tag1 = str_c("mRNA", " | ", as.character(tag1)))
tp2 = tx %>% rename(tag1=p.pDE) %>%
    mutate(tag2=ifelse(tag1=='non_DE',as.character(p.hDE),as.character(p.Dom))) %>%
    count(tag1, tag2) %>% filter(!is.na(tag1)) %>%
    mutate(tag1 = str_c("protein", " | ", as.character(tag1)))
tp = rbind(tp1, tp2) %>%
    mutate(tag1 = factor(tag1), tag2 = factor(tag2, levels=tags))
p1 = cmp_proportion(tp, xangle=10, oneline=T, pal='Set3', legend.pos='left',
                    legend.dir='v',
                    legend.title='DE btw. B & M:')
fo = file.path(dirw, "ctip.mp.hde.pdf")
ggsave(p1, file=fo, width=4, height=5)
#}}}

#{{{ hDE & Dom
tags = tx %>% distinct(p.pDE, p.Dom, p.hDE) %>% filter(!is.na(p.pDE)) %>%
    arrange(p.pDE, p.Dom, p.hDE) %>%
    mutate(tag = ifelse(p.pDE=='non_DE', as.character(p.hDE), as.character(p.Dom))) %>%
    mutate(tag = str_c(as.character(p.pDE), tag, sep=" | ")) %>% pull(tag)
tp = tx %>% filter(!is.na(p.pDE)) %>%
    count(m.pDE, m.hDE, m.Dom, p.pDE, p.hDE, p.Dom) %>%
    mutate(m.hDE=as.character(m.hDE), m.Dom=as.character(m.Dom)) %>%
    mutate(p.hDE=as.character(p.hDE), p.Dom=as.character(p.Dom)) %>%
    mutate(tag1 = ifelse(m.pDE=='non_DE', m.hDE, m.Dom)) %>%
    mutate(tag2 = ifelse(p.pDE=='non_DE', p.hDE, p.Dom)) %>%
    mutate(tag1 = str_c(m.pDE, tag1, sep=' | ')) %>%
    mutate(tag2 = str_c(p.pDE, tag2, sep=' | ')) %>%
    mutate(tag1 = factor(tag1, levels=tags)) %>%
    mutate(tag2 = factor(tag2, levels=tags)) %>%
    select(tag1, tag2, n) %>%
    filter(!is.na(tag1), !is.na(tag2)) %>%
    mutate(tcol=ifelse(n>2000, 'white','black'))

p1 = ggplot(tp, aes(tag1, tag2)) +
    geom_tile(aes(fill=n), color='white') +
    geom_text(aes(label = number(n)), color=tp$tcol, size=2.5) +
    scale_x_discrete(name = 'mRNA abundance', expand=c(0,0)) +
    scale_y_discrete(name = 'protein abundance', expand=c(0,0)) +
    scale_fill_viridis(direction=-1) +
    otheme(xtitle=T, ytitle=T, xtext=T, ytext=T, xtick=T, ytick=T) +
    theme(axis.text.x = element_text(angle=20,hjust=1,vjust=1)) +
    guides(fill=F)
fo = file.path(dirw, "ctip.hde.pdf")
ggsave(p1, file=fo, width=6, height=6)
#}}}
#}}}



