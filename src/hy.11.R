source('functions.R')
diri = file.path(dird, '20_mass_raw')
dirw = file.path(dird, '21_mass')
gcfg = read_genome_conf()
t_syn = read_syn(gcfg) %>% arrange(gid, -ftype) %>% group_by(gid) %>% slice(1) %>% ungroup()
fh = file.path(dird, '01_exp_design', '61.BR5.meta.tsv')
th = read_tsv(fh)

#{{{ read leaf table
sids = th %>% filter(Tissue=='leaf_V2', Rep>1) %>%
    mutate(Genotype=factor(Genotype, levels=gts10)) %>% arrange(Rep,Genotype) %>%
    pull(SampleID)
fi = file.path(diri, 'hybrid.xlsx')
ti = read_xlsx(fi, sheet=2)
colnames(ti)[1:50] = str_c(rep(str_c('batch', 1:5), each=10),
    c(rep(sids[1:10],2), rep(sids[11:20],2), sids[21:30]), sep="_")
ti = ti %>% select(group=Group,gid=Accession,starts_with("batch")) %>%
    group_by(group) %>% mutate(gidx = 1:n()) %>%
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
    select(group, leader, tid, gid, iso, note, everything())
t_pro %>% count(leader)
t_pro %>% count(leader, note)
#}}}

t_pro_cnt = t_pro %>% filter(leader, note=='shared') %>%
    select(-group, -leader, -tid, -iso, -note) %>%
    gather(SampleID, ReadCount, -gid) %>% filter(str_detect(SampleID, '^batch1'))
th2 = t_pro_cnt %>% distinct(SampleID) %>%
    separate(SampleID, c("batch",'sid'), by='_', remove=F) %>%
    inner_join(th, by=c('sid'='SampleID'))

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

tm = res$tm
tm = t_pro_cnt %>% rename(CPM=ReadCount)
#{{{ PCA
tw = tm %>% select(SampleID, gid, CPM) %>% mutate(CPM=asinh(CPM)) %>% spread(SampleID, CPM)
t_exp = tm %>% filter(CPM>=1) %>% count(gid) %>% rename(n.exp=n)
gids = t_exp %>% filter(n.exp >= (ncol(tw)-1) * .1) %>% pull(gid)
tt = tw %>% filter(gid %in% gids)
tt[is.na(tt)] = 0
e = tt %>% select(-gid)
dim(e)
pca <- prcomp(e, center = T, scale. = T)
x = pca['rotation'][[1]]
y = summary(pca)$importance
y[,1:5]

xlab = sprintf("PC1 (%.01f%%)", y[2,1]*100)
ylab = sprintf("PC2 (%.01f%%)", y[2,2]*100)
tp = as_tibble(x[,1:5]) %>%
    add_column(SampleID = rownames(x)) %>%
    left_join(th2, by = 'SampleID')
p_pca = ggplot(tp, aes(PC1, PC2)) +
    geom_point(aes(color=Genotype,shape=batch), size=2) +
    geom_text_repel(aes(label=Genotype), size=2.5) +
    scale_x_continuous(name = xlab) +
    scale_y_continuous(name = ylab) +
    scale_shape_manual(values = 0+c(0:5)) +
    scale_color_aaas() +
    otheme(legend.pos='top.right', legend.dir='v', legend.title=T,
           xtitle=T, ytitle=T,
           margin = c(.2,.2,.2,.2))
fp = file.path(dirw, "leaf.pca.pdf")
ggsave(p_pca, filename = fp, width=8, height=8)
#}}}

tm = t_pro_cnt %>% rename(CPM=ReadCount)
#{{{ tSNE
require(Rtsne)
tw = tm %>% select(SampleID, gid, CPM) %>% mutate(CPM=asinh(CPM)) %>%
    spread(SampleID, CPM)
t_exp = tm %>% filter(CPM>=1) %>% count(gid) %>% rename(n.exp=n)
gids = t_exp %>% filter(n.exp >= (ncol(tw)-1) * .1) %>% pull(gid)
tt = tw %>% filter(gid %in% gids)
tt[is.na(tt)] = 0
dim(tt)
tsne <- Rtsne(t(as.matrix(tt[-1])), dims=2, verbose=T, perplexity=3,
              pca = T, max_iter = 1500)

tp = as_tibble(tsne$Y) %>%
    add_column(SampleID = colnames(tt)[-1]) %>%
    left_join(th2, by = 'SampleID')
x.max=max(tp$V1)
p_tsne = ggplot(tp, aes(x=V1,y=V2)) +
    geom_point(aes(color=Genotype,shape=batch), size=2) +
    geom_text_repel(aes(label=Genotype), size=2.5) +
    scale_x_continuous(name = 'tSNE-1') +
    scale_y_continuous(name = 'tSNE-2') +
    scale_shape_manual(values = 15+c(0:5)) +
    scale_color_aaas() +
    otheme(legend.pos='top.right', legend.dir='v', legend.title=T,
           xtitle=T, ytitle=T,
           margin = c(.2,.2,.2,.2)) +
    guides(fill=F)
fp = file.path(dirw, "leaf.tsne.pdf")
ggsave(p_tsne, filename = fp, width=8, height=8)
#}}}

