source("br.fun.r")
sid = 'me99b'
#sid = 'me99b.m'
dirw = file.path(dirp, ifelse(sid == 'me99b.m', '41_qc_m', "41_qc"))
genome = ifelse(sid == 'me99b.m', 'Mo17', 'B73')
x = load(file.path(dirg, genome, '55.rda'))

#{{{ # btw-replicate correlation (before sample removal)
fi = file.path(dirw, "03.rep.RData")
x = load(fi)
#
tw = tm %>% select(SampleID, gid, CPM) %>% spread(SampleID, CPM)
t_exp = tm %>% group_by(gid) %>% summarise(n.exp = sum(CPM>=1))
gids = t_exp %>% filter(n.exp >= (ncol(tw)-1) * .7) %>% pull(gid)
e = tw %>% filter(gid %in% gids) %>% select(-gid)
dim(e)

ths = th %>% distinct(Tissue, Genotype)
tp = tibble()
tpx = tibble()
for (i in 1:nrow(ths)) {
    tiss = ths$Tissue[i]
    geno = ths$Genotype[i]
    th2 = th %>% filter(Tissue == tiss, Genotype == geno, Treatment <= 3)
    sids = th2$SampleID
    if(length(sids) == 1) next
    #sids = sids[!sids %in% sids_rm]
    #stopifnot(length(sids) == 3)
    sid_pairs = combn(sids, m = 2)
    tc = as_tibble(t(sid_pairs))
    colnames(tc) = c('sid1', 'sid2')
    tc = tc %>% 
        left_join(th2[,c(1,4)], by = c('sid1' = "SampleID")) %>%
        left_join(th2[,c(1,4)], by = c('sid2' = "SampleID")) %>%
        add_column(tissue = tiss, genotype = geno)
    colnames(tc)[3:4] = c("repx", "repy")
    tp1 = th2 %>%
        transmute(repx = 1:nrow(th2), repy = 1:nrow(th2),
                  tissue = tiss, genotype = geno, 
                  lab = SampleID)
    tpx = rbind(tpx, tp1)
    for (cor.opt in c('pearson', 'spearman')) {
        pccs = c()
        for (j in 1:nrow(tc)) {
            pcc = cor(e[,tc$sid1[j]], e[,tc$sid2[j]], method = cor.opt)
            pccs = c(pccs, pcc)
        }
        td = tc %>% add_column(cor.opt = cor.opt, pcc = pccs) %>%
            mutate(lab = sprintf("%.02f", pcc))
        if(cor.opt == 'pearson') {
            tmp = td$repx; td$repx = td$repy; td$repy = tmp
        }
        tp = rbind(tp, td)
    }
}

p1 = ggplot(tp) +
    geom_tile(mapping = aes(x = repx, y = repy, fill = pcc), width=.8, height=.8) +
    geom_text(aes(x = repx, y = repy, label = lab), size = 3) +
    geom_text(data = tpx, aes(x = repx, y = repy, label = lab), size = 2) +
    #geom_point(data = data.frame(x=1:3,y=1:3), aes(x=x,y=y), shape=7, size=6) +
    scale_x_continuous(position = "top") +
    scale_y_reverse() +
    scale_fill_viridis(name = "Pearson (upper right) and Spearman (lower left) Correlation Coefficient", begin = 0.3) +
    facet_wrap(tissue~genotype, ncol = 8, strip.position = 'top') +
    theme_bw() +
    theme(strip.background = element_blank(), strip.text = element_text(size = 8)) + 
    theme(strip.switch.pad.grid = unit(0, 'lines'), strip.switch.pad.wrap = unit(0, 'lines')) +
    theme(legend.position = c(.5,1), legend.direction = "horizontal", legend.justification = c(.5,-.8), legend.background = element_blank()) +
    theme(legend.key.size = unit(1, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size = 8)) +
    theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(2.5,.5,.5,.5), "lines")) +
    theme(panel.grid = element_blank()) +
    theme(axis.title = element_blank(), axis.text = element_blank())
fp = file.path(dirw, "06.raw.cpm.cor.pdf")
ggsave(p1, filename = fp, width = 10, height = 12)
#}}}

#{{{ # sample removal
sids_rm = c(
    'BR207',
    'BR230', # lib size too small
#    'BR007',
#    'BR147', # low correlation w. reps but saved by edgeR norm.
#    'BR184',
#    'BR032',
#    'BR026',
#    'BR083',
    'BR235'
)

fi = file.path(dirp, "03.collect/32.RData")
x = load(fi)
x
tx = t_raw %>% filter(gid %in% gids, !SampleID %in% sids_rm)

res = readcount_norm(tx, t_gs)
th = res$th; tm = res$tm

t_byrep = tm
fo = file.path(dirw, "10.byrep.RData")
save(th, t_byrep, file = fo)
#}}}

#{{{ read data
diri = '~/projects/maize.expression/data/11_qc'
fi = file.path(diri, sid, '20.rc.norm.rda')
y = load(fi)
th = th %>%
    mutate(Genotype = ifelse(Genotype == 'B73xMo17', 'BxM', Genotype)) %>%
    mutate(Genotype = ifelse(Genotype == 'Mo17xB73', 'MxB', Genotype)) %>%
    filter(Genotype %in% gts)
tl = tl %>% filter(SampleID %in% th$SampleID)
tm = tm %>% filter(SampleID %in% th$SampleID)
tiss = unique(th$Tissue); genos = unique(th$Genotype); treas = unique(th$Treatment)
reps = unique(th$Replicate)
#}}}

#{{{ prepare for hclust and pca 
ths = th %>% distinct(Tissue, Genotype)
tw = tm %>% select(SampleID, gid, CPM) %>% spread(SampleID, CPM)
t_exp = tm %>% group_by(gid) %>% summarise(n.exp = sum(CPM>=1))
gids = t_exp %>% filter(n.exp >= (ncol(tw)-1) * .8) %>% pull(gid)
e = tw %>% filter(gid %in% gids) %>% select(-gid)
dim(e)
#}}}

#{{{ btw-replicate correlation
th %>% 
    dplyr::count(Tissue, Genotype) %>%
    spread(Genotype, n) %>% print(n=23)

ths = th %>% distinct(Tissue, Genotype)
tp = tibble()
tpx = tibble()
for (i in 1:nrow(ths)) {
    tiss = ths$Tissue[i]
    geno = ths$Genotype[i]
    th2 = th %>% filter(Tissue == tiss, Genotype == geno, Replicate <= 3) %>% select(-Treatment)
    sids = th2$SampleID
    if(length(sids) == 1) next
    #sids = sids[!sids %in% sids_rm]
    #stopifnot(length(sids) == 3)
    sid_pairs = combn(sids, m = 2)
    tc = as_tibble(t(sid_pairs))
    colnames(tc) = c('sid1', 'sid2')
    tc = tc %>% 
        left_join(th2[,c(1,4)], by = c('sid1' = "SampleID")) %>%
        left_join(th2[,c(1,4)], by = c('sid2' = "SampleID")) %>%
        add_column(tissue = tiss, genotype = geno)
    colnames(tc)[3:4] = c("repx", "repy")
    tp1 = th2 %>%
        transmute(repx = Replicate, repy = Replicate,
                  tissue = tiss, genotype = geno, 
                  lab = SampleID)
    tpx = rbind(tpx, tp1)
    for (cor.opt in c('pearson', 'spearman')) {
        pccs = c()
        for (j in 1:nrow(tc)) {
            pcc = cor(e[,tc$sid1[j]], e[,tc$sid2[j]], method = cor.opt)
            pccs = c(pccs, pcc)
        }
        td = tc %>% add_column(cor.opt = cor.opt, pcc = pccs) %>%
            mutate(lab = sprintf("%.02f", pcc))
        if(cor.opt == 'pearson') {
            tmp = td$repx; td$repx = td$repy; td$repy = tmp
        }
        tp = rbind(tp, td)
    }
}

rep_label <- function(l) {
    l = sprintf("Rep%g", l)
    parse(text = l)
}
p1 = ggplot(tp) +
    geom_tile(mapping = aes(x = repx, y = repy, fill = pcc), width=.8, height=.8) +
    geom_text(aes(x = repx, y = repy, label = lab), size = 3) +
    geom_text(data = tpx, aes(x = repx, y = repy, label = lab), size = 2.5) +
    #geom_point(data = data.frame(x=1:3,y=1:3), aes(x=x,y=y), shape=7, size=6) +
    scale_x_continuous(position = "top") +
    scale_y_reverse(labels = rep_label) +
    scale_fill_viridis(name = "Pearson (upper right) and Spearman (lower left) Correlation Coefficient", begin = 0.3) +
    facet_wrap(tissue~genotype, ncol = 8, strip.position = 'top') +
    theme_bw() +
    #theme(strip.background = element_blank()) +
    theme(strip.text = element_text(size = 8, lineheight = unit(.8, 'lines'), margin = margin(.1,0,.1,0,'lines'))) + 
    theme(legend.position = c(.5,1), legend.direction = "horizontal", legend.justification = c(.5,-.4), legend.background = element_blank()) +
    theme(legend.key.size = unit(1, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size = 8)) +
    theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(2.5,.5,.5,.5), "lines")) +
    theme(panel.grid = element_blank()) +
    theme(axis.title = element_blank()) +
    theme(axis.text.y = element_text(size=8), axis.text.x = element_blank())
fp = sprintf("%s/07.cpm.cor.pdf", dirw)
ggsave(p1, filename = fp, width = 10, height = 12)

tp1 = tp %>% mutate(sid = sid1) %>% select(sid, cor.opt, pcc)
tp2 = tp %>% mutate(sid = sid2) %>% select(sid, cor.opt, pcc)
to = tp1 %>% bind_rows(tp2) %>%
    group_by(sid, cor.opt) %>%
    summarise(pcc = mean(pcc)) %>% ungroup() %>%
    spread(cor.opt, pcc)
fo = file.path(dirw, '07.cpm.cor.tsv')
write_tsv(to, fo)
#}}}

#{{{ hclust 
cor_opt = "pearson"
#cor_opt = "spearman"
hc_opt = "ward.D"
edist <- as.dist(1-cor(e, method = cor_opt))
ehc <- hclust(edist, method = hc_opt)
tree = as.phylo(ehc)
lnames = ehc$labels[ehc$order]
#
tp = th %>% mutate(taxa = SampleID, lab = SampleID) 
if(length(tiss)>1) tp = tp %>% mutate(lab = sprintf("%s %s", lab, Tissue), lab)
if(length(genos)>1) tp = tp %>% mutate(lab = sprintf("%s %s", lab, Genotype), lab)
if(length(treas)>1) tp = tp %>% mutate(lab = sprintf("%s %s", lab, Treatment), lab)
if(length(reps)>1) tp = tp %>% mutate(lab = sprintf("%s %s", lab, Replicate), lab)
t_hc = tp %>% select(taxa, everything())

thp = mutate(th, 
    lab = sprintf("%s %s %s %s", SampleID, Tissue, Genotype, Replicate))
p1 = ggtree(tree) + 
    #geom_tiplab(size = 4, color = 'black', offset = 0.04) +
    ggplot2::xlim(0, 15) + 
    theme_tree2()
p1 = p1 %<+% thp + 
    geom_tiplab(aes(label = lab), size = 4, offset = 0.04) + 
    #geom_text(aes(color = as.character(gt_ok), label = gt), size = 4, nudge_x = 6, hjust = 0) + 
    scale_color_manual(values = c("black", "royalblue", "tomato"))
fo = sprintf("%s/08.hclust.pdf", dirw, cor_opt, hc_opt)
ggsave(p1, filename = fo, width = 12, height = 30)
#}}}

#{{{ PCA
tw = tm %>% select(SampleID, gid, rFPKM) %>% spread(SampleID, rFPKM)
t_exp = tm %>% group_by(gid) %>% summarise(n.exp = sum(rFPKM>=1))
gids = t_exp %>% filter(n.exp >= (ncol(tw)-1) * .5) %>% pull(gid)
e = tw %>% filter(gid %in% gids) %>% select(-gid)
dim(e)
pca <- prcomp(asinh(e), center = F, scale. = F)
x = pca['rotation'][[1]]
y = summary(pca)$importance
y[,1:5]

xlab = sprintf("PC1 (%.01f%%)", y[2,1]*100)
ylab = sprintf("PC2 (%.01f%%)", y[2,2]*100)
tismap = LETTERS[1:length(tissues23)]
names(tismap) = tissues23
tp = as_tibble(x[,1:5]) %>%
    add_column(SampleID = rownames(x)) %>%
    inner_join(th, by = 'SampleID') %>% 
    mutate(lab = tismap[Tissue])
cols17 = c(brewer.pal(8, 'Dark2'), brewer.pal(9, 'Set1'))
cols4 = c(brewer.pal(4, 'Set1'))
p1 = ggplot(tp) +
    geom_point(aes(x = PC1, y = PC2, shape = Tissue, color = Genotype), size = 3) +
    scale_x_continuous(name = xlab) +
    scale_y_continuous(name = ylab) +
    scale_shape_manual(values = as.character(tismap)) +
    scale_color_manual(values = cols4) +
    theme_bw() +
    theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(.2,.2,.2,.2), "lines")) +
    theme(legend.position = 'right', legend.direction = "vertical", legend.justification = c(0,.5)) +
    theme(legend.title = element_blank(), legend.key.size = unit(1, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size = 9)) +
    guides(shape = guide_legend(ncol = 1, byrow = T)) 
fp = sprintf("%s/09.pca.pdf", dirw)
ggsave(p1, filename = fp, width = 11, height = 10)
#}}}

#{{{ heatmap
tps = th %>% select(-Treatment, -paired) %>%
    mutate(Tissue = factor(Tissue, levels = tissues23),
           Genotype = factor(Genotype, levels = gts)) %>%
    arrange(Tissue, Genotype) %>%
    mutate(x = length(Tissue):1)
mat = 1 - as.matrix(edist)
tp = mat %>% as.data.frame() %>% 
    rownames_to_column(var='s1') %>% as_tibble() %>%
    gather(s2, sim, -s1) %>%
    inner_join(tps[,c('SampleID','x')], by = c('s1'='SampleID')) %>%
    inner_join(tps[,c('SampleID','x')], by = c('s2'='SampleID')) %>%
    rename(x = x.x, y = x.y)
tpx = tps %>% group_by(Tissue) %>%
    summarise(xmin = min(x), xmax = max(x), x = median(x),
              lab = Tissue[1]) %>% ungroup()
cols3 = pal_npg()(3)
names(cols3) = gts
tpg = tps %>% group_by(Tissue, Genotype) %>%
    summarise(xmin = min(x), xmax = max(x), x = median(x),
              lab = Genotype[1]) %>%
    ungroup() %>%
    mutate(col.gt = cols3[Genotype]) %>%
    mutate(lab = ifelse(lab=='BxM', 'F1', str_sub(lab,1,1)))
xt = -7; xg = -.5

p = ggplot(tp) +
    geom_tile(aes(x = x, y = y, fill = sim)) +
    geom_segment(data = tpx, mapping = aes(x=xt,xend=xt,y=xmin,yend=xmax), size = 3) +
    geom_segment(data = tpg, mapping = aes(x=xg,xend=xg,y=xmin,yend=xmax), color = tpg$col.gt, size = 1) +
    geom_text(data=tpg, mapping=aes(x=xg-3.5, y = x, label = lab), color = tpg$col.gt, size = 2, hjust = 0) +
    scale_x_continuous(expand = c(0,0)) + 
    scale_y_continuous(breaks = tpx$x, labels = tpx$lab, expand = c(0,0)) +
    scale_fill_gradientn(colors = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)) +
    #scale_fill_viridis() +
    otheme(legend.pos = 'top.center.out', legend.dir = 'h',
           ytext = T) +
    theme(plot.margin = unit(c(2,.2,.2,.2), "lines")) +
    theme(panel.border = element_blank()) +
    theme(legend.key.size = unit(.8, 'lines'), legend.text = element_text(size = 7))
fo = sprintf("%s/08.heatmap.pdf", dirw)
ggsave(p, file = fo, width = 9, height = 8.5)
#}}}

#{{{ save to 10.rda
tmm = tm %>% inner_join(th, by = 'SampleID') %>%
    group_by(Tissue, Genotype, gid) %>%
    summarise(ReadCount = sum(ReadCount), nRC = sum(nRC),
              rCPM = mean(rCPM), rFPKM = mean(rFPKM), 
              CPM = mean(CPM), FPKM = mean(FPKM)) %>% ungroup()
nrow(tmm)/69
fo = file.path(dirw, '10.rda')
save(th, tl, tm, tmm, file = fo)
#}}}

#{{{ hclust to order tissues
tms = tm %>% filter(Genotype == 'B73') %>%
    group_by(gid) %>%
    summarise(ntissue = sum(CPM >= 1)) %>%
    ungroup() %>%
    filter(ntissue >= .8 * length(unique(tm4$Tissue)))
tmw = tm %>% filter(Genotype == 'B73', gid %in% tms$gid) %>%
    ungroup() %>%
    select(Tissue, gid, CPM) %>%
    mutate(CPM = asinh(CPM)) %>%
    spread(Tissue, CPM)
hcl = hclust(dist(t(tmw[,-1])), method = "ward.D")
tissues = hcl$labels[hcl$order]
tissues
#}}}

#{{{ #tissue distance matrix
tw = tm %>% filter(Genotype == 'B73') %>%
    select(Tissue, gid, FPKM) %>% 
    group_by(Tissue, gid) %>%
    summarise(FPKM = mean(FPKM)) %>% 
    ungroup() %>% spread(Tissue, FPKM)
e_raw = tw[,-1]
n_exp = apply(e_raw, 1, myfunc <- function(x) sum(x>=1))
e = e_raw[n_exp >= ncol(e_raw) * 0.6,]
dim(e)
e.c.dist <- as.dist(1-cor(e, method = cor_opt))
e.c.hc <- hclust(e.c.dist, method = hc_opt)
hc = e.c.hc
tree = as.phylo(e.c.hc)
p1 = ggtree(tree) + 
    geom_tiplab(size = 4, color = 'black', offset = 0.01) +
    ggplot2::xlim(0, 1.5) + 
    theme_tree2()
fo = sprintf("%s/51.tissue.hc.pdf", dirw)
ggsave(p1, filename = fo, width = 8, height = 8)

tdt = as.tibble(as.matrix(e.c.dist)) %>%
    mutate(t1 = labels(e.c.dist)) %>%
    gather(t2, distance, -t1)
fo = sprintf("%s/50.tissue.dist.RData", dirw)
save(tdt, file = fo)
#}}}

#{{{# LDA
require(MASS)

tissue_map = tm$Tissue; names(tissue_map) = tm$SampleID
geno_map = tm$Genotype; names(geno_map) = tm$SampleID

e1 = ti[,-1]
n_noexp = apply(e1, 1, myfunc <- function(x) sum(x<1))
idxs = which(n_noexp < 30)
length(idxs)

tl = asinh(t(ti[idxs,-1]))
colnames(tl) = ti$gid[idxs]
tl = data.frame(Tissue = tissue_map[rownames(tl)], Genotype = geno_map[rownames(tl)], tl, stringsAsFactors = F)
tl[1:5,1:5]

r <- lda(formula = Tissue ~ ., data = tl[,-2])#, prior = c(1,1,1)/3)
r$svd^2/sum(r$svd^2)

r$prior
r$counts
#r$means
#r$scaling
prop.lda = r$svd^2/sum(r$svd^2)
prop.lda

plda <- predict(object = r, newdata = tl)
tp = data.frame(Tissue = tl$Tissue, Genotype = tl$Genotype, lda = plda$x)
tp$Tissue = factor(tp$Tissue, levels = unique(tp$Tissue))
tp$Genotype = factor(tp$Genotype, levels = unique(tp$Genotype))

xlab = sprintf("LD1 (%.01f%%)", prop.lda[1]*100)
ylab = sprintf("LD2 (%.01f%%)", prop.lda[2]*100)
cols = c(brewer.pal(8, 'Dark2'), brewer.pal(9, 'Set1'))

p1 <- ggplot(tp) + 
    geom_point(aes(x = lda.LD1, y = lda.LD3, shape = Genotype, color = Tissue)) +
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
fp = sprintf("%s/01.sample.lda.pdf", dirw)
ggsave(p1, filename = fp, width = 8, height = 7)
#}}}

