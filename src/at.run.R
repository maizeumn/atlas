source("functions.R")
sid = 'me99b'
dirw = file.path(dirp, ifelse(sid == 'me99b.m', '41_qc_m', "41_qc"))
genome = ifelse(sid == 'me99b.m', 'Mo17', 'B73')
x = load(file.path(dirg, genome, '55.rda'))

#{{{ read data
diri = '~/projects/rnaseq/data'
fh = sprintf("%s/05_read_list/%s.c.tsv", diri, sid)
th = read_tsv(fh)
fi = file.path(diri, '08_raw_output', sid, 'cpm.rds')
res = readRDS(fi)
tl=res$tl; tm = res$tm
#
th = th %>%
    mutate(Tissue=sprintf("%s|%s",Tissue,Treatment)) %>%
    mutate(Genotype = ifelse(Genotype == 'B73xMo17', 'BxM', Genotype)) %>%
    mutate(Genotype = ifelse(Genotype == 'Mo17xB73', 'MxB', Genotype)) %>%
    filter(Genotype %in% gts) %>%
    select(SampleID,Tissue,Genotype)
tl = tl %>% filter(SampleID %in% th$SampleID)
tm = tm %>% filter(SampleID %in% th$SampleID)
tiss = unique(th$Tissue); genos = unique(th$Genotype)
#}}}

require(MASS)

gids = tm %>% group_by(gid) %>%
    summarise(n.sam = sum(CPM>=1)) %>%
    ungroup() %>% filter(n.sam / nrow(th) >= .7) %>% pull(gid)
length(gids)
#
tw = tm %>% filter(gid %in% gids) %>%
    dplyr::select(gid, SampleID, CPM) %>%
    mutate(CPM = asinh(CPM)) %>%
    spread(gid, CPM) %>%
    inner_join(th, by = 'SampleID') %>%
    dplyr::select(SampleID,Genotype,Tissue,everything())
dim(tw)
tw[1:5,1:5]

r <- lda(formula = Tissue ~ ., data = tw[,-c(1:2)])#, prior = c(1,1,1)/3)
fo = file.path(dirw, 'lda.rds')
saveRDS(r, file=fo)

