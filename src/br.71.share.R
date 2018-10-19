#{{{ head
source("br.fun.r")
dirw = file.path(dirp, "49.coop")
fi = file.path(dirw, "01.master.RData")
x = load(fi)
fd = file.path(dirw, "03.sharing.RData")
x = load(fd)
tis.prop = 0.5
#}}}

#{{{ tissue-specific DE/SPE
dirw = file.path(dirp, "71_output")

to1 = tsh_d %>% filter(ctag == 'pDE', n.tis.tot >= 10) %>%
    transmute(gid = gid, tag = as.character(tsTag)) %>%
    mutate(tag = ifelse(tag=='No data', 'not', tag)) %>%
    mutate(tag = sprintf("%s DE", tag))
to1 %>% count(tag)

to2 = tsh_d %>% filter(ctag == 'pDE', n.tis.tot >= 10) %>%
    transmute(gid = gid, tag = as.character(tag)) %>%
    mutate(tag = ifelse(tag=='non-pDE', 'not DE', tag))
to2 %>% count(tag)

to3 = tsh_d %>% filter(ctag == 'SPE', n.tis.tot >= 10) %>%
    transmute(gid = gid, tag = as.character(tsTag)) %>%
    mutate(tag = ifelse(tag=='No data', 'not', tag)) %>%
    mutate(tag = sprintf("%s SPE", tag))
to3 %>% count(tag)

to4 = tsh_d %>% filter(ctag == 'SPE', n.tis.tot >= 10) %>%
    transmute(gid = gid, tag = as.character(tag)) %>%
    mutate(tag = ifelse(tag=='non-SPE', 'not SPE', tag))
to4 %>% count(tag)

to = to1 %>% 
    inner_join(to2, by = 'gid') %>%
    inner_join(to3, by = 'gid') %>%
    inner_join(to4, by = 'gid')
colnames(to)[2:5] = sprintf("tag%d", 1:4)
fo = file.path(dirw, '01.DE.SPE.specificity.tsv')
write_tsv(to, fo)
#}}}

#{{{
fi = file.path(dirp, "41.qc/10.rep.merged.RData")
x = load(fi)

to1 = t_rep_merged %>% filter(Genotype == 'B73') %>%
    select(Tissue, gid, CPM) %>%
    spread(Tissue, CPM)
dim(to1)
to2 = t_rep_merged %>% filter(Genotype == 'B73') %>%
    select(Tissue, gid, FPKM) %>%
    spread(Tissue, FPKM)
dim(to2)

diro = file.path(dirp, "71_output")
fo1 = file.path(diro, 'CPM_B73.tsv')
write_tsv(to1, fo1)
fo2 = file.path(diro, 'FPKM_B73.tsv')
write_tsv(to2, fo2)
#}}}

#{{{ share w. Briggs group
dirp = '~/projects/maize.expression' 
study = 'briggs'
dirw = file.path(dirp, study, 'data')
diri = file.path(dirp, study, 'data/raw/multiqc_data')
fh = file.path(dirw, '01.reads.tsv')
th = read_tsv(fh)
gts = c("B73", "Mo17", "B73xMo17")
tissues = sort(unique(th$Tissue))
th = th %>% filter(Genotype %in% gts, ! SampleID %in% c('BR207', 'BR230', "BR235"))

fi = file.path(dirw, '20.rc.norm.rda')
x = load(fi)
x

diro = '~/projects/briggs/data/71_output'
fo = file.path(diro, 'samples.tsv')
write_tsv(th, fo)

to = tm %>% select(gid, SampleID, ReadCount) %>%
    spread(SampleID, ReadCount)
fo = file.path(diro, 'raw_readcount.tsv')
write_tsv(to, fo)

to = tm %>% select(gid, SampleID, CPM) %>%
    spread(SampleID, CPM)
fo = file.path(diro, 'cpm.tsv')
write_tsv(to, fo)

to = tm %>% select(gid, SampleID, FPKM) %>%
    spread(SampleID, FPKM)
fo = file.path(diro, 'fpkm.tsv')
write_tsv(to, fo)
#}}}

#{{{ TF expression w. Erika
dirw = file.path(dirp, "71_output")
fi = file.path(dirw, 'tf_v4.txt')
tids = read_tsv(fi, col_names = F) %>% pull(X1)

to = tm %>% filter(gid %in% tids) %>%
    select(Tissue, gid, B73, Mo17, BxM) %>%
    gather(Genotype, CPM, -Tissue, -gid) %>%
    select(gid, Tissue, Genotype, CPM)
fo = file.path(dirw, "tf_cpm.tsv")
write_tsv(to, fo)
#}}}

