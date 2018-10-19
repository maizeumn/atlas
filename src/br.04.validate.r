source('br.fun.r')
require(ape)

dirw = file.path(dirp, "03.collect")

fi = file.path(dirw, "10.stat.RData")
x = load(fi)

#{{{ FPM 
tls = t_raw %>%
    group_by(sid) %>%
    summarise(total_ReadCount = sum(ReadCount))
tl = t_raw %>%
    left_join(tls, by = c('sid' = 'sid')) %>%
    mutate(FPM = ReadCount / total_ReadCount * 1000000) %>%
    select(c(sid, gid, FPM))

#genotyping
tt = t_raw %>%
    filter(nref + nalt >= 30) %>%
    group_by(sid) %>%
    summarise(propb = median(nref/(nref+nalt))) %>%
    mutate(gt = ifelse(propb <= 0.05, "Mo17", 
        ifelse(propb >= 0.95, "B73",
        ifelse(propb >= 0.45 & propb <= 0.55, "B73xMo17",
        ifelse(propb >= 0.30 & propb <= 0.40, "Mo17xB73",
        ifelse(propb >= 0.65 & propb <= 0.75, "B73xMo17", 'unknown')))))) %>%
    replace_na(list(gt = 'unknown'))
table(tt$gt, useNA='ifany')

# read exp design
fd = file.path(dirp, "01.exp.design/exp.design.tsv")
td = read_tsv(fd, col_names = c("SampleID","tissue","genotype","rep","pool","index")) %>%
    separate(index, c("set", "index"), sep = "_") %>%
    mutate(pool = as.integer(substr(pool, 5, nchar(pool))), 
           index = as.integer(substr(index, 6, nchar(index)))) %>%
    mutate(sid = sprintf("pool%02d_index%02d", pool, index)) %>%
    select(c(sid, SampleID, tissue, genotype))

tp = td %>%
    left_join(tt, by = c("sid" = "sid")) %>%
    replace_na(list(gt = 'unknown')) %>%
    mutate(gt_ok = ifelse(gt == genotype, '1', '3'),
           lab = sprintf("%s %s %s", SampleID, tissue, genotype),
           batch = ifelse(SampleID <= 'BR165', '1', '2'))
tp %>% filter(tissue == 'embryo_27DAP', genotype == 'Mo17xB73')
tp$gt_ok[tp$tissue == 'embryo_27DAP' & tp$genotype == 'Mo17xB73'] = 1

# plot hclust tree for inspection
tw = spread(tl, sid, FPM)
e1 = tw[,-1]
dim(e1)

n_exp = apply(e1, 1, myfunc <- function(x) sum(x>=1))
e = e1[n_exp >= ncol(e1) * 0.7, ]
dim(e)

 
cor_opt = "pearson"
#cor_opt = "spearman"
hc_opt = "ward.D"
plot_title = sprintf("dist: %s\nhclust: %s", cor_opt, hc_opt)
e.c.dist <- as.dist(1-cor(e, method = cor_opt))
e.c.hc <- hclust(e.c.dist, method = hc_opt)
hc = e.c.hc
tree = as.phylo(e.c.hc)

p1 = ggtree(tree) + 
    #geom_tiplab(size = 4, color = 'black', offset = 0.04) +
    ggplot2::xlim(0, 21) + 
    theme_tree2()
p1 = p1 %<+% tp + 
    geom_tiplab(aes(color = batch, label = lab), size = 4, offset = 0.04) + 
    geom_text(aes(color = as.character(gt_ok), label = gt), size = 4, nudge_x = 6, hjust = 0) + 
    scale_color_manual(values = c("black", "royalblue", "tomato"))
fo = sprintf("%s/21.validate.pdf", dirw)
ggsave(p1, filename = fo, width = 12, height = 30)
#}}}

#{{{ batch 1 sample mix-up correction
fd = file.path(dirp, "01.exp.design/exp.design.tsv")
td = read_tsv(fd, col_names = c("SampleID","tissue","genotype","rep","pool","index")) %>%
    separate(index, c("set", "index"), sep = "_") %>%
    mutate(pool = as.integer(substr(pool, 5, nchar(pool))), 
           index = as.integer(substr(index, 6, nchar(index)))) %>%
    mutate(sid = sprintf("pool%02d_index%02d", pool, index)) %>%
    select(c(sid, SampleID, tissue, genotype))
colnames(td) = c("sid", "SampleID", "Tissue", "Genotype")

sms1 = c("BR003", "BR006", "BR032")
sms2 = c("BR004", "BR007", "BR029")
for (i in 1:length(sms1)) {
    sm1 = sms1[i]; sm2 = sms2[i]
    idx1 = which(td$SampleID == sm1)
    idx2 = which(td$SampleID == sm2)
    gt1 = td$Genotype[idx1]; gt2 = td$Genotype[idx2]
    td$Genotype[idx1] = gt2; td$Genotype[idx2] = gt1
}

to = as_tibble(cbind(td, Treatment = NA))
to$Tissue = factor(to$Tissue, levels = unique(to$Tissue))
to$Genotype = factor(to$Genotype, levels = unique(to$Genotype))
to = to[order(to$Tissue, to$Genotype, to$SampleID),]

tos = unique(to[,c("Tissue",'Genotype')])
for (i in 1:nrow(tos)) {
    tiss = tos$Tissue[i]; gt = tos$Genotype[i]
    idxs = which(to$Tissue == tiss & to$Genotype == gt)
    to$Treatment[idxs] = 1:length(idxs)
}
fo = file.path(dirw, "31.label.correct.tsv")
write_tsv(to, fo, col_names = T)
#}}}

#{{{ re-generate data table using corrected labels
fl = file.path(dirw, "31.label.correct.tsv")
tl = read_tsv(fl)
x

to = tl %>% inner_join(t_mapping, by = 'sid') %>% select(-sid)
stopifnot(nrow(to) == nrow(t_mapping))
t_mapping = to

to = tl %>% inner_join(t_raw, by = 'sid') %>% select(-sid)
stopifnot(nrow(to) == nrow(t_raw))
t_raw = to

fo = file.path(dirw, "32.RData")
save(t_mapping, t_raw, file = fo)
#}}}

#{{{# batch 2 sample identification
fi = file.path(dirw, "11.read.correct.tsv")
ti = read.table(fi, sep = "\t", header = T, stringsAsFactors = F)[,1:5]

fd = file.path(dirw, "exp.design.tsv")
td = read.table(fd, sep = "\t", header = F, stringsAsFactors = F)
colnames(td) = c("SampleID", "Tissue", "Genotype", "Treatment", "Pool", "Index")
td = cbind(td, cood = sprintf("%s_%s", tolower(td$Pool), substr(td$Index, nchar(td$Index)-1, nchar(td$Index))))

f_gt = file.path(dirw, "../46.ase", "13.gt.tsv")
t_gt = read_tsv(f_gt)

ti2 = merge(ti, td[,c("SampleID","Pool","Index")], by = 'SampleID')
ti3 = merge(ti2, t_gt, by = 'SampleID')
ti4 = ti3 %>%
    filter(SampleID > 'BR165') %>%
    arrange(Tissue, Genotype, Treatment)


tdic = c(
    "BR197" = 'seedlingmeristem_11DAS',
    "BR198" = 'seedlingmeristem_11DAS',
    'BR215' = "seedlingroot_11DAS",
    'BR216' = 'seedlingroot_11DAS',
    'BR187' = 'radicle_root',
    'BR227' = 'radicle_root',
    'BR212' = 'radicle_root',
    'BR214' = 'radicle_root',
    'BR179' = 'seedlingleaf_11DAS',
    'BR180' = 'seedlingleaf_11DAS'
)
df_tdic = data.frame(SampleID = names(tdic), tiss = as.character(tdic), stringsAsFactors = F)
ti5 = merge(ti4, df_tdic, by = 'SampleID', all.x = T)
ti5$tiss[is.na(ti5$tiss)] = ti5$Tissue[is.na(ti5$tiss)]
ti5 = ti5 %>% 
    mutate(tiss_ok = ifelse(Tissue == tiss, '1', '3')) %>%
    arrange(tiss, gt)

fo = file.path(dirw, "15.batch2.tsv")
write.table(ti5, fo, sep = "\t", row.names = F, col.names = T, quote = F)
#}}}

#{{{# batch 2 sample label correction
fi = file.path(dirw, "11.read.correct.tsv")
ti = read.table(fi, sep = "\t", header = T, stringsAsFactors = F)

f_gt = file.path(dirw, "../46.ase", "13.gt.tsv")
t_gt = read_tsv(f_gt)

tdic = c(
    "BR197" = 'seedlingmeristem_11DAS',
    "BR198" = 'seedlingmeristem_11DAS',
    'BR215' = "seedlingroot_11DAS",
    'BR216' = 'seedlingroot_11DAS',
    'BR187' = 'radicle_root',
    'BR227' = 'radicle_root',
    'BR212' = 'radicle_root',
    'BR214' = 'radicle_root',
    'BR179' = 'seedlingleaf_11DAS',
    'BR180' = 'seedlingleaf_11DAS'
)
df_tdic = data.frame(SampleID = names(tdic), tiss = as.character(tdic), stringsAsFactors = F)

ti2 = merge(ti, t_gt, by = 'SampleID')
ti3 = merge(ti2, df_tdic, by = 'SampleID', all.x = T)
ti3$tiss[is.na(ti3$tiss)] = ti3$Tissue[is.na(ti3$tiss)]

to = transmute(ti3, SampleID=SampleID, Species=Species, Tissue=tiss, Genotype=gt, 
    Treatment=Treatment, ReadFile1=ReadFile1, ReadFile2=ReadFile2)
to$Genotype = factor(to$Genotype, levels = unique(to$Genotype))
to$Tissue = factor(to$Tissue, levels = unique(ti3$Tissue))
to = arrange(to, Tissue, Genotype, SampleID)

tos = unique(to[,c("Tissue",'Genotype')])
for (i in 1:nrow(tos)) {
    tiss = tos$Tissue[i]; gt = tos$Genotype[i]
    idxs = which(to$Tissue == tiss & to$Genotype == gt)
    to$Treatment[idxs] = 1:length(idxs)
}

fo = file.path(dirw, "12.read.batch2correct.tsv")
write.table(to, fo, sep = "\t", row.names = F, col.names = T, quote = F)
#}}}

