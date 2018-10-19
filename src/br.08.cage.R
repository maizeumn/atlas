source('br.fun.r')
dirw = file.path(dirp, "03.collect")

fh = file.path(dirw, "31.label.correct.tsv")
th = read_tsv(fh) %>%
    mutate(Replicate = Treatment) %>%
    select(SampleID, Tissue, Genotype, Replicate)

fi = file.path(dirw, "51.rna.quantity.tsv")
ti = read_tsv(fi)
colnames(ti) = c("SampleID", "total_conc", "total_mass_ng", "total_mass", 
                 "total_vol", "box_position", "adapter")
sum(ti$SampleID %in% th$SampleID)

fp = file.path(dirp, '41.qc/07.cpm.cor.tsv')
tp = read_tsv(fp) %>%
    transmute(SampleID = sid, PCC = pearson,  SPC = spearman)

to = ti %>% mutate(total_vol = round(total_mass_ng/total_conc)) %>%
    filter(!is.na(adapter)) %>%
    select(-adapter) %>%
    mutate(required_vol = 3000 / total_conc) %>%
    mutate(taken_vol = ifelse(required_vol >= 16.67, 16.67, required_vol)) %>%
    mutate(taken_mass = taken_vol * total_conc / 1000,
           left_vol = total_vol - taken_vol,
           left_mass = total_mass - taken_mass) %>%
    inner_join(th, by = 'SampleID') %>%
    select(SampleID, Tissue, Genotype, Replicate, everything()) %>%
    left_join(tp, by = 'SampleID')

fo = file.path(dirw, "54.rna.quantity.tsv")
write_tsv(to, fo, na = '')

