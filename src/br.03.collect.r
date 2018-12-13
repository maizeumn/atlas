source("functions.R")
dirw = file.path(dirp, "03_collect")

#{{{ check sequence paths, create read file for running python pipeline
fi = file.path(dirw, "../01_exp_design/exp.design.tsv")
td = read_tsv(fi, col_names = c("SampleID", "Tissue", "Genotype", "Treatment", "pool", "index"))

fi = file.path(dirw, "sequence_paths.txt")
seqpaths = read.table(fi)$V1
tq = tibble()
for (seqpath in seqpaths) {
    fnames = list.files(seqpath)
    idxs = which(endsWith(fnames, ".fastq") | endsWith(fnames, ".fastq.gz"))
    fnames = fnames[idxs]
    x = strsplit(fnames, split = "_")
    pool = sapply(x, "[", 2)
    pool = as.integer(substr(pool, 5, nchar(pool)))
    index = sapply(x, "[", 4)
    index = as.integer(substr(index, 6, nchar(index)))
    pair = sapply(x, "[", 6)
    fname = sprintf("%s/%s", seqpath, fnames)
    sid = sprintf("pool%02d_index%02d", pool, index)
    tos = tibble(fname=fname, pair=pair, sid=sid)
    tq = rbind(tq, tos)
}
tq = spread(tq, pair, fname)
sum(is.na(tq))

td2 = td %>% separate(index, c("set", "index"), sep = "_") %>%
    mutate(poolidx = sprintf("pool%02d_%s", as.integer(str_sub(pool,5,-1)), str_to_lower(index))) %>%
    select(-pool, -set, -index) %>%
    left_join(tq, by = c("poolidx" = "sid"))
td2 %>% filter(is.na(R1))

td3 = td2 %>% select(-poolidx) %>%
    mutate(R1e = file.exists(R1), R2e = file.exists(R2),
           gz1 = (str_sub(R1, -3, -1) == ".gz" ),
           gz2 = (str_sub(R2, -3, -1) == ".gz"))
sum(!td3$R1e + !td3$R2e)
sum(!td3$gz1)
sum(!td3$gz2)
stopifnot(identical(td3$gz1, td3$gz2))

th = td3 %>% mutate(gz = gz1) %>% select(-R1e, -R2e, -gz1, -gz2)
fo = file.path(dirw, "01.reads.tsv")
write_tsv(th, fo)
#}}}


