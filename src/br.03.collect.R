source("functions.R")
dirw = file.path(dird, "03_collect")

#{{{ check sequence paths, create read file for running python pipeline
fi = file.path(dirw, "../01_exp_design/exp.design.tsv")
td = read_tsv(fi, col_names = c("SampleID", "Tissue", "Genotype", "Treatment", "pool", "index"))

fi = file.path(dirw, "sequence_paths.txt")
seqpaths = read_tsv(fi, col_names='fn')$fn

tq1 = tibble(seqdir = seqpaths) %>%
    mutate(fname = map(seqdir, list.files)) %>%
    unnest() %>% filter(str_detect(fname, "\\.fastq(\\.gz)?$")) %>%
    mutate(fpath = file.path(seqdir, fname)) %>%
    separate(fname, c('pre1','pool','nsample','index','s','pair','suf'), sep='_', remove=F) %>%
    mutate(sampleName = sprintf("%s.%s.%s.%s", pre1,pool,nsample,index)) %>%
    mutate(pool = as.integer(str_sub(pool, 5))) %>%
    mutate(index = as.integer(str_sub(index, 6))) %>%
    mutate(sid = sprintf("pool%02d_index%02d", pool, index)) %>%
    select(fpath, pair, sid, sampleName) %>%
    spread(pair, fpath)
#
tq2 = tibble(seqdir = seqpaths) %>%
    mutate(data = map(seqdir, read_msi_fastqc)) %>%
    unnest()
#
tq = tq1 %>% inner_join(tq2, by = 'sampleName') %>% select(-sampleName)

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

to = td3 %>% mutate(gz = gz1, avgLength=avgLength*2) %>%
    select(SampleID, Tissue, Genotype, Treatment, R1, R2, gz, spots, avgLength)
fo = file.path(dirw, "01.reads.tsv")
write_tsv(to, fo)
#}}}


