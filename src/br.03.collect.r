source('br.fun.r')
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

#{{{ collect read mapping stats
diri = "~/scratch/briggs"
fh = file.path(diri, "03.hisat.tsv")
t_mapping = read_tsv(fh)[, -c(2:3,5:8,12)]

# collect raw read count and ASE
diri = "~/scratch/briggs"
gids = NA
to = tibble()
for (i in 1:nrow(t_mapping)) {
    sid = t_mapping$sid[i]

    fh1 = sprintf("%s/26.htseq/%s.txt", diri, sid)
    stopifnot(file.exists(fh1))
    th1 = read_tsv(fh1, col_names = F)
    colnames(th1) = c("gid", "rc")
    ngene = nrow(th1) - 5
    gids0 = th1$gid[-c(ngene+1:5)]
    vals1 = th1$rc[-c(ngene+1:5)]

    fh2 = sprintf("%s/26.htseq/%s.as.txt", diri, sid)
    stopifnot(file.exists(fh2))
    th2 = read_tsv(fh2, col_names = F)
    colnames(th2) = c("gid", "rc")
    stopifnot(ngene == nrow(th2) - 5)
    stopifnot(gids0 == th2$gid[-c(ngene+1:5)])
    vals2 = th2$rc[-c(ngene+1:5)]

    if(i == 1) {
        gids = gids0
    } else {
        stopifnot(identical(gids, gids0))
    }

    fi = sprintf("%s/28.ase/%s.tsv", diri, sid)
    stopifnot(file.exists(fi))
    ti = read_tsv(fi)

    to1 = tibble(sid = sid, gid = gids, ReadCount = vals1, ASReadCount = vals2) %>%
        left_join(ti) %>%
        replace_na(list(nref=0, nalt=0, ncft=0)) %>%
        mutate(nref=as.integer(nref), nalt=as.integer(nalt), ncft=as.integer(ncft))
    to = rbind(to, to1)
    cat(sid, "\n")
}
dim(to)
t_raw = to

fo = file.path(dirw, "10.stat.RData")
save(t_mapping, t_raw, file = fo)
#}}}

