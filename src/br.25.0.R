source("functions.R")
sid = 'me99b'
#sid = 'me99b.m'
dirw = file.path(dirp, ifelse(sid == 'me99b.m', '49_coop_m', "49_coop"))
genome = ifelse(sid == 'me99b.m', 'Mo17', 'B73')
x = load(file.path(dirg, genome, '55.rda'))
diri = file.path(dirp, ifelse(sid == 'me99b.m', '42_de_m', "42_de"))
fm = file.path(diri, '11.de.dom.rda')
y = load(fm)
fi = file.path(dirw, "../44_ase/11.ase.rda")
x = load(fi)
fd = file.path(dirg, 'B73', 'v32', 'v32.gid.txt')
gids = read_tsv(fd, col_names = 'gid') %>% pull(gid)

#{{{ merge data and build master table
pDEs = c("DE_B", "DE_M", "non_DE")
tm1 = dd %>% 
    mutate(pDE = ifelse(is.na(tag.mb), NA,
                 ifelse(tag.mb == -1, 'DE_B',
                 ifelse(tag.mb == 1, 'DE_M', 'non_DE')))) %>%
    mutate(pDE = factor(pDE, levels = pDEs)) %>%
    mutate(silent = ifelse(is.na(pDE), 1, 0)) %>%
    select(Tissue, gid, B73, Mo17, BxM, silent, log2mb,log2fm, pDE,hDE, Dom,DoA)
tm1 %>% dplyr::count(pDE, Dom, hDE, is.na(DoA)) %>% print(n=24)
tm1 %>% filter(!is.na(pDE), is.nan(DoA) | is.infinite(DoA))

taglst = list(
    pDE = levels(tm1$pDE),
    hDE = levels(tm1$hDE),
    Dom = levels(tm1$Dom),
    Reg1 = levels(ase$Reg1),
    Reg2 = levels(ase$Reg2)
)
tm3 = tm1 %>% 
    left_join(ase, by = c("Tissue", 'gid')) 
nrow(tm3)/23
tm3 %>% dplyr::count(pDE, Reg1) %>% spread(Reg1, n)
tm3 %>% dplyr::count(hDE, Reg2) %>% spread(Reg2, n)
tm3 = tm3 %>% 
    mutate(Reg1 = as.character(Reg1),
           Reg2 = as.character(Reg2)) %>%
    mutate(Reg1 = ifelse(!is.na(pDE) & pDE != 'non_DE' & abs(log2mb)>=1, Reg1, NA),
           Reg2 = ifelse(!is.na(pDE) & pDE == 'non_DE', Reg2, NA)) %>%
    mutate(Reg1 = factor(Reg1, levels = taglst$Reg1),
           Reg2 = factor(Reg2, levels = taglst$Reg2))
tm3 %>% dplyr::count(pDE, Reg1) %>% spread(pDE, n)
tm3 %>% dplyr::count(pDE, Reg2) %>% spread(pDE, n)

tm4 = tm3 %>% 
    mutate(MP = (B73 + Mo17) / 2,
           SPE = ifelse(silent == 1, NA,
                 ifelse(B73>=1 & Mo17<0.1 & pDE=='DE_B', 'SPE_B',
                 ifelse(B73<0.1 & Mo17>=1 & pDE=='DE_M', 'SPE_M', 'non_SPE'))),
           HC = ifelse(is.na(SPE), NA,
                ifelse(SPE=='SPE_B' & (BxM>=1 | BxM >= MP), 'HC_B',
                ifelse(SPE=='SPE_M' & (BxM>=1 | BxM >= MP), 'HC_M', 'non_HC')))) %>%
    select(-MP)
tm4 %>% dplyr::count(Tissue, SPE) %>% spread(SPE, n) %>% print(n=23)
tm = tm4 %>% mutate(Tissue = factor(Tissue, levels = tissues23)) 
if(genome == 'B73') 
    tm = tm %>% filter(gid %in% gids)
#}}}

#{{{ sharing
gids = tm %>% distinct(gid) %>% pull(gid)
ctags1 = c('pDE', 'hDE', 'SPE', 'HC')
ctags2 = c('Dom', 'Reg1', 'Reg2')
ctags = c('silent', ctags1, ctags2)
tags_rm = c(1, 'non_DE', 'PL', 'non_SPE', 'non_HC')

#{{{ num genes of different categories per tissue
tn = tibble()
for (ctag in ctags) {
    tn1 = tm %>% mutate(tag = eval(parse(text=ctag))) %>%
        filter(!is.na(tag)) %>%
        group_by(Tissue, tag) %>% 
        summarise(ngene = n()) %>% ungroup() %>%
        filter(! tag %in% tags_rm) %>%
        mutate(ctag = ctag)
    tn = rbind(tn, tn1)
}
t_num = tn %>% mutate(ctag = factor(ctag, levels = ctags))
#}}}

# among tissue sharing
#{{{ expression
tsh_e = tm %>%
    group_by(gid, silent) %>%
    summarise(n.tis = n()) %>%
    ungroup()
tsh_es = tsh_e %>%
    group_by(gid) %>%
    summarise(n.tis.tot = sum(n.tis)) %>%
    ungroup()
tsh_es %>% count(n.tis.tot)
etags = c('Silent', 'Tissue specific', 'Intermediate frequency', 'Constitutive')
tsh_e = tsh_e %>% filter(silent == 0) %>%
    right_join(tsh_es, by = 'gid') %>%
    replace_na(list(n.tis = 0)) %>%
    mutate(prop.tis = n.tis / n.tis.tot,
           etag = ifelse(prop.tis == 0, etags[1],
                  ifelse(prop.tis <= 0.2, etags[2],
                  ifelse(prop.tis < 0.8, etags[3], etags[4])))) %>%
    mutate(etag = factor(etag, levels = etags)) %>%
    select(gid, n.tis, prop.tis, etag)
#}}} 

#{{{ DE/SPE/HC
mixlst = list(
    pDE = 'DE_B/DE_M',
    hDE = 'BP/AP',
    SPE = 'SPE_B/SPE_M',
    HC = 'HC_B/HC_M'
)
tsTags = c('No data', 'Tissue specific', 'Intermediate frequency', 'Constitutive')
tr = tibble()
for (ctag in ctags1) {
    tr1 = tm %>% mutate(tag = eval(parse(text=ctag))) %>%
        filter(!is.na(tag)) %>%
        group_by(gid, tag) %>%
        summarise(n.tis = n()) %>%
        ungroup()
    trs = tr1 %>%
        group_by(gid) %>%
        summarise(n.tis.tot = sum(n.tis)) %>% ungroup()
    tag_fill = sprintf("non-%s", ctag)
    tr1 = tr1 %>% filter(! tag %in% tags_rm) %>%
        arrange(gid, n.tis) %>%
        group_by(gid) %>%
        summarise(n.tis = sum(n.tis),
                  tag = ifelse(length(unique(tag)) == 1, 
                               sprintf("consis. %s", tag[1]), 
                               sprintf("mix of %s", mixlst[[ctag]]))) %>%
        ungroup() %>%
        right_join(trs, by = 'gid') %>%
        replace_na(list(n.tis = 0, tag = tag_fill)) %>%
        mutate(prop.tis = n.tis / n.tis.tot,
               tsTag = ifelse(prop.tis == 0, tsTags[1],
                       ifelse(prop.tis <= 0.2, tsTags[2],
                       ifelse(prop.tis < 0.8, tsTags[3], tsTags[4]))),
               ctag = ctag) %>%
        select(ctag, gid, tag, tsTag, n.tis, n.tis.tot, prop.tis)
    tr = rbind(tr, tr1)
}
tr = tr %>%
    mutate(ctag = factor(ctag, levels = unique(ctag)),
           tsTag = factor(tsTag, levels = tsTags))
levels(tr$ctag)
tsh_d = tr
#}}}

#{{{ Reg/Reg2
tsTags = c('Tissue specific', 'Intermediate frequency', 'Constitutive')
tr = tibble()
for (ctag in ctags2) {
    tr1 = tm %>% mutate(tag = eval(parse(text=ctag))) %>%
        filter(Tissue %in% tissues20) %>%
        mutate(Tissue = factor(Tissue, levels = tissues20)) %>%
        mutate(tag = as.character(tag)) %>%
        filter(!is.na(tag)) %>%
        group_by(gid, tag) %>%
        summarise(n.tis = n()) %>%
        ungroup()
    trs = tr1 %>%
        group_by(gid) %>%
        summarise(n.tis.tot = sum(n.tis)) %>% ungroup()
    tr1 = tr1 %>%
        right_join(trs, by = 'gid') %>%
        replace_na(list(n.tis = 0)) %>%
        mutate(prop.tis = n.tis / n.tis.tot,
               tsTag = ifelse(prop.tis <= 0.2, tsTags[1],
                       ifelse(prop.tis < 0.8, tsTags[2], tsTags[3])),
               ctag = ctag) %>%
        select(ctag, gid, tag, tsTag, n.tis, n.tis.tot, prop.tis)
    tr = rbind(tr, tr1)
}
tr = tr %>%
    mutate(ctag = factor(ctag, levels = unique(ctag)), 
           tsTag = factor(tsTag, levels = tsTags))
unique(tr$tag)
levels(tr$ctag)
tsh_r = tr
#}}}
#}}}

fo = file.path(dirw, "01.master.rda")
save(tm, taglst, t_num, tsh_e, tsh_d, tsh_r, file = fo)

