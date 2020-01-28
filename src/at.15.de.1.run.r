source("br.fun.r")
sid = 'me99b'
#sid = 'me99b.m'
dirw = file.path(dirp, ifelse(sid == 'me99b.m', '42_de_m', "42_de"))
genome = ifelse(sid == 'me99b.m', 'Mo17', 'B73')
x = load(file.path(dirg, genome, '55.rda'))
diri = file.path(dirp, ifelse(sid == 'me99b.m', '41_qc_m', "41_qc"))
fm = file.path(diri, '10.rda')
y = load(fm)
tm = tm %>% inner_join(th[,1:2], by = 'SampleID') %>%
    group_by(Tissue) %>% nest()
th = th %>% inner_join(tl, by = 'SampleID') %>%
    select(-paired, -Treatment) %>%
    group_by(Tissue) %>% nest()
#dirw = file.path(dird, "42_de_old")

#{{{ run DESeq2 and edgeR
#res = run_de_test(tm$data[[1]], th$data[[1]])
t_de = tm %>% inner_join(th, by = 'Tissue') %>%
    mutate(res = map2(data.x, data.y, run_de_test)) %>%
    mutate(deseq = map(res, 'deseq'), edger = map(res, 'edger')) %>%
    select(Tissue, deseq, edger)
fo = file.path(dirw, '10.rda')
save(t_de, file = fo)
#}}}

fi = file.path(dirw, '10.rda')
x = load(fi)

#{{{ call DE
ta = t_de %>%
    mutate(data = map2(deseq, edger, .f = function(d1, d2)
        tibble(gid=d1$gid,
               p.ds = d1$padj.mb, lfc.ds = d1$log2mb,
               p.eg = d2$padj.mb, lfc.eg = d2$log2mb))) %>%
    select(Tissue, data) %>%
    unnest()
    #mutate(logp.ds=-log(p.ds), logp.eg=-log(p.eg)) %>%

pvals = c('**','*','n.s.')
tp0 = ta %>%
    mutate(deseq.p = ifelse(p.ds<.01, '**', ifelse(p.ds<.05, '*', 'n.s.')),
           edger.p = ifelse(p.eg<.01, '**', ifelse(p.ds<.05, '*', 'n.s.')),
           deseq.lfc = ifelse(abs(lfc.ds) >= 1, 'logFC>1', 'logFC<1'),
           edger.lfc = ifelse(abs(lfc.eg) >= 1, 'logFC>1', 'logFC<1'))
tp1 = tp0 %>% select(Tissue, pval = deseq.p, lfc = deseq.lfc) %>%
    mutate(tool = 'deseq')
tp2 = tp0 %>% select(Tissue, pval = edger.p, lfc = edger.lfc) %>%
    mutate(tool = 'edger')
tp = tp1 %>% bind_rows(tp2) %>% dplyr::count(Tissue, tool, pval, lfc) %>%
    mutate(pval = factor(pval, levels = !!pvals))
p = ggplot(tp) +
    geom_bar(aes(x = pval, y = n, fill = lfc), stat='identity',
             position=position_stack(reverse=T),width=.8) + 
    facet_wrap(Tissue~tool, ncol = 10) +
    scale_fill_d3() +#labels = c('logFC<1','logFC>1')) +
    otheme(xtext=T, ytext=T) +
    theme(legend.position = 'top')
fo = file.path(dirw, 'test2.pdf')
ggsave(p, file = fo, width = 8, height = 8)
#}}}

#{{{ compare DESeq2 and edgeR
tags = c("DESeq2+edgeR", "edgeR", "DESeq2", "non-DE")
tp = ta %>%
    mutate(tag.ds = p.ds<.01, tag.eg = p.eg<.01) %>%
    dplyr::count(Tissue, tag.ds, tag.eg) %>%
    mutate(tag = ifelse(tag.ds, ifelse(tag.eg, 'DESeq2+edgeR', 'DESeq2'),
                        ifelse(tag.eg, 'edgeR', 'non-DE'))) %>%
    mutate(tag = factor(tag, levels = tags))
#
p = ggplot(tp) +
    geom_bar(aes(x = Tissue, y = n, fill = tag), stat='identity',
             position=position_stack(reverse=T),width=.8) +
    coord_flip() +
    scale_fill_npg() +
    otheme(xtext=T, ytext=T) +
    theme(legend.position = 'top') +
    theme_classic()
fo = file.path(dirw, 'test.pdf')
ggsave(p, file = fo, width = 6, height = 6)
#}}}

stopifnot(1 > 2)
#{{{ # edgeR ase - not working
run_edgeR_ase <- function(tc, tw) {
    #{{{
    gids = tw$gid
    twd = data.frame(tw[,-1])
    rownames(twd) = gids
    y = DGEList(counts = twd, 
                lib.size = tc$lib.size, 
                norm.factors = tc$norm.factors)
    design = model.matrix(~0 + parent + parent:allele, data = tc)
    #colnames(design) = levels(tc$Genotype)
    #y = estimateDisp(y, design)
    y = estimateGLMCommonDisp(y, design, verbose = T)
    y = estimateGLMTrendedDisp(y, design)
    y = estimateGLMTagwiseDisp(y, design)
    fit = glmFit(y, design)
    
    lrt = glmLRT(fit, contrast = c(0, 0, 1, -1))
    stopifnot(identical(gids, rownames(lrt$table)))
    tags = decideTestsDGE(lrt, adjust.method = "BH", p.value = .05)
    summary(tags)
    tr1 = lrt$table %>% as_tibble() %>%
        bind_cols(tw) %>% arrange(PValue) %>% print(n = 10, width = 2000)
    
    lrt = glmLRT(fit, contrast = c(0, 0, -1, 0))
    stopifnot(identical(gids, rownames(lrt$table)))
    tags = decideTestsDGE(lrt, adjust.method = "BH", p.value = .05)
    summary(tags)
    tr1 = lrt$table %>% as_tibble() %>%
        bind_cols(tw) %>% arrange(PValue) %>% print(n = 10, width = 2000)
    #}}}
}

    ta1 = tm1 %>% filter(Genotype %in% c("B73", "Mo17")) %>%
        mutate(sid = sprintf("%s.%s.Parent", SampleID, Genotype)) %>%
        select(sid, gid, ReadCount) %>%
        spread(sid, ReadCount)
    ta2 = tm1 %>% filter(Genotype == 'BxM') %>%
        mutate(B73 = nref, Mo17 = nalt) %>%
        select(SampleID, gid, B73, Mo17) %>%
        gather(gt, ReadCount, -gid, -SampleID) %>%
        mutate(sid = sprintf("%s.%s.Hybrid", SampleID, gt)) %>%
        select(sid, gid, ReadCount) %>%
        spread(sid, ReadCount)
    ta = ta1 %>% inner_join(ta2, by = 'gid')
    res = run_edgeR(th1, tw)
    tr = res$tr; t_nf = res$t_nf; t_cpm = res$t_cpm
    tr5 = run_edgeR_dom(th1, tr, t_cpm)
    tr = tr %>% left_join(tr5, by = 'gid')
    tah = tibble(sid = colnames(ta)[-1]) %>%
        separate(sid, c("SampleID", "allele", "parent"), sep = "[.]") %>%
        left_join(t_nf, by = 'SampleID')
    #tt = run_edgeR_ase(tah, ta)
    #tr2 = run_deseq2(tc, tw)
    #tx = tr %>% inner_join(tr2, by = c('Tissue', 'gid')) 
    #table(tx[,c("detag", "tag.mb")])
    t_de = rbind(t_de, tr)

doms = c("BLP", "LP", "PD_L", "MP", "PD_H", "HP", "AHP")
td = tr %>%
    filter(tag.mb != 0) %>%
    mutate(tag.lp = ifelse(cpm.b < cpm.m, tag.hb, tag.hm),
           tag.hp = ifelse(cpm.b < cpm.m, tag.hm, tag.hb)) %>%
    mutate(Dom = "MP") %>%
    mutate(Dom = ifelse(tag.nonadd == -1 & tag.lp == -1, 'BLP', Dom)) %>%
    mutate(Dom = ifelse(tag.nonadd == -1 & tag.lp == 0, 'LP', Dom)) %>%
    mutate(Dom = ifelse(tag.nonadd == -1 & tag.lp == 1, 'PD_L', Dom)) %>%
    mutate(Dom = ifelse(tag.nonadd == 1 & tag.hp == -1, 'PD_H', Dom)) %>%
    mutate(Dom = ifelse(tag.nonadd == 1 & tag.hp == 0, 'HP', Dom)) %>%
    mutate(Dom = ifelse(tag.nonadd == 1 & tag.hp == 1, 'AHP', Dom)) %>%
    mutate(Dom = factor(Dom, levels = doms)) %>% 
    mutate(cpm.mp = (cpm.b + cpm.m)/2,
           doa = (cpm.h - cpm.mp)/(pmax(cpm.b, cpm.m) - cpm.mp)) %>%
    mutate(doa = ifelse(doa < -3, -3, doa)) %>%
    mutate(doa = ifelse(doa > 3, 3, doa)) %>%
    transmute(gid = gid, pDE = tag.mb, Dom = Dom, DoA = doa)
table(td$Dom)
#}}}

