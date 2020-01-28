source("br.fun.r")
require(DESeq2)
#require(lmtest)
require(edgeR)
dirw = file.path(dird, "42_de_old")

fi = file.path(dird, '41_qc/10.rda')
x = load(fi)
tm = tm %>% inner_join(th[,1:2], by = 'SampleID') %>%
    group_by(Tissue) %>% nest()
th = th %>% inner_join(tl, by = 'SampleID') %>%
    select(-paired, -Treatment) %>%
    group_by(Tissue) %>% nest()

#{{{ # proof of concept for NB
mu = 20
x = 0:(mu*4)
y1 = dpois(x, lambda = mu)
tp = tibble(model = 'poisson', x = x, y = y1)
for (prob in c(0.1,0.2,0.5,0.7,0.9)) {
    size = mu * prob / (1-prob)
    y = dnbinom(x, size=size, prob=prob)
    tp1 = tibble(model = sprintf("nbinom: size=%.1f, prob=%g", size, prob), x=x, y=y)
    tp = rbind(tp, tp1)
}
tp$model = factor(tp$model, levels = unique(tp$model))
p <- ggplot(tp) +
    geom_line(aes(x=x,y=y,color=model)) +
    geom_vline(xintercept = 20, size = 0.3) +
    theme(legend.position = 'top', legend.direction = 'vertical') +
    scale_color_brewer(palette = "Set1")
fo = file.path(dirw, 'nb.pdf')
ggsave(p, filename = fo, width = 6, height = 6)
#}}}

tissue = "ear_v14"
#{{{ identify DEGs using DESeq2
tc = th %>% filter(Tissue == tissue, Genotype %in% gts) %>% arrange(SampleID)
tcd = data.frame(tc)
rownames(tcd) = tc$SampleID
tw = tm %>% 
    filter(Tissue == tissue, Genotype %in% gts) %>%
    select(SampleID, gid, ReadCount) %>%
    spread(SampleID, ReadCount) %>%
    mutate(totalRC = rowSums(.[grep("BR", names(.))], na.rm = T)) %>%
    filter(totalRC >= 10) %>%
    select(-totalRC)
twd = data.frame(tw[,-1])
rownames(twd) = tw$gid
stopifnot(identical(tcd$SampleID, colnames(twd)))

dds = DESeqDataSetFromMatrix(countData=twd, colData=tcd, design = ~ Genotype)
#dds = estimateSizeFactors(dds1)
dds = DESeq(dds, fitType = 'parametric')
disp = dispersions(dds)

res = results(dds, contrast = c("Genotype", "Mo17", "B73"), pAdjustMethod = "fdr")
td = as_tibble(data.frame(res)) %>%
    add_column(disp = disp) %>%
    mutate(tissue = tissue,
           gid = rownames(res),
           log2MB = log2FoldChange,
           DE_B = padj < .05 & log2MB < -1,
           DE_M = padj < .05 & log2MB > 1) %>%
    select(tissue, gid, disp, DE_B, DE_M, log2MB, pvalue, padj) %>%
    replace_na(list(DE_B = F, DE_M = F)) %>%
    mutate(DE = ifelse(padj<.05, "DE", "non-DE"))
#}}}

#{{{ identify DEGs using mle+dnbinom+lrtest
tx = as_tibble(counts(dds, normalized = T)) %>%
    mutate(gid = tw$gid) %>%
    gather(SampleID, nRC, -gid) %>%
    left_join(tc[,c("SampleID","Genotype")], by = 'SampleID') %>%
    group_by(gid, Genotype) %>%
    summarise(nRC = list(nRC)) %>%
    ungroup() %>%
    spread(Genotype, nRC) %>%
    mutate(disp = disp)
tx2 = as_tibble(counts(dds, normalized = T)) %>%
    mutate(gid = tw$gid) %>%
    gather(SampleID, nRC, -gid) %>%
    left_join(tc[,c("SampleID","Genotype")], by = 'SampleID') %>%
    group_by(gid, Genotype) %>%
    summarise(nRC = mean(nRC)) %>%
    ungroup() %>%
    spread(Genotype, nRC) %>%
    mutate(disp = disp)

get_bic <- function(i, tt) {
    #{{{
    xb = unlist(tt$B73[i]); xm = unlist(tt$Mo17[i])
    size = 1/tt$disp[i]
    bicd = NA; bicn = NA; mode = NA
    probb.s = size / (size + mean(xb))
    probm.s = size / (size + mean(xm))
    xb = round(xb); xm = round(xm)
        LLd <- function(probb, probm) {
            if(probb > 0 & probb < 1 & probm > 0 & probm < 1)
                -sum(dnbinom(xb, size, prob = probb, log = T) + 
                     dnbinom(xm, size, prob = probm, log = T))
            else 100
        }
        LLn <- function(probb) {
            if(probb > 0 & probb < 1)
                -sum(dnbinom(xb, size, prob = probb, log = T) + 
                     dnbinom(xm, size, prob = probb, log = T))
            else 100
        }
        fitd = mle(LLd, start = list(probb = probb.s, probm = probm.s), 
              method = "L-BFGS-B", lower = c(1e-5,1e-5), upper = c(1-1e-5,1-1e-5), 
              nobs = length(xb)+length(xm))
        fitn = mle(LLn, start = list(probb = probb.s), 
              method = "L-BFGS-B", lower = c(1e-5), upper = c(1-1e-5), 
              nobs = length(xb)+length(xm))
        #coef(fitc)
        lrt = lrtest(fitd, fitn)
        pval = lrt[2,5]
        bic = BIC(fitd, fitn)
        bicd = bic$BIC[1]; bicn = bic$BIC[2]
        tb = as_tibble(bic) %>% 
            add_column(mode = c('DE', 'non-DE')) %>% arrange(BIC)
    c('bicd'=bicd, 'bicn'=bicn, 'mode'=tb$mode[1], 'pval'=pval)
    #}}}
}

require(BiocParallel)
bpparam <- MulticoreParam()

y = bplapply(1:nrow(tx), get_bic, tx, BPPARAM = bpparam)

tb = do.call(rbind.data.frame, y) %>% as_tibble()
colnames(tb) = c("bicd", "bicn", "mode", 'pval')
tb = tb %>% mutate(bicd = as.numeric(bicd), 
                   bicn = as.numeric(bicn),
                   pval = as.numeric(pval),
                   padj = p.adjust(pval, "fdr"),
                   DE = ifelse(mode=='DE' & padj < 0.05, "DE", "non-DE"))
stopifnot(nrow(tb) == nrow(tx))

table(tibble(DE_DESeq2 = td$DE, DE_DIY = tb$DE))
#}}}

#{{{ use GLM to identify additive/dominant pattern
t_de = tm %>% inner_join(th, by = 'Tissue') %>%
    mutate(res = map2(data.x, data.y, run_de_test)) %>%
    mutate(deseq = map(res, 'deseq'), edger = map(res, 'edger')) %>%
    select(Tissue, deseq, edger)
fo = file.path(dirw, '10.rda')
save(t_de, file = fo)
#}}}

#{{{
fi = file.path(dirw, '10.rda')
x = load(fi)

tx = t_de %>% select(Tissue, deseq) %>% unnest() %>%
    filter(padj.mb < .01) %>%
    mutate(padj.hh = ifelse(log2mb<0, padj.hb, padj.hm),
           log2hh = ifelse(log2mb<0, log2hb, log2hm),
           padj.hl = ifelse(log2mb<0, padj.hm, padj.hb),
           log2hl = ifelse(log2mb<0, padj.hm, padj.hb)) %>%
    mutate(dom = ifelse(padj.fm<.01, 
        ifelse(log2hm<0,
               ifelse(padj.hl<.01, ifelse(log2hl<0, 'BLP','PD_L'), 'LP'),
               ifelse(padj.hh<.01, ifelse(log2hh>0, 'AHP','PD_H'), 'HP')),
        'MP'))
doms= c("BLP","LP","PD_L","MP","PD_H","HP","AHP")
tp = tx %>% dplyr::count(Tissue, dom) %>%
    mutate(Tissue=factor(Tissue, levels = tissues23)) %>%
    mutate(dom = factor(dom, levels = doms))

tpx = tp %>% group_by(Tissue) %>% summarise(n = sum(n)) %>% ungroup() %>%
    mutate(lab = sprintf("%s (%d)", Tissue, n))
p = ggplot(tp, aes(Tissue, n, fill=dom)) +
    geom_bar(stat='identity',position='stack') +
    scale_x_discrete(breaks = tpx$Tissue, labels = tpx$lab, expand=c(0,0)) +
    scale_fill_simpsons() +
    coord_flip() +
    otheme(xtext=T,ytext=T)
fo = file.path(dirw, 'test2.pdf')
ggsave(p, file=fo, width=6, height=6)

tp %>% group_by(Tissue, dom) %>% 
    summarise(n = n(), q5=quantile(log2hm,.05), q25=quantile(log2hm,.25),
              q50=quantile(log2hm,.5), q75=quantile(log2hm,.75),
              q95=quantile(log2hm,.95))

p = ggplot(tp) +
    #geom_point(aes(x=asinh(mp), y=asinh(BxM), color=log(padj))) +
    geom_density(aes(x = log2HM, fill = dom), alpha=.5) +
    scale_x_continuous(limits = c(-3,3)) +
    scale_fill_aaas() +
    scale_color_viridis()
ggsave(p, file=file.path(dirw,'test.pdf'), width=8,height=8)
#}}}

#{{{ glm test
y = tm %>% filter(Tissue == tissue) %>% select(-Tissue) %>% unnest() %>%
    inner_join(vh[,c('SampleID','Genotype','ac','hybrid')], by = 'SampleID') 

y0 = y %>% filter(gid == vm$gid[1])
fit0 = glm.nb(ReadCount ~ 1, data=y0)
fit1 = glm.nb(ReadCount ~ Genotype, data=y0)
fit2 = glm.nb(ReadCount ~ ac, data=y0)

anova(fit0, fit1, fit2, test = 'Chisq')
#}}}

