source("br.fun.r")
require(DESeq2)
#require(lmtest)
dirw = file.path(dirp, "45.doa")

fi = file.path(dirp, '41.qc/10.byrep.RData')
x = load(fi)
x
tm = t_byrep
th = tm %>% distinct(SampleID, Tissue, Genotype, Treatment)
gts = c("B73", "Mo17", "BxM")
tissues = levels(tm$Tissue)
fi = file.path(dirp, "42_de/10.de.dom.rda")
x = load(fi)
t_de = t_de %>% select(-hDE, -hDE2)

#{{{ assign inheritance pattern using mle+dnbinom+lrtest
tw = tm %>% 
    filter(Genotype %in% gts) %>%
    select(SampleID, Tissue, Genotype, gid, ReadCount) %>%
    left_join(t_sf, by = 'SampleID') %>%
    mutate(ReadCount = round(ReadCount/sizeFactor)) %>%
    select(-sizeFactor)
tt = tw %>% group_by(Tissue, Genotype, gid) %>%
    summarise(nRC = list(ReadCount)) %>%
    ungroup() %>%
    spread(Genotype, nRC) %>%
    inner_join(t_de, by = c("Tissue",'gid')) %>%
    filter(pDE != 'non-DE')
tx = tw %>% group_by(Tissue, Genotype, gid) %>%
    summarise(nRC = mean(ReadCount)) %>%
    ungroup() %>%
    spread(Genotype, nRC) %>%
    inner_join(t_de, by = c("Tissue",'gid')) %>%
    filter(pDE != 'non-DE')
 
get_bic <- function(i, tt) {
    #{{{
    xb = unlist(tt$B73[i]); xm = unlist(tt$Mo17[i]); xh = unlist(tt$BxM[i])
    size = 1/tt$disp[i]
    bica = NA; bicd = NA; bicr = NA; bico = NA; mode = NA
    probb.s = size / (size + mean(xb))
    probm.s = size / (size + mean(xm))
    probh.s = size / (size + mean(xh))
    xb = round(xb); xm = round(xm); xh = round(xh)
    LLa <- function(probb, probm) {
        probh = 2*probb*probm/(probb+probm)
        if(probb > 0 & probb < 1 & probm > 0 & probm < 1)
            -sum(dnbinom(xb, size, prob = probb, log = T)) + 
            -sum(dnbinom(xm, size, prob = probm, log = T)) +
            -sum(dnbinom(xh, size, prob = probh, log = T))
        else 100
    }
    LLd <- function(probb, probm) {
        probh = min(probb, probm)
        if(probb > 0 & probb < 1 & probm > 0 & probm < 1)
            -sum(dnbinom(xb, size, prob = probb, log = T)) + 
            -sum(dnbinom(xm, size, prob = probm, log = T)) +
            -sum(dnbinom(xh, size, prob = probh, log = T))
        else 100
    }
    LLr <- function(probb, probm) {
        probh = max(probb, probm)
        if(probb > 0 & probb < 1 & probm > 0 & probm < 1)
            -sum(dnbinom(xb, size, prob = probb, log = T)) +
            -sum(dnbinom(xm, size, prob = probm, log = T)) +
            -sum(dnbinom(xh, size, prob = probh, log = T))
        else 100
    }
    LLo <- function(probb, probm, probh) {
        if(probb > 0 & probb < 1 & probm > 0 & probm < 1 & probh > 0 & probh < 1)
            -sum(dnbinom(xb, size, prob = probb, log = T)) +
            -sum(dnbinom(xm, size, prob = probm, log = T)) +
            -sum(dnbinom(xh, size, prob = probh, log = T))
        else 100
    }
    nobs = length(xb) + length(xm) + length(xh)
    fita = mle(LLa, start = list(probb = probb.s, probm = probm.s),
          method = "L-BFGS-B", lower = rep(1e-5,2), upper = rep(1-1e-5,3),
          nobs = nobs)
    fitd = mle(LLd, start = list(probb = probb.s, probm = probm.s),
          method = "L-BFGS-B", lower = rep(1e-5,2), upper = rep(1-1e-5,3),
          nobs = nobs)
    fitr = mle(LLr, start = list(probb = probb.s, probm = probm.s),
          method = "L-BFGS-B", lower = rep(1e-5,2), upper = rep(1-1e-5,3),
          nobs = nobs)
    fito = mle(LLo, start = list(probb = probb.s, probm = probm.s, probh = probh.s),
          method = "L-BFGS-B", lower = rep(1e-5,3), upper = rep(1-1e-5,3),
          nobs = nobs)
    #coef(fitc)
    bic = BIC(fita, fitd, fitr, fito)
    bica = bic$BIC[1]; bicd = bic$BIC[2]; bicr = bic$BIC[3]; bico = bic$BIC[4]
    tb = as_tibble(bic) %>%
        add_column(mode = c('add','dom','rec','other')) %>% arrange(BIC)
    c('bica'=bica,'bicd'=bicd,'bicr'=bicr,'bico'=bico,'mode'=tb$mode[1])
    #}}}
}

require(BiocParallel)
bpparam <- MulticoreParam()

y = bplapply(1:nrow(tx), get_bic, tx, BPPARAM = bpparam)
tb = do.call(rbind.data.frame, y) %>% as_tibble()
colnames(tb) = c("bica", "bicd", "bicr", "bico", "mode")
tb = tb %>%
    transmute(bica = as.numeric(bica),
              bicd = as.numeric(bicd),
              bicr = as.numeric(bicr),
              bico = as.numeric(bico))
stopifnot(nrow(tb) == nrow(tx))
table(tb$Dom)
t_dom = tx %>% bind_cols(tb) %>% select(-mode)

fo = file.path(dirw, "09.dom.rda")
save(t_dom, file = fo)
#}}}


