require(DESeq2)
require(VGAM)
require(gamlss)
require(BiocParallel)
source("br.fun.r")
bpparam <- MulticoreParam()
#bplog(bpparam) <- TRUE

diri = file.path(dirp, "41_qc")
fm = file.path(diri, '10.rda')
y = load(fm)

dirw = file.path(dird, '44_ase')
diri = '~/projects/maize.expression/data/08_raw_output'
fi = file.path(diri, 'me99b', 'ase.tsv')
ti = read_tsv(fi)

if(1==2) {
#{{{ ASE stats / QC
ti1 = ti %>% inner_join(th[,1:3], by = 'SampleID') %>%
    mutate(ntc = n0 + n1 + ncft, nt = n0 + n1,
           pcft = ncft / ntc, pref = n0 / nt) %>%
    filter(nt >= 10, pcft <= .05)

tps = th %>% mutate(Tissue = factor(Tissue, levels = tissues23)) %>%
    arrange(Tissue, Genotype, Replicate) %>%
    mutate(i = 1:length(SampleID)) 
tpx = tps %>% group_by(Tissue) %>%
    summarise(xmed = mean(i), xmin = min(i), xmax = max(i)) %>% ungroup()
tps = tps %>% mutate(i = factor(i, levels = 1:length(SampleID)))
tp = ti1 %>% inner_join(tps[,c('SampleID','i')], by = 'SampleID')
tps = tp %>% count(SampleID, i, Tissue, Genotype)
p = ggplot(tp) +
    geom_boxplot(aes(x = i, y = pref, color = Genotype), outlier.shape = NA, width = .7) +
    scale_x_discrete(name = 'num. genes', breaks = tps$i, labels = tps$n) +
    scale_y_continuous(name = 'Proportion reads w. B73 allele') +
    facet_wrap(.~Tissue, scale = 'free', ncol = 3) + 
    scale_color_aaas() +
    otheme(xtitle = T, xtext = T, ytitle = T, ytext = T, 
           ygrid = T, xticks = T, yticks = T,
           legend.pos = 'bottom.right') +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 8))
fo = file.path(dirw, '03.prop.ref.pdf')
ggsave(p, file = fo, width = 8, height = 10)

tis = ti1 %>%
    group_by(SampleID, Tissue, Genotype) %>%
    summarise(n.gene = n(),
              pcft.q25 = quantile(pcft, .25),
              pcft.q50 = quantile(pcft, .5),
              pcft.q75 = quantile(pcft, .75),
              pref.q25 = quantile(pref, .25, na.rm = T),
              pref.q50 = quantile(pref, .5, na.rm = T),
              pref.q75 = quantile(pref, .75, na.rm = T))
tis %>% 
    filter(Genotype == 'BxM') %>% 
    select(SampleID, Tissue, n.gene, pcft.q25, pcft.q50, pcft.q75) %>%
    print(n=69)
fo = file.path(dirw, '05.ase.stat.rda')
save(tis, file = fo)
#}}}
}

#{{{ MLE approach to estimate cis/trans
ta = ti %>% inner_join(th[,1:3], by = 'SampleID') %>%
    filter(Genotype == 'BxM') %>%
    mutate(ntc = n0 + n1 + ncft, nt = n0 + n1,
           pcft = ncft / ntc, pref = n0 / nt) %>%
    filter(pcft <= .05, nt >= 10) %>%
    select(SampleID, Tissue, gid, n0, n1, ncft) %>%
    group_by(Tissue, gid) %>%
    summarise(ref = sum(n0), alt = sum(n1), cft = sum(ncft),
           refl = list(n0), altl = list(n1), cftl = list(ncft)) %>%
    ungroup()
ta %>% count(Tissue) %>% print(n=23)

tt = tmm %>%
    select(Tissue, gid, Genotype, nRC) %>%
    spread(Genotype, nRC) %>%
    inner_join(ta, by = c('Tissue','gid')) %>%
    filter(B73+Mo17 >= 20, ref + alt >= 10)
    #mutate(ref = ifelse(ref+alt <= 1000, ref, round(ref*1000/(ref+alt))),
    #       alt = ifelse(ref+alt <= 1000, alt, round(alt*1000/(ref+alt))),
    #       cft = ifelse(ref+alt <= 1000, cft, round(cft*1000/(ref+alt)))) 
tt %>% count(Tissue) %>% print(n=23)

get_bic <- function(i, tt) {
    #{{{
    prob = tt$B73[i] / (tt$B73[i] + tt$Mo17[i])
    x = unlist(tt$refl[i]); size = unlist(tt$refl[i]) + unlist(tt$altl[i])
    prob.h = sum(x) / sum(size)
    bicc = NA; bict = NA; bicct = NA; reg = NA
    if(prob == 0)  
        reg = ifelse(prob.h < 0.05, 'cis',
              ifelse(prob.h > 0.45 & prob.h < 0.55, 'trans', 'cis+trans'))
    else if(prob == 1)
        reg = ifelse(prob.h > 0.95, 'cis',
              ifelse(prob.h > 0.45 & prob.h < 0.55, 'trans', 'cis+trans'))
    else {
        LL <- function(prob, rho) {
            if(prob > 0 & prob < 1 & rho >= 0 & rho < 1)
                -sum(VGAM::dbetabinom(x, size, prob = prob, rho = rho, log = TRUE))
            else 100
        }
        fitc = mle(LL, start = list(prob = prob, rho = 0), 
              method = "L-BFGS-B", lower = c(1e-5,1e-5), upper = c(1-1e-5,1-1e-5), 
              fixed = list(prob = prob),
              nobs = length(x))
        fitt = mle(LL, start = list(prob = prob, rho = 0), 
              method = "L-BFGS-B", lower = c(1e-5,1e-5), upper = c(1-1e-5,1-1e-5), 
              fixed = list(prob = 0.5),
              nobs = length(x))
        fit = mle(LL, start = list(prob = prob, rho = 0), 
              method = "L-BFGS-B", lower = c(1e-5,1e-5), upper = c(1-1e-5,1-1e-5), 
              nobs = length(x))
        #coef(fitc)
        bic = BIC(fitc, fitt, fit)
        bicc = bic$BIC[1]; bict = bic$BIC[2]; bicct = bic$BIC[3]
        tb = as_tibble(bic) %>% 
            add_column(reg = c('cis', 'trans','cis+trans')) %>% arrange(BIC)
        reg = tb$reg[1]
    }
    c('cis'=bicc, 'trans'=bict, 'cis+trans'=bicct, 'reg'=reg)
    #}}}
}

y = bplapply(1:nrow(tt), get_bic, tt, BPPARAM = bpparam)

tb = do.call(rbind.data.frame, y) %>% as_tibble()
colnames(tb) = c("bic.c", "bic.t", "bic.ct", "reg")
tb = tb %>% mutate(bic.c = as.numeric(bic.c), 
                   bic.t = as.numeric(bic.t),
                   bic.ct = as.numeric(bic.ct))
stopifnot(nrow(tb) == nrow(tt))

ttb = cbind(tt, tb) %>%
    select(-refl, -altl, -cftl) %>%
    as_tibble() %>%
    mutate(Tissue = factor(Tissue, levels = tissues23))

reg_levels = c("cis only", "trans only", "cis+trans", "cis B73", "cis Mo17", "trans B73", "trans Mo17")
t_ase = ttb %>%
    mutate(prop.p = B73/(B73+Mo17), prop.h = ref/(ref+alt)) %>%
    mutate(reg = ifelse(reg == 'cis', 'cis only',
                 ifelse(reg == 'trans', 'trans only',
                 ifelse(reg == 'cis+trans',
                 ifelse((prop.h<prop.p & prop.h>=.5)|(prop.h>prop.p & prop.h<=.5), 'cis+trans',
                 ifelse((prop.h<prop.p & prop.h<.5 & prop.p>.5), "trans Mo17", 
                 ifelse((prop.h>prop.p & prop.h>.5 & prop.p<.5), "trans B73",
                 ifelse((prop.h<prop.p & prop.p<.5), "cis Mo17",
                 ifelse((prop.h>prop.p & prop.p>.5), "cis B73", "unc"))))), 
                 'unc')))) %>%
    mutate(reg = factor(reg, levels = reg_levels))
table(t_ase$reg)

fo = sprintf("%s/10.mle.rda", dirw)
save(t_ase, file = fo)
#}}}

if(1==2) {
#{{{ proof of concept for beta-binomial
size = 50
prob = 0.7
x = 0:50
y1 = dbinom(x, size=size, prob=prob)
tp = tibble(model = 'binom', x = x, y = y1)
for (rho in c(0.01,0.05,0.1)) {
    y = dbetabinom(x, size=size, prob=prob, rho=rho)
    tp1 = tibble(model = sprintf("beta-binom: rho=%g", rho), x=x, y=y)
    tp = rbind(tp, tp1)
}

tp$model = factor(tp$model, levels = unique(tp$model))
p <- ggplot(tp) +
    geom_line(aes(x=x,y=y,color=model)) +
    geom_vline(xintercept = 35, size = 0.3) +
    scale_color_brewer(palette = "Set1")
fo = file.path(dirw, 'bb.pdf')
ggsave(p, filename = fo, width = 7, height = 5)
#}}}
}
