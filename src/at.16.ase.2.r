source("br.fun.r")
dirw = file.path(dird, "44_ase")

if(1 == 2) {
#{{{# Indel mapping bias
sids = c("BR006", "BR054", "BR151")
opts = c("snp+indel", "snp")

for (sid in sids) {
for (opt in opts) {
    fi = sprintf("%s/%s.5.%s.bed", diri, sid, opt)
    ti = read.table(fi, header = F, sep = "\t", as.is = T)
    colnames(ti) = c("chr", "beg", "end", "rid", "nref", "nalt", "nunk", "nerr")
    ti = ti[ti$nref+ti$nalt>0, ]

    n0 = nrow(ti)

    n1 = sum(ti$nref + ti$nalt == 1)
    n11 = sum(ti$nref == 1 & ti$nalt == 0)
    n12 = sum(ti$nref == 0 & ti$nalt == 1)

    n2 = sum(ti$nref + ti$nalt >= 2)
    n21 = sum(ti$nref > 1 & ti$nalt == 0)
    n22 = sum(ti$nref == 0 & ti$nalt > 1)
    n23 = sum(ti$nref > 0 & ti$nalt > 0)


    lab = tl$label[tl$SampleID==sid]
    outputs = c(
        '',
        sprintf("%s [%s] against %s", sid, lab, opt),
        sprintf("%8d total reads with >=1 variants:", n0),
        sprintf("%8d: only 1 variants", n1),
        sprintf("   %8d [%5.02f%%] B73 allele", n11, n11/n1*100),
        sprintf("   %8d [%5.02f%%] Mo17 allele", n12, n12/n1*100),
        sprintf("%8d: >=2 variants", n2),
        sprintf("   %8d [%5.02f%%] with all B73 allele", n21, n21/n2*100),
        sprintf("   %8d [%5.02f%%] with all Mo17 allele", n22, n22/n2*100),
        sprintf("   %8d [%5.02f%%] with both B and M alleles", n23, n23/n2*100),
        ''
    )
    cat(paste(outputs, collapse = "\n"))
}
}
#}}}
}

#{{{ process cis/trans categories
fi = file.path(dirw, "../42_de/11.de.dom.rda")
y = load(fi)
td = dd %>% select(Tissue, gid, pDE)
    #mutate(Tissue = factor(Tissue, levels = tissues23)) %>%
fi = file.path(dirw, "10.mle.rda")
x = load(fi)
t_ase = t_ase %>% mutate(Tissue = as.character(Tissue))

regs1 = c("cis only",
          "cis + trans: opp. dir.",
          "cis + trans: same dir.",
          "trans only",
          "unexpected bias")
regs2 = c("conserved", "B73 biased", "Mo17 biased")
ta = t_ase %>%
    inner_join(dd, by = c("Tissue", 'gid')) %>%
    filter(!is.na(pDE)) %>%
    mutate(Reg1 = as.character(reg),
           Reg2 = as.character(reg)) %>%
    mutate(Reg1 = ifelse(pDE == 'non_DE', NA, Reg1)) %>%
    mutate(Reg1 = ifelse(Reg1 %in% c("cis B73", "cis Mo17"), 'cis + trans: same dir.', Reg1)) %>%
    mutate(Reg1 = ifelse(Reg1 %in% c("cis+trans"), 'cis + trans: opp. dir.', Reg1)) %>%
    mutate(Reg1 = ifelse(Reg1 %in% c("trans B73", "trans Mo17"), 'unexpected bias', Reg1)) %>%
    mutate(Reg2 = ifelse(pDE == 'non_DE',
                  ifelse(Reg2 %in% c("cis B73", "trans B73"), "B73 biased",
                  ifelse(Reg2 %in% c("cis Mo17", "trans Mo17"), "Mo17 biased",
                        "conserved")), NA)) %>%
    mutate(Reg1 = factor(Reg1, levels = regs1),
           Reg2 = factor(Reg2, levels = regs2))
ta %>% count(pDE, Reg1) %>% spread(pDE, n)
ta %>% count(pDE, Reg2) %>% spread(pDE, n)

ase = ta %>% select(Tissue, gid, prop.p, prop.h, Reg1, Reg2, bic.c, bic.t, bic.ct)
fo = sprintf("%s/11.ase.rda", dirw)
save(ase, file = fo)
#}}}

