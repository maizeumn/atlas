source('br.fun.r')
sid = 'me99b'
#sid = 'me99b.m'
dirw = file.path(dirp, ifelse(sid == 'me99b.m', '41_qc_m', "41_qc"))

diri = '~/projects/maize.expression/data/11_qc'
fi = file.path(diri, sid, '20.rc.norm.rda')
fi = file.path(diri, sid, '10.mapping.stat.tsv')
ti = read_tsv(fi)

types = c("FailedQC", "Unmap", "Map_LowQual", "Unassigned", "Assigned", 'Total')
tx = ti %>%
    mutate(Genotype = ifelse(Genotype == 'B73xMo17', 'BxM', Genotype)) %>%
    mutate(Genotype = ifelse(Genotype == 'Mo17xB73', 'MxB', Genotype)) %>%
    filter(Genotype %in% gts) %>%
    mutate(Total = total,
           FailedQC = dropped,
           Unmap = pair_unmap + unpair_unmap,
           Map_LowQual = pair_map-pair_map_hq + pair_orphan-pair_orphan_hq + unpair_map-unpair_map_hq,
           Unassigned = pair_map_hq+pair_orphan_hq+unpair_map_hq - Assigned,
           Assigned = Assigned) %>%
    select(Tissue, Genotype, Treatment, FailedQC, Unmap, Map_LowQual, Unassigned, Assigned, Total) %>%
    gather(type, ReadPairCount, -Tissue, -Genotype, -Treatment) %>%
    group_by(Tissue, Genotype, type) %>%
    summarise(ReadPairCount = sum(ReadPairCount)) %>% ungroup() %>%
    mutate(Tissue = factor(Tissue, levels = tissues23)) %>%
    mutate(Genotype = factor(Genotype, levels = gts)) %>%
    mutate(type = factor(type, levels = rev(types)))

#{{{ read mapping plot
tps = tx %>% filter(type == 'Total') %>% select(-type) %>%
    rename(total=ReadPairCount)
tp = tx %>% filter(type != 'Total') %>%
    inner_join(tps, by = c('Tissue','Genotype')) %>%
    mutate(prop = ReadPairCount/total) %>%
    mutate(lab = sprintf("%.02f", ReadPairCount/total)) %>%
    mutate(lab = str_replace(lab, "^0[\\.]", "."))
#
cols5 = pal_simpsons()(5)[c(1,2,4,5,3)]
pr = ggplot(tp) +
    geom_bar(aes(x = Genotype, y = ReadPairCount/1000000, fill = type), stat='identity', position = position_stack(reverse = T), width = .7) +
    scale_fill_manual(values = cols5) +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_continuous(name='Read Pairs (/million)', expand = expand_scale(mult=c(0,.02))) +
    coord_flip() +
    facet_grid(Tissue ~ ., switch = 'y') +
    otheme(legend.pos = 'none', xgrid = T, ygrid = F, 
           xtitle = T, ytitle = F, xtext = T, ytext = T) +
    theme(strip.background = element_blank(), strip.placement = "outside") +
    theme(strip.text.y = element_text(size = 8, angle = 180, hjust=1)) +
    theme(legend.position = c(.5,1), legend.direction = "horizontal", legend.justification = c(.5,0), legend.background = element_blank()) +
    theme(legend.title = element_blank(), legend.key.size = unit(.6, 'lines'), legend.key.width = unit(.6, 'lines'), legend.text = element_text(size = 8)) +
    theme(axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8)) +
    theme(plot.margin = unit(c(1.5,.5,.5,.5), "lines"))
fr = file.path(dirw, "01.readmapping.pdf")
ggsave(pr, filename=fr, width=6, height=8)
#}}}

#{{{ mapping prop. plot
tp = tx %>% inner_join(tps, by = c('Tissue', 'Genotype')) %>%
    mutate(lab = ifelse(type == 'Total', sprintf("%.01fM", total/1000000),
                        sprintf("%.0f%%", ReadPairCount/total*100))) %>%
    mutate(lab = str_replace(lab, "^0[\\.]", "."))
#
cols6 = pal_simpsons()(6)[c(6,1,2,4,5,3)]
pr = ggplot(tp) +
    geom_tile(aes(x = type, y = 1, fill = type), alpha = .5) +
    geom_text(aes(x = type, y = 1, label = lab), size = 2.8) +
    scale_fill_manual(values = cols6) +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    facet_grid(Tissue ~ Genotype, switch = 'y') +
    otheme(legend.pos = 'none', xgrid = F, ygrid = F, 
           xtitle = F, ytitle = F, xtext = F, ytext = F) +
    theme(panel.border = element_blank(), panel.spacing = unit(.1,'lines')) +
    theme(strip.background = element_blank(), strip.placement = "outside") +
    theme(strip.text.y = element_text(size = 8, angle = 180, hjust=1)) +
    theme(legend.position = c(.5,1), legend.direction = "horizontal", legend.justification = c(.5,-.6), legend.background = element_blank()) +
    theme(legend.title = element_blank(), legend.key.size = unit(.6, 'lines'), legend.key.width = unit(.6, 'lines'), legend.text = element_text(size = 8)) +
    guides(fill = guide_legend(nrow = 1)) +
    theme(plot.margin = unit(c(1.5,.5,.5,.5), "lines"))
fr = file.path(dirw, "02.readmapping.prop.pdf")
ggsave(pr, filename=fr, width=8, height=8)
#}} 
#}}}

#{{{ save to 00.table.rda
tissues6 = c('coleoptile_tip', 'radicle_root', 'embryo_imbibedseed', 'seedlingleaf_11DAS', 'seedlingroot_11DAS', 'seedlingmeristem_11DAS')
tt = ti %>%
    mutate(Genotype = ifelse(Genotype == 'B73xMo17', 'BxM', Genotype)) %>%
    mutate(Genotype = ifelse(Genotype == 'Mo17xB73', 'MxB', Genotype)) %>%
    filter(Genotype %in% gts) %>%
    mutate(Condition = ifelse(Tissue %in% tissues6, 'Growth chamber', 'Field')) %>%
    mutate(TotalReadPair = total,
           TrimmedReadPair = surviving,
           MappingRate = pair_map/surviving,
           UniqueMappingRate = pair_map_hq/surviving,
           AssignedRate = Assigned/surviving) %>%
    mutate(MappingRate = sprintf("%.1f%%", MappingRate*100),
           UniqueMappingRate = sprintf("%.1f%%", UniqueMappingRate*100),
           AssignedRate = sprintf("%.1f%%", AssignedRate*100)) %>%
    select(SampleID, Tissue, Genotype, Replicate, Condition,
           TotalReadPair, TrimmedReadPair,
           MappingRate, UniqueMappingRate, AssignedRate)
ft = file.path(dirw, "00.table.rda")
save(tt, file = ft)
#}}}

