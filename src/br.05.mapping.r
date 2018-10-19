source('br.fun.r')
sid = 'me99b'
#sid = 'me99b.m'
dirw = file.path(dirp, ifelse(sid == 'me99b.m', '41_qc_m', "41_qc"))

#{{{ read mapping stats
fi = file.path(dirp, "03.collect/32.RData")
x = load(fi)

ti = t_mapping
ti = ti %>% 
    mutate(Genotype = ifelse(Genotype == 'B73xMo17', 'BxM', Genotype)) %>%
    mutate(Genotype = ifelse(Genotype == 'Mo17xB73', 'MxB', Genotype)) %>%
    filter(Genotype %in% c("B73", "Mo17", "BxM"))
ti %>% group_by(Genotype) %>%
    summarise(mapping_rate = sum(Pair_Map_Hq)/sum(Pair))
tissues = unique(ti$Tissue)
gts = unique(ti$Genotype)
ti = ti %>%
    mutate(Tissue = factor(Tissue, levels = tissues),
           Genotype = factor(Genotype, levels = gts),
           Treatment = sprintf("Rep%d", Treatment))
tp = ti %>% mutate(Unmap = Pair_Unmap,
                   Orphan = Pair_Orphan,
                   Map_LowQual = Pair_Map - Pair_Map_Hq,
                   Map_HighQual = Pair_Map_Hq) %>%
    select(Tissue, Genotype, Treatment, Unmap, Orphan, Map_LowQual, Map_HighQual) %>%
    gather(type, ReadPairCount, -Tissue, -Genotype, -Treatment) %>%
    group_by(Tissue, Genotype, type) %>%
    summarise(ReadPairCount = sum(ReadPairCount)) %>% ungroup()
tps = tp %>%
    distinct(Tissue, Genotype) %>%
    mutate(lab = sprintf("%s: %s", Tissue, Genotype),
           x = 1:length(Tissue))

types = rev(c("Unmap", "Orphan", "Map_LowQual", "Map_HighQual"))
cols4 = pal_simpsons()(4)[c(1,2,4,3)]
pr = ggplot(tp) +
    geom_bar(aes(x = Genotype, y = ReadPairCount/1000000, fill = type), stat='identity', position = position_stack(reverse = T), width = 0.7) +
    scale_fill_manual(values = cols4) +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_continuous(name='Read Pairs (/million)', expand = c(0,0)) +
    coord_flip() +
    theme_bw() +
    facet_grid(Tissue ~ ., switch = 'y') +
    theme(axis.ticks.y = element_blank()) +
    theme(strip.background = element_blank(), strip.placement = "outside") +
    theme(strip.text.y = element_text(size = 8, angle = 180, hjust=1)) +
    theme(legend.position = c(.5,1), legend.direction = "horizontal", legend.justification = c(.5,0), legend.background = element_blank()) +
    theme(legend.title = element_blank(), legend.key.size = unit(.6, 'lines'), legend.key.width = unit(.6, 'lines'), legend.text = element_text(size = 8)) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank()) +
    theme(plot.margin = unit(c(1.5,.5,.5,.5), "lines")) +
    theme(axis.title.x = element_text(size = 9)) +
    theme(axis.title.y = element_blank()) +
    theme(axis.text.x = element_text(size = 8)) +
    theme(axis.text.y = element_text(size = 8)) 
fr = file.path(dirw, "01.readmapping.pdf")
ggsave(pr, filename=fr, width=6, height=8)

t_sam = t_mapping %>% 
    transmute(SampleID = SampleID,
              Tissue = Tissue, 
              Genotype = Genotype,
              Replicate = Treatment,
              TotalReadPair = ReadPairCount,
              TrimmedReadPair = RetainedReadPairCount,
              MappingRate = sprintf("%.01f%%", Pair_Map/Pair*100),
              UniqueMappingRate = sprintf("%.01f%%", Pair_Map_Hq/Pair*100))

tp = ti %>% mutate(Unmap = Pair_Unmap,
                   Orphan = Pair_Orphan,
                   Map_LowQual = Pair_Map - Pair_Map_Hq,
                   Map_HighQual = Pair_Map_Hq) %>%
    select(Tissue, Genotype, Treatment, Unmap, Orphan, Map_LowQual, Map_HighQual) %>%
    gather(type, ReadPairCount, -Tissue, -Genotype, -Treatment) 
tt1 = tp %>% filter(Treatment != 'Rep4') %>%
    group_by(Tissue, Genotype, Treatment) %>%
    filter(Genotype != 'MxB') %>%
    summarise(ReadPairCount = sum(ReadPairCount)) %>%
    ungroup() %>%
    mutate(Genotype = factor(Genotype, levels = c("B73", "Mo17", "BxM"))) %>%
    arrange(Tissue, Genotype) %>%
    mutate(GenoCond = sprintf("%s_%s", Genotype, Treatment),
           ReadPairCount = ReadPairCount / 1000000)
tts = tt1 %>% distinct(Genotype, Treatment, GenoCond) %>%
    arrange(Genotype, Treatment)
tt = tt1 %>% select(Tissue, GenoCond, ReadPairCount) %>%
    mutate(GenoCond = factor(GenoCond, levels = tts$GenoCond)) %>%
    spread(GenoCond, ReadPairCount)

cnames = c(colnames(tt)[1], tts$Treatment)
ngt = tts %>% count(Genotype)
gtheader = as.integer(c(1, ngt$n))
names(gtheader) = c(' ', as.character(ngt$Genotype))
ft = file.path(dirw, "00.table.RData")
save(t_sam, tt, cnames, ngt, gtheader, file = ft)
#}}}

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
tt = ti %>%
    mutate(Genotype = ifelse(Genotype == 'B73xMo17', 'BxM', Genotype)) %>%
    mutate(Genotype = ifelse(Genotype == 'Mo17xB73', 'MxB', Genotype)) %>%
    filter(Genotype %in% gts) %>%
    mutate(TotalReadPair = total,
           TrimmedReadPair = surviving,
           MappingRate = pair_map/surviving,
           UniqueMappingRate = pair_map_hq/surviving,
           AssignedRate = Assigned/surviving) %>%
    mutate(MappingRate = sprintf("%.1f%%", MappingRate*100), 
           UniqueMappingRate = sprintf("%.1f%%", UniqueMappingRate*100), 
           AssignedRate = sprintf("%.1f%%", AssignedRate*100)) %>%
    select(SampleID, Tissue, Genotype, Replicate, 
           TotalReadPair, TrimmedReadPair,
           MappingRate, UniqueMappingRate)
ft = file.path(dirw, "00.table.rda")
save(tt, file = ft)
#}}}

