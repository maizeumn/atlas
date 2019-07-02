source("functions.R")
dirw = file.path(dird, "01_exp_design")

#{{{ 23 dev atlas * 3 genotypes
fi = file.path(dirw, "01.tsv")
ti = read_tsv(fi)
ti

#{{{ create reverse mapping for grinding tubebox
boxnames = c("BR_P01_Boxx", "BR_P02_Box3", "BR_P03_Box6", "BR_R01_Box1", "BR_R02_Box2", "BR_R03_Box4", "BR_R04_Box5", "BR_R05_Box7")

to = data.frame()
for (boxname in boxnames) {
    fi = sprintf("%s/11_tubebox_grinding/%s.tsv", dirw, boxname)
    ti = read.table(fi, sep="\t", header = F, as.is = T)
    ti = cbind(rownum = rownames(ti), ti)
    tos = reshape(ti, direction = 'long', idvar = 'rownum', varying = colnames(ti)[2:ncol(ti)], timevar = "colnum", v.names = "sid")
    ps = strsplit(boxname, split = "_")[[1]]
    tos = cbind(bname = ps[2], bname2 = ps[3], tos)
    to = rbind(to, tos)
}
to = to[to$sid != '' & !is.na(to$sid),]
lmap = LETTERS[1:9]
to = cbind(to, coord = sprintf("%s_%s_%s%d", to$bname, to$bname2, lmap[as.numeric(to$rownum)], as.numeric(to$colnum)))

t1 = to[to$bname %in% c("P01", "P02", "P03"),]
t1$sid = as.integer(t1$sid)
stopifnot(identical(1:nsam, sort(t1$sid)))

t2 = to[!to$bname %in% c("P01", "P02", "P03"),]
lens = sapply(t2$sid, nchar)
sid = substr(t2$sid, 1, lens-1)
tag = substr(t2$sid, lens, lens)
t2$sid = as.integer(sid)
t21 = t2[tag == 'A',]
t22 = t2[tag == 'B',]
stopifnot(identical(1:nsam, sort(t21$sid)))
stopifnot(identical(1:nsam, sort(t22$sid)))

tr = merge(t01, t1[,c('sid','coord')], by.x = 'Sample', by.y = 'sid')
colnames(tr)[ncol(tr)] = "Storage.Protein"
tr = merge(tr, t21[,c('sid','coord')], by.x = 'Sample', by.y = 'sid')
colnames(tr)[ncol(tr)] = "Storage.RNA.A"
tr = merge(tr, t22[,c('sid','coord')], by.x = 'Sample', by.y = 'sid')
colnames(tr)[ncol(tr)] = "Storage.RNA.B"

tr = tr[match(tr$Sample, t01$Sample),]
fo = file.path(dirw, "12.tubebox.grinding.tsv")
#write.table(tr[,-c(2:5)], fo, sep = "\t", row.names = F, col.names = F, quote = F)
#}}}

#{{{ multiplex design
f33 = file.path(dirw, "33.multiplex.done.tsv")
t33 = read.table(f33, sep = "\t", as.is = T, header = F)
t33 = t33[!is.na(t33$V2),]
t33$V2 = sprintf("SetB_Index%02d", t33$V2)

seta = c(2,4:7,12:16,18,19)
setb = c(1,3,8:11,20:23,25,27)
namesa = sprintf("SetA_Index%02d", seta)
namesb = sprintf("SetB_Index%02d", setb)
idxnames = c(rbind(namesa, namesb))

poolmap = c(
	"blade_v12" = 1,
	"auricle_v12" = 1,
	"sheath_v12" = 2,
	"internode_v12" = 2,
	"tassel_v12" = 3,
	"ear_v14" = 3,
	"silk_0DAP" = 4,
	"spikelets_0DAP" = 4,
	"husk_0DAP" = 5,
	"tasselstem_0DAP" = 5,
	"floret_0DAP" = 6,
	"flagleaf_0DAP" = 7,
	"root_0DAP" = 7,
	"kernel_14DAP" = 6,
	"endosperm_14DAP" = 8,
	"embryo_27DAP" = 8,
	"endosperm_27DAP" = 9
)
tr = cbind(t01, pool = NA, idxname = NA)
for (tissue in unique(tr$Tissue)) {
	stopifnot(tissue %in% names(poolmap))
	idxs = which(tr$Tissue == tissue)
	tr$pool[idxs] = sprintf("Pool%d", poolmap[tissue])
}

for (i in 1:nrow(t33)) {
	tr$idxname[tr$Sample == t33$V1[i]] = t33$V2[i]
}
for (pool in unique(tr$pool)) {
	if(pool != "Pool7") {
		idxs = which(tr$pool == pool)
		tr$idxname[idxs] = idxnames[1:length(idxs)]
	} else {
		idxs = which(tr$pool == pool & is.na(tr$idxname))
		idxnames_used = tr$idxname[tr$pool == pool & !is.na(tr$idxname)]
		idxnames_unused = idxnames[! idxnames %in% idxnames_used]
		tr$idxname[idxs] = idxnames_unused[1:length(idxs)]
	}
}
ddply(tr, .(pool), summarise, nidx = length(unique(idxname)))
tr = tr[match(tr$Sample, t01$Sample),]
fo = file.path(dirw, "39.pool.tsv")
write_tsv(tr[,-c(2:5)], fo)
#}}}

#{{{ multiplex for additional 54 samples
seta = c(2,4:7,12:16,18,19)
setb = c(1,3,8:11,20:23,25,27)
namesa = sprintf("SetA_Index%02d", seta)
namesb = sprintf("SetB_Index%02d", setb)
idxnames = c(rbind(namesa, namesb))

poolmap = c(
	"coleoptile_tip" = 10,
	"radicle_root" = 10,
	"embryo_imbibedseed" = 11,
	"seedlingleaf_11DAS" = 11,
	"seedlingroot_11DAS" = 12,
	"seedlingmeristem_11DAS" = 12
)
tr = cbind(t01[t01$SampleID > "BR165",], pool = NA, idxname = NA)
for (tissue in unique(tr$Tissue)) {
	stopifnot(tissue %in% names(poolmap))
	idxs = which(tr$Tissue == tissue)
	tr$pool[idxs] = sprintf("Pool%d", poolmap[tissue])
}
for (pool in unique(tr$pool)) {
	idxs = which(tr$pool == pool)
	tr$idxname[idxs] = idxnames[1:length(idxs)]
}
ddply(tr, .(pool), summarise, nidx = length(unique(idxname)))
fo = file.path(dirw, "39.pool2.tsv")
#write_tsv(tr, fo)
#}}}
#}}}

#{{{ hybrid panel sheet
ti = crossing(Genotype=gts10, Tissue=tissues3, Rep=1:4) %>%
    mutate(Genotype=factor(Genotype, levels=gts10)) %>%
    mutate(Tissue=factor(Tissue, levels=tissues3)) %>%
    arrange(Tissue, Genotype, Rep) %>%
    mutate(SampleID = sprintf("BR%03d", 500+(1:length(Tissue)))) %>%
    select(SampleID, everything())
fo = file.path(dirw, '61.meta.tsv')
write_tsv(ti, fo)
#}}}



