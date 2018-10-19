require(dplyr)
require(plyr)
require(GenomicRanges)
require(randomForest)
source("/home/springer/zhoux379/source/genie3/GENIE3_R/genie3.R")

dirw = file.path(Sys.getenv("misc2"), "grn23", "62.genie3")

## run genie3
if (FALSE) {
require(randomForest)
source("/home/springer/zhoux379/source/genie3/GENIE3_R/genie3.R")

dirw = file.path(Sys.getenv("misc2"), "grn23", "62.genie3")

fi = file.path(dirw, "01.matrix.tsv")
fr = '/home/springer/zhoux379/data/genome/Zmays_v4/TF/11.TF.txt'
rids = scan(fr, what = character())

expr.matrix <- read.expr.matrix(fi, form="rows.are.samples")
gids = rownames(expr.matrix)
rids = rids[rids %in% gids]

weight.matrix <- get.weight.matrix(expr.matrix, input.idx=rids)
link.list <- get.link.list(weight.matrix, report.max=10000)

save(weight.matrix, link.list, file = file.path(dirw, "10.genie3.rda"))
}

##
fi = file.path(dirw, "../37.rpkm.filtered.tsv")
ti = read.table(fi, header = T, sep = "\t", as.is = T)
gids = ti$gid

to = ti
rownames(to) = ti$gid
to = to[,-1]
fo = file.path(dirw, "01.matrix.tsv")
write.table(t(to), fo, sep = "\t", row.names = F, col.names = T, quote = F)

fr = '/home/springer/zhoux379/data/genome/Zmays_v4/TF/11.TF.txt'
rids = scan(fr, what = character())
rids = rids[rids %in% gids]
write(rids, file.path(dirw, "11.TF.txt"))

##
fx = file.path(dirw, "10.genie3.rda")
x = load(fx)
x

tl = link.list
colnames(tl) = c("rid", "tid", "im")
tl = tl[order(tl$im, decreasing = T),]

## read aracne links
fa = file.path(dirw, "../61.aracne/output/network.txt")
ta = read.table(fa, header = T, sep = "\t", as.is = T)
colnames(ta) = c("rid", "tid", "MI", "pvalue")
ta = ta[order(ta$MI, decreasing = T),]

tas = ta[1:10000,]
tls = tl[1:10000,]
apairs = sprintf("%s-%s", tas$rid, tas$tid)
gpairs = sprintf("%s-%s", tls$rid, tls$tid)
sum(apairs %in% gpairs)


ta2 = ddply(tas, .(rid), summarise, ntarget = length(unique(tid)))
ta2 = ta2[order(ta2$ntarget, decreasing = T), ]
ta2[1:10,]

tl2 = ddply(tl, .(rid), summarise, ntarget = length(unique(tid)))
tl2 = tl2[order(tl2$ntarget, decreasing = T),]
tl2[1:10,]

gids = tl2$rid[1:100]
nt_a = c(); nt_g = c(); nt_common = c()
for (gid in gids) {
	tids_a = unique(tas$tid[tas$rid == gid])
	tids_g = unique(tls$tid[tls$rid == gid])
	nt_a = c(nt_a, length(tids_a))
	nt_g = c(nt_g, length(tids_g))
	nt_common = c(nt_common, sum(tids_a %in% tids_g))
}
to = data.frame(rid = gids, aracne_targets = nt_a, genie3_targets = nt_g, targets_in_common = nt_common, stringsAsFactors = F)
