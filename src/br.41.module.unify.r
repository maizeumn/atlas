source("br.fun.R")

dirw = file.path(Sys.getenv("misc2"), "briggs")

## read GO
fg = '/home/springer/zhoux379/data/genome/Zmays_v4/61.interpro/15.tsv'
tg = read.table(fg, header = F, as.is = T, sep = "\t", quote = "")
colnames(tg) = c("gid", "goid", "goterm")
#tg = tg[tg$gid %in% gids,]
tf = unique(tg[,c("gid","goterm")])
colnames(tf)[2] = "funcat"
tfs = ddply(tf, .(funcat), summarise, size = length(gid))
tfs = tfs[order(tfs$size, decreasing = T), ]
tfs = tfs[tfs$size >= 5,]
tfs = cbind(tfs, mid = sprintf("GO.%03d", 1:nrow(tfs)))

tg2 = tf[tf$funcat %in% tfs$funcat,]
tg3 = merge(tg2, tfs[,c('funcat','mid')], by = 'funcat')
t_go = tg3[order(tg3$mid, tg3$gid), c(3,2,1)]
t_go = cbind(t_go, opt = 'GO')

## read CornCyc
fc = '/home/springer/zhoux379/data/genome/Zmays_v4/corncyc/01.tsv'
tc = read.table(fc, header = F, sep = "\t", quote = "")
gidlst = strsplit(tc$V7, split = ",")
lens = sapply(gidlst, length)
tc = data.frame(pl = rep(tc$V1, lens), 
	p1 = rep(tc$V2, lens), p2 = rep(tc$V3, lens), 
	p3 = rep(tc$V4, lens), p4 = rep(tc$V5, lens), 
	ps = rep(tc$V6, lens), gid = unlist(gidlst), stringsAsFactors = F)
table(tc$p1)

tc = tc[,c('pl','gid')]
colnames(tc)[1] = 'funcat'
tcs = ddply(tc, .(funcat), summarise, size = length(funcat))
tcs = tcs[order(tcs$size, decreasing = T), ]
tcs = tcs[tcs$size >= 5,]
tcs = cbind(tcs, mid = sprintf("CornCyc.%03d", 1:nrow(tcs)))

tc2 = tc[tc$funcat %in% tcs$funcat,]
tc3 = merge(tc2, tcs[,c('funcat','mid')], by = 'funcat')
t_corncyc = tc3[order(tc3$mid, tc3$gid), c(3,2,1)]
t_corncyc = cbind(t_corncyc, opt = "CornCyc")

### read mcl clusters
tma = data.frame()
for (gt in gts) {
	fm = sprintf("%s/51.camoco/03.modules.%s.tsv", dirw, gt)
	tm = read.table(fm, header = T, sep = "\t", as.is = T)
	colnames(tm) = c('funcat', 'gid')
	
	tms = ddply(tm, .(funcat), summarise, size = length(gid))
	tms = tms[order(tms$size, decreasing = T), ]
	tms = tms[tms$size >= 5,]
	tms = cbind(tms, mid = sprintf("MCL.%s.%03d", gt, 1:nrow(tms)))
	
	tm2 = tm[tm$funcat %in% tms$funcat,]
	tm3 = merge(tm2, tms[,c('funcat','mid')], by = 'funcat')
	tm3 = cbind(tm3, opt = sprintf("MCL.%s", gt))
	tma = rbind(tma, tm3[order(tm3$mid, tm3$gid), c(3,2,1,4)])
}
table(unique(tma[,-2])$opt)

### unify modules
to = rbind(t_go, t_corncyc, tma)
fo = file.path(dirw, "59.allmodules.tsv")
write.table(to, fo, sep = "\t", row.names = F, col.names = T, quote = F)

###
fm = file.path(dirw, "59.allmodules.tsv")
tm = read.table(fm, header = T, sep = "\t", quote = "")

grp = group_by(tm, opt, mid)
tm1 = summarise(grp, size = n())
grp = group_by(tm1, opt)
tm2 = summarise(grp, num_groups=n(), total_genes=sum(size), min_size=min(size), max_size=max(size), mean_size=mean(size), median_size=median(size))

fo = file.path(dirw, "59.modules", "10.sum.tsv")
write.table(tm2, fo, sep = "\t", row.names = F, col.names = T, quote = F)
