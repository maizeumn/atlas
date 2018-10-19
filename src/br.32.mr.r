require(plyr)
require(ape)
require(dplyr)
require(tidyr)
require(ggplot2)

options(stringsAsFactors = FALSE)

dirw = file.path(Sys.getenv("misc2"), "grn23", "54.mr")
dirw = '/home/springer/zhoux379/scratch/briggs2/54.mr'

fi = file.path(dirw, "../36.long.filtered.tsv")
ti = read.table(fi, sep = "\t", header = T, as.is = T)

gts = c("B73", "Mo17", "B73xMo17")
gt = 'B73'
tiw = spread(ti[ti$Genotype == gt, -c(3,4)], Tissue, fpkm)

expr = t(as.matrix(tiw[,-1]))
colnames(expr) = tiw[,1]
gids = tiw[,1]

pcc.matrix = cor(asinh(expr), method = 'pearson')
pcc = pcc.matrix[lower.tri(pcc.matrix)]

ii = which(lower.tri(pcc.matrix))
colidx = as.integer((ii - 1) / nrow(pcc.matrix)) + 1
rowidx = ii - (colidx - 1) * nrow(pcc.matrix)

ng = nrow(tiw)
coex <- matrix(rep(-1, ng*ng), nrow=ng)
coex[lower.tri(coex)] = pcc
coex = t(coex)
coex[lower.tri(coex)] = pcc

rankmatrix = t(apply(-coex, 1, rank))
mr = sapply(1:length(pcc), myfunc <- function(i) 
	sqrt(rankmatrix[rowidx[i], colidx[i]] * rankmatrix[colidx[i], rowidx[i]])
)

nets = sprintf("N%d", 1:5)
params = c(5, 10, 25, 50, 100)
for (i in 1:length(nets)) {
    net = nets[i]
    ew = exp((-(mr-1)/params[i]))
    fe = sprintf("%s/02.edgeweight.%s.rda", dirw, net)
    save(ew, file = fe)

    idxs = which(ew >= 0.01)
    to = data.frame(g1 = gids[rowidx[idxs]], g2 = gids[colidx[idxs]], ew = ew[idxs])
    fo = sprintf("%s/03.edges.%s.txt", dirw, net)
    write.table(to, fo, sep = " ", row.names = F, col.names = F, quote = F)

    fc = sprintf("%s/04.clusterone.%s.txt", dirw, net)
    cmd = sprintf("java -jar $src/cluster_ine-1.0.jar -f edge_list -F plain %s > %s", fo, fc)
    system(cmd)

    fm = sprintf("%s/05.modules.%s.tsv", dirw, net)
    cmd = sprintf("mcl2tsv.py %s %s", fc, fm)
    system(cmd)
}

# 
