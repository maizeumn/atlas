require(dplyr)
require(plyr)
require(GenomicRanges)

options(stringsAsFactors = FALSE)

dirw = file.path(Sys.getenv("misc2"), "grn23", "61.aracne")

### run aracne
module load java
java -Xmx8G -jar $src/aracne/Aracne.jar -e matrix.txt -o output --tfs 11.TF.txt --pvalue 1E-8 --seed 1 --calculateThreshold

java -Xmx8G -jar $src/aracne/Aracne.jar -e matrix.txt -o output --tfs 11.TF.txt --pvalue 1E-8 --seed 1

rm cmds.sh
for i in {1..100}
do
echo "java -Xmx2580M -jar $src/aracne/Aracne.jar -e matrix.txt -o output --tfs 11.TF.txt --pvalue 1E-8 --seed $i" >> cmds.sh
done

parallel --citation -j 24 < cmds.sh

java -Xmx8G -jar $src/aracne/Aracne.jar -o output --consolidate

###
fi = file.path(dirw, "output/network.txt")
ti = read.table(fi, header = T, sep = "\t", as.is = T)

length(unique(ti$Regulator))

ti2 = ddply(ti, .(Regulator), summarise, ntarget = length(unique(Target)))
ti2[order(ti2$ntarget, decreasing = T), ][1:10,]