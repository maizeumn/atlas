RNA-Seq Manuscript Figures
================
October 26, 2018





































![](../data/49_coop/11.expr.pdf)

Figure 1. Summary of experimental design and data collected for this
study. (A) The datasets and analyses that were performed are summarized.
Genes are divided according to whether they show evidence for
differential expression (DE) between the two parents, B73 and Mo17 and
the subset of genes that exhibit single parent expression (SPE) was then
identified. The additivity of expression for all DE genes was assessed
and classified. A subset of DE genes that include sequence polymorphisms
and are expressed at sufficient levels were used to assess and classify
cis/trans regulatory variation. (B) t-SNE visualization of all 204
samples. Log2 transformed CPM values of a set of 17,433 genes expressed
(CPM \>= 1) in at least 143 (70%) samples were used in this analysis.
Principle component analysis (PCA) followed by nearest-neighbor
graph-baed clustering (t-SNE) was performed to display all samples in 2D
space. The color indicates tissue and the shape indicates genotype. (C)
The numbers of genes that are detected (Counts per Million (CPM) \>1 in
at least two samples in each tissue) in 0-23 tissues are shown.
Different colors represent the proportion out of all 23 tissues where
each gene is expressed: not expressed in any tissue (“Silent”),
expressed in less than 20% of tissues (Tissue specific), expressed in
20-80% tissues (Intermediate frequency) and those expressed in more than
80% tissues (Constitutive). (D) The proportion genes in each expression
category (defined in panel B) that are non-syntenic (relative to other
grasses including sorghum and rice) or lack any known domains
(Interproscan, see methods) were compared to the background gene set
(all genes).

![](../data/49_coop/25.DE.pdf)

Figure 2. Analysis of developmental dynamics of differential expression.
(A) The number of DE genes for each tissue is indicated by the square
symbols with the genotype exhibit higher expression indicated by color
(blue - B73 and orange - Mo17). The number of genes with single parent
expression (SPE - DE genes with expression \<0.1 CPM (Counts per
Million) in one parent) for each tissue is shown by the circle. (B) The
number of DE genes that are detected in 1-23 tissues is shown. The color
indicates which genotypes is more highly expressed as in (A) with pink
indicating genes for which some tissues exhibit higher expression for
B73 and other tissues with higher expressed for Mo17. (C) The numbers of
SPE genes that are detected in 1-23 tissues. (D) The set of DE/SPE genes
that have detectable expression in at least 10 tissues were classified
according to whether the DE/SPE pattern is observed in less than 20% of
expressed tissues (red, “tissue-specific), in 20-80% of expressed
tissues (blue, “intermediate frequency”) or more than 80% of expressed
tissues (green, “constitutive”). (E-F) Characterization of DE/SPE
patterns relative to presence of expression. Genes were first binned by
how many tissues it is expressed in (y-axis) with numbers in the
parentheses indicating the sample size (corresponding to Figure 1B).
Genes within each bin were further binned by the proportion of expressed
tissues showing DE/SPE patterns. Genes showing DE/SPE in none of the
expressed tissues were colored in gray, while genes showing DE/SPE in at
least one tissue were color coded by the proportion of tissues
exhibiting DE/SPE (G) The proportion of genes that exhibit DE or SPE
patterns that are non-syntenic, lack any known domains (Interproscan) or
lack any homologs in public databases was determined and compared for
all genes, all expressed genes (genes expressed in at least one out of
total 23 tissues), non-DE genes (genes expressed in at least 10 tissues
but not showing DE in any tissue), DE genes (genes expressed in at least
10 tissues and DE in at least one tissue) and SPE genes (genes expressed
in at least 10 tissues and SPE in at least one tissue). For the set of
DE/SPE genes we also show the frequencies for the subset of genes with
constitutive (DE/SPE in more than 80% expressed tissues) patterns.

![](../data/49_coop/25.Dom.pdf)

Figure 3. Classification of non-additive expression patterns. (A) The
distribution of scaled difference (sometimes referred to as
dominance/additivity (d/a) values) is shown. For each gene the scaled
difference value is calculated as the (F1-MP)/(HP-MP) such that a value
of -1 would indicate expression at the same level of the low parent, a
value of 0 indicates mid-parental expression and a value of 1 indicates
high parent expression levels. The d/a distributions are shown for all
DE genes (All DE), DE genes that are below 2-fold change (DE1-2), DE
genes that are 2-4 fold change between parents (DE2-4), DE genes that
are 4-8 fold change between parents (DE4-8), DE genes that are above 8
fold change between parents (DE8+). In (B) the scaled difference
distributions are shown for all DE genes compared to SPE and non-SPE
genes. (C) The proportion of genes showing different additivity patterns
(BLP: below low-parent, LP: low-parent, MP: mid-parent, HP: high-parent
and AHP: above high-parent) was determined for aforementioned gene sets.
The numbers in parentheses are total number of DE instances (passing
given thresholds) summed across 20 non-seed tissues. (D) The set of
genes that are DE in at least five tissues were used to examine the
prevalence of additivity patterns across development. Among these genes
there are 836 examples with AHP expression in at least one tissue, 3,756
HP examples, 12,957 MP examples, 4,942 LP examples and 625 BLP examples.
For each gene we determined the proportion of tissues exhibiting the
additivity pattern of interest and use a box-plot to visualize the
distribution of values.

![](../data/49_coop/25.Reg.pdf)

Figure 4. Analysis of biased allelic expression patterns and regulatory
variation classification across tissues. (A-B) Scatterplot showing the
parental B73 allele ratio (x-axis, CPMB/(CPMB/CPMM)) and hybrid B73
allele ratio (y-axis) for DE genes (A) and non-DE genes (B) in maize
root tissue. The colors represent different regulatory variation
classifications determined for each gene (see methods). (C) The
pie-chart (left) shows the proportion of all differentially expressed
genes (between the two parents) that were assigned to different
regulatory mechanisms across all tissues. The plots on the right show
the enrichment or depletion (as fold change relative to background
proportion from the left pie-chart) for subsets of genes for each
regulatory variation classification. The subsets of genes include
different levels of fold change in expression (DE2-4, DE4-8, DE8+ and
SPE), different additivity patterns (BLP, LP, MP, HP, AHP) and genes
that were characterized by previous eQTL study (in Shoot Apex Meristem)
to be regulated by only cis-eQLT(s), only trans-eQTL(s) or by both
cis-eQTL(s) and trans-eQTL(s). For each subset of genes the proportion
of each regulatory classification was compared to background proportion
(left pie-chart) with the ratio plotted as bars along x-axis. P-values
were determined using hypergeometric test (lower.tail = FALSE for
enrichment and lower.tail = TRUE for depletion) and labelled as "\*" (P
\< 0.01) or "\*\*" (P \< 0.001). (D-E) For genes that are DE and have
allele-specific expression data in at least five tissues we assessed the
consistency of regulation variation classifications. In (D) the subset
of genes that are classified into each pattern in at least one tissue
were used to assess the proportion of tissues that show this pattern. In
(E) all genes classified into a specific pattern in at least one tissue
were used to assess whether that classification was tissue specific
(showing that pattern in less than 20% of DE tissues), intermediate
frequency (20%-80% DE tissues) or constitutive (more than 80% of DE
tissues).

    ## Warning in kable_styling(., latex_options = c("repeat_header",
    ## "hold_position"), : Please specify format in kable. kableExtra can
    ## customize either HTML or LaTeX outputs. See https://haozhu233.github.io/
    ## kableExtra/ for details.

| Sample | Tissue                  | Genotype | Rep | Condition      | TotalPairs | TrimmedPairs | MappingRate | UniqueMappingRate |
| :----- | :---------------------- | :------- | --: | :------------- | ---------: | -----------: | :---------- | :---------------- |
| BR001  | blade\_v12              | B73      |   1 | Field          |  9,456,463 |    9,367,415 | 98.4%       | 91.6%             |
| BR002  | blade\_v12              | B73      |   2 | Field          | 10,614,657 |   10,530,486 | 98.5%       | 92.9%             |
| BR003  | blade\_v12              | Mo17     |   1 | Field          |  9,456,508 |    9,374,216 | 92.2%       | 87.1%             |
| BR004  | blade\_v12              | B73      |   3 | Field          | 10,983,851 |   10,894,598 | 98.4%       | 92.0%             |
| BR005  | blade\_v12              | Mo17     |   2 | Field          |  9,878,003 |    9,791,778 | 92.5%       | 86.6%             |
| BR006  | blade\_v12              | BxM      |   1 | Field          | 13,109,236 |   13,005,210 | 95.2%       | 83.8%             |
| BR007  | blade\_v12              | Mo17     |   3 | Field          | 10,136,440 |   10,049,207 | 94.2%       | 41.5%             |
| BR008  | blade\_v12              | BxM      |   2 | Field          |  9,958,057 |    9,762,863 | 91.8%       | 81.6%             |
| BR009  | blade\_v12              | BxM      |   3 | Field          | 10,368,573 |   10,283,237 | 95.1%       | 88.1%             |
| BR010  | auricle\_v12            | B73      |   1 | Field          | 12,161,908 |   12,062,687 | 97.4%       | 92.3%             |
| BR011  | auricle\_v12            | B73      |   2 | Field          |  8,756,804 |    8,681,487 | 97.7%       | 91.8%             |
| BR012  | auricle\_v12            | B73      |   3 | Field          |  9,205,753 |    9,115,520 | 98.4%       | 93.5%             |
| BR013  | auricle\_v12            | Mo17     |   1 | Field          | 12,629,990 |   12,534,713 | 92.4%       | 68.5%             |
| BR014  | auricle\_v12            | Mo17     |   2 | Field          | 12,106,840 |   12,015,331 | 91.9%       | 85.6%             |
| BR015  | auricle\_v12            | Mo17     |   3 | Field          | 12,789,833 |   12,688,046 | 91.9%       | 86.2%             |
| BR016  | auricle\_v12            | BxM      |   1 | Field          | 13,001,085 |   12,894,656 | 95.2%       | 90.2%             |
| BR017  | auricle\_v12            | BxM      |   2 | Field          | 10,302,708 |   10,221,028 | 95.0%       | 89.9%             |
| BR018  | auricle\_v12            | BxM      |   3 | Field          | 10,717,426 |   10,612,110 | 93.9%       | 88.7%             |
| BR019  | sheath\_v12             | B73      |   1 | Field          | 11,363,325 |   11,252,448 | 98.3%       | 92.6%             |
| BR020  | sheath\_v12             | B73      |   2 | Field          | 11,653,190 |   11,560,292 | 98.5%       | 92.6%             |
| BR021  | sheath\_v12             | B73      |   3 | Field          | 10,397,926 |   10,313,536 | 98.3%       | 91.6%             |
| BR022  | sheath\_v12             | Mo17     |   1 | Field          | 10,300,996 |   10,223,999 | 91.0%       | 84.9%             |
| BR023  | sheath\_v12             | Mo17     |   2 | Field          | 10,610,559 |   10,504,881 | 90.6%       | 84.8%             |
| BR024  | sheath\_v12             | Mo17     |   3 | Field          | 10,818,591 |   10,732,498 | 91.1%       | 85.2%             |
| BR025  | sheath\_v12             | BxM      |   1 | Field          | 10,460,515 |   10,317,769 | 89.3%       | 80.0%             |
| BR026  | sheath\_v12             | BxM      |   2 | Field          | 13,454,079 |   13,205,917 | 90.2%       | 84.6%             |
| BR027  | sheath\_v12             | BxM      |   3 | Field          | 10,947,765 |   10,796,345 | 93.0%       | 85.9%             |
| BR028  | internode\_v12          | B73      |   1 | Field          | 11,704,246 |   11,566,394 | 97.2%       | 92.1%             |
| BR032  | internode\_v12          | Mo17     |   1 | Field          | 18,897,584 |   18,677,705 | 90.6%       | 84.7%             |
| BR030  | internode\_v12          | B73      |   2 | Field          | 11,511,692 |   11,382,708 | 98.1%       | 91.9%             |
| BR031  | internode\_v12          | Mo17     |   2 | Field          | 11,675,275 |   11,554,947 | 90.9%       | 85.9%             |
| BR029  | internode\_v12          | B73      |   3 | Field          |  7,280,863 |    7,207,320 | 97.4%       | 90.0%             |
| BR033  | internode\_v12          | Mo17     |   3 | Field          | 11,763,511 |   11,618,868 | 90.2%       | 85.3%             |
| BR034  | internode\_v12          | BxM      |   1 | Field          | 11,446,013 |   11,309,805 | 94.4%       | 89.6%             |
| BR035  | internode\_v12          | BxM      |   2 | Field          | 10,813,338 |   10,680,875 | 94.0%       | 88.2%             |
| BR036  | internode\_v12          | BxM      |   3 | Field          | 12,185,937 |   12,056,501 | 94.8%       | 89.3%             |
| BR037  | tassel\_v12             | B73      |   1 | Field          | 11,869,598 |   11,767,998 | 98.2%       | 87.0%             |
| BR038  | tassel\_v12             | B73      |   2 | Field          | 13,057,367 |   12,963,185 | 98.2%       | 80.5%             |
| BR039  | tassel\_v12             | B73      |   3 | Field          | 11,971,434 |   11,876,376 | 98.4%       | 86.7%             |
| BR040  | tassel\_v12             | Mo17     |   1 | Field          | 10,198,649 |   10,119,317 | 92.8%       | 78.9%             |
| BR041  | tassel\_v12             | Mo17     |   2 | Field          |  9,541,534 |    9,458,818 | 92.5%       | 86.3%             |
| BR042  | tassel\_v12             | Mo17     |   3 | Field          | 12,122,651 |   12,028,360 | 92.8%       | 73.7%             |
| BR043  | tassel\_v12             | BxM      |   1 | Field          | 11,577,660 |   11,490,867 | 94.8%       | 87.3%             |
| BR044  | tassel\_v12             | BxM      |   2 | Field          |  9,209,055 |    9,080,007 | 94.3%       | 81.8%             |
| BR045  | tassel\_v12             | BxM      |   3 | Field          | 11,345,207 |   11,254,362 | 94.8%       | 87.8%             |
| BR046  | ear\_v14                | B73      |   1 | Field          | 10,279,754 |   10,201,193 | 97.1%       | 77.0%             |
| BR047  | ear\_v14                | B73      |   2 | Field          | 10,892,889 |   10,799,410 | 98.2%       | 88.4%             |
| BR048  | ear\_v14                | B73      |   3 | Field          | 10,423,251 |   10,340,695 | 98.2%       | 87.5%             |
| BR049  | ear\_v14                | Mo17     |   1 | Field          | 11,158,750 |   11,087,463 | 93.4%       | 88.3%             |
| BR050  | ear\_v14                | Mo17     |   2 | Field          | 10,615,140 |   10,544,171 | 92.6%       | 87.5%             |
| BR051  | ear\_v14                | Mo17     |   3 | Field          | 10,130,338 |   10,058,913 | 93.5%       | 88.5%             |
| BR052  | ear\_v14                | BxM      |   1 | Field          |  9,779,316 |    9,705,723 | 96.0%       | 90.8%             |
| BR053  | ear\_v14                | BxM      |   2 | Field          | 10,765,471 |   10,691,005 | 96.1%       | 90.7%             |
| BR054  | ear\_v14                | BxM      |   3 | Field          | 10,866,113 |   10,786,414 | 96.0%       | 90.8%             |
| BR055  | silk\_0DAP              | B73      |   1 | Field          | 10,460,214 |   10,354,824 | 98.3%       | 92.9%             |
| BR056  | silk\_0DAP              | B73      |   2 | Field          | 10,094,098 |   10,015,873 | 98.3%       | 87.9%             |
| BR057  | silk\_0DAP              | B73      |   3 | Field          |  9,351,362 |    9,259,488 | 98.2%       | 92.8%             |
| BR058  | silk\_0DAP              | Mo17     |   1 | Field          |  9,154,814 |    9,077,113 | 91.2%       | 84.3%             |
| BR059  | silk\_0DAP              | Mo17     |   2 | Field          |  9,271,623 |    9,184,908 | 91.2%       | 85.6%             |
| BR060  | silk\_0DAP              | Mo17     |   3 | Field          |  8,610,999 |    8,539,340 | 91.1%       | 86.0%             |
| BR061  | silk\_0DAP              | BxM      |   1 | Field          | 12,378,233 |   12,274,018 | 94.4%       | 88.3%             |
| BR062  | silk\_0DAP              | BxM      |   2 | Field          | 11,117,019 |   10,918,321 | 93.2%       | 87.6%             |
| BR063  | silk\_0DAP              | BxM      |   3 | Field          | 12,831,768 |   12,714,257 | 95.0%       | 89.4%             |
| BR064  | spikelets\_0DAP         | B73      |   1 | Field          | 12,862,409 |   12,756,383 | 97.9%       | 92.1%             |
| BR065  | spikelets\_0DAP         | B73      |   2 | Field          | 13,140,565 |   13,016,984 | 98.5%       | 92.8%             |
| BR066  | spikelets\_0DAP         | B73      |   3 | Field          | 12,366,130 |   12,261,404 | 98.5%       | 92.9%             |
| BR067  | spikelets\_0DAP         | Mo17     |   1 | Field          | 12,302,995 |   12,214,346 | 92.0%       | 85.9%             |
| BR068  | spikelets\_0DAP         | Mo17     |   2 | Field          | 12,364,997 |   12,266,006 | 92.3%       | 85.0%             |
| BR069  | spikelets\_0DAP         | Mo17     |   3 | Field          | 12,363,396 |   12,261,701 | 92.5%       | 86.7%             |
| BR070  | spikelets\_0DAP         | BxM      |   1 | Field          | 10,329,056 |   10,246,209 | 95.4%       | 89.4%             |
| BR071  | spikelets\_0DAP         | BxM      |   2 | Field          | 12,278,485 |   12,184,457 | 95.4%       | 89.8%             |
| BR072  | spikelets\_0DAP         | BxM      |   3 | Field          | 12,426,693 |   12,323,355 | 95.7%       | 89.9%             |
| BR073  | husk\_0DAP              | B73      |   1 | Field          | 12,468,897 |   12,338,745 | 98.3%       | 92.5%             |
| BR074  | husk\_0DAP              | B73      |   2 | Field          | 11,486,069 |   11,389,989 | 98.3%       | 92.2%             |
| BR075  | husk\_0DAP              | B73      |   3 | Field          | 11,705,643 |   11,598,042 | 98.3%       | 92.1%             |
| BR076  | husk\_0DAP              | Mo17     |   1 | Field          | 12,818,182 |   12,713,872 | 92.1%       | 83.4%             |
| BR077  | husk\_0DAP              | Mo17     |   2 | Field          | 12,761,502 |   12,635,684 | 91.4%       | 85.5%             |
| BR078  | husk\_0DAP              | Mo17     |   3 | Field          | 12,826,901 |   12,717,414 | 91.9%       | 85.6%             |
| BR079  | husk\_0DAP              | BxM      |   1 | Field          | 12,503,981 |   12,395,438 | 94.7%       | 89.1%             |
| BR080  | husk\_0DAP              | BxM      |   2 | Field          | 11,607,643 |   11,406,667 | 93.7%       | 88.0%             |
| BR081  | husk\_0DAP              | BxM      |   3 | Field          | 13,790,861 |   13,666,823 | 95.3%       | 89.8%             |
| BR082  | tasselstem\_0DAP        | B73      |   1 | Field          | 12,451,564 |   12,349,251 | 97.8%       | 92.1%             |
| BR083  | tasselstem\_0DAP        | B73      |   2 | Field          | 12,221,050 |   12,112,847 | 98.3%       | 92.7%             |
| BR084  | tasselstem\_0DAP        | B73      |   3 | Field          | 11,567,848 |   11,471,842 | 98.0%       | 88.0%             |
| BR085  | tasselstem\_0DAP        | Mo17     |   1 | Field          |  9,540,165 |    9,463,815 | 90.5%       | 83.1%             |
| BR086  | tasselstem\_0DAP        | Mo17     |   2 | Field          |  9,089,631 |    9,011,444 | 91.1%       | 85.8%             |
| BR087  | tasselstem\_0DAP        | Mo17     |   3 | Field          |  9,244,929 |    9,157,904 | 91.2%       | 85.7%             |
| BR088  | tasselstem\_0DAP        | BxM      |   1 | Field          |  8,876,713 |    8,797,378 | 94.5%       | 88.0%             |
| BR089  | tasselstem\_0DAP        | BxM      |   2 | Field          |  9,021,520 |    8,943,627 | 93.9%       | 88.0%             |
| BR090  | tasselstem\_0DAP        | BxM      |   3 | Field          |  9,072,250 |    8,993,652 | 94.2%       | 87.0%             |
| BR091  | floret\_0DAP            | B73      |   1 | Field          | 10,046,688 |    9,917,126 | 98.1%       | 89.5%             |
| BR092  | floret\_0DAP            | B73      |   2 | Field          | 10,344,979 |   10,247,053 | 98.3%       | 88.6%             |
| BR093  | floret\_0DAP            | B73      |   3 | Field          |  9,586,452 |    9,469,113 | 98.1%       | 89.7%             |
| BR094  | floret\_0DAP            | Mo17     |   1 | Field          |  9,591,540 |    9,496,979 | 90.9%       | 83.0%             |
| BR095  | floret\_0DAP            | Mo17     |   2 | Field          | 10,171,616 |   10,050,093 | 90.1%       | 83.2%             |
| BR096  | floret\_0DAP            | Mo17     |   3 | Field          |  8,968,665 |    8,870,021 | 91.6%       | 84.9%             |
| BR097  | floret\_0DAP            | BxM      |   1 | Field          |  8,971,648 |    8,871,487 | 95.1%       | 89.1%             |
| BR098  | floret\_0DAP            | BxM      |   2 | Field          |  9,125,085 |    8,871,603 | 93.3%       | 86.3%             |
| BR099  | floret\_0DAP            | BxM      |   3 | Field          | 10,114,073 |   10,000,848 | 95.0%       | 86.9%             |
| BR118  | kernel\_14DAP           | B73      |   1 | Field          | 10,936,328 |   10,843,315 | 98.3%       | 87.3%             |
| BR119  | kernel\_14DAP           | B73      |   2 | Field          | 11,056,044 |   10,956,542 | 98.4%       | 89.0%             |
| BR120  | kernel\_14DAP           | B73      |   3 | Field          | 10,542,981 |   10,452,024 | 98.3%       | 86.7%             |
| BR121  | kernel\_14DAP           | Mo17     |   1 | Field          | 13,288,857 |   13,186,732 | 94.6%       | 85.2%             |
| BR122  | kernel\_14DAP           | Mo17     |   2 | Field          | 13,244,937 |   13,134,240 | 94.5%       | 86.1%             |
| BR123  | kernel\_14DAP           | Mo17     |   3 | Field          | 11,612,393 |   11,504,790 | 94.3%       | 83.0%             |
| BR124  | kernel\_14DAP           | BxM      |   1 | Field          | 12,287,598 |   12,188,660 | 97.5%       | 82.4%             |
| BR125  | kernel\_14DAP           | BxM      |   2 | Field          | 11,390,752 |   11,303,494 | 97.4%       | 82.3%             |
| BR126  | kernel\_14DAP           | BxM      |   3 | Field          | 10,802,003 |   10,718,433 | 97.5%       | 82.8%             |
| BR100  | flagleaf\_0DAP          | B73      |   1 | Field          |  8,870,907 |    8,779,729 | 98.3%       | 90.6%             |
| BR101  | flagleaf\_0DAP          | B73      |   2 | Field          |  8,927,886 |    8,852,968 | 98.5%       | 92.2%             |
| BR102  | flagleaf\_0DAP          | B73      |   3 | Field          |  9,310,919 |    9,227,730 | 98.6%       | 92.3%             |
| BR103  | flagleaf\_0DAP          | Mo17     |   1 | Field          | 12,454,393 |   12,351,934 | 92.0%       | 87.2%             |
| BR104  | flagleaf\_0DAP          | Mo17     |   2 | Field          | 11,124,928 |   11,016,585 | 92.4%       | 86.5%             |
| BR105  | flagleaf\_0DAP          | Mo17     |   3 | Field          |  9,139,066 |    9,063,499 | 92.0%       | 86.9%             |
| BR106  | flagleaf\_0DAP          | BxM      |   1 | Field          |  8,237,380 |    8,161,554 | 95.1%       | 89.4%             |
| BR107  | flagleaf\_0DAP          | BxM      |   2 | Field          | 11,375,396 |   11,265,249 | 95.3%       | 89.9%             |
| BR108  | flagleaf\_0DAP          | BxM      |   3 | Field          |  8,739,031 |    8,651,672 | 95.2%       | 89.0%             |
| BR109  | root\_0DAP              | B73      |   1 | Field          | 12,780,458 |   12,681,574 | 98.7%       | 94.0%             |
| BR110  | root\_0DAP              | B73      |   2 | Field          | 11,071,413 |   10,726,484 | 95.8%       | 91.2%             |
| BR111  | root\_0DAP              | B73      |   3 | Field          | 14,352,722 |   14,232,937 | 98.7%       | 94.2%             |
| BR112  | root\_0DAP              | Mo17     |   1 | Field          | 13,432,400 |   13,337,049 | 91.9%       | 87.2%             |
| BR113  | root\_0DAP              | Mo17     |   2 | Field          | 11,241,291 |   11,148,033 | 91.9%       | 87.7%             |
| BR114  | root\_0DAP              | Mo17     |   3 | Field          | 11,748,397 |   11,657,233 | 91.8%       | 87.4%             |
| BR115  | root\_0DAP              | BxM      |   1 | Field          | 12,447,495 |   12,322,988 | 94.7%       | 90.3%             |
| BR116  | root\_0DAP              | BxM      |   2 | Field          | 13,501,613 |   13,391,660 | 95.0%       | 90.7%             |
| BR117  | root\_0DAP              | BxM      |   3 | Field          | 12,014,758 |   11,917,490 | 95.3%       | 90.9%             |
| BR130  | endosperm\_14DAP        | B73      |   1 | Field          |  7,971,614 |    7,908,501 | 98.4%       | 87.7%             |
| BR131  | endosperm\_14DAP        | B73      |   2 | Field          |  8,354,751 |    8,297,059 | 98.6%       | 87.4%             |
| BR132  | endosperm\_14DAP        | B73      |   3 | Field          |  8,106,969 |    8,043,877 | 98.4%       | 87.1%             |
| BR133  | endosperm\_14DAP        | Mo17     |   1 | Field          |  7,996,705 |    7,944,390 | 95.1%       | 86.0%             |
| BR134  | endosperm\_14DAP        | Mo17     |   2 | Field          |  7,111,069 |    7,056,236 | 95.0%       | 86.1%             |
| BR135  | endosperm\_14DAP        | Mo17     |   3 | Field          |  7,053,101 |    7,003,061 | 94.9%       | 85.4%             |
| BR136  | endosperm\_14DAP        | BxM      |   1 | Field          |  7,342,993 |    7,293,685 | 97.3%       | 80.3%             |
| BR137  | endosperm\_14DAP        | BxM      |   2 | Field          |  7,408,587 |    7,282,363 | 95.5%       | 82.6%             |
| BR138  | endosperm\_14DAP        | BxM      |   3 | Field          |  8,295,102 |    8,242,013 | 97.5%       | 81.5%             |
| BR142  | embryo\_27DAP           | B73      |   1 | Field          |  7,372,483 |    7,316,640 | 98.5%       | 86.8%             |
| BR143  | embryo\_27DAP           | B73      |   2 | Field          |  7,042,198 |    6,983,724 | 98.7%       | 90.0%             |
| BR144  | embryo\_27DAP           | B73      |   3 | Field          |  6,553,609 |    6,498,382 | 98.3%       | 87.1%             |
| BR145  | embryo\_27DAP           | Mo17     |   1 | Field          |  7,077,358 |    7,002,686 | 91.6%       | 83.7%             |
| BR146  | embryo\_27DAP           | Mo17     |   2 | Field          |  7,533,408 |    7,451,111 | 91.3%       | 83.5%             |
| BR147  | embryo\_27DAP           | Mo17     |   3 | Field          |  8,261,322 |    8,162,520 | 92.3%       | 47.1%             |
| BR148  | embryo\_27DAP           | BxM      |   1 | Field          |  7,602,835 |    7,530,655 | 94.5%       | 84.1%             |
| BR149  | embryo\_27DAP           | BxM      |   2 | Field          |  8,456,705 |    8,336,995 | 94.2%       | 85.4%             |
| BR150  | embryo\_27DAP           | BxM      |   3 | Field          |  8,081,181 |    7,993,598 | 94.9%       | 83.0%             |
| BR154  | endosperm\_27DAP        | B73      |   1 | Field          | 17,736,115 |   17,600,317 | 98.3%       | 68.2%             |
| BR155  | endosperm\_27DAP        | B73      |   2 | Field          | 17,169,716 |   17,057,377 | 98.3%       | 68.4%             |
| BR156  | endosperm\_27DAP        | B73      |   3 | Field          | 18,730,409 |   18,611,675 | 98.4%       | 69.2%             |
| BR157  | endosperm\_27DAP        | Mo17     |   1 | Field          | 14,649,158 |   14,549,188 | 96.2%       | 66.5%             |
| BR158  | endosperm\_27DAP        | Mo17     |   2 | Field          | 15,811,476 |   15,702,263 | 96.4%       | 66.1%             |
| BR159  | endosperm\_27DAP        | Mo17     |   3 | Field          | 13,053,684 |   12,967,048 | 96.6%       | 67.5%             |
| BR160  | endosperm\_27DAP        | BxM      |   1 | Field          | 14,814,247 |   14,719,774 | 97.8%       | 71.6%             |
| BR161  | endosperm\_27DAP        | BxM      |   2 | Field          | 14,413,309 |   14,133,430 | 95.7%       | 67.7%             |
| BR162  | endosperm\_27DAP        | BxM      |   3 | Field          | 14,883,939 |   14,783,781 | 97.7%       | 68.9%             |
| BR166  | coleoptile\_tip         | B73      |   1 | Growth chamber | 12,267,403 |   12,126,234 | 98.0%       | 89.7%             |
| BR167  | coleoptile\_tip         | B73      |   2 | Growth chamber | 11,485,549 |   11,384,771 | 98.4%       | 93.3%             |
| BR168  | coleoptile\_tip         | B73      |   3 | Growth chamber | 12,399,486 |   12,290,630 | 98.5%       | 91.7%             |
| BR169  | coleoptile\_tip         | Mo17     |   1 | Growth chamber | 11,947,303 |   11,849,179 | 92.4%       | 87.4%             |
| BR170  | coleoptile\_tip         | Mo17     |   2 | Growth chamber | 12,242,354 |   12,126,048 | 92.5%       | 87.4%             |
| BR171  | coleoptile\_tip         | Mo17     |   3 | Growth chamber | 13,259,951 |   13,138,932 | 92.4%       | 86.6%             |
| BR172  | coleoptile\_tip         | BxM      |   1 | Growth chamber | 12,391,702 |   12,284,694 | 95.3%       | 90.0%             |
| BR173  | coleoptile\_tip         | BxM      |   2 | Growth chamber | 12,048,829 |   11,802,713 | 94.0%       | 88.8%             |
| BR174  | coleoptile\_tip         | BxM      |   3 | Growth chamber | 12,187,226 |   12,068,545 | 95.1%       | 73.1%             |
| BR175  | radicle\_root           | B73      |   1 | Growth chamber | 12,295,498 |   12,171,020 | 97.0%       | 92.8%             |
| BR176  | radicle\_root           | B73      |   2 | Growth chamber | 11,534,982 |   11,416,060 | 98.1%       | 93.0%             |
| BR177  | radicle\_root           | B73      |   3 | Growth chamber | 10,646,209 |   10,526,034 | 97.8%       | 92.8%             |
| BR178  | radicle\_root           | Mo17     |   1 | Growth chamber | 12,066,875 |   11,965,649 | 90.8%       | 84.2%             |
| BR179  | radicle\_root           | Mo17     |   2 | Growth chamber | 11,661,817 |   11,558,386 | 90.2%       | 85.7%             |
| BR180  | radicle\_root           | Mo17     |   3 | Growth chamber | 10,959,894 |   10,841,662 | 90.1%       | 85.5%             |
| BR181  | radicle\_root           | BxM      |   1 | Growth chamber | 13,122,821 |   13,004,621 | 94.2%       | 89.9%             |
| BR182  | radicle\_root           | BxM      |   2 | Growth chamber | 12,468,832 |   12,358,623 | 94.8%       | 90.2%             |
| BR183  | radicle\_root           | BxM      |   3 | Growth chamber | 11,587,217 |   11,485,773 | 94.1%       | 89.5%             |
| BR184  | embryo\_imbibedseed     | B73      |   1 | Growth chamber | 16,521,544 |   16,009,522 | 98.4%       | 86.6%             |
| BR242  | embryo\_imbibedseed     | B73      |   2 | Growth chamber | 18,100,445 |   17,634,613 | 98.4%       | 87.7%             |
| BR243  | embryo\_imbibedseed     | B73      |   3 | Growth chamber | 16,539,617 |   16,051,406 | 98.4%       | 88.7%             |
| BR245  | embryo\_imbibedseed     | Mo17     |   1 | Growth chamber | 15,060,983 |   14,635,237 | 91.8%       | 78.2%             |
| BR188  | embryo\_imbibedseed     | Mo17     |   2 | Growth chamber | 19,136,137 |   18,583,167 | 91.9%       | 77.8%             |
| BR227  | embryo\_imbibedseed     | Mo17     |   3 | Growth chamber | 10,606,785 |   10,249,652 | 92.2%       | 62.6%             |
| BR191  | embryo\_imbibedseed     | BxM      |   1 | Growth chamber | 12,528,211 |   11,753,391 | 95.2%       | 83.8%             |
| BR192  | embryo\_imbibedseed     | BxM      |   2 | Growth chamber | 16,939,295 |   16,506,855 | 95.8%       | 84.8%             |
| BR193  | seedlingleaf\_11DAS     | B73      |   1 | Growth chamber | 12,931,365 |   12,562,955 | 98.4%       | 90.7%             |
| BR194  | seedlingleaf\_11DAS     | B73      |   2 | Growth chamber | 12,993,398 |   12,476,399 | 98.7%       | 91.5%             |
| BR195  | seedlingleaf\_11DAS     | B73      |   3 | Growth chamber | 12,944,343 |   12,533,704 | 98.7%       | 91.8%             |
| BR196  | seedlingleaf\_11DAS     | Mo17     |   1 | Growth chamber | 13,154,046 |   12,759,277 | 91.3%       | 83.2%             |
| BR197  | seedlingleaf\_11DAS     | Mo17     |   2 | Growth chamber | 14,456,660 |   14,047,653 | 91.8%       | 85.8%             |
| BR198  | seedlingleaf\_11DAS     | Mo17     |   3 | Growth chamber | 12,627,543 |   12,111,783 | 91.5%       | 82.8%             |
| BR199  | seedlingleaf\_11DAS     | BxM      |   1 | Growth chamber | 13,369,929 |   12,922,403 | 95.2%       | 87.0%             |
| BR200  | seedlingleaf\_11DAS     | BxM      |   2 | Growth chamber | 12,269,095 |   11,890,258 | 94.9%       | 88.5%             |
| BR201  | seedlingleaf\_11DAS     | BxM      |   3 | Growth chamber | 16,931,135 |   16,451,260 | 95.4%       | 86.3%             |
| BR202  | seedlingroot\_11DAS     | B73      |   1 | Growth chamber | 14,286,923 |   13,837,005 | 98.3%       | 92.4%             |
| BR203  | seedlingroot\_11DAS     | B73      |   2 | Growth chamber | 13,799,806 |   13,493,069 | 98.5%       | 93.1%             |
| BR204  | seedlingroot\_11DAS     | B73      |   3 | Growth chamber | 12,973,928 |   12,661,827 | 98.6%       | 93.2%             |
| BR187  | embryo\_imbibedseed     | Mo17     |   4 | Growth chamber | 12,018,006 |   11,651,375 | 90.8%       | 74.2%             |
| BR206  | seedlingroot\_11DAS     | Mo17     |   1 | Growth chamber | 11,969,669 |   11,569,090 | 90.5%       | 84.8%             |
| BR208  | seedlingroot\_11DAS     | BxM      |   1 | Growth chamber | 14,247,885 |   13,919,656 | 95.2%       | 88.8%             |
| BR209  | seedlingroot\_11DAS     | BxM      |   2 | Growth chamber | 12,311,137 |   11,773,916 | 95.0%       | 89.2%             |
| BR211  | seedlingmeristem\_11DAS | B73      |   1 | Growth chamber | 15,300,678 |   14,943,561 | 98.7%       | 93.6%             |
| BR212  | seedlingmeristem\_11DAS | B73      |   2 | Growth chamber | 13,228,641 |   12,892,735 | 98.9%       | 93.1%             |
| BR213  | seedlingmeristem\_11DAS | B73      |   3 | Growth chamber | 14,690,786 |   14,301,943 | 98.9%       | 93.1%             |
| BR214  | seedlingmeristem\_11DAS | Mo17     |   1 | Growth chamber | 17,091,762 |   16,707,591 | 93.7%       | 88.4%             |
| BR215  | seedlingmeristem\_11DAS | Mo17     |   2 | Growth chamber | 14,961,784 |   14,610,501 | 93.7%       | 88.6%             |
| BR216  | seedlingmeristem\_11DAS | Mo17     |   3 | Growth chamber | 14,085,672 |   13,755,241 | 93.9%       | 88.9%             |
| BR217  | seedlingmeristem\_11DAS | BxM      |   1 | Growth chamber | 15,456,137 |   15,088,536 | 96.4%       | 91.4%             |
| BR218  | seedlingmeristem\_11DAS | BxM      |   2 | Growth chamber | 14,791,193 |   14,459,872 | 96.5%       | 89.3%             |
| BR219  | seedlingmeristem\_11DAS | BxM      |   3 | Growth chamber | 13,255,400 |   12,905,133 | 96.4%       | 91.0%             |

    ## Warning in kable_styling(., latex_options = c("repeat_header",
    ## "scale_down", : Please specify format in kable. kableExtra can customize
    ## either HTML or LaTeX outputs. See https://haozhu233.github.io/kableExtra/
    ## for details.

    ## Warning in column_spec(., 1, bold = T): Please specify format in kable.
    ## kableExtra can customize either HTML or LaTeX outputs. See https://
    ## haozhu233.github.io/kableExtra/ for details.

    ## Warning in collapse_rows(., columns = 1:2, latex_hline = "major", valign
    ## = "top", : Please specify format in kable. kableExtra can customize either
    ## HTML or LaTeX outputs. See https://haozhu233.github.io/kableExtra/ for
    ## details.

|                            | Gene set     | GO           | P-value  | GO name                                                           |
| :------------------------- | :----------- | :----------- | :------- | :---------------------------------------------------------------- |
| All genes                  | Constitutive | <GO:0005739> | 3.08e-22 | mitochondrion                                                     |
| All genes                  | Constitutive | <GO:0022625> | 1.45e-11 | cytosolic large ribosomal subunit                                 |
| All genes                  | Constitutive | <GO:0006412> | 4.17e-11 | translation                                                       |
| All genes                  | Constitutive | <GO:0005774> | 7.24e-10 | vacuolar membrane                                                 |
| All genes                  | Constitutive | <GO:0009793> | 2.54e-09 | embryo development ending in seed dormancy                        |
| All genes                  | Constitutive | <GO:0000398> | 3.25e-06 | mRNA splicing, via spliceosome                                    |
| All genes                  | Constitutive | <GO:0005743> | 4.24e-05 | mitochondrial inner membrane                                      |
| All genes                  | Silent       | <GO:0005654> | 4.47e-09 | nucleoplasm                                                       |
| All genes                  | Silent       | <GO:0043161> | 3.49e-06 | proteasome-mediated ubiquitin-dependent protein catabolic process |
| All genes                  | Silent       | <GO:0005737> | 2.47e-05 | cytoplasm                                                         |
| non-DEGs btw. B73 and Mo17 | Above-Parent | <GO:0022625> | 1.36e-08 | cytosolic large ribosomal subunit                                 |
| non-DEGs btw. B73 and Mo17 | Above-Parent | <GO:0042256> | 1.43e-05 | mature ribosome assembly                                          |
| non-DEGs btw. B73 and Mo17 | Above-Parent | <GO:0002181> | 4.30e-05 | cytoplasmic translation                                           |
| non-DEGs btw. B73 and Mo17 | Above-Parent | <GO:0009941> | 1.85e-04 | chloroplast envelope                                              |
| non-DEGs btw. B73 and Mo17 | Above-Parent | <GO:0010287> | 2.02e-04 | plastoglobule                                                     |
| non-DEGs btw. B73 and Mo17 | Above-Parent | <GO:0009535> | 1.93e-03 | chloroplast thylakoid membrane                                    |
| non-DEGs btw. B73 and Mo17 | Above-Parent | <GO:0005774> | 2.09e-03 | vacuolar membrane                                                 |
| non-DEGs btw. B73 and Mo17 | Above-Parent | <GO:0006412> | 4.39e-03 | translation                                                       |
| non-DEGs btw. B73 and Mo17 | Below-Parent | <GO:0022625> | 1.09e-22 | cytosolic large ribosomal subunit                                 |
| non-DEGs btw. B73 and Mo17 | Below-Parent | <GO:0042256> | 2.63e-20 | mature ribosome assembly                                          |
| non-DEGs btw. B73 and Mo17 | Below-Parent | <GO:0002181> | 4.38e-13 | cytoplasmic translation                                           |
| non-DEGs btw. B73 and Mo17 | Below-Parent | <GO:0022627> | 5.96e-13 | cytosolic small ribosomal subunit                                 |
| non-DEGs btw. B73 and Mo17 | Below-Parent | <GO:0006412> | 6.76e-13 | translation                                                       |
| non-DEGs btw. B73 and Mo17 | Below-Parent | <GO:0048046> | 4.04e-12 | apoplast                                                          |
| non-DEGs btw. B73 and Mo17 | Below-Parent | <GO:0005794> | 9.05e-12 | Golgi apparatus                                                   |
| non-DEGs btw. B73 and Mo17 | Below-Parent | <GO:0042788> | 8.01e-09 | polysomal ribosome                                                |
| non-DEGs btw. B73 and Mo17 | Below-Parent | <GO:0005774> | 4.04e-08 | vacuolar membrane                                                 |
| non-DEGs btw. B73 and Mo17 | Below-Parent | <GO:0005618> | 1.88e-06 | cell wall                                                         |
| non-DEGs btw. B73 and Mo17 | Below-Parent | <GO:0009506> | 3.10e-05 | plasmodesma                                                       |
| non-DEGs btw. B73 and Mo17 | Below-Parent | <GO:0046686> | 7.22e-05 | response to cadmium ion                                           |
| non-DEGs btw. B73 and Mo17 | Below-Parent | <GO:0000028> | 4.78e-04 | ribosomal small subunit assembly                                  |
| non-DEGs btw. B73 and Mo17 | Below-Parent | <GO:0009735> | 1.23e-03 | response to cytokinin                                             |
| non-DEGs btw. B73 and Mo17 | Below-Parent | <GO:0005886> | 1.23e-03 | plasma membrane                                                   |
| non-DEGs btw. B73 and Mo17 | Below-Parent | <GO:0009409> | 2.49e-03 | response to cold                                                  |
| non-DEGs btw. B73 and Mo17 | Below-Parent | <GO:0009651> | 3.26e-03 | response to salt stress                                           |
| DEGs btw. B73 and Mo17     | HP/AHP       | <GO:0009535> | 7.76e-31 | chloroplast thylakoid membrane                                    |
| DEGs btw. B73 and Mo17     | HP/AHP       | <GO:0009941> | 9.82e-19 | chloroplast envelope                                              |
| DEGs btw. B73 and Mo17     | HP/AHP       | <GO:0009570> | 3.78e-16 | chloroplast stroma                                                |
| DEGs btw. B73 and Mo17     | HP/AHP       | <GO:0009768> | 6.89e-08 | photosynthesis, light harvesting in photosystem I                 |
| DEGs btw. B73 and Mo17     | HP/AHP       | <GO:0009507> | 2.61e-06 | chloroplast                                                       |
| DEGs btw. B73 and Mo17     | HP/AHP       | <GO:0009773> | 3.32e-06 | photosynthetic electron transport in photosystem I                |
| DEGs btw. B73 and Mo17     | HP/AHP       | <GO:0010287> | 6.58e-06 | plastoglobule                                                     |
| DEGs btw. B73 and Mo17     | HP/AHP       | <GO:0031977> | 1.49e-05 | thylakoid lumen                                                   |
| DEGs btw. B73 and Mo17     | HP/AHP       | <GO:0016036> | 2.32e-05 | cellular response to phosphate starvation                         |
| DEGs btw. B73 and Mo17     | HP/AHP       | <GO:0010027> | 9.02e-05 | thylakoid membrane organization                                   |
| DEGs btw. B73 and Mo17     | LP/BLP       | <GO:0005886> | 6.65e-08 | plasma membrane                                                   |
