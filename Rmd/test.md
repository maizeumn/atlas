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
