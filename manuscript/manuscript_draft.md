Introduction

Since ..., linkage mapping in line crosses such as intercross and backcross populations have been a way to identify chromosomal regions associated with complex traits. In standard interval mapping (and related methods, such as Haley—Knott regression or composite interval mapping (?)), one uses genotypes on a relatively sparse genetic map to infer the probability of the genotype in a grid along the genome.

Since ... linkage mapping has largely been superseded by genome-wide association studies, which rely on mapping genetic variation from a population sample. Genome-wide association studies of population samples are driven by linkage disequilibrium between causative variants and markers, and thus have potentially much finer genomic resolution than a linkage study.

In highly structured populations, such as livestock and plant populations, the independent sampling assumption is problematic, and various methods to adjust for population structure have been suggested. Linear mixed models provide a flexible modelling framework to adjust for population structure.

When genetic mapping is to be applied in an intercross population, one is faced with the choice between linkage mapping and genome-wide association methods. There are many examples in the literature of such intercross resource populations that have been analysed by linear mixed model genome-wide association studies. A decade earlier, the same population would probably have been analysed by linkage mapping.

In this paper, we simulated an F<sub>2</sub> resource population based on two divergent breeds, as well as an advanced intercross, and analysed the populations with line cross linkage mapping and linear mixed model genome-wide association, using a higher marker density. We find that, ...

# Methods

We evaluated the properties of linkage mapping and linear mixed model genome-wide associations in line cross populations by simulating an F2 resource population and an advanced intercross based on two divergent breeds. The simulations consisted of:

1.  Coalescent simulation to generate founder individuals sampled from two divergent breeds

2.  Sampling of causal variants for a quantitative trait to be mapped

3.  Sampling of neutral markers for a sparse and a dense genetic map

4.  Random mating between founders of both breeds to form an F1 generation of 100 individuals

5.  For the F2 scenario, random mating of F1 individuals to form an F2 of 500 individuals

6.  For the advanced intercross, seven additional generations of random mating to form seven generations of 100 individuals, and an eight generation of 500 individuals

7.  Linkage mapping and genome-wide association in the last generation

We performed simulations with AlphaSimR (REF).

## Coalescent simulation

## Quantitative trait and causal variants

For each trait, we selected randomly one causative variant per chromosome.

## Genetic maps

For the linkage mapping, we used a sparser genetic map of 100 markers per chromosome, totalling 1000 markers.

For the genome-wide association, we used a denser genetic map of 1000 markers per chromosome, totalling 10,000 markers.

For the advanced intercross scenario, we dropped any markers that had fixed during random mating. On average, XXX markers on the sparse map and XXX markers on the dense map fixed during the generation of the advanced intercross.

For the linkage mapping method, we used only fully informative markers (that is, markers that were differentially homozygous between all the founding individuals). On average, this reduced the number of markers on the sparse map to XXX.

## Linkage mapping

We performed linkage mapping with standard interval mapping as implemented in R/qtl (REF).

We assessed significance by performing a permutation test with 100 permutations.

We formed confidence intervals around the QTL peaks with the 1.8 LOD drop method (REF), expanded to the closest marker.

We estimated effect sizes by running fitqtl, including all detected QTL, and extracting only the additive effect.

## Genome-wide association

We performed genome-wide association with linear mixed models as implemented in rrBLUP (REF) using the EMMAX model (REF).

We assessed significance by adjusting the p-values with Benjamini-Hochberg FDR correction (REF), and using a 5% FDR.

We formed linkage disequilibrium intervals around the GWAS hits by including all markers that had R<sup>2</sup> &gt; 0.2 with the significant marker.

We estimated the effect sizes of the GWAS by running the mixed.solve model of the rrBLUP package using the significant markers. When there was more than one significant marker, we averaged the estimates.

## Metrics

We evaluated the performance of linkage mapping and genome-wide association by ...

To make the comparison between linkage mapping and genome-wide association fairer to the GWAS method, we only counted one false positive per chromosome, for a maximum of ten. This is because the single-QTL scan we employed for linkage mapping only models one QTL per chromosome.

# Results

Figure 1. Random examples of Manhattan plots and linkage plots for four simulated datasets. Lines indicate thresholds used, and the red points around the significant hits indicate the 1.8 LOD drop interval, and the R<sup>2</sup> &gt; 0.2 linkage disequilibrium interval, respectively.

Figure 2. Scatterplots between the number of false positives (x-axis), and true positives (y-axis) for each method and scenario.

Figure 3. Relationship between true effect size and estimated effect size for each method and scenario.

Table 1. Number of true positives, false positives, average length of confidence intervals, and average inflation factor for each method and scenario.

# Discussion

Despite being less than fashionable, linkage mapping remains an appropriate way to analyse quantitative traits in line cross populations.

Author contributions

MJ and DW conceived the study. MJ performed the simulations. MJ and DW wrote the paper.
