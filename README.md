---
marp: false
---

# Mapping patterns of mutational constraint in GPCR drug targets
## Overview
* gnomadIC (gnomAD Inference of Constraint) is a package for estimation of constraint from the gnomAD dataset, using custom annotations and aggregation groups. The constraint computation pipeline is written in [Hail 0.2](https://hail.is) and is based on the workflow from the gnomAD flagship manuscript [Karczewski et al., 2019](https://www.biorxiv.org/content/10.1101/531210v2)
* Analysis can be run using a set of chosen target genes using constraint_analysis.py --targets data/targets/targets.txt. You can either keep standard settings (aggegrating missense variants using provided PolyPhen2 labels) or add your own annotations. Custom annotations must be linked to HGVSP ids for variants of your choice. These can then be used to modify the annotation column before aggregation.
* This repo also contains notebooks required to perform the analyses detailed in our manuscript (Boardman et al, in preparation). 

## Methods
### Selection of GPCR genes 
We selected 393 non-olfactory protein coding GPCR genes in humans using the GPCRdb and IUPhar GuideToPharmacology databases by taking the full list and removing olfactory receptors and pseudogenes (Supplementary Table 1). We annotated genes with the receptor family (from GPCRdb) (Kooistra et al., 2021). Additionally, we extracted 1000 control genes at random.
### Extraction and annotation of genetic variants
We extracted the frequencies of all observed single-nucleotide variants (SNVs) in the identified GPCR and control genes from the gnomAD 2.1.1 exomes dataset hosted on Google Cloud. In addition, we extracted an enumerated table of all possible SNVs in these genes from the same source. Predicted impacts of these SNVs from the EnsEMBL VEP framework were also downloaded; possible impacts were 'non-coding', 'synonymous', 'missense', or 'loss-of-function'. Predicted loss-of-function (pLoF) variants consist of gained stops, frameshifts or splice site mutations; the LOFTEE tool was used to filter these into high and low-confidence pLoF variants. Missense variants predicted by PolyPhen2 as 'probably damaging' were classed as 'predicted pathogenic missense' (pPM); other missense variants were classed as 'predicted benign missense' (pBM). 
### Estimation of mutational constraint
We then calculated estimates of constraint using custom Python scripts and the Hail framework for genome-scale data processing (Karczewski, 2020); (Hail 0.2.). Observed polymorphisms with an allele fraction > 0.001 were removed. Using linear models based on nucleotide context, methylation level and coverage trained by Karczewski et al, we predicted the mutation rate for each possible variant. We then summed this for each type of mutation in each gene to obtain the expected frequency in a sample of our size. We then calculated the observed/expected ratio. We obtained Poisson confidence intervals for this ratio by evaluating the cumulative density function of a Poisson distribution with mean $n_{exp}c$ (Eq. 3). We refer to the constraint ratio in our sample as $n_{obs}/n_{exp}$ and to the upper bound of the confidence interval as $OEUF$ (observed-expected upper bound fraction). We also calculate the log-likelihood that a gene is unconstrained, which we refer to as
$$ \frac{dN}{dS}=  \frac{\text{number of non-synonymous substitutions}}{\text{number of synonymous substitutions}} (1)$$ 
$$ \frac{n_{obs}}{n_{exp}}(\text{gene, impact})=\frac{\text{number of observed variants}}{\text{number of expected variants}}  (2) $$
$$ n_{obs} \sim Poisson(n_{exp}c) $$
$$ \text{OEUF}=x \text{ such that }  P(C<x)=0.95 (3) $$
$$ \text{Log-likelihood unconstrained} =\text{log} [P(C>1)] =  -\text{log}[1-∫_0^1p(c)  dc]  (4) $$

### GPCR genetic disease genes
We extracted 51 GPCR genes which cause inherited diseases when mutated from the Online Encyclopaedia of Mendelian Inheritance in Man (OMIM) database {, 2021 #564}. We then queried the literature supporting this database for the inheritance pattern of the disease (autosomal dominant, autosomal recessive, or X-linked) and the functional effect of disease mutations on receptor signalling (loss-of-function or gain-of-function) (Supplementary Table 3). 
### GPCR drug target genes
Secondly, we identified GPCR genes coding for 110 proteins targeted by FDA-approved drugs, using the known activities of drugs recorded in DrugBank {Wishart, 2017 #565} and compared these to results from a recent paper by Congreve et al (Congreve et al., 2020) (Supplementary Table 2). The drug mode-of-action was obtained from DrugBank and classified as either activating (full agonist, partial agonist, positive allosteric modulator) or inactivating (antagonist, inverse agonist, negative allosteric modulator). In addition, 27 targets were associated with drugs which had other impacts on the protein, including ‘binder’, ‘potentiator’, ’inducer’ or ‘downregulator’; these were all classified as ‘other’. 
### GPCR genes with mouse phenotypes
Finally, phenotypes associated with homozygous deletion of 171 different GPCR genes in mice were annotated using data from the International Mouse Phenotyping Consortium (IMPC) (Supplementary Table 4) {Dickinson, 2016 #566}. For each gene, 14 homozygous knockout mice were phenotyped using a standardised pipeline. Genes for which knockout is associated with preweaning lethality of either complete or incomplete penetrance were classed as ‘lethal’, while genes associated with any other phenotype were classed as ‘non-lethal’. 
### Analysis and plotting
Associations between OEUF and phenotypic annotations were tested using the pandas, statsmodels and sklearn packages in python. Relationships were visualised using matplotlib and seaborn. Structures were visualised using Mol*.

## Analysis of constraint in GPCRs
* A set of 398 non-olfactory GPCR genes in humans was selected using the GPCRdb and IUPhar Guide To Pharmacology databases. 1000 control genes were selected at random from a list of protein-coding human genes. Lists of gene symbols for each can be found in `results/selected_gpcrs.txt` and `results/selected_controls.txt`
* Variants were extracted and counted by category. 30564 pPM variants and 3278 pLoF variants were observed.


<!--- 
Table of genes included in analysis (GPCRs and controls)
Table of observed & expected variant count by gene
-->

* Synonymous and benign variants were consistent with no significant constraint. By constrast, ploF and pPM variants were significantly constrained in a number of genes. Constraint is stronger in pLoF but can be more precisely estimated in pPM variants. Synonymous and benign variants are not subject to much constraint and are peaked around 1, with a few genes forming a shoulder (are these hyper-mutated?). Damaging missense mutations are more constrained on average; the peak is shifted to the left and is just below 1. Some genes are highly constrained against pLoF mutations, but for others it is imporssible to get robust estimates of constraint against pLoF mutations. 
* 
<!--- 
Table of constraint, upper/lower bounds, CI width by variant type
Figure shwoing distributions of point estimate, confidence interval, evidence for constraint
Table of genes w/evidence for constraint
-->

![](plots/oeuf_density_by_variant_class.png)
![](plots/logP_density_by_variant_class.png)

<!---
Figure showing variance of point estimate between transcripts of the same gene; the same with only functional transcripts included (>70% of canonical length)

--->

![](plots/canonical_CDS_length_by_gene.png)


* The association between constraint and disease gene status was invesigated; gene sets were taken from our adapted version of the gene sets repository. ROC curves were plotted
* 66 GPCR genes which are associated with heritable traits in humans due to loss and gain of function mutations were annotated using OMIM data and Schoeneberg et al [3], including 15 GPCR genes associated with genetic diseases due to gain-of-function mutations and 51 GPCR genes were associated with genetic diseases due to loss-of-function mutations. 10 genes are associated with diseases due to both gain-of-function and loss-of-function. For gain-of-function diseases, no significant difference was observed, while for loss-of-function diseases a significant difference was observed.
![](plots/constraint_by_disease_roc.png)
* 98 GPCR genes which are the targets of FDA-approved drugs were annotated using data from Congreve et al [4], including 75 targets of activating drugs (partial/full agonist or PAM) and 61 targets of inactivating drugs (antagonist, inverse agonist or NAM). 35 GPCR genes are the targets of both activating and inactivating drugs.
![](plots/constraint_by_drug_target_roc.png)
* The association between constraint and disease annotations was tested by training linear models on the ranks of the pLoF constraint and on the value of the damaging missense constraint. Loss-of-function and gain-of-function diseases were treated separately;  A similar analysis was carried out for genes which are the targets of FDA-approved drugs.

![](plots/aminergic_constraint_scatterplot.png)
This shows that the different aminergic families overlap in constraint, and that the two constraint metrics are strongly correlated for both

## References

[1] K. J. Karczewski et al., ‘The mutational constraint spectrum quantified from variation in 141,456 humans’, Nature, 581, pp. 434–443, 2020, doi :10.1038/s41586-020-2308-7
[2] E. V. Minikel et al., ‘Evaluating drug targets through human loss-of-function genetic variation’, Nature, 581, pp. 459–464, 2020, doi: 10.1038/s41586-020-2267-z
[3] T. Schöneberg and I. Liebscher, ‘Mutations in G Protein–Coupled Receptors: Mechanisms, Pathophysiology and Potential Therapeutic Approaches’, Pharmacological Reviews, 73, pp. 89-119, 2021, doi: https://doi.org/10.1124/pharmrev.120.000011
[4] M. Congreve, C. de Graaf, N. A. Swain, and C. G. Tate, ‘Impact of GPCR Structures on Drug Discovery’, Cell, 181, pp. 81–91, 2020, doi: 10.1016/j.cell.2020.03.003

## Directory and scripts

* notebooks/gnomad_extraction_BigQuery.ipynb 
    * purpose =  extraction of population variants in GPCR gene regions with annotated effects for all possible transcripts
    * input = data/Ensembl_Grch37_gpcr_gene_regions.csv
    * product = data/gnomad_v2.1.1_gpcr_variants_unfiltered.csv
* notebooks/gnomad_variant_filtering.ipynb
    * purpose = filtering of population variants that affect canonical transcripts & splitting into synonymous, missense, pLoF & modifiers
    * input = data/gnomad_v2.1.1_gpcr_variants_unfiltered.csv
    * outputs = data/gnomad_v2.1.1_gpcr_variants_(missense|plof|synonymous|modifier).csv
* notebooks/gnomad_missense_analysis.ipynb
    * purpose = estimation and plotting of constraint against missense variants
    * input = data/gnomad_v2.1.1_gpcr_variants_missense.csv, data/gnomad_v2.1.1_gpcr_variants_synonymous.csv
    * output = plots  

## To do
* 
* Harmonize README with slidedeck and paper
* Test hail logging 
* Organise gene sets
* 