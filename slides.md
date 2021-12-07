---
marp: true
theme: default
math: katex
---

26/11/20
# Updates on GPCR constraint 
* Paper draft & timeline
* Drug target information
* IMPC phenotyping
* Isoforms
---
# Drug target information
* For up-to-date information extracted approved and investigational drugs from DrugBank
* Slight changes to drug data from Congreve 2019
---
# IMPC phenotyping
* Lethal = 'preweaning lethality'. Allow incomplete/complete penetrance
* Apparent focus on understudied receptors: ~40 orphan receptors, 8 chemokines, few aminergics
* Therefore highly relevant to focus of paper on understudied receptor. Where constraint and orphan receptors agree seems likely that these receptors are of high importance.
---
# Isoforms
* Majority (62%) of GPCRs have only one isoform, and majority of isoforms differ only at C & N terminal truncations
* Therefore expressed non-canonical isoforms are unlikely to have different effects to the canonical
* Might be worth looking at some GPCRs that do show functional alternative splicing e.g. CNR1
---
# Follow up from Wednesday meeting
- Clarification of mutation nomenclature
- Comparison of GPCRs to control genes
- Equations for estimation of constraint
- Transcript-specific constraint effects
---

# Clarification of mutation nomenclature

- Constraint inferred against 2 classes of mutations
- pLoF (putative loss-of-function) (predicted by LOFTEE algorithm)
    - stop gained, frameshift, splice site
    - rare (often < 10/gene in GPCRs) 
    - high confidence 
- damaging missense (predicted by PolyPhen)
    - single amino acid change, predicted pathogenic 
    - common (often >50/gene in GPCRs)
    - average 85% AUC for DMS data 
- How should I refer to these carefully: loss of function vs. altered function?


---

# Constraint estimation equations

$$
\text{Mutations} \sim Poisson(n_{exp}\lambda) \\
n_{exp}\lambda \approx n_{obs} \pm \sqrt{n_{obs}}\\
\text{Constraint } \lambda\approx n_{obs}/n_{exp} \pm \sqrt{n_{obs}/n_{exp}^{2}}\\
\text{upper bound} = x \text{ such that } P(\lambda > x) = 0.05\\
\text{p value for null hypothesis} = P(\lambda > 1)
$$

---
# Weight of evidence for constraint
![height:7cm](plots/oeuf_density_by_variant_class.png)
![height:7cm](plots/logP_density_by_variant_class.png)

---

# GPCRs have slightly shorter CDS than average

![h:8.5cm](plots/canonical_CDS_length_by_gene.png)

GPCR mean | GPCR median | non-GPCR mean | non-GPCR median
--- | --- | --- | --- 
1525|	1155 | 1844	| 1419


---

# Canonical and non-canonical GPCR transcripts have similar damaging missense constraint distributions
![height:10cm](plots/mis_pphen_constraint_by_canonical_status.png)

Seems unlikely that including non-canonical transcripts would alter results 




