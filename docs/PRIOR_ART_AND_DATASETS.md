# Prior art and datasets

All records below were checked against PubMed on 2026-11-26. DOIs are in `evidence/literature.yaml`.

## Prior-art comparison

| Capability | Existing work | What PersistTrace adds |
|---|---|---|
| Osimertinib DTP transcriptomes in four lines | Criscione et al., npj Precision Oncology 2022 (source of GSE193258). Clustered genes by trend, found an EMT cluster that persisted through washout and a cell-cycle cluster that tracked regrowth. | Gene-level separation of acute and persister effects with uncertainty, meta-analysis across lines, explicit withdrawal classes and half-lives, held-out-line testing |
| ML persister signature | Chang et al., Molecular Medicine Reports 2026. Four genes, knockdown impaired DTP formation. | Uses their genes as perturbation evidence; does not claim to be first |
| TROP2 in persisters | Baldacci et al., Cancer Discovery 2025 (knockout: no effect on formation); Liao et al., Cancer Cell 2026 (contributes to maintenance) | Shows the disagreement and the differing outcomes in one passport |
| Persister heterogeneity | Oren et al., Nature 2021: cycling and non-cycling persisters from distinct lineages in PC9 | Programme time courses; heterogeneity left to the single-cell phase |
| Mutagenesis in persisters | Isozaki et al., Nature 2023: APOBEC3A deletion delays resistance | Recorded with outcome "resistance evolution", kept apart from DTP formation |
| Single-cell resistance states | Kashima et al., Cancer Research 2021: CD74 in the drug-tolerant state | Recorded as functional evidence in its stated context |
| Replication stress | Gini et al., bioRxiv 2025: CDK4/6 overexpression, replication stress, resistance | Recorded as preprint-level context for the replication panel |

## Dataset tiers

**Tier A, analysed.** GSE193258, bulk RNA-seq, 60 samples. Counts from NCBI uniform processing (GRCh38.p13). Sample labels parsed from the series matrix and validated against the 4 x 5 x 3 design; the run stops if they do not match.

**Tier B, planned, not run.**
* GSE193259: the multi-omic SuperSeries (ATAC-seq, ChIP-seq). Assay counts are not independent biological samples.
* GSE243562: single-cell data from patient-derived xenografts, untreated and residual tumour.
* GSE275199: PC9 osimertinib arm only. The H358 sotorasib arm is a different drug and cancer and must not be pooled.
* GSE302284: six samples with mixed tissue and model contexts. Limited clinical context, not a validation cohort.

**Tier C, to be checked.** The roadmap named GSE335846 and GSE335847 (cell-cycle sorted RNA-seq and TrAEL-seq under prolonged osimertinib). These accessions and their linked paper could not be verified during this build because GEO pages were not reachable from the build environment and PubMed returned no matching paper. They stay out of every claim until checked by hand.

## Inclusion and exclusion

Included: samples whose cell line and condition parse unambiguously. Excluded: none so far. Any QC outlier (mean within-line Spearman correlation more than 3 SD below its line) is flagged on the site but kept, because removal rules were not preregistered.
