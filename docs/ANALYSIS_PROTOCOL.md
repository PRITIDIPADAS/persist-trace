# Analysis protocol

Frozen on 2026-11-26, before the pipeline had been run on the real data. The code that defines every hypothesis, threshold and comparison below was committed (commit 0986f30) before the first real run started. Earlier exploration used only synthetic data. Nothing here is backdated.

## Primary question

In EGFR-mutant NSCLC cell lines treated with osimertinib, which expression changes belong to the persister state rather than the acute drug response, and how much of each change remains after the drug is withdrawn?

## Data

GSE193258: PC9, NCI-H1975, HCC827, HCC2935. Conditions: DMSO control, 24 h osimertinib (500 nM), 21-day persisters, 24 h washout, long washout (3 to 10 days, chosen per line by regrowth). Three biological replicates each, 60 samples. Counts from NCBI's uniform GRCh38 processing. Inclusion requires every line and condition to have exactly three samples; otherwise the run stops.

## Primary analysis

1. Per line, DESeq2 (pydeseq2), design ~condition, Wald tests for seven contrasts. Genes kept if they reach 10 counts in at least 3 samples of that line.
2. Decomposition per line: beta is the Huber-weighted slope of P (persister vs control) on A (acute vs control) across genes with base mean above 20. S = P minus beta A, with SE = sqrt(SE_P squared + beta squared SE_A squared).
3. DerSimonian-Laird random-effects meta-analysis of S and P across the four lines; BH FDR.
4. Reversibility classes per line at FDR 0.05 and |log2FC| at least 1 (washout persistence uses half that effect size), consensus when at least 3 of 4 lines agree.

## Hypotheses

* **H1.** In leave-one-cell-line-out testing, a persister-specific signature (top 100 genes by S, anchored to significant P) separates persister from acute samples in the held-out line better than a naive treated-versus-control signature. Metric: mean AUROC across folds. Margin: 0.02.
* **H1b.** Persister-specific signature genes replicate (|z| above 2, same sign) in the held-out line at least as often as naive genes.
* **H2.** Reproducible persister genes split into more than one withdrawal class.
* **H3.** At least one transcription factor shows significant persister-specific activity (CollecTRI, ULM, FDR 0.05). Supporting only.
* **H4.** At least one candidate is labelled marker-without-support or contested.
* **H5.** Heterogeneity of persister states. Declared untestable with bulk data.
* **H6.** Exploratory rank comparison for genes with published perturbation records.

## Secondary and exploratory

Classifier benchmark (nearest centroid, L2 logistic, elastic net, random forest, XGBoost) on three feature sets. Hallmark programme scores. Robustness grid over FDR {0.01, 0.05, 0.1} and |log2FC| {0.5, 1, 1.5}; leave-one-line-out signature stability; OLS versus robust beta.

## Interpretation rules

A result is called reproducible only if it holds in at least 3 of 4 lines. Verdicts are written by code: "supported", "no clear difference" or "contradicted". Any post-hoc analysis added after the first real run is labelled post hoc in the results and in this file.

## Post-hoc changes

None yet.
