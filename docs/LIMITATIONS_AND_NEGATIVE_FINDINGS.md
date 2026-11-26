# Limitations and negative findings

## Limits of the data

* Four cell lines, three biological replicates per condition. Leave-one-line-out testing has four folds, so every accuracy figure is uncertain.
* Bulk RNA-seq averages over the whole dish. It cannot show whether a few cycling persisters sit inside a mostly dormant population, and it cannot tell reversal apart from outgrowth of other cells.
* Long washouts were chosen by regrowth and differ by line (H1975 3 days, HCC827 4 days, PC9 7 days, HCC2935 10 days, read from the GEO sample records). Half-lives assume simple exponential decay through two time points.
* No lineage tracing, no new wet-lab work, no clinical data.

## Limits of the methods

* The decomposition assumes one genome-wide slope per line. Genes whose acute and persister responses are related in a non-linear way are not modelled.
* The standard error of S ignores the covariance that comes from shared controls. This widens it.
* Transcription-factor activity is inferred from the expression of target genes. It is supporting evidence, not independent evidence.
* The evidence library was curated by hand from abstracts and, where noted, full text. It is incomplete, and newer papers will change some labels.
* Gene sets and regulons are downloaded at run time. Their versions are recorded on the site.

## Limits that were removed in the second build

* "Bulk data cannot test heterogeneity": PC9 single cells (GSE275199) and PDX tumour cells (GSE243562) are now analysed.
* "No replication-stress data": the cell-cycle sorted persister dataset GSE335846 is now analysed.
* "Regulator analysis depends on a live download": CollecTRI and the hallmark sets are pinned in the repository.

## Negative and null results

* H1b, contradicted: persister-specific genes replicate in the held-out line less often than naive genes (56% against 89%). The specific programme varies more between lines than the shared stress does.
* Removing stress did not improve five-condition classifiers (best balanced accuracy about 0.77 either way).
* In PDX residual tumours the naive and stress-removed signatures transfer equally well (AUROC about 0.70); removing stress gives no extra edge in vivo.
* In PC9 single cells the naive signature separates exposed from unexposed cells slightly better (0.89 against 0.85).
* Only 4 genes reverse reproducibly after washout in 3 or more lines; reversal is mostly a PC9 property.
* TROP2 is not reproducibly enriched across the four lines (2 of 4 up, I² 0.99).
* Gene-level lists are threshold-sensitive (leave-one-line-out Jaccard 0.17 to 0.33).
