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

## Negative and null results

These are filled from the real run in `site/data/results.json` and shown on the Validation panel, including any hypothesis marked "no clear difference" or "contradicted". They are not removed if they are unflattering.
