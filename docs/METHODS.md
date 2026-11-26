# Methods

## Counts and labels

NCBI-generated raw counts for GSE193258 (GRCh38.p13) and the matching gene annotation are downloaded with their SHA-256 checksums recorded in `build/manifests/file_manifest.csv`. Sample labels come from the GEO series matrix. Cell line and condition are parsed from titles and characteristics, and the run stops unless every line and condition has exactly three samples.

## Differential expression

Each cell line is fitted separately with pydeseq2 (design `~condition`, Cooks refit on). Seven Wald contrasts are extracted: acute, persister, persister vs acute, and both washouts against control and against persisters.

## Decomposition

For line l and gene g, P = log2FC(persister vs control) and A = log2FC(acute vs control).

    P_g = beta_l * A_g + S_g

beta_l is a Huber-weighted slope through the origin, fitted on informative genes (base mean above 20, standard errors below 1.5). The standard error of S uses sqrt(SE_P^2 + beta^2 SE_A^2). Because P and A share the control samples, their errors are positively correlated and this formula overstates the variance of S, which is the conservative direction.

Across lines, S and P are combined with a DerSimonian-Laird random-effects model. tau^2 and I^2 are reported. A gene counts as agreeing in a line if |z| > 2 with the same sign.

## Withdrawal and memory

Memory at washout time t is m(t) = log2FC(washout vs control) / log2FC(persister vs control), computed only for genes with a significant persister change. Where the long-washout day T is in the metadata, k is the least-squares solution of log m = -k t through m(0) = 1 and the points t = 1 day and t = T:

    k = -(1 * log m(1) + T * log m(T)) / (1 + T^2),   half-life = ln 2 / k

Genes whose memory is not between 0.02 and 1.5 at both points are left unfitted, because exponential decay does not describe them.

## Programmes and regulators

MSigDB hallmark sets and CollecTRI regulons are downloaded through decoupler (Enrichr fallback for hallmark sets). Scores are the t-statistic of a univariate linear model of the gene-level meta z on set membership or regulon weight (the ULM method), with BH correction.

## Benchmark

Leave one cell line out. Expression is log-CPM centred within each line using only its own sample means (no labels). Signatures and the top 2000 variable genes are chosen from the three training lines only. Models: nearest centroid, L2 logistic, elastic-net logistic, random forest, XGBoost. Metrics: balanced accuracy, macro F1, persister/acute confusion rate, bootstrap 95% intervals over the 60 held-out predictions.

## Switching model

    d/dt [S, P] = [[r_S - a, b], [a, r_P - b]] [S, P]

with drug-dependent rates. The dominant eigenvalue (tr + sqrt(tr^2 - 4 det)) / 2 sets long-run growth. Parameters are illustrative; GSE193258 contains no cell counts to fit them. The model and the memory fit are implemented a second time in `matlab/` and compared to the Python results in CI.
