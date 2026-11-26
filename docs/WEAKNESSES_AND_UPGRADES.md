# Weaknesses in the original plan, and what changed

The original PERSIST-TRACE roadmap was ambitious and mostly sound. Reviewed against the literature and the actual data, it had problems that a sharp judge would find in a few minutes. This file lists them and what the build does about each one.

## 1. Scientific weaknesses

**The central question needs single-cell data, but the core dataset is bulk.**
"Are persisters one population or several states?" cannot be answered from GSE193258, which is bulk RNA-seq. Bulk averages hide subpopulations. The build keeps the question, answers what bulk data can answer (stress versus persistence, reversibility, cross-line reproducibility), and marks heterogeneity as hypothesis H5, "not testable with bulk RNA-seq". The single-cell datasets are listed as the next step rather than claimed.

**The TROP2 "negative control" is now contested.**
The plan used TROP2 as the textbook marker that is not a driver. That came from Baldacci et al. (Cancer Discovery 2025): TACSTD2 knockout in PC9 did not change DTP formation. In June 2026 Liao et al. (Cancer Cell) reported that TROP2 contributes to DTP maintenance. Those are different outcomes (formation versus maintenance) in different models, so both can hold. The passport now labels TROP2 "Contested mechanism" and shows both records with their outcomes. This turned a weakness into the best demo moment of the project.

**Closest prior art was not engaged.**
Chang et al. (Molecular Medicine Reports 2026) already used machine-learning feature selection to build a four-gene persister signature (BTG1, ING4, PNRC1, CALCOCO1) and validated it by knockdown. The project cannot claim to be first at ML persister signatures. Those four genes are now in the evidence library and shown on the site.

**"Reversibility" was undefined at the population level.**
A bulk signal returning toward baseline fits reversal, but it fits outgrowth of other cells just as well. The build uses explicit, population-level classes and says so on the page.

**Long washouts are not the same length.**
The source paper chose 3 to 10 days per line based on regrowth. The pipeline reads the duration from GEO metadata when present and only fits decay half-lives where it has a real number.

**Perturbation evidence for persisters is sparse.**
DepMap dependency is ordinary cell-line fitness, not persister fitness, and LINCS signatures come from other contexts. The evidence library therefore uses hand-checked published experiments and records what each one measured (formation, maintenance, vulnerability, apoptosis, resistance). A knockdown that only changes acute apoptosis is never promoted to "dependency".

**Foundation models do not fit the primary data.**
scGPT and Geneformer are single-cell models. Running them on 60 bulk samples would add a buzzword and no evidence. The benchmark uses strong simple baselines instead and keeps foundation models for the single-cell phase, where they can be compared fairly.

## 2. Statistical weaknesses

* Four cell lines means four held-out folds. Any accuracy number has wide uncertainty. Bootstrap intervals are shown and described as optimistic.
* Three replicates per condition per line. DESeq2 with per-line fits and a random-effects meta-analysis across lines handles this honestly; heterogeneity (I squared) is reported for every gene.
* The persister and acute contrasts share the same control samples, so their errors are correlated. The decomposition ignores the covariance, which makes the error on S larger, not smaller. That choice is stated.
* Inferred transcription-factor activity is computed from target-gene expression, so it is not independent evidence. H3 says this.
* Literature-based ranking checks (H6) use a tiny, biased reference set and are labelled exploratory.

## 3. Engineering weaknesses

* The plan listed about forty tools. Most were never going to run. The build runs a small set end to end and tests it.
* The website could not run heavy analyses live. All site data come from a recorded pipeline run, with checksums of every input file.
* No data reach the site unless they come from the real run. Synthetic data exist only for tests and carry a red banner if ever loaded.

## 4. Presentation weaknesses

* "PERSIST-TRACE" contained a hyphen. The project name is now PersistTrace. The repository URL keeps its original slug because it already existed.
* The proposed images were third-party photos from news and company sites. They cannot be reused. Every visual on the site is generated from the data or drawn in code.
* The pitch leaned on scale ("hundreds of models"). Judges from other fields reward one clear question, one clean test and one honest surprise. The site is built around that.

## 5. What makes the build distinctive

1. **Persistence-specific signal decomposition.** P = beta A + S, estimated per line with a robust slope, then meta-analysed. Simple enough to explain on a whiteboard, and it can be checked.
2. **Transcriptional memory.** The fraction of each persister change still present after washout, plus a half-life where the data allow.
3. **Evidence passports with six separate dimensions** and written rules for the final label. No composite score.
4. **A contested case shown openly** (TROP2), with both papers and their outcomes.
5. **A preregistered benchmark** on held-out cell lines, with the hypotheses committed to the repository before the first real run.
6. **Cross-language check.** The memory fit and the switching model are implemented again in Octave/MATLAB and compared in CI.
