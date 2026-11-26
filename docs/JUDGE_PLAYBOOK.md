# Judge playbook

## The 15-second version

"When lung cancer is treated with osimertinib, a few cells survive without mutating. I built a pipeline that separates their real survival programme from ordinary drug stress, measures how much of it fades when the drug stops, and labels each gene by how strong the evidence is. It runs on public data and anyone can rerun it."

## The 60-second version

1. **Problem.** Targeted drugs shrink EGFR-mutant lung tumours, but drug-tolerant persister cells survive and can seed relapse. They change state, not DNA.
2. **Gap.** Comparing treated with untreated cells mixes day-one stress with the persister state. A gene that looks important may just be stressed.
3. **Method.** 60 RNA-seq samples, four cell lines, five conditions. For every gene, the persister change is split into the part explained by acute stress (beta times A) and the persister-specific remainder S.
4. **Test.** Signatures are built on three cell lines and tested on the fourth, against a naive signature. The hypotheses were committed before the first real run.
5. **Honesty.** Every gene gets a passport. TROP2 shows why: two strong papers disagree about it, and the passport shows both.

## Two-minute live demo (works offline once the page has loaded)

| Time | Screen | Say |
|---|---|---|
| 0:00 | Hero animation | "Every dot is a real experiment. Watch the four cell lines walk from control to persister and back." |
| 0:20 | Decomposition | "The yellow line is stress. Coloured points changed in persisters for reasons stress does not explain." |
| 0:45 | Reversibility, type a gene | "This is how much the cell remembers 24 hours and about a week after the drug is gone." |
| 1:10 | Escape graph | "Solid lines are measured or published. Dashed lines are computer inference. They never look the same." |
| 1:30 | Passport, TROP2 | "A 2025 knockout said TROP2 is a passenger. A 2026 paper says it helps persisters survive. Different outcomes, different models. My system calls it contested instead of picking a side." |
| 1:50 | Validation arena | "Here is the preregistered test on an unseen cell line, including where my idea did not win." |

## Angles for judges outside cancer biology

* **Physicist or engineer:** the switching model is a 2x2 linear system; relapse is decided by the sign of the largest eigenvalue. The memory half-life is a first-order decay fit.
* **Statistician:** random-effects meta-analysis across cell lines, held-out-line validation, bootstrap intervals, robustness grid, no cell-level pseudoreplication.
* **Computer scientist:** leakage-safe splits, baselines before complex models, deterministic pipeline, CI with tests, cross-language check in Octave.
* **Clinician:** the claim boundary. Nothing here is treatment advice. Cell lines are not patients.
* **Chemist:** osimertinib binds EGFR covalently; persisters survive without the T790M or C797S style mutations, which is why a non-genetic explanation is needed.

## Questions to expect

**What is new here?** The combination: an explicit stress-versus-persistence decomposition, population-level reversibility with memory half-lives, and passports that keep markers apart from mechanisms. Persister signatures from machine learning already exist (Chang et al. 2026), so that is not the claim.

**How do you know it is not just stress?** The decomposition removes the part of each change predicted by the acute response, and H1 tests whether that helps separate persisters from acute samples in a cell line the model never saw. The answer is reported whichever way it goes.

**Does this prove a gene causes resistance?** No. Causality needs perturbation experiments. The passport only says "perturbation-supported" when a published experiment changed a persister outcome.

**Why not use a big AI model?** The core data are 60 bulk samples. Strong simple models are the honest baseline. Foundation models belong in the single-cell phase.

**What if your hypothesis failed?** Then the site says "contradicted" in red. A failed prediction is still a result.

**What would you do with a lab?** Test the top persister-specific genes that have no perturbation record yet, with knockdowns measured for both DTP formation and DTP maintenance, since TROP2 shows those can differ.

## Board layout (if a poster is required)

Left column: the hook question, one sentence on persisters, the five-condition design. Centre: the decomposition plot and the memory histogram. Right: the TROP2 passport and the H1 result. Bottom strip: QR code to the live site and the claim boundary sentence.
