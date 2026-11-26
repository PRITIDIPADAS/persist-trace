# Claim boundaries and evidence standard

| Evidence | What it allows us to say |
|---|---|
| Differential expression in one line | The gene is associated with that condition in that line |
| Same direction in 3 or more of 4 lines | The association is reproducible across these models |
| S significant after removing beta A | The change is not explained by the acute response |
| Return toward baseline after washout | The population-level pattern fits reversibility |
| Programme enrichment | A gene set is over-represented; nothing about cause |
| Inferred TF activity | Candidate regulatory relationship, not independent of expression |
| Published perturbation changing a persister outcome | Functional support in the cited model and outcome |
| Published perturbation with no effect | Evidence against necessity for that outcome |
| Published studies disagreeing | Contested; both shown |

## Passport labels and their rules

Rules are applied in this order. The first that matches wins.

1. **Contested mechanism**: published perturbations of persister outcomes disagree.
2. **Perturbation-supported dependency**: a published perturbation changed a persister outcome (formation, maintenance or vulnerability) and expression is at least partly enriched.
3. **Marker without functional support**: enriched, but a perturbation found no effect on persister formation.
4. **Candidate regulator**: the gene is a transcription factor whose inferred activity changes significantly.
5. **DTP-associated marker**: reproducibly enriched and not explained by acute stress.
6. **Shared stress response**: most of the change is already present after 24 h.
7. **Insufficient evidence**: everything else.

## Words this project does not use about its own results

Cure, causes resistance, clinically validated, predicts patient outcome, first AI to. A test in the repository fails if those claims appear in the site or README without a negation.

## Scope

Cell-line models of EGFR-mutant NSCLC under osimertinib. No diagnosis, no prescribing, no patient-level prediction, no replacement for laboratory validation.
